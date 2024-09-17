# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import itertools
import numpy as np

import qiime2


def cc(
    ctx, diversity_measure, metadata, compare, distance_to, time_column,
    subject_column, reference_column=None, control_column=None,
    filter_missing_references=False, baseline_timepoint=None, where=None,
    against_group=None, alternative='two-sided', p_val_approx='auto'
):

    raincloud_plot = ctx.get_action('stats', 'plot_rainclouds')
    group_timepoints = ctx.get_action('fmt', 'group_timepoints')

    results = []
    time_dist, ref_dist = group_timepoints(
                          diversity_measure=diversity_measure,
                          metadata=metadata, distance_to=distance_to,
                          time_column=time_column,
                          reference_column=reference_column,
                          subject_column=subject_column,
                          control_column=control_column,
                          filter_missing_references=filter_missing_references,
                          where=where,
                          baseline_timepoint=baseline_timepoint,
                          group_column=False)

    if compare == 'reference' or compare == 'all-pairwise':
        mann_whitney_u = ctx.get_action('stats', 'mann_whitney_u')
        stats = mann_whitney_u(distribution=ref_dist, compare=compare,
                               reference_group=against_group,
                               against_each=time_dist, alternative=alternative,
                               p_val_approx=p_val_approx)

    else:
        wilcoxon_srt = ctx.get_action('stats', 'wilcoxon_srt')
        stats = wilcoxon_srt(distribution=time_dist, compare=compare,
                             baseline_group=against_group,
                             alternative=alternative,
                             p_val_approx=p_val_approx)

    results += stats
    results += raincloud_plot(data=time_dist, stats=stats[0])

    return tuple(results)


def group_timepoints(
       diversity_measure: pd.Series, metadata: qiime2.Metadata,
        distance_to: str, time_column: str, reference_column: str = None,
        group_column: str = False,
        subject_column: str = False, control_column: str = None,
        filter_missing_references: bool = False,
        baseline_timepoint: str = None,
        where: str = None) -> (pd.DataFrame, pd.DataFrame):

    if isinstance(diversity_measure.index, pd.MultiIndex):
        diversity_measure.index = _sort_multi_index(diversity_measure.index)

    (is_beta, used_references, time_col, subject_col, group_col,
     used_controls) = \
        _data_filtering(diversity_measure=diversity_measure,
                        metadata=metadata, distance_to=distance_to,
                        time_column=time_column,
                        reference_column=reference_column,
                        group_column=group_column,
                        subject_column=subject_column,
                        control_column=control_column,
                        filter_missing_references=filter_missing_references,
                        baseline_timepoint=baseline_timepoint, where=where)

    original_measure_name = diversity_measure.name
    diversity_measure.name = 'measure'
    diversity_measure.index.name = 'id'

    ordered_df = _ordered_dists(diversity_measure, is_beta, used_references,
                                time_col, subject_col=subject_col,
                                group_col=group_col)

    id_annotation = {
        'title': used_references.index.name,
        'description': '...'
    }
    # id, measure, group, [subject]
    ordered_df['id'].attrs.update(id_annotation)
    ordered_df['measure'].attrs.update({
        'title': ('Distance to %s' % used_references.name)
        if is_beta else original_measure_name,
        'description': '...'
    })
    ordered_df['group'].attrs.update({
        'title': time_col.name,
        'description': '...'
    })
    if group_col is not None:
        ordered_df['class'].attrs.update({
            'title': "selected metadata column",
            'description': '...'
        })

        ordered_df['level'].attrs.update({
            'title': group_col.name,
            'description': '...'
        })
    if subject_col is not None:
        ordered_df['subject'].attrs.update({
            'title': subject_col.name,
            'description': '...'
        })

    independent_df = _independent_dists(diversity_measure, metadata,
                                        used_references, is_beta,
                                        used_controls)

    # id, measure, group, [A, B]
    if is_beta:
        independent_df['id'].attrs.update({
            'title': 'Pairwise Comparison',
            'description': 'The pairwise comparisons within a group,'
                           ' seperated by "..". Use column A and B for easier'
                           ' parsing.'
        })
    else:
        independent_df['id'].attrs.update(id_annotation)

    independent_df['measure'].attrs.update({
        'title': 'distance' if is_beta else original_measure_name,
        'description': 'Pairwise distance between A and B' if is_beta else
                       original_measure_name
    })
    independent_df['group'].attrs.update({
        'title': used_references.name if used_controls is None else
        '%s or %s' % (used_references.name, used_controls.name),
        'description': '...'
    })

    if is_beta:
        independent_df['A'].attrs.update(id_annotation)
        independent_df['B'].attrs.update(id_annotation)

    return ordered_df, independent_df


# HELPER FUNCTION FOR DATA FILTERING
def _data_filtering(diversity_measure: pd.Series, metadata: qiime2.Metadata,
                    distance_to: str, time_column: str,
                    reference_column: str = None,
                    group_column: str = None,
                    subject_column: str = False, control_column: str = None,
                    filter_missing_references: bool = False,
                    baseline_timepoint: str = None, where: str = None):

    if diversity_measure.empty:
        raise ValueError('Empty diversity measure detected.'
                         ' Please make sure your diversity measure'
                         ' contains data.')

    if isinstance(diversity_measure.index, pd.MultiIndex):
        is_beta = True
        ids_with_data = set(itertools.chain.from_iterable(
            diversity_measure.index))
    else:
        is_beta = False
        ids_with_data = diversity_measure.index

    metadata = metadata.filter_ids(ids_to_keep=ids_with_data)

    if where is not None:
        metadata = (metadata
                    .filter_ids(ids_to_keep=metadata
                                .get_ids(where=where))
                    )

    if distance_to == "donor" and baseline_timepoint is not None:
        raise ValueError("`donor` was provided to the `distance_to` parameter"
                         " and a value was provided to `baseline_timepoint`."
                         " These values can not be passed in together.")
    elif distance_to == "donor" and reference_column is None:
        raise ValueError("`donor` was provided to the `distance_to` parameter"
                         " and a `reference_column` was not provided. Please"
                         " provide a `reference_column` if you are"
                         " investigating distance to donor")
    elif distance_to == "baseline" and reference_column is not None:
        raise ValueError("`baseline` was provided to the `distance_to`"
                         " parameter and a value was provided to"
                         " `reference_column`. These values can not be passed"
                         " in together.")
    elif distance_to == "baseline" and baseline_timepoint is None:
        raise ValueError("`baseline` was provided to the `distance_to`"
                         " parameter and a `baseline_timepoint` was not"
                         " provided. Please provide a `baseline_timepoint`"
                         " if you are investigating distance to baseline")

    time_col = _get_series_from_col(
        md=metadata, col_name=time_column,
        param_name='time_column',
        expected_type=qiime2.NumericMetadataColumn)
    if distance_to == 'donor':
        reference_col = _get_series_from_col(
            md=metadata, col_name=reference_column,
            param_name='reference_column',
            expected_type=qiime2.CategoricalMetadataColumn)
        used_references = reference_col[~time_col.isna()]
    elif distance_to == 'baseline':
        used_references = \
            _get_to_baseline_ref(time_col=time_col, time_column=time_column,
                                 baseline_timepoint=baseline_timepoint,
                                 subject_column=subject_column,
                                 metadata=metadata)
    if used_references.isna().any():
        if filter_missing_references:
            used_references = used_references.dropna()
        else:
            nan_references = used_references.index[used_references.isna()]
            raise KeyError('Missing references for the associated sample data.'
                           ' Please make sure that all samples with a'
                           ' timepoint value have an associated reference.'
                           ' IDs where missing references were found:'
                           ' %s' % (tuple(nan_references),))

    available_references = (used_references.isin(ids_with_data))
    if not available_references.all():
        if filter_missing_references:
            used_references = used_references[available_references]
        else:
            raise KeyError('References included in the metadata are missing'
                           ' from the diversity measure. Please make sure all'
                           ' references included in the metadata are also'
                           ' present in the diversity measure.'
                           ' Missing references: %s'
                           % list(used_references[~available_references]
                                  .unique()))

    if used_references.empty:
        raise KeyError('No references were found within the diversity metric.')

    subject_col = None
    if subject_column:
        subject_col = _get_series_from_col(
            md=metadata, col_name=subject_column,
            param_name='subject_column',
            expected_type=qiime2.CategoricalMetadataColumn)

    group_col = None
    if group_column:
        group_col = _get_series_from_col(
            md=metadata, col_name=group_column,
            param_name='group_column',
            expected_type=qiime2.CategoricalMetadataColumn)

    used_controls = None
    if control_column is not None:
        control_col = _get_series_from_col(md=metadata,
                                           col_name=control_column,
                                           param_name='control_column')
        used_controls = control_col[~control_col.isna()]

    return (is_beta, used_references, time_col, subject_col,
            group_col, used_controls)


# HELPER FUNCTION FOR DATA Filtering
def _get_series_from_col(md, col_name, param_name, expected_type=None,
                         drop_missing_values=False):
    try:
        column = md.get_column(col_name)
    except ValueError as e:
        raise ValueError("There was an issue with the argument for %r. %s"
                         % (param_name, e)) from e

    if expected_type is not None and not isinstance(column, expected_type):
        if type(expected_type) is tuple:
            exp = tuple(e.type for e in expected_type)
        else:
            exp = expected_type.type

        raise ValueError("Provided column for %r is %r, not %r."
                         % (param_name, column.type, exp))

    if drop_missing_values:
        column = column.drop_missing_values()

    return column.to_series()


# HELPER FUNCTION FOR sorting a multi-index (for dist matrix and metadata)
def _sort_multi_index(index):
    sorted_levels = list(map(sorted, index))
    sorted_multi = pd.MultiIndex.from_tuples(sorted_levels)
    return sorted_multi


# HELPER FUNCTION FOR Dists1D[Ordered | NestedOrdered, Matched | Independent]
def _ordered_dists(diversity_measure: pd.Series, is_beta,
                   used_references, time_col, subject_col, group_col):
    if is_beta:
        idx = pd.MultiIndex.from_frame(
            used_references.to_frame().reset_index())
        idx = _sort_multi_index(idx)
        idx.names = ['id', 'reference']
        diversity_measure.index.names = ['id', 'reference']
    else:
        idx = used_references.index
        idx.name = 'id'

    try:
        sliced_df = (diversity_measure[idx]
                     .to_frame()
                     .reset_index()
                     .set_index('id')
                     )
    except KeyError:
        raise KeyError('Pairwise comparisons were unsuccessful. Please double'
                       ' check that your chosen reference column contains'
                       ' values that are also present in the ID column for'
                       ' the associated metadata.')

    if is_beta:
        sliced_df.index = used_references.index
        sliced_df.index.name = 'id'

    ordinal_df = sliced_df[['measure']]
    ordinal_df['group'] = time_col

    if subject_col is not None:
        ordinal_df['subject'] = subject_col

    if group_col is not None:
        ordinal_df['class'] = group_col.name
        ordinal_df['level'] = group_col

    return ordinal_df.reset_index()


# HELPER FUNCTION FOR GroupDists[Unordered, Independent]
def _independent_dists(diversity_measure, metadata,
                       used_references, is_beta, used_controls):
    unique_references = sorted(used_references.unique())

    if is_beta:
        if len(unique_references) > 1:
            ref_idx = pd.MultiIndex.from_tuples(
                itertools.combinations(unique_references, 2))
            ref_idx.names = ['A', 'B']

        else:
            ref_idx = pd.MultiIndex(levels=[[], []],
                                    codes=[[], []],
                                    names=['A', 'B'])

        diversity_measure.index.names = ['A', 'B']

        if used_controls is not None:
            grouped_md = (metadata
                          .to_dataframe()
                          .loc[used_controls.index]
                          .groupby(used_controls)
                          )
            ctrl_list = list()
            for group_id, grouped_ctrls in grouped_md:
                if len(grouped_ctrls.index) < 2:
                    continue
                ctrl_combos = list(
                    itertools.combinations(
                        grouped_ctrls.index, 2)
                )
                ctrl_idx = pd.MultiIndex.from_tuples(ctrl_combos)
                ctrl_series = pd.Series(group_id, index=ctrl_idx)
                ctrl_list.append(ctrl_series)

            if len(ctrl_list) >= 1:
                ctrl_series = pd.concat(ctrl_list)
                ctrl_series.index.names = ['A', 'B']

            else:
                ctrl_series = pd.Series()

            ctrl_series.name = 'group'

    else:
        ref_idx = list(unique_references)
        if used_controls is not None:
            ctrl_series = used_controls
            ctrl_series.index.name = 'id'

    try:
        nominal_df = diversity_measure[ref_idx].to_frame().reset_index()
    except KeyError:
        raise KeyError('Pairwise comparisons were unsuccessful. Please double'
                       ' check that your chosen reference column contains'
                       ' values that are also present in the ID column for'
                       ' the associated metadata.')

    nominal_df['group'] = 'reference'

    if used_controls is not None:
        ctrl_group = diversity_measure[ctrl_series.index].to_frame()
        ctrl_group['group'] = ctrl_series
        ctrl_group = ctrl_group.reset_index()
        nominal_df = pd.concat([nominal_df, ctrl_group])
        nominal_df = nominal_df.reset_index(drop=True)

    if 'A' in nominal_df.columns:
        if not nominal_df.empty:
            nominal_df['id'] = nominal_df[['A', 'B']].agg('..'.join, axis=1)
        else:
            nominal_df['id'] = []
        nominal_df = nominal_df[['id', 'measure', 'group', 'A', 'B']]

    return nominal_df


# Helper Function For to-baseline datafiltering
def _get_to_baseline_ref(time_col, baseline_timepoint, time_column,
                         subject_column, metadata):

    temp_baseline_ref = []
    reference_list = []
    baseline_ref_df = pd.DataFrame()
    # All valid FMT samples have to have a time column
    metadata = metadata.to_dataframe()[~time_col.isna()]
    if float(baseline_timepoint) not in metadata[time_column].values:
        raise AssertionError('The provided baseline timepoint'
                             f' {baseline_timepoint} was not'
                             f' found in `metadata` '
                             f' column {time_column}.')
    for sub, samples in metadata.groupby([subject_column]):
        reference = \
            samples[samples[
                time_column] == float(baseline_timepoint)].index.to_list()
        if len(reference) > 1:
            raise ValueError('More than one baseline sample was found per'
                             ' subject. Only one baseline sample can be'
                             ' used as a reference. Please group baseline'
                             ' replicates.')
        elif len(reference) == 0:
            # If there is no baseline for a subject,
            # This will either drop with filter-missing-references or
            # or error and say that they need to pass
            # filter-missing-references
            reference = [np.nan]
        temp_baseline_ref = temp_baseline_ref + samples.index.to_list()
        reference_list = \
            reference_list + (reference * len(samples.index.to_list()))
    # I dont see any way that this hits because of my above assertion but
    # I think its a good check so I am leavig it in.
    if len(reference_list) == 0:
        raise AssertionError('No baseline samples'
                             ' were found in the metadata.'
                             ' Please confirm that a valid'
                             ' baseline timepoint was given.')
    if pd.Series(reference_list).isnull().all():
        raise AssertionError('No baseline samples'
                             ' were connected via subject.'
                             ' Please confirm that all valid'
                             ' baseline timepoint where all baseline samples'
                             ' have a corresponding subject')
    baseline_ref_df['sample_name'] = temp_baseline_ref
    baseline_ref_df['relevant_baseline'] = reference_list
    baseline_ref_df = \
        baseline_ref_df[~baseline_ref_df['sample_name'].isin(
            reference_list)].set_index('sample_name')
    reference_col = _get_series_from_col(
        md=qiime2.Metadata(baseline_ref_df), col_name='relevant_baseline',
        param_name='reference_column',
        expected_type=qiime2.CategoricalMetadataColumn)
    # this is so the variables for distance to donor and distance to
    # baseline have the same variable name
    used_references = reference_col
    return used_references
