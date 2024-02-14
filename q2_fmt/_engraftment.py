# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import itertools

import qiime2


def engraftment(
    ctx, diversity_measure, metadata, compare, time_column,
    reference_column, subject_column, control_column=None,
    filter_missing_references=False, where=None, against_group=None,
    alternative='two-sided', p_val_approx='auto'
):

    raincloud_plot = ctx.get_action('stats', 'plot_rainclouds')
    group_timepoints = ctx.get_action('fmt', 'group_timepoints')

    results = []

    time_dist, ref_dist = group_timepoints(diversity_measure, metadata,
                                           time_column, reference_column,
                                           subject_column, control_column,
                                           filter_missing_references, where)

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
        time_column: str, reference_column: str, subject_column: str = False,
        control_column: str = None, filter_missing_references: bool = False,
        where: str = None) -> (pd.DataFrame, pd.DataFrame):

    if isinstance(diversity_measure.index, pd.MultiIndex):
        diversity_measure.index = _sort_multi_index(diversity_measure.index)

    is_beta, used_references, time_col, subject_col, used_controls = \
        _data_filtering(diversity_measure, metadata, time_column,
                        reference_column, subject_column, control_column,
                        filter_missing_references, where)

    original_measure_name = diversity_measure.name
    diversity_measure.name = 'measure'
    diversity_measure.index.name = 'id'
             
    ordered_df = _ordered_dists(diversity_measure, is_beta, used_references,
                                time_col, subject_col)

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


# TODO: It seems like if I am baseline or donor is reference I would need to
# change code here But I think that makes the most sense and seems easy! 
# HELPER FUNCTION FOR DATA FILTERING
def _data_filtering(diversity_measure: pd.Series, metadata: qiime2.Metadata,
                    time_column: str, reference_column: str, distance_to: str,
                    subject_column: str = False, control_column: str = None,
                    filter_missing_references: bool = False,
                    where: str = None, baseline_timepoint: int = None):
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
        raise ValueError("'donor' was provided to the distance_to parameter"
                         " and a value was provided to baseline_timepoint."
                         " These values can not be passed in together.")
    elif distance_to == "donor" and reference_column is None:
        raise ValueError("'donor' was provided to the distance_to parameter"
                         " and a reference_column was not provided. Please"
                         " provide a reference_column if you are investigating"
                         " distance to donor")
    elif distance_to == "baseline" and reference_column is not None:
        raise ValueError("'baseline' was provide to the distance_to parameter"
                         " and a value was provided to reference_column."
                         " These values can not be passed in together.")
    elif distance_to == "baseline" and baseline_timepoint is None: 
        raise ValueError("'baseline' was provided to the distance_to parameter"
                         " and a baseline_timepoint was not provided. Please"
                         " provide a baseline_timepoint if you are" 
                         " investigating distance to baseline")
    elif distance_to == "baseline" and baseline_timepoint is not None: 
        ##todo I need to get the baseline reference into a series 

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

    time_col = _get_series_from_col(
        md=metadata, col_name=time_column,
        param_name='time_column',
        expected_type=qiime2.NumericMetadataColumn)
    if distance_to == 'donor':
        reference_col = _get_series_from_col(
            md=metadata, col_name=reference_column,
            param_name='reference_column',
            expected_type=qiime2.CategoricalMetadataColumn)

    # TODO: only applicable to donor as the reference column 
    used_references = reference_col[~time_col.isna()]

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

    used_controls = None
    if control_column is not None:
        control_col = _get_series_from_col(md=metadata,
                                           col_name=control_column,
                                           param_name='control_column')
        used_controls = control_col[~control_col.isna()]

    return is_beta, used_references, time_col, subject_col, used_controls


# HELPER FUNCTION FOR sorting a multi-index (for dist matrix and metadata)
def _sort_multi_index(index):
    sorted_levels = list(map(sorted, index))
    sorted_multi = pd.MultiIndex.from_tuples(sorted_levels)
    return sorted_multi


# HELPER FUNCTION FOR GroupDists[Ordered, Matched | Independent]
def _ordered_dists(diversity_measure: pd.Series, is_beta,
                   used_references, time_col, subject_col):
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
