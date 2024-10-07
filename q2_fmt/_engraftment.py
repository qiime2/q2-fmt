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

from q2_fmt._util import (_get_to_baseline_ref, _get_series_from_col,
                          _sort_multi_index, _ordered_dists,
                          _independent_dists)


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
            _get_to_baseline_ref(time_col=time_col,
                                 time_column_name=time_column,
                                 baseline_timepoint=baseline_timepoint,
                                 subject_column_name=subject_column,
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
