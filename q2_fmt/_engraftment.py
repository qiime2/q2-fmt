# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import itertools

import qiime2
from q2_fmt._stats import mann_whitney_u, wilcoxon_srt
from q2_fmt._visualizer import plot_rainclouds

#TODO: add in control comparison hypothesis
# Questions from Greg's study
# Increase in alpha div from baseline
# Difference in alpha at a single timepoint (week vs. control)
# Beta: is distance to donor significantly *lower* (not just different) than it was at baseline (week 0, pre-FMT)
# Beta[donor, subject]_w < Beta[donor, subject]_0


# takes the same inputs as group_timepoints and then runs the corresponding outputs
# (either time dist or ref dist) through either wilcoxon or mann_whitney and then
# through plot_rainclouds - final output includes stats.qza and raincloud_plot.qzv
def engraftment(
    diversity_measure: pd.Series, metadata: qiime2.Metadata,
    hypothesis: str, time_column: str, reference_column: str,
    subject_column: str = False, control_column: str = None,
    filter_missing_references: bool = False, where: str = None,
    baseline_group: str=None, reference_group: str=None,
    against_each: pd.DataFrame=None,
    p_val_approx: str='auto') -> (pd.DataFrame):
    # should the function return (pd.DataFrame, None) for the viz?

    time_dist, ref_dist = group_timepoints(diversity_measure, metadata,
    time_column, reference_column, subject_column, control_column,
    filter_missing_references, where)

    if hypothesis == 'reference' or 'all-pairwise':
        stats = mann_whitney_u(distribution=ref_dist, hypothesis=hypothesis,
                               reference_group=reference_group,
                               against_each=against_each, p_val_approx=p_val_approx)

    if hypothesis == 'baseline' or 'consecutive':
        stats = wilcoxon_srt(distribution=time_dist, hypothesis=hypothesis,
                             baseline_group=baseline_group, p_val_approx=p_val_approx)

    viz = plot_rainclouds(data=stats)

    return stats, viz

def group_timepoints(
        diversity_measure: pd.Series, metadata: qiime2.Metadata,
        time_column: str, reference_column: str, subject_column: str = False,
        control_column: str = None, filter_missing_references: bool = False,
        where: str = None) -> (pd.DataFrame, pd.DataFrame):

    if isinstance(diversity_measure.index, pd.MultiIndex):
        diversity_measure.index = _sort_multi_index(diversity_measure.index)

    is_beta, used_references, time_col, subject_col, used_controls = \
        _data_filtering(diversity_measure, metadata, time_column, reference_column,
                        subject_column, control_column, filter_missing_references,
                        where)

    original_measure_name = diversity_measure.name
    diversity_measure.name = 'measure'
    diversity_measure.index.name = 'id'

    ordered_df = _ordered_dists(diversity_measure, is_beta, used_references,
                                time_col, subject_col)

    id_annotation = {
        'unit': used_references.index.name,
        'description': '...'
    }
    # id, measure, group, [subject]
    ordered_df['id'].attrs.update(id_annotation)
    ordered_df['measure'].attrs.update({
        'unit': ('Distance to %s' % used_references.name)
                if is_beta else original_measure_name,
        'description': '...'
    })
    ordered_df['group'].attrs.update({
        'unit': time_col.name,
        'description': '...'
    })
    if subject_col is not None:
        ordered_df['subject'].attrs.update({
            'unit': subject_col.name,
            'description': '...'
        })


    independent_df = _independent_dists(diversity_measure, metadata,
                                        used_references, is_beta, used_controls)

    # id, measure, group, [A, B]
    if is_beta:
        independent_df['id'].attrs.update({
            'unit': 'Pairwise Comparison',
            'description': 'The pairwise comparisons within a group,'
                           ' seperated by "..". Use column A and B for easier'
                           ' parsing.'
        })
    else:
        independent_df['id'].attrs.update(id_annotation)

    independent_df['measure'].attrs.update({
        'unit': 'distance' if is_beta else original_measure_name,
        'description': 'Pairwise distance between A and B' if is_beta else
                       original_measure_name
    })
    independent_df['group'].attrs.update({
        'unit': used_references.name if used_controls is None else
                '%s or %s' % (used_references.name, used_controls.name),
        'description': '...'
    })
    if is_beta:
        independent_df['A'].attrs.update(id_annotation)
        independent_df['B'].attrs.update(id_annotation)

    return ordered_df, independent_df

# HELPER FUNCTION FOR DATA FILTERING
def _data_filtering(diversity_measure: pd.Series, metadata: qiime2.Metadata,
        time_column: str, reference_column: str, subject_column: str = False,
        control_column: str = None, filter_missing_references: bool = False,
        where: str = None):

    if diversity_measure.empty:
        raise ValueError('Empty diversity measure detected.'
                         ' Please make sure your diversity measure contains data.')

    if isinstance(diversity_measure.index, pd.MultiIndex):
        is_beta = True
        ids_with_data = set(itertools.chain.from_iterable(
            diversity_measure.index))
    else:
        is_beta = False
        ids_with_data = diversity_measure.index

    metadata = metadata.filter_ids(ids_to_keep=ids_with_data)

    if where is not None:
        metadata = metadata.filter_ids(ids_to_keep=metadata.get_ids(where=where))

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

    time_col = _get_series_from_col(md=metadata, col_name=time_column, param_name='time_column',
                                    expected_type=qiime2.NumericMetadataColumn)

    reference_col = _get_series_from_col(md=metadata, col_name=reference_column, param_name='reference_column',
                                         expected_type=qiime2.CategoricalMetadataColumn)

    used_references = reference_col[~time_col.isna()]

    if used_references.isna().any():
        if filter_missing_references:
            used_references = used_references.dropna()
        else:
            nan_references = used_references.index[used_references.isna()]
            raise KeyError('Missing references for the associated sample data. Please make sure'
                        ' that all samples with a timepoint value have an associated reference.'
                        ' IDs where missing references were found: %s' % (tuple(nan_references),))

    available_references = (used_references.isin(ids_with_data))
    if not available_references.all():
        if filter_missing_references:
            used_references = used_references[available_references]
        else:
            raise KeyError('References included in the metadata are missing from the diversity measure.'
                        ' Please make sure all references included in the metadata are also present'
                        ' in the diversity measure. Missing references: %s'
                        % list(used_references[~available_references].unique())
            )

    if used_references.empty:
        raise KeyError('No references were found within the diversity metric.')

    subject_col = None
    if subject_column:
            subject_col = _get_series_from_col(md=metadata, col_name=subject_column, param_name='subject_column',
                                               expected_type=qiime2.CategoricalMetadataColumn)

    used_controls = None
    if control_column is not None:
        control_col = _get_series_from_col(md=metadata, col_name=control_column, param_name='control_column')
        used_controls = control_col[~control_col.isna()]

    return is_beta, used_references, time_col, subject_col, used_controls

# HELPER FUNCTION FOR sorting a multi-index (for dist matrix and metadata)
def _sort_multi_index(index):
    sorted_levels = list(map(sorted, index))
    sorted_multi = pd.MultiIndex.from_tuples(sorted_levels)
    return sorted_multi

# HELPER FUNCTION FOR GroupDists[Ordered, Matched | Independent]
def _ordered_dists(diversity_measure: pd.Series, is_beta, used_references, time_col, subject_col):
    if is_beta:
        idx = pd.MultiIndex.from_frame(
            used_references.to_frame().reset_index())
        idx = _sort_multi_index(idx)
        idx.names = ['id', 'reference']
    else:
        idx = used_references.index
        idx.name = 'id'

    try:
        sliced_df = diversity_measure[idx].to_frame().reset_index().set_index('id')
    except KeyError:
        raise KeyError('Pairwise comparisons were unsuccessful. Please double check that your'
        ' chosen reference column contains values that are also present in the ID column for'
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
def _independent_dists(diversity_measure, metadata, used_references, is_beta, used_controls):
    unique_references = sorted(used_references.unique())

    if is_beta:
        try:
            ref_idx = pd.MultiIndex.from_tuples(
                itertools.combinations(unique_references, 2))
        except TypeError:
            raise TypeError('Single reference value detected. More than one unique reference must be'
                            ' provided for successful grouping.')

        ref_idx.names = ['A', 'B']

        if used_controls is not None:
            grouped_md = metadata.to_dataframe().loc[used_controls.index].groupby(used_controls)
            ctrl_list = list()
            for group_id, grouped_ctrls in grouped_md:
                if len(grouped_ctrls.index) < 2:
                    continue
                ctrl_combos = list(itertools.combinations(grouped_ctrls.index, 2))
                ctrl_idx = pd.MultiIndex.from_tuples(ctrl_combos)
                ctrl_series = pd.Series(group_id, index=ctrl_idx)
                ctrl_list.append(ctrl_series)

            try:
                ctrl_series = pd.concat(ctrl_list)
            except ValueError:
                raise ValueError('One or less controls detected. When including controls in your data,'
                                 ' please include more than one for successful grouping.')

            ctrl_series.name = 'group'
            ctrl_series.index.names = ['A', 'B']

    else:
        ref_idx = list(unique_references)
        if used_controls is not None:
            ctrl_series = used_controls
            ctrl_series.index.name = 'id'

    try:
        nominal_df = diversity_measure[ref_idx].to_frame().reset_index()
    except KeyError:
        raise KeyError('Pairwise comparisons were unsuccessful. Please double check that your'
        ' chosen reference column contains values that are also present in the ID column for'
        ' the associated metadata.')

    nominal_df['group'] = 'reference'

    if used_controls is not None:
        ctrl_group = diversity_measure[ctrl_series.index].to_frame()
        ctrl_group['group'] = ctrl_series
        ctrl_group = ctrl_group.reset_index()
        nominal_df = pd.concat([nominal_df, ctrl_group])
        nominal_df = nominal_df.reset_index(drop=True)

    if 'A' in nominal_df.columns:
        nominal_df['id'] = nominal_df[['A', 'B']].agg('..'.join, axis=1)
        nominal_df = nominal_df[['id', 'measure', 'group', 'A', 'B']]

    return nominal_df
