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


def group_timepoints(
        diversity_measure: pd.Series, metadata: qiime2.Metadata,
        time_column: str, reference_column: str, subject_column: str = None,
        control_column: str = None) -> (pd.DataFrame, pd.DataFrame):

    # Data filtering
    if isinstance(diversity_measure.index, pd.MultiIndex):
        is_beta = True
        ids_with_data = set(itertools.chain.from_iterable(
            diversity_measure.index))
    else:
        is_beta = False
        ids_with_data = diversity_measure.index

    metadata = metadata.filter_ids(ids_to_keep=ids_with_data)

    try:
        time_col = metadata.get_column(time_column)
    except:
        raise ValueError('time_column provided not present within the metadata')

    if not isinstance(time_col, qiime2.NumericMetadataColumn):
        raise TypeError('Non-numeric characters detected in time_column.')
    else:
        time_col = time_col.to_series()

    try:
        reference_col = metadata.get_column(reference_column)
    except:
        raise ValueError('reference_column provided not present within the metadata')
    reference_col = reference_col.to_series()
    used_references = reference_col[~time_col.isna()]

    if subject_column is not None:
        try:
            subject_col = metadata.get_column(subject_column)
        except:
            raise ValueError('subject_column provided not present within the metadata')
        subject_col = subject_col.to_series()

    if control_column is not None:
        try:
            control_col = metadata.get_column(control_column)
        except:
            raise ValueError('control_column provided not present within the metadata')
        control_col = control_col.to_series()
        used_controls = control_col[~control_col.isna()]

    diversity_measure.name = 'measure'
    diversity_measure.index.name = 'id'

    # GroupDists[Gordinal]
    if is_beta:
        idx = pd.MultiIndex.from_frame(
            used_references.to_frame().reset_index())
        idx.names = ['id', 'reference']
    else:
        idx = used_references.index
        idx.name = 'id'

    sliced_df = diversity_measure[idx].to_frame().reset_index().set_index('id')
    ordinal_df = sliced_df[['measure']]
    ordinal_df['group'] = time_col
    if subject_column is not None:
        ordinal_df['subject'] = subject_col

    # GroupDists[Gnominal]
    unique_references = used_references.unique()

    if is_beta:
        ref_idx = pd.MultiIndex.from_tuples(
            itertools.combinations(unique_references, 2))
        ref_idx.names = ['A', 'B']

        if control_column is not None:
            grouped_md = metadata.to_dataframe().loc[used_controls.index].groupby(used_controls)
            ctrl_list = list()
            for group_id, grouped_ctrls in grouped_md:
                if len(grouped_ctrls.index) < 2:
                    continue
                ctrl_combos = list(itertools.combinations(grouped_ctrls.index, 2))
                ctrl_idx = pd.MultiIndex.from_tuples(ctrl_combos)
                ctrl_series = pd.Series(group_id, index=ctrl_idx)
                ctrl_list.append(ctrl_series)

            ctrl_series = pd.concat(ctrl_list)
            ctrl_series.name = 'group'
            ctrl_series.index.names = ['A', 'B']

    else:
        ref_idx = list(unique_references)
        if control_column is not None:
            ctrl_series = used_controls
            ctrl_series.index.name = 'id'

    nominal_df = diversity_measure[ref_idx].to_frame().reset_index()
    nominal_df['group'] = 'reference'

    if control_column is not None:
        ctrl_group = diversity_measure[ctrl_series.index].to_frame()
        ctrl_group['group'] = ctrl_series
        ctrl_group = ctrl_group.reset_index()
        nominal_df = pd.concat([nominal_df, ctrl_group])
        nominal_df = nominal_df.reset_index(drop=True)

    if 'A' in nominal_df.columns:
        nominal_df['id'] = nominal_df[['A', 'B']].agg('..'.join, axis=1)
        nominal_df = nominal_df[['id', 'measure', 'group', 'A', 'B']]

    nominal_df = nominal_df.set_index('id')

    return ordinal_df, nominal_df
