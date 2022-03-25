# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import itertools
import typing
import skbio

import qiime2

# Dataframe does something basic
def _dataframe_adds_blank_column(input_dataframe, column_name):
    input_dataframe[column_name] = ''

    return input_dataframe

def dataframe_adds_blank_column(dataframe: pd.DataFrame, column_name: str) -> (pd.DataFrame):
    return _dataframe_adds_blank_column(input_dataframe=dataframe, column_name=column_name)

# Group timepoints - combines beta div & metadata based on groups & specified metric
def _group_timepoints(dist_matrix, metadata, group_column):
    # Instantiating the output dataframe that contains the distances & their associated group metric
    output_df = pd.DataFrame(
        columns=['distance_measure', '{}'.format(group_column)]
    )
    # Creating the grouped metadata based on the specified group metric
    grouped_md = metadata.groupby(metadata[group_column])

    # Creating all possible sample pairings for the distance measures within each group
    for group_value, grouped_samples in grouped_md:
        group_combos = itertools.combinations(grouped_samples.index, 2)

        # Iterating through each sample pairing
        for combo in group_combos:
            # Constructing a DF that contains each sample pairing & associated group value
            row = pd.DataFrame([
                [dist_matrix[combo[0]][combo[1]], group_value]
                ],
                columns=['distance measure', '{}'.format(group_column)],
                index=[frozenset(combo)],
                )
            # Appending each row that contains the sample pairing & group value to the output DF
            output_df = output_df.append(row)

# def group_timepoints(dist_matrix: pd.DataFrame, metadata: pd.DataFrame, column_name: str) -> (pd.DataFrame):
#     return _group_timepoints(dist_matrix=dist_matrix, metadata=metadata, column_name=column_name)


def group_timepoints(
        diversity_measure: typing.Union[skbio.DistanceMatrix, pd.Series],
        metadata: qiime2.Metadata, time_column: str, reference_column: str,
        subject_column: str = None, control_column: str = None
        ) -> (pd.DataFrame, pd.DataFrame):

    # Data filtering
    if isinstance(diversity_measure, skbio.DistanceMatrix):
        is_beta = True
        ids_with_data = list(diversity_measure.ids)
        diversity_measure = diversity_measure.to_series()
    else:
        is_beta = False
        ids_with_data = list(diversity_measure.index)

    metadata = metadata.filter_ids(ids_to_keep=ids_with_data)
    time_col = metadata.get_column(time_column).to_series()
    reference_col = metadata.get_column(reference_column).to_series()
    if subject_column is not None:
        subject_col = metadata.get_column(subject_column).to_series()
    if control_column is not None:
        control_col = metadata.get_column(control_column).to_series()

    used_references = reference_col[~time_col.isna()]

    # GroupDists[Ordinal]
    diversity_measure.name = 'measure'
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

    # GroupDists[Nominal]
    nominal_df = pd.DataFrame()

    return ordinal_df, nominal_df
