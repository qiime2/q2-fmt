# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from q2_fmt._format import TSVFileFormat

# TSV to dataframe
def _tsv_to_dataframe(tsv_filepath):
    df = pd.read_csv(tsv_filepath, sep='\t', skip_blank_lines=True,
                     header=None, dtype=object)

    df.set_index(df.columns[0], drop=True, append=False, inplace=True)

    return df

def tsv_to_dataframe(tsv_filepath: TSVFileFormat) -> (pd.DataFrame):
    return _tsv_to_dataframe(tsv_filepath=tsv_filepath)

# Dataframe does something basic
def _dataframe_adds_blank_column(input_dataframe, column_name):
    input_dataframe = tsv_to_dataframe(tsv_filepath=TSVFileFormat)
    input_dataframe[column_name] = ''

    return input_dataframe

def dataframe_adds_blank_column(dataframe: pd.DataFrame) -> (pd.DataFrame):
    return _dataframe_adds_blank_column(input_dataframe=dataframe)
