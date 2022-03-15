# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

# Dataframe does something basic
def _dataframe_adds_blank_column(input_dataframe, column_name):
    input_dataframe[column_name] = ''

    return input_dataframe

def dataframe_adds_blank_column(dataframe: pd.DataFrame, column_name: str) -> (pd.DataFrame):
    return _dataframe_adds_blank_column(input_dataframe=dataframe, column_name=column_name)

# Group timepoints - combines beta div & metadata based on groups & specified metric
