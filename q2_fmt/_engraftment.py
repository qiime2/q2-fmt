# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
# import qiime2

from q2_fmt._format import TSVFileFormat

def _tsv_to_dataframe(tsv_filepath):
    df = pd.read_csv(tsv_filepath, sep='\t', skip_blank_lines=True,
                     header=None, dtype=object)

    df.set_index(df.columns[0], drop=True, append=False, inplace=True)

    return df

def tsv_to_dataframe(tsv_filepath: TSVFileFormat) -> (pd.DataFrame):
    return _tsv_to_dataframe(tsv_filepath=tsv_filepath)
