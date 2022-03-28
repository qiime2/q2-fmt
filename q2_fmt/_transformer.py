# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os import sep

import pandas as pd

import qiime2
import skbio
from q2_types.distance_matrix import LSMatFormat
from q2_fmt.plugin_setup import plugin
from . import TSVFileFormat

def _tsv_format_to_dataframe(filepath, has_header=None):
    """Read a TSV file into a dataframe.

    Parameters
    ----------
    filepath : str
        The TSV file to be read.
    has_header : bool, optional
        If `None`, autodetect the header: only `TBD - PLACEHOLDER` is
        recognized, optionally followed by other columns. If `True`, the file
        must have the expected header described above otherwise an error is
        raised. If `False`, the file is read without assuming a header.

    Returns
    -------
    pd.DataFrame
        Dataframe containing parsed contents of the TSV file.

    """
    df = pd.read_csv(filepath, sep='\t', skip_blank_lines=True,
                     header=None, dtype=object)

    df.set_index(df.columns[0], drop=True, append=False, inplace=True)

    return df

@plugin.register_transformer
def _1(ff: TSVFileFormat) -> pd.DataFrame:
    return _tsv_format_to_dataframe(str(ff), has_header=None)


@plugin.register_transformer
def _2(obj: pd.DataFrame) -> TSVFileFormat:
    ff = TSVFileFormat()
    obj.to_csv(str(ff), sep='\t')
    return ff


@plugin.register_transformer
def _3(ff: LSMatFormat) -> pd.Series:
    dm = skbio.DistanceMatrix.read(str(ff), format='lsmat', verify=False)
    return dm.to_series()
