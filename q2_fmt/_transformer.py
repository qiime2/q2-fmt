# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

import skbio
from q2_types.distance_matrix import LSMatFormat
from q2_fmt.plugin_setup import plugin
from . import RecordTSVFileFormat, AnnotatedTSVDirFmt


@plugin.register_transformer
def _1(ff: RecordTSVFileFormat) -> pd.DataFrame:
    df = pd.read_csv(str(ff), sep='\t', skip_blank_lines=True, header=0)
    return df


@plugin.register_transformer
def _2(obj: pd.DataFrame) -> RecordTSVFileFormat:
    ff = RecordTSVFileFormat()
    obj.to_csv(str(ff), sep='\t', index=False)
    return ff


@plugin.register_transformer
def _3(ff: LSMatFormat) -> pd.Series:
    dm = skbio.DistanceMatrix.read(str(ff), format='lsmat', verify=False)
    return dm.to_series()

@plugin.register_transformer
def _4(df: AnnotatedTSVDirFmt) -> pd.DataFrame:
    data = df.data.view(pd.DataFrame)
    metadata = df.metadata.view(pd.DataFrame)
    metadata = metadata.set_index('column')

    for column in data.columns:
        # not sure what the semantics are, so do our best
        data[column].attrs.update(metadata.loc[column].to_dict())

    return data

@plugin.register_transformer
def _5(obj: pd.DataFrame) -> AnnotatedTSVDirFmt:
    metadata = []
    for col in obj.columns:
        metadata.append(obj[col].attrs)

    metadata_df = pd.DataFrame(metadata, index=obj.columns.copy())
    metadata_df.index.name = 'column'
    metadata_df = metadata_df.reset_index()

    dir_fmt = AnnotatedTSVDirFmt()

    dir_fmt.data.write_data(obj, pd.DataFrame)
    dir_fmt.metadata.write_data(metadata_df, pd.DataFrame)

    return dir_fmt
