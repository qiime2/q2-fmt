# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import frictionless as fls
from frictionless import Resource

from q2_fmt.plugin_setup import plugin
from ._format import (NDJSONFileFormat, DataResourceSchemaFileFormat,
                      TabularDataResourceDirFmt)


@plugin.register_transformer
def _1(obj: pd.DataFrame) -> NDJSONFileFormat:
    ff = NDJSONFileFormat()
    obj.to_json(str(ff), lines=True, orient='records')
    return ff


@plugin.register_transformer
def _2(obj: DataResourceSchemaFileFormat) -> fls.Resource:
    return fls.Resource(str(obj))


@plugin.register_transformer
def _3(df: TabularDataResourceDirFmt) -> pd.DataFrame:
    path = df.data.view(NDJSONFileFormat)
    data = pd.read_json(str(path), lines=True)
    resource = df.metadata.view(fls.Resource)

    for field in resource.schema.fields:
        data[field['name']].attrs = field.to_dict()

    return data


@plugin.register_transformer
def _4(obj: pd.DataFrame) -> TabularDataResourceDirFmt:
    metadata_resource = Resource()

    for col in obj.columns:
        series = obj[col]
        dtype = series.convert_dtypes().dtype
        metadata = series.attrs.copy()

        if pd.api.types.is_float_dtype(dtype):
            schema_dtype = 'number'
        elif pd.api.types.is_integer_dtype(dtype):
            schema_dtype = 'integer'
        else:
            schema_dtype = 'string'

        metadata['name'] = col
        metadata['type'] = schema_dtype
        metadata_resource.schema.add_field(source=metadata)

    metadata_resource.data = 'data.ndjson'
    metadata_resource.format = 'ndjson'

    dir_fmt = TabularDataResourceDirFmt()

    dir_fmt.data.write_data(obj, pd.DataFrame)
    metadata_resource.to_json(str(dir_fmt.path/'dataresource.json'))

    return dir_fmt
