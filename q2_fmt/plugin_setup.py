# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from tokenize import String
import pandas as pd
import importlib

import qiime2.plugin
# from q2_types.sample_data import SampleData

import q2_fmt
from q2_fmt import TSVFileFormat
from q2_fmt._engraftment import dataframe_adds_blank_column

plugin = qiime2.plugin.Plugin(name='fmt',
                version=q2_fmt.__version__,
                website='https://github.com/qiime2/q2-fmt',
                package='q2_diversity',
                description='This QIIME 2 plugin supports FMT analyses.',
                short_description='Plugin for analyzing FMT data.')

plugin.register_formats(TSVFileFormat)
# plugin.register_semantic_types(ModelTests)
# plugin.register_semantic_type_to_format(
#    SampleData[ModelTests], TSVFileDirFmt)

plugin.methods.register_function(
    function=dataframe_adds_blank_column,
    inputs={'input_dataframe': pd.DataFrame},
    parameters={'column_name': str},
    outputs={'output_dataframe': pd.DataFrame},
    input_descriptions={
        'input_dataframe': ('The original dataframe to be modified.')
    },
    parameter_descriptions={
        'column_name': ('The name of the blank column to be added to the dataframe.')
    },
    output_descriptions={
        'dataframe': ('The resulting dataframe.')
    },
    name='Modifies a dataframe with a specified blank column',
    description=('This method adds a named blank column to an existing dataframe.')
)

importlib.import_module('q2_fmt._transformer')
