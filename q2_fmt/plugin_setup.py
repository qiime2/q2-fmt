# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

import qiime2.plugin
from q2_types.sample_data import SampleData

import q2_fmt
from q2_fmt import TSVFileFormat, ModelTests
from q2_fmt._engraftment import dataframe_adds_blank_column
from q2_fmt._format import TSVFileDirFmt
from q2_fmt._visualizer import hello_world
from q2_fmt._type import GroupDist, Ordinal, Nominal

plugin = qiime2.plugin.Plugin(name='fmt',
                version=q2_fmt.__version__,
                website='https://github.com/qiime2/q2-fmt',
                package='q2_diversity',
                description='This QIIME 2 plugin supports FMT analyses.',
                short_description='Plugin for analyzing FMT data.')

plugin.register_formats(TSVFileFormat, TSVFileDirFmt)
plugin.register_semantic_types(ModelTests, GroupDist, Ordinal, Nominal)
plugin.register_semantic_type_to_format(
    SampleData[ModelTests], TSVFileDirFmt, GroupDist[Ordinal | Nominal])

plugin.methods.register_function(
    function=dataframe_adds_blank_column,
    inputs={'dataframe': SampleData[ModelTests]},
    parameters={'column_name': qiime2.plugin.Str},
    outputs=[('output_dataframe', SampleData[ModelTests])],
    input_descriptions={
        'dataframe': 'The original dataframe to be modified.'
    },
    parameter_descriptions={
        'column_name': 'The name of the blank column to be added to the dataframe.'
    },
    output_descriptions={
        'output_dataframe': 'The resulting dataframe.'
    },
    name='Modifies a dataframe with a specified blank column',
    description='This method adds a named blank column to an existing dataframe.'
)

plugin.visualizers.register_function(
    function=hello_world,
    inputs={},
    parameters={},
    input_descriptions={},
    parameter_descriptions={},
    name='Hello World Viz',
    description='Placeholder visualizer that outputs hello world.'
)

importlib.import_module('q2_fmt._transformer')
