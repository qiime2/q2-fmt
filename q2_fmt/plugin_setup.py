# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

import qiime2.plugin
# from q2_types.sample_data import SampleData

import q2_fmt
from q2_fmt import TSVFileFormat

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
    function=q2_fmt.tsv_to_dataframe,
    inputs={'tsv_filepath': TSVFileFormat},
    parameters={},
    outputs={'dataframe': pd.DataFrame},
    input_descriptions={
        'tsv_filepath': ('The TSV file to be transformed into a dataframe.')
    },
    output_descriptions={
        'dataframe': ('The resulting dataframe.')
    },
    name='Transform a TSV file into a dataframe',
    description=('This method transforms a TSV file into a dataframe.')
)
