# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

from qiime2.plugin import Str, Plugin, Metadata, TypeMap, Bool, Choices
from q2_types.sample_data import SampleData, AlphaDiversity
from q2_types.distance_matrix import DistanceMatrix

import q2_fmt
from q2_fmt import TSVFileFormat, ModelTests
from q2_fmt._engraftment import group_timepoints
from q2_fmt._stats import mann_whitney_u, wilcoxon_srt
from q2_fmt._format import TSVFileDirFmt
from q2_fmt._visualizer import hello_world
from q2_fmt._type import GroupDist, Matched, Independent, Ordered, Unordered

plugin = Plugin(name='fmt',
                version=q2_fmt.__version__,
                website='https://github.com/qiime2/q2-fmt',
                package='q2_fmt',
                description='This QIIME 2 plugin supports FMT analyses.',
                short_description='Plugin for analyzing FMT data.')

plugin.register_formats(TSVFileFormat, TSVFileDirFmt)
plugin.register_semantic_types(ModelTests, GroupDist, Matched, Independent,
                               Ordered, Unordered)
plugin.register_semantic_type_to_format(
    GroupDist[Ordered | Unordered, Matched | Independent], TSVFileDirFmt)

T_subject, T_dependence = TypeMap({
    Bool % Choices(False): Independent,
    Str: Matched
})

plugin.methods.register_function(
    function=group_timepoints,
    inputs={'diversity_measure': DistanceMatrix | SampleData[AlphaDiversity] },
    parameters={'metadata': Metadata, 'time_column': Str,
                'reference_column': Str, 'subject_column': T_subject, 'control_column': Str},
    outputs=[('timepoint_dists', GroupDist[Ordered, T_dependence]),
             ('reference_dists', GroupDist[Unordered, Independent])],
    parameter_descriptions={
        'metadata': 'The sample metadata.',
        'time_column': 'The column within the `metadata` that the `diversity_measure` should be grouped by.'
                       ' This column should contain simple integer values.',
        'control_column': 'The column within the `metadata` that contains any relevant control group IDs.'
                          ' Actual treatment samples should not contain any value within this column.',
        'reference_column': 'The column within the `metadata` that contains the sample to use as a reference'
                            ' for a given beta `diversity_measure`.'
                            ' For example, this may be the relevant donor sample to compare against.',
        'subject_column': 'The column within the `metadata` that contains the subject ID to be tracked against timepoints.',
    },
    output_descriptions={
        'timepoint_dists': 'The distributions for the `diversity_measure`, grouped by the selected `time_column`.'
                           ' May also contain subject IDs, if `subject_column` is provided in the `metadata`.',
        'reference_dists': 'The inter-group reference and inter-group control (when provided) distributions.'
                           ' When `diversity_measure` is DistanceMatrix, the inter-group calculations will be all pairwise comparisons within a group.'
                           ' Otherwise, these are just the per-sample measurements of alpha-diversity.'
    },
    name='',
    description=''
)

plugin.methods.register_function(
    function=mann_whitney_u,
    inputs={'distribution': GroupDist[Unordered | Ordered, Independent],
            'against_each': GroupDist[Unordered | Ordered, Matched | Independent]},
    parameters={'hypothesis': Str % Choices('reference', 'all-pairwise'),
                'reference_group': Str},
    outputs=[('stats', ModelTests)],
    parameter_descriptions={
        'hypothesis': 'The hypothesis that will be used to analyze the input `distribution`.'
                      ' Either `reference` or `all-pairwise` must be selected.',
        'reference_group': 'If `reference` is the selected hypothesis, this is the column that will be used'
                           ' to compare all samples against.',
    },
    output_descriptions={
        'stats': 'The Mann-Whitney U distribution for either the `reference` or `all-pairwise` hypothesis.',
    },
    name='Mann-Whitney U Test',
    description=''
)

plugin.methods.register_function(
    function=wilcoxon_srt,
    inputs={'distribution': GroupDist[Ordered, Matched]},
    parameters={'hypothesis': Str % Choices('baseline', 'consecutive'),
                'baseline_group': Str},
    outputs=[('stats', ModelTests)],
    parameter_descriptions={
        'hypothesis': 'The hypothesis that will be used to analyze the input `distribution`.'
                      ' Either `baseline` or `consecutive` must be selected.',
        'baseline_group': 'If `baseline` is the selected hypothesis, this is the column that will be used'
                          ' to compare all samples against.',
    },
    output_descriptions={
        'stats': 'The Wilcoxon SRT distribution for either the `baseline` or `consecutive` hypothesis.',
    },
    name='Wilcoxon Signed Rank Test',
    description=''
)

# Dummy Visualizer
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
