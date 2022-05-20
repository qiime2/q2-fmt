# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

from qiime2.plugin import (Str, Plugin, Metadata, TypeMap,
                           Bool, Choices, Visualization)
from q2_types.sample_data import SampleData, AlphaDiversity
from q2_types.distance_matrix import DistanceMatrix

import q2_fmt
from q2_fmt import RecordTSVFileFormat
from q2_fmt._engraftment import group_timepoints
from q2_fmt._stats import mann_whitney_u, wilcoxon_srt
from q2_fmt._format import AnnotatedTSVDirFmt
from q2_fmt._visualizer import plot_rainclouds
from q2_fmt._type import (GroupDist, Matched, Independent, Ordered,
                          Unordered, StatsTable, Pairwise)
import q2_fmt._examples as ex

plugin = Plugin(name='fmt',
                version=q2_fmt.__version__,
                website='https://github.com/qiime2/q2-fmt',
                package='q2_fmt',
                description='This QIIME 2 plugin supports FMT analyses.',
                short_description='Plugin for analyzing FMT data.')

plugin.register_formats(RecordTSVFileFormat, AnnotatedTSVDirFmt)
plugin.register_semantic_types(StatsTable, Pairwise, GroupDist, Matched,
                               Independent, Ordered, Unordered)
plugin.register_semantic_type_to_format(
    GroupDist[Ordered | Unordered,
              Matched | Independent] | StatsTable[Pairwise], AnnotatedTSVDirFmt
    )

T_subject, T_dependence = TypeMap({
    Bool % Choices(False): Independent,
    Str: Matched
})

T_engraft_subject, T_compare, _ = TypeMap({
    (Bool % Choices(False),
        Str % Choices("reference", "all-pairwise")): Visualization,
    (Str, Str % Choices("baseline", "consecutive")): Visualization
})

plugin.pipelines.register_function(
    function=q2_fmt.engraftment,
    inputs={'diversity_measure': DistanceMatrix | SampleData[AlphaDiversity]},
    parameters={'metadata': Metadata,
                'compare': T_compare,
                'time_column': Str, 'reference_column': Str,
                'subject_column': T_engraft_subject, 'control_column': Str,
                'filter_missing_references': Bool, 'where': Str,
                'against_group': Str,
                'alternative': Str % Choices('two-sided', 'greater', 'less'),
                'p_val_approx': Str % Choices('auto', 'exact', 'asymptotic')},
    outputs=[
        ('stats', StatsTable[Pairwise]),
        ('raincloud_plot', Visualization)
    ],
    input_descriptions={
        'diversity_measure': '',
    },
    parameter_descriptions={
        'metadata': 'The sample `metadata`.',
        'compare': 'The type of comparison that will be used to analyze the'
                   ' input `diversity_measure`.'
                   ' The "baseline" comparison defines Group A as the'
                   ' timepoint provided to `against_group` (sourced from'
                   ' `time_column`), and Group B as all other timepoints'
                   ' contained in `time_column`. The "reference" comparison'
                   ' defines Group A as the reference/control provided to'
                   ' `against_group` (sourced from either `reference_column`'
                   ' or `control_column`), and Group B as all timepoints'
                   ' contained in `time_column`. The "consecutive" comparison'
                   ' defines Group A as "timepoint n", and Group B as'
                   ' "timepoint n+1" (both sourced from `time_column`).'
                   ' The "all-pairwise" comparison defines Group A as all'
                   ' groups in `reference_column` and `control_column`, and'
                   ' Group B is all timepoints in `time_column`.',
        'time_column': 'The column within the `metadata` that the'
                       ' `diversity_measure` should be grouped by.'
                       ' This column should contain simple integer values.',
        'control_column': 'The column within the `metadata` that contains any'
                          ' relevant control group IDs.'
                          ' Actual treatment samples should not contain any'
                          ' value within this column.',
        'reference_column': 'The column within the `metadata` that contains'
                            ' the sample to use as a reference for a given'
                            ' beta `diversity_measure`. For example, this'
                            ' may be the relevant donor sample to compare'
                            ' against.',
        'subject_column': 'The column within the `metadata` that contains the'
                          ' subject ID to be tracked against timepoints.',
        'filter_missing_references': 'Filter out references contained within'
                                     ' the metadata that are not present in'
                                     ' the diversity measure.'
                                     ' Default behavior is to raise an error.',
        'where': 'Additional filtering for the associated `metadata` file.'
                 ' This can be used to filter by a subset of the `metadata`,'
                 ' such as a specific value in one of the `metadata` columns.',
        'against_group': 'Based on the selected comparison, this is the column'
                         ' that will be used to compare all samples against.',
        'alternative': 'The "two-sided" alternative hypothesis is that the'
                       ' median of Group A does not equal the median of Group'
                       ' B. The "greater" alternative hypothesis is that the'
                       ' median of group A is greater than the median of Group'
                       ' B. The "less" alternative hypothesis is that the'
                       ' median of group A is less than the median of Group'
                       ' B.',
        'p_val_approx': '"exact" will calculate an exact p-value'
                        ' for distributions, "asymptotic" will use a normal'
                        ' distribution, and "auto" will use either "exact"'
                        ' when one of the groups has less than 8 observations'
                        ' and there are no ties, otherwise "asymptotic".'
    },
    output_descriptions={
        'stats': 'Either the Mann-Whitney U or Wilcoxon SRT table'
                 ' for the chosen comparison.',
        'raincloud_plot': 'Raincloud plot for the computed significance test'
                          ' (either Mann-Whitney U or Wilxocon SRT) from the'
                          ' grouped diversity data and selected comparison.',
    },
    name='Engraftment Pipeline for FMT Analysis',
    description='Chloe will write this.',
    examples={
        'engraftment_baseline': ex.engraftment_baseline
    }
)

plugin.methods.register_function(
    function=group_timepoints,
    inputs={'diversity_measure': DistanceMatrix | SampleData[AlphaDiversity]},
    parameters={'metadata': Metadata, 'time_column': Str,
                'reference_column': Str, 'subject_column': T_subject,
                'control_column': Str, 'filter_missing_references': Bool,
                'where': Str},
    outputs=[('timepoint_dists', GroupDist[Ordered, T_dependence]),
             ('reference_dists', GroupDist[Unordered, Independent])],
    parameter_descriptions={
        'metadata': 'The sample metadata.',
        'time_column': 'The column within the `metadata` that the'
                       ' `diversity_measure` should be grouped by.'
                       ' This column should contain simple integer values.',
        'control_column': 'The column within the `metadata` that contains any'
                          ' relevant control group IDs.'
                          ' Actual treatment samples should not contain any'
                          ' value within this column.',
        'reference_column': 'The column within the `metadata` that contains'
                            ' the sample to use as a reference'
                            ' for a given beta `diversity_measure`.'
                            ' For example, this may be the relevant donor'
                            ' sample to compare against.',
        'subject_column': 'The column within the `metadata` that contains the'
                          ' subject ID to be tracked against timepoints.',
        'filter_missing_references': 'Filter out references contained within'
                                     ' the metadata that are not present'
                                     ' in the diversity measure.'
                                     ' Default behavior is to raise an error.',
        'where': 'Additional filtering for the associated `metadata` file.'
                 ' This can be used to filter by a subset of the `metadata`,'
                 ' such as a specific value in one of the `metadata` columns.',
    },
    output_descriptions={
        'timepoint_dists': 'The distributions for the `diversity_measure`,'
                           ' grouped by the selected `time_column`.'
                           ' May also contain subject IDs, if `subject_column`'
                           ' is provided in the `metadata`.',
        'reference_dists': 'The inter-group reference and inter-group control'
                           ' (when provided) distributions.'
                           ' When `diversity_measure` is DistanceMatrix, the'
                           ' inter-group calculations will be all pairwise'
                           ' comparisons within a group.'
                           ' Otherwise, these are just the per-sample'
                           ' measurements of alpha-diversity.'
    },
    name='',
    description='',
    examples={
        'group_timepoints_alpha_ind': ex.group_timepoints_alpha_independent,
        'group_timepoints_beta': ex.group_timepoints_beta
    }
)

plugin.methods.register_function(
    function=mann_whitney_u,
    inputs={'distribution': GroupDist[Unordered | Ordered, Independent],
            'against_each': GroupDist[Unordered | Ordered,
                                      Matched | Independent]},
    parameters={'compare': Str % Choices('reference', 'all-pairwise'),
                'reference_group': Str,
                'alternative': Str % Choices('two-sided', 'greater', 'less'),
                'p_val_approx': Str % Choices('auto', 'exact', 'asymptotic')},
    outputs=[('stats', StatsTable[Pairwise])],
    parameter_descriptions={
        'compare': 'The comparison that will be used to analyze the input'
                   ' `distribution`. Either "reference" or "all-pairwise"'
                   ' must be selected. The "reference" comparison defines'
                   ' Group A as the reference/control provided to'
                   ' `reference_group` (sourced from either `reference_column`'
                   ' or `control_column`), and Group B as all other groups.'
                   ' The "all-pairwise" comparison compares all groups to'
                   ' all other groups. If `against_each` is used, this will'
                   ' define Group B.',
        'reference_group': 'If "reference" is the selected comparison, this'
                           ' is the column that will be used to compare all'
                           ' other groups against.',
        'alternative': 'The "two-sided" alternative hypothesis is that the'
                       ' median of Group A does not equal the median of Group'
                       ' B. The "greater" alternative hypothesis is that the'
                       ' median of group A is greater than the median of Group'
                       ' B. The "less" alternative hypothesis is that the'
                       ' median of group A is less than the median of Group'
                       ' B.',
        'p_val_approx': '"exact" will calculate an exact p-value for'
                        ' distributions, "asymptotic" will use a normal'
                        ' distribution, and "auto" will use either "exact"'
                        ' when one of the groups has less than 8 observations'
                        ' and there are no ties, otherwise "asymptotic".'
    },
    output_descriptions={
        'stats': 'The Mann-Whitney U table for either the "reference"'
                 ' or "all-pairwise" comparison.',
    },
    name='Mann-Whitney U Test',
    description='',
    examples={
        'mann_whitney_pairwise': ex.mann_whitney_pairwise
    }
)

plugin.methods.register_function(
    function=wilcoxon_srt,
    inputs={'distribution': GroupDist[Ordered, Matched]},
    parameters={'compare': Str % Choices('baseline', 'consecutive'),
                'baseline_group': Str,
                'alternative': Str % Choices('two-sided', 'greater', 'less'),
                'p_val_approx': Str % Choices('auto', 'exact', 'asymptotic')},
    outputs=[('stats', StatsTable[Pairwise])],
    parameter_descriptions={
        'compare': 'The type of comparison that will be used to analyze the'
                   ' input `distribution`.'
                   ' The "baseline" comparison defines Group A as the'
                   ' timepoint provided to `baseline_group` (sourced from'
                   ' `time_column`), and Group B as all other timepoints'
                   ' contained in `time_column`. The "consecutive" comparison'
                   ' defines Group A as "timepoint n", and Group B as'
                   ' "timepoint n+1" (both sourced from `time_column`).',
        'baseline_group': 'If "baseline" is the selected comparison, this is'
                          ' the column that will be used to compare all'
                          ' other groups against.',
        'alternative': 'The "two-sided" alternative hypothesis is that the'
                       ' median of Group A does not equal the median of Group'
                       ' B. The "greater" alternative hypothesis is that the'
                       ' median of group A is greater than the median of Group'
                       ' B. The "less" alternative hypothesis is that the'
                       ' median of group A is less than the median of Group'
                       ' B.',
        'p_val_approx': '"exact" will calculate an exact p-value for'
                        ' distributions of up to 25 (inclusive) measurements,'
                        ' "asymptotic" will use a normal distribution,'
                        ' and "auto" will use either "exact" or "approx"'
                        ' depending on size.'
    },
    output_descriptions={
        'stats': 'The Wilcoxon SRT table for either the "baseline"'
                 ' or "consecutive" comparison.',
    },
    name='Wilcoxon Signed Rank Test',
    description='',
    examples={
        'wilcoxon_baseline0': ex.wilcoxon_baseline0
    }
)

plugin.visualizers.register_function(
    function=plot_rainclouds,
    inputs={
        'data': GroupDist[Ordered, Matched],
        'stats': StatsTable[Pairwise]
    },
    parameters={},
    input_descriptions={
        'data': 'The group distributions to plot.',
        'stats': 'Statistical tests to display.'
    },
    parameter_descriptions={},
    name='Raincloud plots',
    description='Plot raincloud distributions for each group.'
)

importlib.import_module('q2_fmt._transformer')
importlib.import_module('q2_fmt._validator')
