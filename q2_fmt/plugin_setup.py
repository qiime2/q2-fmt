# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

from qiime2.plugin import (Str, Plugin, Metadata, TypeMap,
                           Bool, Choices, Visualization, Properties, Citations,
                           List)
from q2_types.sample_data import SampleData, AlphaDiversity
from q2_types.distance_matrix import DistanceMatrix

import q2_fmt
from q2_types.feature_table import (
    FeatureTable, Frequency, RelativeFrequency, PresenceAbsence)
from q2_stats._type import (Dist1D, Matched, Independent, Ordered,
                            Unordered, StatsTable, Pairwise, NestedOrdered)
import q2_fmt._examples as ex

citations = Citations.load('citations.bib', package='q2_fmt')

plugin = Plugin(name='fmt',
                version=q2_fmt.__version__,
                website='https://github.com/qiime2/q2-fmt',
                package='q2_fmt',
                description='This QIIME 2 plugin supports FMT analyses.',
                short_description='Plugin for analyzing FMT data.')

T_group, T_nested = TypeMap({
    Bool % Choices(False): Ordered,
    Str: NestedOrdered
})

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
                'distance_to': Str,
                'time_column': Str, 'reference_column': Str,
                'subject_column': T_engraft_subject, 'control_column': Str,
                'filter_missing_references': Bool, 'baseline_timepoint': Str,
                'where': Str, 'against_group': Str,
                'alternative': Str % Choices('two-sided', 'greater', 'less'),
                'p_val_approx': Str % Choices('auto', 'exact', 'asymptotic')},
    outputs=[
        ('stats', StatsTable[Pairwise]),
        ('raincloud_plot', Visualization)
    ],
    input_descriptions={
        'diversity_measure': 'A diveristy measure to compare against donated'
                             ' microbiome',
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
        'distance_to': 'what reference type to calculate distance to against',
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
        'baseline_timepoint': 'If `baseline` is selected for `distance_to`,'
                              ' the timepoint to use as baseline reference.',
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
    name='Engraftment Pipeline for FMT Diversity Community'
         ' Coalescence Analysis',
    description='Applies group timepoints to Fecal Microbiota Transplant data'
         ' and then applies statistical tests (Wilcoxon signed rank test or'
         ' Mann–Whitney U test) on alpha or beta diversity metrics to assess'
         ' Community Coalescense of the FMT.',
    examples={
        'engraftment_baseline': ex.engraftment_baseline
    }
)

plugin.methods.register_function(
    function=q2_fmt.group_timepoints,
    inputs={'diversity_measure': DistanceMatrix | SampleData[AlphaDiversity]},
    parameters={'metadata': Metadata, 'distance_to': Str, 'time_column': Str,
                'reference_column': Str, 'subject_column': T_subject,
                'group_column': T_group, 'control_column': Str,
                'filter_missing_references': Bool, 'where': Str,
                'baseline_timepoint': Str, 'where': Str},
    outputs=[('timepoint_dists', Dist1D[Ordered, T_dependence]),
             ('reference_dists', Dist1D[Unordered, Independent])],
    input_descriptions={
        'diversity_measure': 'diversity metric to put in a long form data'
                             ' structure'
    },
    parameter_descriptions={
        'metadata': 'The sample metadata.',
        'distance_to': 'what reference type to calculate distance to against',
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
        'group_column': 'The column within the metadata that contains'
                        ' information about groups (ex: treatment group)'
                        ' in order to compare engraftment between groups',
        'filter_missing_references': 'Filter out references contained within'
                                     ' the metadata that are not present'
                                     ' in the diversity measure.'
                                     ' Default behavior is to raise an error.',
        'where': 'Additional filtering for the associated `metadata` file.'
                 ' This can be used to filter by a subset of the `metadata`,'
                 ' such as a specific value in one of the `metadata` columns.',
        'baseline_timepoint': 'If `baseline` is selected for `distance_to`,'
                              ' the timepoint to use as baseline reference.',
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
    name='Prep Method for organizing metadata and diversity metrics'
         ' for appropriate statisics ',
    description='Puts diversity metric data in a long form data structure'
                ' that includes relevant metadata information',
    examples={
        'group_timepoints_alpha_ind': ex.group_timepoints_alpha_independent,
        'group_timepoints_beta': ex.group_timepoints_beta
    }
)

plugin.pipelines.register_function(
    function=q2_fmt.peds,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency |
                                  PresenceAbsence]},
    parameters={'metadata': Metadata,
                'peds_metric': Str % Choices('feature', 'sample'),
                'time_column': Str, 'reference_column': Str,
                'subject_column': Str,
                'filter_missing_references': Bool,
                'drop_incomplete_subjects': Bool,
                'drop_incomplete_timepoint': List[Str],
                'level_delimiter': Str},
    outputs=[('peds_heatmap', Visualization)],
    input_descriptions={'table': 'A feature table to calculate PEDS on'},
    parameter_descriptions={
        'metadata': 'The sample metadata.',
        'peds_metric': 'PEDS metric to run.',
        'time_column': 'The column within the `metadata` that the'
                       ' `table` should be grouped by. This column'
                       ' should contain simple integer values.',
        'reference_column': 'The column within the `metadata` that contains'
                            ' the sample to use as a reference'
                            ' for a given `table`.'
                            ' For example, this may be the relevant donor'
                            ' sample to compare against.',
        'subject_column': 'The column within the `metadata` that contains the'
                          ' subject ID to be tracked against timepoints.',
        'filter_missing_references': 'Filter out references contained within'
                                     ' the metadata that are not present'
                                     ' in the table.'
                                     ' Default behavior is to raise an error.',
        'drop_incomplete_subjects': 'Filter out subjects that do not have'
                                    ' a sample at every timepoint.'
                                    ' Default behavior is to raise an error.',
        'drop_incomplete_timepoint': 'Filter out a list of provided timepoint.'
                                     ' This will be performed before'
                                     ' drop_incomplete_subjects if the'
                                     ' drop_incomplete_subjects parameter is'
                                     ' passed.',
        'level_delimiter': 'delimiter to split taxonomic string on'},
    output_descriptions={'peds_heatmap': 'PEDS heatmap visualization'},
    name='PEDS pipeline to calculate feature or sample PEDS',
    description='Runs a pipeline to calculate sample or feature PEDS,'
                '  and generate the relevant heatmap'
)

plugin.visualizers.register_function(
    function=q2_fmt.peds_heatmap,
    inputs={'data': Dist1D[Ordered, Matched] % Properties("peds")},
    input_descriptions={'data': 'PEDS output to plot'},
    parameters={'level_delimiter': Str},
    parameter_descriptions={
                            'level_delimiter': 'delimiter to split taxonomic'
                                               ' string on'},
    name='PEDS Heatmap',
    description='Plot heatmap for PEDS value over time'
)

plugin.methods.register_function(
    function=q2_fmt.sample_peds,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency |
                                  PresenceAbsence]},
    parameters={'metadata': Metadata, 'time_column': Str,
                'reference_column': Str, 'subject_column': Str,
                'filter_missing_references': Bool,
                'drop_incomplete_subjects': Bool,
                'drop_incomplete_timepoint': List[Str]},
    outputs=[('peds_dists', Dist1D[Ordered, Matched] % Properties("peds"))],
    input_descriptions={'table': 'A feature table to calculate PEDS on'},
    parameter_descriptions={
        'metadata': 'The sample metadata.',
        'time_column': 'The column within the `metadata` that the'
                       ' `table` should be grouped by. This column'
                       ' should contain simple integer values.',
        'reference_column': 'The column within the `metadata` that contains'
                            ' the sample to use as a reference'
                            ' for a given `table`.'
                            ' For example, this may be the relevant donor'
                            ' sample to compare against.',
        'subject_column': 'The column within the `metadata` that contains the'
                          ' subject ID to be tracked against timepoints.',
        'filter_missing_references': 'Filter out references contained within'
                                     ' the metadata that are not present'
                                     ' in the table.'
                                     ' Default behavior is to raise an error.',
        'drop_incomplete_subjects': 'Filter out subjects that do not have'
                                    ' a sample at every timepoint.'
                                    ' Default behavior is to raise an error.',
        'drop_incomplete_timepoint': 'Filter out a list of provided timepoint.'
                                     ' This will be performed before'
                                     ' drop_incomplete_subjects if the'
                                     ' drop_incomplete_subjects parameter is'
                                     ' passed.'
    },
    output_descriptions={
        'peds_dists': 'The distributions for the PEDS measure,'
                      ' grouped by the selected `time_column`.'
                      ' also contains the Numerator and Denominator for'
                      ' PEDS calulations. May also contain subject IDs,'
                      ' if `subject_column` is provided in the `metadata`.'
    },
    name='Proportional Engraftment of Donor Strains (Features) in each'
         ' recipient sample',
    description='Calculates percentage of microbes that where found in the '
    'donated material that are found in the recipient.',
    citations=[citations['aggarwala_precise_2021']],
    examples={
        'peds_methods': ex.peds_method
    }
)
plugin.methods.register_function(
    function=q2_fmt.feature_peds,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency |
                                  PresenceAbsence]},
    parameters={'metadata': Metadata, 'time_column': Str,
                'reference_column': Str, 'subject_column': Str,
                'filter_missing_references': Bool},
    outputs=[('peds_dists', Dist1D[Ordered, Matched] % Properties("peds"))],
    input_descriptions={'table': 'A feature table to calculate PEDS on.'},
    parameter_descriptions={
        'metadata': 'The sample metadata.',
        'time_column': 'The column within the `metadata` that the'
                       ' `table` should be grouped by. This column'
                       ' should contain simple integer values.',
        'reference_column': 'The column within the `metadata` that contains'
                            ' the sample to use as a reference'
                            ' for a given `table`.'
                            ' For example, this may be the relevant donor'
                            ' sample to compare against.',
        'subject_column': 'The column within the `metadata` that contains the'
                          ' subject ID to be tracked against timepoints.',
        'filter_missing_references': 'Filter out references contained within'
                                     ' the metadata that are not present'
                                     ' in the table.'
                                     ' Default behavior is to raise an error.',
    },
    output_descriptions={
        'peds_dists': 'The distributions for the PEDS measure,'
                      ' grouped by the selected `time_column`.'
                      ' also contains the Numerator and Denominator for'
                      ' PEDS calulations. May also contain subject IDs,'
                      ' if `subject_column` is provided in the `metadata`.'
    },
    name='Porportional Engraftment of Donor Strains per feature',
    description='Calculates how many recipients recieved a given'
                ' donated material feature ',
    examples={
        'peds_methods': ex.feature_peds_method
    }
)

importlib.import_module('q2_fmt._transformer')
