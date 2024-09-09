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
                           Int, Range, List)
from q2_types.sample_data import SampleData, AlphaDiversity
from q2_types.distance_matrix import DistanceMatrix

import q2_fmt
import q2_fmt._peds
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

# Description Variables
metadata = 'The sample `metadata`.'
time_column = ('The column within `metadata` that the `table` should be'
               ' grouped by. This column should contain integer values.')
reference_column = ('The column in `metadata` indicating the sample to use as'
                    ' the reference. For example, this may list the relevant'
                    ' donor sample for each recipient sample.')
subject_column = ('The column within `metadata` that contains the subject'
                  ' ID for each sample.')
filter_missing_references = ('Filter out references contained within the'
                             ' metadata that are not present in the input'
                             ' feature table. Default behavior is to raise an'
                             ' error if some reference ids are in the metadata'
                             ' but not the feature table.')
drop_incomplete_subjects = ('Filter out subjects that do not have a sample at'
                            ' every timepoint. Default behavior is to raise an'
                            ' error if any subject is missing a timepoint.')
drop_incomplete_timepoints = ('Filter out multiple specified timepoints.'
                              ' This is useful for removing frequently missing'
                              ' timepoints which cause many subjects to be'
                              ' dropped. Default behavior is to raise an error'
                              ' if any subject is missing a timepoint.')
level_delimiter = 'delimiter to split taxonomic label on'
control_column = ('The column within `metadata` that contains any relevant'
                  ' control group IDs. Actual treatment samples should not'
                  ' contain any value within this column.')
distance_to = 'What reference type to calculate distance to against'
baseline_timepoint = ('If `baseline` is selected for `distance_to`, the'
                      ' timepoint to use as baseline reference.')
where = ('Additional filtering for the associated `metadata` file. This can be'
         ' used to filter by a subset of `metadata`, such as a specific value'
         ' in one of `metadata` columns.')

per_subject_stats = ('Table describing significance of PEDS scores compared to'
                     ' mismatched donor-recipient pairs on a per-subject'
                     ' basis.')
global_stats = ('Table describing significance of PEDS scores across all'
                ' subjects.')
peds_table = 'The feature table to calculate PEDS on.'
peds_dists = ('The distributions for the PEDS measure, grouped by the selected'
              ' `time_column`. Also contains the numerator and denominator for'
              ' PEDS calulations. May also contain subject IDs, if'
              ' `subject_column` is  provided in `metadata`.')


plugin.pipelines.register_function(
    function=q2_fmt.cc,
    inputs={'diversity_measure': DistanceMatrix | SampleData[AlphaDiversity]},
    parameters={'metadata': Metadata,
                'compare': T_compare,
                'distance_to': Str % Choices('baseline', 'donor'),
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
        'metadata': metadata,
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
        'distance_to': distance_to,
        'time_column': time_column,
        'control_column': control_column,
        'reference_column': reference_column,
        'subject_column': subject_column,
        'filter_missing_references': filter_missing_references,
        'baseline_timepoint': baseline_timepoint,
        'where': where,
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
                          ' (either Mann-Whitney U or Wilcoxon SRT) from the'
                          ' grouped diversity data and selected comparison.',
    },
    name='Pipeline for FMT Diversity Community'
         ' Coalescence Analysis',
    description='Applies group timepoints to Fecal Microbiota Transplant data'
         ' and then applies statistical tests (Wilcoxon signed rank test or'
         ' Mann–Whitney U test) on alpha or beta diversity metrics to assess'
         ' Community Coalescense (CC) of the FMT.',
    examples={
        'cc_baseline': ex.cc_baseline
    }
)

plugin.methods.register_function(
    function=q2_fmt.group_timepoints,
    inputs={'diversity_measure': DistanceMatrix | SampleData[AlphaDiversity]},
    parameters={'metadata': Metadata,
                'distance_to': Str % Choices('baseline', 'donor'),
                'time_column': Str,
                'reference_column': Str, 'subject_column': T_subject,
                'group_column': T_group, 'control_column': Str,
                'filter_missing_references': Bool, 'where': Str,
                'baseline_timepoint': Str, 'where': Str},
    outputs=[('timepoint_dists', Dist1D[T_nested, T_dependence]),
             ('reference_dists', Dist1D[Unordered, Independent])],
    input_descriptions={
        'diversity_measure': 'diversity metric to put in a long form data'
                             ' structure'
    },
    parameter_descriptions={
        'metadata': metadata,
        'distance_to': distance_to,
        'time_column': time_column,
        'control_column': control_column,
        'reference_column': reference_column,
        'subject_column': subject_column,
        'group_column': 'The column within the metadata that contains'
                        ' information about groups (ex: treatment group)'
                        ' in order to compare engraftment between groups',
        'filter_missing_references': filter_missing_references,
        'where': where,
        'baseline_timepoint': baseline_timepoint,

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
                'drop_incomplete_timepoints': List[Str],
                'level_delimiter': Str},
    outputs=[('peds_heatmap', Visualization)],
    input_descriptions={'table': peds_table},
    parameter_descriptions={
        'metadata': metadata,
        'peds_metric': 'PEDS metric to run.',
        'time_column': time_column,
        'reference_column': reference_column,
        'subject_column': subject_column,
        'filter_missing_references': filter_missing_references,
        'drop_incomplete_subjects': drop_incomplete_subjects,
        'drop_incomplete_timepoints': drop_incomplete_timepoints,
        'level_delimiter': level_delimiter},
    output_descriptions={'peds_heatmap': 'PEDS heatmap visualization'},
    name='PEDS pipeline to calculate feature or sample PEDS',
    description='Runs a pipeline to calculate sample or feature PEDS,'
                '  and generate the relevant heatmap'
)

plugin.visualizers.register_function(
    function=q2_fmt.peds_heatmap,
    inputs={'data': Dist1D[Ordered, Matched] % Properties("peds") |
            Dist1D[Ordered, Matched] % Properties("pprs"),
            'per_subject_stats': StatsTable[Pairwise],
            'global_stats': StatsTable[Pairwise]},
    input_descriptions={'data': 'PEDS output to plot',
                        'per_subject_stats': per_subject_stats,
                        'global_stats': global_stats},
    parameters={'level_delimiter': Str},
    parameter_descriptions={
                            'level_delimiter': level_delimiter},
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
                'drop_incomplete_timepoints': List[Str]},
    outputs=[('peds_dists', Dist1D[Ordered, Matched] % Properties("peds"))],
    input_descriptions={'table': peds_table},
    parameter_descriptions={
        'metadata': metadata,
        'time_column': time_column,
        'reference_column': reference_column,
        'subject_column': subject_column,
        'filter_missing_references': filter_missing_references,
        'drop_incomplete_subjects': drop_incomplete_subjects,
        'drop_incomplete_timepoints': drop_incomplete_timepoints
    },
    output_descriptions={
        'peds_dists': peds_dists
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
    input_descriptions={'table': peds_table},
    parameter_descriptions={
        'metadata': metadata,
        'time_column': time_column,
        'reference_column': reference_column,
        'subject_column': subject_column,
        'filter_missing_references': filter_missing_references,
    },
    output_descriptions={
        'peds_dists': peds_dists
    },
    name='Porportional Engraftment of Donor Strains per feature',
    description='Calculates how many recipients recieved a given'
                ' donated material feature ',
    examples={
        'peds_methods': ex.feature_peds_method
    }
)

plugin.methods.register_function(
    function=q2_fmt._peds.sample_pprs,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency |
                                  PresenceAbsence]},
    parameters={'metadata': Metadata, 'time_column': Str,
                'baseline_timepoint': Str, 'subject_column': Str,
                'filter_missing_references': Bool,
                'drop_incomplete_subjects': Bool,
                'drop_incomplete_timepoint': List[Str]},
    outputs=[('peds_dists', Dist1D[Ordered, Matched] % Properties("peds"))],
    input_descriptions={
        'table': peds_table},
    parameter_descriptions={
        'metadata': metadata,
        'time_column': time_column,
        'baseline_timepoint': baseline_timepoint,
        'subject_column': subject_column,
        'filter_missing_references': filter_missing_references,
        'drop_incomplete_subjects': drop_incomplete_subjects,
        'drop_incomplete_timepoints': drop_incomplete_timepoints
    },
    output_descriptions={
        'peds_dists': peds_dists
    },
    name='',
    description='',
    citations=[citations['aggarwala_precise_2021']],
    examples={
        'peds_methods': ex.peds_method}
)

plugin.methods.register_function(
    function=q2_fmt.peds_simulation,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency |
                                  PresenceAbsence]},
    parameters={'metadata': Metadata,
                'time_column': Str,
                'reference_column': Str,
                'subject_column': T_subject,
                'filter_missing_references': Bool,
                'drop_incomplete_subjects': Bool,
                'drop_incomplete_timepoints': List[Str],
                'num_iterations': Int % Range(99, None)},
    outputs=[('per_subject_stats', StatsTable[Pairwise]),
             ('global_stats', StatsTable[Pairwise])],
    parameter_descriptions={
        'metadata': metadata,
        'time_column': time_column,
        'reference_column': reference_column,
        'subject_column': subject_column,
        'filter_missing_references': filter_missing_references,
        'drop_incomplete_subjects': drop_incomplete_subjects,
        'drop_incomplete_timepoints': drop_incomplete_timepoints,
        'num_iterations': 'The number of iterations to run the Monte Carlo'
                          ' simulation on'
    },
    output_descriptions={
        'per_subject_stats': per_subject_stats,
        'global_stats': global_stats
    },
    name='PEDS Monte Carlo simulation',
    description='A Monte Carlo simulation that randomizes the relationships'
                ' between donors and recipients, to test whether the PEDS'
                ' score between a recipient and their actual donor is'
                ' significantly higher than PEDS scores between other'
                ' recipients paired with random donors. This is intended to'
                ' only work in studies where there are distinct donors,'
                ' and will yield insignificant results if there are too'
                ' few donors. This method attempts to quantify the extent to'
                ' which indicator features that are unique to a given donor'
                ' are transferred to their recipients, as opposed to features'
                ' that are not indicative of any specific donor. Note: '
                ' PEDS Monte Carlo simulation may have dependency issues'
                ' between samples and the simulated background'
                ' distribution that can make the test'
                ' overly conservative. This can be fixed by filtering down to'
                ' a single timepoint before running this method. Additionally'
                ' if there are many baseline timepoints, the global test may'
                ' be too conservative and this can be addressed by filtering'
                ' out baseline samples prior to running this method.',
    citations=[citations['aggarwala_precise_2021'],
               citations['stouffer_1949_american'],
               citations['Benjamini_fdr_1995']],
    examples={
        'peds_methods': ex.simulation_peds_method
    }
)

importlib.import_module('q2_fmt._transformer')
