# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from skbio.stats.distance import DistanceMatrix

from qiime2.plugin.testing import TestPluginBase
from qiime2 import Metadata

from q2_fmt._engraftment import group_timepoints
from q2_fmt._stats import wilcoxon_srt, mann_whitney_u
from q2_fmt._examples import faithpd_timedist_factory, faithpd_refdist_factory


class TestBase(TestPluginBase):
    package='q2_fmt.tests'

    def setUp(self):
        super().setUp()

        self.md_beta = Metadata.load(self.get_data_path('sample_metadata_donors.tsv'))
        self.md_alpha = Metadata.load(self.get_data_path('sample_metadata_alpha_div.tsv'))

        self.dm = DistanceMatrix.read(self.get_data_path('dist_matrix_donors.tsv')).to_series()
        self.alpha = pd.read_csv(self.get_data_path('alpha_div.tsv'), sep='\t', index_col=0, squeeze=True)

        self.faithpd_timedist = faithpd_timedist_factory().view(pd.DataFrame)
        self.faithpd_refdist = faithpd_refdist_factory().view(pd.DataFrame)


class ErrorMixins:
    def test_with_time_column_input_not_in_metadata(self):
        with self.assertRaisesRegex(ValueError, 'time_column.*foo.*metadata'):
            group_timepoints(diversity_measure=self.div,
                             metadata=self.md,
                             time_column='foo',
                             reference_column='relevant_donor',
                             control_column='control')

    def test_with_reference_column_input_not_in_metadata(self):
        with self.assertRaisesRegex(ValueError, 'reference_column.*foo.*metadata'):
            group_timepoints(diversity_measure=self.div,
                             metadata=self.md,
                             time_column='days_post_transplant',
                             reference_column='foo',
                             control_column='control')

    def test_with_control_column_input_not_in_metadata(self):
        with self.assertRaisesRegex(ValueError, 'control_column.*foo.*metadata'):
            group_timepoints(diversity_measure=self.div,
                             metadata=self.md,
                             time_column='days_post_transplant',
                             reference_column='relevant_donor',
                             control_column='foo')

    def test_with_non_numeric_time_column(self):
        with self.assertRaisesRegex(ValueError, 'time_column.*categorical.*numeric'):
            group_timepoints(diversity_measure=self.div,
                             metadata=self.md,
                             time_column='non_numeric_time_column',
                             reference_column='relevant_donor',
                             control_column='control')

class TestAlphaErrors(TestBase, ErrorMixins):
    def setUp(self):
        super().setUp()

        self.div = self.alpha
        self.md = self.md_alpha

class TestBetaErrors(TestBase, ErrorMixins):
    def setUp(self):
        super().setUp()

        self.div = self.dm
        self.md = self.md_beta

class TestGroupTimepoints(TestBase):
    # Beta Diversity (Distance Matrix) Test Cases
    def test_beta_dists_with_donors_and_controls(self):
        exp_time_df = pd.DataFrame({
            'id': ['sampleA', 'sampleB', 'sampleC', 'sampleD', 'sampleE'],
            'measure': [0.45, 0.40, 0.28, 0.78, 0.66],
            'group': [7.0, 7.0, 9.0, 11.0, 11.0]
        })

        exp_ref_df = pd.DataFrame({
            'id': ['donor1..donor2', 'donor1..donor3', 'donor2..donor3',
                   'sampleB..sampleC', 'sampleB..sampleD', 'sampleC..sampleD'],
            'measure': [0.24, 0.41, 0.74, 0.37, 0.44, 0.31],
            'group': ['reference', 'reference', 'reference', 'control1', 'control1', 'control1'],
            'A': ['donor1', 'donor1', 'donor2', 'sampleB', 'sampleB', 'sampleC'],
            'B': ['donor2', 'donor3', 'donor3', 'sampleC', 'sampleD', 'sampleD']
        })

        time_df, ref_df = group_timepoints(diversity_measure=self.dm,
                                           metadata=self.md_beta,
                                           time_column='days_post_transplant',
                                           reference_column='relevant_donor',
                                           control_column='control')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_beta_dists_with_donors_controls_and_subjects(self):
        exp_time_df = pd.DataFrame({
            'id': ['sampleA', 'sampleB', 'sampleC', 'sampleD', 'sampleE'],
            'measure': [0.45, 0.40, 0.28, 0.78, 0.66],
            'group': [7.0, 7.0, 9.0, 11.0, 11.0],
            'subject': ['subject1', 'subject1', 'subject1', 'subject2', 'subject2']
        })

        exp_ref_df = pd.DataFrame({
            'id': ['donor1..donor2', 'donor1..donor3', 'donor2..donor3',
                   'sampleB..sampleC', 'sampleB..sampleD', 'sampleC..sampleD'],
            'measure': [0.24, 0.41, 0.74, 0.37, 0.44, 0.31],
            'group': ['reference', 'reference', 'reference', 'control1', 'control1', 'control1'],
            'A': ['donor1', 'donor1', 'donor2', 'sampleB', 'sampleB', 'sampleC'],
            'B': ['donor2', 'donor3', 'donor3', 'sampleC', 'sampleD', 'sampleD']
        })

        time_df, ref_df = group_timepoints(diversity_measure=self.dm,
                                           metadata=self.md_beta,
                                           time_column='days_post_transplant',
                                           reference_column='relevant_donor',
                                           control_column='control',
                                           subject_column='subject')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_beta_dists_with_same_donor_for_all_samples(self):
        with self.assertRaisesRegex(TypeError, 'Single reference value detected'):
            group_timepoints(diversity_measure=self.dm,
                             metadata=self.md_beta,
                             time_column='days_post_transplant',
                             reference_column='relevant_donor_all',
                             control_column='control')

    def test_beta_dists_with_one_donor_and_controls(self):
        with self.assertRaisesRegex(KeyError, 'Missing references for the associated sample data'):
            group_timepoints(diversity_measure=self.dm,
                           metadata=self.md_beta,
                           time_column='days_post_transplant',
                           reference_column='single_donor',
                           control_column='control')

    def test_beta_dists_with_donors_and_one_control(self):
        with self.assertRaisesRegex(ValueError, 'One or less controls detected'):
            group_timepoints(diversity_measure=self.dm,
                           metadata=self.md_beta,
                           time_column='days_post_transplant',
                           reference_column='relevant_donor',
                           control_column='single_control')

    def test_beta_dists_with_donors_no_controls(self):
        exp_time_df = pd.DataFrame({
            'id': ['sampleA', 'sampleB', 'sampleC', 'sampleD', 'sampleE'],
            'measure': [0.45, 0.40, 0.28, 0.78, 0.66],
            'group': [7.0, 7.0, 9.0, 11.0, 11.0]
        })

        exp_ref_df = pd.DataFrame({
            'id': ['donor1..donor2', 'donor1..donor3', 'donor2..donor3'],
            'measure': [0.24, 0.41, 0.74],
            'group': ['reference', 'reference', 'reference'],
            'A': ['donor1', 'donor1', 'donor2'],
            'B': ['donor2', 'donor3', 'donor3']
        })

        time_df, ref_df = group_timepoints(diversity_measure=self.dm,
                                           metadata=self.md_beta,
                                           time_column='days_post_transplant',
                                           reference_column='relevant_donor')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_beta_dists_no_donors_with_controls(self):
        with self.assertRaisesRegex(TypeError, r"group_timepoints\(\) missing 1 required positional argument: "
                                    "'reference_column'"):
            group_timepoints(diversity_measure=self.dm,
                             metadata=self.md_beta,
                             time_column='days_post_transplant',
                             control_column='control')

    def test_beta_dists_with_invalid_ref_column(self):
        with self.assertRaisesRegex(KeyError, 'References included in the metadata are missing'
                                    ' from the diversity measure.*foo.*bar.*baz'):
            group_timepoints(diversity_measure=self.dm,
                             metadata=self.md_beta,
                             time_column='days_post_transplant',
                             reference_column='invalid_ref_control',
                             control_column='control')

    def test_beta_dists_with_empty_diversity_series(self):
        empty_beta_series = pd.Series()

        with self.assertRaisesRegex(ValueError, 'Empty diversity measure detected'):
            group_timepoints(diversity_measure=empty_beta_series,
                             metadata=self.md_beta,
                             time_column='days_post_transplant',
                             reference_column='relevant_donor',
                             control_column='control')

    def test_beta_dists_with_extra_samples_in_metadata_not_in_diversity(self):
        extra_md = Metadata.load(self.get_data_path('sample_metadata_donors_missing.tsv'))

        exp_time_df = pd.DataFrame({
            'id': ['sampleA', 'sampleB', 'sampleC', 'sampleD', 'sampleE'],
            'measure': [0.45, 0.40, 0.28, 0.78, 0.66],
            'group': [7.0, 7.0, 9.0, 11.0, 11.0]
        })

        exp_ref_df = pd.DataFrame({
            'id': ['donor1..donor2', 'donor1..donor3', 'donor2..donor3',
                   'sampleB..sampleC', 'sampleB..sampleD', 'sampleC..sampleD'],
            'measure': [0.24, 0.41, 0.74, 0.37, 0.44, 0.31],
            'group': ['reference', 'reference', 'reference', 'control1', 'control1', 'control1'],
            'A': ['donor1', 'donor1', 'donor2', 'sampleB', 'sampleB', 'sampleC'],
            'B': ['donor2', 'donor3', 'donor3', 'sampleC', 'sampleD', 'sampleD']
        })

        time_df, ref_df = group_timepoints(diversity_measure=self.dm,
                                           metadata=extra_md,
                                           time_column='days_post_transplant',
                                           reference_column='relevant_donor',
                                           control_column='control')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_beta_dists_with_extra_samples_in_diversity_not_in_metadata(self):
        extra_dm = DistanceMatrix.read(self.get_data_path('dist_matrix_donors_missing.tsv')).to_series()

        with self.assertRaisesRegex(ValueError, 'The following IDs are not present in the metadata'):
            group_timepoints(diversity_measure=extra_dm,
                             metadata=self.md_beta,
                             time_column='days_post_transplant',
                             reference_column='relevant_donor',
                             control_column='control')

    # Alpha Diversity (Series) Test Cases
    def test_alpha_dists_with_donors_controls(self):
        exp_time_df = pd.DataFrame({
            'id': ['sampleA', 'sampleB', 'sampleC', 'sampleD',
                   'sampleE', 'sampleF', 'sampleG'],
            'measure': [24, 37, 15, 6, 44, 17, 29],
            'group': [7.0, 7.0, 9.0, 11.0, 11.0, 9.0, 7.0],
        })

        exp_ref_df = pd.DataFrame({
            'id': ['donor1', 'donor2', 'donor3', 'donor4',
                   'sampleC', 'sampleD', 'sampleE', 'sampleF'],
            'measure': [32, 51, 3, 19, 15, 6, 44, 17],
            'group': ['reference', 'reference', 'reference', 'reference',
                      'control1', 'control1', 'control2', 'control2']
        })

        time_df, ref_df = group_timepoints(diversity_measure=self.alpha,
                                           metadata=self.md_alpha,
                                           time_column='days_post_transplant',
                                           reference_column='relevant_donor',
                                           control_column='control')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_alpha_dists_with_donors_controls_and_subjects(self):
        exp_time_df = pd.DataFrame({
            'id': ['sampleA', 'sampleB', 'sampleC', 'sampleD',
                   'sampleE', 'sampleF', 'sampleG'],
            'measure': [24, 37, 15, 6, 44, 17, 29],
            'group': [7.0, 7.0, 9.0, 11.0, 11.0, 9.0, 7.0],
            'subject': ['subject1', 'subject1', 'subject2',
                        'subject1', 'subject2', 'subject2', 'subject1']
        })

        exp_ref_df = pd.DataFrame({
            'id': ['donor1', 'donor2', 'donor3', 'donor4',
                   'sampleC', 'sampleD', 'sampleE', 'sampleF'],
            'measure': [32, 51, 3, 19, 15, 6, 44, 17],
            'group': ['reference', 'reference', 'reference', 'reference',
                      'control1', 'control1', 'control2', 'control2']
        })

        time_df, ref_df = group_timepoints(diversity_measure=self.alpha,
                                           metadata=self.md_alpha,
                                           time_column='days_post_transplant',
                                           reference_column='relevant_donor',
                                           control_column='control',
                                           subject_column='subject')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_alpha_dists_with_same_donor_for_all_samples(self):
        exp_time_df = pd.DataFrame({
            'id': ['sampleA', 'sampleB', 'sampleC', 'sampleD',
                   'sampleE', 'sampleF', 'sampleG'],
            'measure': [24, 37, 15, 6, 44, 17, 29],
            'group': [7.0, 7.0, 9.0, 11.0, 11.0, 9.0, 7.0],
        })

        exp_ref_df = pd.DataFrame({
            'id': ['donor1', 'sampleC', 'sampleD', 'sampleE', 'sampleF'],
            'measure': [32, 15, 6, 44, 17],
            'group': ['reference', 'control1', 'control1', 'control2', 'control2']
        })

        time_df, ref_df = group_timepoints(diversity_measure=self.alpha,
                                           metadata=self.md_alpha,
                                           time_column='days_post_transplant',
                                           reference_column='relevant_donor_all',
                                           control_column='control')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_alpha_dists_with_one_donor_and_controls(self):
        with self.assertRaisesRegex(KeyError, 'Missing references for the associated sample data'):
            group_timepoints(diversity_measure=self.alpha,
                             metadata=self.md_alpha,
                             time_column='days_post_transplant',
                             reference_column='single_donor',
                             control_column='control')

    def test_alpha_dists_with_donors_and_one_control(self):
        exp_time_df = pd.DataFrame({
            'id': ['sampleA', 'sampleB', 'sampleC', 'sampleD',
                   'sampleE', 'sampleF', 'sampleG'],
            'measure': [24, 37, 15, 6, 44, 17, 29],
            'group': [7.0, 7.0, 9.0, 11.0, 11.0, 9.0, 7.0]
        })

        exp_ref_df = pd.DataFrame({
            'id': ['donor1', 'donor2', 'donor3', 'donor4', 'sampleB'],
            'measure': [32, 51, 3, 19, 37],
            'group': ['reference', 'reference', 'reference', 'reference', 'control1']
        })

        time_df, ref_df = group_timepoints(diversity_measure=self.alpha,
                             metadata=self.md_alpha,
                             time_column='days_post_transplant',
                             reference_column='relevant_donor',
                             control_column='single_control')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_alpha_dists_with_donors_no_controls(self):
        exp_time_df = pd.DataFrame({
            'id': ['sampleA', 'sampleB', 'sampleC', 'sampleD',
                   'sampleE', 'sampleF', 'sampleG'],
            'measure': [24, 37, 15, 6, 44, 17, 29],
            'group': [7.0, 7.0, 9.0, 11.0, 11.0, 9.0, 7.0],
        })

        exp_ref_df = pd.DataFrame({
            'id': ['donor1', 'donor2', 'donor3', 'donor4'],
            'measure': [32, 51, 3, 19],
            'group': ['reference', 'reference', 'reference', 'reference']
        })

        time_df, ref_df = group_timepoints(diversity_measure=self.alpha,
                                           metadata=self.md_alpha,
                                           time_column='days_post_transplant',
                                           reference_column='relevant_donor')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_alpha_dists_no_donors_with_controls(self):
        with self.assertRaisesRegex(TypeError, r"group_timepoints\(\) missing 1 required positional argument: "
                                    "'reference_column'"):
            group_timepoints(diversity_measure=self.alpha,
                             metadata=self.md_alpha,
                             time_column='days_post_transplant',
                             control_column='control')

    def test_alpha_dists_with_invalid_ref_column(self):
        with self.assertRaisesRegex(KeyError, 'References included in the metadata are missing'
                                    ' from the diversity measure.*foo.*bar.*baz'):
            group_timepoints(diversity_measure=self.alpha,
                             metadata=self.md_alpha,
                             time_column='days_post_transplant',
                             reference_column='invalid_ref_control',
                             control_column='control')

    def test_alpha_dists_with_empty_diversity_series(self):
        empty_alpha_series = pd.Series()

        with self.assertRaisesRegex(ValueError, 'Empty diversity measure detected'):
            group_timepoints(diversity_measure=empty_alpha_series,
                             metadata=self.md_alpha,
                             time_column='days_post_transplant',
                             reference_column='relevant_donor',
                             control_column='control')

    def test_alpha_dists_with_extra_samples_in_metadata_not_in_diversity(self):
        extra_md = Metadata.load(self.get_data_path('sample_metadata_alpha_div_missing.tsv'))

        exp_time_df = pd.DataFrame({
            'id': ['sampleA', 'sampleB', 'sampleC', 'sampleD',
                   'sampleE', 'sampleF', 'sampleG'],
            'measure': [24, 37, 15, 6, 44, 17, 29],
            'group': [7.0, 7.0, 9.0, 11.0, 11.0, 9.0, 7.0],
        })

        exp_ref_df = pd.DataFrame({
            'id': ['donor1', 'donor2', 'donor3', 'donor4',
                   'sampleC', 'sampleD', 'sampleE', 'sampleF'],
            'measure': [32, 51, 3, 19, 15, 6, 44, 17],
            'group': ['reference', 'reference', 'reference', 'reference',
                      'control1', 'control1', 'control2', 'control2']
        })

        time_df, ref_df = group_timepoints(diversity_measure=self.alpha,
                                           metadata=extra_md,
                                           time_column='days_post_transplant',
                                           reference_column='relevant_donor',
                                           control_column='control')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_alpha_dists_with_extra_samples_in_diversity_not_in_metadata(self):
        extra_alpha = pd.read_csv(self.get_data_path('alpha_div_missing.tsv'), sep='\t', index_col=0, squeeze=True)

        with self.assertRaisesRegex(ValueError, 'The following IDs are not present in the metadata'):
            group_timepoints(diversity_measure=extra_alpha,
                             metadata=self.md_alpha,
                             time_column='days_post_transplant',
                             reference_column='relevant_donor',
                             control_column='control')

    def test_examples(self):
        self.execute_examples()

class TestStats(TestBase):
    # Wilcoxon SRT test cases

    # Data in the exp_stats_data dataframes were pulled from Greg Caporaso's
    # Autism study repo on github, which can be found here:
    # https://github.com/caporaso-lab/autism-fmt1/blob/18-month-followup/16S/engraftment.ipynb
    def test_wilcoxon_with_faith_pd_baseline0_asymptotic(self):
        exp_stats_data = pd.DataFrame({
            'A:group': [0.0, 0.0, 0.0, 0.0],
            'A:n': [18, 18, 18, 18],
            'A:measure': [9.54973486, 9.54973486, 9.54973486, 9.54973486],
            'B:group': [3, 10, 18, 100],
            'B:n': [17, 18, 18, 16],
            'B:measure': [9.592979726, 10.9817719, 11.39392352, 12.97286672],
            'n': [17, 18, 18, 16],
            'test-statistic': [70.0, 20.0, 4.0, 5.0],
            'p-value': [0.758312374, 0.004337022, 0.000386178, 0.001123379],
            'q-value': [0.758312374, 0.005782696, 0.00154471, 0.002246758]
        })

        stats_data = wilcoxon_srt(distribution=self.faithpd_timedist, hypothesis='baseline',
                                  baseline_group='0', p_val_approx='asymptotic')

        pd.testing.assert_frame_equal(stats_data, exp_stats_data)

    def test_wilcoxon_with_faith_pd_consecutive_asymptotic(self):
        exp_stats_data = pd.DataFrame({
            'A:group': [0, 3, 10, 18],
            'A:n': [18, 17, 18, 18],
            'A:measure': [9.54973486, 9.592979726, 10.9817719, 11.393923515],
            'B:group': [3, 10, 18, 100],
            'B:n': [17, 18, 18, 16],
            'B:measure': [9.592980, 10.981772, 11.393924, 12.972867],
            'n': [17, 17, 18, 16],
            'test-statistic': [70.0, 26.0, 83.0, 24.0],
            'p-value': [0.758312, 0.016822, 0.913301, 0.022895],
            'q-value': [1.000000, 0.067288, 0.913301, 0.045790]
        })

        stats_data = wilcoxon_srt(distribution=self.faithpd_timedist,
                                  hypothesis='consecutive', p_val_approx='asymptotic')

        pd.testing.assert_frame_equal(stats_data, exp_stats_data)

    def test_wilcoxon_consecutive_hypothesis_with_baseline_group(self):
        with self.assertRaisesRegex(ValueError, "`consecutive` was selected as the hypothesis,"
                                    " but a `baseline_group` was added."):
            wilcoxon_srt(distribution=self.faithpd_timedist,
                         hypothesis='consecutive', baseline_group='reference')

    def test_wilcoxon_invalid_hypothesis(self):
        with self.assertRaisesRegex(ValueError, "Invalid hypothesis. Please either choose"
                                    " `baseline` or `consecutive` as your hypothesis."):
            wilcoxon_srt(distribution=self.faithpd_timedist, hypothesis='foo')

    def test_wilcoxon_invalid_baseline_group(self):
        with self.assertRaisesRegex(ValueError, "'foo' was not found as a group"
                                    " within the distribution."):
            wilcoxon_srt(distribution=self.faithpd_timedist,
                         hypothesis='baseline', baseline_group='foo')

    # Mann-Whitney U test cases

    # Data in the exp_stats_data dataframes were calculated 'by hand' in a jupyter notebook
    # using the same data, manually organized into groups and subsequently compared using
    # scipy.stats.mannwhitneyu to calculate the test-statistic and p-values
    # Notebook can be found here, for reference:
    # https://gist.github.com/lizgehret/c9add7b451e5e91b1017a2a963276bff
    def test_mann_whitney_pairwise_against_each(self):
        exp_stats_data = pd.DataFrame({
            'A:group': ['control', 'control', 'control', 'control', 'control',
                        'reference', 'reference', 'reference', 'reference', 'reference'],
            'A:n': [23, 23, 23, 23, 23, 5, 5, 5, 5, 5],
            'A:measure': [11.64962736, 11.64962736, 11.64962736, 11.64962736, 11.64962736,
                          10.24883918, 10.24883918, 10.24883918, 10.24883918, 10.24883918],
            'B:group': [0, 3, 10, 18, 100, 0, 3, 10, 18, 100],
            'B:n': [18, 17, 18, 18, 16, 18, 17, 18, 18, 16],
            'B:measure': [9.54973486, 9.592979726, 10.9817719, 11.39392352, 12.97286672,
                          9.54973486, 9.592979726, 10.9817719, 11.39392352, 12.97286672],
            'n': [41, 40, 41, 41, 39, 23, 22, 23, 23, 21],
            'test-statistic': [282.0, 260.0, 194.0, 190.0, 104.0,
                               49.0, 43.0, 20.0, 14.0, 6.0],
            'p-value': [0.050330911733538534, 0.07994303215567311, 0.7426248650660427,
                        0.6646800940267454, 0.02321456407322841, 0.7941892150565809,
                        1.0, 0.06783185968744732, 0.023005953105134484,
                        0.0056718704407604376],
            'q-value': [0.12582728, 0.13323839, 0.92828108, 0.94954299, 0.07738188,
                        0.88243246, 1.0, 0.13566372, 0.11502977, 0.0567187],
        })

        stats_data = mann_whitney_u(distribution=self.faithpd_refdist,
                                    against_each=self.faithpd_timedist,
                                    hypothesis='all-pairwise',
                                    p_val_approx='asymptotic')

        pd.testing.assert_frame_equal(stats_data, exp_stats_data)

    def test_mann_whitney_reference(self):
        exp_stats_data = pd.DataFrame({
            'A:group': ['reference'],
            'A:n': [5],
            'A:measure': [10.2488392],
            'B:group': ['control'],
            'B:n': [23],
            'B:measure': [11.6496274],
            'n': [28],
            'test-statistic': [37.0],
            'p-value': [0.23025583],
            'q-value': [0.23025583],
        })

        stats_data = mann_whitney_u(distribution=self.faithpd_refdist,
                                    hypothesis='reference',
                                    reference_group='reference',
                                    p_val_approx='asymptotic')

        pd.testing.assert_frame_equal(stats_data, exp_stats_data)

    def test_mann_whitney_all_pairwise_hypothesis_with_reference_group(self):
        with self.assertRaisesRegex(ValueError, "`all-pairwise` was selected as the"
                                    " hypothesis, but a `reference_group` was added."):
            mann_whitney_u(distribution=self.faithpd_refdist,
                           hypothesis='all-pairwise',
                           reference_group='reference')

    def test_mann_whitney_invalid_hypothesis(self):
        with self.assertRaisesRegex(ValueError, "Invalid hypothesis. Please either"
                         " choose `reference` or `all-pairwise` as your hypothesis."):
            mann_whitney_u(distribution=self.faithpd_refdist,
                           hypothesis='foo')

    def test_mann_whitney_invalid_reference_group(self):
        with self.assertRaisesRegex(ValueError, "'foo' was not found as a group"
                                    " within the distribution."):
            mann_whitney_u(distribution=self.faithpd_refdist,
                           hypothesis='reference', reference_group='foo')
