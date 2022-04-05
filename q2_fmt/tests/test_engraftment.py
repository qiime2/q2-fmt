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

class TestGroupTimepoints(TestPluginBase):
    package = 'q2_fmt.tests'

    def setUp(self):
        super().setUp()

        self.md_beta = Metadata.load(self.get_data_path('sample_metadata_donors.tsv'))
        self.md_alpha = Metadata.load(self.get_data_path('sample_metadata_alpha_div.tsv'))

        self.dm = DistanceMatrix.read(self.get_data_path('dist_matrix_donors.tsv')).to_series()
        self.alpha = pd.read_csv(self.get_data_path('alpha_div.tsv'), sep='\t', index_col=0, squeeze=True)

    # Beta Diversity (Distance Matrix) Test Cases
    def test_beta_dists_with_donors_and_controls(self):
        exp_time_df = pd.DataFrame({
            'id': ['sampleA', 'sampleB', 'sampleC', 'sampleD', 'sampleE'],
            'measure': [0.45, 0.40, 0.28, 0.78, 0.66],
            'group': [7.0, 7.0, 9.0, 11.0, 11.0]
        }).set_index('id')

        exp_ref_df = pd.DataFrame({
            'id': ['donor1..donor2', 'donor1..donor3', 'donor2..donor3',
                   'sampleB..sampleC', 'sampleB..sampleD', 'sampleC..sampleD'],
            'measure': [0.24, 0.41, 0.74, 0.37, 0.44, 0.31],
            'group': ['reference', 'reference', 'reference', 'control1', 'control1', 'control1'],
            'A': ['donor1', 'donor1', 'donor2', 'sampleB', 'sampleB', 'sampleC'],
            'B': ['donor2', 'donor3', 'donor3', 'sampleC', 'sampleD', 'sampleD']
        }).set_index('id')

        time_df, ref_df = group_timepoints(diversity_measure=self.dm,
                                           metadata=self.md_beta,
                                           time_column='days_post_transplant',
                                           reference_column='relevant_donor',
                                           control_column='control')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

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
        }).set_index('id')

        exp_ref_df = pd.DataFrame({
            'id': ['donor1..donor2', 'donor1..donor3', 'donor2..donor3'],
            'measure': [0.24, 0.41, 0.74],
            'group': ['reference', 'reference', 'reference'],
            'A': ['donor1', 'donor1', 'donor2'],
            'B': ['donor2', 'donor3', 'donor3']
        }).set_index('id')

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

    def test_beta_dists_with_non_numeric_time_column(self):
        with self.assertRaisesRegex(TypeError, 'Non-numeric characters detected in time_column'):
            group_timepoints(diversity_measure=self.dm,
                             metadata=self.md_beta,
                             time_column='non_numeric_time_column',
                             reference_column='relevant_donor',
                             control_column='control')

    def test_beta_dists_with_time_column_input_not_in_metadata(self):
        with self.assertRaisesRegex(ValueError, 'time_column provided not present within the metadata'):
            group_timepoints(diversity_measure=self.dm,
                             metadata=self.md_beta,
                             time_column='foo',
                             reference_column='relevant_donor',
                             control_column='control')

    def test_beta_dists_with_reference_column_input_not_in_metadata(self):
        with self.assertRaisesRegex(ValueError, 'reference_column provided not present within the metadata'):
            group_timepoints(diversity_measure=self.dm,
                             metadata=self.md_beta,
                             time_column='days_post_transplant',
                             reference_column='foo',
                             control_column='control')

    def test_beta_dists_with_control_column_input_not_in_metadata(self):
        with self.assertRaisesRegex(ValueError, 'control_column provided not present within the metadata'):
            group_timepoints(diversity_measure=self.dm,
                             metadata=self.md_beta,
                             time_column='days_post_transplant',
                             reference_column='relevant_donor',
                             control_column='foo')

    def test_beta_dists_with_invalid_ref_column(self):
        with self.assertRaisesRegex(KeyError, 'Pairwise comparisons were unsuccessful'):
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

    # Alpha Diversity (Series) Test Cases
    def test_alpha_dists_with_donors_controls(self):
        exp_time_df = pd.DataFrame({
            'id': ['sampleA', 'sampleB', 'sampleC', 'sampleD',
                   'sampleE', 'sampleF', 'sampleG'],
            'measure': [24, 37, 15, 6, 44, 17, 29],
            'group': [7.0, 7.0, 9.0, 11.0, 11.0, 9.0, 7.0],
        }).set_index('id')

        exp_ref_df = pd.DataFrame({
            'id': ['donor1', 'donor2', 'donor4', 'donor3',
                   'sampleC', 'sampleD', 'sampleE', 'sampleF'],
            'measure': [32, 51, 19, 3, 15, 6, 44, 17],
            'group': ['reference', 'reference', 'reference', 'reference',
                      'control1', 'control1', 'control2', 'control2']
        }).set_index('id')

        time_df, ref_df = group_timepoints(diversity_measure=self.alpha,
                                           metadata=self.md_alpha,
                                           time_column='days_post_transplant',
                                           reference_column='relevant_donor',
                                           control_column='control')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_alpha_dists_with_one_donor_and_controls(self):
        with self.assertRaisesRegex(ValueError, 'Missing references for the associated sample data'):
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
        }).set_index('id')

        exp_ref_df = pd.DataFrame({
            'id': ['donor1', 'donor2', 'donor4', 'donor3', 'sampleB'],
            'measure': [32, 51, 19, 3, 37],
            'group': ['reference', 'reference', 'reference', 'reference', 'control1']
        }).set_index('id')

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
        }).set_index('id')

        exp_ref_df = pd.DataFrame({
            'id': ['donor1', 'donor2', 'donor4', 'donor3'],
            'measure': [32, 51, 19, 3],
            'group': ['reference', 'reference', 'reference', 'reference']
        }).set_index('id')

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

    def test_alpha_dists_with_non_numeric_time_column(self):
        with self.assertRaisesRegex(TypeError, 'Non-numeric characters detected in time_column'):
            group_timepoints(diversity_measure=self.alpha,
                             metadata=self.md_alpha,
                             time_column='non_numeric_time_column',
                             reference_column='relevant_donor',
                             control_column='control')

    def test_alpha_dists_with_time_column_input_not_in_metadata(self):
        with self.assertRaisesRegex(ValueError, 'time_column provided not present within the metadata'):
            group_timepoints(diversity_measure=self.alpha,
                             metadata=self.md_alpha,
                             time_column='foo',
                             reference_column='relevant_donor',
                             control_column='control')

    def test_alpha_dists_with_reference_column_input_not_in_metadata(self):
        with self.assertRaisesRegex(ValueError, 'reference_column provided not present within the metadata'):
            group_timepoints(diversity_measure=self.alpha,
                             metadata=self.md_alpha,
                             time_column='days_post_transplant',
                             reference_column='foo',
                             control_column='control')

    def test_alpha_dists_with_control_column_input_not_in_metadata(self):
        with self.assertRaisesRegex(ValueError, 'control_column provided not present within the metadata'):
            group_timepoints(diversity_measure=self.alpha,
                             metadata=self.md_alpha,
                             time_column='days_post_transplant',
                             reference_column='relevant_donor',
                             control_column='foo')

    def test_alpha_dists_with_invalid_ref_column(self):
        with self.assertRaisesRegex(KeyError, 'Pairwise comparisons were unsuccessful'):
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

#TODO: edge cases

# when samples are present in diversity but not metadata (error) & vice versa (ignored)
# beta div - what happens when all (or a majority of) samples are associated with a single donor class (comparwisesons of a single donor)
