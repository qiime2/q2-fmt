# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from importlib_metadata import metadata
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
    def test_beta_dists_with_donors_controls(self):
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

    def test_beta_dists_with_non_numeric_time_column(self):
        with self.assertRaisesRegex(TypeError, 'Non-numeric characters detected in time_column'):
            group_timepoints(diversity_measure=self.dm,
                             metadata=self.md_beta,
                             time_column='non_numeric_time_column',
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

    def test_alpha_dists_with_non_numeric_time_column(self):
        with self.assertRaisesRegex(TypeError, 'Non-numeric characters detected in time_column'):
            group_timepoints(diversity_measure=self.alpha,
                             metadata=self.md_alpha,
                             time_column='non_numeric_time_column',
                             reference_column='relevant_donor',
                             control_column='control')

#TODO: edge cases

# no subject column
# no control
# control only 1 vs. more
# how can NaN behavior fail
# sometimes ref won't refer to samples - what if those values aren't in the table? (present in both id and ref cols)
# column param input not present in metadata
# ref/control param provided doesn't contain relevant data
# when samples are present in diversity but not metadata (error) & vice versa (ignored)
# when diversity series is empty (error)
