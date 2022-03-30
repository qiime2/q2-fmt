# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from numpy import float64
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
        self.alpha = pd.read_csv(self.get_data_path('alpha_div.tsv'), sep='\t', index_col=0)

    def test_beta_dists_with_donors_controls(self):
        exp_time_df = pd.DataFrame(
            [['sampleA', '0.45', '7.0'],
            ['sampleB', '0.40', '7.0'],
            ['sampleC', '0.28', '9.0'],
            ['sampleD', '0.78', '11.0'],
            ['sampleE', '0.66', '11.0']],
            columns=['id', 'measure', 'group'],
            dtype=float64).set_index('id')

        exp_ref_df = pd.DataFrame(
            [['donor1..donor2', '0.24', 'reference', 'donor1', 'donor2'],
            ['donor1..donor3', '0.41', 'reference', 'donor1', 'donor3'],
            ['donor2..donor3', '0.74', 'reference', 'donor2', 'donor3'],
            ['sampleB..sampleC', '0.37', 'control1', 'sampleB', 'sampleC'],
            ['sampleB..sampleD', '0.44', 'control1', 'sampleB', 'sampleD'],
            ['sampleC..sampleD', '0.31', 'control1', 'sampleC', 'sampleD']],
            columns=['id', 'measure', 'group', 'A', 'B'],
            dtype=float64).set_index('id')

        time_df, ref_df = group_timepoints(diversity_measure=self.dm,
                               metadata=self.md_beta,
                               time_column='days_post_transplant',
                               reference_column='relevant_donor',
                               control_column='control')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_alpha_dists_with_donors_controls(self):
        exp_time_df = pd.DataFrame(
            [['sampleA', '24', '7'],
            ['sampleB', '37', '7'],
            ['sampleC', '15', '9'],
            ['sampleD', '6', '11'],
            ['sampleE', '44', '11'],
            ['sampleF', '17', '9'],
            ['sampleG', '29', '7']],
            columns=['id', 'measure', 'group'],
            dtype=float64).set_index('id')

        exp_ref_df = pd.DataFrame(
            [['donor1', '32', 'reference'],
            ['donor2', '51', 'reference'],
            ['donor4', '19', 'reference'],
            ['donor3', '3', 'reference']],
            columns=['id', 'measure', 'group'],
            dtype=float64).set_index('id')

        print(self.dm)
        print('.......')
        print(self.alpha)
        print('.......')
        print(self.md_alpha.to_dataframe())

        time_df, ref_df = group_timepoints(diversity_measure=self.alpha,
                               metadata=self.md_alpha,
                               time_column='days_post_transplant',
                               reference_column='relevant_donor',
                               control_column='control')

        print(time_df)
        print('..........')
        print(ref_df)

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)
