# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import pandas as pd
from skbio.stats.distance import DistanceMatrix

from qiime2.plugin.testing import TestPluginBase
from qiime2 import Metadata

from q2_fmt._engraftment import group_timepoints

class TestGroupTimepoints(TestPluginBase):
    package = 'q2_fmt.tests'

    # def dataframe_adds_blank_column(self):
    #     exp = pd.DataFrame(
    #             [['L1S8', 'AGCTGACTAGTC', 'gut', '2008', ''],
    #             ['L1S57', 'ACACACTATGGC', 'gut', '2009', '']],
    #             columns=['sample-id', 'barcode-sequence', 'body-site', 'year', 'ziggy'],
    #             dtype=object)

    #     obs = (TSVFileFormat('data/test_metadata.tsv'))

    #     self.assertEqual(exp, obs)

    def setUp(self):
        super().setUp()

        self.md_beta = Metadata.load(self.get_data_path('sample_metadata_donors.tsv'))
        self.md_alpha = Metadata.load(self.get_data_path('sample_metadata_alpha_div.tsv'))

        self.dm = DistanceMatrix.read(self.get_data_path('dist_matrix_donors.tsv'))
        self.alpha = pd.read_csv(self.get_data_path('alpha_div.tsv'),
                                 index_col=False, header=0).transpose()[0]

    def test_beta_timepoint_dist_with_donors_controls(self):
        exp_dm = pd.DataFrame(
            [['sampleA', '0.45', '7'],
            ['sampleB', '0.41', '7'],
            ['sampleC', '0.28', '9'],
            ['sampleD', '0.78', '11'],
            ['sampleE', '0.66', '11']],
            columns=['id', 'measure', 'group'],
            dtype=object)

        obs = group_timepoints(diversity_measure=self.dm,
                               metadata=self.md_beta,
                               time_column='days_post_transplant',
                               reference_column='relevant_donor',
                               control_column='control')

        print(obs)

    # def test_beta_ref_dist_with_donors_controls(self):
    #     exp = pd.DataFrame(
    #         [['sampleA', '0.450200535', '7'],
    #         ['sampleB', '0.409489887', '7'],
    #         ['sampleC', '0.28264371', '9'],
    #         ['sampleD', '0.78092451', '11'],
    #         ['sampleE', '0.66315908', '11']],
    #         columns=['id', 'measure', 'group'],
    #         dtype=object)

    #     obs = group_timepoints(diversity_measure=self.dm,
    #                            metadata=self.md_beta,
    #                            time_column='days_post_transplant',
    #                            reference_column='relevant_donor',
    #                            control_column='control')

    # def test_alpha_timepoint_dist_with_donors_controls(self):
    #     obs = group_timepoints(diversity_measure=self.alpha,
    #                            metadata=self.md_alpha,
    #                            time_column='days_post_transplant',
    #                            reference_column='relevant_donor',
    #                            control_column='control')

    # def test_alpha_ref_dist_with_donors_controls(self):
    #     obs = group_timepoints(diversity_measure=self.alpha,
    #                            metadata=self.md_alpha,
    #                            time_column='days_post_transplant',
    #                            reference_column='relevant_donor',
    #                            control_column='control')
