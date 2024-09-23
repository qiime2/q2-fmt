# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import numpy as np
from skbio.stats.distance import DistanceMatrix
from scipy.stats import false_discovery_control

from qiime2.plugin.testing import TestPluginBase
from qiime2 import Metadata

from q2_fmt._engraftment import group_timepoints, _get_to_baseline_ref
from q2_fmt._peds import (_compute_peds, sample_peds,
                          _filter_associated_reference,
                          _check_reference_column, _check_for_time_column,
                          _check_subject_column, _check_column_type,
                          _drop_incomplete_timepoints, feature_peds,
                          _check_column_missing, _rename_features,
                          peds_simulation, _create_mismatched_pairs,
                          _simulate_uniform_distro, _create_sim_masking,
                          _mask_recipient, _create_duplicated_recip_table,
                          _per_subject_stats, _global_stats, _peds_sim_stats,
                          sample_pprs)


class TestBase(TestPluginBase):
    package = 'q2_fmt.tests'

    def setUp(self):
        super().setUp()

        self.md_beta = Metadata.load(self.get_data_path(
                           'sample_metadata_donors.tsv'))
        self.md_alpha = Metadata.load(self.get_data_path(
                            'sample_metadata_alpha_div.tsv'))

        self.dm = DistanceMatrix.read(self.get_data_path(
                      'dist_matrix_donors.tsv')).to_series()
        self.alpha = pd.read_csv(self.get_data_path('alpha_div.tsv'),
                                 sep='\t', index_col=0).squeeze('columns')


class ErrorMixins:
    def test_with_time_column_input_not_in_metadata(self):
        with self.assertRaisesRegex(ValueError,
                                    'time_column.*foo.*metadata'):
            group_timepoints(diversity_measure=self.div,
                             metadata=self.md,
                             distance_to='donor',
                             time_column='foo',
                             reference_column='relevant_donor',
                             control_column='control')

    def test_with_reference_column_input_not_in_metadata(self):
        with self.assertRaisesRegex(ValueError,
                                    'reference_column.*foo.*metadata'):
            group_timepoints(diversity_measure=self.div,
                             metadata=self.md,
                             distance_to='donor',
                             time_column='days_post_transplant',
                             reference_column='foo',
                             control_column='control')

    def test_with_control_column_input_not_in_metadata(self):
        with self.assertRaisesRegex(ValueError,
                                    'control_column.*foo.*metadata'):
            group_timepoints(diversity_measure=self.div,
                             metadata=self.md,
                             distance_to='donor',
                             time_column='days_post_transplant',
                             reference_column='relevant_donor',
                             control_column='foo')

    def test_with_non_numeric_time_column(self):
        with self.assertRaisesRegex(ValueError,
                                    'time_column.*categorical.*numeric'):
            group_timepoints(diversity_measure=self.div,
                             metadata=self.md,
                             distance_to='donor',
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
            'group': ['reference', 'reference', 'reference',
                      'control1', 'control1', 'control1'],
            'A': ['donor1', 'donor1', 'donor2',
                  'sampleB', 'sampleB', 'sampleC'],
            'B': ['donor2', 'donor3', 'donor3',
                  'sampleC', 'sampleD', 'sampleD']
        })

        time_df, ref_df = group_timepoints(diversity_measure=self.dm,
                                           metadata=self.md_beta,
                                           distance_to='donor',
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
            'subject': ['subject1', 'subject2',
                        'subject1', 'subject1', 'subject2']
        })

        exp_ref_df = pd.DataFrame({
            'id': ['donor1..donor2', 'donor1..donor3', 'donor2..donor3',
                   'sampleB..sampleC', 'sampleB..sampleD', 'sampleC..sampleD'],
            'measure': [0.24, 0.41, 0.74, 0.37, 0.44, 0.31],
            'group': ['reference', 'reference', 'reference',
                      'control1', 'control1', 'control1'],
            'A': ['donor1', 'donor1', 'donor2',
                  'sampleB', 'sampleB', 'sampleC'],
            'B': ['donor2', 'donor3', 'donor3',
                  'sampleC', 'sampleD', 'sampleD']
        })

        time_df, ref_df = group_timepoints(diversity_measure=self.dm,
                                           metadata=self.md_beta,
                                           distance_to='donor',
                                           time_column='days_post_transplant',
                                           reference_column='relevant_donor',
                                           control_column='control',
                                           subject_column='subject')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_beta_dists_with_same_donor_for_all_samples(self):
        _, ref_df = group_timepoints(diversity_measure=self.dm,
                                     metadata=self.md_beta,
                                     distance_to='donor',
                                     time_column='days_post_transplant',
                                     reference_column='relevant_donor_all')

        self.assertTrue(ref_df.empty)

    def test_beta_dists_with_one_donor_and_controls(self):
        with self.assertRaisesRegex(KeyError,
                                    'Missing references for the associated'
                                    ' sample data'):
            group_timepoints(diversity_measure=self.dm,
                             metadata=self.md_beta,
                             distance_to='donor',
                             time_column='days_post_transplant',
                             reference_column='single_donor',
                             control_column='control')

    def test_beta_dists_with_donors_and_one_control(self):
        exp_ref_df = pd.DataFrame({
            'id': ['donor1..donor2', 'donor1..donor3', 'donor2..donor3'],
            'measure': [0.24, 0.41, 0.74],
            'group': ['reference', 'reference', 'reference'],
            'A': ['donor1', 'donor1', 'donor2'],
            'B': ['donor2', 'donor3', 'donor3']
        })

        _, ref_df = group_timepoints(diversity_measure=self.dm,
                                     metadata=self.md_beta,
                                     distance_to='donor',
                                     time_column='days_post_transplant',
                                     reference_column='relevant_donor',
                                     control_column='single_control')

        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

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
                                           distance_to='donor',
                                           time_column='days_post_transplant',
                                           reference_column='relevant_donor')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_beta_dists_no_donors_with_controls(self):
        with self.assertRaisesRegex(ValueError, "`donor` was provided to the"
                                    " `distance_to` parameter and a"
                                    " `reference_column` was not provided."
                                    " Please provide a `reference_column` *"):
            group_timepoints(diversity_measure=self.dm,
                             metadata=self.md_beta,
                             distance_to='donor',
                             time_column='days_post_transplant',
                             control_column='control')

    def test_beta_dists_with_invalid_ref_column(self):
        with self.assertRaisesRegex(KeyError, 'References included in the'
                                    ' metadata are missing from the diversity'
                                    ' measure.*foo.*bar.*baz'):
            group_timepoints(diversity_measure=self.dm,
                             metadata=self.md_beta,
                             distance_to='donor',
                             time_column='days_post_transplant',
                             reference_column='invalid_ref_control',
                             control_column='control')

    def test_beta_dists_with_empty_diversity_series(self):
        empty_beta_series = pd.Series()

        with self.assertRaisesRegex(ValueError,
                                    'Empty diversity measure detected'):
            group_timepoints(diversity_measure=empty_beta_series,
                             metadata=self.md_beta,
                             distance_to='donor',
                             time_column='days_post_transplant',
                             reference_column='relevant_donor',
                             control_column='control')

    def test_beta_dists_with_extra_samples_in_metadata_not_in_diversity(self):
        extra_md = Metadata.load(self.get_data_path(
                       'sample_metadata_donors_missing.tsv'))

        exp_time_df = pd.DataFrame({
            'id': ['sampleA', 'sampleB', 'sampleC', 'sampleD', 'sampleE'],
            'measure': [0.45, 0.40, 0.28, 0.78, 0.66],
            'group': [7.0, 7.0, 9.0, 11.0, 11.0]
        })

        exp_ref_df = pd.DataFrame({
            'id': ['donor1..donor2', 'donor1..donor3', 'donor2..donor3',
                   'sampleB..sampleC', 'sampleB..sampleD', 'sampleC..sampleD'],
            'measure': [0.24, 0.41, 0.74, 0.37, 0.44, 0.31],
            'group': ['reference', 'reference', 'reference',
                      'control1', 'control1', 'control1'],
            'A': ['donor1', 'donor1', 'donor2',
                  'sampleB', 'sampleB', 'sampleC'],
            'B': ['donor2', 'donor3', 'donor3',
                  'sampleC', 'sampleD', 'sampleD']
        })

        time_df, ref_df = group_timepoints(diversity_measure=self.dm,
                                           metadata=extra_md,
                                           distance_to='donor',
                                           time_column='days_post_transplant',
                                           reference_column='relevant_donor',
                                           control_column='control')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_beta_dists_with_extra_samples_in_diversity_not_in_metadata(self):
        extra_dm = DistanceMatrix.read(self.get_data_path(
                       'dist_matrix_donors_missing.tsv')).to_series()

        with self.assertRaisesRegex(ValueError,
                                    'The following IDs are not present'
                                    ' in the metadata'):
            group_timepoints(diversity_measure=extra_dm,
                             metadata=self.md_beta,
                             distance_to='donor',
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
                                           distance_to='donor',
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
                                           distance_to='donor',
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
            'group': ['reference', 'control1',
                      'control1', 'control2', 'control2']
        })

        time_df, ref_df = group_timepoints(
            diversity_measure=self.alpha, metadata=self.md_alpha,
            time_column='days_post_transplant',
            distance_to='donor',
            reference_column='relevant_donor_all',
            control_column='control')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_alpha_dists_with_one_donor_and_controls(self):
        with self.assertRaisesRegex(KeyError,
                                    'Missing references for the associated'
                                    ' sample data'):
            group_timepoints(diversity_measure=self.alpha,
                             metadata=self.md_alpha,
                             distance_to='donor',
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
            'group': ['reference', 'reference', 'reference',
                      'reference', 'control1']
        })

        time_df, ref_df = group_timepoints(diversity_measure=self.alpha,
                                           metadata=self.md_alpha,
                                           distance_to='donor',
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
                                           distance_to='donor',
                                           time_column='days_post_transplant',
                                           reference_column='relevant_donor')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_alpha_dists_no_donors_with_controls(self):
        with self.assertRaisesRegex(ValueError, "`donor` was provided to the"
                                    " `distance_to` parameter and a"
                                    " `reference_column` was not provided."
                                    " Please provide a `reference_column` *"):
            group_timepoints(diversity_measure=self.alpha,
                             metadata=self.md_alpha,
                             distance_to='donor',
                             time_column='days_post_transplant',
                             control_column='control')

    def test_alpha_dists_with_invalid_ref_column(self):
        with self.assertRaisesRegex(KeyError, 'References included in the'
                                    ' metadata are missing from the diversity'
                                    ' measure.*foo.*bar.*baz'):
            group_timepoints(diversity_measure=self.alpha,
                             metadata=self.md_alpha,
                             distance_to='donor',
                             time_column='days_post_transplant',
                             reference_column='invalid_ref_control',
                             control_column='control')

    def test_alpha_dists_with_empty_diversity_series(self):
        empty_alpha_series = pd.Series()

        with self.assertRaisesRegex(ValueError,
                                    'Empty diversity measure detected'):
            group_timepoints(diversity_measure=empty_alpha_series,
                             metadata=self.md_alpha,
                             distance_to='donor',
                             time_column='days_post_transplant',
                             reference_column='relevant_donor',
                             control_column='control')

    def test_alpha_dists_with_extra_samples_in_metadata_not_in_diversity(self):
        extra_md = Metadata.load(self.get_data_path(
                       'sample_metadata_alpha_div_missing.tsv'))

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
                                           distance_to='donor',
                                           time_column='days_post_transplant',
                                           reference_column='relevant_donor',
                                           control_column='control')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_alpha_dists_with_extra_samples_in_diversity_not_in_metadata(self):
        extra_alpha = pd.read_csv(self.get_data_path('alpha_div_missing.tsv'),
                                  sep='\t', index_col=0).squeeze('columns')

        with self.assertRaisesRegex(ValueError, 'The following IDs are not'
                                    ' present in the metadata'):
            group_timepoints(diversity_measure=extra_alpha,
                             metadata=self.md_alpha,
                             distance_to='donor',
                             time_column='days_post_transplant',
                             reference_column='relevant_donor',
                             control_column='control')

    def test_d2_donor_reference_col_baseline_tp(self):
        with self.assertRaisesRegex(ValueError, "`donor` was provided to the"
                                    " `distance_to` parameter and a value was"
                                    " provided to `baseline_timepoint`. These"
                                    " values can not be passed in together."):
            group_timepoints(diversity_measure=self.alpha,
                             metadata=self.md_alpha, distance_to='donor',
                             reference_column='relevant_donor',
                             baseline_timepoint="1",
                             time_column='days_post_transplant')

    def test_d2_donor_no_reference_col(self):
        with self.assertRaisesRegex(ValueError, "`donor` was provided to the"
                                    " `distance_to` parameter and a"
                                    " `reference_column` was not provided."
                                    " Please provide a `reference_column` *"):
            group_timepoints(diversity_measure=self.alpha,
                             metadata=self.md_alpha, distance_to='donor',
                             time_column='days_post_transplant')

    def test_d2_baseline_baseline_tp_ref_col(self):
        with self.assertRaisesRegex(ValueError, "`baseline` was provided to"
                                    " the `distance_to` parameter and a value"
                                    " was provided to `reference_column`. *"):
            group_timepoints(diversity_measure=self.alpha,
                             metadata=self.md_alpha, distance_to='baseline',
                             reference_column='relevant_donor',
                             baseline_timepoint="1",
                             time_column='days_post_transplant')

    def test_d2_baseline_no_baseline_tp(self):
        with self.assertRaisesRegex(ValueError, "`baseline` was provided to"
                                    " the `distance_to` parameter and a"
                                    " `baseline_timepoint` was not provided."):
            group_timepoints(diversity_measure=self.alpha,
                             metadata=self.md_alpha, distance_to='baseline',
                             time_column='days_post_transplant')

    def test_d2_baseline_alpha(self):
        md_baseline = Metadata.load(self.get_data_path(
                       'sample_metadata_a_div_baseline.txt'))

        exp_time_df = pd.DataFrame({
            'id': ['sampleD', 'sampleC', 'sampleE', 'sampleF'],
            'measure': [6, 15, 44, 17,],
            'group': [11.0, 9.0, 11.0, 9.0,],
            'subject': ['subject1', 'subject2', 'subject2', 'subject2']
        })

        exp_ref_df = pd.DataFrame({
            'id': ['sampleA', 'sampleB',
                   'sampleC', 'sampleD', 'sampleE', 'sampleF'],
            'measure': [24, 37, 15, 6, 44, 17],
            'group': ['reference', 'reference',
                      'control1', 'control1', 'control2', 'control2']
        })

        time_df, ref_df = group_timepoints(diversity_measure=self.alpha,
                                           metadata=md_baseline,
                                           distance_to='baseline',
                                           time_column='days_post_transplant',
                                           subject_column='subject',
                                           baseline_timepoint=7.0,
                                           control_column='control')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_d2_baseline_beta_dists(self):
        exp_time_df = pd.DataFrame({
            'id': ['sampleC', 'sampleD', 'sampleE'],
            'measure': [0.52, 0.51, 0.36],
            'group': [9.0, 11.0, 11.0],
            'subject': ['subject1', 'subject1', 'subject2']
        })

        exp_ref_df = pd.DataFrame({
            'id': ['sampleA..sampleB'],
            'measure': [0.48],
            'group': ['reference'],
            'A': ['sampleA'],
            'B': ['sampleB']
        })

        time_df, ref_df = group_timepoints(diversity_measure=self.dm,
                                           metadata=self.md_beta,
                                           distance_to='baseline',
                                           subject_column='subject',
                                           time_column='days_post_transplant',
                                           baseline_timepoint=7)

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_no_baseline_duplicates(self):
        with self.assertRaisesRegex(ValueError, "More than one baseline sample"
                                    " was found per subject.*"):
            group_timepoints(diversity_measure=self.alpha,
                             metadata=self.md_alpha, distance_to='baseline',
                             baseline_timepoint="7",
                             time_column='days_post_transplant',
                             subject_column='subject')

    def test_d2_baseline_alpha_missing_reference(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'sample4'],
            'subject': ['sub1', 'sub1', 'sub2', 'sub2'],
            'group': [1, 2, 2, 3]}).set_index('id')
        md_baseline = Metadata(metadata_df)

        obs_feature = pd.Series(data=[1, 0, 1, 0],
                                index=['sample1', 'sample2',
                                       'sample3', 'sample4'])
        with self.assertRaisesRegex(KeyError,
                                    'Missing references for the associated'
                                    ' sample data'):
            group_timepoints(diversity_measure=obs_feature,
                             metadata=md_baseline,
                             distance_to='baseline',
                             time_column='group',
                             subject_column='subject',
                             baseline_timepoint=1)

    def test_d2_baseline_alpha_filt_missing_reference(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'sample4'],
            'subject': ['sub1', 'sub1', 'sub2', 'sub2'],
            'group': [1, 2, 2, 3]}).set_index('id')
        md_baseline = Metadata(metadata_df)

        obs_feature = pd.Series(data=[1, 0, 1, 0],
                                index=['sample1', 'sample2',
                                       'sample3', 'sample4'])

        exp_time_df = pd.DataFrame({
            'id': ['sample2'],
            'measure': [0],
            'group': [2.0],
            'subject': ['sub1',]
        })

        exp_ref_df = pd.DataFrame({
            'id': ['sample1'],
            'measure': [1],
            'group': ['reference']
        })
        time_df, ref_df = group_timepoints(diversity_measure=obs_feature,
                                           metadata=md_baseline,
                                           distance_to='baseline',
                                           time_column='group',
                                           subject_column='subject',
                                           baseline_timepoint=1,
                                           filter_missing_references=True)

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_d2_baseline_alpha_drop_na_tp(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'sample4'],
            'subject': ['sub1', 'sub1', 'sub2', 'sub2'],
            'group': [1, 2, np.nan, np.nan]}).set_index('id')
        md_baseline = Metadata(metadata_df)

        obs_feature = pd.Series(data=[1, 0, 1, 0],
                                index=['sample1', 'sample2',
                                       'sample3', 'sample4'])

        exp_time_df = pd.DataFrame({
            'id': ['sample2'],
            'measure': [0],
            'group': [2.0],
            'subject': ['sub1',]
        })

        exp_ref_df = pd.DataFrame({
            'id': ['sample1'],
            'measure': [1],
            'group': ['reference']
        })
        time_df, ref_df = group_timepoints(diversity_measure=obs_feature,
                                           metadata=md_baseline,
                                           distance_to='baseline',
                                           time_column='group',
                                           subject_column='subject',
                                           baseline_timepoint=1,
                                           filter_missing_references=True)

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_d2_baseline_alpha_invalid_tp(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'sample4'],
            'subject': ['sub1', 'sub1', 'sub2', 'sub2'],
            'group': [1, 2, 2, 3]}).set_index('id')
        md_baseline = Metadata(metadata_df)

        obs_feature = pd.Series(data=[1, 0, 1, 0],
                                index=['sample1', 'sample2',
                                       'sample3', 'sample4'])
        with self.assertRaisesRegex(AssertionError,
                                    'The provided .* group.'):
            group_timepoints(diversity_measure=obs_feature,
                             metadata=md_baseline,
                             distance_to='baseline',
                             time_column='group',
                             subject_column='subject',
                             baseline_timepoint=7)

    def test_examples(self):
        self.execute_examples()

    def test_baseline_reference(self):
        time_col = pd.Series([1, 2, 3, 1, 2, 3], index=['sample1', 'sample2',
                                                        'sample3', 'sample4',
                                                        'sample5', 'sample6'])
        baseline_timepoint = "1"
        time_column = "group"
        subject_column = "subject"
        metadata = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'sample5', 'sample6'],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', 'sub2', 'sub2'],
            'group': [1, 2, 3, 1, 2, 3]}).set_index('id')
        ref = _get_to_baseline_ref(time_col, baseline_timepoint, time_column,
                                   subject_column, Metadata(metadata))

        exp_ref = pd.Series(['sample1', 'sample1', 'sample4', 'sample4'],
                            index=['sample2', 'sample3', 'sample5',
                                   'sample6'])
        exp_ref.index.name = 'sample_name'
        exp_ref.name = 'relevant_baseline'

        pd.testing.assert_series_equal(ref, exp_ref)


class TestPeds(TestBase):
    def test_get_donor(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor2', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub2', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 1, 2, np.nan,
                      np.nan]}).set_index('id')
        reference_series = metadata_df['Ref'].dropna()
        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Feature1': [1, 0, 1, 1, 1, 1],
            'Feature2': [1, 1, 1, 1, 1, 1],
            'Feature3': [0, 0, 1, 1, 1, 1]}).set_index('id')
        peds_df = pd.DataFrame(columns=['id', 'measure',
                                        'transfered_donor_features',
                                        'total_donor_features', 'donor',
                                        'subject', 'group'])
        peds_df = _compute_peds(peds_df=peds_df, peds_type="Sample",
                                peds_time=np.nan,
                                reference_series=reference_series,
                                table=table_df, metadata=metadata_df,
                                time_column="group", reference_column="Ref",
                                subject_column="subject")
        peds_df = peds_df.set_index("id")
        donor = peds_df.at["sample1", "donor"]
        self.assertEqual(donor, "donor1")

    def test_get_subject(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor2', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub2', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 1, 2, np.nan,
                      np.nan]}).set_index('id')
        reference_series = metadata_df['Ref'].dropna()
        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Feature1': [1, 0, 1, 1, 1, 1],
            'Feature3': [1, 1, 1, 1, 1, 1]}).set_index('id')
        peds_df = pd.DataFrame(columns=['id', 'measure',
                                        'transfered_donor_features',
                                        'total_donor_features', 'donor',
                                        'subject',
                                        'group'])
        peds_df = _compute_peds(peds_df=peds_df, peds_type="Sample",
                                peds_time=np.nan,
                                reference_series=reference_series,
                                table=table_df, metadata=metadata_df,
                                time_column="group", reference_column="Ref",
                                subject_column="subject")
        peds_df = peds_df.set_index("id")
        subject = peds_df.at["sample1", "subject"]
        self.assertEqual(subject, "sub1")

    def test_timepoint(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor2', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub2', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 1, 2, np.nan,
                      np.nan]}).set_index('id')
        reference_series = metadata_df['Ref'].dropna()
        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Feature1': [1, 0, 1, 1, 1, 1],
            'Feature2': [1, 1, 1, 1, 1, 1],
            'Feature3': [0, 0, 1, 1, 1, 1]}).set_index('id')
        peds_df = pd.DataFrame(columns=['id', 'measure',
                                        'transfered_donor_features',
                                        'total_donor_features', 'donor',
                                        'subject', 'group'])
        peds_df = _compute_peds(peds_df=peds_df, peds_type="Sample",
                                peds_time=np.nan,
                                reference_series=reference_series,
                                table=table_df, metadata=metadata_df,
                                time_column="group", reference_column="Ref",
                                subject_column="subject")
        peds_df = peds_df.set_index("id")
        tp = peds_df.at["sample3", "group"]
        self.assertEqual(tp, 1)

    def test_no_donors(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': [np.nan, np.nan, np.nan, np.nan,
                    np.nan, np.nan],
            'subject': ['sub1', 'sub1', 'sub2', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 1, 2, np.nan,
                      np.nan]}).set_index('id')
        reference_series = metadata_df['Ref']
        with self.assertRaisesRegex(KeyError, 'Missing references for'
                                    ' the associated sample data. Please make'
                                    ' sure that all samples with a timepoint'
                                    ' value have an associated reference.'
                                    ' IDs where missing references were found'
                                    ':.*'):
            _filter_associated_reference(reference_series=reference_series,
                                         metadata=metadata_df,
                                         time_column="group",
                                         filter_missing_references=False,
                                         reference_column="Ref")

    def test_incomplete_timepoints(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 3, 2, np.nan,
                      np.nan]}).set_index('id')
        metadata = Metadata(metadata_df)
        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Feature1': [1, 0, 1, 1, 1, 1],
            'Feature2': [1, 1, 1, 1, 1, 1]}).set_index('id')
        with self.assertRaisesRegex(ValueError, 'Missing timepoints for'
                                    ' associated subjects. Please make sure'
                                    ' that all subjects have all timepoints.'
                                    ' You can drop these subjects by using the'
                                    ' drop_incomplete_subjects parameter or'
                                    ' drop any timepoints that have large'
                                    ' numbers of subjects missing by using the'
                                    ' drop_incomplete_timepoints parameter. .*'
                                    '[\'sub2\']'):
            sample_peds(table=table_df, metadata=metadata,
                        time_column="group",
                        reference_column="Ref",
                        subject_column="subject")

    def test_incomplete_timepoints_with_flag(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub2', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 2, 3, np.nan,
                      np.nan]}).set_index('id')
        metadata = Metadata(metadata_df)
        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Feature1': [1, 0, 1, 1, 1, 1],
            'Feature2': [1, 1, 1, 1, 1, 1],
            'Feature3': [0, 0, 1, 1, 1, 1]}).set_index('id')
        sample_peds_df = sample_peds(table=table_df, metadata=metadata,
                                     time_column="group",
                                     reference_column="Ref",
                                     subject_column="subject",
                                     drop_incomplete_timepoints=[1, 3])

        exp_peds_df = pd.DataFrame({
            'id': ['sample2', 'sample3'],
            'measure': [0.333333, 1],
            'transfered_donor_features': [1, 3],
            'total_donor_features': [3, 3],
            'donor': ["donor1", "donor1"],
            'subject': ["sub1", "sub2"],
            'group': [2.0, 2.0]
            })
        pd.testing.assert_frame_equal(sample_peds_df, exp_peds_df)

    def test_incorrect_reference_column_name(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 3, 2, np.nan,
                      np.nan]}).set_index('id')
        with self.assertRaisesRegex(KeyError, ".*the provided"
                                    " `--p-reference-column`: `R` in the"
                                    " metadata"):
            _check_reference_column(metadata_df, "R")

    def test_incorrect_group_column_name(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 3, 2, np.nan,
                      np.nan]}).set_index('id')
        with self.assertRaisesRegex(KeyError,
                                    ".*the provided `--p-time-column`: `time`"
                                    " in the metadata"):
            _check_for_time_column(metadata_df, 'time')

    def test_incorrect_subject_column_name(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 3, 2, np.nan,
                      np.nan]}).set_index('id')
        with self.assertRaisesRegex(KeyError, ".*the provided"
                                    " `--p-subject-column`: `sub` in the"
                                    " metadata"):
            _check_subject_column(metadata_df, 'sub')

    def test_no_feature_overlap(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 3, 2, np.nan,
                      np.nan]}).set_index('id')
        metadata = Metadata(metadata_df)
        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Feature1': [0, 0, 1, 1, 1, 1],
            'Feature2': [0, 1, 1, 1, 1, 1],
            'Feature3': [0, 0, 1, 1, 1, 1]}).set_index('id')
        sample_peds_df = sample_peds(table=table_df, metadata=metadata,
                                     time_column="group",
                                     reference_column="Ref",
                                     subject_column="subject",
                                     drop_incomplete_subjects=True)
        exp_peds_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3'],
            'measure': [0, 0.333333, 1],
            'transfered_donor_features': [0, 1, 3],
            'total_donor_features': [3, 3, 3],
            'donor': ["donor1", "donor1", "donor1"],
            'subject': ["sub1", "sub1", "sub1"],
            'group': [1.0, 2.0, 3.0]
            })
        pd.testing.assert_frame_equal(sample_peds_df, exp_peds_df)

    def test_feature_overlap(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 3, 2, np.nan,
                      np.nan]}).set_index('id')
        metadata = Metadata(metadata_df)
        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Feature1': [0, 0, 1, 1, 1, 1],
            'Feature2': [0, 1, 1, 1, 1, 1],
            'Feature3': [0, 0, 1, 1, 1, 1]}).set_index('id')
        sample_peds_df = sample_peds(table=table_df, metadata=metadata,
                                     time_column="group",
                                     reference_column="Ref",
                                     subject_column="subject",
                                     drop_incomplete_subjects=True)
        TDFs1 = sample_peds_df.set_index("id").at['sample1',
                                                  'transfered_donor_features']
        TDFs2 = sample_peds_df.set_index("id").at['sample2',
                                                  'transfered_donor_features']
        TDFs3 = sample_peds_df.set_index("id").at['sample3',
                                                  'transfered_donor_features']
        self.assertEqual(TDFs2, 1)
        self.assertEqual(TDFs1, 0)
        self.assertEqual(TDFs3, 3)

    def test_peds_calc(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 3, 2, np.nan,
                      np.nan]}).set_index('id')
        metadata = Metadata(metadata_df)
        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Feature1': [0, 0, 1, 1, 1, 1],
            'Feature2': [0, 1, 1, 1, 1, 1],
            'Feature3': [0, 0, 1, 1, 1, 1]}).set_index('id')
        sample_peds_df = sample_peds(table=table_df, metadata=metadata,
                                     time_column="group",
                                     reference_column="Ref",
                                     subject_column="subject",
                                     drop_incomplete_subjects=True)
        TDFs1 = sample_peds_df.set_index("id").at['sample1',
                                                  'measure']
        TDFs2 = sample_peds_df.set_index("id").at['sample2',
                                                  'measure']
        TDFs3 = sample_peds_df.set_index("id").at['sample3',
                                                  'measure']
        self.assertEqual(TDFs2, 1/3)
        self.assertEqual(TDFs1, 0)
        self.assertEqual(TDFs3, 1)

# this test doesn't make sense anymore because of the refactor
    """ def test_no_feature_in_donor(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 3, 2, np.nan,
                      np.nan]}).set_index('id')
        metadata = Metadata(metadata_df)
        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Feature1': [1, 0, 1, 1, 0, 1],
            'Feature2': [1, 1, 1, 1, 0, 1],
            'Feature3': [0, 0, 1, 1, 0, 1]}).set_index('id')
        with self.assertRaisesRegex(ValueError, "Donor Sample donor1.*in it."):
            sample_peds(table=table_df, metadata=metadata,
                        time_column="group",
                        reference_column="Ref",
                        subject_column="subject",
                        drop_incomplete_subjects=True)
 """
    def test_unique_subjects_in_timepoints(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 2, 1, np.nan,
                      np.nan]}).set_index('id')
        metadata = Metadata(metadata_df)
        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Feature1': [1, 0, 1, 1, 1, 1],
            'Feature2': [1, 1, 1, 1, 1, 1],
            'Feature3': [0, 0, 1, 1, 1, 1]}).set_index('id')
        with self.assertRaisesRegex(ValueError, 'There is more than one'
                                    ' occurrence of.*Subject sub1.*[1,2,2]'):
            sample_peds(table=table_df, metadata=metadata,
                        time_column="group",
                        reference_column="Ref",
                        subject_column="subject",
                        drop_incomplete_subjects=True)

    def test_feature_peds_calc(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'donor1'],
            'Ref': ['donor1', 'donor1', 'donor1', np.nan],
            'subject': ['sub1', 'sub1', 'sub1', np.nan],
            'group': [1, 1, 1, np.nan]}).set_index('id')
        metadata = Metadata(metadata_df)
        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'donor1'],
            'Feature1': [0, 0, 1, 1],
            'Feature2': [0, 1, 1, 1],
            'Feature3': [0, 0, 1, 0]}).set_index('id')
        feature_peds_df = feature_peds(table=table_df, metadata=metadata,
                                       time_column="group",
                                       reference_column="Ref",
                                       subject_column="subject")
        TDFs1 = feature_peds_df.set_index("id").at['Feature1',
                                                   'measure']
        TDFs2 = feature_peds_df.set_index("id").at['Feature2',
                                                   'measure']
        self.assertEqual(TDFs1, 1/3)
        self.assertEqual(TDFs2, 2/3)

    def test_sample_id_match(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 2, 1, np.nan,
                      np.nan]}).set_index('id')
        metadata = Metadata(metadata_df)
        table_df = pd.DataFrame({
            'id': ['s1', 's2', 's3', 's4',
                   'd1', 'd2'],
            'Feature1': [1, 0, 1, 1, 1, 1],
            'Feature2': [1, 1, 1, 1, 1, 1],
            'Feature3': [0, 0, 1, 1, 1, 1]}).set_index('id')
        with self.assertRaisesRegex(ValueError, "The following IDs are not"
                                    " present in the metadata: 'd1', 'd2',"
                                    " 's1', 's2', 's3', 's4'"):
            feature_peds(table=table_df, metadata=metadata,
                         time_column="group",
                         reference_column="Ref",
                         subject_column="subject")

    def test_column_type_nonnumeric(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', np.nan,
                        np.nan],
            'group': ["t1", "t2", "t3", "t2", np.nan,
                      np.nan]}).set_index('id')
        metadata_obj = Metadata(metadata_df)
        column_properties = metadata_obj.columns
        with self.assertRaisesRegex(AssertionError, ".*Column with non-numeric"
                                    " values that was selected: `group`"):
            _check_column_type(column_properties, 'time', 'group', 'numeric')

    def test_column_type_noncategorical(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', np.nan,
                    np.nan],
            'subject': [1, 1, 1, 2, np.nan,
                        np.nan],
            'group': [1, 2, 3, 2, np.nan,
                      np.nan]}).set_index('id')
        metadata_obj = Metadata(metadata_df)
        column_properties = metadata_obj.columns
        with self.assertRaisesRegex(AssertionError, ".*Column with"
                                    " non-categorical values that was"
                                    " selected: `subject`"):
            _check_column_type(column_properties, 'subject', 'subject',
                               'categorical')

    def test_reference_series_not_in_table(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['1', '1', '2', '2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub2', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 1, 2, np.nan,
                      np.nan]}).set_index('id')
        reference_series = metadata_df['Ref'].dropna()
        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Feature1': [1, 0, 1, 1, 1, 1],
            'Feature2': [1, 1, 1, 1, 1, 1],
            'Feature3': [0, 0, 1, 1, 1, 1]}).set_index('id')
        peds_df = pd.DataFrame(columns=['id', 'measure',
                                        'transfered_donor_features',
                                        'total_donor_features', 'donor',
                                        'subject', 'group'])
        with self.assertRaisesRegex(AssertionError, ".*['1' '2'].*"):
            _compute_peds(peds_df=peds_df, peds_type="Sample",
                          peds_time=np.nan,
                          reference_series=reference_series,
                          table=table_df, metadata=metadata_df,
                          time_column="group", reference_column="Ref",
                          subject_column="subject")

    def test_column_name_is_ID(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 3, 2, np.nan,
                      np.nan]}).set_index('id')
        with self.assertRaisesRegex(KeyError, ".*`--p-subject-column` can not"
                                    " be the same as the index of"
                                    " metadata: `id`"):
            _check_column_missing(metadata_df, 'id', 'subject', KeyError)

    def test_drop_incomplete_timepoints(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 3, 2, np.nan,
                      np.nan]}).set_index('id')
        metadata_df = _drop_incomplete_timepoints(metadata_df, "group", [3])
        self.assertEqual(metadata_df["group"].unique()[0], float(1))
        self.assertEqual(metadata_df["group"].unique()[1], float(2))

    def test_drop_incomplete_timepoints_list(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 3, 2, np.nan,
                      np.nan]}).set_index('id')
        metadata_df = _drop_incomplete_timepoints(metadata_df, "group", [3, 2])
        self.assertEqual(metadata_df["group"].dropna().unique(), [float(1)])

    def test_rename_features_with_delim(self):
        metadata_df = pd.DataFrame({
                'id': ['sample1', 'sample2', 'sample3',
                       'donor1'],
                'Ref': ['donor1', 'donor1', 'donor1', np.nan],
                'subject': ['sub1', 'sub1', 'sub1', np.nan],
                'group': [1, 1, 1, np.nan]}).set_index('id')
        metadata = Metadata(metadata_df)
        table_df = pd.DataFrame({
                'id': ['sample1', 'sample2', 'sample3',
                       'donor1'],
                'Feature;1': [0, 0, 1, 1],
                'Feature;2': [0, 1, 1, 1],
                'Feature;3': [0, 0, 1, 0]}).set_index('id')
        feature_peds_df = feature_peds(table=table_df, metadata=metadata,
                                       time_column="group",
                                       reference_column="Ref",
                                       subject_column="subject")
        _rename_features(data=feature_peds_df, level_delimiter=";")
        Fs1 = feature_peds_df.set_index("id").at['Feature 1',
                                                 'subject']
        Fs2 = feature_peds_df.set_index("id").at['Feature 2',
                                                 'subject']
        self.assertEqual("1", Fs1)
        self.assertEqual("2", Fs2)

    def test_rename_features_with_no_delim(self):
        metadata_df = pd.DataFrame({
                'id': ['sample1', 'sample2', 'sample3',
                       'donor1'],
                'Ref': ['donor1', 'donor1', 'donor1', np.nan],
                'subject': ['sub1', 'sub1', 'sub1', np.nan],
                'group': [1, 1, 1, np.nan]}).set_index('id')
        metadata = Metadata(metadata_df)
        table_df = pd.DataFrame({
                'id': ['sample1', 'sample2', 'sample3',
                       'donor1'],
                'Feature1': [0, 0, 1, 1],
                'Feature2': [0, 1, 1, 1],
                'Feature3': [0, 0, 1, 0]}).set_index('id')
        feature_peds_df = feature_peds(table=table_df, metadata=metadata,
                                       time_column="group",
                                       reference_column="Ref",
                                       subject_column="subject")
        _rename_features(data=feature_peds_df, level_delimiter=None)
        Fs1 = feature_peds_df.set_index("id").at['Feature1',
                                                 'subject']
        Fs2 = feature_peds_df.set_index("id").at['Feature2',
                                                 'subject']
        self.assertEqual("Feature1", Fs1)
        self.assertEqual("Feature2", Fs2)

    def test_rename_features_with_wrong_delim(self):
        metadata_df = pd.DataFrame({
                'id': ['sample1', 'sample2', 'sample3',
                       'donor1'],
                'Ref': ['donor1', 'donor1', 'donor1', np.nan],
                'subject': ['sub1', 'sub1', 'sub1', np.nan],
                'group': [1, 1, 1, np.nan]}).set_index('id')
        metadata = Metadata(metadata_df)
        table_df = pd.DataFrame({
                'id': ['sample1', 'sample2', 'sample3',
                       'donor1'],
                'Feature;1': [0, 0, 1, 1],
                'Feature;2': [0, 1, 1, 1],
                'Feature;3': [0, 0, 1, 0]}).set_index('id')
        feature_peds_df = feature_peds(table=table_df, metadata=metadata,
                                       time_column="group",
                                       reference_column="Ref",
                                       subject_column="subject")
        _rename_features(data=feature_peds_df, level_delimiter=":")
        Fs1 = feature_peds_df.set_index("id").at['Feature;1',
                                                 'subject']
        Fs2 = feature_peds_df.set_index("id").at['Feature;2',
                                                 'subject']
        self.assertEqual("Feature;1", Fs1)
        self.assertEqual("Feature;2", Fs2)

    def test_rename_features_with_blank_label(self):
        metadata_df = pd.DataFrame({
                'id': ['sample1', 'sample2', 'sample3',
                       'donor1'],
                'Ref': ['donor1', 'donor1', 'donor1', np.nan],
                'subject': ['sub1', 'sub1', 'sub1', np.nan],
                'group': [1, 1, 1, np.nan]}).set_index('id')
        metadata = Metadata(metadata_df)
        table_df = pd.DataFrame({
                'id': ['sample1', 'sample2', 'sample3',
                       'donor1'],
                'Feature;1;__': [0, 0, 1, 1],
                'Feature;2': [0, 1, 1, 1],
                'Feature;3': [0, 0, 1, 0]}).set_index('id')
        feature_peds_df = feature_peds(table=table_df, metadata=metadata,
                                       time_column="group",
                                       reference_column="Ref",
                                       subject_column="subject")
        _rename_features(data=feature_peds_df, level_delimiter=";")
        Fs1 = feature_peds_df.set_index("id").at['Feature 1 __',
                                                 'subject']
        Fs2 = feature_peds_df.set_index("id").at['Feature 2',
                                                 'subject']
        self.assertEqual("1", Fs1)
        self.assertEqual("2", Fs2)

    def test_peds_nan_tp(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, np.nan, 2, np.nan,
                      np.nan]}).set_index('id')
        metadata = Metadata(metadata_df)
        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Feature1': [0, 0, 1, 1, 1, 1],
            'Feature2': [0, 1, 1, 1, 1, 1],
            'Feature3': [0, 0, 1, 1, 1, 1]}).set_index('id')
        sample_peds_df = sample_peds(table=table_df, metadata=metadata,
                                     time_column="group",
                                     reference_column="Ref",
                                     subject_column="subject",
                                     drop_incomplete_subjects=True)
        obs_samples = sample_peds_df['id'].to_list()
        exp_sample = ['sample1', 'sample2']
        self.assertEqual(obs_samples, exp_sample)

    def test_peds_no_donor_in_table(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor2', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub2', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 1, 2, np.nan,
                      np.nan]}).set_index('id')
        metadata = Metadata(metadata_df)
        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1'],
            'Feature1': [0, 0, 1, 1, 1],
            'Feature2': [0, 1, 1, 1, 1],
            'Feature3': [0, 0, 1, 1, 1]}).set_index('id')
        with self.assertRaisesRegex(KeyError, "References included in the"
                                    " metadata are missing from the feature"
                                    " table.*"):
            sample_peds(table=table_df, metadata=metadata,
                        time_column="group",
                        reference_column="Ref",
                        subject_column="subject",
                        drop_incomplete_subjects=True)

    def test_peds_no_donor_in_table_flag(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor2', 'donor2', np.nan,
                    np.nan],
            'subject': ['sub1', 'sub1', 'sub2', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 1, 2, np.nan,
                      np.nan]}).set_index('id')
        metadata = Metadata(metadata_df)
        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1'],
            'Feature1': [0, 0, 1, 1, 1],
            'Feature2': [0, 1, 1, 1, 1],
            'Feature3': [0, 0, 1, 1, 1]}).set_index('id')
        sample_peds_df = sample_peds(table=table_df, metadata=metadata,
                                     time_column="group",
                                     reference_column="Ref",
                                     subject_column="subject",
                                     drop_incomplete_subjects=True, 
                                     filter_missing_references=True)
        obs_samples = sample_peds_df['id'].to_list()
        exp_sample = ['sample1', 'sample2']
        self.assertEqual(obs_samples, exp_sample)

    def test_pprs(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'sample5', 'sample6'],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', 'sub2', 'sub2'],
            'group': [1, 2, 3, 1, 2, 3]}).set_index('id')
        metadata = Metadata(metadata_df)
        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'sample5', 'sample6'],
            'Feature1': [1, 0, 1, 0, 0, 0],
            'Feature2': [0, 0, 1, 1, 0, 1]}).set_index('id')
        sample_pprs_df = sample_pprs(table=table_df, metadata=metadata,
                                     time_column="group",
                                     subject_column="subject",
                                     baseline_timepoint="1")

        exp_pprs_df = pd.DataFrame({
            'id': ['sample2', 'sample3',  'sample5', 'sample6'],
            'measure': [0.0, 1.0, 0.0, 1.0],
            'transfered_baseline_features': [0, 1, 0, 1],
            'total_baseline_features': [1, 1, 1, 1],
            'baseline': ["sample1", "sample1", "sample4", "sample4"],
            'subject': ["sub1", "sub1", "sub2", "sub2"],
            'group': [2.0, 3.0, 2.0, 3.0]
            })
        pd.testing.assert_frame_equal(sample_pprs_df, exp_pprs_df)

    def test_pprs_incomplete_timepoints_with_flag(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'pre1', 'pre2'],
            'subject': ['sub1', 'sub1', 'sub2', 'sub2', 'sub1',
                        'sub2'],
            'group': [1, 2, 2, 3, 0,
                      0]}).set_index('id')
        metadata = Metadata(metadata_df)
        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'pre1', 'pre2'],
            'Feature1': [1, 0, 1, 1, 1, 1],
            'Feature2': [1, 1, 1, 1, 1, 1],
            'Feature3': [0, 0, 1, 1, 1, 1]}).set_index('id')
        sample_pprs_df = sample_pprs(table=table_df, metadata=metadata,
                                     time_column="group",
                                     baseline_timepoint="0",
                                     subject_column="subject",
                                     drop_incomplete_timepoints=[1, 3])

        exp_pprs_df = pd.DataFrame({
            'id': ['sample2', 'sample3'],
            'measure': [0.333333, 1],
            'transfered_baseline_features': [1, 3],
            'total_baseline_features': [3, 3],
            'baseline': ["pre1", "pre2"],
            'subject': ["sub1", "sub2"],
            'group': [2.0, 2.0]
            })
        pd.testing.assert_frame_equal(sample_pprs_df, exp_pprs_df)

    def test_pprs_baseline_sub_incomplete_timepoints_with_flag(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'pre1', 'pre2'],
            'subject': ['sub1', 'sub1', 'sub2', 'sub2', np.nan,
                        np.nan],
            'group': [1, 2, 2, 3, 0,
                      0]}).set_index('id')
        metadata = Metadata(metadata_df)
        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'pre1', 'pre2'],
            'Feature1': [1, 0, 1, 1, 1, 1],
            'Feature2': [1, 1, 1, 1, 1, 1],
            'Feature3': [0, 0, 1, 1, 1, 1]}).set_index('id')
        with self.assertRaisesRegex(AssertionError, "No baseline samples"
                                    " were connected via subject. .*"):
            sample_pprs(table=table_df, metadata=metadata,
                        time_column="group",
                        baseline_timepoint="0",
                        subject_column="subject",
                        drop_incomplete_timepoints=[1, 3])


class TestSim(TestBase):
    def test_high_donor_overlap(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'donor1', 'donor2', 'donor3'],
            'Ref': ['donor1', 'donor2', 'donor3', np.nan, np.nan,
                    np.nan],
            'subject': ['sub1', 'sub2', 'sub3', np.nan, np.nan,
                        np.nan],
            'group': [1, 1, 1, np.nan, np.nan,
                      np.nan],
            "Location": [np.nan, np.nan,
                         np.nan, 'test', 'test',
                         'test']}).set_index('id')

        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'donor1', 'donor2', 'donor3'],
            'Feature1': [1, 0, 0, 1, 0, 0],
            'Feature2': [0, 1, 0, 0, 1, 0],
            'Feature3': [0, 0, 1, 0, 0, 1]}).set_index('id')
        metadata = Metadata(metadata_df)

        stats, _ = peds_simulation(metadata=metadata,
                                   table=table_df,
                                   time_column="group",
                                   reference_column="Ref",
                                   subject_column="subject",
                                   num_iterations=999)
        real_median = np.median(stats["A:measure"].values)
        fake_median = np.median(stats["B:measure"].values)
        self.assertGreater(real_median, fake_median)

    def test_low_donor_overlap(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'donor1', 'donor2', 'donor3'],
            'Ref': ['donor1', 'donor2', 'donor3', np.nan, np.nan,
                    np.nan],
            'subject': ['sub1', 'sub2', 'sub3', np.nan, np.nan,
                        np.nan],
            'group': [1, 1, 1, np.nan, np.nan,
                      np.nan],
            "Location": [np.nan, np.nan,
                         np.nan, 'test', 'test',
                         'test']}).set_index('id')

        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'donor1', 'donor2', 'donor3'],
            'Feature1': [1, 0, 0, 0, 1, 1],
            'Feature2': [0, 1, 0, 1, 0, 1],
            'Feature3': [0, 0, 1, 1, 1, 0]}).set_index('id')
        metadata = Metadata(metadata_df)

        stats, _ = peds_simulation(metadata=metadata,
                                   table=table_df,
                                   time_column="group",
                                   reference_column="Ref",
                                   subject_column="subject",
                                   num_iterations=999)

        real_median = np.median(stats["A:measure"].values)
        fake_median = np.median(stats["B:measure"].values)
        self.assertGreater(fake_median, real_median)

    def test_single_donor(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'donor1'],
            'Ref': ['donor1', 'donor1', 'donor1', np.nan],
            'subject': ['sub1', 'sub2', 'sub3', np.nan],
            'group': [1, 1, 1, np.nan],
            "Location": [np.nan, np.nan,
                         np.nan, 'test']}).set_index('id')

        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'donor1'],
            'Feature1': [1, 0, 0, 0],
            'Feature2': [0, 1, 0, 1],
            'Feature3': [0, 0, 1, 1]}).set_index('id')
        metadata = Metadata(metadata_df)

        with self.assertRaisesRegex(AssertionError, "There is only one"
                                    " donated microbiome in your data. *"):
            peds_simulation(metadata=metadata,
                            table=table_df,
                            time_column="group",
                            reference_column="Ref",
                            subject_column="subject",
                            num_iterations=999)

    def test_create_mismatched_pairs(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'donor1', 'donor2', 'donor3'],
            'Ref': ['donor1', 'donor2', 'donor3', np.nan, np.nan,
                    np.nan],
            'subject': ['sub1', 'sub2', 'sub3', np.nan, np.nan,
                        np.nan],
            'group': [1, 1, 1, np.nan, np.nan,
                      np.nan],
            "Location": [np.nan, np.nan,
                         np.nan, 'test', 'test',
                         'test']}).set_index('id')

        recip_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3'],
            'Feature1': [1, 0, 0],
            'Feature2': [0, 1, 0],
            'Feature3': [0, 0, 1]}).set_index('id')

        used_references = pd.Series(data=['donor1', 'donor2', 'donor3'],
                                    index=['sample1', 'sample2', 'sample3'],
                                    name='Ref')
        used_references.index.name = "id"
        mismatched_df = _create_mismatched_pairs(recip_df, metadata_df,
                                                 used_references,
                                                 reference_column='Ref')
        exp_mismatched_df = pd.DataFrame({'id': ["sample1", "sample1",
                                                 "sample2", "sample2",
                                                 "sample3", "sample3"],
                                          "Ref": ["donor2", "donor3",
                                                  "donor1", "donor3",
                                                  "donor1", "donor2"]}
                                         ).set_index('id')
        pd.testing.assert_frame_equal(mismatched_df, exp_mismatched_df)

    def test_mask_recipient(self):
        donor_mask = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
        recip_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3'],
            'Feature1': [1, 0, 0],
            'Feature2': [0, 1, 0],
            'Feature3': [0, 0, 1]}).set_index('id')
        recip_mask = _mask_recipient(donor_mask, recip_df)
        exp_r_mask = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        np.testing.assert_array_equal(recip_mask, exp_r_mask)

    def test_mask_recipient_donor_one(self):
        donor_mask = np.array([[0, 1, 0], [0, 1, 0], [0, 1, 0]])
        recip_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3'],
            'Feature1': [1, 0, 0],
            'Feature2': [0, 1, 0],
            'Feature3': [0, 0, 1]}).set_index('id')
        recip_mask = _mask_recipient(donor_mask, recip_df)
        exp_r_mask = [[0, 0, 0], [0, 1, 0], [0, 0, 0]]
        np.testing.assert_array_equal(recip_mask, exp_r_mask)

    def test_simulate_uniform_distro(self):
        # Note: This tests has a VERY small chance to have intermit failures
        # if by random chance 1, 2, or 3 are not selected by random.choice.
        mismatch_peds = [1, 2, 3]

        iterations = 999

        mismatchpairs_df = _simulate_uniform_distro(mismatch_peds,
                                                    iterations)
        self.assertIn(1, mismatchpairs_df)
        self.assertIn(2, mismatchpairs_df)
        self.assertIn(3, mismatchpairs_df)
        self.assertEquals(mismatchpairs_df.size, iterations)

    def test_one_iter_simulate_uniform_distro(self):
        mismatch_peds = [0, 0, 0, 0, 0, 0]

        iterations = 1

        mismatchpairs_df = _simulate_uniform_distro(mismatch_peds,
                                                    iterations)
        self.assertEquals(mismatchpairs_df.size, iterations)

    def test_create_sim_masking(self):

        mismatched_df = pd.DataFrame({'id': ["sample1", "sample1",
                                             "sample2", "sample2",
                                             "sample3", "sample3"],
                                      "Ref": ["donor2", "donor3",
                                              "donor1", "donor3",
                                              "donor1", "donor2"]}
                                     ).set_index('id')

        donor_df = pd.DataFrame({
            'id': ['donor1', 'donor2', 'donor3'],
            'Feature1': [1, 0, 0],
            'Feature2': [0, 1, 0],
            'Feature3': [0, 0, 1]}).set_index('id')

        exp_mask = [[0, 1, 0],
                    [0, 0, 1],
                    [1, 0, 0],
                    [0, 0, 1],
                    [1, 0, 0],
                    [0, 1, 0]]

        donor_mask = _create_sim_masking(mismatched_df, donor_df,
                                         reference_column='Ref')
        np.testing.assert_array_equal(donor_mask, exp_mask)

    def test_create_one_donor_sim_masking(self):

        mismatched_df = pd.DataFrame({'id': ["sample1",
                                             "sample2",
                                             "sample3"],
                                      "Ref": ["donor2",
                                              "donor2",
                                              "donor2"]}
                                     ).set_index('id')

        donor_df = pd.DataFrame({
            'id': ['donor2'],
            'Feature1': [1],
            'Feature2': [0],
            'Feature3': [0]}).set_index('id')

        exp_mask = [[1, 0, 0],
                    [1, 0, 0],
                    [1, 0, 0]]

        donor_mask = _create_sim_masking(mismatched_df, donor_df,
                                         reference_column='Ref')
        np.testing.assert_array_equal(donor_mask, exp_mask)

    def test_create_duplicated_table(self):
        recip_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3'],
            'Feature1': [1, 0, 0],
            'Feature2': [0, 1, 0],
            'Feature3': [0, 0, 1]}).set_index('id')

        mismatched_df = pd.DataFrame({'id': ["sample1", "sample1",
                                             "sample2", "sample2",
                                             "sample3", "sample3"],
                                      "Ref": ["donor2", "donor3",
                                              "donor1", "donor3",
                                              "donor1", "donor2"]}
                                     ).set_index('id')

        duplicated_recip_table = _create_duplicated_recip_table(mismatched_df,
                                                                recip_df)

        exp_d_r_table = pd.DataFrame({
            'id': ['sample1', 'sample1', 'sample2', 'sample2',
                   'sample3', 'sample3'],
            'Feature1': [1, 1, 0, 0, 0, 0],
            'Feature2': [0, 0, 1, 1, 0, 0],
            'Feature3': [0, 0, 0, 0, 1, 1]}).set_index('id')

        pd.testing.assert_frame_equal(duplicated_recip_table,
                                      exp_d_r_table)

    def test_create_no_duplicated_table(self):
        recip_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3'],
            'Feature1': [1, 0, 0],
            'Feature2': [0, 1, 0],
            'Feature3': [0, 0, 1]}).set_index('id')

        mismatched_df = pd.DataFrame({'id': ["sample1",
                                             "sample2",
                                             "sample3"],
                                      "Ref": ["donor2",
                                              "donor2",
                                              "donor2"]}
                                     ).set_index('id')

        duplicated_recip_table = _create_duplicated_recip_table(mismatched_df,
                                                                recip_df)
        pd.testing.assert_frame_equal(duplicated_recip_table,
                                      recip_df)

    def test_per_subject_stats_labels(self):
        mismatched_peds = [0, 0, 0, 0]
        actual_temp = pd.Series(data=[1, 1, 1, 1],
                                index=["sample1", "sample2", "sample3",
                                       "sample4"])
        iterations = 10

        p_s_stats = _per_subject_stats(mismatched_peds,
                                       actual_temp, iterations)

        exp_column_names = ["A:group", "A:n", "A:measure",
                            "B:group", "B:n", "B:measure", "n",
                            "test-statistic", "p-value", "q-value"]
        np.testing.assert_array_equal(p_s_stats.columns.values,
                                      exp_column_names)

    def test_per_subject_stats(self):
        mismatched_peds = [0, 0, 0, 0]
        actual_temp = pd.Series(data=[1, 1, 1, 1],
                                index=["sample1", "sample2", "sample3",
                                       "sample4"])
        iterations = 10

        p_s_stats = _per_subject_stats(mismatched_peds,
                                       actual_temp, iterations)

        exp_test_stats = pd.Series([10, 10, 10, 10])

        np.testing.assert_array_equal(p_s_stats["test-statistic"].values,
                                      exp_test_stats.values)

    def test_per_subject_stats_q(self):
        mismatched_peds = [0, 0, 0, 0]
        actual_temp = pd.Series(data=[1, 0, 1, 0],
                                index=["sample1", "sample2", "sample3",
                                       "sample4"])
        iterations = 10

        p_s_stats = _per_subject_stats(mismatched_peds,
                                       actual_temp, iterations)

        exp_p = ([1/11, 11/11, 1/11, 11/11])

        exp_q = false_discovery_control(ps=exp_p, method='bh')

        np.testing.assert_array_equal(p_s_stats["q-value"].values,
                                      exp_q)

    def test_global_stats_label(self):
        p_series = pd.Series(data=[0.001, 0.001, 0.001, 0.001])

        p = _global_stats(p_series)

        exp_labels = ["Measure", "n", "test-statistic", "p-value", "q-value"]
        np.testing.assert_array_equal(p.columns.values,
                                      exp_labels)

    def test_global_stats(self):
        p_series = pd.Series(data=[0.001, 0.001, 0.001, 0.001])

        p = _global_stats(p_series)

        np.testing.assert_array_equal(p["n"].values, [4])

    def test_peds_sim_stats_good_match(self):
        value = 1
        peds_iters = pd.Series(data=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                               index=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        num_iterations = 10

        count_gte, count_less, per_subject_p = _peds_sim_stats(value,
                                                               peds_iters,
                                                               num_iterations)
        exp_count_gte = 0
        exp_count_less = 10
        exp_per_subject_p = (1/11)

        self.assertEqual(count_gte, exp_count_gte)
        self.assertEqual(count_less, exp_count_less)
        self.assertEqual(per_subject_p, exp_per_subject_p)

    def test_peds_sim_stats_bad_match(self):
        value = 0
        peds_iters = pd.Series(data=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                               index=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        num_iterations = 10

        count_gte, count_less, per_subject_p = _peds_sim_stats(value,
                                                               peds_iters,
                                                               num_iterations)
        exp_count_gte = 10
        exp_count_less = 0
        exp_per_subject_p = (11/11)

        self.assertEqual(count_gte, exp_count_gte)
        self.assertEqual(count_less, exp_count_less)
        self.assertEqual(per_subject_p, exp_per_subject_p)

    def test_peds_sim_stats_equal_match(self):
        value = .5
        peds_iters = pd.Series(data=[.5, .5, .5, .5, .5, .5, .5, .5, .5, .5],
                               index=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        num_iterations = 10

        count_gte, count_less, per_subject_p = _peds_sim_stats(value,
                                                               peds_iters,
                                                               num_iterations)
        exp_count_gte = 10
        exp_count_less = 0
        exp_per_subject_p = (11/11)

        self.assertEqual(count_gte, exp_count_gte)
        self.assertEqual(count_less, exp_count_less)
        self.assertEqual(per_subject_p, exp_per_subject_p)

    def test_peds_sim_stats_50_percent_bad_match(self):
        value = .5
        peds_iters = pd.Series(data=[.5, 1, .5, 1, .5, 1, .5, 1, .5, 1],
                               index=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        num_iterations = 10

        count_gte, count_less, per_subject_p = _peds_sim_stats(value,
                                                               peds_iters,
                                                               num_iterations)
        exp_count_gte = 10
        exp_count_less = 0
        exp_per_subject_p = (11/11)

        self.assertEqual(count_gte, exp_count_gte)
        self.assertEqual(count_less, exp_count_less)
        self.assertEqual(per_subject_p, exp_per_subject_p)

    def test_peds_sim_stats_50_good_match(self):
        value = .5
        peds_iters = pd.Series(data=[.5, 0, .5, 0, .5, 0, .5, 0, .5, 0],
                               index=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        num_iterations = 10

        count_gte, count_less, per_subject_p = _peds_sim_stats(value,
                                                               peds_iters,
                                                               num_iterations)
        exp_count_gte = 5
        exp_count_less = 5
        exp_per_subject_p = (6/11)

        self.assertEqual(count_gte, exp_count_gte)
        self.assertEqual(count_less, exp_count_less)
        self.assertEqual(per_subject_p, exp_per_subject_p)

    def test_peds_sim_stats_99_iters(self):
        value = .5
        peds_iters = pd.Series(data=[0, .5, 0, .5, 0, .5, 0, .5, 0, .5,
                                     0, .5, 0, .5, 0, .5, 0, .5, 0, .5, 0, .5,
                                     0, .5, 0, .5, 0, .5, 0, .5, 0, .5, 0, .5,
                                     0, .5, 0, .5, 0, .5, 0, .5, 0, .5, 0, .5,
                                     0, .5, 0, .5, 0, .5, 0, .5, 0, .5, 0, .5,
                                     0, .5, 0, .5, 0, .5, 0, .5, 0, .5, 0, .5,
                                     0, .5, 0, .5, 0, .5, 0, .5, 0, .5, 0, .5,
                                     0, .5, 0, .5, 0, .5, 0, .5, 0, .5, 0, .5,
                                     0, .5, 0, .5, 0], index=list(range(99)))
        num_iterations = 99

        count_gte, count_less, per_subject_p = _peds_sim_stats(value,
                                                               peds_iters,
                                                               num_iterations)
        exp_count_gte = 49
        exp_count_less = 50
        exp_per_subject_p = (50/100)

        self.assertEqual(count_gte, exp_count_gte)
        self.assertEqual(count_less, exp_count_less)
        self.assertEqual(per_subject_p, exp_per_subject_p)
