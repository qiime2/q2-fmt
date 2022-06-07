# ----------------------------------------------------------------------------
# Copyright (c) 2022-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from skbio.stats.distance import DistanceMatrix

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin import ValidationError
from qiime2 import Metadata

from q2_fmt._engraftment import group_timepoints
from q2_fmt._validator import (validate_all_dist_columns_present,
                               validate_unique_subjects_within_group)


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
                                 sep='\t', index_col=0, squeeze=True)


class ErrorMixins:
    def test_with_time_column_input_not_in_metadata(self):
        with self.assertRaisesRegex(ValueError,
                                    'time_column.*foo.*metadata'):
            group_timepoints(diversity_measure=self.div,
                             metadata=self.md,
                             time_column='foo',
                             reference_column='relevant_donor',
                             control_column='control')

    def test_with_reference_column_input_not_in_metadata(self):
        with self.assertRaisesRegex(ValueError,
                                    'reference_column.*foo.*metadata'):
            group_timepoints(diversity_measure=self.div,
                             metadata=self.md,
                             time_column='days_post_transplant',
                             reference_column='foo',
                             control_column='control')

    def test_with_control_column_input_not_in_metadata(self):
        with self.assertRaisesRegex(ValueError,
                                    'control_column.*foo.*metadata'):
            group_timepoints(diversity_measure=self.div,
                             metadata=self.md,
                             time_column='days_post_transplant',
                             reference_column='relevant_donor',
                             control_column='foo')

    def test_with_non_numeric_time_column(self):
        with self.assertRaisesRegex(ValueError,
                                    'time_column.*categorical.*numeric'):
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
            'group': ['reference', 'reference', 'reference',
                      'control1', 'control1', 'control1'],
            'A': ['donor1', 'donor1', 'donor2',
                  'sampleB', 'sampleB', 'sampleC'],
            'B': ['donor2', 'donor3', 'donor3',
                  'sampleC', 'sampleD', 'sampleD']
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
                                           time_column='days_post_transplant',
                                           reference_column='relevant_donor',
                                           control_column='control',
                                           subject_column='subject')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_beta_dists_with_same_donor_for_all_samples(self):
        _, ref_df = group_timepoints(diversity_measure=self.dm,
                                     metadata=self.md_beta,
                                     time_column='days_post_transplant',
                                     reference_column='relevant_donor_all')

        self.assertTrue(ref_df.empty)

    def test_beta_dists_with_one_donor_and_controls(self):
        with self.assertRaisesRegex(KeyError,
                                    'Missing references for the associated'
                                    ' sample data'):
            group_timepoints(diversity_measure=self.dm,
                             metadata=self.md_beta,
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
                                           time_column='days_post_transplant',
                                           reference_column='relevant_donor')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_beta_dists_no_donors_with_controls(self):
        with self.assertRaisesRegex(
            TypeError,
            r"group_timepoints\(\) missing 1 required positional argument: "
                "'reference_column'"):
            group_timepoints(diversity_measure=self.dm,
                             metadata=self.md_beta,
                             time_column='days_post_transplant',
                             control_column='control')

    def test_beta_dists_with_invalid_ref_column(self):
        with self.assertRaisesRegex(KeyError, 'References included in the'
                                    ' metadata are missing from the diversity'
                                    ' measure.*foo.*bar.*baz'):
            group_timepoints(diversity_measure=self.dm,
                             metadata=self.md_beta,
                             time_column='days_post_transplant',
                             reference_column='invalid_ref_control',
                             control_column='control')

    def test_beta_dists_with_empty_diversity_series(self):
        empty_beta_series = pd.Series()

        with self.assertRaisesRegex(ValueError,
                                    'Empty diversity measure detected'):
            group_timepoints(diversity_measure=empty_beta_series,
                             metadata=self.md_beta,
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
            'group': ['reference', 'control1',
                      'control1', 'control2', 'control2']
        })

        time_df, ref_df = group_timepoints(
            diversity_measure=self.alpha, metadata=self.md_alpha,
            time_column='days_post_transplant',
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
        with self.assertRaisesRegex(
            TypeError,
            r"group_timepoints\(\) missing 1 required positional argument: "
                "'reference_column'"):
            group_timepoints(diversity_measure=self.alpha,
                             metadata=self.md_alpha,
                             time_column='days_post_transplant',
                             control_column='control')

    def test_alpha_dists_with_invalid_ref_column(self):
        with self.assertRaisesRegex(KeyError, 'References included in the'
                                    ' metadata are missing from the diversity'
                                    ' measure.*foo.*bar.*baz'):
            group_timepoints(diversity_measure=self.alpha,
                             metadata=self.md_alpha,
                             time_column='days_post_transplant',
                             reference_column='invalid_ref_control',
                             control_column='control')

    def test_alpha_dists_with_empty_diversity_series(self):
        empty_alpha_series = pd.Series()

        with self.assertRaisesRegex(ValueError,
                                    'Empty diversity measure detected'):
            group_timepoints(diversity_measure=empty_alpha_series,
                             metadata=self.md_alpha,
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
                                           time_column='days_post_transplant',
                                           reference_column='relevant_donor',
                                           control_column='control')

        pd.testing.assert_frame_equal(time_df, exp_time_df)
        pd.testing.assert_frame_equal(ref_df, exp_ref_df)

    def test_alpha_dists_with_extra_samples_in_diversity_not_in_metadata(self):
        extra_alpha = pd.read_csv(self.get_data_path('alpha_div_missing.tsv'),
                                  sep='\t', index_col=0, squeeze=True)

        with self.assertRaisesRegex(ValueError, 'The following IDs are not'
                                    ' present in the metadata'):
            group_timepoints(diversity_measure=extra_alpha,
                             metadata=self.md_alpha,
                             time_column='days_post_transplant',
                             reference_column='relevant_donor',
                             control_column='control')

    def test_examples(self):
        self.execute_examples()


class TestValidators(TestBase):
    def test_validators_missing_columns_in_dist(self):
        with self.assertRaisesRegex(ValidationError, '"group" not found'
                                    ' in distribution.'):
            df = pd.DataFrame({
                'id': ['S340445', 'S892825', 'S460691'],
                'measure': [7.662921088, 8.431734297, 8.513263823]
            })
            validate_all_dist_columns_present(df, level=min)

    def test_validators_unique_subjects_not_duplicated_per_group(self):
        with self.assertRaisesRegex(ValidationError, 'Unique subject found'
                                    ' more than once within an individual'
                                    ' group.*0.*P26'):
            df = pd.DataFrame({
                'id': ['S116625', 'S813956'],
                'measure': [7.662921088, 8.431734297],
                'group': [0, 0],
                'subject': ['P26', 'P26']
            })
            validate_unique_subjects_within_group(df, level=min)
