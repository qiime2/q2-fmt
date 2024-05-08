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

from qiime2.plugin.testing import TestPluginBase
from qiime2 import Metadata

from q2_fmt._engraftment import group_timepoints
from q2_fmt._peds import (_compute_peds, sample_peds,
                          _filter_associated_reference,
                          _check_reference_column, _check_for_time_column,
                          _check_subject_column, _check_column_type,
                          _drop_incomplete_timepoints, feature_peds,
                          _check_column_missing, _rename_features, peds_bootstrap)



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


class TestPeds(TestBase):
    def test_get_donor(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor2', 'donor2', float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub1', 'sub2', 'sub2', float("Nan"),
                        float("Nan")],
            'group': [1, 2, 1, 2, float("Nan"),
                      float("Nan")]}).set_index('id')
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
            'Ref': ['donor1', 'donor1', 'donor2', 'donor2', float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub1', 'sub2', 'sub2', float("Nan"),
                        float("Nan")],
            'group': [1, 2, 1, 2, float("Nan"),
                      float("Nan")]}).set_index('id')
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
            'Ref': ['donor1', 'donor1', 'donor2', 'donor2', float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub1', 'sub2', 'sub2', float("Nan"),
                        float("Nan")],
            'group': [1, 2, 1, 2, float("Nan"),
                      float("Nan")]}).set_index('id')
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
            'Ref': [float("Nan"), float("Nan"), float("Nan"), float("Nan"),
                    float("Nan"), float("Nan")],
            'subject': ['sub1', 'sub1', 'sub2', 'sub2', float("Nan"),
                        float("Nan")],
            'group': [1, 2, 1, 2, float("Nan"),
                      float("Nan")]}).set_index('id')
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
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', float("Nan"),
                        float("Nan")],
            'group': [1, 2, 3, 2, float("Nan"),
                      float("Nan")]}).set_index('id')
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
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', float("Nan"),
                        float("Nan")],
            'group': [1, 2, 3, 2, float("Nan"),
                      float("Nan")]}).set_index('id')
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
                                     drop_incomplete_subjects=True)

        exp_peds_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3'],
            'measure': [0.666667, 0.333333, 1],
            'transfered_donor_features': [2, 1, 3],
            'total_donor_features': [3, 3, 3],
            'donor': ["donor1", "donor1", "donor1"],
            'subject': ["sub1", "sub1", "sub1"],
            'group': [1.0, 2.0, 3.0]
            })
        pd.testing.assert_frame_equal(sample_peds_df, exp_peds_df)

    def test_incorrect_reference_column_name(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', float("Nan"),
                        float("Nan")],
            'group': [1, 2, 3, 2, float("Nan"),
                      float("Nan")]}).set_index('id')
        with self.assertRaisesRegex(KeyError, ".*the provided"
                                    " `--p-reference-column`: `R` in the"
                                    " metadata"):
            _check_reference_column(metadata_df, "R")

    def test_incorrect_group_column_name(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', float("Nan"),
                        float("Nan")],
            'group': [1, 2, 3, 2, float("Nan"),
                      float("Nan")]}).set_index('id')
        with self.assertRaisesRegex(KeyError,
                                    ".*the provided `--p-time-column`: `time`"
                                    " in the metadata"):
            _check_for_time_column(metadata_df, 'time')

    def test_incorrect_subject_column_name(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', float("Nan"),
                        float("Nan")],
            'group': [1, 2, 3, 2, float("Nan"),
                      float("Nan")]}).set_index('id')
        with self.assertRaisesRegex(KeyError, ".*the provided"
                                    " `--p-subject-column`: `sub` in the"
                                    " metadata"):
            _check_subject_column(metadata_df, 'sub')

    def test_no_feature_overlap(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', float("Nan"),
                        float("Nan")],
            'group': [1, 2, 3, 2, float("Nan"),
                      float("Nan")]}).set_index('id')
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
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', float("Nan"),
                        float("Nan")],
            'group': [1, 2, 3, 2, float("Nan"),
                      float("Nan")]}).set_index('id')
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
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', float("Nan"),
                        float("Nan")],
            'group': [1, 2, 3, 2, float("Nan"),
                      float("Nan")]}).set_index('id')
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
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', float("Nan"),
                        float("Nan")],
            'group': [1, 2, 3, 2, float("Nan"),
                      float("Nan")]}).set_index('id')
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
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', float("Nan"),
                        float("Nan")],
            'group': [1, 2, 2, 1, float("Nan"),
                      float("Nan")]}).set_index('id')
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
            'Ref': ['donor1', 'donor1', 'donor1', float("Nan")],
            'subject': ['sub1', 'sub1', 'sub1', float("Nan")],
            'group': [1, 1, 1, float("Nan")]}).set_index('id')
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
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', float("Nan"),
                        float("Nan")],
            'group': [1, 2, 2, 1, float("Nan"),
                      float("Nan")]}).set_index('id')
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
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', float("Nan"),
                        float("Nan")],
            'group': ["t1", "t2", "t3", "t2", float("Nan"),
                      float("Nan")]}).set_index('id')
        metadata_obj = Metadata(metadata_df)
        column_properties = metadata_obj.columns
        with self.assertRaisesRegex(AssertionError, ".*Column with non-numeric"
                                    " values that was selected: `group`"):
            _check_column_type(column_properties, 'time', 'group', 'numeric')

    def test_column_type_noncategorical(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', float("Nan"),
                    float("Nan")],
            'subject': [1, 1, 1, 2, float("Nan"),
                        float("Nan")],
            'group': [1, 2, 3, 2, float("Nan"),
                      float("Nan")]}).set_index('id')
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
            'Ref': ['1', '1', '2', '2', float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub1', 'sub2', 'sub2', float("Nan"),
                        float("Nan")],
            'group': [1, 2, 1, 2, float("Nan"),
                      float("Nan")]}).set_index('id')
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
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', float("Nan"),
                        float("Nan")],
            'group': [1, 2, 3, 2, float("Nan"),
                      float("Nan")]}).set_index('id')
        with self.assertRaisesRegex(KeyError, ".*`--p-subject-column` can not"
                                    " be the same as the index of"
                                    " metadata: `id`"):
            _check_column_missing(metadata_df, 'id', 'subject', KeyError)

    def test_drop_incomplete_timepoints(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', float("Nan"),
                        float("Nan")],
            'group': [1, 2, 3, 2, float("Nan"),
                      float("Nan")]}).set_index('id')
        metadata_df = _drop_incomplete_timepoints(metadata_df, "group", [3])
        self.assertEqual(metadata_df["group"].unique()[0], float(1))
        self.assertEqual(metadata_df["group"].unique()[1], float(2))

    def test_drop_incomplete_timepoints_list(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'donor1', 'donor2'],
            'Ref': ['donor1', 'donor1', 'donor1', 'donor2', float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', float("Nan"),
                        float("Nan")],
            'group': [1, 2, 3, 2, float("Nan"),
                      float("Nan")]}).set_index('id')
        metadata_df = _drop_incomplete_timepoints(metadata_df, "group", [3, 2])
        self.assertEqual(metadata_df["group"].dropna().unique(), [float(1)])

    def test_rename_features_with_delim(self):
        metadata_df = pd.DataFrame({
                'id': ['sample1', 'sample2', 'sample3',
                       'donor1'],
                'Ref': ['donor1', 'donor1', 'donor1', float("Nan")],
                'subject': ['sub1', 'sub1', 'sub1', float("Nan")],
                'group': [1, 1, 1, float("Nan")]}).set_index('id')
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
                'Ref': ['donor1', 'donor1', 'donor1', float("Nan")],
                'subject': ['sub1', 'sub1', 'sub1', float("Nan")],
                'group': [1, 1, 1, float("Nan")]}).set_index('id')
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
                'Ref': ['donor1', 'donor1', 'donor1', float("Nan")],
                'subject': ['sub1', 'sub1', 'sub1', float("Nan")],
                'group': [1, 1, 1, float("Nan")]}).set_index('id')
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
                'Ref': ['donor1', 'donor1', 'donor1', float("Nan")],
                'subject': ['sub1', 'sub1', 'sub1', float("Nan")],
                'group': [1, 1, 1, float("Nan")]}).set_index('id')
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
        print(feature_peds_df)
        Fs1 = feature_peds_df.set_index("id").at['Feature 1 __',
                                                 'subject']
        Fs2 = feature_peds_df.set_index("id").at['Feature 2',
                                                 'subject']
        self.assertEqual("1", Fs1)
        self.assertEqual("2", Fs2)


class TestBoot(TestBase):
    def test_high_donor_overlap(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'donor1', 'donor2', 'donor3'],
            'Ref': ['donor1', 'donor2', 'donor3', float("Nan"), float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub2', 'sub3', float("Nan"), float("Nan"),
                        float("Nan")],
            'group': [1, 1, 1, float("Nan"), float("Nan"),
                      float("Nan")],
            "Location": [float("Nan"), float("Nan"),
                         float("Nan"), 'test', 'test',
                         'test']}).set_index('id')

        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'donor1', 'donor2', 'donor3'],
            'Feature1': [1, 0, 0, 1, 0, 0],
            'Feature2': [0, 1, 0, 0, 1, 0],
            'Feature3': [0, 0, 1, 0, 0, 1]}).set_index('id')

        p, real_mean, fake_mean = peds_bootstrap(metadata_df=metadata_df,
                                                 table=table_df,
                                                 time_column="group",
                                                 reference_column="Ref",
                                                 subject_column="subject",
                                                 iters=999)
        self.assertGreater(real_mean, fake_mean)

    def test_low_donor_overlap(self):
        metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'donor1', 'donor2', 'donor3'],
            'Ref': ['donor1', 'donor2', 'donor3', float("Nan"), float("Nan"),
                    float("Nan")],
            'subject': ['sub1', 'sub2', 'sub3', float("Nan"), float("Nan"),
                        float("Nan")],
            'group': [1, 1, 1, float("Nan"), float("Nan"),
                      float("Nan")],
            "Location": [float("Nan"), float("Nan"),
                         float("Nan"), 'test', 'test',
                         'test']}).set_index('id')

        table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'donor1', 'donor2', 'donor3'],
            'Feature1': [1, 0, 0, 0, 0, 1],
            'Feature2': [0, 1, 0, 1, 0, 0],
            'Feature3': [0, 0, 1, 0, 1, 0]}).set_index('id')

        p, real_median, fake_median = peds_bootstrap(metadata_df=metadata_df,
                                                     table=table_df,
                                                     time_column="group",
                                                     reference_column="Ref",
                                                     subject_column="subject",
                                                     iters=999)
        self.assertGreater(fake_median, real_median)

