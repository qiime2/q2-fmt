# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

import qiime2


def _get_data_from_tests(path):
    return importlib.resources.files('q2_fmt.tests.data') / path


def alpha_md_factory():
    return qiime2.Metadata.load(
        _get_data_from_tests('sample_metadata_alpha_div.tsv'))


def beta_md_factory():
    return qiime2.Metadata.load(
        _get_data_from_tests('sample_metadata_donors.tsv'))


def alpha_div_factory():
    return qiime2.Artifact.import_data(
        'SampleData[AlphaDiversity]', _get_data_from_tests('alpha_div.tsv'))


def beta_div_factory():
    return qiime2.Artifact.import_data(
        'DistanceMatrix', _get_data_from_tests('dist_matrix_donors.tsv'))


def faithpd_md_factory():
    return qiime2.Metadata.load(
        _get_data_from_tests('metadata-faithpd.tsv')
    )


def faithpd_div_factory():
    return qiime2.Artifact.import_data(
        'SampleData[AlphaDiversity]', _get_data_from_tests('faithpd.tsv')
    )


def feature_table_factory():
    return qiime2.Artifact.import_data(
        'FeatureTable[Frequency]', _get_data_from_tests('feature-table.biom')
    )


def pedf_md_factory():
    return qiime2.Metadata.load(
        _get_data_from_tests('md-peds-usage.txt')
    )


def pedf_dist_factory():
    return qiime2.Artifact.import_data(
        "Dist1D[Ordered, Matched] % Properties('pedf')",
        _get_data_from_tests('peds_dist.table.jsonl')
    )


def group_timepoints_alpha_independent(use):
    alpha = use.init_artifact('alpha', alpha_div_factory)
    metadata = use.init_metadata('metadata', alpha_md_factory)

    timepoints, references = use.action(
        use.UsageAction('fmt', 'group_timepoints'),
        use.UsageInputs(
            diversity_measure=alpha,
            distance_to='donor',
            metadata=metadata,
            time_column='days_post_transplant',
            reference_column='relevant_donor',
            subject_column=False,
            group_column=False
        ),
        use.UsageOutputNames(
            timepoint_dists='timepoint_dists',
            reference_dists='reference_dists'
        )
    )

    timepoints.assert_output_type('Dist1D[Ordered, Independent]')
    references.assert_output_type('Dist1D[Unordered, Independent]')


def group_timepoints_beta(use):
    beta = use.init_artifact('beta', beta_div_factory)
    metadata = use.init_metadata('metadata', beta_md_factory)

    timepoints, references = use.action(
        use.UsageAction('fmt', 'group_timepoints'),
        use.UsageInputs(
            diversity_measure=beta,
            metadata=metadata,
            distance_to='donor',
            time_column='days_post_transplant',
            reference_column='relevant_donor',
            subject_column='subject',
            group_column=False
        ),
        use.UsageOutputNames(
            timepoint_dists='timepoint_dists',
            reference_dists='reference_dists'
        )
    )

    timepoints.assert_output_type('Dist1D[Ordered, Matched]')
    references.assert_output_type('Dist1D[Unordered, Independent]')


# Engraftment example using faith PD, baseline0 comparison
def cc_baseline(use):
    md = use.init_metadata('md', faithpd_md_factory)
    div_measure = use.init_artifact('div_measure', faithpd_div_factory)

    stats_table, raincloud = use.action(
        use.UsageAction('fmt', 'cc'),
        use.UsageInputs(
            diversity_measure=div_measure,
            metadata=md,
            distance_to='donor',
            compare='baseline',
            time_column='week',
            reference_column='InitialDonorSampleID',
            subject_column='SubjectID',
            where='SampleType="stool"',
            filter_missing_references=True,
            against_group='0',
            p_val_approx='asymptotic'
        ),
        use.UsageOutputNames(
            stats='stats',
            raincloud_plot='raincloud_plot'
        )
    )

    stats_table.assert_output_type('StatsTable[Pairwise]')
    raincloud.assert_output_type('Visualization')


def pedf_method(use):
    md = use.init_metadata('md', pedf_md_factory)
    table = use.init_artifact('table', feature_table_factory)

    pedf_group_dists, = use.action(
        use.UsageAction('fmt', 'pedf'),
        use.UsageInputs(
            table=table,
            metadata=md,
            time_column='time_point',
            reference_column='Donor',
            subject_column='SubjectID'
        ),
        use.UsageOutputNames(
            pedf_dists='pedf_dist'
        )

    )

    pedf_group_dists.assert_output_type("Dist1D[Ordered, Matched] %"
                                        " Properties('pedf')")


def prdf_method(use):
    md = use.init_metadata('md', pedf_md_factory)
    table = use.init_artifact('table', feature_table_factory)

    prdf_group_dists, = use.action(
        use.UsageAction('fmt', 'prdf'),
        use.UsageInputs(
            table=table,
            metadata=md,
            time_column='time_point',
            reference_column='Donor',
            subject_column='SubjectID'
        ),
        use.UsageOutputNames(
            prdf_dists='prdf_dist'
        )

    )

    prdf_group_dists.assert_output_type("Dist1D[Ordered, Matched] %"
                                        " Properties('prdf')")


def heatmap(use):
    pedf_dist = use.init_artifact('pedf_dist', pedf_dist_factory)

    pedf_heatmap, = use.action(
        use.UsageAction('fmt', 'heatmap'),
        use.UsageInputs(
            data=pedf_dist,
        ),
        use.UsageOutputNames(
            visualization='heatmap',

        )
    )

    pedf_heatmap.assert_output_type('Visualization')


def perm_pedf_method(use):
    md = use.init_metadata('md', pedf_md_factory)
    table = use.init_artifact('table', feature_table_factory)

    actual_sample_pedf, pedf_stats, global_stats = use.action(
        use.UsageAction('fmt', 'pedf_permutation_test'),
        use.UsageInputs(
            table=table,
            metadata=md,
            time_column='time_point',
            reference_column='Donor',
            subject_column='SubjectID',
            num_resamples=99,
            sampling_depth=400
        ),
        use.UsageOutputNames(
            actual_sample_pedf='actual_sample_pedf',
            per_subject_stats='per_subject_stats',
            global_stats='global_stats'
        )

    )

    pedf_stats.assert_output_type("StatsTable[Pairwise]")
    global_stats.assert_output_type("StatsTable[Pairwise]")


def pprf_method(use):
    md = use.init_metadata('md', pedf_md_factory)
    table = use.init_artifact('table', feature_table_factory)

    pprf_group_dists, = use.action(
        use.UsageAction('fmt', 'pprf'),
        use.UsageInputs(
            table=table,
            metadata=md,
            time_column='time_point',
            baseline_timepoint='1',
            subject_column='SubjectID',
            filter_missing_references=False
        ),
        use.UsageOutputNames(
            pprf_dists='pprf_dist'
        )

    )

    pprf_group_dists.assert_output_type("Dist1D[Ordered, Matched] %"
                                        " Properties('pprf')")


def detect_donor_indicators_method(use):
    md = use.init_metadata('md', pedf_md_factory)
    table = use.init_artifact('table', feature_table_factory)

    differentials, da_barplot = use.action(
        use.UsageAction('fmt', 'detect_donor_indicators'),
        use.UsageInputs(
            table=table,
            metadata=md,
            time_column='time_point',
            reference_column='Donor',
            baseline_timepoint='1',
        ),
        use.UsageOutputNames(
            differentials='differentials',
            da_barplot='da_barplot'
        )

    )

    differentials.assert_output_type("FeatureData[DifferentialAbundance]")
    da_barplot.assert_output_type("Visualization")
