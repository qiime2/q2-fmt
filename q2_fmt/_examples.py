# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources

import qiime2


def _get_data_from_tests(path):
    return pkg_resources.resource_filename('q2_fmt.tests',
                                           os.path.join('data', path))


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


def faithpd_timedist_factory():
    return qiime2.Artifact.import_data(
        'GroupDist[Ordered, Matched]', _get_data_from_tests('faithpd_timedist')
    )


def faithpd_refdist_factory():
    return qiime2.Artifact.import_data(
        'GroupDist[Unordered, Independent]',
        _get_data_from_tests('faithpd_refdist')
    )


def faithpd_md_factory():
    return qiime2.Metadata.load(
        _get_data_from_tests('metadata-faithpd.tsv')
    )


def faithpd_div_factory():
    return qiime2.Artifact.import_data(
        'SampleData[AlphaDiversity]', _get_data_from_tests('faithpd.tsv')
    )


def group_timepoints_alpha_independent(use):
    alpha = use.init_artifact('alpha', alpha_div_factory)
    metadata = use.init_metadata('metadata', alpha_md_factory)

    timepoints, references = use.action(
        use.UsageAction('fmt', 'group_timepoints'),
        use.UsageInputs(
            diversity_measure=alpha,
            metadata=metadata,
            time_column='days_post_transplant',
            reference_column='relevant_donor',
            subject_column=False
        ),
        use.UsageOutputNames(
            timepoint_dists='timepoint_dists',
            reference_dists='reference_dists'
        )
    )

    timepoints.assert_output_type('GroupDist[Ordered, Independent]')
    references.assert_output_type('GroupDist[Unordered, Independent]')


def group_timepoints_beta(use):
    beta = use.init_artifact('beta', beta_div_factory)
    metadata = use.init_metadata('metadata', beta_md_factory)

    timepoints, references = use.action(
        use.UsageAction('fmt', 'group_timepoints'),
        use.UsageInputs(
            diversity_measure=beta,
            metadata=metadata,
            time_column='days_post_transplant',
            reference_column='relevant_donor',
            subject_column='subject',
        ),
        use.UsageOutputNames(
            timepoint_dists='timepoint_dists',
            reference_dists='reference_dists'
        )
    )

    timepoints.assert_output_type('GroupDist[Ordered, Matched]')
    references.assert_output_type('GroupDist[Unordered, Independent]')


def wilcoxon_baseline0(use):
    timedist = use.init_artifact('timedist', faithpd_timedist_factory)

    stats_table, = use.action(
        use.UsageAction('fmt', 'wilcoxon_srt'),
        use.UsageInputs(
            distribution=timedist,
            hypothesis='baseline',
            baseline_group='0',
            p_val_approx='asymptotic',
        ),
        use.UsageOutputNames(
            stats='stats'
        )
    )

    stats_table.assert_output_type('StatsTable[Pairwise]')


def mann_whitney_pairwise(use):
    timedist = use.init_artifact('timedist', faithpd_timedist_factory)
    refdist = use.init_artifact('refdist', faithpd_refdist_factory)

    stats_table, = use.action(
        use.UsageAction('fmt', 'mann_whitney_u'),
        use.UsageInputs(
            distribution=refdist,
            hypothesis='all-pairwise',
            against_each=timedist,
            p_val_approx='asymptotic',
        ),
        use.UsageOutputNames(
            stats='stats'
        )
    )

    stats_table.assert_output_type('StatsTable[Pairwise]')


# Engraftment example using faith PD, baseline0 hypothesis
def engraftment_baseline(use):
    md = use.init_metadata('md', faithpd_md_factory)
    div_measure = use.init_artifact('div_measure', faithpd_div_factory)

    stats_table, raincloud = use.action(
        use.UsageAction('fmt', 'engraftment'),
        use.UsageInputs(
            diversity_measure=div_measure,
            metadata=md,
            hypothesis='baseline',
            time_column='week',
            reference_column='InitialDonorSampleID',
            subject_column='SubjectID',
            where='SampleType="stool"',
            filter_missing_references=True,
            against_group='0',
            p_val_approx='asymptotic',
        ),
        use.UsageOutputNames(
            stats='stats',
            raincloud_plot='raincloud_plot'
        )
    )

    stats_table.assert_output_type('StatsTable[Pairwise]')
    raincloud.assert_output_type('Visualization')
