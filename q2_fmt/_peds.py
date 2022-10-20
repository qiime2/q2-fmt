# ----------------------------------------------------------------------------
# Copyright (c) 2022-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2
import pandas as pd


def sample_peds(table: pd.DataFrame, metadata: qiime2.Metadata,
                time_column: str, reference_column: str, subject_column: str,
                filter_missing_references: bool = False,
                drop_incomplete_subjects: bool = False) -> (pd.DataFrame):
    # TODO: Make incomplete samples possible move this to heatmap
    metadata = metadata.to_dataframe()

    try:
        num_timepoints = metadata[time_column].dropna().unique().size
    except Exception as e:
        raise KeyError('There was an error finding %s in the metadata'
                       % time_column) from e
    try:
        metadata[time_column].dropna().astype(int)
    except Exception as e:
        raise KeyError('%s must be numeric' % time_column) from e
    try:
        subject_series = metadata[subject_column]
    except Exception as e:
        raise KeyError('There was an error finding %s in the metadata'
                       % subject_column) from e
    used_subject_series = subject_series[~metadata[time_column].isna()]
    subject_occurrence_df = (used_subject_series.value_counts()
                             .to_frame())
    if (subject_occurrence_df[subject_column] != num_timepoints).any():
        if drop_incomplete_subjects:
            subject_to_keep = (subject_occurrence_df
                               .loc[subject_occurrence_df[subject_column]
                                    == num_timepoints].index)
            metadata = metadata[metadata[subject_column].isin(subject_to_keep)]
        elif (subject_occurrence_df[subject_column] < num_timepoints).any():
            incomplete_subjects = (subject_occurrence_df
                                   .loc[subject_occurrence_df[subject_column]
                                        < num_timepoints].index).to_list()
            raise KeyError('Missing timepoints for associated subjects.'
                           ' Please make sure that all subjects have all'
                           ' timepoints or use drop_incomplete_subjects'
                           ' parameter. The incomplete subjects were %s'
                           % incomplete_subjects)
        elif (subject_occurrence_df[subject_column] > num_timepoints).any():
            raise KeyError('There are more occurrences of subjects than'
                           'timepoints. Make sure that all of your relevant'
                           'samples have associated timepoints.')

    try:
        reference_series = metadata[reference_column]
    except Exception as e:
        raise KeyError('There was an error finding %s in the metadata'
                       % reference_column) from e

    used_references = reference_series[~metadata[time_column].isna()]

    if used_references.isna().any():
        if filter_missing_references:
            used_references = used_references.dropna()
        else:
            nan_references = used_references.index[used_references.isna()]
            raise KeyError('Missing references for the associated sample data.'
                           ' Please make sure that all samples with a'
                           ' timepoint value have an associated reference.'
                           ' IDs where missing references were found:'
                           ' %s' % (tuple(nan_references),))

    peds_df = _compute_peds(used_references, table, metadata, time_column,
                            reference_column, subject_column)
    return peds_df


def _compute_peds(reference_series: pd.Series, table: pd.Series,
                  metadata: qiime2.Metadata, time_column: str,
                  reference_column: str,
                  subject_column: str) -> (pd.DataFrame):

    peds_series_list = []
    for sample in reference_series.index:
        donor = metadata.loc[sample][reference_column]
        subject = metadata.loc[sample][subject_column]
        timepoint = metadata.loc[sample][time_column]

        donor_present_list = _get_observed_features(table, donor)
        donor_num_present = _count_observed_features(donor_present_list)
        if donor_num_present == 0:
            raise ValueError('Donor Sample %s has no features in it.' % donor)
        sample_present_list = _get_observed_features(table, sample)

        intersect = (donor_present_list & sample_present_list).sum()
        intersect_num_present = _count_observed_features(intersect)

        peds = (intersect_num_present/donor_num_present)
        peds_series_list.append((sample, peds, intersect_num_present,
                                 donor_num_present, donor, subject, timepoint))
    peds_df = pd.DataFrame(peds_series_list,
                           columns=['id', 'measure',
                                    'transfered_donor_features',
                                    'total_donor_features', 'donor', 'subject',
                                    'group'])

    # use title for correcting ugly names
    peds_df['id'].attrs.update({
        'title': reference_series.index.name,
        'description': 'Sample IDs'
    })
    peds_df['measure'].attrs.update({
        'title': "PEDS",
        'description': 'Proportional Engraftment of Donor Strains '
    })
    peds_df['group'].attrs.update({
        'title': time_column,
        'description': 'Time'
    })
    peds_df["subject"].attrs.update({
        'title': subject_column,
        'description': 'ID to link samples across time'
    })
    peds_df["transfered_donor_features"].attrs.update({
        'title': "Transfered Donor Features",
        'description': '...'
    })
    peds_df['total_donor_features'].attrs.update({
        'title': "Total Donor Features",
        'description': '...'
    })
    peds_df['donor'].attrs.update({
        'title': reference_column,
        'description': 'Donor'
    })
    return peds_df


def _get_observed_features(table, id):
    try:
        present = table.loc[id] > 0
    except Exception as e:
        raise KeyError('There was an error finding the sample %s in'
                       ' your feature table' % id) from e
    return present


def _count_observed_features(present_list):
    num_present = present_list.sum()
    return num_present
