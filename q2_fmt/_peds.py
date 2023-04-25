# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2
import pandas as pd
import numpy as np
import warnings


def sample_peds(table: pd.DataFrame, metadata: qiime2.Metadata,
                time_column: str, reference_column: str, subject_column: str,
                filter_missing_references: bool = False,
                drop_incomplete_subjects: bool = False) -> (pd.DataFrame):
    ids_with_data = table.index
    metadata = metadata.filter_ids(ids_to_keep=ids_with_data)
    # TODO: Make incomplete samples possible move this to heatmap
    metadata = metadata.to_dataframe()
    num_timepoints = _check_for_time_column(metadata, time_column)
    _check_time_column_numeric(metadata, time_column)
    reference_series = _check_reference_column(metadata, reference_column)
    # return things that should be removed
    metadata, used_references = \
        _check_associated_reference(reference_series, metadata, time_column,
                                    filter_missing_references,
                                    reference_column)
    subject_series = _check_subject_column(metadata, subject_column)
    _check_duplicate_subject_timepoint(subject_series, metadata,
                                        subject_column, time_column)
    # return things that should be removed
    metadata, used_references = \
        _check_subjects_in_all_timepoint(subject_series, num_timepoints,
                                         drop_incomplete_subjects, metadata,
                                         subject_column, used_references)

    peds_df = pd.DataFrame(columns=['id', 'measure',
                                    'transfered_donor_features',
                                    'total_donor_features', 'donor', 'subject',
                                    'group'])
    peds_df = _compute_peds(peds_df=peds_df, peds_type="Sample",
                            peds_time=np.nan, reference_series=used_references,
                            table=table, metadata=metadata,
                            time_column=time_column,
                            subject_column=subject_column,
                            reference_column=reference_column)
    return peds_df


def feature_peds(table: pd.DataFrame, metadata: qiime2.Metadata,
                 time_column: str, reference_column: str, subject_column: str,
                 filter_missing_references: bool = False) -> (pd.DataFrame):
    ids_with_data = table.index
    metadata = metadata.filter_ids(ids_to_keep=ids_with_data)
    metadata = metadata.to_dataframe()

    _ = _check_for_time_column(metadata, time_column)
    _check_time_column_numeric(metadata, time_column)
    reference_series = _check_reference_column(metadata, reference_column)
    metadata, used_references = \
        _check_associated_reference(reference_series, metadata, time_column,
                                    filter_missing_references,
                                    reference_column)
    peds_df = pd.DataFrame(columns=['id', 'measure', 'recipients with feature',
                                    'all possible recipients with feature',
                                    'group', 'subject'])
    for time, time_metadata in metadata.groupby(time_column):
        peds_df = _compute_peds(peds_df=peds_df, peds_type="Feature",
                                peds_time=time,
                                reference_series=used_references, table=table,
                                metadata=time_metadata,
                                time_column=time_column,
                                subject_column=subject_column,
                                reference_column=reference_column)
    return peds_df


def _compute_peds(peds_df: pd.Series, peds_type: str, peds_time: int,
                  reference_series: pd.Series, table: pd.Series,
                  metadata: qiime2.Metadata, time_column: str,
                  subject_column: str,
                  reference_column: str) -> (pd.DataFrame):
    table = table > 0
    donordf = table[table.index.isin(reference_series)]
    recipdf = _create_recipient_table(reference_series, metadata, table)

    donormask = _create_masking(time_metadata=metadata, donordf=donordf,
                                recipdf=recipdf,
                                reference_column=reference_column)
    maskedrecip = donormask & recipdf
    if peds_type == "Sample":
        num_sum = np.sum(maskedrecip, axis=1)
        donor_sum = np.sum(donormask, axis=1)
        for count, sample in enumerate(recipdf.index):
            sample_row = metadata.loc[sample]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                peds = num_sum[count] / donor_sum[count]

            peds_df.loc[len(peds_df)] = [sample, peds, num_sum[count],
                                         donor_sum[count],
                                         sample_row[reference_column],
                                         sample_row[subject_column],
                                         sample_row[time_column]]
        peds_df['id'].attrs.update({
            'title': reference_series.index.name,
            'description': 'Sample IDs'
        })
        peds_df['measure'].attrs.update({
            'title': "PEDS",
            'description': 'Proportional Engraftment of Donor Strains'
        })
        peds_df['group'].attrs.update({
            'title': time_column,
            'description': 'Time'
        })
        peds_df["subject"].attrs.update({
            'title': subject_column,
            'description': 'Subject IDs linking samples across time'
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

    elif peds_type == "Feature":
        num_sum = np.sum(maskedrecip, axis=0)
        donor_sum = np.sum(donormask, axis=0)
        for count, feature in enumerate(recipdf.columns):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                peds = num_sum[count] / donor_sum[count]
            peds_df.loc[len(peds_df)] = [feature, peds, num_sum[count],
                                         donor_sum[count], peds_time, feature]
            peds_df = peds_df.dropna()
        peds_df['id'].attrs.update({
            'title': "FeatureID",
            'description': ''
        })
        peds_df['measure'].attrs.update({
            'title': "PEDS",
            'description': 'Proportional Engraftment of Donor Strains'
        })
        peds_df['group'].attrs.update({
            'title': time_column,
            'description': 'Time'
        })
        peds_df['subject'].attrs.update({
            'title': "Feature ID",
            'description': ''
        })
    else:
        raise KeyError('There was an error finding which PEDS methods to use')
    return peds_df


# Filtering methods
def _check_for_time_column(metadata, time_column):
    try:
        num_timepoints = metadata[time_column].dropna().unique().size
    except Exception as e:
        raise KeyError('There was an error finding %s in the metadata'
                       % time_column) from e
    return num_timepoints


def _check_time_column_numeric(metadata, time_column):
    try:
        metadata[time_column].dropna().astype(int)
    except Exception as e:
        raise ValueError('%s must be numeric' % time_column) from e


def _check_reference_column(metadata, reference_column):
    try:
        reference_series = metadata[reference_column]
    except Exception as e:
        raise KeyError('There was an error finding %s in the metadata'
                       % reference_column) from e
    return reference_series


def _check_associated_reference(reference_series, metadata, time_column,
                                filter_missing_references, reference_column):

    used_references = reference_series[~metadata[time_column].isna()]

    if used_references.isna().any():
        if filter_missing_references:
            metadata = metadata.dropna(subset=[reference_column])
            used_references = used_references.dropna()
        else:
            nan_references = used_references.index[used_references.isna()]
            raise KeyError('Missing references for the associated sample data.'
                           ' Please make sure that all samples with a'
                           ' timepoint value have an associated reference.'
                           ' IDs where missing references were found:'
                           ' %s' % (tuple(nan_references),))
    return metadata, used_references


def _check_subject_column(metadata, subject_column):
    try:
        subject_series = metadata[subject_column]
    except Exception as e:
        raise KeyError('There was an error finding %s in the metadata'
                       % subject_column) from e
    return subject_series


def _check_duplicate_subject_timepoint(subject_series, metadata,
                                        subject_column, time_column):
    for subject in subject_series:
        subject_df = metadata[metadata[subject_column] == subject]
        if not subject_df[time_column].is_unique:
            timepoint_list = subject_df[time_column].to_list()
            raise ValueError('There is more than one occurrence of a subject'
                             ' in a timepoint. All subjects must occur only'
                             ' once per timepoint. Subject %s appears in '
                             ' timepoints: %s' % (subject, timepoint_list))


def _check_subjects_in_all_timepoints(subject_series, num_timepoints,
                                     drop_incomplete_subjects, metadata,
                                     subject_column, used_references):

    subject_occurrence_series = (subject_series.value_counts())
    if (subject_occurrence_series < num_timepoints).any():
        if drop_incomplete_subjects:
            subject_to_keep = (subject_occurrence_series[
                                subject_occurrence_series ==
                                num_timepoints].index)
            metadata = metadata[metadata[subject_column].isin(subject_to_keep)]
            used_references = used_references.filter(axis=0,
                                                     items=metadata.index)
        else:
            incomplete_subjects = (subject_occurrence_series[
                                    subject_occurrence_series
                                    != num_timepoints].index).to_list()
            raise ValueError('Missing timepoints for associated subjects.'
                             ' Please make sure that all subjects have all'
                             ' timepoints or use drop_incomplete_subjects'
                             ' parameter. The incomplete subjects were %s'
                             % incomplete_subjects)
    return metadata, used_references


# PEDS calculation methods
def _create_recipient_table(reference_series, time_metadata, tabledf):
    subset_reference_series = \
        reference_series[reference_series.index.isin(time_metadata.index)]
    recipdf = tabledf[tabledf.index.isin(subset_reference_series.index)]
    return recipdf


def _create_masking(time_metadata, donordf, recipdf, reference_column):
    donor_index_masking = []
    for sample in recipdf.index:
        donor = time_metadata.loc[sample, reference_column]
        donor_index_masking.append(donordf.index.get_loc(donor))
    donordf = donordf.to_numpy()
    donormask = donordf[donor_index_masking]
    donormask = donormask.astype(int)
    recipdf = recipdf.to_numpy()
    recipdf = recipdf.astype(int)
    return donormask


def _mask_recipient(donormask, recipdf):
    maskedrecip = donormask & recipdf
    return maskedrecip


# Decommissioned Methods
def _compute_feature_peds(reference_series: pd.Series, table: pd.Series,
                          metadata: qiime2.Metadata) -> (pd.DataFrame):
    table = table > 0
    pedsdf = pd.DataFrame(columns=["FeatureID", "Timepoint", "Numerator",
                                   "Denominator", "PEDS"])
    for time, time_metadata in metadata.groupby("TPO"):
        donordf = table[table.index.isin(reference_series)]
        recipdf = _create_recipient_table(time_metadata, table)
        donormask = _create_masking(time_metadata, donordf, recipdf)
        maskedrecip = _mask_recipient(donormask, recipdf)
        num_sum = np.sum(maskedrecip, axis=0)
        donor_sum = np.sum(donormask, axis=0)

        for count, feature in enumerate(recipdf.columns):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                peds = num_sum[count] / donor_sum[count]
            pedsdf.loc[len(pedsdf)] = [feature, time, num_sum[count],
                                       donor_sum[count], peds]
    return pedsdf


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


def _compute_sample_peds(reference_series: pd.Series, table: pd.Series,
                         metadata: qiime2.Metadata, time_column: str,
                         reference_column: str,
                         subject_column: str) -> (pd.DataFrame):

    peds_series_list = []
    for sample in reference_series.index:
        sample_row = metadata.loc[sample]
        donor = sample_row[reference_column]
        subject = sample_row[subject_column]
        timepoint = sample_row[time_column]

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
