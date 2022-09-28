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
    except Exception:
        raise KeyError('There was and error finding %s in the metadata'
                       % time_column)
    try:
        metadata[time_column].dropna().astype(int)
    except Exception:
        raise KeyError('%s must be numeric' % time_column)
    try:
        subject_occcurance_df = (metadata[subject_column].value_counts()
                                 .to_frame())
    except Exception:
        raise KeyError('There was and error finding %s in the metadata'
                       % subject_column)
    if (subject_occcurance_df[subject_column] != num_timepoints).any():
        if drop_incomplete_subjects:
            subject_to_keep = (subject_occcurance_df
                               .loc[subject_occcurance_df[subject_column]
                                    == num_timepoints].index)
            metadata = metadata[metadata[subject_column].isin(subject_to_keep)]
        else:
            incomplete_subjects = (subject_occcurance_df
                                   .loc[subject_occcurance_df[subject_column]
                                        < num_timepoints].index).to_list()
            raise KeyError('Missing timepoints for associated subjects.'
                           ' Please make sure that all subjects have all'
                           ' timepoints or use drop_incomplete_subjects'
                           ' parameter. The incomplete subjects were %s'
                           % incomplete_subjects)

    try:
        reference_series = metadata[reference_column]
    except Exception:
        raise KeyError('There was and error finding %s in the metadata'
                       % reference_column)

    if reference_series.isna().any():
        if filter_missing_references:
            reference_series = reference_series.dropna()
        else:
            nan_references = reference_series.index[reference_series.isna()]
            raise KeyError('Missing references for the associated sample data.'
                           ' Please make sure that all samples with a'
                           ' timepoint value have an associated reference.'
                           ' IDs where missing references were found:'
                           ' %s' % (tuple(nan_references),))

    peds_df = _compute_peds(reference_series, table, metadata, time_column,
                            reference_column, subject_column)
    return peds_df


def _compute_peds(reference_series: pd.Series, table: pd.Series,
                  metadata: qiime2.Metadata, time_column: str,
                  reference_column: str,
                  subject_column: str) -> (pd.DataFrame):

    PEDSserieslist = []
    for sample in reference_series.index:
        donor = metadata.loc[sample][reference_column]
        subject = metadata.loc[sample][subject_column]
        timepoint = metadata.loc[sample][time_column]

        donorPresentList = _get_observed_features(table, donor)
        donorNumPresent = _count_observed_features(donorPresentList)

        samplePresentList = _get_observed_features(table, sample)

        intersect = (donorPresentList & samplePresentList).sum()
        intersectNumPresent = _count_observed_features(intersect)

        peds = (intersectNumPresent/donorNumPresent)
        (PEDSserieslist.append((sample, peds, intersectNumPresent,
                                donorNumPresent, donor, subject, timepoint)))
    PEDSdf = pd.DataFrame(PEDSserieslist,
                          columns=['id', 'measure',
                                   'Transfered_Donor_Features',
                                   'Total_Donor_Features', 'donor', 'subject',
                                   'group'])

    # use title for correcting ugly names
    PEDSdf['id'].attrs.update({
        'title': reference_series.index.name,
        'description': 'Sample IDS'
    })
    PEDSdf['measure'].attrs.update({
        'title': "PEDS",
        'description': 'Proportional Engraftment of Donor Strains '
    })
    PEDSdf['group'].attrs.update({
        'title': time_column,
        'description': 'Time'
    })
    PEDSdf["subject"].attrs.update({
        'title': subject_column,
        'description': 'ID to link samples across time'
    })
    PEDSdf["Transfered_Donor_Features"].attrs.update({
        'title': "Transfered Donor Features",
        'description': '...'
    })
    PEDSdf['Total_Donor_Features'].attrs.update({
        'title': "Total Donor Features",
        'description': '...'
    })
    PEDSdf['donor'].attrs.update({
        'title': reference_column,
        'description': 'donor'
    })
    return PEDSdf


def _get_observed_features(table, id):
    try:
        present = table.loc[id] > 0
    except Exception:
        raise KeyError('There was an error finding the sample %s in'
                       ' your feature table' % id)
    return present


def _count_observed_features(presentList):
    numPresent = presentList.sum()
    return numPresent
