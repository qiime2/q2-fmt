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
from scipy.stats import mannwhitneyu
from collections import Counter
import os
import pkg_resources
import jinja2
import json

from q2_fmt._util import json_replace


def peds(ctx, table, metadata, peds_metric, time_column, reference_column,
         subject_column, filter_missing_references=False,
         drop_incomplete_subjects=False, drop_incomplete_timepoint=None,
         level_delimiter=None):

    peds_heatmap = ctx.get_action('fmt', 'peds_heatmap')

    results = []

    if peds_metric == 'sample':
        sample_peds = ctx.get_action('fmt', 'sample_peds')
        peds_dist = sample_peds(
            table=table, metadata=metadata, time_column=time_column,
            subject_column=subject_column, reference_column=reference_column,
            drop_incomplete_subjects=drop_incomplete_subjects,
            drop_incomplete_timepoint=drop_incomplete_timepoint,
            filter_missing_references=filter_missing_references)

    else:
        if drop_incomplete_subjects or drop_incomplete_timepoint:
            warnings.warn('Feature PEDS was selected as the PEDS metric, which'
                          ' does not accept `drop_incomplete_subjects` or'
                          ' `drop_incomplete_timepoint` as parameters. One'
                          ' (or both) of these parameters were detected in'
                          ' your input, and will be ignored.')

        feature_peds = ctx.get_action('fmt', 'feature_peds')
        peds_dist = feature_peds(
            table=table, metadata=metadata, time_column=time_column,
            subject_column=subject_column, reference_column=reference_column,
            filter_missing_references=filter_missing_references)
    results += peds_heatmap(data=peds_dist[0], level_delimiter=level_delimiter)

    return tuple(results)


def peds_heatmap(output_dir: str, data: pd.DataFrame,
                 level_delimiter: str = None):
    _rename_features(data=data, level_delimiter=level_delimiter)

    J_ENV = jinja2.Environment(
        loader=jinja2.PackageLoader('q2_fmt', 'assets')
    )

    x_label = "group"
    y_label = "subject"
    gradient = "measure"

    x_label_name = data[x_label].attrs['title']
    y_label_name = data[y_label].attrs['title']
    measure_name = data[gradient].attrs['title']
    title = f'{measure_name} of {y_label_name} across {x_label_name}'

    index = J_ENV.get_template('index.html')
    data = json.loads(data.to_json(orient='records'))
    spec_fp = pkg_resources.resource_filename(
        'q2_fmt', os.path.join('assets', 'spec.json')
    )
    with open(spec_fp) as fh:
        json_obj = json.load(fh)

    order = {"order": "ascending"}

    full_spec = json_replace(json_obj, data=data, x_label=x_label,
                             x_label_name=x_label_name,
                             y_label=y_label, y_label_name=y_label_name,
                             title=title, measure=gradient,
                             measure_name=measure_name, order=order)

    with open(os.path.join(output_dir, "index.html"), "w") as fh:
        spec_string = json.dumps(full_spec)
        fh.write(index.render(spec=spec_string))


def sample_peds(table: pd.DataFrame, metadata: qiime2.Metadata,
                time_column: str, reference_column: str, subject_column: str,
                filter_missing_references: bool = False,
                drop_incomplete_subjects: bool = False,
                drop_incomplete_timepoint: list = None) -> (pd.DataFrame):
    ids_with_data = table.index
    metadata = metadata.filter_ids(ids_to_keep=ids_with_data)
    column_properties = metadata.columns
    # TODO: Make incomplete samples possible move this to heatmap
    metadata = metadata.to_dataframe()
    num_timepoints = _check_for_time_column(metadata, time_column)
    _check_column_type(column_properties, "time",
                       time_column, "numeric")
    reference_series = _check_reference_column(metadata, reference_column)
    _check_column_type(column_properties, "reference",
                       reference_column, "categorical")
    # return things that should be removed
    metadata, used_references = \
        _filter_associated_reference(reference_series, metadata, time_column,
                                     filter_missing_references,
                                     reference_column)
    subject_series = _check_subject_column(metadata, subject_column)
    _check_column_type(column_properties, "subject",
                       subject_column, "categorical")
    _check_duplicate_subject_timepoint(subject_series, metadata,
                                       subject_column, time_column)
    if drop_incomplete_timepoint is not None:
        metadata = _drop_incomplete_timepoints(metadata, time_column,
                                               drop_incomplete_timepoint)
        table.filter(items=metadata.index)
    # return things that should be removed
    metadata, used_references = \
        _check_subjects_in_all_timepoints(subject_series, num_timepoints,
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
    column_properties = metadata.columns
    metadata = metadata.to_dataframe()

    _ = _check_for_time_column(metadata, time_column)
    _check_column_type(column_properties, "time",
                       time_column, "numeric")
    reference_series = _check_reference_column(metadata, reference_column)
    _check_column_type(column_properties, "reference",
                       reference_column, "categorical")
    metadata, used_references = \
        _filter_associated_reference(reference_series, metadata, time_column,
                                     filter_missing_references,
                                     reference_column)
    _ = _check_subject_column(metadata, subject_column)
    _check_column_type(column_properties, "subject",
                       subject_column, "categorical")
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
    reference_overlap = reference_series.isin(table.index)
    try:
        assert all(reference_overlap)
    except AssertionError as e:
        missing_ref = reference_series[~reference_overlap].unique()
        raise AssertionError('Reference IDs: %s provided were not found in'
                             ' the feature table. Please confirm that all'
                             ' values in reference column are present in the'
                             ' feature table' % missing_ref) from e
    donor_df = table[table.index.isin(reference_series)]
    recip_df = _create_recipient_table(reference_series, metadata, table)

    donormask = _create_masking(time_metadata=metadata, donor_df=donor_df,
                                recip_df=recip_df,
                                reference_column=reference_column)
    maskedrecip = donormask & recip_df
    if peds_type == "Sample":
        num_sum = np.sum(maskedrecip, axis=1)
        donor_sum = np.sum(donormask, axis=1)
        for count, sample in enumerate(recip_df.index):
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
        for count, feature in enumerate(recip_df.columns):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                peds = num_sum[count] / donor_sum[count]
            peds_df.loc[len(peds_df)] = [feature, peds, num_sum[count],
                                         donor_sum[count], peds_time, feature]
            peds_df = peds_df.dropna()

        peds_df['id'].attrs.update({
            'title': "Feature ID",
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


# prep method
def _rename_features(level_delimiter, data: pd.DataFrame):
    if ("recipients with feature" in data.columns and
            level_delimiter is not None):
        group_name = data["group"].attrs['title']
        subject_name = data["subject"].attrs['title']
        measure_name = data["measure"].attrs['title']

        y_labels = []
        seen = Counter()
        subject_seen = []
        for i, sub in enumerate(data['subject']):
            if level_delimiter in sub:
                fields = [field for field in sub.split(level_delimiter)
                          if not field.endswith('__')]
            else:
                # This is necessary to handle a case where the delimiter
                # isn't found but the sub ends with __. In that case, sub would
                # be completely thrown out.
                fields = [sub]
            subject_seen.append(sub)
            most_specific = fields[-1]
            if most_specific in seen and sub not in subject_seen:
                y_labels.append(f"{seen[most_specific]}: {most_specific} *")
            else:
                y_labels.append(most_specific)
            seen[most_specific] += 1
        data['subject'] = y_labels

        data['id'] = [i.replace(level_delimiter, ' ') for i in data['id']]

        # currently attrs get deleted with df is changed. right now the best
        # way to solve this is by saving them as temp and saving them at the
        # end

        data['subject'].attrs.update({'title': subject_name,
                                      'description': ''})
        data['group'].attrs.update({'title': group_name,
                                    'description': ''})
        data['measure'].attrs.update({'title': measure_name,
                                      'description': ''})


# Filtering methods
def _check_for_time_column(metadata, time_column):
    try:
        num_timepoints = metadata[time_column].dropna().unique().size
    except Exception as e:
        _check_column_missing(metadata, time_column, "time", e)
    return num_timepoints


def _check_reference_column(metadata, reference_column):
    try:
        reference_series = metadata[reference_column]
    except Exception as e:
        _check_column_missing(metadata, reference_column, "reference", e)
    return reference_series


def _filter_associated_reference(reference_series, metadata, time_column,
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
        _check_column_missing(metadata, subject_column, "subject", e)
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


def _drop_incomplete_timepoints(metadata, time_column,
                                drop_incomplete_timepoint):
    for time in drop_incomplete_timepoint:
        try:
            assert (float(time)
                    in metadata[time_column].unique())
        except AssertionError as e:
            raise AssertionError('The provided incomplete timepoint `%s` was'
                                 ' not found in the metadata. Please check'
                                 ' that the incomplete timepoint provided is'
                                 ' in your provided --p-time-column: `%s`'
                                 % (time, time_column)) from e
        metadata = metadata[metadata[time_column] !=
                            float(time)]
    return metadata


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
                             ' timepoints. You can drop these subjects by'
                             ' using the drop_incomplete_subjects parameter or'
                             ' drop any timepoints that have large numbers'
                             ' of subjects missing by using the'
                             ' drop_incomplete_timepoints parameter. The'
                             ' incomplete subjects were %s'
                             % incomplete_subjects)
    return metadata, used_references


def _check_column_type(column_properties, parameter_type, column, column_type):
    try:
        assert column_properties[column].type == column_type
    except AssertionError as e:
        raise AssertionError('Non-%s values found in `--p-%s-column`.'
                             ' Please make sure the column selected contains'
                             ' the correct MetadataColumn type. Column with'
                             ' non-%s values that was'
                             ' selected: `%s`' % (column_type, parameter_type,
                                                  column_type, column)) from e


def _check_column_missing(metadata, column, column_param, e):
    if column == metadata.index.name:
        raise KeyError('The `--p-%s-column` input provided was the'
                       ' same as the index of the metadata.'
                       ' `--p-%s-column` can not be the same as the'
                       ' index of metadata:'
                       ' `%s`' % (column_param, column_param, column)) from e
    else:
        raise KeyError('There was an error finding the provided'
                       ' `--p-%s-column`: `%s` in the metadata'
                       % (column_param, column)) from e


# PEDS calculation methods
def _create_recipient_table(reference_series, time_metadata, table_df):
    subset_reference_series = \
        reference_series[reference_series.index.isin(time_metadata.index)]
    recip_df = table_df[table_df.index.isin(subset_reference_series.index)]
    return recip_df


def _create_masking(time_metadata, donor_df, recip_df, reference_column):
    donor_index_masking = []
    for sample in recip_df.index:
        donor = time_metadata.loc[sample, reference_column]
        donor_index_masking.append(donor_df.index.get_loc(donor))
    donor_df = donor_df.to_numpy()
    donor_mask = donor_df[donor_index_masking]
    donor_mask = donor_mask.astype(int)
    recipdf = recip_df.to_numpy()
    recipdf = recipdf.astype(int)
    return donor_mask


def _mask_recipient(donor_mask, recip_df):
    maskedrecip = donor_mask & recip_df
    return maskedrecip


def peds_bootstrap(table: pd.DataFrame, metadata: qiime2.Metadata,
                   time_column: str, reference_column: str,
                   subject_column: str,
                   filter_missing_references: bool = False,
                   drop_incomplete_subjects: bool = False,
                   bootstrap_replicates: int = 999):
    metadata_df = metadata.to_dataframe
    # TODO: Grab Donor in a more logic way
    donor = metadata_df[metadata.index.isin(
        _check_reference_column(metadata=metadata_df,
                                reference_column=reference_column))]
    recipient = metadata_df.loc[metadata_df[reference_column].notnull()]
    fake_donor = []
    for i in range(0, bootstrap_replicates+1):
        if i == 0:
            peds, = sample_peds(table=table, metadata=metadata,
                                time_column=time_column,
                                reference_column=reference_column,
                                subject_column=subject_column)
            peds = peds.view(pd.DataFrame).set_index("id")
            real_temp = peds["measure"].to_list()
        else:
            shifted_list = recipient[reference_column].sample(frac=1).to_list()
            recipient.loc[:, reference_column] = shifted_list
            metadata_df = pd.concat([donor, recipient])
            metadata = qiime2.Metadata(metadata_df)
            peds, = sample_peds(table=table, metadata=metadata,
                                time_column=time_column,
                                reference_column=reference_column,
                                subject_column=subject_column)
            peds = peds.view(pd.DataFrame).set_index("id")
            fake_donor = fake_donor + peds["measure"].to_list()
    s, p = mannwhitneyu(real_temp, fake_donor, alternative='greater')
