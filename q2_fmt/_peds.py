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
import random
import itertools
import warnings
from scipy.stats import false_discovery_control, combine_pvalues
from collections import Counter
import os
import pkg_resources
import jinja2
import json

from qiime2 import Metadata
from q2_fmt._engraftment import _get_to_baseline_ref
from q2_fmt._util import json_replace
from q2_stats._visualizer import _make_stats


def peds(ctx, table, metadata, peds_metric, time_column, reference_column,
         subject_column, filter_missing_references=False,
         drop_incomplete_subjects=False, drop_incomplete_timepoints=None,
         level_delimiter=None):

    heatmap = ctx.get_action('fmt', 'heatmap')

    results = []

    if peds_metric == 'sample':
        sample_peds = ctx.get_action('fmt', 'sample_peds')
        peds_dist = sample_peds(
            table=table, metadata=metadata, time_column=time_column,
            subject_column=subject_column, reference_column=reference_column,
            drop_incomplete_subjects=drop_incomplete_subjects,
            drop_incomplete_timepoints=drop_incomplete_timepoints,
            filter_missing_references=filter_missing_references)

    else:
        if drop_incomplete_subjects or drop_incomplete_timepoints:
            warnings.warn('Feature PEDS was selected as the PEDS metric, which'
                          ' does not accept `drop_incomplete_subjects` or'
                          ' `drop_incomplete_timepoints` as parameters. One'
                          ' (or both) of these parameters were detected in'
                          ' your input, and will be ignored.')

        feature_peds = ctx.get_action('fmt', 'feature_peds')
        peds_dist = feature_peds(
            table=table, metadata=metadata, time_column=time_column,
            subject_column=subject_column, reference_column=reference_column,
            filter_missing_references=filter_missing_references)
    results += heatmap(data=peds_dist[0],
                       level_delimiter=level_delimiter)

    return tuple(results)


def heatmap(output_dir: str, data: pd.DataFrame,
            level_delimiter: str = None,
            per_subject_stats: pd.DataFrame = None,
            global_stats: pd.DataFrame = None):
    try:
        assert "baseline" not in data.columns or (global_stats is None and
                                                  per_subject_stats is None)
    except AssertionError as e:
        raise AssertionError("The input data provided was created with"
                             " `fmt sample_pprs`. This is not compatible with"
                             " statistics created from `fmt peds-simulation`"
                             " because they are created from separate"
                             " references (i.e. baseline and donor)") from e

    _rename_features(data=data, level_delimiter=level_delimiter)
    gstats = None
    table1 = None
    psstats = None
    if global_stats is not None:
        gstats = global_stats.to_html(index=False)
    if per_subject_stats is not None:
        table1, psstats = _make_stats(per_subject_stats)
    J_ENV = jinja2.Environment(
        loader=jinja2.PackageLoader('q2_fmt', 'assets')
    )

    x_label = "group"
    y_label = "subject"
    gradient = "measure"
    if "all possible recipients with feature" in data.columns:
        n_label = "all possible recipients with feature"
    elif "total_donor_features" in data.columns:
        n_label = "total_donor_features"
    elif "total_baseline_features" in data.columns:
        n_label = "total_baseline_features"
    data_denom = "datum['%s']" % n_label

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
                             measure_name=measure_name, order=order,
                             n_label=n_label, data_denom=data_denom)

    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        spec_string = json.dumps(full_spec)
        fh.write(index.render(spec=spec_string,
                              persubjectstats=psstats,
                              globalstats=gstats,
                              table1=table1))


def sample_peds(table: pd.DataFrame, metadata: qiime2.Metadata,
                time_column: str, reference_column: str, subject_column: str,
                filter_missing_references: bool = False,
                drop_incomplete_subjects: bool = False,
                drop_incomplete_timepoints: list = None) -> (pd.DataFrame):

    ids_with_data = table.index
    metadata = metadata.filter_ids(ids_to_keep=ids_with_data)
    column_properties = metadata.columns
    metadata = metadata.to_dataframe()
    if drop_incomplete_timepoints is not None:
        metadata = _drop_incomplete_timepoints(metadata, time_column,
                                               drop_incomplete_timepoints)
        table.filter(items=metadata.index)
    # TODO: Make incomplete samples possible move this to heatmap
    num_timepoints, time_col = _check_for_time_column(metadata, time_column)
    _check_column_type(column_properties, "time",
                       time_column, "numeric")
    metadata = metadata.filter(items=time_col.index, axis=0)
    reference_series = _check_reference_column(metadata, reference_column)
    _check_column_type(column_properties, "reference",
                       reference_column, "categorical")
    # return things that should be removed
    metadata, used_references = \
        _filter_associated_reference(reference_series, metadata, time_column,
                                     filter_missing_references,
                                     reference_column, ids_with_data)
    subject_series = _check_subject_column(metadata, subject_column)
    _check_column_type(column_properties, "subject",
                       subject_column, "categorical")
    _check_duplicate_subject_timepoint(subject_series, metadata,
                                       subject_column, time_column)

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

    _, _ = _check_for_time_column(metadata, time_column)
    _check_column_type(column_properties, "time",
                       time_column, "numeric")
    reference_series = _check_reference_column(metadata, reference_column)
    _check_column_type(column_properties, "reference",
                       reference_column, "categorical")
    metadata, used_references = \
        _filter_associated_reference(reference_series, metadata, time_column,
                                     filter_missing_references,
                                     reference_column, ids_with_data)
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
                  reference_column: str = None) -> (pd.DataFrame):
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
    if peds_type == "Sample" or peds_type == "PPRS":
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
        if peds_type == "PPRS":
            transfered = "transfered_baseline_features"
            total = 'total_baseline_features'
            ref = 'baseline'
            measure_description = ('Proportional Persistence of Recipient'
                                   ' Strains')
        elif peds_type == "Sample":
            transfered = "transfered_donor_features"
            total = 'total_donor_features'
            ref = 'donor'
            measure_description = 'Proportional Engraftment of Donor Strains'

        peds_df['id'].attrs.update({
            'title': reference_series.index.name,
            'description': 'Sample IDs'
        })
        peds_df['measure'].attrs.update({
            'title': peds_type,
            'description': measure_description
        })
        peds_df['group'].attrs.update({
            'title': time_column,
            'description': 'Time'
        })
        peds_df["subject"].attrs.update({
            'title': subject_column,
            'description': 'Subject IDs linking samples across time'
        })
        peds_df[transfered].attrs.update({
            'title': "Transfered Reference Features",
            'description': '...'
        })
        peds_df[total].attrs.update({
            'title': "Total Reference Features",
            'description': '...'
        })
        peds_df[ref].attrs.update({
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
        time_col = metadata[time_column].dropna()
        num_timepoints = time_col.unique().size
    except Exception as e:
        _check_column_missing(metadata, time_column, "time", e)
    return num_timepoints, time_col


def _check_reference_column(metadata, reference_column):
    try:
        reference_series = metadata[reference_column]
    except Exception as e:
        _check_column_missing(metadata, reference_column, "reference", e)
    return reference_series


def _filter_associated_reference(reference_series, metadata, time_column,
                                 filter_missing_references, reference_column,
                                 ids_with_data):
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
    available_references = (used_references.isin(ids_with_data))
    if not available_references.all():
        if filter_missing_references:
            used_references = used_references[available_references]
        else:
            raise KeyError('References included in the metadata are missing'
                           ' from the feature table. Please make sure all'
                           ' references included in the metadata are also'
                           ' present in the table.'
                           ' Missing references: %s'
                           % list(used_references[~available_references]
                                  .unique()))

    used_references = used_references[available_references]
    metadata = metadata.filter(items=used_references.index, axis=0)
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
                                drop_incomplete_timepoints):
    for time in drop_incomplete_timepoints:
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
    return donor_mask


def _mask_recipient(donor_mask, recip_df):
    maskedrecip = donor_mask & recip_df
    return maskedrecip


def sample_pprs(table: pd.DataFrame, metadata: qiime2.Metadata,
                time_column: str, baseline_timepoint: str, subject_column: str,
                filter_missing_references: bool = False,
                drop_incomplete_subjects: bool = False,
                drop_incomplete_timepoints: list = None) -> (pd.DataFrame):
    ids_with_data = table.index
    metadata = metadata.filter_ids(ids_to_keep=ids_with_data)
    column_properties = metadata.columns
    # TODO: Make incomplete samples possible move this to heatmap
    metadata = metadata.to_dataframe()
    if drop_incomplete_timepoints:
        metadata = _drop_incomplete_timepoints(metadata, time_column,
                                               drop_incomplete_timepoints)
        table.filter(items=metadata.index)
    num_timepoints, time_col = _check_for_time_column(metadata, time_column)
    _check_column_type(column_properties, 'time',
                       time_column, 'numeric')
    metadata = metadata.filter(items=time_col.index, axis=0)

    used_references =\
        _get_to_baseline_ref(time_col=metadata[time_column],
                             baseline_timepoint=baseline_timepoint,
                             time_column=time_column,
                             subject_column=subject_column,
                             metadata=Metadata(metadata))

    subject_series = _check_subject_column(metadata, subject_column)
    _check_column_type(column_properties, 'subject',
                       subject_column, 'categorical')
    _check_duplicate_subject_timepoint(subject_series, metadata,
                                       subject_column, time_column)
    # return things that should be removed
    metadata, used_references = \
        _check_subjects_in_all_timepoints(subject_series, num_timepoints,
                                          drop_incomplete_subjects, metadata,
                                          subject_column, used_references)

    peds_df = pd.DataFrame(columns=['id', 'measure',
                                    'transfered_baseline_features',
                                    'total_baseline_features', 'baseline',
                                    'subject', 'group'])
    baseline_metadata = metadata.join(used_references)
    peds_df = _compute_peds(peds_df=peds_df, peds_type='PPRS',
                            peds_time=np.nan, reference_series=used_references,
                            table=table, metadata=baseline_metadata,
                            time_column=time_column,
                            subject_column=subject_column,
                            reference_column=used_references.name)
    return peds_df


def peds_simulation(table: pd.DataFrame, metadata: qiime2.Metadata,
                    time_column: str, reference_column: str,
                    subject_column: str,
                    filter_missing_references: bool = False,
                    drop_incomplete_subjects: bool = False,
                    drop_incomplete_timepoints: list = None,
                    num_iterations: int = 999) -> (pd.DataFrame, pd.DataFrame):

    ids_with_data = table.index
    metadata = metadata.filter_ids(ids_to_keep=ids_with_data)

    metadata_df = metadata.to_dataframe()
    reference_series = _check_reference_column(metadata_df,
                                               reference_column)

    (metadata_df,
     used_references) = _filter_associated_reference(reference_series,
                                                     metadata_df, time_column,
                                                     filter_missing_references,
                                                     reference_column,
                                                     ids_with_data)

    if len(used_references.unique()) == 1:
        raise AssertionError("There is only one donated microbiome in your"
                             " data. A Monte Carlo simulation shuffles"
                             " your donated microbiome and recipient pairing"
                             " and needs more than one donated microbiome"
                             " to successfully shuffle.")
    recipient = metadata_df.loc[metadata_df[reference_column].notnull()]

    if len(recipient.index.unique()) == 1:
        raise AssertionError("There is only one recipient in your"
                             " data. A Monte Carlo simulation shuffles"
                             " your donated microbiome and recipient pairing"
                             " and needs more than one recipient"
                             " to successfully shuffle.")

    peds = sample_peds(
           table=table, metadata=metadata,
           time_column=time_column,
           reference_column=reference_column,
           subject_column=subject_column,
           filter_missing_references=filter_missing_references,
           drop_incomplete_subjects=drop_incomplete_subjects,
           drop_incomplete_timepoints=drop_incomplete_timepoints
           ).set_index("id")
    actual_peds = peds["measure"]

    # Mismatch simulation:
    table = table > 0
    recip_df = _create_recipient_table(used_references, metadata_df, table)
    donor_df = table[table.index.isin(used_references)]
    mismatched_df = \
        _create_mismatched_pairs(recip_df,
                                 metadata_df,
                                 used_references,
                                 reference_column)
    duplicated_recip_table = _create_duplicated_recip_table(mismatched_df,
                                                            recip_df)
    donor_mask = _create_sim_masking(mismatched_df, donor_df, reference_column)
    recip_mask = _mask_recipient(donor_mask, duplicated_recip_table)
    # Numerator for PEDS Calc. (Number of Donor features in the Recipient)
    num_engrafted_donor_features = np.sum(recip_mask.values, axis=1)
    # Denominator for PEDS Calc. (Number of unique features in the Donor)
    num_donor_features = np.sum(donor_mask, axis=1)
    # This ignores warnings that come from dividing by 0.
    # mismatched_peds will be Nan if the denominator is 0 and thats reasonable.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mismatched_peds = num_engrafted_donor_features/num_donor_features
    per_sub_stats = _per_subject_stats(mismatched_peds,
                                       actual_peds, num_iterations)
    global_stats = _global_stats(per_sub_stats['p-value'])
    return per_sub_stats, global_stats


# peds_simulation helper functions

def _create_mismatched_pairs(recip_df, metadata, used_references,
                             reference_column):
    """Creates a Dataframe of all the incorrect Donor-Recipient pairs

    Creates a list of tuples of all donor-recipient pairs then filters out
    any real donor-recipient pairs. The result is a Dataframe that includes
    only mismatched donor-recipient pairs.

    Parameters
    ----------
    recip_df: pd.DataFrame
        A feature table of FMT recipients.
    metadata: pd.DataFrame
        The sample metadata for the study.
    used_references: pd.Series
        A series with recipients as the index and the associated reference
        as the values.
    reference_column: Str
        Name of the reference column in the Sample Metadata.

    Returns
    -------
    mismatched_df: pd.DataFrame
        A DataFrame containing all mismatched pairs of donors and recipients.

    Examples
    --------
    >>> recip_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3'],
            'Feature1': [1, 0, 0],
            'Feature2': [0, 1, 0],
            'Feature3': [0, 0, 1]}).set_index('id')

    >>> metadata_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'donor1', 'donor2', 'donor3'],
            'Ref': ['donor1', 'donor2', 'donor3',
                    np.nan, np.nan, np.nan],
            'subject': ['sub1', 'sub2', 'sub3',
                        np.nan, np.nan, np.nan],
            'group': [1, 1, 1,
                      np.nan, np.nan, np.nan],
            'Location': [np.nan, np.nan,np.nan,
                         'test', 'test','test']}).set_index('id')

    >>> used_references = pd.Series(data=['donor1', 'donor2', 'donor3'],
                                    index=pd.Index(['sample1', 'sample2',
                                                    'sample3'], name='id'),
                                    name='Ref')

    >>> reference_column = "Ref"

    >>> _create_mismatched_pairs(recip_df, metadata_df, used_references,
                                 reference_column)
    pd.DataFrame({'id': ["sample1", "sample1", "sample2", "sample2",
                         "sample3", "sample3"],
                  "Ref": ["donor2", "donor3", "donor1", "donor3",
                          "donor1", "donor2"]}).set_index('id')
    """
    matched_pairs = list(zip(used_references.index, used_references))
    donors = metadata[reference_column].dropna().unique()
    # Generates all donor recipient pairs and then removes
    # matched donor and recipient pairs
    filtered = \
        [pair
         for pair in itertools.product(recip_df.index,
                                       donors)if pair not in matched_pairs]
    idx, values = zip(*filtered)
    mismatched_df = pd.DataFrame({"id": idx,
                                  reference_column: values}).set_index("id")
    return mismatched_df


def _create_duplicated_recip_table(mismatched_df, recip_df):
    """Creates a recipient feature table that is the same dimensions as
    the mismatched_df

    Creates a recipient table that duplicates feature information so that the
    recip_df is the same dimension as the mismatched_df. This makes numpy array
    math easier later on.

    Parameters
    ----------
    mismatched_df: pd.DataFrame
        A Dataframe that contains recipient samples as the index and mismatched
        donors as the values. A recipient sample will appear as many times as
        there are mismatched donors to pair with.
    recip_df: pd.DataFrame
        A feature table of FMT recipients.

    Returns
    -------
    duplicated_table: pd.DataFrame
        A recipent feature table where feature information is duplicated so
        that recipient is in the same diminsion as mismatched_df

    Examples
    --------
    >>> mismatched_df = pd.DataFrame({
            'id': ["sample1", "sample1", "sample2", "sample2",
                   "sample3", "sample3"],
            "Ref": ["donor2", "donor3", "donor1", "donor3",
                    "donor1", "donor2"]}).set_index('id')

    >>> recip_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3'],
            'Feature1': [1, 0, 0],
            'Feature2': [0, 1, 0],
            'Feature3': [0, 0, 1]}).set_index('id')

    >>> _create_duplicated_recip_table(mismatched_df, recip_df)

    pd.DataFrame({
            'id': ["sample1", "sample1", "sample2", "sample2",
                   "sample3", "sample3"],
            'Feature1': [1, 1, 0, 0, 0, 0],
            'Feature2': [0, 0, 1, 1, 0, 0],
            'Feature3': [0, 0, 0, 0, 1, 1]}).set_index('id')
    """
    duplicated_table = \
        mismatched_df.merge(recip_df, left_index=True, right_index=True,
                            how='left').drop(mismatched_df.columns, axis=1)
    return duplicated_table


def _create_sim_masking(mismatched_df, donor_df, reference_column):
    """Create a donor mask to mask recipient features that aren't in the donor.

    Creates a Donor Numpy array that duplicates donor samples to match up
    with duplicated_table. This will mask recipient feature that aren't
    in the donor.

    Parameters
    ----------
    mismatched_df: pd.DataFrame
        A Dataframe that contains recipient samples as the index and mismatched
        donors as the values. A recipient sample will appear as many times as
        there are mismatched donors to pair with.
    donor_df: pd.DataFrame
        A feature table of FMT donors.
    reference_column: str
        Name of the reference column in the Sample Metadata.

    Returns
    -------
    donor_mask: ndarray
        a numpy array of donor feature information with the same order as
        duplicated_recip_table. This allows for easy np array math.

    Examples
    --------
    >>> mismatched_df = pd.DataFrame({
                'id': ["sample1", "sample1", "sample2", "sample2",
                       "sample3", "sample3"],
                "Ref": ["donor2", "donor3", "donor1", "donor3",
                        "donor1", "donor2"]}
                                ).set_index('id')
    >>> donor_df = pd.DataFrame({
              'id': ['donor1', 'donor2', 'donor3'],
              'Feature1': [1, 0, 0],
              'Feature2': [0, 1, 0],
              'Feature3': [0, 0, 1]}).set_index('id')

    >>> _create_sim_masking(mismatched_df, donor_df, reference_column)

    ndarray[[0, 1, 0],
            [0, 0, 1],
            [1, 0, 0],
            [0, 0, 1],
            [1, 0, 0],
            [0, 1, 0]]
    """
    donors = mismatched_df[reference_column]
    donor_index_masking = []
    for donor in donors:
        donor_index_masking.append(donor_df.index.get_loc(donor))
    donor_df_np = donor_df.to_numpy()
    donor_mask = donor_df_np[donor_index_masking]
    return donor_mask


def _simulate_uniform_distro(mismatched_peds, k):
    """Randomly samples the mismatched PEDS values.

    Creates a uniform distribution of mismatched PEDS values by randomly
    sampling `num_iterations` times with replacement.

    Parameters
    ----------
    mismatched_peds: list
        A list that contains all mismatched PEDS values.
    k: int
        Number of iterations(`k`) to run simulations (Number of times to
        randomly sample mismatched_peds)

    Returns
    -------
    peds_iters: pd.Series
        a pd.Series with all num_iterations of mismatched PEDS Values. This
        will later be compared to an actual PEDS value.

    Examples
    --------
    >>> mismatched_peds = [0, 0, 0,0 ]
    >>> num_iterations = 10

    >>> _simulate_uniform_distro(mismatched_peds, num_iterations)

    pd.Series(data = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              index = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    """
    peds_iters = random.choices(mismatched_peds, k=k)
    # Tranforming to series for easy series math in _per_subject_stats()
    peds_iters = pd.Series(peds_iters)
    return peds_iters


def _peds_sim_stats(value, peds_iters, num_iterations):
    """Calculates test statistics, and p-value.

    Calculates tests statistics (`count_gte` and ` count_less`) and p-value
    for PEDS Monte Carlo Simulation

    Parameters
    ----------
    value: float
        A actual PEDS value to compare against
    peds_iters: pd.Series
        a pd.Series with all num_iterations of mismatched PEDS Values. This
        will be compared to the actual PEDS value.
    num_iterations: int
        Number of iterations to run simulations (Number of times to
        randomly sample mismatched_peds)

    Returns
    -------
    count_gte: int
        Count of mismatched PEDS values that were greater than the actual PEDS
        value
    count_less: int
        Count of mismatched PEDS values that were less than the actual PEDS
        value. This is calculated by ``num_interations - count_gte``
    per_subject_p: float
        The p-value associated with the above test stats.

    Examples
    --------
    >>> value = 1
    >>> peds_iters = pd.Series(data = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                               index = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    >>> num_iterations = 10

    >>> _peds_sim_stats()

    count_gte = 0
    count_less = 10
    per_subject_p = (1/11) (Note: the 10 iterations is not enough to get a
                            significant p-value)
    """
    gte_series = peds_iters >= value
    count_gte = gte_series.sum()
    count_less = num_iterations-count_gte
    # adding 1 here because you can mathmatically can get p-value of 0 from a
    # Monte Carlo Simulation
    per_subject_p = (count_gte + 1)/(num_iterations+1)
    return count_gte, count_less, per_subject_p


def _per_subject_stats(mismatched_peds, actual_peds,
                       num_iterations):
    """Creates per subject PEDS stats

    Creates per subject PEDS stats by sampling mismatch PEDS values 'iteration'
    number of times. The list of mismatched PEDS values is compared to the
    actual PEDS value.

    Parameters
    ----------
    mismatched_peds: list
        A list that contains all mismatched PEDS values.
    actual_peds: pd.Series
        A Series containing Sample IDs as the index and actual PEDS values for
        the sample as the value.
    num_iterations: int
        Number of num_iterations to run simulations (Number of times to
        randomly sample mismatched_peds)

    Returns
    -------
    per_subject_stats: pd.DataFrame
        a Dataframe containing test statistics, p-values, q-values and
        sample size info. Note: this is the DataFrame that is returned
        to the user as the per-subject-stats.qza

    Examples
    --------
    >>> mismatched_peds = [0, 0, 0, 0]
    >>> actual_peds = pd.Series([data = [1, 1, 1, 1],
                                 index = ["sample1", "sample2",
                                          "sample3", "sample4"]])
    >>> num_iterations = 10

    >>> _per_subject_stats(mismatched_peds, actual_peds,
                       num_iterations):

    pd.Dataframe({"A:group": ["sample1", "sample2", "sample3", "sample4"],
                   "A:n": [1, 1, 1, 1],
                   "A:measure" : [1, 1, 1, 1],
                   "B:group" : ["shuffled recipients", "shuffled recipients",
                                "shuffled recipients", "shuffled recipients"],
                   "B:n": [4, 4, 4, 4],
                   "B:measure": [0, 0, 0, 0],
                   "n": [10, 10, 10, 10],
                   "test-statistic": [10, 10, 10, 10],
                   "p-value": [0.001, 0.001, 0.001, 0.001],
                   "q-value": [0.004, 0.002, 0.00133, 0.001]})
    """
    peds_iters_means = []
    count_less_list = []
    per_subject_p_list = []
    for value in actual_peds:
        peds_iters = _simulate_uniform_distro(mismatched_peds, num_iterations)
        peds_iters_means.append(peds_iters.mean())
        _, count_less, per_subject_p = _peds_sim_stats(value, peds_iters,
                                                       num_iterations)
        count_less_list.append(count_less)
        per_subject_p_list.append(per_subject_p)

    per_subject_q = false_discovery_control(ps=per_subject_p_list, method='bh')
    per_sub_stats = pd.DataFrame({'A:group': actual_peds.index,
                                 'A:n': 1,
                                  'A:measure': actual_peds.values,
                                  'B:group': "shuffled recipients",
                                  'B:n': len(mismatched_peds),
                                  'B:measure': peds_iters_means,
                                  'n': num_iterations,
                                  'test-statistic': count_less_list,
                                  'p-value': per_subject_p_list,
                                  'q-value': per_subject_q})
    per_sub_stats['A:group'].attrs.update({'title': 'actual_values',
                                          'description': 'PEDS values'
                                           ' calculated with actual recipient'
                                           ' and donated microbiome pairing'})
    per_sub_stats['B:group'].attrs.update({'title': 'shuffled_values',
                                          'description': 'PEDS values'
                                           ' calculated with shuffled'
                                           ' recipient and donated microbiome'
                                           ' pairings'})
    n = {'title': 'count', 'description': 'Number of recipients and donated'
         ' microbiome pairings'}
    per_sub_stats['A:n'].attrs.update(n)
    per_sub_stats['B:n'].attrs.update(n)
    measure = {
        'title': 'Mean PEDS Value',
        'description': 'Mean PEDS Value'
    }
    per_sub_stats['A:measure'].attrs.update({'title': 'PEDS Value',
                                            'description': 'PEDS Value'
                                             })
    per_sub_stats['B:measure'].attrs.update(measure)
    per_sub_stats['n'].attrs.update(dict(title='count', description='Number of'
                                    'comparisons'))
    per_sub_stats['test-statistic'].attrs.update(
        dict(title='Iteration',
             description='Number of num_iterations that agree with the'
                         ' ALTERNATIVE hypothesis'))
    per_sub_stats['p-value'].attrs.update(dict(title='one-tailed',
                                          description='one-tail p-value'))
    per_sub_stats['q-value'].attrs.update(
        dict(title='Benjamini–Hochberg', description='FDR corrections using'
             ' Benjamini–Hochberg procedure'))
    return per_sub_stats


def _global_stats(p_series):
    """Creates global PEDS stats

    Uses Stouffers Method on the per-subject-values to create a global PEDS
    statistic.

    Parameters
    ----------
    p_series: pd.Series
        Series of all per-subject-stats

    Returns
    -------
    global_stats: pd.DataFrame
        a Dataframe containing global stat statistics, p-values, and q-values.
        Note: this is the DataFrame that is returned
        to the user as the global-stats.qza

    Examples
    --------
    >>> p_series = pd.Series([data = [0.001, 0.001, 0.001, 0.001])

    >>> _global_stats(p_series)

    pd.Dataframe({"measure": ["p-values"],
                   "n": [4],
                   "test-statistic": [20.93],
                   "p-value": [0.001]
                   "q-value": [NaN]})
    """
    stats, p = combine_pvalues(p_series, method="stouffer")
    global_stats = pd.DataFrame([["p-values", p_series.size, stats, p,
                                  np.nan]],
                                columns=['Measure', 'n', 'test-statistic',
                                         'p-value', 'q-value'])
    global_stats['Measure'].attrs.update({'title': ('p-value'),
                                          'description': ("p-value")})
    global_stats['test-statistic'].attrs.update(dict(title="Stouffer's",
                                                description="Stouffer's Z"
                                                " score method"))
    global_stats['p-value'].attrs.update(dict(title='one-tail p-value',
                                         description='one-tail p-value'))
    global_stats['q-value'].attrs.update(
        dict(title='Benjamini–Hochberg', description='FDR corrections using'
             'Benjamini–Hochberg procedure'))
    return global_stats
