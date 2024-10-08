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

import os
import pkg_resources
import jinja2
import json

from qiime2 import Metadata
from q2_fmt._util import (json_replace, _rename_features,
                          _check_reference_column, _check_for_time_column,
                          _check_subject_column, _filter_associated_reference,
                          _check_duplicate_subject_timepoint,
                          _drop_incomplete_timepoints,
                          _drop_incomplete_subjects,
                          _check_column_type,
                          _create_recipient_table,
                          _create_masking, _create_mismatched_pairs,
                          _create_duplicated_recip_table, _create_sim_masking,
                          _per_subject_stats, _global_stats, _mask_recipient,
                          _get_to_baseline_ref, _create_used_references)
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
            filter_missing_references=filter_missing_references)

    else:
        feature_peds = ctx.get_action('fmt', 'feature_peds')
        peds_dist = feature_peds(
            table=table, metadata=metadata, time_column=time_column,
            subject_column=subject_column, reference_column=reference_column,
            filter_missing_references=filter_missing_references)
    results += heatmap(data=peds_dist[0],
                       level_delimiter=level_delimiter,
                       drop_incomplete_subjects=drop_incomplete_subjects,
                       drop_incomplete_timepoints=drop_incomplete_timepoints)

    return tuple(results)


def heatmap(output_dir: str, data: pd.DataFrame,
            level_delimiter: str = None,
            per_subject_stats: pd.DataFrame = None,
            global_stats: pd.DataFrame = None,
            drop_incomplete_timepoints: list = None,
            drop_incomplete_subjects: bool = None):
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
        if drop_incomplete_subjects or drop_incomplete_timepoints:
            warnings.warn('Feature PEDS was selected as the PEDS metric, which'
                          ' does not accept `drop_incomplete_subjects` or'
                          ' `drop_incomplete_timepoints` as parameters. One'
                          ' (or both) of these parameters were detected in'
                          ' your input, and will be ignored.')
            drop_incomplete_timepoints = None
            drop_incomplete_subjects = False
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

    data = _drop_incomplete_timepoints(data, drop_incomplete_timepoints)
    data = _drop_incomplete_subjects(data, drop_incomplete_subjects)
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
                filter_missing_references: bool = False) -> (pd.DataFrame):

    # making sure that samples exist in the table
    ids_with_data = table.index
    metadata = metadata.filter_ids(ids_to_keep=ids_with_data)
    column_properties = metadata.columns
    metadata_df = metadata.to_dataframe()

    time_col = _check_for_time_column(metadata_df, time_column)
    _check_column_type(column_properties, "time",
                       time_column, "numeric")
    metadata_df = metadata_df.filter(items=time_col.index, axis=0)
    reference_series = _check_reference_column(metadata_df, reference_column)
    _check_column_type(column_properties, "reference",
                       reference_column, "categorical")
    used_references = _create_used_references(reference_series, metadata_df,
                                              time_column)
    # return things that should be removed
    metadata_df, used_references = \
        _filter_associated_reference(used_references, metadata_df,
                                     filter_missing_references, ids_with_data)
    subject_series = _check_subject_column(metadata_df, subject_column)
    _check_column_type(column_properties, "subject",
                       subject_column, "categorical")
    _check_duplicate_subject_timepoint(subject_series, metadata_df,
                                       subject_column, time_column)

    peds_df = pd.DataFrame(columns=['id', 'measure',
                                    'transfered_donor_features',
                                    'total_donor_features', 'donor', 'subject',
                                    'group'])
    peds_df = _compute_peds(peds_df=peds_df, peds_type="Sample",
                            peds_time=np.nan, reference_series=used_references,
                            table=table, metadata=metadata_df,
                            time_column=time_column,
                            subject_column=subject_column,
                            reference_column=reference_column)
    return peds_df


def feature_peds(table: pd.DataFrame, metadata: qiime2.Metadata,
                 time_column: str, reference_column: str, subject_column: str,
                 filter_missing_references: bool = False) -> (pd.DataFrame):
    # making sure that samples exist in the table
    ids_with_data = table.index
    metadata = metadata.filter_ids(ids_to_keep=ids_with_data)
    column_properties = metadata.columns
    metadata_df = metadata.to_dataframe()

    _check_for_time_column(metadata_df, time_column)
    _check_column_type(column_properties, "time",
                       time_column, "numeric")
    reference_series = _check_reference_column(metadata_df, reference_column)
    _check_column_type(column_properties, "reference",
                       reference_column, "categorical")
    used_references = _create_used_references(reference_series, metadata_df,
                                              time_column)
    metadata_df, used_references = \
        _filter_associated_reference(used_references, metadata_df,
                                     filter_missing_references, ids_with_data)
    _check_subject_column(metadata_df, subject_column)
    _check_column_type(column_properties, "subject",
                       subject_column, "categorical")
    peds_df = pd.DataFrame(columns=['id', 'measure', 'recipients with feature',
                                    'all possible recipients with feature',
                                    'group', 'subject'])
    for time, time_metadata in metadata_df.groupby(time_column):
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

    donormask = _create_masking(metadata_df=metadata, donor_df=donor_df,
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
            peds_type = "Sample PEDS"
            transfered = "transfered_donor_features"
            total = 'total_donor_features'
            ref = 'donor'
            measure_description = 'Proportional Engraftment of Donor Strains'

        peds_df['id'].attrs.update({
            'title': metadata.index.name,
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
            'title': "Feature PEDS",
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


def sample_pprs(table: pd.DataFrame, metadata: qiime2.Metadata,
                time_column: str, baseline_timepoint: str, subject_column: str,
                filter_missing_references: bool) -> (pd.DataFrame):
    # making sure that samples exist in the table
    ids_with_data = table.index
    metadata = metadata.filter_ids(ids_to_keep=ids_with_data)
    # This retains types from the QIIME2 Metadata.
    column_properties = metadata.columns
    metadata_df = metadata.to_dataframe()

    time_col = _check_for_time_column(metadata_df, time_column)
    _check_column_type(column_properties, 'time',
                       time_column, 'numeric')
    metadata_df = metadata_df.filter(items=time_col.index, axis=0)

    used_references =\
        _get_to_baseline_ref(time_col=metadata_df[time_column],
                             baseline_timepoint=baseline_timepoint,
                             time_column_name=time_column,
                             subject_column_name=subject_column,
                             metadata=Metadata(metadata_df))
    metadata_df, used_references = \
        _filter_associated_reference(used_references, metadata_df,
                                     filter_missing_references, ids_with_data)
    subject_series = _check_subject_column(metadata_df, subject_column)
    _check_column_type(column_properties, 'subject',
                       subject_column, 'categorical')
    _check_duplicate_subject_timepoint(subject_series, metadata_df,
                                       subject_column, time_column)

    peds_df = pd.DataFrame(columns=['id', 'measure',
                                    'transfered_baseline_features',
                                    'total_baseline_features', 'baseline',
                                    'subject', 'group'])
    baseline_metadata = metadata_df.join(used_references)
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
                    num_iterations: int = 999) -> (pd.DataFrame, pd.DataFrame):

    ids_with_data = table.index
    metadata = metadata.filter_ids(ids_to_keep=ids_with_data)

    metadata_df = metadata.to_dataframe()
    reference_series = _check_reference_column(metadata_df,
                                               reference_column)
    used_references = _create_used_references(reference_series, metadata_df,
                                              time_column)
    (metadata_df,
     used_references) = _filter_associated_reference(used_references,
                                                     metadata_df,
                                                     filter_missing_references,
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
           filter_missing_references=filter_missing_references
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
