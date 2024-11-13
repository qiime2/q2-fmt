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

from q2_fmt._util import (_rename_features,
                          _check_reference_column, _check_for_time_column,
                          _check_subject_column, _filter_associated_reference,
                          _check_duplicate_subject_timepoint,
                          _drop_incomplete_timepoints,
                          _drop_incomplete_subjects,
                          _check_column_type,
                          _create_recipient_table,
                          _create_masking, _create_mismatched_pairs,
                          _create_duplicated_tables,
                          _per_subject_stats, _global_stats, _mask_recipient,
                          _get_to_baseline_ref, _create_used_references,
                          _median, _subsample, _check_rarefaction_parameters)

from q2_stats.util import json_replace
from q2_stats.plots.raincloud import _make_stats


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
                             " `fmt pprf`. This is not compatible with"
                             " statistics created from `fmt pedf-simulation`"
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
            warnings.warn('PRDF was selected as the proportion metric, which'
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


def pedf(table: pd.DataFrame, metadata: qiime2.Metadata,
         time_column: str, reference_column: str, subject_column: str,
         filter_missing_references: bool = False,
         sampling_depth: int = None,
         num_resamples: int = 0) -> (pd.DataFrame):

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
    _check_rarefaction_parameters(num_resamples=num_resamples,
                                  sampling_depth=sampling_depth)

    num_list = []
    denom_list = []

    if num_resamples == 0:
        # This will allow the next block of code to run one with
        # no samping depth
        num_resamples = 1
    for x in range(num_resamples):
        pedf_df = pd.DataFrame(columns=['id',
                                        'transfered_donor_features',
                                        'total_donor_features',
                                        'donor', 'subject',
                                        'group'])
        if sampling_depth:
            table = _subsample(table, sampling_depth)
        pedf_df = _compute_proportion(df=pedf_df, type="Sample",
                                      time=np.nan,
                                      reference_series=used_references,
                                      table=table, metadata=metadata_df,
                                      time_column=time_column,
                                      subject_column=subject_column,
                                      reference_column=reference_column)
        # Set the index for pd matching when concating and summing
        pedf_df = pedf_df.set_index('id')
        num_list.append(pedf_df['transfered_donor_features'])
        denom_list.append(pedf_df['total_donor_features'])

    numerator_df = pd.concat(num_list, axis=1)
    denominator_df = pd.concat(denom_list, axis=1)

    median_numerator_series = _median(numerator_df)
    median_denominator_series = _median(denominator_df)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pedf = median_numerator_series/median_denominator_series

    pedf_df['measure'] = pedf
    pedf_df['transfered_donor_features'] = median_numerator_series
    pedf_df['total_donor_features'] = median_denominator_series
    pedf_df = pedf_df.reset_index()

    pedf_df['id'].attrs.update({
        'title': metadata_df.index.name,
        'description': 'Sample IDs'
    })
    pedf_df['measure'].attrs.update({
        'title': "Sample PEDF",
        'description': 'Proportional Engraftment of Donor Strains'
    })
    pedf_df['group'].attrs.update({
        'title': time_column,
        'description': 'Time'
    })
    pedf_df['subject'].attrs.update({
        'title': subject_column,
        'description': 'Subject IDs linking samples across time'
    })
    pedf_df['transfered_donor_features'].attrs.update({
        'title': "Transfered Reference Features",
        'description': '...'
    })
    pedf_df['total_donor_features'].attrs.update({
        'title': "Total Reference Features",
        'description': '...'
    })
    pedf_df['donor'].attrs.update({
        'title': reference_column,
        'description': 'Donor'
    })
    return pedf_df


def prdf(table: pd.DataFrame, metadata: qiime2.Metadata,
         time_column: str, reference_column: str, subject_column: str,
         filter_missing_references: bool = False,
         sampling_depth: int = None,
         num_resamples: int = 0) -> (pd.DataFrame):
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
    _check_rarefaction_parameters(num_resamples=num_resamples,
                                  sampling_depth=sampling_depth)

    num_list = []
    denom_list = []
    if num_resamples == 0:
        # This will allow the next block of code to run one with
        # no samping depth
        num_resamples = 1
    for x in range(num_resamples):
        prdf_df =\
            pd.DataFrame(columns=['id', 'recipients with feature',
                                  'all possible recipients with feature',
                                  'group', 'subject'])
        if sampling_depth:
            table = _subsample(table, sampling_depth)

        for time, time_metadata in metadata_df.groupby(time_column):
            prdf_df = _compute_proportion(df=prdf_df, type="Feature",
                                          time=time,
                                          reference_series=used_references,
                                          table=table,
                                          metadata=time_metadata,
                                          time_column=time_column,
                                          subject_column=subject_column,
                                          reference_column=reference_column)
        prdf_df = prdf_df.set_index('id')
        num_list.append(prdf_df['recipients with feature'])
        denom_list.append(prdf_df['all possible recipients with feature'])

    numerator_df = pd.concat(num_list, axis=1)
    denominator_df = pd.concat(denom_list, axis=1)
    median_numerator_series = _median(numerator_df)
    median_denominator_series = _median(denominator_df)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        prdf = median_numerator_series/median_denominator_series

    prdf_df['measure'] = prdf.values
    prdf_df['recipients with feature'] = median_numerator_series.values
    prdf_df['all possible recipients with feature'] = \
        median_denominator_series.values
    prdf_df = prdf_df.reset_index()

    prdf_df['id'].attrs.update({
        'title': "Feature ID",
        'description': ''
    })
    prdf_df['measure'].attrs.update({
        'title': "Feature PEDS",
        'description': 'Proportional Engraftment of Donor Strains'
    })
    prdf_df['group'].attrs.update({
        'title': time_column,
        'description': 'Time'
    })
    prdf_df['subject'].attrs.update({
        'title': "Feature ID",
        'description': ''
    })
    return prdf_df


def _compute_proportion(df: pd.Series, type: str, time: int,
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
    if type == "Sample" or type == "PPRS":
        num_sum = np.sum(maskedrecip, axis=1)
        donor_sum = np.sum(donormask, axis=1)
        for count, sample in enumerate(recip_df.index):
            sample_row = metadata.loc[sample]
            df.loc[len(df)] = [sample, num_sum[count],
                               donor_sum[count],
                               sample_row[reference_column],
                               sample_row[subject_column],
                               sample_row[time_column]]

    elif type == "Feature":
        num_sum = np.sum(maskedrecip, axis=0)
        donor_sum = np.sum(donormask, axis=0)
        for count, feature in enumerate(recip_df.columns):
            df.loc[len(df)] = [feature, num_sum[count],
                               donor_sum[count], time, feature]
            df = df.dropna()
    else:
        raise KeyError('There was an error finding which PEDS methods to use')
    return df


def pprf(table: pd.DataFrame, metadata: qiime2.Metadata,
         time_column: str, baseline_timepoint: str, subject_column: str,
         filter_missing_references: bool, sampling_depth: int = None,
         num_resamples: int = 0) -> (pd.DataFrame):
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
    baseline_metadata = metadata_df.join(used_references)

    num_list = []
    denom_list = []
    if num_resamples == 0:
        # This will allow the next block of code to run one with
        # no samping depth
        num_resamples = 1
    for x in range(num_resamples):
        pprf_df = pd.DataFrame(columns=['id',
                                        'transfered_baseline_features',
                                        'total_baseline_features',
                                        'baseline', 'subject',
                                        'group'])
        if sampling_depth:
            table = _subsample(table, sampling_depth)
        pprf_df = _compute_proportion(df=pprf_df, type='PPRS',
                                      time=np.nan,
                                      reference_series=used_references,
                                      table=table, metadata=baseline_metadata,
                                      time_column=time_column,
                                      subject_column=subject_column,
                                      reference_column=used_references.name)
        # Set the index for pd matching when concating and summing
        pprf_df = pprf_df.set_index('id')

        num_list.append(pprf_df['transfered_baseline_features'])
        denom_list.append(pprf_df['total_baseline_features'])

    numerator_df = pd.concat(num_list, axis=1)
    denominator_df = pd.concat(denom_list, axis=1)

    median_numerator_series = _median(numerator_df)
    median_denominator_series = _median(denominator_df)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pprf = median_numerator_series/median_denominator_series

    pprf_df['measure'] = pprf
    pprf_df['transfered_baseline_features'] = median_numerator_series
    pprf_df['total_baseline_features'] = median_denominator_series
    pprf_df = pprf_df.reset_index()

    pprf_df['id'].attrs.update({
        'title': metadata_df.index.name,
        'description': 'Sample IDs'
    })
    pprf_df['measure'].attrs.update({
        'title': 'PPRS',
        'description': 'Proportional Persistence of Recipient Strains'
    })
    pprf_df['group'].attrs.update({
        'title': time_column,
        'description': 'Time'
    })
    pprf_df['subject'].attrs.update({
        'title': subject_column,
        'description': 'Subject IDs linking samples across time'
    })
    pprf_df['transfered_baseline_features'].attrs.update({
        'title': "Transfered Reference Features",
        'description': '...'
    })
    pprf_df['total_baseline_features'].attrs.update({
        'title': "Total Reference Features",
        'description': '...'
    })
    pprf_df['baseline'].attrs.update({
        'title': used_references.name,
        'description': 'recipeint baseline'
    })
    return pprf_df


def pedf_permutation_test(table: pd.DataFrame, metadata: qiime2.Metadata,
                          time_column: str, reference_column: str,
                          subject_column: str, sampling_depth: int,
                          filter_missing_references: bool = False,
                          pedf_rarefaction: bool = True,
                          num_resamples: int = 999,
                          ) -> (pd.DataFrame, pd.DataFrame, pd.DataFrame):

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

    if not pedf_rarefaction:
        actual_pedf_num_resamples = 1
    else:
        actual_pedf_num_resamples = num_resamples

    pedf_df = pedf(
                   table=table, metadata=metadata,
                   time_column=time_column,
                   reference_column=reference_column,
                   subject_column=subject_column,
                   filter_missing_references=filter_missing_references,
                   num_resamples=actual_pedf_num_resamples,
                   sampling_depth=sampling_depth)

    actual_pedf = pedf_df[['id', 'measure']].set_index('id')['measure']

    # Mismatch simulation:
    recip_df = _create_recipient_table(used_references, metadata_df, table)
    donor_df = table[table.index.isin(used_references)]
    mismatched_series = \
        _create_mismatched_pairs(recip_df,
                                 metadata_df,
                                 used_references,
                                 reference_column)
    simulated_mismatched_series = mismatched_series.sample(n=num_resamples,
                                                           axis=0,
                                                           replace=True)
    simulated_recip_table, simulated_donor_table =\
        _create_duplicated_tables(simulated_mismatched_series, recip_df,
                                  donor_df)
    # concating or recip and donor tables so column number stays the same after
    # subsampling
    simulated_table = pd.concat([simulated_recip_table,
                                 simulated_donor_table])
    rarefied_simulated_table = _subsample(table=simulated_table,
                                          sampling_depth=sampling_depth)

    rarefied_recip_simulated_table =\
        rarefied_simulated_table.filter(items=simulated_recip_table.index,
                                        axis=0)

    rarefied_donor_simulated_table =\
        rarefied_simulated_table.filter(items=simulated_donor_table.index,
                                        axis=0)

    rarefied_donor_simulated_table = rarefied_donor_simulated_table > 0
    rarefied_recip_simulated_table = rarefied_recip_simulated_table > 0
    donor_mask = rarefied_donor_simulated_table.to_numpy()
    recip_mask = _mask_recipient(donor_mask, rarefied_recip_simulated_table)
    # Numerator for PEDS Calc. (Number of Donor features in the Recipient)
    num_engrafted_donor_features = np.sum(recip_mask.values, axis=1)
    # Denominator for PEDS Calc. (Number of unique features in the Donor)
    num_donor_features = np.sum(donor_mask, axis=1)
    # This ignores warnings that come from dividing by 0.
    # mismatched_pedf will be Nan if the denominator is 0 and thats reasonable.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mismatched_pedf = num_engrafted_donor_features/num_donor_features
    per_sub_stats = _per_subject_stats(mismatched_pedf,
                                       actual_pedf)
    global_stats = _global_stats(per_sub_stats['p-value'])

    return pedf_df, per_sub_stats, global_stats
