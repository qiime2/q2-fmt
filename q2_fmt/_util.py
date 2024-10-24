# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
from collections import Counter

from scipy.stats import false_discovery_control, combine_pvalues

import random
import itertools
import biom

import numpy as np
import qiime2


# PEDS/PPRS Data Prep Helper Methods
def _check_column_missing(metadata, column_name, param_name, e):
    """Checks if column is in the metdata

    Checks that column is in the metdata and creates helpful error messages.

    Parameters
    ----------
    metadata: pd.Dataframe
        Study `Metadata`

    column_name: str
       `Metadata` column name

    param_name: str
        qiime2 parameter that above column name was provided for. This allows
        for better erroring.

    e: str
        traceback from try accept that this method is called from


    Examples
    --------
    >>> metadata = pd.Dataframe({'id': ['sample1', 'sample2'],
                   'time': [1, 1]}).set_index('id')

    >>> column_name = 'time'

    >>> param_name = 'group'

    >>> e = '[Tracback Here]'

    >>> _check_column_missing(metadata, column_name, param_name, e)
    """
    if column_name == metadata.index.name:
        raise KeyError('The `--p-%s-column` input provided was the'
                       ' same as the index of the metadata.'
                       ' `--p-%s-column` can not be the same as the'
                       ' index of metadata:'
                       ' `%s`' % (param_name, param_name,
                                  column_name)) from e
    else:
        raise KeyError('There was an error finding the provided'
                       ' `--p-%s-column`: `%s` in the metadata'
                       % (param_name, column_name)) from e


def _check_reference_column(metadata, reference_column_name):
    """Checks if refernce column exists and returns it as a series if it
    exists.

    Checks if refernce colum exists and returns it as a series if it exists.
    Errors using helper method `_check_column_missing`

    Parameters
    ----------
    metadata: pd.Dataframe
        Study `Metadata`

    reference_column_name: str
       reference column name inside `Metadata`

    Returns
    ------
    reference_series: pd.series
        Series with recipient sample id's as the index and
        their donor or "reference" as the value of the series.

    Examples
    --------
    >>> metadata = pd.Dataframe({'id': ['sample1', 'sample2'],
                   'reference': ['donor1', 'donor1'],
                   'time': [1,2]}).set_index('id')

    >>> column_name = 'reference'

    >>> _check_reference_column(metadata, reference_column_name)

    pd.Series(data=['donor1', 'donor1'],
              index=pd.Index(['sample1','sample2'], name='id'),
              name='reference')


    """
    try:
        reference_series = metadata[reference_column_name]
    except Exception as e:
        _check_column_missing(metadata, reference_column_name, 'reference', e)
    return reference_series


def _check_for_time_column(metadata, time_column_name):
    """Checks if time column exists and returns it as a series if it exists.

    Checks if time column exists and returns it as a series if it exists.
    Errors using helper method `_check_column_missing`.

    Parameters
    ----------
    metadata: pd.Dataframe
        Study `Metadata`

    time_column_name: str
       time column name inside `Metadata`

    Returns
    -------
    num_timepoints: int
        number of timepoints in the metadata

    time_col: pd.Series
        pandas series with recpient samples as the index and
        timepoints as the values.

    Examples
    --------
    >>> metadata = pd.Dataframe({'id': ['sample1', 'sample2'],
                   'reference': ['donor1', 'donor1'],
                   'time': [1,2]}).set_index('id')

    >>> column_name = 'time'

    >>> _check_for_time_column(metadata, time_column_name)

    time_series = pd.Series(data=['1', '2'],
                         index=pd.Index(['sample1','sample2'], name='id'),
                         name='reference')

    """
    try:
        time_series = metadata[time_column_name].dropna()
    except Exception as e:
        _check_column_missing(metadata, time_column_name, 'time', e)
    return time_series


def _check_subject_column(metadata, subject_column_name):
    """Checks if subject column exists and returns it as a series if it exists.

    Checks if subjuct column exists and returns it as a series if it exists.
    Errors using helper method `_check_column_missing`

    Parameters
    ----------
    metadata: pd.Dataframe
        Study `Metadata`

    subject_column_name: str
        subject column name inside `Metadata`

    Returns
    -------
    subject_series: pd.series
        Series with recipient sample id's as the index and
        their subject as the value of the series.

    Examples
    --------
    >>> metadata = pd.Dataframe({'id': ['sample1', 'sample2'],
                   'reference': ['donor1', 'donor1'],
                   'time': [1,2],
                   'subject': ['sub1','sub1']}).set_index('id')

    >>> subject_column_name = 'subject'

    >>> _check_subject_column(metadata, subject_column_name)

    pd.Series(data=['sub1', 'sub1'],
              index=pd.Index(['sample1','sample2'], name='id'),
              name='reference')

    """
    try:
        subject_series = metadata[subject_column_name]
    except Exception as e:
        _check_column_missing(metadata, subject_column_name, "subject", e)
    return subject_series


def _create_used_references(reference_series, metadata_df, time_column_name):
    """creates used references pd.Series.

    creates a used reference pd.Series by checking which recipients have time
    values

    Parameters
    ----------
    reference_series: pd.series
        Series with recipient sample id's as the index
        and their donor or "reference" as the value of the series.

    metadata_df: pd.dataframe
        Study `Metadata`

    time_column_name: str
        time column name inside `Metadata`

    Returns
    --------
    used_references: pd.Series
        Filtered Series with recipient sample id's as the index and
        their donor or "reference" as the value of the series. Filtered to
        exclude missing references if filter_missing_references is true.

    Examples
    --------
    >>> reference_series = pd.Series(data=['donor1', 'donor1'],
                                    index=pd.Index(['sample1',
                                                    'sample2'], name='id'),
                                    name='reference')

    >>> metadata_df = pd.DataFrame({'id': ['sample1', 'sample2'],
                   'reference': ['donor1', 'donor1'],
                   'time': [1,np.NaN]}).set_index('id')
    >>> time_column_name = 'time'

    >>>  _create_used_references(reference_series, metadata_df,
                                 time_column_name)

        pd.Series(data=['donor1'],
                  index=pd.Index(['sample1'], name='id'),
                  name='reference')
    """
    used_references = reference_series[~metadata_df[time_column_name].isna()]
    return used_references


def _filter_associated_reference(used_references, metadata_df,
                                 filter_missing_references,
                                 ids_with_data):
    """Errors on/Filters references that are missing in metadata and or
    in the associated feature-table

    Errors if there are samples in metadata that have a timepoint but no donor.
    OR if there are references mentioned in the metadata that are not found in
    the associated.
    If filter_missing_references is true this method will filter out both cases
    listed above.
    Returns filtered metadata and a list of refernces to be used in downstream
    analysis

    Parameters
    ----------
    used_references: pd.series
        Series with recipient sample id's (with timepoint values) as the index
        and their donor or "reference" as the value of the series.

    metadata_df: pd.dataframe
        Study `Metadata`

    time_column_name: str
        time column name inside `Metadata`

    filter_missing_references: bool
        Boolean for whether or not to filter

    id_with_data: pd.Index
        list(index) of all the sample ids in `Metadata` that are found in the
        associated feature table.

    Returns
    --------
    metadata_df: pd.dataframe
        Study `Metadata` that is filtered to exclude missing references if
        filter_missing_references is true.

    used_references: pd.Series
        Filtered Series with recipient sample id's as the index and
        their donor or "reference" as the value of the series. Filtered to
        exclude missing references if filter_missing_references is true.

    Examples
    --------
    >>> used_references = pd.Series(data=['donor1', np.Nan],
                                    index=pd.Index(['sample1',
                                                    'sample2'], name='id'),
                                    name='reference')

    # Note: sample2 should be filtered out because it has no
    reference.
    >>> metadata_df = pd.Dataframe({'id': ['sample1', 'sample2'],
                   'reference': ['Donor1', np.Nan ],
                   'time': [1,2]}).set_index('id')

    >>> filter_missing_references = True

    >>> id_with_data = Index(['sample1', 'sample2', 'donor1'])

    >>>  _filter_associated_reference(used_references, metadata_df,
                                      filter_missing_references,
                                      ids_with_data)
    metadata_df =  pd.Dataframe({'id': ['sample1'],
                   'reference': ['Donor1'],
                   'time': [1]}).set_index('id')

    used_references = pd.Series(data=['donor1'],
                                index=pd.Index(['sample1'], name='id'),
                                name='reference')
    """
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
    metadata_df = metadata_df.filter(items=used_references.index, axis=0)
    return metadata_df, used_references


def _check_duplicate_subject_timepoint(subject_series, metadata,
                                       subject_column_name, time_column_name):
    """Checks if recipients in a subject group occur more than once in a
    timepoint.

    Errors if there is more than one sample in a timepoint within a subject
    group

    Parameters
    ----------
    subject_series: pd.Series
        Series with recipient sample id's as the index and
        their subject as the value of the series.

    metadata: pd.Dataframe
        Study `Metadata`

    subject_column_name: str
        subject column name inside `Metadata`

    time_column_name: str
        time column name inside `Metadata`


    Examples
    --------
    >>> pd.Series(data=['sub1', 'sub1'],
                  index=pd.Index(['sample1','sample2'], name='id'),
                  name='subject')

    >>> metadata = pd.Dataframe({'id': ['sample1', 'sample2'],
                   'reference': ['Donor1', 'Donor1'],
                   'time': [1,2],
                   'subject': ['sub1','sub1']}).set_index('id')

    >>> subject_column_name = 'subject'

    >>> time_column_name = 'time'

    >>>  _check_duplicate_subject_timepoint(subject_series, metadata,
                                            subject_column_name,
                                            time_column_name)
    """
    for subject in subject_series:
        subject_df = metadata[metadata[subject_column_name] == subject]
        if not subject_df[time_column_name].is_unique:
            timepoint_list = subject_df[time_column_name].to_list()
            raise ValueError('There is more than one occurrence of a subject'
                             ' in a timepoint. All subjects must occur only'
                             ' once per timepoint. Subject %s appears in '
                             ' timepoints: %s' % (subject, timepoint_list))


def _check_column_type(column_properties, param_name, column_name,
                       column_type):
    """Checks if metadata column type

    Checks that column is the appropriate type. For example: Time columns are
    expected to be numeric.

    Parameters
    ----------
   column_properties: OrderedDict
        Dictionary with all columns names as the key and their type as a value.

    param_name: str
        qiime2 parameter that above column name was provided for. This allows
        for better erroring.

    column_name: str
       `Metadata` column name

    column_type: str
        string of the type to check that the column is.
        (Categorical or Numeric)


    Examples
    --------
    >>> column_properties =  OrderedDict([('SubjectID',
                                   ColumnProperties(type='categorical',
                                   missing_scheme='blank')),
                                 ('time_point',
                                   ColumnProperties(type='numeric',
                                   missing_scheme='blank')),
                                 ('Donor', ColumnProperties(type='categorical',
                                   missing_scheme='blank'))])


    >>> param_name = 'group'

    >>> column_name = 'time_point'

    >>>  column_type = "numeric"

    >>> _check_column_type(column_properties, param_name, column_name,
                           column_type)

    """
    try:
        assert column_properties[column_name].type == column_type
    except AssertionError as e:
        raise AssertionError('Non-%s values found in `--p-%s-column`.'
                             ' Please make sure the column selected contains'
                             ' the correct MetadataColumn type. Column with'
                             ' non-%s values that was'
                             ' selected: `%s`' % (column_type, param_name,
                                                  column_type,
                                                  column_name)) from e


# Porportion Calculations Methods (for PPRS and PEDS)
def _create_recipient_table(reference_series, metadata_df, table_df):
    """Creates a feature table with just the recipient samples.

   This method takes a pd.series (reference_series) which contains the
   recipient and donor pairs, study `Metadata` which has been filtered down
   to relevant recipients and the study feature table and returns a feature
   table that is filtered down to relevant recipients.

    Parameters
    ----------
    reference_series: pd.Series
        A series with recipients as the index and the associated reference
        as the values.
    metadata_df: pd.DataFrame
        The sample metadata for the study filtered down to recipients with
        relevant timepoints.
    table_df: pd.DataFrame
        the study feature table.

    Returns
    -------
    recip_df: pd.DataFrame
        A DataFrame feature table of FMT recipients.

    Examples
    --------
    >>> reference_series = pd.Series(data=['donor1', 'donor2', 'donor3'],
                                     index=pd.Index(['sample1', 'sample2',
                                                    'sample3'], name='id'),
                                     name='Ref')

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

    >>> table_df = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',  'donor1', 'donor2',
                   'donor3'],
            'Feature1': [1, 0, 0, 1, 0, 0],
            'Feature2': [0, 1, 0, 0, 1, 0],
            'Feature3': [0, 0, 10, 0, 1]}).set_index('id')

    >>> _create_recipient_table(recip_df, metadata_df, used_references,
                                 reference_column)
        pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3'],
            'Feature1': [1, 0, 0],
            'Feature2': [0, 1, 0],
            'Feature3': [0, 0, 1]}).set_index('id')
    """
    subset_reference_series = \
        reference_series[reference_series.index.isin(metadata_df.index)]
    recip_df = table_df[table_df.index.isin(subset_reference_series.index)]
    return recip_df


def _create_masking(metadata_df, donor_df, recip_df, reference_column):
    """Create a donor mask to mask recipient features that aren't in the donor.

    Creates a Donor Numpy array that duplicates donor samples to match up
    with table. This will mask recipient feature that aren't
    in the donor.

    Parameters
    ----------
    metadata_df: pd.DataFrame
        Study `Metadata`
    donor_df: pd.DataFrame
        A feature table of FMT donors.
    recip_df: pd.DataFrame
        A feature table of recipient donors.
    reference_column: str
        Name of the reference column in the Sample Metadata.

    Returns
    -------
    donor_mask: ndarray
        a numpy array of donor feature information with the same order as
        duplicated_recip_table. This allows for easy np array math.

    Examples
    --------
    >>> metadata_df = pd.DataFrame({
              'id': ['donor1', 'donor2', 'donor3', 's1', 's2', 's3'],
              'reference': [np.Nan, np.Nan, np.Nan, 'donor1', 'donor2',
                            'donor3'],
              'time': [np.Nan, np.Nan, np.Nan, 1, 1, 1],
              'Feature3':[np.Nan, np.Nan, np.Nan, 'sub1', 'sub2',
                            'sub3'],}).set_index('id')

    >>> donor_df = pd.DataFrame({
              'id': ['donor1', 'donor2', 'donor3'],
              'Feature1': [1, 0, 0],
              'Feature2': [0, 1, 0],
              'Feature3': [0, 0, 1]}).set_index('id')

    >>> recip_df = pd.DataFrame({
              'id': ['s1', 's2', 's3'],
              'Feature1': [1, 0, 0],
              'Feature2': [0, 1, 0],
              'Feature3': [0, 0, 1]}).set_index('id')

    >>> reference_column = 'reference'

    >>> _create_masking(metadata_df, donor_df, recip_df, reference_column)

    ndarray[[1, 0, 0],
           [0, 1, 0],
           [0, 0, 1]]
    """
    donor_index_masking = []
    for sample in recip_df.index:
        donor = metadata_df.loc[sample, reference_column]
        donor_index_masking.append(donor_df.index.get_loc(donor))
    donor_df = donor_df.to_numpy()
    donor_mask = donor_df[donor_index_masking]
    donor_mask = donor_mask.astype(int)
    return donor_mask


def _mask_recipient(donor_mask, recip_df):
    """Uses np array math to multiple the donor_mask on to the recipient.

    Uses np array math to multiple the donor_mask on to the recipient.
    This will mask recipient feature that aren't
    in the donor.

    Parameters
    ----------
    donor_mask: ndarray
        numpy array of donor feature table
    recip_df: pd.DataFrame
        A feature table of recipient donors.

    Returns
    -------
    masked_recip: ndarray
        a numpy array of donor feature in recipients

    Examples
    --------

    >>> donor_df = ndarray[[1, 0, 0],
                          [0, 1, 0],
                          [0, 0, 1]]

    >>> recip_df = pd.DataFrame({
              'id': ['s1', 's2', 's3'],
              'Feature1': [1, 0, 0],
              'Feature2': [0, 1, 0],
              'Feature3': [0, 0, 1]}).set_index('id')

    >>> _mask_recipient(donor_mask, recip_df)

    ndarray[[1, 0, 0],
           [0, 1, 0],
           [0, 0, 1]]
    """

    masked_recip = donor_mask & recip_df
    return masked_recip


def _median(df):
    """Calculates the median of rows in a dataframe.

    takes a pd.DataFrame and calculates the median for each sample(row) in the
    dataframe and returns a series with median values per sample

    Parameters
    ----------
    df: pd.DataFrame
        pd.DataFrame with values for each iteration of rarefaction

    Returns
    -------
    median_series: pd.Series
        pd.Series with median values from rarefactions

    Examples
    --------

    >>> df = pd.DataFrame({
              'id': ['s1', 's2', 's3'],
              'transfered_donor_features': [1, 1, 1],
              'transfered_donor_features0': [0, 1, 0.5],
              'transfered_donor_features1': [0.75, 0.50, .5]}).set_index('id')

    >>> _median(df)

        pd.Series(data=[.75, 1, 0.5],
                  index=pd.Index(['s1', 's2', 's3'], name='id'))
    """
    median_series = df.median(axis=1)
    return median_series


def _subsample(table, sampling_depth):
    """ subsamples feature table

    takes a pd.DataFrame feature table and transforms it to a biom table and
    uses bioms subsample function and returns a rarefied feature-table

    Parameters
    ----------
    df: pd.DataFrame
        the study feature-table
    sampling_depth: int
        number of observation to subsample each sample(row)

    Returns
    -------
    table: pd.DataFrame
        pd.DataFrame rarified feature table

    Examples
    --------

    >>> table = pd.DataFrame({
              'id': ['s1', 's2', 's3'],
              'feature1': [10, 10, 10],
              'feature2': [2, 1, 5],
              'feature3': [1, 5, 1]}).set_index('id')

    sampling_depth = 7

    >>> _subsample(table, sampling_depth)

       pd.DataFrame({
              'id': ['s1', 's2', 's3'],
              'feature1': [3, 3, 1],
              'feature2': [2, 1, 4],
              'feature3': [1, 5, 1]}).set_index('id')
"""

    table_biom = biom.Table(table.T.values,
                            sample_ids=table.T.columns.to_list(),
                            observation_ids=table.T.index.to_list())

    subsampled_table_biom = table_biom.subsample(sampling_depth,
                                                 axis='sample',
                                                 by_id=False)
    table = subsampled_table_biom.to_dataframe(True).T
    return table


def _check_rarefaction_parameters(num_resamples, sampling_depth):
    """ Checks if both rarefaction parameters have values or are None/0

    Checks that if sampling depth is none, num_resamples is 0 or that both have
    values.

    Parameters
    ----------
    num_resamples: int
        Number of iterations to preform rarefaction
    sampling_depth: int
        Number of oberservations to subsample each sample to.


    Examples
    --------
    >>> num_resamples = 2

    >>> sampling_depth = 100

    >>> _check_rarefaction_parameters(num_resamples, sampling_depth)
    """
    if (num_resamples == 0) != (sampling_depth is None):
        raise AssertionError('`num_resamples` and `sampling depth` parameters'
                             ' must be passed in together. In order to run'
                             ' rarefaction, a `sampling_depth`'
                             ' (how many observations to resample to)'
                             ' and `num_resamples` (how many iterations of'
                             ' resampling) must be provided. To rarefy, please'
                             ' provided `1` as the value for`num_resamples`.'
                             ' If you would not like to rarefy or preform'
                             ' rarefaction, please provide `0` as the value'
                             ' for `num_resamples` and `none` as the value for'
                             ' `sampling_depth`')


# Heatmap Helper Methods
def _rename_features(level_delimiter, data: pd.DataFrame):
    """Splits Feature Labels using a provided level_delimiter.

    Uses a level delimiter to split a feature label in place.
    This method checks to see if the provide `data` is a `feature-peds` data
    structure and will not split the label if it is not a `feature-peds` data
    structure. It expects the label to split to be in the index of the data.

    Parameters
    ----------
    level_delimiter: str
        Character to split index string on.

    data: pd.dataframe
        PEDS dataframe to split label on

    Examples
    --------
    >>> data = pd.Dataframe({'id': ['Genus1;species1', 'Genus2;species2'],
                   'measure': [.5, .75],
                   'group': [1, 1],
                   'recipients with feature': [10, 15]
                   'all possible recipients with feature': [20, 20]
                   'subject': ['Genus1;species1', 'Genus2;species2']})

    >>> level_delimiter = ';'

    Note there is no return but this is what data will look like at the end:
    >>> _rename_features(level_delimiter, data: pd.DataFrame)
        pd.Dataframe({'id': [species1, species2],
                   'measure': [.5, .75],
                   'group': [1, 1],
                   'recipients with feature': [10, 15]
                   'all possible recipients with feature': [20, 20]
                   'subject': [Genus1;species1, Genus2;species2]})
    """

    if ('recipients with feature' in data.columns and
            level_delimiter is not None):
        temp_group_name = data['group'].attrs['title']
        temp_subject_name = data['subject'].attrs['title']
        temp_measure_name = data['measure'].attrs['title']

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
                y_labels.append(f'{seen[most_specific]}: {most_specific} *')
            else:
                y_labels.append(most_specific)
            seen[most_specific] += 1
        data['subject'] = y_labels

        data['id'] = [i.replace(level_delimiter, ' ') for i in data['id']]

        # currently attrs get deleted with df is changed. right now the best
        # way to solve this is by saving them as temp and saving them at the
        # end

        data['subject'].attrs.update({'title': temp_subject_name,
                                      'description': ''})
        data['group'].attrs.update({'title': temp_group_name,
                                    'description': ''})
        data['measure'].attrs.update({'title': temp_measure_name,
                                      'description': ''})


def _drop_incomplete_timepoints(data, drop_incomplete_timepoints):
    """Filters Dist1D using a list of timepoints to drop

    Takes a Dist1D (viewable as a pandas Dataframe) and a list of timepoints to
    drop and filters the Dist1D. This should help create beautiful heatmaps!

    Parameters
    ----------
    data: pd.Dataframe
        The Dist1D produced by a portion method(feature, or sample PEDS or
        PPRS) to filter.
    drop_incomplete_timepoints: list[Str]
        The list of timepoints to filter on.

    Returns
    -------
    data: pd.DataFrame
        A filtered Dist1D for the heatmap

    Examples
    --------
    >>> data = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'donor1', 'donor2'],
            'measure': [.4, .5, 1,
                        np.nan, np.nan],
            'subject': ['sub1', 'sub1', 'sub2',
                        np.nan, np.nan],
            'group': [1, 2, 1,
                      np.nan, np.nan],
            'transfered_donor_features': [4, 5, 10,
                                          np.nan, np.nan],
            'total_donor_features': [10, 10, 10,
                                     np.nan, np.nan]}).set_index('id')

    >>> drop_incomplete_timepoints = ['2']

    >>> _drop_incomplete_timepoints(data, drop_incomplete_timepoints)

        pd.DataFrame({
            'id': ['sample1', 'sample3',
                   'donor1', 'donor2'],
            'measure': [.4, 1,
                        np.nan, np.nan],
            'subject': ['sub1', 'sub2',
                        np.nan, np.nan],
            'group': [1, 1,
                      np.nan, np.nan],
            'transfered_donor_features': [4, 10,
                                          np.nan, np.nan],
            'total_donor_features': [10, 10,
                                     np.nan, np.nan]}).set_index('id')
    """
    if drop_incomplete_timepoints:
        for time in drop_incomplete_timepoints:
            data = data[data['group'] != float(time)]
    return data


def _drop_incomplete_subjects(data, drop_incomplete_subjects):
    """Filters Dist1D so that there are no incomplete subjects

    Takes a Dist1D (viewable as a pandas Dataframe) and a bool of whether or
    not to drop incomplete subjects and filters the Dist1D.
    This should help create beautiful heatmaps!

    Parameters
    ----------
    data: pd.Dataframe
        The Dist1D produced by a portion method(feature, or sample PEDS or
        PPRS) to filter.
    drop_incomplete_subjects: Bool
        Whether or not to  drop incomplete subjects.

    Returns
    -------
    data: pd.DataFrame
        A filtered Dist1D for the heatmap

    Examples
    --------
    >>> data = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3',
                   'donor1', 'donor2'],
            'measure': [.4, .5, 1,
                        np.nan, np.nan],
            'subject': ['sub1', 'sub1', 'sub2',
                        np.nan, np.nan],
            'group': [1, 2, 1,
                      np.nan, np.nan],
            'transfered_donor_features': [4, 5, 10,
                                          np.nan, np.nan],
            'total_donor_features': [10, 10, 10,
                                     np.nan, np.nan]}).set_index('id')

    >>> drop_incomplete_subjects = True

    >>> _drop_incomplete_subjects(data, drop_incomplete_timepoints)
        pd.DataFrame({
            'id': ['sample1', 'sample2',
                   'donor1', 'donor2'],
            'measure': [.4, .5,
                        np.nan, np.nan],
            'subject': ['sub1', 'sub1', 'sub2',
                        np.nan, np.nan],
            'group': [1, 2,
                      np.nan, np.nan],
            'transfered_donor_features': [4, 5,
                                          np.nan, np.nan],
            'total_donor_features': [10, 10,
                                     np.nan, np.nan]}).set_index('id')


    """
    num_timepoints = data['group'].unique().size
    subject_occurrence_series = (data['subject'].value_counts())
    if (subject_occurrence_series < num_timepoints).any():
        if drop_incomplete_subjects:
            subject_to_keep = (subject_occurrence_series[
                                subject_occurrence_series ==
                                num_timepoints].index)
            data = data[data['subject'].isin(subject_to_keep)]
    return data


# Helper Method for Heatmap Spec
# TODO: this will be imported from stats soon
def json_replace(json_obj, **values):
    """
    Search for elements of `{"{{REPLACE_PARAM}}": "some_key"}` and replace
    with the result of `values["some_key"]`.
    """
    if type(json_obj) is dict and list(json_obj) == ["{{REPLACE_PARAM}}"]:
        param_name = json_obj["{{REPLACE_PARAM}}"]
        return values[param_name]

    if type(json_obj) is list:
        return [json_replace(x, **values) for x in json_obj]

    elif type(json_obj) is dict:
        return {key: json_replace(value, **values)
                for key, value in json_obj.items()}

    else:
        return json_obj


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


# CC/PPRS Helper Function
def _get_to_baseline_ref(time_col, baseline_timepoint, time_column_name,
                         subject_column_name, metadata):
    """Creates a reference column where the baseline samples are the reference
    and returns a pd.series

    Matches recipient samples to their baseline sample and returns a
    pd.series


    Parameters
    ----------
    time_column_name: str
        time column name inside `Metadata`
    baseline_timepoint: str
        string of the baseline timepoint
    time_column: pd.Series
        pandas series with recpient samples as the index and
        timepoints as the values.
    subject_column_name: str
        subject column name inside `Metadata`
    metadata: qiime2.Metadata
        Study `Metadata`

    Returns
    -------
    used_references: pd.Series
        Filtered Series with recipient sample id's as the index and
        their "reference" (baseline) as the value of the series.

    Examples
    --------
    >>> time_col = pd.Series([1, 2, 3, 1, 2, 3], index=['sample1', 'sample2',
                                                        'sample3', 'sample4',
                                                        'sample5', 'sample6'])
    >>> baseline_timepoint = "1"
    >>> time_column_name = "group"
    >>> subject_column_name = "subject"
    >>> metadata = pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'sample5', 'sample6'],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', 'sub2', 'sub2'],
            'group': [1, 2, 3, 1, 2, 3]}).set_index('id')

    >>> _get_to_baseline_ref(time_col, baseline_timepoint, time_column,
                                   subject_column, Metadata(metadata))
    pd.Series(['sample1', 'sample1', 'sample4', 'sample4'],
                            index=['sample2', 'sample3', 'sample5',
                                   'sample6'])
    """

    temp_baseline_ref = []
    reference_list = []
    baseline_ref_df = pd.DataFrame()
    # All valid FMT samples have to have a time column
    metadata = metadata.to_dataframe()[~time_col.isna()]
    if float(baseline_timepoint) not in metadata[time_column_name].values:
        raise AssertionError('The provided baseline timepoint'
                             f' {baseline_timepoint} was not'
                             f' found in `metadata` '
                             f' column {time_column_name}.')
    for sub, samples in metadata.groupby([subject_column_name]):
        reference = \
            samples[samples[
                time_column_name] == float(baseline_timepoint)].index.to_list()
        if len(reference) > 1:
            raise ValueError('More than one baseline sample was found per'
                             ' subject. Only one baseline sample can be'
                             ' used as a reference. Please group baseline'
                             ' replicates.')
        elif len(reference) == 0:
            # If there is no baseline for a subject,
            # This will either drop with filter-missing-references or
            # or error and say that they need to pass
            # filter-missing-references
            reference = [np.nan]
        temp_baseline_ref = temp_baseline_ref + samples.index.to_list()
        reference_list = \
            reference_list + (reference * len(samples.index.to_list()))
    # I dont see any way that this hits because of my above assertion but
    # I think its a good check so I am leavig it in.
    if len(reference_list) == 0:
        raise AssertionError('No baseline samples'
                             ' were found in the metadata.'
                             ' Please confirm that a valid'
                             ' baseline timepoint was given.')
    if pd.Series(reference_list).isnull().all():
        raise AssertionError('No baseline samples'
                             ' were connected via subject.'
                             ' Confirm that all subjects have a'
                             ' baseline timepoint')
    baseline_ref_df['sample_name'] = temp_baseline_ref
    baseline_ref_df['relevant_baseline'] = reference_list
    baseline_ref_df = \
        baseline_ref_df[~baseline_ref_df['sample_name'].isin(
            reference_list)].set_index('sample_name')
    reference_col = _get_series_from_col(
        md=qiime2.Metadata(baseline_ref_df), col_name='relevant_baseline',
        param_name='reference_column',
        expected_type=qiime2.CategoricalMetadataColumn)
    # this is so the variables for distance to donor and distance to
    # baseline have the same variable name
    used_references = reference_col
    return used_references


# HELPER FUNCTION FOR DATA Filtering in CC and _get_to_baseline_ref
def _get_series_from_col(md, col_name, param_name, expected_type=None,
                         drop_missing_values=False):
    """Creates a pd.Series from qiime2.Metadata column.

    Creates a pd.Series from qiime2.Metadata column. This functions also check
    if the metadata column provide is on the right type (Either Categorical
    or Numeric). This functions also drops missing values if needed.

    Parameters
    ----------
    md: qiim2.Metadata
        Study `Metadata`
    col_name: str
        subject column name inside `Metadata`
    expected_type: (None| qiime2.NumericMetadataColumn |
                    qiime2.CategoricalMetadataColumn)
        qiime2 metadata type to check column type against
    drop_missing_samples: bool
        boolean for filtering Nans

    Returns
    -------
    pd.Series
        Filtered Series with recipient sample id's as the index and
        metadata column of interest as the value of the series.

    Examples
    --------
    >>> md = Metadata(pd.DataFrame({
            'id': ['sample1', 'sample2', 'sample3', 'sample4',
                   'sample5', 'sample6'],
            'subject': ['sub1', 'sub1', 'sub1', 'sub2', 'sub2', 'sub2'],
            'group': [1, 2, 3, 1, 2, 3]}).set_index('id'))
    >>> column_name = "group"
    >>>  expected_type = qiime2.NumericMetadataColumn
    >>> drop_missing_samples = False


    >>> _get_series_from_col(md, col_name, param_name, expected_type,
                             drop_missing_samples)
    pd.Series([1, 2, 3, 1, 2, 3],
              index=['sample1', 'sample2', 'sample3', 'sample4',
                                   'sample5', 'sample6'])
    """
    try:
        column = md.get_column(col_name)
    except ValueError as e:
        raise ValueError("There was an issue with the argument for %r. %s"
                         % (param_name, e)) from e

    if expected_type is not None and not isinstance(column, expected_type):
        if type(expected_type) is tuple:
            exp = tuple(e.type for e in expected_type)
        else:
            exp = expected_type.type

        raise ValueError("Provided column for %r is %r, not %r."
                         % (param_name, column.type, exp))

    if drop_missing_values:
        column = column.drop_missing_values()

    return column.to_series()


# HELPER FUNCTION FOR sorting a multi-index (for distance matrix and metadata)
def _sort_multi_index(index):
    """sorts a multi-index (for distance matrix and metadata) dataframe

    Takes an multi-index as input and returns the sorted multi-index

    Parameters
    ----------
    index: pd.MultiIndex
        Multi-index from a `Distance Matrix`


    Returns
    -------
    sorted_multi: pd.MultiIndex
        Sorted Multi-index from a `Distance Matrix`


    Examples
    --------
    >>> index = MultiIndex([('bar', 'one'),
            ('bar', 'two'),
            ('baz', 'one'),
            ('baz', 'two'),
            ('foo', 'one'),
            ('foo', 'two'),
            ('qux', 'one'),
            ('qux', 'two')])

    >>> _sort_multi_index(index)

    MultiIndex([('bar', 'one'),
            ('bar', 'two'),
            ('baz', 'one'),
            ('baz', 'two'),
            ('foo', 'one'),
            ('foo', 'two'),
            ('one', 'qux'),
            ('qux', 'two')])

    """

    sorted_levels = list(map(sorted, index))
    sorted_multi = pd.MultiIndex.from_tuples(sorted_levels)
    return sorted_multi


# HELPER FUNCTION FOR Dists1D[Ordered | NestedOrdered, Matched | Independent]
def _ordered_dists(diversity_measure: pd.Series, is_beta,
                   used_references, time_col, subject_col, group_col):
    """Creates the output dataframe of group timepoints.

    Creates  Dists1D[Ordered | NestedOrdered, Matched | Independent] from
    diversity metrics and pd.Series containing reference info, time info,
    subject info and group info.

    Parameters
    ----------
    diversity_measure: pd.Series
        Diveristy metric that group-timepoints is preparing with metadata
    is_beta: bool
        boolean to discribe if the data is a beta distance matix or a
        alpha diverity metrix
    used_references: pd.Series
        Filtered Series with recipient sample id's as the index and
        their donor or "reference" as the value of the series. Filtered to
        exclude missing references if filter_missing_references is true.
    time_col: pd.Series
        pandas series with recpient samples as the index and
        timepoints as the values.
    subject_col: pd.series
        pandas series with recipient sample id's as the index and
        their subject as the value of the series.
    group_col: pd.Series
        pandas series with recipient sample id's as the index and
        their group as the value of the series.

    Returns
    -------
    ordinal_df: pd.Dataframe
        Dataframe output of  Dists1D[Ordered | NestedOrdered, Matched |
        Independent]


    Examples
    --------
    >>> diversity_measure = pd.Series(data=[24, 37,
                                        15, 6,
                                        44, 17,
                                        29, 32,
                                        51, 3],
                                 index=pd.Index(['sampleA', 'sampleB',
                                                 'sampleC', 'sampleD',
                                                 'sampleE', 'sampleF',
                                                 'sampleG', 'donor1',
                                                 'donor2', 'donor3'],
                                                 name='id'),
                                 name='measure')
    >>> is_beta = False
    >>> used_references = pd.Series(data=['donor1', 'donor2',
                                      'donor1', 'donor1',
                                      'donor2', 'donor2',
                                      'donor3'],
                                 index=pd.Index(['sampleA', 'sampleB',
                                                 'sampleC', 'sampleD',
                                                 'sampleE', 'sampleF',
                                                 'sampleG', 'donor1',
                                                 'donor2', 'donor3'],
                                                 name='id'),
                                 name='relevant_donor')
    >>> time_col = pd.Series(data=[7.0, 7.0,
                               9.0, 11.0,
                               11.0, 9.0,
                               7.0, Np.Nan,
                               Np.Nan, Np.Nan],
                         index=pd.Index(['sampleA', 'sampleB',
                                         'sampleC', 'sampleD',
                                         'sampleE', 'sampleF',
                                         'sampleG', 'donor1',
                                         'donor2', 'donor3'],
                                         name='id'),
                         name='days_post_transplant')
    >>> subject_col = pd.Series(data=['subject1', 'subject2',
                               'subject1', 'subject1',
                               'subject2', 'subject2',
                               'subject3', Np.Nan,
                               Np.Nan, Np.Nan],
                         index=pd.Index(['sampleA', 'sampleB',
                                         'sampleC', 'sampleD',
                                         'sampleE', 'sampleF',
                                         'sampleG', 'donor1',
                                         'donor2', 'donor3',
                                         'donor4'], name='id'),
                         name='subject')
    >>> group_col = None

    >>> _ordered_dists(diversity_measure: pd.Series, is_beta,
                   used_references, time_col, subject_col, group_col):
    pd.Dataframe({"id": ['sampleA', 'sampleB','sampleC', 'sampleD',
                         'sampleE', 'sampleF','sampleG']
                   "measure" : [24, 37, 15, 6,
                                44, 17, 29],
                   "group" : [7.0, 7.0,9.0, 11.0,
                              11.0, 9.0, 7.0],
                   "subject" : ['subject1', 'subject2', 'subject1', 'subject1',
                               'subject2', 'subject2', 'subject3']})
    """
    if is_beta:
        idx = pd.MultiIndex.from_frame(
            used_references.to_frame().reset_index())
        idx = _sort_multi_index(idx)
        idx.names = ['id', 'reference']
        diversity_measure.index.names = ['id', 'reference']
    else:
        idx = used_references.index
        idx.name = 'id'

    try:
        sliced_df = (diversity_measure[idx]
                     .to_frame()
                     .reset_index()
                     .set_index('id')
                     )
    except KeyError:
        raise KeyError('Pairwise comparisons were unsuccessful. Please double'
                       ' check that your chosen reference column contains'
                       ' values that are also present in the ID column for'
                       ' the associated metadata.')

    if is_beta:
        sliced_df.index = used_references.index
        sliced_df.index.name = 'id'

    ordinal_df = sliced_df[['measure']]
    ordinal_df['group'] = time_col

    if subject_col is not None:
        ordinal_df['subject'] = subject_col

    if group_col is not None:
        ordinal_df['class'] = group_col.name
        ordinal_df['level'] = group_col

    return ordinal_df.reset_index()


# HELPER FUNCTION FOR Dists[Unordered, Independent]
def _independent_dists(diversity_measure, metadata,
                       used_references, is_beta, used_controls):
    """Creates the independent output dataframe of group timepoints.

    Creates  Dists[Unordered, Independent] from diversity metrics, metadata
    and pd.Series containing reference and control info.

    Parameters
    ----------
    diversity_measure: pd.Series
        Diveristy metric that group-timepoints is preparing with metadata
    is_beta: bool
        boolean to discribe if the data is a beta distance matix or a
        alpha diverity metrix
    used_references: pd.Series
        Filtered Series with recipient sample id's as the index and
        their donor or "reference" as the value of the series. Filtered to
        exclude missing references if filter_missing_references is true.
    metadata: qiime2.Metadata
        Study `Metadata`
    used_controls: pd.Series
        pandas series with recipient sample id's as the index and
        their relevant control as the value of the series.

    Returns
    -------
    nominal_df: pd.Dataframe
        Dataframe output of  Dists1D[Ordered | NestedOrdered, Matched |
        Independent]


    Examples
    --------
    >>> diversity_measure = pd.Series(data=[24, 37,
                                        15, 6,
                                        44, 17,
                                        29, 32,
                                        51, 3],
                                 index=pd.Index(['sampleA', 'sampleB',
                                                 'sampleC', 'sampleD',
                                                 'sampleE', 'sampleF',
                                                 'sampleG', 'donor1',
                                                 'donor2', 'donor3'],
                                                 name='id'),
                                 name='measure')
    >>> is_beta = False
    >>> used_references = pd.Series(data=['donor1', 'donor2',
                                      'donor1', 'donor1',
                                      'donor2', 'donor2',
                                      'donor3', np.Nan,
                                       np.Nan, np.Nan],
                                 index=pd.Index(['sampleA', 'sampleB',
                                                 'sampleC', 'sampleD',
                                                 'sampleE', 'sampleF',
                                                 'sampleG', 'donor1',
                                                 'donor2', 'donor3'],
                                                 name='id'),
                                 name='relevant_donor')
    >>> metadata = Metadata(pd.Dataframe({"id": ['sampleA', 'sampleB',
                                                 'sampleC', 'sampleD',
                                                 'sampleE', 'sampleF',
                                                 'sampleG', 'donor1',
                                                 'donor2', 'donor3']
                                "reference" : ['donor1', 'donor2',
                                                'donor1', 'donor1',
                                                'donor2', 'donor2',
                                                'donor3', np.Nan,
                                                np.Nan, np.Nan]]}))
    >>> used_controls = None

    >>> _independent_dists(diversity_measure: pd.Series, is_beta,
                   used_references, time_col, subject_col, group_col):

    pd.Dataframe({"id": ['donor1', 'donor2', 'donor3', 'donor4',
                         'sampleE', 'sampleF','sampleG']
                   "measure" : [32, 51, 3, 19],
                   "group" : ['reference', 'reference','reference',
                              'reference']})
    """
    unique_references = sorted(used_references.unique())

    if is_beta:
        if len(unique_references) > 1:
            ref_idx = pd.MultiIndex.from_tuples(
                itertools.combinations(unique_references, 2))
            ref_idx.names = ['A', 'B']

        else:
            ref_idx = pd.MultiIndex(levels=[[], []],
                                    codes=[[], []],
                                    names=['A', 'B'])

        diversity_measure.index.names = ['A', 'B']

        if used_controls is not None:
            grouped_md = (metadata
                          .to_dataframe()
                          .loc[used_controls.index]
                          .groupby(used_controls)
                          )
            ctrl_list = list()
            for group_id, grouped_ctrls in grouped_md:
                if len(grouped_ctrls.index) < 2:
                    continue
                ctrl_combos = list(
                    itertools.combinations(
                        grouped_ctrls.index, 2)
                )
                ctrl_idx = pd.MultiIndex.from_tuples(ctrl_combos)
                ctrl_series = pd.Series(group_id, index=ctrl_idx)
                ctrl_list.append(ctrl_series)

            if len(ctrl_list) >= 1:
                ctrl_series = pd.concat(ctrl_list)
                ctrl_series.index.names = ['A', 'B']

            else:
                ctrl_series = pd.Series()

            ctrl_series.name = 'group'

    else:
        ref_idx = list(unique_references)
        if used_controls is not None:
            ctrl_series = used_controls
            ctrl_series.index.name = 'id'

    try:
        nominal_df = diversity_measure[ref_idx].to_frame().reset_index()
    except KeyError:
        raise KeyError('Pairwise comparisons were unsuccessful. Please double'
                       ' check that your chosen reference column contains'
                       ' values that are also present in the ID column for'
                       ' the associated metadata.')

    nominal_df['group'] = 'reference'

    if used_controls is not None:
        ctrl_group = diversity_measure[ctrl_series.index].to_frame()
        ctrl_group['group'] = ctrl_series
        ctrl_group = ctrl_group.reset_index()
        nominal_df = pd.concat([nominal_df, ctrl_group])
        nominal_df = nominal_df.reset_index(drop=True)

    if 'A' in nominal_df.columns:
        if not nominal_df.empty:
            nominal_df['id'] = nominal_df[['A', 'B']].agg('..'.join, axis=1)
        else:
            nominal_df['id'] = []
        nominal_df = nominal_df[['id', 'measure', 'group', 'A', 'B']]

    return nominal_df
