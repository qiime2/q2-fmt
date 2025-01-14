# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
from qiime2 import Metadata

from q2_fmt._util import (_check_for_time_column, _check_reference_column,
                          _check_column_type, _check_subject_column,
                          _check_duplicate_subject_timepoint,
                          _create_used_references,
                          _filter_associated_reference
                          )


def detect_donor_indicators(ctx, table, reference_column, time_column,
                            baseline_timepoint, metadata,
                            level_delimiter=None):
    filter = ctx.get_action('feature_table', 'filter_samples')
    ancombc = ctx.get_action('composition', 'ancombc')
    da_barplot = ctx.get_action('composition', 'da_barplot')
    results = []

    _check_for_time_column(metadata.to_dataframe(),
                           time_column)
    _check_reference_column(metadata.to_dataframe(),
                            reference_column)

    ids_to_keep = get_baseline_donor_md(metadata=metadata,
                                        reference_column=reference_column,
                                        time_column=time_column,
                                        baseline_timepoint=baseline_timepoint)

    filtered_table, = filter(table=table,
                             metadata=Metadata(ids_to_keep))
    dataloaf, = ancombc(table=filtered_table, metadata=Metadata(ids_to_keep),
                        reference_levels=["type::donor"], formula='type')
    results.append(dataloaf)
    viz, = da_barplot(data=dataloaf, significance_threshold=0.05,
                      level_delimiter=level_delimiter)
    results.append(viz)
    return tuple(results)


def get_baseline_donor_md(metadata, reference_column, time_column,
                          baseline_timepoint):
    """Creates a metadata for differentiating baseline and donor
    ----------
    metadata: pd.Dataframe
        Study `Metadata`
    reference_column: str
       name of reference column in `Metadata` column
    time_column: str
       name of reference column in `Metadata` column
    baseline_timepoint: str
        timepoint that represents baseline
    Examples
    --------
    >>> metadata = pd.DataFrame({'id': ['sample1', 'sample2', 'donor1'],
                   'reference': ['donor1', 'donor1', np.nan],
                   'time': [1, 2, np.nan],
                   'subject': ['sub1','sub1', np.nan]}).set_index('id')
    >>> time_column = 'time'
    >>> reference_column = 'reference'
    >>>  baseline_timepoint = '1'
    >>> get_baseline_donor_md(metadata, reference_column, time_column,
                              baseline_timepoint)
        pd.DataFrame({'id': ['sample1', 'donor1'],
                   'reference': ['donor1', np.nan],
                   'time': [1,np.nan],
                   'subject': ['sub1',np.nan]}).set_index('id')
    """
    md_df = metadata.to_dataframe()
    ids_to_keep =\
        pd.Series(index=md_df[reference_column].dropna().unique(),
                  data='donor', name='type')
    ids_to_keep =\
        pd.concat([ids_to_keep,
                   pd.Series(index=md_df[md_df[time_column] ==
                                         float(baseline_timepoint)
                                         ].index.to_list(),
                             data='baseline', name=type)])

    ids_to_keep = ids_to_keep.to_frame()
    ids_to_keep.index.name = 'id'
    ids_to_keep = ids_to_keep.rename({0: "type"}, axis=1)
    return ids_to_keep


def indicator_tracking_prep(
    table: pd.DataFrame, metadata: Metadata, time_column: str,
    reference_column: str,
    subject_column: str, indicator_id: str,
    filter_missing_references: bool = False
) -> (pd.DataFrame):
    # making sure that samples exist in the table
    ids_with_data = table.index
    metadata = metadata.filter_ids(ids_to_keep=ids_with_data)
    column_properties = metadata.columns
    metadata_df = metadata.to_dataframe()

    time_col = _check_for_time_column(metadata_df, time_column)
    _check_column_type(column_properties, "time",
                       time_column, "numeric")
    metadata_df = metadata_df.filter(items=time_col.index, axis=0)
    subject_series = _check_subject_column(metadata_df, subject_column)
    _check_column_type(column_properties, "subject",
                       subject_column, "categorical")
    _check_duplicate_subject_timepoint(subject_series, metadata_df,
                                       subject_column, time_column)
    reference_series = _check_reference_column(metadata_df, reference_column)
    _check_column_type(column_properties, "reference",
                       reference_column, "categorical")
    used_references = _create_used_references(reference_series, metadata_df,
                                              time_column)
    # return things that should be removed
    metadata_df, used_references = \
        _filter_associated_reference(used_references, metadata_df,
                                     filter_missing_references, ids_with_data)
    try:
        measure = table[indicator_id]
    except KeyError:
        raise KeyError(f'{indicator_id} was not found in feature-table.'
                       ' Please check input feature-table and confirm that'
                       ' the feature of interest is in the feature table.'
                       ' This is commonly caused because the provided'
                       ' feature-table was not collapsed but a taxon string'
                       ' was provided as the indicator name')
    ordinal_dist = pd.DataFrame(data={'measure': measure,
                                      'group': time_col,
                                      'subject': subject_series},
                                index=used_references.index)
    return ordinal_dist.reset_index()
