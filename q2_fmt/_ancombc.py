# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
from qiime2 import Metadata

from q2_fmt._util import _check_for_time_column, _check_reference_column


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
