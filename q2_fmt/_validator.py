# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from qiime2.plugin import ValidationError
from q2_fmt.plugin_setup import plugin
from q2_fmt import GroupDist, Ordered, Unordered, Matched


@plugin.register_validator(GroupDist[Ordered | Unordered, Matched])
def validate_unique_subjects_within_group(data: pd.DataFrame, level):
    if 'subject' not in data.columns:
        raise ValidationError('Subject column not found in data.')

    for group_id, group_df in data.groupby('group'):
        if group_df['subject'].duplicated().any():
            dupes = list(group_df['subject'][group_df['subject'].duplicated()])
            raise ValidationError(
                'Unique subject found more than once within an individual'
                ' group. Group(s) where duplicated subject was found:'
                f' {group_id}, Duplicated subjects:, {dupes}')
