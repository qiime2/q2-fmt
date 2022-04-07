# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import scipy.stats

# C : T_n
# dist_a = GroupDist[Ordered | Unordered, Independent]
# dist_optional = GroupDist[Ordered | Unordered, Matched | Independent]

# hypothesis: {'reference', 'all-pairwise'}

# (dist, hypothesis, alpha=0.05, reference=None, against_matched=None, previous_tests=list[teststatics])

def mann_whitney_u(distribution: pd.DataFrame, hypothesis: str,
                   reference_group: str=None,
                   against_each: pd.DataFrame=None) -> pd.DataFrame:
    # hypothesis = 'reference'
    #   check if reference is set
    #   if against_each is set, then it will be reference vs against_each
    #   else it will be reference vs other's in distribution
    #
    # hypothesis = 'all-pairwise'
    #   reference should be None
    #   if against_each is set, then each group in distribution vs against_each
    #   else it will be distribution vs distribution where group_i != group_j
    return pd.DataFrame()

# T_baseline : T_n
# T_n-1 : T_n
# GroupDist[Ordered, Matched]

# (dist, hypothesis, alpha=0.05, baseline=None)
def wilcoxon_srt(distribution: pd.DataFrame, hypothesis: str,
                 baseline_group: str=None) -> pd.DataFrame:
    # hypothesis = 'baseline'
    #   check if baseline is set
    #   compare each group in distribution (except baseline) against baseline
    #
    # hypothesis = 'consecutive'
    #   baseline should be None
    #   compare each group against the last one
    return pd.DataFrame()


# output for both:
# baseline
# a   b   statistic   p-val   q-val
# 0   1   ....
# 0   2   ....
# 0   3   ....
#
# consecutive
# 0   1   ...
# 1   2   ...
# 2   3   ...

# control
# controla   0
# controla   1
# controla   2
# controla   3
