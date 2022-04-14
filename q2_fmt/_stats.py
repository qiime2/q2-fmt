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
                 baseline_group: str=None, p_val_approx: str='auto') -> pd.DataFrame:

    if hypothesis == 'baseline':
        comparisons = _comp_baseline(distribution, baseline_group)
    elif hypothesis == 'consecutive':
        if baseline_group is not None:
            raise Exception()
        comparisons = _comp_consecutive(distribution)
    else:
        raise ValueError()

    table = []
    for comp_a, comp_b in comparisons:
        group_a = distribution[distribution['group'] == comp_a]
        group_b = distribution[distribution['group'] == comp_b]

        group_a = group_a.set_index('subject')['measure']
        group_b = group_b.set_index('subject')['measure']
        row = _compare_wilcoxon(group_a, group_b, p_val_approx)
        row['A:group'] = comp_a
        row['B:group'] = comp_b
        table.append(row)

    df = pd.DataFrame(table)

    df['q-value'] = _fdr_correction(df['p-value'])

    df = df[['A:group', 'A:n', 'A:measure', 'B:group', 'B:n', 'B:measure',
             'n', 'test-statistic', 'p-value', 'q-value']]

    return df

def _comp_baseline(distribution, baseline_group):
    baseline_group = float(baseline_group)
    group = distribution['group']
    if baseline_group is None:
        raise ValueError()

    if not (group == baseline_group).any():
        print(group)
        print(group == baseline_group)
        raise ValueError()

    for comp_b in group[group != baseline_group].unique():
        yield (baseline_group, comp_b)


def _comp_consecutive(distribution):
    group = distribution['group']
    timepoints = list(sorted(group.unique()))
    yield from zip(timepoints, timepoints[1:])


def _compare_wilcoxon(group_a, group_b, p_val_approx) -> dict:
    comp = pd.merge(group_a.to_frame(), group_b.to_frame(), how='outer',
                    left_index=True, right_index=True)
    filtered = comp.dropna()

    results = {
        'A:n': len(group_a),
        'B:n': len(group_b),
        'A:measure': group_a.median(),
        'B:measure': group_b.median(),
        'n': len(filtered.index),
    }

    stat, p_val = scipy.stats.wilcoxon(
        filtered.iloc[:, 0], filtered.iloc[:, 1],
        nan_policy='raise', mode=p_val_approx, alternative='two-sided')

    results['test-statistic'] = stat
    results['p-value'] = p_val

    return results

def _fdr_correction(series):
    return series

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
