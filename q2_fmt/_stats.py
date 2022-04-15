# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import itertools
import pandas as pd
import scipy.stats


def mann_whitney_u(distribution: pd.DataFrame, hypothesis: str,
                   reference_group: str=None,
                   against_each: pd.DataFrame=None,
                   p_val_approx: str='auto') -> pd.DataFrame:
    dists = [distribution]
    if against_each is not None:
        dists.append(against_each)

    if hypothesis == 'reference':
        comparisons = _comp_reference(distribution, reference_group,
                                      against_each=against_each)
    elif hypothesis == 'all-pairwise':
        if reference_group is not None:
            raise ValueError("")
        comparisons = _comp_all_pairwise(distribution,
                                         against_each=against_each)
    else:
        raise ValueError()

    table = []
    for (idx_a, comp_a), (idx_b, comp_b) in comparisons:
        a_dist = dists[idx_a]
        b_dist = dists[idx_b]

        group_a = a_dist[a_dist['group'] == comp_a]['measure']
        group_b = b_dist[b_dist['group'] == comp_b]['measure']

        row = _compare_mannwhitneyu(group_a, group_b, p_val_approx)
        row['A:group'] = comp_a
        row['B:group'] = comp_b
        table.append(row)

    df = pd.DataFrame(table)

    df['q-value'] = _fdr_correction(df['p-value'])

    df = df[['A:group', 'A:n', 'A:measure', 'B:group', 'B:n', 'B:measure',
             'n', 'test-statistic', 'p-value', 'q-value']]

    return df


def _comp_reference(distribution, reference_group, against_each=None):
    group = distribution['group']
    reference_group = _get_reference_from_column(group, reference_group,
                                                 'reference_group')

    if against_each is None:
        for other in group[group != reference_group].unique():
            yield ((0, reference_group), (0, other))
    else:
        for other in against_each['group'].unique():
            yield ((0, reference_group), (1, other))


def _comp_all_pairwise(distribution, against_each=None):
    if against_each is None:
        for (comp_a, comp_b) in itertools.combinations(
                distribution['group'].unique(), 2):
            yield ((0, comp_a), (0, comp_b))
    else:
        for (comp_a, comp_b) in itertools.product(
                distribution['group'].unique(),
                against_each['group'].unique()):
            yield ((0, comp_a), (1, comp_b))


def _compare_mannwhitneyu(group_a, group_b, p_val_approx):
    stat, p_val = scipy.stats.mannwhitneyu(
        group_a, group_b, method=p_val_approx, nan_policy='raise')

    return {
        'A:n': len(group_a),
        'B:n': len(group_b),
        'A:measure': group_a.median(),
        'B:measure': group_b.median(),
        'n': len(group_a) + len(group_b),
        'test-statistic': stat,
        'p-value': p_val
    }


def wilcoxon_srt(distribution: pd.DataFrame, hypothesis: str,
                 baseline_group: str=None, p_val_approx: str='auto') -> pd.DataFrame:

    if hypothesis == 'baseline':
        comparisons = _comp_baseline(distribution, baseline_group)
    elif hypothesis == 'consecutive':
        if baseline_group is not None:
            raise ValueError()
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
    group = distribution['group']
    baseline_group = _get_reference_from_column(group, baseline_group,
                                                'baseline_group')

    for comp_b in group[group != baseline_group].unique():
        yield (baseline_group, comp_b)


def _comp_consecutive(distribution):
    group = distribution['group']
    timepoints = list(sorted(group.unique()))
    yield from zip(timepoints, timepoints[1:])


def _compare_wilcoxon(group_a, group_b, p_val_approx) -> dict:
    if p_value_approx == 'asymptotic':
        # wilcoxon differs from mannwhitneyu in arg value, but does the same
        # test using a normal dist instead of the permutational dist so
        # normalize the naming in Q2
        p_value_approx = 'approx'
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


def _get_reference_from_column(series, reference_value, param_name):
    if reference_value is None:
        raise ValueError("%s must be provided." % param_name)

    if (series == reference_value).any():
        return reference_value
    elif (series == float(reference_value)).any():
        return float(reference_value)
    else:
        raise ValueError("%r was not found as a group within the distribution."
                         % reference_value)
