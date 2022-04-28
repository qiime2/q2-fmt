# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources
import jinja2
import json
import pandas as pd


def plot_rainclouds(output_dir: str, data: pd.DataFrame,
                    stats: pd.DataFrame=None):
    table1 = None
    if stats is not None:
        table1, stats = _make_stats(stats)

    J_ENV = jinja2.Environment(
        loader=jinja2.PackageLoader('q2_fmt', 'assets')
    )

    x_label = data['measure'].attrs['unit']
    y_label = data['group'].attrs['unit']
    subject_unit = data['subject'].attrs['unit']
    title = f'{x_label} of {subject_unit} across {y_label}'
    figure1 = (
        f'Raincloud plots showing the distribution of subjects\''
        f' measure of {x_label} across {y_label}. Kernel density estimation'
        f' performed using a bandwidth calculated by Scott\'s method. Boxplots'
        f' show the min and max of the data (whiskers) as well as the first,'
        f' second (median), and third quartiles (box). '
        f' Points and connecting lines represent individual subjects'
        f' with a consistent jitter added across groups such that slopes'
        f' across adjacent groups are visually comparable between subjects.')

    index = J_ENV.get_template('index.html')
    data = json.loads(data.to_json(orient='records'))

    spec_fp = pkg_resources.resource_filename(
        'q2_fmt', os.path.join('assets', 'spec.json'))
    with open(spec_fp) as fh:
        json_obj = json.load(fh)

    full_spec = json_replace(json_obj,
                             data=data, x_label=x_label, y_label=y_label,
                             title=title)

    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        spec_string = json.dumps(full_spec)
        fh.write(index.render(spec=spec_string, stats=stats,
                              figure1=figure1, table1=table1))


def json_replace(json_obj, **values):
    """
    Search for elements of `{"{{REPLACE_PARAM}}": "some_key"}` and replace
    with the result of `values["some_key"]`.
    """
    if type(json_obj) is list:
        return [json_replace(x, **values) for x in json_obj]
    elif type(json_obj) is dict:
        new = {}
        for key, value in json_obj.items():
            if type(value) is dict and list(value) == ["{{REPLACE_PARAM}}"]:
                param_name = value["{{REPLACE_PARAM}}"]
                new[key] = values[param_name]
            else:
                new[key] = json_replace(value, **values)
        return new
    else:
        return json_obj


def _make_stats(stats):
    method = stats['test-statistic'].attrs['unit']
    group_unit = (stats['A:group'].attrs['unit']
                  + ' vs ' + stats['B:group'].attrs['unit'])
    pval_method = stats['p-value'].attrs['unit']
    qval_method = stats['q-value'].attrs['unit']
    table1 = (f'{method} tests between groups ({group_unit}), with'
              f' {pval_method} p-value calculations and {qval_method}'
              f' correction for multiple comparisons (q-value).')
    df = pd.DataFrame(index=stats.index)
    group_a = _make_group_col('A', stats)
    df[group_a.name] = group_a
    group_b = _make_group_col('B', stats)
    df[group_b.name] = group_b
    df['A'] = stats['A:measure']
    df['B'] = stats['B:measure']
    df = df.merge(stats.iloc[:, 6:], left_index=True, right_index=True)
    df.columns = pd.MultiIndex.from_tuples([
        ('Group A', stats['A:group'].attrs['unit']),
        ('Group B', stats['B:group'].attrs['unit']),
        ('A', stats['A:measure'].attrs['unit']),
        ('B', stats['B:measure'].attrs['unit']),
        ('', 'n'),
        ('', 'test-statistic'),
        ('', 'p-value'),
        ('', 'q-value'),
    ])
    html = df.to_html(index=False)

    return table1, html


def _make_group_col(prefix, df):
    group_series = df[prefix + ':group']
    group_n = df[prefix + ':n']

    if (group_series.dtype == float
            and group_series.apply(float.is_integer).all()):
        group_series = group_series.astype(int)

    group_series = group_series.apply(str)
    group_n = " (n=" + group_n.apply(str) + ")"

    series = group_series + group_n
    series.name = f"Group " + prefix
    return series
