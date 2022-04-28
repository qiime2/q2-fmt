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


def plot_rainclouds(output_dir: str, data: pd.DataFrame):
    J_ENV = jinja2.Environment(
        loader=jinja2.PackageLoader('q2_fmt', 'assets')
    )

    x_label = data['measure'].attrs['unit']
    y_label = data['group'].attrs['unit']
    title = f'{x_label} across {y_label}'

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
        fh.write(index.render(spec=spec_string))


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
