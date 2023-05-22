# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
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
from collections import Counter


def plot_heatmap(output_dir: str, data: pd.DataFrame,
                 level_delimiter: str = None):
    J_ENV = jinja2.Environment(
        loader=jinja2.PackageLoader('q2_fmt', 'assets')
    )

    x_label = data['group'].attrs['title']
    measure = data['measure'].attrs['title']
    subject_title_temp = data['subject'].attrs['title']
    if ("recipients with feature" in data.columns and
            level_delimiter is not None):
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
                # be complete thrown out.
                fields = [sub]
            subject_seen.append(sub)
            most_specific = fields[-1]
            if most_specific in seen and sub not in subject_seen:
                y_labels.append(f"{seen[most_specific]}: {most_specific} *")
            else:
                y_labels.append(most_specific)
            seen[most_specific] += 1
        data['y_label'] = y_labels

        data['id'] = [i.replace(level_delimiter, ' ') for i in data['id']]
        data['subject'] = [i.replace(level_delimiter, ' ')
                           for i in data['subject']]

        data['y_label'].attrs.update({
                'title': "Feature ID",
                'description': ''})
    else:
        data['y_label'] = data["subject"]
        data['y_label'].attrs.update({
                'title': subject_title_temp,
                'description': ''})

    y_label = data['y_label'].attrs['title']
    title = f'{measure} of {y_label} across {x_label}'

    index = J_ENV.get_template('index.html')
    data = json.loads(data.to_json(orient='records'))

    spec_fp = pkg_resources.resource_filename(
        'q2_fmt', os.path.join('assets', "spec.json")
    )
    with open(spec_fp) as fh:
        json_obj = json.load(fh)

    full_spec = json_replace(json_obj, data=data, x_label=x_label,
                             y_label=y_label, title=title, measure=measure)

    with open(os.path.join(output_dir, "index.html"), "w") as fh:
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
