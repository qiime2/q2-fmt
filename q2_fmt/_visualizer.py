# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from dataclasses import Field
import os
import pkg_resources
import jinja2
import json
import pandas as pd

from q2_fmt import GroupDist

def hello_world(output_dir: str):
    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write('hello world')

    template = pkg_resources.resource_filename('q2_fmt', os.path.join('assets', 'index.html'))

def plot_rainclouds(output_dir: str, data: pd.DataFrame):
    J_ENV = jinja2.Environment(
        loader=jinja2.PackageLoader('q2_fmt', 'assets')
    )

    index = J_ENV.get_template('index.html')
    data = json.loads(data.to_json(orient='records'))

    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write(index.render(spec=json.dumps(_spec_(data))))

def _spec_(data):
    return {
    "$schema": "https://vega.github.io/schema/vega/v5.json",
    "description": "A basic bar chart example, with value labels shown upon mouse hover.",
    "width": 400,
    "height": 200,
    "padding": 5,

    "data": [
        {
        "name": "table",
        "values": data
        },
        {
        "name": "median",
        "source": "table",
        "transform": [
            {"type": "aggregate", "groupby": ["subject"],
             "fields": ["measure"], "ops": ["median"]},
        ]
        }
    ],

    "signals": [
        {
        "name": "tooltip",
        "value": {},
        "on": [
            {"events": "rect:mouseover", "update": "datum"},
            {"events": "rect:mouseout",  "update": "{}"}
        ]
        }
    ],

    "scales": [
        {
        "name": "xscale",
        "type": "band",
        "domain": {"data": "table", "field": "subject"},
        "range": "width",
        "padding": 0.05,
        "round": True
        },
        {
        "name": "yscale",
        "domain": {"data": "table", "field": "measure"},
        "nice": True,
        "range": "height"
        }
    ],

    "axes": [
        { "orient": "bottom", "scale": "xscale" },
        { "orient": "left", "scale": "yscale" }
    ],

    "marks": [
        {
        "type": "rect",
        "from": {"data":"median"},
        "encode": {
            "enter": {
            "x": {"scale": "xscale", "field": "subject"},
            "width": {"scale": "xscale", "band": 1},
            "y": {"scale": "yscale", "field": "median_measure"},
            "y2": {"scale": "yscale", "value": 0}
            },
            "update": {
            "fill": {"value": "steelblue"}
            },
            "hover": {
            "fill": {"value": "magenta"}
            }
        }
        },
        {
        "type": "text",
        "encode": {
            "enter": {
            "align": {"value": "center"},
            "baseline": {"value": "bottom"},
            "fill": {"value": "#333"}
            },
            "update": {
            "x": {"scale": "xscale", "signal": "tooltip.subject", "band": 0.5},
            "y": {"scale": "yscale", "signal": "tooltip.median_measure", "offset": -2},
            "text": {"signal": "tooltip.median_measure"},
            "fillOpacity": [
                {"test": "datum === tooltip", "value": 0},
                {"value": 1}
            ]
            }
        }
        }
    ]
    }
