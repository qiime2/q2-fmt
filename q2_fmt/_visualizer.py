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

def hello_world(output_dir: str):
    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write('hello world')

    template = pkg_resources.resource_filename('q2_fmt', os.path.join('assets', 'index.html'))

def plot_rainclouds(output_dir: str):
    J_ENV = jinja2.Environment(
        loader=jinja2.PackageLoader('q2_fmt', 'assets')
    )

    index = J_ENV.get_template('index.html')

    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write(index.render(spec=json.dumps(_spec_())))

def _spec_():
    return {
    "$schema": "https://vega.github.io/schema/vega/v5.json",
    "description": "A basic bar chart example, with value labels shown upon mouse hover.",
    "width": 400,
    "height": 200,
    "padding": 5,

    "data": [
        {
        "name": "table",
        "values": [
            {"category": "A", "amount": 28},
            {"category": "B", "amount": 55},
            {"category": "C", "amount": 43},
            {"category": "D", "amount": 91},
            {"category": "E", "amount": 81},
            {"category": "F", "amount": 53},
            {"category": "G", "amount": 19},
            {"category": "H", "amount": 87}
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
        "domain": {"data": "table", "field": "category"},
        "range": "width",
        "padding": 0.05,
        "round": True
        },
        {
        "name": "yscale",
        "domain": {"data": "table", "field": "amount"},
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
        "from": {"data":"table"},
        "encode": {
            "enter": {
            "x": {"scale": "xscale", "field": "category"},
            "width": {"scale": "xscale", "band": 1},
            "y": {"scale": "yscale", "field": "amount"},
            "y2": {"scale": "yscale", "value": 0}
            },
            "update": {
            "fill": {"value": "steelblue"}
            },
            "hover": {
            "fill": {"value": "red"}
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
            "x": {"scale": "xscale", "signal": "tooltip.category", "band": 0.5},
            "y": {"scale": "yscale", "signal": "tooltip.amount", "offset": -2},
            "text": {"signal": "tooltip.amount"},
            "fillOpacity": [
                {"test": "datum === tooltip", "value": 0},
                {"value": 1}
            ]
            }
        }
        }
    ]
    }
