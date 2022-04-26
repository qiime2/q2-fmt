# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources

def hello_world(output_dir: str):
    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write('hello world')

    template = pkg_resources.resource_filename('q2_fmt', os.path.join('assets', 'index.html.jinja2'))

def plot_rainclouds(output_dir: str):
    pass
