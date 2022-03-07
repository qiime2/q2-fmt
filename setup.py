# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

import versioneer

setup(
    name='q2-fmt',
    version=versioneer.get_version(),
    packages=find_packages(),
    package_data={},
    author='Liz Gehret',
    author_email='elizabeth.gehret@nau.edu',
    description='QIIME 2 Plugin used for FMT analyses.',
    license='BSD-3-Clause',
    url='https://github.com/qiime2/q2-fmt',
    zip_safe=False,
    entry_points={
        'qiime2.plugins': ['q2-fmt=q2_fmt.plugin_setup:plugin']
    }
)
