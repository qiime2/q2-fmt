# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import pandas as pd

from qiime2.plugin.testing import TestPluginBase

from q2_fmt import tsv_to_dataframe

class TestTsvToDataframe(TestPluginBase):
    package = 'q2_fmt.tests'

    def test_tsv_to_dataframe_method(self):
        exp = pd.DataFrame(
            [['L1S8', 'AGCTGACTAGTC', 'gut', '2008'],
            ['L1S57', 'ACACACTATGGC', 'gut', '2009']],
            columns=['sample-id', 'barcode-sequence', 'body-site', 'year'],
            dtype=object)

        obs = tsv_to_dataframe('data/test_metadata.tsv')

        self.assertEqual(exp, obs)
