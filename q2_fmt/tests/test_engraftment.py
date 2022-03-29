# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from qiime2.plugin.testing import TestPluginBase

from q2_fmt._format import TSVFileFormat

class TestEngraftment(TestPluginBase):
    pass
    # package = 'q2_fmt.tests'

    # def dataframe_adds_blank_column(self):
    #     exp = pd.DataFrame(
    #             [['L1S8', 'AGCTGACTAGTC', 'gut', '2008', ''],
    #             ['L1S57', 'ACACACTATGGC', 'gut', '2009', '']],
    #             columns=['sample-id', 'barcode-sequence', 'body-site', 'year', 'ziggy'],
    #             dtype=object)

    #     obs = (TSVFileFormat('data/test_metadata.tsv'))

    #     self.assertEqual(exp, obs)
