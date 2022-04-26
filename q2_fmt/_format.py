# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import model

import pandas as pd


class RecordTSVFileFormat(model.TextFileFormat):
    """Format for TSV file.

    first line is headers


    More to be added on this later.

    """
    def _validate_(self, level):
        pass


class AnnotatedTSVDirFmt(model.DirectoryFormat):
    data = model.File('data.tsv', format=RecordTSVFileFormat)
    metadata = model.File('metadata.tsv', format=RecordTSVFileFormat)

    def _validate_(self, level='min'):
        data = self.data.view(pd.DataFrame)
        metadata = self.metadata.view(pd.DataFrame)

        if list(data.columns) != list(metadata['column']):
            raise model.ValidationError(
                'The metadata TSV does not completely describe the data TSV'
                ' columns.')

        if metadata.columns[0] != 'column':
            raise model.ValidationError('The metadata TSV does not start with'
                                        ' "column" on the header line.')
