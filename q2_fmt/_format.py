# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import model


class TSVFileFormat(model.TextFileFormat):
    """Format for TSV file.

    More to be added on this later.

    """


TSVFileDirFmt = model.SingleFileDirectoryFormat(
    'TSVFileFormat', 'file.tsv', TSVFileFormat)
