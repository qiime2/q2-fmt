# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._version import get_versions
from ._format import TSVFileFormat, AnnotatedTSVDirFmt
from ._type import (ModelTests, GroupDist,
                    Ordered, Unordered, Matched, Independent)

__version__ = get_versions()['version']
del get_versions

__all__ = ['TSVFileFormat', 'AnnotatedTSVDirFmt', 'ModelTests',
           'GroupDist', 'Ordered', 'Unordered', 'Matched', 'Independent']
