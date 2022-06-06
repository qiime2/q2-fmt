# ----------------------------------------------------------------------------
# Copyright (c) 2022-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import SemanticType

StatsTable = SemanticType('StatsTable', field_names=['kind'])

Pairwise = SemanticType('Pairwise', variant_of=StatsTable.field['kind'])

GroupDist = SemanticType('GroupDist', field_names=['order', 'dependence'])

Ordered = SemanticType('Ordered', variant_of=GroupDist.field['order'])
Unordered = SemanticType('Unordered', variant_of=GroupDist.field['order'])

Matched = SemanticType("Matched", variant_of=GroupDist.field['dependence'])
Independent = SemanticType("Independent",
                           variant_of=GroupDist.field['dependence'])
