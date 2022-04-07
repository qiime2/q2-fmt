# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import SemanticType
from q2_types.sample_data import SampleData

ModelTests = SemanticType('ModelTests', variant_of=SampleData.field['type'])

GroupDist = SemanticType('GroupDist', field_names=['order', 'dependence'])

Ordered = SemanticType('Ordered', variant_of=GroupDist.field['order'])
Unordered = SemanticType('Unordered', variant_of=GroupDist.field['order'])

Matched = SemanticType("Matched", variant_of=GroupDist.field['dependence'])
Independent = SemanticType("Independent",
                           variant_of=GroupDist.field['dependence'])


GroupDist[Ordered, Matched]  # <- Wilcoxon
GroupDist[Ordered, Independent]  # <- U
# GroupDist[Unordered, Matched]  # <- Wilcoxon (+ a group to compare against)
GroupDist[Unordered, Independent]  # <- U
