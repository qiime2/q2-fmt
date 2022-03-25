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

GroupDist = SemanticType('GroupDist', field_names='content')

Gordinal =  SemanticType('Gordinal', variant_of=GroupDist.field['content'])
Gnominal = SemanticType('Gnominal', variant_of=GroupDist.field['content'])
