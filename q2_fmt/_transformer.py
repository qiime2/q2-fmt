# ----------------------------------------------------------------------------
# Copyright (c) 2022-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

import skbio
from q2_types.distance_matrix import LSMatFormat
from q2_fmt.plugin_setup import plugin


@plugin.register_transformer
def _1(ff: LSMatFormat) -> pd.Series:
    dm = skbio.DistanceMatrix.read(str(ff), format='lsmat', verify=False)
    return dm.to_series()
