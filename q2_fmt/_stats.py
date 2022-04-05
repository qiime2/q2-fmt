# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# TODO: Add in significance testing for group_timepoints()
def is_significant():
    pass

# C: T_n
def mann_whitney_u(dist_a, dist_b, group):
    pass

# T_0: T_n
# T_n-1: T_n
def wilcoxon_srt(time_dist, hypothesis):
    pass


# output for both:
# baseline
# a   b   statistic   p-val   q-val
# 0   1   ....
# 0   2   ....
# 0   3   ....
#
# slope
# 0   1   ...
# 1   2   ...
# 2   3   ...

# control
# controla   0
# controla   1
# controla   2
# controla   3
