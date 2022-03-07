# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import Plugin

import q2_fmt

# from q2_fmt.format_types import (
    # Formats go here )

plugin = Plugin(name='fmt',
                version=q2_fmt.__version__,
                website='https://github.com/qiime2/q2-fmt',
                package='q2_diversity',
                description='This QIIME 2 plugin supports FMT analyses.',
                short_description='Plugin for analyzing FMT data.')

# plugin.register_formats(
#     # Formats go here
# )

# plugin.register_semantic_types(
#     # Semantic types go here
# )
# plugin.register_semantic_type_to_format(
#     # Semantic types and formats go here
# )
