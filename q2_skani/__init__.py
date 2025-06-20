# flake8: noqa
# ----------------------------------------------------------------------------
# Copyright (c) 2025, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from .skani import compare_seqs

__all__ = ["compare_seqs"]

try:
    from ._version import __version__
except ModuleNotFoundError:
    __version__ = "0.0.0+notfound"
