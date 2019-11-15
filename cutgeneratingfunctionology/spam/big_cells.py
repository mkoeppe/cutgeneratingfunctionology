r"""
Utilities for programming with big cells -- public interface

This module overrides globals such as ``min`` and ``sorted``.
"""

from __future__ import division, print_function, absolute_import

from cutgeneratingfunctionology.spam.big_cells_impl import *
from cutgeneratingfunctionology.spam.big_cells_impl import big_cells_min as min
from cutgeneratingfunctionology.spam.big_cells_impl import big_cells_sorted as sorted
