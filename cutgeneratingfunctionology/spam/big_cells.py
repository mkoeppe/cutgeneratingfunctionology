r"""
Utilities for programming with big cells -- public interface

Importing symbols from this module overrides globals such as ``min`` and ``sorted``.
"""

from __future__ import division, print_function, absolute_import

from cutgeneratingfunctionology.spam.big_cells_impl import *
from cutgeneratingfunctionology.spam.big_cells_impl import big_cells_min as min
from cutgeneratingfunctionology.spam.big_cells_impl import big_cells_sorted as sorted

def bigcellify_module(module):
    """
    Modify ``module`` by overriding some of its globals by big cells functions.

    WARNING:

    This is, of course, dangerous.

    EXAMPLES::

        sage: print("Danger ahead!"); from six.moves import reload_module; from cutgeneratingfunctionology import *; from cutgeneratingfunctionology.spam import *; reload_module(big_cells); reload_module(igp); big_cells.bigcellify_module(igp); from cutgeneratingfunctionology.igp import *
        Danger ahead!...
        sage: igp.sorted
        <...big_cells_sorted...>
        sage: sorted
        <...big_cells_sorted...>
    """
    module.min = min
    module.sorted = sorted
