r"""
Basic semialgebraic sets using the Maple interface
"""

from __future__ import division, print_function, absolute_import

from cutgeneratingfunctionology.spam.basic_semialgebraic import BasicSemialgebraicSet_base

class BasicSemialgebraicSet_maple(BasicSemialgebraicSet_base):

    """
    Basic semialgebraic set using the Maple interface.

    Maple references:

     - https://www.maplesoft.com/support/help/Maple/view.aspx?path=SolveTools/SemiAlgebraic
     - https://www.maplesoft.com/support/help/Maple/view.aspx?path=RegularChains
     - https://www.maplesoft.com/support/help/Maple/view.aspx?path=RegularChains%2fSemiAlgebraicSetTools

    EXAMPLES::

        sage: import sage.interfaces.maple; from sage.interfaces.maple import Maple  # not tested - for interactive use
        sage: maple = sage.interfaces.maple.maple = Maple(server='logic.math.ucdavis.edu')   # not tested - for interactive use
        sage: from cutgeneratingfunctionology.spam.semialgebraic_maple import *
        sage: bsa = BasicSemialgebraicSet_maple(QQ, 2)    # optional - maple

    """

    def __init__(self, base_ring, ambient_dim, maple=None):
        if maple is None:
            from sage.interfaces.maple import maple
        self._maple = maple
        maple.with_package("SolveTools")
        ## self._semialgebraic = maple.SemiAlgebraic(sys, vars)
