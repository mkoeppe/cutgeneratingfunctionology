from __future__ import division, print_function, absolute_import

from cutgeneratingfunctionology.spam.basic_semialgebraic import BasicSemialgebraicSet_base
from itertools import chain

class BasicSemialgebraicSet_formal_relint(BasicSemialgebraicSet_base):

    r"""
    Represent the formal relative interior (see method ``formal_relint``) of
    another basic semialgebraic set ``upstairs_bsa``.
    """

    def __init__(self, upstairs_bsa):
        base_ring = upstairs_bsa.base_ring()
        ambient_dim = upstairs_bsa.ambient_dim()
        poly_ring = upstairs_bsa.poly_ring()
        super(BasicSemialgebraicSet_formal_relint, self).__init__(base_ring, ambient_dim, poly_ring=poly_ring)
        self._upstairs_bsa = upstairs_bsa

    def _repr_(self):
        return 'BasicSemialgebraicSet_formal_relint({})'.format(self.upstairs())

    def upstairs(self):
        return self._upstairs_bsa

    def eq_poly(self):
        return self.upstairs().eq_poly()

    def le_poly(self):
        return []

    def lt_poly(self):
        return chain(self.upstairs().lt_poly(), self.upstairs().le_poly())
