from __future__ import division, print_function, absolute_import

from cutgeneratingfunctionology.spam.basic_semialgebraic import BasicSemialgebraicSet_base
from itertools import chain

class BasicSemialgebraicSet_intersection(BasicSemialgebraicSet_base):

    r"""
    Represent the intersection of finitely many basic semialgebraic sets.
    See method ``intersection``.
    """

    def __init__(self, bsa_list):
        bsa_list = list(bsa_list)
        base_ring = bsa_list[0].base_ring()
        ambient_dim = bsa_list[0].ambient_dim()
        poly_ring = bsa_list[0].poly_ring()
        if not all(bsa.base_ring() == base_ring
                   and bsa.ambient_dim() == ambient_dim
                   and bsa.poly_ring() == poly_ring
                   for bsa in bsa_list):
            raise ValueError("all sets in the intersection must have the same base_ring, ambient_dim, poly_ring")
        super(BasicSemialgebraicSet_intersection, self).__init__(base_ring, ambient_dim, poly_ring=poly_ring)
        self._bsa_list = bsa_list

    def _repr_(self):
        return 'BasicSemialgebraicSet_intersection({})'.format(self._bsa_list)

    def eq_poly(self):
        return set(chain(*[bsa.eq_poly() for bsa in self._bsa_list]))

    def le_poly(self):
        return set(chain(*[bsa.le_poly() for bsa in self._bsa_list]))

    def lt_poly(self):
        return set(chain(*[bsa.lt_poly() for bsa in self._bsa_list]))

    def linear_function_upper_bound(self, form):
        return min(bsa.linear_function_upper_bound(form) for bsa in self._bsa_list)

    def is_polynomial_constraint_valid(self, lhs, op):
        return any(bsa.is_polynomial_constraint_valid(lhs, op) for bsa in self._bsa_list)
