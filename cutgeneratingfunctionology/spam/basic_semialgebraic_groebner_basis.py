from __future__ import division, print_function, absolute_import

from cutgeneratingfunctionology.spam.basic_semialgebraic import BasicSemialgebraicSet_base
from itertools import chain
from sage.rings.ideal import Ideal
import operator

class BasicSemialgebraicSet_groebner_basis(BasicSemialgebraicSet_base):

    r"""
    Represent the equations of ``upstairs_bsa`` by its reduced Groebner basis B, and the inequalites by the reduced polynomials modulo B.

    EXAMPLES::

        sage: Q.<x,y,z> = QQ[]
        sage: veronese = BasicSemialgebraicSet_veronese(poly_ring=Q)
        sage: veronese.add_polynomial_constraint(x^2 + y*z - 11, operator.lt)
        sage: veronese.add_polynomial_constraint(2*y - z, operator.eq)
        sage: veronese.add_polynomial_constraint(y + z + 3, operator.eq)
        sage: bsa_eq_lt_le = BasicSemialgebraicSet_eq_lt_le_sets.from_bsa(veronese); bsa_eq_lt_le
        BasicSemialgebraicSet_eq_lt_le_sets(eq=[z + 2, y + 1], lt=[x^2 + y*z - 11], le=[])
        sage: bsa_gb = BasicSemialgebraicSet_groebner_basis(veronese); bsa_gb
        BasicSemialgebraicSet_groebner_basis(eq=[y + 1, z + 2], lt=[x^2 - 9], le=[])
        sage: bsa_gb.add_polynomial_constraint(x-2, operator.eq); bsa_gb
        BasicSemialgebraicSet_groebner_basis(eq=[x - 2, y + 1, z + 2], lt=[], le=[])
        sage: list(bsa_gb.upstairs().lt_poly())
        [x^2 + y*z - 11]
        sage: list(bsa_gb.upstairs().eq_poly())
        [x - 2, z + 2, y + 1]
        sage: bsa_gb1 = BasicSemialgebraicSet_groebner_basis.from_bsa(bsa_eq_lt_le, poly_ring=bsa_eq_lt_le.poly_ring().change_ring(order='invlex')); bsa_gb1
        BasicSemialgebraicSet_groebner_basis(eq=[z + 2, y + 1], lt=[x^2 - 9], le=[])
        sage: bsa_gb2 = BasicSemialgebraicSet_groebner_basis(bsa_eq_lt_le); bsa_gb2
        BasicSemialgebraicSet_groebner_basis(eq=[y + 1, z + 2], lt=[x^2 - 9], le=[])
        sage: bsa_gb2.add_polynomial_constraint(x^3-8, operator.eq); bsa_gb2
        BasicSemialgebraicSet_groebner_basis(eq=[x^3 - 8, y + 1, z + 2], lt=[x^2 - 9], le=[]
        sage: bsa_gb2.add_polynomial_constraint(x^5-32, operator.eq); bsa_gb2
        BasicSemialgebraicSet_groebner_basis(eq=[x - 2, y + 1, z + 2], lt=[], le=[])
        sage: bsa_gb2.add_polynomial_constraint(x^2-1, operator.lt)
        sage: bsa_gb2.is_empty()
        True
        sage: bsa_gb2.upstairs().is_empty()
        Traceback (most recent call last):
        ...
        NotImplementedError...
    """

    def __init__(self, upstairs_bsa, base_ring=None, poly_ring=None):
        if base_ring is None:
            base_ring = upstairs_bsa.base_ring()
        ambient_dim = upstairs_bsa.ambient_dim()
        if poly_ring is None:
            poly_ring = upstairs_bsa.poly_ring()
        super(BasicSemialgebraicSet_groebner_basis, self).__init__(base_ring, ambient_dim, poly_ring=poly_ring)
        self._upstairs_bsa = upstairs_bsa
        self._ideal = poly_ring.ideal([poly_ring(l) for l in self.upstairs().eq_poly()])
        self._gb = self._ideal.groebner_basis('libsingular:groebner')
        #FIXME: should we take the radical ideal and its gb?

    def _repr_(self):
        return 'BasicSemialgebraicSet_groebner_basis(eq={}, lt={}, le={})'.format(list(self.eq_poly()), list(self.lt_poly()), list(self.le_poly()))

    def upstairs(self):
        return self._upstairs_bsa

    def eq_poly(self):
        for l in self.groebner_basis():
            if not (l in self.base_ring() and l == 0):
                yield l

    def le_poly(self):
        for l in self.upstairs().le_poly():
            l_red = self.poly_ring()(l).reduce(self.groebner_basis())
            if not (l_red in self.base_ring() and l_red <= 0):
                yield l_red

    def lt_poly(self):
        for l in self.upstairs().lt_poly():
            l_red = self.poly_ring()(l).reduce(self.groebner_basis())
            if not (l_red in self.base_ring() and l_red < 0):
                yield l_red

    def add_polynomial_constraint(self, lhs, op):
        self.upstairs().add_polynomial_constraint(lhs, op)
        if op == operator.eq:
            # FIXME: make it incremental
            #self._ideal = self.poly_ring().ideal(self.upstairs().eq_poly())
            gens = list(self.groebner_basis()) + [self.poly_ring()(l).reduce(self.groebner_basis()) for l in self.upstairs().eq_poly()]
            self._ideal = self.poly_ring().ideal(gens)
            self._gb = self._ideal.groebner_basis('libsingular:groebner')

    def is_polynomial_constraint_valid(self, lhs, op):
        """lazy implementation. First try the same method as in class BasicSemialgebraicSet_eq_lt_le_sets, then try the method of self.upstairs(). """
        if op == operator.lt:
            if lhs in self.lt_poly():
                return True
        elif op == operator.gt:
            if -lhs in self.lt_poly():
                return True
        elif op == operator.eq:
            if (lhs in self.eq_poly()) or (-lhs in self.eq_poly()):
                return True
        elif op == operator.le:
            if (lhs in self.le_poly()) or (lhs in self.lt_poly()):
                return True
        elif op == operator.ge:
            if (-lhs in self.le_poly()) or (-lhs in self.lt_poly()):
                return True
        else:
            raise ValueError("{} is not a supported operator".format(op))
        return self.upstairs().is_polynomial_constraint_valid(lhs, op)

    def ideal(self):
        return self._ideal

    def groebner_basis(self):
        return self._gb

    @classmethod
    def from_bsa(cls, bsa, base_ring=None, poly_ring=None):
        if bsa.__class__ == cls and (base_ring is None or base_ring == bsa.base_ring()) and (poly_ring is None or poly_ring == bsa.poly_ring()):
            return bsa
        return cls(bsa, base_ring=base_ring, poly_ring=poly_ring)
