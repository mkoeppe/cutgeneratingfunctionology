from __future__ import division, print_function, absolute_import

from cutgeneratingfunctionology.spam.semialgebraic_predicate import BasicSemialgebraicSet_predicate
from cutgeneratingfunctionology.spam.basic_semialgebraic import BasicSemialgebraicSet_base

class BasicSemialgebraicSet_local(BasicSemialgebraicSet_predicate):
    """
    EXAMPLES::

        sage: from cutgeneratingfunctionology.spam.basic_semialgebraic_local import BasicSemialgebraicSet_local
        sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import BasicSemialgebraicSet_eq_lt_le_sets, BasicSemialgebraicSet_veronese
        sage: import logging; logging.disable(logging.INFO)             # Suppress output in automatic tests.

    This class does not do any computation and does not implement the methods
    ``eq_poly`` etc.  Only when a new basic semialgebraic set is requested to
    be created using ``from_bsa(self)``, inequalities valid for a "local patch"
    (heuristic component) of the basic semialgebraic set ``self.upstairs()``
    are computed.

    In the following example, ``BasicSemialgebraicSet_eq_lt_le_sets.from_bsa(self)``
    is called.  Because this class does not remove any redundancy, the
    resulting inequality description in this example is non-minimal::

        sage: test_point = [1, 1/2]
        sage: P.<x,y> = QQ[]
        sage: polynomial_map = [2*y, y]
        sage: bsa = BasicSemialgebraicSet_eq_lt_le_sets(eq=[-x + 2*y], lt=[3*x - 4, y^2 - x], le=[]).section(polynomial_map)
        sage: sorted(bsa.lt_poly())
        [6*y - 4, y^2 - 2*y]
        sage: bsa_local = BasicSemialgebraicSet_local(bsa, test_point)
        sage: sorted(BasicSemialgebraicSet_eq_lt_le_sets.from_bsa(bsa_local).lt_poly())
        [-y, y - 2, 6*y - 4]
        sage: bsa_local_veronese = BasicSemialgebraicSet_veronese.from_bsa(bsa_local)
        sage: sorted(bsa_local_veronese.lt_poly())
        [-y, 3*y - 2]

    Interesting examples: {y = x^2 and 36*x^2+y^2-10*x*y-50*x-y+24<=0} has two components.
    The example has been chosen carefully so that it factors over QQ to avoid implementation restrictions::

        sage: P.<x,y> = QQ[]
        sage: pt1 = [3/2, 9/4]; pt2 = [7/2, 49/4]
        sage: bsa = BasicSemialgebraicSet_eq_lt_le_sets(eq=[y - x^2], lt=[], le=[36*x^2+y^2-10*x*y-50*x-y+24]).section([x, x^2])
        sage: sorted(bsa.le_poly())
        [x^4 - 10*x^3 + 35*x^2 - 50*x + 24]
        sage: bsa_local_1 = BasicSemialgebraicSet_local(bsa, pt1)
        sage: sorted(BasicSemialgebraicSet_eq_lt_le_sets.from_bsa(bsa_local_1).le_poly())
        [-x + 1, x - 4, x - 3, x - 2]
        sage: bsa_local_veronese_1 = BasicSemialgebraicSet_veronese.from_bsa(bsa_local_1)
        sage: sorted(bsa_local_veronese_1.le_poly())
        [-x + 1, x - 2]
        sage: bsa_local_2 = BasicSemialgebraicSet_local(bsa, pt2)
        sage: bsa_local_veronese_2 = BasicSemialgebraicSet_veronese.from_bsa(bsa_local_2)
        sage: sorted(bsa_local_veronese_2.le_poly())
        [-x + 3, x - 4]

    Another example: {y = x^2 and 1/9*(x-5)^2 + (y-10)^2 <= 25}, which has two components, but it does not factor over QQ. We have to use AA, but then we cannot use PPL to implement the upstairs BSA.  We resort to upstairs_bsa_class='mip'.  However, that offers only nonstrict inequalities and does not guarantee minimal representations:  Whether redundant inequalities are removed depends on the insertion order. Simplification in bsa_local_veronese_1.le_poly() is as expected. No simplification is done in bsa_local_veronese_2.le_poly()::

        sage: P.<x,y> = AA[] # not using QQ because the polynomial does not factor over QQ>
        sage: pt1 = [3, 9]; pt2 = [-3, 9]
        sage: bsa = BasicSemialgebraicSet_eq_lt_le_sets(eq=[y - x^2], lt=[], le=[1/9*(x-5)^2 + (y-10)^2 - 25]).section([x, x^2])
        sage: sorted(bsa.le_poly())[0].factor()
        (x - 3.871152417778312?) * (x - 2.254871046366155?) * (x + 2.376196051287999?) * (x + 3.749827412856468?)
        sage: bsa_local_1 = BasicSemialgebraicSet_local(bsa, pt1)
        sage: sorted(BasicSemialgebraicSet_eq_lt_le_sets.from_bsa(bsa_local_1).le_poly())
        [-x - 3.749827412856468?,
         -x - 2.376196051287999?,
         -x + 2.254871046366155?,
         x - 3.871152417778312?]
        sage: bsa_local_veronese_1 = BasicSemialgebraicSet_veronese.from_bsa(bsa_local_1, upstairs_bsa_class='mip')
        sage: sorted(bsa_local_veronese_1.le_poly())
        [-x + 2.254871046366154, x - 3.871152417778312]
        sage: bsa_local_2 = BasicSemialgebraicSet_local(bsa, pt2)
        sage: bsa_local_veronese_2 = BasicSemialgebraicSet_veronese.from_bsa(bsa_local_2, upstairs_bsa_class='mip')
        sage: sorted(bsa_local_veronese_2.le_poly()) # expect [-x - 3.749827412856468, x + 2.376196051287999]
        [-x - 3.749827412856468,
        x - 3.871152417778312,
        x - 2.254871046366154,
        x + 2.376196051287999]
    """
    def __init__(self, bsa, test_point, base_ring=None, ambient_dim=None, poly_ring=None):
        assert (test_point in bsa)
        if base_ring is None:
            base_ring = bsa.base_ring()
        if ambient_dim is None:
            ambient_dim = bsa.ambient_dim()
        if poly_ring is None:
            poly_ring = bsa.poly_ring()
        super(BasicSemialgebraicSet_local, self).__init__(bsa.__contains__, true_point=test_point, base_ring=base_ring, ambient_dim=ambient_dim, poly_ring=poly_ring, allow_refinement=True)

    def test_point(self):
        return self._true_point

    @classmethod
    def from_bsa(cls, bsa, test_point):
        if bsa.__class__ == cls and bsa.test_point() == test_point:
            return bsa
        return cls(bsa, test_point, poly_ring=bsa.poly_ring())

class BasicSemialgebraicSet_local_with_test_point(BasicSemialgebraicSet_base):
    """
    A BasicSemialgebraicSet_local_with_test_point has an "upstairs" to which it records the factored inequalities, and a testpoint. That upstairs can be any bsa, and would be a veronese instance for implementing K. We factor polynomial using testpoint in add_polynomial_constraint of the new class, where it adds the factors to the upstairs. Its lt_poly etc. delegate to upstairs. Later, create a new K._local_bsa from test point and set its upstairs =  K._factor_bsa. Simplifying the description with the new local BSA: you just use from_bsa to make the local BSA and then you query it.
    """
    def __init__(self, upstairs_bsa, test_point, base_ring=None, ambient_dim=None, poly_ring=None):
        assert (test_point in bsa)
        if base_ring is None:
            base_ring = bsa.base_ring()
        if ambient_dim is None:
            ambient_dim = bsa.ambient_dim()
        if poly_ring is None:
            poly_ring = bsa.poly_ring()
        super(BasicSemialgebraicSet_local, self).__init__(bsa.__contains__, true_point=test_point, base_ring=base_ring, ambient_dim=ambient_dim, poly_ring=poly_ring, allow_refinement=True)

    def _repr_(self):
        return 'BasicSemialgebraicSet_local_with_test_point(eq={}, lt={}, le={})'.format(list(self.eq_poly()), list(self.lt_poly()), list(self.le_poly()))

    def upstairs(self):
        return self._upstairs_bsa

    def test_point(self):
        return self._true_point

    @classmethod
    def from_bsa(cls, bsa, test_point):
        if bsa.__class__ == cls and bsa.test_point() == test_point:
            return bsa
        return cls(bsa, test_point, poly_ring=bsa.poly_ring())

    def eq_poly(self):
        for l in self.groebner_basis():
            if not (l in self.base_ring() and l == 0):
                yield l

    def le_poly(self):
        for l in self.upstairs().le_poly():
            l_red = l.reduce(self.groebner_basis())
            if not (l_red in self.base_ring() and l_red <= 0):
                yield l_red

    def lt_poly(self):
        for l in self.upstairs().lt_poly():
            l_red = l.reduce(self.groebner_basis())
            if not (l_red in self.base_ring() and l_red < 0):
                yield l_red

    def add_polynomial_constraint(self, lhs, op):
        self.upstairs().add_polynomial_constraint(lhs, op)
        if op == operator.eq:
            # FIXME: make it incremental
            #self._ideal = self.poly_ring().ideal(self.upstairs().eq_poly())
            gens = list(self.groebner_basis()) + [l.reduce(self.groebner_basis()) for l in self.upstairs().eq_poly()]
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
        return self.upstairs().is_polynomial_constraint_valid()

    @classmethod
    def from_bsa(cls, bsa, base_ring=None, poly_ring=None):
        pass
