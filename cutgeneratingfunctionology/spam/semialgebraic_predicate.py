r"""
Semialgebraic and basic semialgebraic sets defined by a Python predicate

Converted to explicit form when requested, using function tracing.

"""

from __future__ import division, print_function, absolute_import

from cutgeneratingfunctionology.spam.basic_semialgebraic import BasicSemialgebraicSet_base

class BasicSemialgebraicSet_predicate(BasicSemialgebraicSet_base):

    """
    A basic semialgebraic set defined by a Python predicate.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.spam.semialgebraic_predicate import BasicSemialgebraicSet_predicate
        sage: def is_in_unit_disk(x): return x[0]^2 + x[1]^2 <= 1
        sage: bsa = BasicSemialgebraicSet_predicate(is_in_unit_disk)
    """

    def __init__(self, predicate, true_point=None, base_ring=None, ambient_dim=None, poly_ring=None):
        r"""
        Construct a basic semialgebraic set of points where ``predicate`` is true.

        If ``true_point`` is not ``None``, it must be a point where ``predicate`` is true.
        """
        # FIXME: Normalize predicate's inputs
        super(BasicSemialgebraicSet_predicate, self).__init__(base_ring=base_ring, ambient_dim=ambient_dim, poly_ring=poly_ring)
        self._predicate = predicate
        self._true_point = true_point

    def __contains__(self, point):
        r"""
        Whether the set contains the ``point`` (vector).
        """
        return self._predicate(point)

    def find_point(self, point):
        r"""
        Find a point in ``self``.
        """
        if self._true_point is not None and self._true_point in self:
            return self._true_point
        return super(BasicSemialgebraicSet_predicate, self).find_point()

    @classmethod
    def from_bsa(cls, bsa):
        if bsa.__class__ == cls:
            return bsa
        return cls(bsa.__contains__, poly_ring=bsa.poly_ring())

    def add_constraints_to(self, bsa):
        from cutgeneratingfunctionology.igp import ParametricRealField
        K = ParametricRealField(values=self._true_point, bsa=bsa,
                                allow_refinement=False, big_cells=True, mutable_values=False)
        assert self._predicate(K.gens())
