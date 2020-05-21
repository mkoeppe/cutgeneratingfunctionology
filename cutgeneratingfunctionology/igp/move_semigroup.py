"""
Move semigroups
"""
from __future__ import division, print_function, absolute_import

from sage.structure.element import Element, MonoidElement
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.homset import Homset
from sage.structure.richcmp import richcmp, op_NE, op_EQ

from .fast_linear import fast_linear_function
from .fast_piecewise import PiecewiseLinearFunction_1d
from .intervals import *

class FunctionalDirectedMove (MonoidElement):

    r"""
    Return a piecewise function to represent a functional directed move
    from a list of domain intervals and the functional directed move.
    """
    def __init__(self, domain_intervals, directed_move, parent=None):
        function = fast_linear_function(directed_move[0], directed_move[1])
        pieces = [ (interval, function) for interval in domain_intervals ]
        self._piecewise = piecewise = PiecewiseLinearFunction_1d(pieces)
        self.directed_move = tuple(directed_move)
        if parent is None:
            from sage.structure.element import get_coercion_model
            from sage.categories.rings import Rings
            cm = get_coercion_model()
            domain = cm.common_parent(directed_move[1], *piecewise.end_points())
            if domain in Rings():
                domain = domain.fraction_field()
            parent = FullMoveSemigroup(domain)
        super(FunctionalDirectedMove, self).__init__(parent)

    def __repr__(self):
        return "<FunctionalDirectedMove %s with domain %s, range %s>" % (self.directed_move, self.intervals(), self.range_intervals())

    def __call__(self, x):
        return self._piecewise(x)

    def __hash__(self):
        return id(self) # FIXME: should be hash((self.intervals(), self.directed_move))

    def _richcmp_(self, other, op):
        # FIXME: Implement restriction partial order
        return richcmp(self._piecewise, other._piecewise, op)

    def sign(self):
        r"""
        Return the sign of the move

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: h = FunctionalDirectedMove([[0.3, 0.4]], (1,0))
            sage: h.sign()
            1
        """
        return self.directed_move[0]

    def is_functional(self):
        return True

    def additive_faces(self, is_backward_translation=False):
        r"""
        Map ``FunctionalDirectedMove`` back to one-dimensional additive face(s) in the 2d-diagram.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: face_hor = Face([[2/5, 3/5],[4/5],[6/5,7/5]])
            sage: face_ver = Face([[4/5],[2/5, 3/5],[6/5,7/5]])
            sage: fdm_f = face_hor.functional_directed_move()
            sage: fdm_b = face_hor.functional_directed_move(is_backward_translation=True)
            sage: fdm_f.additive_faces() == [face_hor, face_ver]
            True
            sage: fdm_b.additive_faces(is_backward_translation=True) == [face_hor, face_ver]
            True
            sage: face_dia = Face([[2/5, 3/5],[2/5, 1/5],[4/5]])
            sage: fdm_d = face_dia.functional_directed_move()
            sage: fdm_d.additive_faces()
            [<Face ([1/5, 2/5], [2/5, 3/5], [4/5])>,
             <Face ([2/5, 3/5], [1/5, 2/5], [4/5])>]
        """
        from . import fractional, Face
        directed_move = self.directed_move
        domain = self.intervals()
        codomain = self.range_intervals()
        n = len(domain)
        faces = []
        if directed_move[0] == 1:
            if not is_backward_translation:
                # forward translation
                t = fractional(directed_move[1])
                for i in range(n):
                    a, b = domain[i][0], domain[i][1]
                    faces += [Face([[a, b],[t],[a+t, b+t]]),Face([[t],[a, b],[a+t, b+t]])]
            else:
                # backward translation
                t = fractional(-directed_move[1])
                for i in range(n):
                    a, b = codomain[i][0], codomain[i][1]
                    faces += [Face([[a, b],[t],[a+t, b+t]]),Face([[t],[a, b],[a+t, b+t]])]
        else:
            # reflection
            r = directed_move[1]
            for i in range(n):
                a, b = domain[i][0], domain[i][1]
                faces.append(Face([[a, b],[r-b, r-a],[r]]))
        return faces

    def __getitem__(self, item):
        return self.directed_move[item]

    def can_apply(self, x):
        r"""
        Determine if self can apply on `x`.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: h = FunctionalDirectedMove([[0.3, 0.4], [0.58, 0.68]], (1,0))
            sage: h.can_apply(0.3)
            True
            sage: h.can_apply(0.2)
            False
        """
        try:
            self(x)
            return True
        except ValueError:
            return False

    def apply_ignoring_domain(self, x):
        r"""
        Apply self on x by ignoring the domain (use modulo 1)

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: h = FunctionalDirectedMove([[0.3, 0.4], [0.58, 0.68]], (1,0))
            sage: h.apply_ignoring_domain(1/10)
            1/10
            sage: h = FunctionalDirectedMove([[0.1, 0.6]], (-1,1))
            sage: h.apply_ignoring_domain(1/2)
            1/2

        """
        from . import fractional, Face
        move_sign = self.sign()
        if move_sign == 1:
            next_x = fractional(x + self.directed_move[1])
        elif move_sign == -1:
            next_x = fractional(self.directed_move[1]-x)
        return next_x

    def apply_to_coho_interval(self, interval, inverse=False):
        # This does not do error checking.  Some code depends on this fact!
        # FIXME: This should be made clear in the name of this function.
        r"""
        Returns a range inverval from a given interval by applying the move.

        If the move sign is 1, the user can take the inverse of the operation,
        i.e `y = x - t_1`.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: h = FunctionalDirectedMove([[0.3, 0.4]], (-1, 1))
            sage: h.apply_to_coho_interval([1/10, 1/2])
            <Int[1/2, 9/10]>
            sage: h = FunctionalDirectedMove([[0.3, 0.4]], (1, 1/10))
            sage: h.apply_to_coho_interval([1/10, 1/2])
            <Int[1/5, 3/5]>
            sage: h.apply_to_coho_interval([1/10, 1/2], inverse=True)
            <Int[0, 2/5]>
        """
        if len(interval) <= 2:
            interval = coho_interval_from_interval(interval) # FIXME: Can be removed if FastPiecewise exclusively uses coho intervals.
        directed_move = self.directed_move
        move_sign = directed_move[0]
        if move_sign == 1:
            if inverse:
                result = closed_or_open_or_halfopen_interval(interval[0] - directed_move[1], interval[1] - directed_move[1], \
                                                             interval.left_closed, interval.right_closed)
            else:
                result = closed_or_open_or_halfopen_interval(interval[0] + directed_move[1], interval[1] + directed_move[1], \
                                                             interval.left_closed, interval.right_closed)
        elif move_sign == -1:
            result = closed_or_open_or_halfopen_interval(directed_move[1] - interval[1], directed_move[1] - interval[0], \
                                                         interval.right_closed, interval.left_closed)
        else:
            raise ValueError("Move not valid: %s" % list(move))
        return result

    def intervals(self):
        return self._piecewise.intervals()

    def range_intervals(self):
        r"""
        Return the range intervals of self.
        """
        return [ self.apply_to_coho_interval(interval) for interval in self.intervals() ]

    def is_identity(self):
        r"""
        Determine whether self is a (restricted) identity function or not.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: h = FunctionalDirectedMove([[0.3, 0.4]], (1, 0))
            sage: h.is_identity()
            True
            sage: h = FunctionalDirectedMove([[0.3, 0.4]], (-1, 1))
            sage: h.is_identity()
            False
        """
        return self.directed_move[0] == 1 and self.directed_move[1] == 0

    def restricting(self, components):
        r"""
        Returns a new move by removing ``self.restricted(component)`` for component in components.
        (The result may have the empty set as its domain.)
        """
        domain = self.intervals()                        # sorted.
        restricting_domain_list = []
        for component in components:
            preimages = [ self.apply_to_coho_interval(interval, inverse=True) for interval in component ]
            preimages.sort(key=coho_interval_left_endpoint_with_epsilon)
            restricting_domain_list.append(list(intersection_of_coho_intervals([component, preimages])))
        new_domain = union_of_coho_intervals_minus_union_of_coho_intervals([domain], restricting_domain_list, remove_closure=True )
        return FunctionalDirectedMove(new_domain, self.directed_move)

    def reduced_by_components(self, covered_components, pts_of_discontinuity=True):
        from . import interval_including_endpoints_if_continuous
        from . import equiv7_mode
        domain_move = [interval_including_endpoints_if_continuous(interval, pts_of_discontinuity, self) for interval in self.intervals()]
        covered_domains = []
        if equiv7_mode:
            for component in covered_components:
                preimages =  [ self.apply_to_coho_interval(interval, inverse=True) for interval in component ]
                preimages.sort(key=coho_interval_left_endpoint_with_epsilon)
                covered_domain = list(intersection_of_coho_intervals([component, preimages]))
                covered_domains.append(covered_domain) #list of open sets in the new version
            domain_component_open = union_of_coho_intervals_minus_union_of_coho_intervals(covered_domains,[])
            domain_component = union_of_coho_intervals_minus_union_of_coho_intervals([[interval_including_endpoints_if_continuous(interval, pts_of_discontinuity, self) for interval in domain_component_open]],[])
            domains = []
            for j in domain_component:
                keep = False
                for i in domain_move:
                    if coho_intervals_intersecting(i, j):
                        # extend to the boundary of covered interval.
                        keep = True
                        break
                if keep:
                    domains.append(j)
            extended_domains = union_of_coho_intervals_minus_union_of_coho_intervals([domain_move, domains], [])
            reduced_domains = []
            for i in extended_domains:
                i_open = open_interval(i[0], i[1])
                keep = True
                for j in domain_component_open:
                    if coho_interval_contained_in_coho_interval(i_open, j):
                        # remove interval i from domain if it is fully covered.
                        keep = False
                        break
                if keep:
                    reduced_domains.append(i)
            open_domains = [open_interval(i[0],i[1]) for i in reduced_domains]
        else:
            for component in covered_components:
                preimages =  [ self.apply_to_coho_interval(interval, inverse=True) for interval in component ]
                preimages.sort(key=coho_interval_left_endpoint_with_epsilon)
                restricted_domain = intersection_of_coho_intervals([component, preimages])
                covered_domain = [interval_including_endpoints_if_continuous(interval, pts_of_discontinuity, self) for interval in restricted_domain]
                covered_domains.append(covered_domain)
            domain_component = union_of_coho_intervals_minus_union_of_coho_intervals(covered_domains,[])
            domain_component_indices = set([])
            domain_move_indices = set(range(len(domain_move)))
            domains = []
            for index_i in list(domain_move_indices):
                i = domain_move[index_i]
                for index_j in range(len(domain_component)):
                    j = domain_component[index_j]
                    if coho_interval_contained_in_coho_interval(i, j):
                        # remove interval i from domain if it is fully covered.
                        domain_move_indices.remove(index_i)
                    elif coho_intervals_intersecting(i, j):
                        # extend to the boundary of covered interval.
                        domain_component_indices.add(index_j)
            domains = [[domain_move[index_i] for index_i in sorted(domain_move_indices)], [domain_component[index_j] for index_j in sorted(domain_component_indices)]]
            union_domains = union_of_coho_intervals_minus_union_of_coho_intervals(domains, [])
            open_domains = [open_interval(i[0],i[1]) for i in union_domains]
        return FunctionalDirectedMove(open_domains, self.directed_move)

    def plot(self, rgbcolor=None, color=None, *args, **kwds):
        from copy import copy
        kwds = copy(kwds)
        kwds['aspect_ratio'] = 1.0
        # ignore discontinuity markers in the moves diagram
        from cutgeneratingfunctionology.igp import show_moves_with_discontinuity_markers, show_translations_and_reflections_by_color
        kwds['discontinuity_markers'] = show_moves_with_discontinuity_markers
        if rgbcolor is None:
            if show_translations_and_reflections_by_color and self.sign() == -1:
                rgbcolor='red'
            else:
                rgbcolor='blue'
        return self._piecewise.plot(self, color=rgbcolor, *args, **kwds)

    def __invert__(self):
        r"""
        Returns the (pseudo-)inverse of ``self``.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: h = FunctionalDirectedMove([[3/5, 4/5]], (1, 1/5))
            sage: ~h
            <FunctionalDirectedMove (1, -1/5) with domain [(4/5, 1)], range [<Int[3/5, 4/5]>]>
            sage: ~~h == h
            True
            sage: h = FunctionalDirectedMove([[1/5, 2/5], [3/5, 4/5]], (-1, 1))
            sage: ~h == h
            True
            sage: h = FunctionalDirectedMove([[3/5, 4/5]], (-1, 1)); h
            <FunctionalDirectedMove (-1, 1) with domain [(3/5, 4/5)], range [<Int[1/5, 2/5]>]>
            sage: ~h
            <FunctionalDirectedMove (-1, 1) with domain [(1/5, 2/5)], range [<Int[3/5, 4/5]>]>
        """
        if self.sign() == 1:
            return FunctionalDirectedMove(self.range_intervals(), (1, -self[1]))
        else:
            return FunctionalDirectedMove(self.range_intervals(), self.directed_move)

    def _mul_(self, other):
        r"""
        Compute the directed move that corresponds to the directed move ``self`` after ``other``.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp.move_semigroup import FunctionalDirectedMove
            sage: FunctionalDirectedMove([(5/10,7/10)],(1, 2/10)) * FunctionalDirectedMove([(2/10,4/10)],(1,2/10))
            <FunctionalDirectedMove (1, 2/5) with domain [(3/10, 2/5)], range [<Int[7/10, 4/5]>]>

        """
        from . import compose_functional_directed_moves
        return compose_functional_directed_moves(self, other)

class FullMoveSemigroup(UniqueRepresentation, Homset):

    r"""
    The inverse semigroup `\Gamma^\subseteq(``domain``)` of restricted translations and reflections.

    Identity is the unrestricted identity translation.

    TESTS::

        sage: from cutgeneratingfunctionology.igp.move_semigroup import FullMoveSemigroup
        sage: TestSuite(FullMoveSemigroup(QQ)).run()
        sage: TestSuite(FullMoveSemigroup(QQ)).run(verbose=True, catch=False)  # not tested - for interactive use

    """

    Element = FunctionalDirectedMove

    # FIXME: Create category of inverse semigroups

    def __init__(self, domain):
        from sage.categories.sets_cat import SetsWithPartialMaps
        domain_cat = domain.category()
        if hasattr(domain_cat, "WithBasis"):   # Adapted from Homset.
            # The above is a lame but fast check that category is a
            # subcategory of Modules(...).
            base_ring = domain.base_ring()
        else:
            base_ring = domain
        cat = SetsWithPartialMaps()
        super(FullMoveSemigroup, self).__init__(domain, domain, cat, base_ring)

    def _element_constructor_(self, x, *args, **kwds):
        if isinstance(x, FunctionalDirectedMove):
            if args or kwds:
                raise TypeError("no extra arguments are allowed")
            return self.element_class(x.intervals(), x.directed_move, parent=self)
        return self.element_class(x, parent=self, *args, **kwds)

    def _an_element_(self):
        from sage.rings.rational_field import QQ
        return self.element_class([(0, QQ('1/2')), (QQ('1/2'), QQ('3/4'))],
                                  (1, QQ('1/4')), parent=self)

    def one(self):
        # FIXME: FastPiecewise cannot handle unbounded intervals.
        # So we have to use a restricted identity as the neutral element.
        return self.element_class([[0, 1]], (1, 0), parent=self)
