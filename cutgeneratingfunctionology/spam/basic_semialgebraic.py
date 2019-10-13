r"""
Mutable basic semialgebraic sets
"""

from __future__ import division, print_function, absolute_import

from sage.structure.element import Element
from sage.modules.free_module_element import vector
from sage.rings.all import QQ, ZZ
from sage.rings.real_double import RDF
from sage.misc.abstract_method import abstract_method

from copy import copy
import operator

try:
    from ppl import Variable, Constraint, Linear_Expression, Constraint_System, NNC_Polyhedron, Poly_Con_Relation, Poly_Gen_Relation, Generator, MIP_Problem, point as ppl_point
except ImportError:
    # old Sage
    from sage.libs.ppl import Variable, Constraint, Linear_Expression, Constraint_System, NNC_Polyhedron, Poly_Con_Relation, Poly_Gen_Relation, Generator, MIP_Problem, point as ppl_point

from sage.numerical.mip import MixedIntegerLinearProgram, MIPVariable, MIPSolverException

class BasicSemialgebraicSet_base(Element):

    """
    Abstract base class of mutable basic semialgebraic sets.
    """

    def __init__(self, ambient_dim):
        """
        Initialize a basic semialgebraic set as the universe in ``ambient_dim``.
        """
        self._ambient_dim = ambient_dim

    def ambient_dim(self):
        return self._ambient_dim

    @abstract_method
    def __contains__(self, point):
        """
        Whether the set contains the ``point`` (vector).
        """

    @abstract_method
    def closure_polyhedron(self):
        """
        If the topological closure is a polyhedron, return it in
        the form of a Sage polyhedron.  Otherwise raise an error.
        """

    def add_space_dimensions_and_embed(self, space_dim_to_add):
        """
        Mutate ``self`` by injecting it into a higher dimensional space.
        """
        self._ambient_dim += space_dim_to_add

    @abstract_method
    def add_linear_constraint(self, lhs, op):
        """
        ``lhs`` should be a vector of length ambient_dim.
        Add the constraint ``lhs`` * x ``op`` 0,
        where ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.
        """

    @abstract_method
    def is_linear_constraint_valid(self, lhs, op):
        """
        Whether the constraint ``lhs`` * x ``op`` 0
        is satisfied for all points of ``self``.
        """

    @abstract_method
    def find_point(self):
        """
        Find a point in ``self``.
        """
        # default implementation could go through self.closure_polyhedron()

## (1) In the first step, we implement the following class.  Everything is linear.
## Rewrite all direct uses of PPL in ParametricRealFieldElement, ParametricRealField
## using method calls to this class.

poly_is_included = Poly_Con_Relation.is_included()
#strictly_intersects = Poly_Con_Relation.strictly_intersects()
point_is_included = Poly_Gen_Relation.subsumes()
#con_saturates = Poly_Con_Relation.saturates()

from sage.arith.functions import lcm

class BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(BasicSemialgebraicSet_base):

    """
    A (possibly half-open) polyhedral basic semialgebraic set,
    represented by a PPL ``NNC_Polyhedron``

    """

    def __init__(self, ambient_dim=None, polyhedron=None):
        r"""
        Initialize a basic semialgebraic set as the universe in
        ``ambient_dim``, or, if ``polyhedron`` (an ``NNC_Polyhedron``,
        which after that belongs to this object) is provided, as
        that.
        """
        if ambient_dim is None:
            if polyhedron is None:
                raise ValueError("at least one of ambient_dim and polyhedron must be provided")
            ambient_dim = polyhedron.space_dimension()
        super(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron, self).__init__(ambient_dim)
        if polyhedron is None:
            self._polyhedron = NNC_Polyhedron(ambient_dim, 'universe')
        else:
            self._polyhedron = polyhedron

    def __copy__(self):
        r"""
        Make a copy of ``self``.

        TESTS:

        Test that it is actually making a copy of the (mutable!) NNC_Polyhedron::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(2)
            sage: P._polyhedron is copy(P)._polyhedron
            False
        """
        return self.__class__(polyhedron=copy(self._polyhedron))

    # override the abstract methods
    def __contains__(self, point):
        """
        Whether the set contains the ``point`` (vector).
        """
        rational_list = [ QQ(x) for x in point ]
        num_list = [x.numerator() for x in rational_list]
        den_list = [x.denominator() for x in rational_list]
        common_den = lcm(den_list)
        coef = [common_den // den_list[i] * num_list[i] for i in range(len(rational_list))]
        pt = ppl_point(Linear_Expression(coef, 0), common_den)
        return self._polyhedron.relation_with(pt).implies(point_is_included)

    def find_point(self):
        """
        Find a point in ``self``.
        """
        def to_point(g):
            den = g.divisor()
            return vector(QQ, ( QQ(x)/den for x in g.coefficients() ))
        def to_vector(g):
            return vector(QQ, ( QQ(x) for x in g.coefficients() ))
        points = [ to_point(g) for g in self._polyhedron.generators()
                   if g.is_point() or g.is_closure_point() ]
        rays = [ to_vector(g) for g in self._polyhedron.generators()
                 if g.is_ray() ]
        if points:
            p = sum(points) / len(points)
            if rays:
                p += sum(rays) / len(rays)
            return p
        raise NotImplementedError("find_test_point implementation cannot handle this case")

    def add_space_dimensions_and_embed(self, space_dim_to_add):
        """
        Mutate ``self`` by injecting it into a higher dimensional space.
        """
        super(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron, self).add_space_dimensions_and_embed(space_dim_to_add)
        self._polyhedron.add_space_dimensions_and_embed(space_dim_to_add)

    def is_linear_constraint_valid(self, linexpr, op):
        """
        Whether the constraint ``lhs`` * x ``op`` 0
        is satisfied for all points of ``self``.
        """
        # REFACTOR: Right now linexpr is assumed to be a PPL Linear_Expression. Change that
        if op == operator.lt:
            constraint_to_add = (linexpr < 0)
        elif op == operator.eq:
            constraint_to_add = (linexpr == 0)
        elif op == operator.le:
            constraint_to_add = (linexpr <= 0)
        else:
            raise ValueError("{} is not a supported operator".format(op))
        return self._polyhedron.relation_with(constraint_to_add).implies(poly_is_included)

    def add_linear_constraint(self, linexpr, op):
        """
        ``lhs`` should be a vector of length ambient_dim.
        Add the constraint ``lhs`` * x ``op`` 0,
        where ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.
        """
        if op == operator.lt:
            constraint_to_add = (linexpr < 0)
        elif op == operator.eq:
            constraint_to_add = (linexpr == 0)
        elif op == operator.le:
            constraint_to_add = (linexpr <= 0)
        else:
            raise ValueError("{} is not a supported operator".format(op))
        self._polyhedron.add_constraint(constraint_to_add)


## (2) Then
class BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram(BasicSemialgebraicSet_base):

    """
    A (possibly half-open) polyhedral basic semialgebraic set,
    represented by a PPL ``NNC_Polyhedron``
    """

    def __init__(self, ambient_dim, solver=None):
        """
        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram(3)

        """
        super(BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram, self).__init__(ambient_dim)
        self._mip = MixedIntegerLinearProgram(solver=solver)

    def __copy__(self):
        raise NotImplementedError()


## (3) Then introduce the following class to simplify the code in parametric.sage
class BasicSemialgebraicSet_eq_lt_le_sets(BasicSemialgebraicSet_base):

    """
    A basic semialgebraic set, represented in a straightforward way
    as 3 finite sets of polynomial constraints `p(x) OP 0`.
    """

    def __init__(self, ambient_dim=None, ring=None, eq=[], lt=[], le=[]):
        """
        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_eq_lt_le_sets(3)

        """
        if ring is None:
            polys = list(eq) + list(lt) + list(le)
            if polys:
                ring = polys[0].parent()
        if (ambient_dim is None) and (ring is not None):
            ambient_dim = ring.ngens()
        super(BasicSemialgebraicSet_eq_lt_le_sets, self).__init__(ambient_dim)
        self._ring = ring
        self._eq = set(eq)
        self._lt = set(lt)
        self._le = set(le)

    def __copy__(self):
        """
        Make a copy of ``self``.

        TESTS:

        Test that it is actually making a copy of the (mutable!) sets::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P = BasicSemialgebraicSet_eq_lt_le_sets(2)
            sage: P._eq is copy(P)._eq
            False
            sage: P._lt is copy(P)._lt
            False
            sage: P._le is copy(P)._le
            False
        """
        return self.__class__(self.ambient_dim(),
                              self._ring,
                              self._eq, self._lt, self._le)

    # override the abstract methods
    def __contains__(self, x):
        return all(f(x) == 0 for f in self._eq) and all(f(x) <= 0 for f in self._le) and all(f(x) < 0 for f in self._lt)

    def __repr__(self):
        return 'BasicSemialgebraicSet_eq_lt_le_sets(eq = {}, lt = {}, le = {})'.format(list(self._eq), list(self._lt), list(self._le))

    ### TODO: Add implementations of more of the methods.


## (4) Later... introduce a class that takes care of the monomial lifting etc.

class BasicSemialgebraicSet_veronese(BasicSemialgebraicSet_base):

    """
    A basic semialgebraic set that delegates to another semialgebraic set
    via a Veronese embedding (RLT).
    """

    pass
