r"""
Mutable basic semialgebraic sets
"""

from __future__ import division, print_function, absolute_import

from sage.structure.element import Element
from sage.modules.free_module_element import vector
from sage.rings.all import QQ, ZZ
from sage.rings.real_double import RDF
from sage.misc.abstract_method import abstract_method
from sage.arith.misc import gcd

from copy import copy
import operator

try:
    from ppl import Variable, Constraint, Linear_Expression, Constraint_System, NNC_Polyhedron, Poly_Con_Relation, Poly_Gen_Relation, Generator, MIP_Problem, point as ppl_point
except ImportError:
    # old Sage
    from sage.libs.ppl import Variable, Constraint, Linear_Expression, Constraint_System, NNC_Polyhedron, Poly_Con_Relation, Poly_Gen_Relation, Generator, MIP_Problem, point as ppl_point

from sage.numerical.mip import MixedIntegerLinearProgram, MIPVariable, MIPSolverException

from sage.structure.sage_object import SageObject
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

class BasicSemialgebraicSet_base(SageObject):    # SageObject until we decide if this should be facade parent, or an Element.

    """
    Abstract base class of mutable basic semialgebraic sets.
    """

    def __init__(self, base_ring, ambient_dim):
        """
        Initialize a basic semialgebraic set as the universe in ``ambient_dim``.

        ``base_ring`` is the ring in which the coefficients of polynomials live.
        """
        super(BasicSemialgebraicSet_base, self).__init__()
        self._ambient_dim = ambient_dim
        self._base_ring = base_ring

    def ambient_dim(self):
        return self._ambient_dim

    def ambient_space(self):
        return self.base_ring() ** self.ambient_dim()

    def base_ring(self):
        return self._base_ring

    def poly_ring(self):
        return PolynomialRing(self.base_ring(), "x", self.ambient_dim())

    @abstract_method
    def eq_poly(self):
        r"""
        Return a list, set, or iterator of the polynomials `f` in equations `f(x) = 0`
        in the description of ``self``.

        Together, ``eq_poly``, ``lt_poly``, and ``le_poly`` describe ``self``.
        """

    @abstract_method
    def lt_poly(self):
        r"""
        Return a list, set, or iterator of the polynomials `f` in strict inequalities `f(x) < 0`
        in the description of ``self``.

        Together, ``eq_poly``, ``lt_poly``, and ``le_poly`` describe ``self``.
        """

    @abstract_method
    def le_poly(self):
        r"""
        Return a list, set, or iterator of the polynomials `f` in inequalities `f(x) \leq 0`
        in the description of ``self``.

        Together, ``eq_poly``, ``lt_poly``, and ``le_poly`` describe ``self``.
        """

    def __contains__(self, point):
        """
        Whether the set contains the ``point`` (vector).
        """
        return all(f(point) == 0 for f in self.eq_poly()) and all(f(point) <= 0 for f in self.le_poly()) and all(f(point) < 0 for f in self.lt_poly())

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
        ## right now lhs needs to be a PPL thing. to be changed, and add rhs.
        """
        ``lhs`` should be a vector of length ambient_dim.
        Add the constraint ``lhs`` * x ``op`` 0,
        where ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.
        """
        # default implementation should call add_polynomial_constraint

    @abstract_method
    def is_linear_constraint_valid(self, lhs, op):
        ## right now lhs needs to be a PPL thing. to be changed, and add rhs.
        """
        Whether the constraint ``lhs`` * x ``op`` 0
        is satisfied for all points of ``self``.
        """
        # default implementation should call is_linear_constraint_valid

    @abstract_method
    def add_polynomial_constraint(self, lhs, op):
        """
        ``lhs`` should be a polynomial.
        Add the constraint ``lhs``(x) ``op`` 0,
        where ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.
        """

    @abstract_method
    def is_polynomial_constraint_valid(self, lhs, op):
        """
        Whether the constraint ``lhs`` * x ``op``
        is satisfied for all points of ``self``.
        """

    @abstract_method
    def find_point(self):
        """
        Find a point in ``self``.
        """
        # default implementation could go through self.closure_polyhedron()

    def section(self, polynomial_map, bsa_class='section'):
        r"""
        Define the semialgebraic set that is a section of ``self``.

        ``polynomial_map`` is a vector, list, or tuple of `n` polynomials
        defining a polynomial map `F` from some R^m to R^n where `n` is
        ``self.ambient_dim()``.

        This defines the section `X = \{ x \in R^m : F(x) \in Y \}`,
        where `Y` is ``self``.

        """
        bsa_section = BasicSemialgebraicSet_section(self, polynomial_map)
        if bsa_class == 'section':
            bsa_class = BasicSemialgebraicSet_section
        if bsa_class == BasicSemialgebraicSet_section:
            return bsa_section
        else:
            return bsa_class.from_bsa(bsa_section)

    def plot(self, alpha=0.5, plot_points=300, slice_value=None, **kwds):
        r"""
        Plot the semialgebraic set or a slice (section) of it.

        - If slice_value is given, plot the slice of the cell according to the parameter values in slice_value that are not None. See examples in ``SemialgebraicComplex.plot()``.
        - plot_points controls the quality of the plotting.

        Plot the slice in (x,y)-space with z=4::

            sage: complex.plot(slice_value=[None, None, 4])          # not tested

        Plot the slice in (y,z)-space with x=4::

            sage: complex.plot(slice_value=[4, None, None])          # not tested
        """
        ## Refactor SemialgebraicComplexComponent.plot and plot2dslice through this method.
        raise NotImplementedError()


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
        base_ring = QQ
        super(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron, self).__init__(base_ring, ambient_dim)
        if polyhedron is None:
            self._polyhedron = NNC_Polyhedron(ambient_dim, 'universe')
        else:
            self._polyhedron = polyhedron

    def poly_ring(self):
        """
        Return the polynomial ring.  Variable names match that of the PPL polyhedron:
        x0, x1, ...

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(2)
            sage: P.poly_ring()
            Multivariate Polynomial Ring in x0, x1 over Rational Field

        """
        # Start at 0 to match the names on the PPL side.
        return PolynomialRing(self.base_ring(), ["x{}".format(i) for i in range(self.ambient_dim())])

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

    def _repr_(self):
        return 'BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(polyhedron={})'.format(self._polyhedron)

    def eq_poly(self):
        r"""
        Return a list of the polynomials `f` in equations `f(x) = 0`
        in the description of ``self``.

        Together, ``eq_poly``, ``lt_poly``, and ``le_poly`` describe ``self``.
        """
        mineq = []
        for c in self._polyhedron.minimized_constraints():
            if c.is_equality():
                coeff = c.coefficients()
                # observe: coeffients in a constraint of NNC_Polyhedron could have gcd != 1. # take care of this.
                gcd_c = gcd(gcd(coeff), c.inhomogeneous_term())
                t = sum(QQ(x)/gcd_c*y for x, y in zip(coeff, self.poly_ring().gens())) + QQ(c.inhomogeneous_term())/gcd_c
                mineq.append(t)
        return mineq

    def lt_poly(self):
        r"""
        Return a list of the polynomials `f` in strict inequalities `f(x) < 0`
        in the description of ``self``.

        Together, ``eq_poly``, ``lt_poly``, and ``le_poly`` describe ``self``.
        """
        minlt = []
        for c in self._polyhedron.minimized_constraints():
            if c.is_strict_inequality():
                coeff = c.coefficients()
                gcd_c = gcd(gcd(coeff), c.inhomogeneous_term())
                # constraint is written with '>', while lt_poly records '<' relation
                t = sum(-QQ(x)/gcd_c*y for x, y in zip(coeff, self.poly_ring().gens())) - QQ(c.inhomogeneous_term())/gcd_c
                minlt.append(t)
        return minlt

    def le_poly(self):
        r"""
        Return a list of the polynomials `f` in inequalities `f(x) \leq 0`
        in the description of ``self``.

        Together, ``eq_poly``, ``lt_poly``, and ``le_poly`` describe ``self``.
        """
        minle = []
        for c in self._polyhedron.minimized_constraints():
            if c.is_nonstrict_inequality():
                coeff = c.coefficients()
                gcd_c = gcd(gcd(coeff), c.inhomogeneous_term())
                # constraint is written with '>=', while lt_poly records '<=' relation
                t = sum(-QQ(x)/gcd_c*y for x, y in zip(coeff, self.poly_ring().gens())) - QQ(c.inhomogeneous_term())/gcd_c
                minle.append(t)
        return minle

    # override the default implementation
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

    # override the abstract methods
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

    def __init__(self, base_ring, ambient_dim, solver=None):
        """
        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram(QQ, 3)

        """
        super(BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram, self).__init__(base_ring, ambient_dim)
        self._mip = MixedIntegerLinearProgram(solver=solver)

    def __copy__(self):
        raise NotImplementedError()


## (3) Then introduce the following class to simplify the code in parametric.sage
class BasicSemialgebraicSet_eq_lt_le_sets(BasicSemialgebraicSet_base):

    """
    A basic semialgebraic set, represented in a straightforward way
    as 3 finite sets of polynomial constraints `p(x) OP 0`.
    """

    def __init__(self, base_ring=None, ambient_dim=None, poly_ring=None, eq=[], lt=[], le=[]):
        """
        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_eq_lt_le_sets(QQ, 3)

        """
        if poly_ring is None:
            polys = list(eq) + list(lt) + list(le)
            if polys:
                poly_ring = polys[0].parent()
        if (ambient_dim is None) and (poly_ring is not None):
            ambient_dim = poly_ring.ngens()
        if base_ring is None and poly_ring is not None:
            base_ring = poly_ring.base_ring()
        if ambient_dim is None or base_ring is None:
            raise ValueError("if eq, lt, and le are all empty, must provide either poly_ring or both of base_ring and ambient_dim")
        if poly_ring is None:
            poly_ring = PolynomialRing(base_ring, "x", ambient_dim)
        super(BasicSemialgebraicSet_eq_lt_le_sets, self).__init__(base_ring, ambient_dim)
        self._poly_ring = poly_ring
        self._eq = set(eq)
        self._lt = set(lt)
        self._le = set(le)

    def poly_ring(self):
        return self._poly_ring

    def __copy__(self):
        """
        Make a copy of ``self``.

        TESTS:

        Test that it is actually making a copy of the (mutable!) sets::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P = BasicSemialgebraicSet_eq_lt_le_sets(QQ, 2)
            sage: P._eq is copy(P)._eq
            False
            sage: P._lt is copy(P)._lt
            False
            sage: P._le is copy(P)._le
            False
        """
        return self.__class__(self.base_ring(), self.ambient_dim(), self.poly_ring(),
                              self._eq, self._lt, self._le)

    # override the abstract methods

    def __repr__(self):
        return 'BasicSemialgebraicSet_eq_lt_le_sets(eq = {}, lt = {}, le = {})'.format(list(self._eq), list(self._lt), list(self._le))

    def eq_poly(self):
        r"""
        Return a list, set, or iterator of the polynomials `f` in equations `f(x) = 0`
        in the description of ``self``.

        Together, ``eq_poly``, ``lt_poly``, and ``le_poly`` describe ``self``.
        """
        return self._eq

    def lt_poly(self):
        r"""
        Return a list, set, or iterator of the polynomials `f` in strict inequalities `f(x) < 0`
        in the description of ``self``.

        Together, ``eq_poly``, ``lt_poly``, and ``le_poly`` describe ``self``.
        """
        return self._lt

    def le_poly(self):
        r"""
        Return a list, set, or iterator of the polynomials `f` in inequalities `f(x) \leq 0`
        in the description of ``self``.

        Together, ``eq_poly``, ``lt_poly``, and ``le_poly`` describe ``self``.
        """
        return self._le

    def add_polynomial_constraint(self, lhs, op):
        """
        ``lhs`` should be a polynomial.
        Add the constraint ``lhs``(x) ``op`` 0,
        where ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.
        """
        if op == operator.lt:
            self._lt.add(lhs)
        elif op == operator.eq:
            self._eq.add(lhs)
        elif op == operator.le:
            self._le.add(lhs)
        else:
            raise ValueError("{} is not a supported operator".format(op))


    ### TODO: Add implementations of more of the methods.


class BasicSemialgebraicSet_section(BasicSemialgebraicSet_base):

    """
    Section of another ``BasicSemialgebraicSet``.
    """

    def __init__(self, upstairs_bsa, polynomial_map, ambient_dim=None):
        """
        EXAMPLES:

        McCormick::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: PY.<x,y,xy> = QQ[]
            sage: bsa = BasicSemialgebraicSet_eq_lt_le_sets(le=[-x, -y, -xy, x-2, y-3, xy-6])
            sage: PX.<u,v> = QQ[]
            sage: F = [u, v, u*v]
            sage: section = bsa.section(F); section
            BasicSemialgebraicSet_section(BasicSemialgebraicSet_eq_lt_le_sets(eq = [], lt = [], le = [-y, -x, xy - 6, y - 3, x - 2, -xy]), polynomial_map=[u, v, u*v])
            sage: section.ambient_dim()
            2
            sage: sorted(section.eq_poly()), sorted(section.le_poly())
            ([], [-v, v - 3, -u, u - 2, -u*v, u*v - 6])

        An interval::

            sage: PZ.<z> = QQ[]
            sage: interval = bsa.section([1, 1, z])
            sage: sorted(f for f in interval.le_poly() if f not in QQ)
            [-z, z - 6]

        A 0-dimensional section of 0-dimensional set::

            sage: Q0 = PolynomialRing(QQ, [])
            sage: point = BasicSemialgebraicSet_eq_lt_le_sets(poly_ring=Q0)
            sage: section = point.section([])
            sage: section.poly_ring()
            Multivariate Polynomial Ring in no variables over Rational Field

        """
        if len(polynomial_map) != upstairs_bsa.ambient_dim():
            raise ValueError("polynomial_map must have the same length as the dimension of the underlying bsa")
        if polynomial_map and ambient_dim is None:
            ambient_dim = polynomial_map[0].parent().ngens()
        if not all(poly.parent().ngens() == ambient_dim for poly in polynomial_map):
            raise ValueError("elements in polynomial_map must come from a polynomial ring with the same number of variables as the ambient dimension")
        base_ring = upstairs_bsa.base_ring()
        super(BasicSemialgebraicSet_section, self).__init__(base_ring, ambient_dim)
        self._upstairs_bsa = upstairs_bsa
        self._polynomial_map = polynomial_map

    def poly_ring(self):
        if not self._polynomial_map:
            return PolynomialRing(self.base_ring(), [])
        else:
            return polynomial_map[0].parent()

    def eq_poly(self):
        for p in self._upstairs_bsa.eq_poly():
            yield p(self._polynomial_map)

    def le_poly(self):
        for p in self._upstairs_bsa.le_poly():
            yield p(self._polynomial_map)

    def lt_poly(self):
        for p in self._upstairs_bsa.lt_poly():
            yield p(self._polynomial_map)

    def _repr_(self):
        return 'BasicSemialgebraicSet_section({}, polynomial_map={})'.format(self._upstairs_bsa, self._polynomial_map)

## (4) Later... introduce a class that takes care of the monomial lifting etc.

class BasicSemialgebraicSet_veronese(BasicSemialgebraicSet_section):

    """
    A basic semialgebraic set that delegates to another semialgebraic set
    via a Veronese embedding (RLT).
    """

    def __init__(self, upstairs_bsa, monomial_list, v_dict, ambient_dim=None):
        """
        EXAMPLES:

        Trivial initialization of the universe in dimension 3::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: upstairs_bsa_ppl = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(ambient_dim=0)
            sage: veronese = BasicSemialgebraicSet_veronese(upstairs_bsa_ppl, [], dict(), ambient_dim=3)
            sage: veronese.ambient_space()
            Vector space of dimension 3 over Rational Field

        Adding initial space dimensions::

            sage: names = ['u', 'v']
            sage: n = len(names)
            sage: P = PolynomialRing(QQ, names)
            sage: monomial_list = list(P.gens())
            sage: v_dict = {P.gens()[i]:i for i in range(n)}
            sage: polyhedron = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(ambient_dim=n)
            sage: veronese = BasicSemialgebraicSet_veronese(polyhedron, monomial_list, v_dict)
            sage: veronese.ambient_space()
            Vector space of dimension 2 over Rational Field

        """
        super(BasicSemialgebraicSet_veronese, self).__init__(upstairs_bsa, monomial_list,
                                                             ambient_dim=ambient_dim)
        self._v_dict = v_dict

    def monomial_list(self):
        """
        A list that maps the index of each generator in ``self.upstairs()``
        to a monomial.
        """
        return self._polynomial_map

    def v_dict(self):
        """
        A dictionary that maps each monomial to the index of its corresponding generator
        in ``self.upstairs()``.
        """
        return self._v_dict()

    def add_polynomial_constraint(self, lhs, op):
        """
        ``lhs`` should be a polynomial.
        Add the constraint ``lhs``(x) ``op`` 0,
        where ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.
        """
        # Like polynomial_to_linexpr.
        # Call self.upstairs().add_space_dimensions_and_embed when needed.
        raise NotImplementedError()
