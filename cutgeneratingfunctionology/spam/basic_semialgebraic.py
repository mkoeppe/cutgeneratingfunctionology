r"""
Mutable basic semialgebraic sets
"""

from __future__ import division, print_function, absolute_import

from sage.structure.element import Element
from sage.modules.free_module_element import vector
from sage.rings.all import QQ, ZZ
from sage.rings.real_double import RDF
from sage.rings.infinity import Infinity
from sage.misc.abstract_method import abstract_method
from sage.arith.misc import gcd

from copy import copy
from itertools import chain
import operator

try:
    from ppl import Variable, Constraint, Linear_Expression, Constraint_System, NNC_Polyhedron, Poly_Con_Relation, Poly_Gen_Relation, Generator, MIP_Problem, point as ppl_point
except ImportError:
    # old Sage
    from sage.libs.ppl import Variable, Constraint, Linear_Expression, Constraint_System, NNC_Polyhedron, Poly_Con_Relation, Poly_Gen_Relation, Generator, MIP_Problem, point as ppl_point

from sage.numerical.mip import MixedIntegerLinearProgram, MIPVariable, MIPSolverException

from sage.structure.sage_object import SageObject
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

def _bsa_class(bsa_class):
    """
    Translate a class nickname to a class.
    """
    if isinstance(bsa_class, type):
        return bsa_class
    elif bsa_class == 'formal_closure':
        return BasicSemialgebraicSet_formal_closure
    elif bsa_class == 'section':
        return BasicSemialgebraicSet_section
    elif bsa_class == 'ppl':
        return BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron
    elif bsa_class == 'mip':
        return BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram
    else:
        raise ValueError("unknown bsa class: {}".format(bsa_class))

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

    # Default implementation. Subclasses are encouraged to provide
    # faster implementations.
    @classmethod
    def from_bsa(cls, bsa, **init_kwds):
        """
        Initialize a basic semialgebraic set of class ``cls`` to be the same
        as ``bsa``.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: upstairs_bsa_ppl = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(ambient_dim=0)
            sage: veronese = BasicSemialgebraicSet_veronese(upstairs_bsa_ppl, [], dict(), ambient_dim=0)
            sage: Q.<x0,x1,x2> = QQ[]
            sage: lhs = 27/113 * x0^2 + x1*x2 + 1/2
            sage: veronese.add_polynomial_constraint(lhs, operator.lt)
            sage: list(veronese.lt_poly())
            [54*x0^2 + 226*x1*x2 + 113]
            sage: bsa_eq_lt_le_sets = BasicSemialgebraicSet_eq_lt_le_sets.from_bsa(veronese)
            sage: list(bsa_eq_lt_le_sets.lt_poly())
            [54*x0^2 + 226*x1*x2 + 113]
            sage: upstairs_bsa_eq_lt_le_sets = BasicSemialgebraicSet_eq_lt_le_sets.from_bsa(upstairs_bsa_ppl)
            sage: list(upstairs_bsa_eq_lt_le_sets.lt_poly())
            [54*x0 + 226*x1 + 113]

        Test that ``BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron`` supports this method::

            sage: upstairs_bsa_ppl_again = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron.from_bsa(upstairs_bsa_eq_lt_le_sets)
            sage: list(upstairs_bsa_ppl_again.lt_poly())
            [54*x0 + 226*x1 + 113]

        Test that ``BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram``
        supports this method and accepts keyword arguments::

            sage: upstairs_bsa_mip = BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram.from_bsa(upstairs_bsa_eq_lt_le_sets, solver='ppl')    # not tested - unfinished implementation

        """
        base_ring = init_kwds.pop('base_ring', bsa.base_ring())
        ambient_dim = bsa.ambient_dim()
        self = cls(base_ring=base_ring, ambient_dim=ambient_dim, **init_kwds)
        for p in bsa.eq_poly():
            self.add_polynomial_constraint(p, operator.eq)
        for p in bsa.lt_poly():
            self.add_polynomial_constraint(p, operator.lt)
        for p in bsa.le_poly():
            self.add_polynomial_constraint(p, operator.le)
        return self

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

    def formal_closure(self, bsa_class='formal_closure'):
        """
        Return the basic semialgebraic set obtained by replacing all strict
        inequalities by <= inequalities.  This is a superset of the topological closure.

        By default, the formal closure is represented by an instance of class
        ``BasicSemialgebraicSet_formal_closure``; use the argument ``bsa_class``
        to choose another class.  See ``_bsa_class`` for the allowed class nicknames.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: bsa = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(ambient_dim=2)
            sage: bsa.add_linear_constraint([1, 1], -3, operator.lt)
            sage: list(bsa.lt_poly()), list(bsa.le_poly())
            ([x0 + x1 - 3], [])
            sage: closure = bsa.formal_closure(); closure
            BasicSemialgebraicSet_formal_closure(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(...))
            sage: list(closure.eq_poly()), list(closure.lt_poly()), list(closure.le_poly())
            ([], [], [x0 + x1 - 3])
            sage: closure = bsa.formal_closure(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron); closure
            BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(polyhedron=A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 1 ray, 1 line)
            sage: list(closure.eq_poly()), list(closure.lt_poly()), list(closure.le_poly())
            ([], [], [x0 + x1 - 3])

        In general, this is only a superset of the topological closure::

            sage: R.<x> = QQ['x']
            sage: bsa = BasicSemialgebraicSet_eq_lt_le_sets(lt={x, -x})
            sage: formal_closure = bsa.formal_closure(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron)
            sage: list(formal_closure.eq_poly())
            [x0]
            sage: bsa_ppl = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron.from_bsa(bsa)
            sage: closure = bsa_ppl.closure()
            sage: list(closure.eq_poly())
            [-1]

        """
        bsa_formal_closure = BasicSemialgebraicSet_formal_closure(self)
        bsa_class = _bsa_class(bsa_class)
        if bsa_class == BasicSemialgebraicSet_formal_closure:
            return bsa_formal_closure
        else:
            return bsa_class.from_bsa(bsa_formal_closure)

    @abstract_method
    def closure(self, bsa_class):
        """
        Return the basic semialgebraic set that is the topological closure
        of ``self``.
        """

    def add_space_dimensions_and_embed(self, space_dim_to_add):
        """
        Mutate ``self`` by injecting it into a higher dimensional space.
        """
        self._ambient_dim += space_dim_to_add

    def add_linear_constraint(self, lhs, cst, op):
        """
        Add the constraint ``lhs`` * x + cst ``op`` 0,
        where ``lhs`` is a vector of length ``ambient_dim`` and
        ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.

        This implementation just calls ``self.add_polynomial_constraint``.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: bsa = BasicSemialgebraicSet_eq_lt_le_sets(base_ring=QQ, ambient_dim=2)
            sage: lhs = vector(QQ, [1, 1])
            sage: bsa.add_linear_constraint(lhs, -2, operator.lt)
            sage: list(bsa.lt_poly())
            [x0 + x1 - 2]
        """
        poly = sum(coeff * gen for coeff, gen in zip(lhs, self.poly_ring().gens())) + cst
        self.add_polynomial_constraint(poly, op)

    def is_linear_constraint_valid(self, lhs, cst, op):
        """
        Whether the constraint ``lhs`` * x + cst ``op`` 0
        is satisfied for all points of ``self``.

        This implementation uses ``linear_function_upper_bound``
        and ``linear_function_lower_bound`` and raises an error
        if the information provided does not suffice to decide
        the validity of the constraint.
        """
        op_known = False
        if op in (operator.lt, operator.le, operator.eq):
            op_known = True
            if not op(self.linear_function_upper_bound(lhs) + cst, 0):
                raise NotImplementedError
        if op in (operator.gt, operator.ge, operator.eq):  # eq appears in both conditions
            op_known = True
            if not op(self.linear_function_lower_bound(lhs) + cst, 0):
                raise NotImplementedError
        if not op_known:
            raise ValueError("{} is not a supported operator".format(op))
        return True

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

    def linear_function_upper_bound(self, form):
        """
        Find an upper bound for ``form`` (a vector) on ``self``.

        The default implementation just returns +oo.
        Subclasses should provide more useful bounds.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: bsa = BasicSemialgebraicSet_base(QQ, 2)
            sage: bsa.linear_function_upper_bound(vector(QQ, [1, 1]))
            +Infinity

        """
        return +Infinity

    def linear_function_lower_bound(self, form):
        """
        Find a lower bound for ``form`` (a vector) on ``self``.

        This implementation delegates to ``linear_function_upper_bound``.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: bsa = BasicSemialgebraicSet_base(QQ, 2)
            sage: bsa.linear_function_lower_bound(vector(QQ, [1, 1]))
            -Infinity

        """
        form = vector(form)
        return -self.linear_function_upper_bound(-form)

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

        By default, the section is represented by an instance of class
        ``BasicSemialgebraicSet_section``; use the argument ``bsa_class``
        to choose another class.  See ``_bsa_class`` for the allowed
        class nicknames.

        """
        bsa_section = BasicSemialgebraicSet_section(self, polynomial_map)
        bsa_class = _bsa_class(bsa_class)
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

class BasicSemialgebraicSet_polyhedral(BasicSemialgebraicSet_base):

    """
    An abstract class of polyhedral basic semialgebraic sets.

    """

    @abstract_method
    def add_linear_constraint(self, lhs, cst, op):
        """
        Add the constraint ``lhs`` * x + cst ``op`` 0,
        where ``lhs`` is a vector of length ``ambient_dim`` and
        ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.

        In subclasses of ``BasicSemialgebraicSet_polyhedral``,
        this method should be defined.
        """

    def add_polynomial_constraint(self, lhs, op):
        """
        Add the constraint ``lhs``(x) ``op`` 0,
        where ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.

        ``lhs`` should be a linear polynomial that can be converted to
        an element of ``self.poly_ring()``.

        This implementation checks that ``lhs`` is a linear polynomial
        and then delegates to ``add_linear_constraint``.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(1)
            sage: f = 2 * P.poly_ring().gen() + 1
            sage: P.add_polynomial_constraint(f, operator.le)
            sage: list(P.le_poly())
            [2*x0 + 1]

        With conversion of polynomials from different rings::

            sage: R.<xyzzy> = QQ[]
            sage: P.add_polynomial_constraint(xyzzy + 3, operator.lt)
            sage: list(P.lt_poly())
            [x0 + 3]

        """
        lhs = self.poly_ring()(lhs)   # convert if necessary
        if lhs.degree() > 1:
            raise ValueError("{} is not a valid linear polynomial.".format(lhs))
        cst = lhs.constant_coefficient()
        try:
            lhs_vector = vector(lhs.coefficient(x) for x in self.poly_ring().gens())
        except AttributeError:
            # univariate polynomials (Polynomial_rational_flint) unfortunately
            # have a different API that does not define "coefficient".
            lhs_vector = [lhs[1]]  # coefficient of the linear term
        self.add_linear_constraint(lhs_vector, cst, op)



## (1) In the first step, we implement the following class.  Everything is linear.
## Rewrite all direct uses of PPL in ParametricRealFieldElement, ParametricRealField
## using method calls to this class.

poly_is_included = Poly_Con_Relation.is_included()
#strictly_intersects = Poly_Con_Relation.strictly_intersects()
point_is_included = Poly_Gen_Relation.subsumes()
#con_saturates = Poly_Con_Relation.saturates()

from sage.arith.functions import lcm

class BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(BasicSemialgebraicSet_polyhedral):

    """
    A (possibly half-open) polyhedral basic semialgebraic set,
    represented by a PPL ``NNC_Polyhedron``

    """

    def __init__(self, ambient_dim=None, polyhedron=None, base_ring=None):
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
        if base_ring is None:
            base_ring = QQ
        if base_ring is not QQ:
            raise ValueError("only base_ring=QQ is supported")
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

    def closure(self, bsa_class='formal_closure'):
        """
        Return the basic semialgebraic set that is the topological closure
        of ``self``.
        """
        # Because our description consists of minimized constraints, the closure is
        # just the formal closure.
        return self.formal_closure(bsa_class=bsa_class)

    def eq_poly(self):
        r"""
        Return a list of the polynomials `f` in equations `f(x) = 0`
        in the description of ``self``.

        Together, ``eq_poly``, ``lt_poly``, and ``le_poly`` describe ``self``.
        """
        for c in self._polyhedron.minimized_constraints():
            if c.is_equality():
                coeff = c.coefficients()
                # observe: coeffients in a constraint of NNC_Polyhedron could have gcd != 1. # take care of this.
                gcd_c = gcd(gcd(coeff), c.inhomogeneous_term())
                t = sum(QQ(x)/gcd_c*y for x, y in zip(coeff, self.poly_ring().gens())) + QQ(c.inhomogeneous_term())/gcd_c
                yield t

    def lt_poly(self):
        r"""
        Return a list of the polynomials `f` in strict inequalities `f(x) < 0`
        in the description of ``self``.

        Together, ``eq_poly``, ``lt_poly``, and ``le_poly`` describe ``self``.
        """
        for c in self._polyhedron.minimized_constraints():
            if c.is_strict_inequality():
                coeff = c.coefficients()
                gcd_c = gcd(gcd(coeff), c.inhomogeneous_term())
                # constraint is written with '>', while lt_poly records '<' relation
                t = sum(-QQ(x)/gcd_c*y for x, y in zip(coeff, self.poly_ring().gens())) - QQ(c.inhomogeneous_term())/gcd_c
                yield t

    def le_poly(self):
        r"""
        Return a list of the polynomials `f` in inequalities `f(x) \leq 0`
        in the description of ``self``.

        Together, ``eq_poly``, ``lt_poly``, and ``le_poly`` describe ``self``.
        """
        for c in self._polyhedron.minimized_constraints():
            if c.is_nonstrict_inequality():
                coeff = c.coefficients()
                gcd_c = gcd(gcd(coeff), c.inhomogeneous_term())
                # constraint is written with '>=', while lt_poly records '<=' relation
                t = sum(-QQ(x)/gcd_c*y for x, y in zip(coeff, self.poly_ring().gens())) - QQ(c.inhomogeneous_term())/gcd_c
                yield t

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

    @staticmethod
    def _ppl_constraint(lhs, cst, op):
        """
        Make a PPL ``Constraint`` ``lhs`` * x + cst ``op`` 0.
        """
        lcd = lcm(lcm(x.denominator() for x in lhs), cst.denominator())
        linexpr = Linear_Expression(lhs * lcd, cst * lcd)
        if op == operator.lt:
            return (linexpr < 0)
        elif op == operator.gt:
            return (linexpr > 0)
        elif op == operator.eq:
            return (linexpr == 0)
        elif op == operator.le:
            return (linexpr <= 0)
        elif op == operator.ge:
            return (linexpr >= 0)
        else:
            raise ValueError("{} is not a supported operator".format(op))

    def is_linear_constraint_valid(self, lhs, cst, op):
        """
        Whether the constraint ``lhs`` * x + cst ``op`` 0
        is satisfied for all points of ``self``.
        """
        constraint = self._ppl_constraint(lhs, cst, op)
        return self._polyhedron.relation_with(constraint).implies(poly_is_included)

    def add_linear_constraint(self, lhs, cst, op):
        """
        ``lhs`` should be a vector of length ambient_dim.
        Add the constraint ``lhs`` * x + cst ``op`` 0,
        where ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.
        """
        constraint = self._ppl_constraint(lhs, cst, op)
        self._polyhedron.add_constraint(constraint)


## class BasicSemialgebraicSet_polyhedral_ppl_MIP_Problem(BasicSemialgebraicSet_base):

##     """
##     A closed polyhedral basic semialgebraic set,
##     represented by a PPL ``MIP_Problem``.
##     """

### do we need the above, or should we just use the below, using solver='ppl'?


## (2) Then
class BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram(BasicSemialgebraicSet_polyhedral):

    """
    A closed polyhedral basic semialgebraic set,
    represented by a Sage ``MixedIntegerLinearProgram``
    with only continuous variables.
    """

    def __init__(self, base_ring, ambient_dim, solver=None):
        """
        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram(QQ, 3)

        """
        super(BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram, self).__init__(base_ring, ambient_dim)
        self._mip = MixedIntegerLinearProgram(solver=solver, maximization=True)

    def __copy__(self):
        raise NotImplementedError()

    def mip(self):
        return self._mip

    def mip_gens(self):
        """
        Return the components of the MIP variable corresponding to the
        space dimensions.
        """
        mip_var = self.mip().default_variable()
        for i in range(self.ambient_dim()):
            yield mip_var[i]

    def closure(self, bsa_class='mip'):
        """
        Return the basic semialgebraic set that is the topological closure
        of ``self``, which is ``self`` itself.
        """
        return self

    def is_linear_constraint_valid(self, lhs, cst, op):
        """
        Whether the constraint ``lhs`` * x + cst ``op`` 0
        is satisfied for all points of ``self``.

        In the current implementation, this raises ``NotImplementedError``
        in unbounded/infeasible situations.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram(QQ, 1, solver='ppl')
            sage: S.add_linear_constraint([1], -1, operator.le)
            sage: S.add_linear_constraint([1], 0, operator.ge)
            sage: S.is_linear_constraint_valid([1], -2, operator.lt)
            True
            sage: S.is_linear_constraint_valid([1], -1, operator.lt)
            False
            sage: S.is_linear_constraint_valid([1], -1, operator.le)
            True
            sage: S.is_linear_constraint_valid([1], 1/2, operator.le)
            False

        """
        try:
            # Use the default implementation, which delegates to
            # ``linear_function_upper_bound`` and
            # ``linear_function_lower_bound``.
            return super(BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram, self).is_linear_constraint_valid(lhs, cst, op)
        except NotImplementedError:
            # As long as ``linear_function_upper_bound`` cannot distinguish
            # infeasible from unbounded, we have a problem here.
            if self.linear_function_upper_bound(lhs) == Infinity or self.linear_function_lower_bound(lhs) == -Infinity:
                raise NotImplementedError
            # The default implementation checked
            # ``linear_function_upper_bound`` and
            # ``linear_function_lower_bound`` already.  These are the sup and
            # inf in our case, and because we are a closed polyhedral set, this
            # suffices to decide validity.
            return False

    def _mip_linear_function(self, form, cst=0):
        """
        Obtain a linear function object representing
        ``form`` * x + ``cst``.
        """
        return self.mip().sum(coeff * gen for coeff, gen in zip(form, self.mip_gens())) + cst

    def linear_function_upper_bound(self, form):
        """
        Find an upper bound for ``form`` (a vector) on ``self``.

        In this implementation, this is done by solving the LP,
        so this upper bound is the supremum.

        However, if ``self`` is empty, this may return +oo
        because we cannot distinguish infeasible from unbounded.

        """
        mip = self.mip()
        objective = self._mip_linear_function(form)
        mip.set_objective(objective)
        try:
            return mip.solve()
        except MIPSolverException as e:
            # FIXME: Here we should really distinguish infeasible and unbounded.
            # Unfortunately there is no solver-independent protocol for this.
            return +Infinity

    def add_linear_constraint(self, lhs, cst, op):
        """
        Add the constraint ``lhs`` * x + cst ``op`` 0,
        where ``lhs`` is a vector of length ``ambient_dim`` and
        ``op`` is one of ``operator.eq``, ``operator.le``, ``operator.ge``.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram(QQ, 1, solver='ppl')
            sage: S.add_linear_constraint([1], -1, operator.lt)
            Traceback (most recent call last):
            ...
            ValueError: strict < is not allowed, use <= instead
            sage: S.add_linear_constraint([1], -1, operator.le)
            sage: S.add_linear_constraint([1], 0, operator.ge)
            sage: S.linear_function_upper_bound([1])
            1
            sage: S.linear_function_lower_bound([1])
            0

        """
        constraint_lhs = self._mip_linear_function(lhs, cst)
        constraint = op(constraint_lhs, 0)
        self.mip().add_constraint(constraint)

    def eq_poly(self):
        """
        Generate the polynomials `f` in equations `f(x) = 0`
        in the description of ``self``.

        Together, ``eq_poly`` and ``le_poly`` describe ``self``.
        (``lt_poly`` gives the empty list because ``self`` is closed.)

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram(QQ, 1, solver='ppl')
            sage: S.add_linear_constraint([1], -1, operator.eq)
            sage: sorted(S.eq_poly())
            [x - 1]
        """
        poly_ring = self.poly_ring()
        for lb, (indices, coeffs), ub in self.mip().constraints():
            if lb == ub:
                form = sum(coeff * poly_ring.gen(i) for i, coeff in zip(indices, coeffs))
                yield form - ub

    def le_poly(self):
        """
        Generate the polynomials `f` in inequalities `f(x) \leq 0`
        in the description of ``self``.

        Together, ``eq_poly`` and ``le_poly`` describe ``self``.
        (``lt_poly`` gives the empty list because ``self`` is closed.)

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram(QQ, 1, solver='ppl')
            sage: S.add_linear_constraint([1], -1, operator.le)
            sage: S.add_linear_constraint([1], 0, operator.ge)
            sage: sorted(S.le_poly())
            [-x, x - 1]
        """
        poly_ring = self.poly_ring()
        for lb, (indices, coeffs), ub in self.mip().constraints():
            if lb != ub:
                form = sum(coeff * poly_ring.gen(i) for i, coeff in zip(indices, coeffs))
                if lb not in (None, -Infinity):
                    yield lb - form
                if ub not in (None, +Infinity):
                    yield form - ub

    def lt_poly(self):
        """
        Return the empty list because ``self`` is closed.

        Together, ``eq_poly`` and ``le_poly`` describe ``self``.
        """
        return []

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
        elif op == operator.gt:
            self._lt.add(-lhs)
        elif op == operator.eq:
            self._eq.add(lhs)
        elif op == operator.le:
            self._le.add(lhs)
        elif op == operator.ge:
            self._le.add(-lhs)
        else:
            raise ValueError("{} is not a supported operator".format(op))


    ### TODO: Add implementations of more of the methods.


class BasicSemialgebraicSet_section(BasicSemialgebraicSet_base):

    """
    Section of another ``BasicSemialgebraicSet``.

    See ``BasicSemialgebraicSet_base.section``.
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

    def __copy__(self):
        """
        Make a copy of ``self``.

        """
        return self.__class__(copy(self.upstairs()), copy(self.polynomial_map()),
                              self.ambient_dim())

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

    def upstairs(self):
        return self._upstairs_bsa

    def polynomial_map(self):
        return self._polynomial_map

## (4) Later... introduce a class that takes care of the monomial lifting etc.

class BasicSemialgebraicSet_veronese(BasicSemialgebraicSet_section):

    """
    A basic semialgebraic set that delegates to another semialgebraic set
    via a Veronese embedding (reformulation-linearization, RLT).

    Expand the polynomial in the standard monomial basis and replace each monomial by a new variable.
    Record monomials in monomial_list and their corresponding variables in v_dict. The resulting
    linear expression in the extended space will be provided as inequality or equation in
    upstairs_bsa, which could, for example be represented by a PPL not-necessarily-closed polyhedron.
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

    def __copy__(self):
        """
        Make a copy of ``self``.

        """
        return self.__class__(copy(self.upstairs()), copy(self.polynomial_map()),
                              copy(self.v_dict()), self.ambient_dim())

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
        return self._v_dict

    def add_polynomial_constraint(self, lhs, op):
        """
        ``lhs`` should be a polynomial.
        Add the constraint ``lhs``(x) ``op`` 0,
        where ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: upstairs_bsa_ppl = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(ambient_dim=0)
            sage: veronese = BasicSemialgebraicSet_veronese(upstairs_bsa_ppl, [], dict(), ambient_dim=0)
            sage: Q.<x0,x1,x2> = QQ[]
            sage: lhs = 27/113 * x0^2 + x1*x2 + 1/2
            sage: veronese.add_polynomial_constraint(lhs, operator.lt)
            sage: list(veronese.lt_poly())
            [54*x0^2 + 226*x1*x2 + 113]
            sage: veronese.monomial_list()
            [x0^2, x1*x2]
            sage: veronese.v_dict()
            {x1*x2: 1, x0^2: 0}
            sage: lhs2 = x0 + 1/3*x1*x2
            sage: veronese.add_polynomial_constraint(lhs2, operator.lt)
            sage: list(veronese.lt_poly())
            [54*x0^2 + 226*x1*x2 + 113, x1*x2 + 3*x0]
            sage: veronese.monomial_list()
            [x0^2, x1*x2, x0]
            sage: veronese.v_dict()
            {x0: 2, x1*x2: 1, x0^2: 0}
        """
        space_dim_to_add = 0
        upstairs_lhs_coeff = [0] * self.upstairs().ambient_dim()
        upstairs_lhs_cst = 0
        for m in lhs.monomials():
            coeffm = lhs.monomial_coefficient(m)
            if m == 1:
                upstairs_lhs_cst = coeffm
            else:
                nv = self.v_dict().get(m, None)
                if nv is None:
                    nv = len(self.monomial_list())
                    self.v_dict()[m] = nv
                    self.monomial_list().append(m)
                    space_dim_to_add += 1
                    upstairs_lhs_coeff.append(coeffm)
                else:
                    upstairs_lhs_coeff[nv] = coeffm
        if space_dim_to_add:
            self.upstairs().add_space_dimensions_and_embed(space_dim_to_add)
        upstairs_lhs = sum(QQ(x)*y for x, y in zip(upstairs_lhs_coeff, self.upstairs().poly_ring().gens())) + QQ(upstairs_lhs_cst)
        self.upstairs().add_polynomial_constraint(upstairs_lhs, op)

class BasicSemialgebraicSet_formal_closure(BasicSemialgebraicSet_base):

    r"""
    Represent the formal closure (see method ``formal_closure``) of
    another basic semialgebraic set ``upstairs_bsa``.
    """

    def __init__(self, upstairs_bsa):
        base_ring = upstairs_bsa.base_ring()
        ambient_dim = upstairs_bsa.ambient_dim()
        super(BasicSemialgebraicSet_formal_closure, self).__init__(base_ring, ambient_dim)
        self._upstairs_bsa = upstairs_bsa

    def _repr_(self):
        return 'BasicSemialgebraicSet_formal_closure({})'.format(self.upstairs())

    def upstairs(self):
        return self._upstairs_bsa

    def eq_poly(self):
        return self.upstairs().eq_poly()

    def le_poly(self):
        return chain(self.upstairs().le_poly(), self.upstairs().lt_poly())

    def lt_poly(self):
        return []
