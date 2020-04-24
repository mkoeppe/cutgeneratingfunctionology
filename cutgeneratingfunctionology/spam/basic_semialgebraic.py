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
from sage.plot.graphics import Graphics
from sage.calculus.var import var
from sage.plot.contour_plot import region_plot, implicit_plot
from sage.matrix.constructor import matrix
from sage.structure.element import get_coercion_model
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing

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
    r"""
    Translate a class nickname to a class.
    """
    if isinstance(bsa_class, type):
        return bsa_class
    elif bsa_class == 'formal_closure':
        from .basic_semialgebraic_formal_closure import BasicSemialgebraicSet_formal_closure
        return BasicSemialgebraicSet_formal_closure
    elif bsa_class == 'formal_relint':
        from .basic_semialgebraic_formal_relint import BasicSemialgebraicSet_formal_relint
        return BasicSemialgebraicSet_formal_relint
    elif bsa_class == 'mathematica':
        from .semialgebraic_mathematica import BasicSemialgebraicSet_mathematica
        return BasicSemialgebraicSet_mathematica
    elif bsa_class == 'section':
        return BasicSemialgebraicSet_section
    elif bsa_class == 'ppl':
        return BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron
    elif bsa_class == 'mip':
        return BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram
    elif bsa_class == 'linear_system':
        from .basic_semialgebraic_linear_system import BasicSemialgebraicSet_polyhedral_linear_system
        return BasicSemialgebraicSet_polyhedral_linear_system
    elif bsa_class == 'intersection':
        from .basic_semialgebraic_intersection import BasicSemialgebraicSet_intersection
        return BasicSemialgebraicSet_intersection
    elif bsa_class == 'predicate':
        from .semialgebraic_predicate import BasicSemialgebraicSet_predicate
        return BasicSemialgebraicSet_predicate
    elif bsa_class == 'veronese':
        return BasicSemialgebraicSet_veronese
    elif bsa_class == 'groebner_basis':
        return BasicSemialgebraicSet_groebner_basis
    else:
        raise ValueError("unknown bsa class: {}".format(bsa_class))

class BasicSemialgebraicSet_base(SageObject):    # SageObject until we decide if this should be facade parent, or an Element.

    r"""
    Abstract base class of mutable basic semialgebraic sets.
    """

    @classmethod
    def _poly_ring_from_options(cls, poly_ring=None, base_ring=None, ambient_dim=None, names=None, **options):
        """
        Interpret the `poly_ring`, `base_ring`, `ambient_dim`, and `names` options.

        Return a tuple of these four, normalized.
        """
        if poly_ring is None:
            if names is not None:
                poly_ring = PolynomialRing(base_ring, names)
            elif ambient_dim is not None and base_ring is not None:
                poly_ring = PolynomialRing(base_ring, [cls._default_name(i) for i in range(ambient_dim)])
        if ambient_dim is None:
            if poly_ring is not None:
                ambient_dim = poly_ring.ngens()
            elif names is not None:
                ambient_dim = len(names)
        if base_ring is None:
            if poly_ring is not None:
                base_ring = poly_ring.base_ring()
        if names is None:
            if poly_ring is not None:
                names = poly_ring.gens()
        return poly_ring, base_ring, ambient_dim, names

    @classmethod
    def _default_name(cls, index):
        return "x{}".format(index)

    def __init__(self, base_ring=None, ambient_dim=None, **options):
        """
        Initialize a basic semialgebraic set as the universe in ``ambient_dim``.

        ``base_ring`` is the ring in which the coefficients of polynomials live.
        """
        poly_ring, base_ring, ambient_dim, names = self._poly_ring_from_options(
            base_ring=base_ring, ambient_dim=ambient_dim, **options)
        if not poly_ring:
            raise ValueError("poly_ring is not determined by the provided input")
        super(BasicSemialgebraicSet_base, self).__init__()
        self._ambient_dim = ambient_dim
        #self._ambient_space = base_ring ** ambient_dim
        self._base_ring = base_ring
        self._poly_ring = poly_ring

    # Default implementation. Subclasses are encouraged to provide
    # faster implementations.
    @classmethod
    def from_bsa(cls, bsa, base_ring=None, poly_ring=None, **init_kwds):
        r"""
        Initialize a basic semialgebraic set of class ``cls`` to be the same
        as ``bsa``.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: upstairs_bsa_ppl = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(ambient_dim=0)
            sage: Q.<x0,x1,x2> = QQ[]
            sage: veronese = BasicSemialgebraicSet_veronese(upstairs_bsa_ppl, [], v_dict=dict(), poly_ring=Q)
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

        Test that ``BasicSemialgebraicSet_veronese`` supports this method without keyword arguments::

            sage: veronese_bsa_again = BasicSemialgebraicSet_veronese.from_bsa(bsa_eq_lt_le_sets)
            sage: list(veronese_bsa_again.lt_poly())
            [54*x0^2 + 226*x1*x2 + 113]

        Test that ``BasicSemialgebraicSet_veronese`` supports this method and accepts keyword arguments::

            sage: new_upstairs_bsa = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(0)
            sage: veronese_bsa_again = BasicSemialgebraicSet_veronese.from_bsa(bsa_eq_lt_le_sets, upstairs_bsa=new_upstairs_bsa, polynomial_map=[], v_dict={}, poly_ring=bsa_eq_lt_le_sets.poly_ring())
            sage: list(veronese_bsa_again.lt_poly())
            [54*x0^2 + 226*x1*x2 + 113]

        If ``bsa`` is already of class ``cls`` and the correct ``base_ring``, it is just returned::

            sage: BasicSemialgebraicSet_eq_lt_le_sets.from_bsa(upstairs_bsa_eq_lt_le_sets) is upstairs_bsa_eq_lt_le_sets
            True

        We can use this method to upgrade a ``base_ring``::

            sage: upstairs_bsa_eq_lt_le_sets_AA = BasicSemialgebraicSet_eq_lt_le_sets.from_bsa(upstairs_bsa_eq_lt_le_sets, base_ring=AA)
            sage: upstairs_bsa_eq_lt_le_sets_AA is upstairs_bsa_eq_lt_le_sets
            False
            sage: upstairs_bsa_eq_lt_le_sets_AA.base_ring()
            Algebraic Real Field

        """
        if bsa.__class__ == cls and (base_ring is None or base_ring == bsa.base_ring()) and (poly_ring is None or poly_ring == bsa.poly_ring()) and ('polynomial_map' not in init_kwds):
            return bsa
        if poly_ring is None:
            if base_ring is not None:
                poly_ring, base_ring, ambient_dim, names = cls._poly_ring_from_options(
                    base_ring=base_ring,
                    ambient_dim=bsa.ambient_dim(),
                    names=bsa.poly_ring().gens())
            else:
                poly_ring = bsa.poly_ring()
        self = cls(poly_ring=poly_ring, **init_kwds)
        bsa.add_constraints_to(self)
        return self

    def add_constraints_to(self, bsa):
        """
        Add the constraints of ``self`` to ``bsa``.
        """
        for p in sorted(self.eq_poly()):
            bsa.add_polynomial_constraint(p, operator.eq)
        for p in sorted(self.lt_poly()):
            bsa.add_polynomial_constraint(p, operator.lt)
        for p in sorted(self.le_poly()):
            bsa.add_polynomial_constraint(p, operator.le)

    def ambient_dim(self):
        return self._ambient_dim

    def linear_forms_space(self):
        return self.base_ring() ** self.ambient_dim()

    def ambient_space(self, field=None):
        if field is None:
            field = QQ
        return field ** self.ambient_dim() 

    def base_ring(self):
        return self._base_ring

    def poly_ring(self):
        """
        Return the polynomial ring.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(2)
            sage: P.poly_ring()
            Multivariate Polynomial Ring in x0, x1 over Rational Field
        """
        return self._poly_ring

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
        r"""
        Whether the set contains the ``point`` (vector).
        """
        return all(f(*point) == 0 for f in self.eq_poly()) and all(f(*point) <= 0 for f in self.le_poly()) and all(f(*point) < 0 for f in self.lt_poly())

    @abstract_method
    def closure_polyhedron(self):
        r"""
        If the topological closure is a polyhedron, return it in
        the form of a Sage polyhedron.  Otherwise raise an error.
        """

    def formal_closure(self, bsa_class='formal_closure'):
        r"""
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
            BasicSemialgebraicSet_formal_closure(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(Constraint_System {-x0-x1+3>0}, names=[x0, x1]))
            sage: list(closure.eq_poly()), list(closure.lt_poly()), list(closure.le_poly())
            ([], [], [x0 + x1 - 3])
            sage: closure = bsa.formal_closure(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron); closure
            BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(Constraint_System {-x0-x1+3>=0}, names=[x0, x1])
            sage: list(closure.eq_poly()), list(closure.lt_poly()), list(closure.le_poly())
            ([], [], [x0 + x1 - 3])

        In general, this is only a superset of the topological closure::

            sage: R.<x> = QQ['x']
            sage: bsa = BasicSemialgebraicSet_eq_lt_le_sets(lt={x, -x})
            sage: formal_closure = bsa.formal_closure(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron)
            sage: list(formal_closure.eq_poly())
            [x]
            sage: bsa_ppl = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron.from_bsa(bsa)
            sage: closure = bsa_ppl.closure()
            sage: list(closure.eq_poly())
            [-1]

        """
        from .basic_semialgebraic_formal_closure import BasicSemialgebraicSet_formal_closure
        bsa_formal_closure = BasicSemialgebraicSet_formal_closure(self)
        bsa_class = _bsa_class(bsa_class)
        return bsa_class.from_bsa(bsa_formal_closure)

    def formal_relint(self, bsa_class='formal_relint'):
        r"""
        Return the basic semialgebraic set obtained by replacing all non-strict
        inequalities by strict inequalities.  This is a subset of the topological relative interior.

        By default, the formal relative interior is represented by an instance of class
        ``BasicSemialgebraicSet_formal_relint``; use the argument ``bsa_class``
        to choose another class.  See ``_bsa_class`` for the allowed class nicknames.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: bsa = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(ambient_dim=2)
            sage: bsa.add_linear_constraint([1, 1], -3, operator.le)
            sage: list(bsa.lt_poly()), list(bsa.le_poly())
            ([], [x0 + x1 - 3])
            sage: relint = bsa.formal_relint(); relint
            BasicSemialgebraicSet_formal_relint(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(Constraint_System {-x0-x1+3>=0}, names=[x0, x1]))
            sage: list(relint.eq_poly()), list(relint.lt_poly()), list(relint.le_poly())
            ([], [x0 + x1 - 3], [])
            sage: relint = bsa.formal_relint(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron); relint
            BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(Constraint_System {-x0-x1+3>0}, names=[x0, x1])
            sage: list(relint.eq_poly()), list(relint.lt_poly()), list(relint.le_poly())
            ([], [x0 + x1 - 3], [])

        In general, this is only a subset of the topological relative interior::

            sage: R.<x> = QQ['x']
            sage: bsa = BasicSemialgebraicSet_eq_lt_le_sets(le={x, -x})
            sage: formal_relint = bsa.formal_relint(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron)
            sage: list(formal_relint.eq_poly())
            [-1]
            sage: bsa_ppl = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron.from_bsa(bsa)
            sage: relint = bsa_ppl.relint()
            sage: list(relint.eq_poly())
            [x]

        """
        from .basic_semialgebraic_formal_relint import BasicSemialgebraicSet_formal_relint
        bsa_formal_relint = BasicSemialgebraicSet_formal_relint(self)
        bsa_class = _bsa_class(bsa_class)
        return bsa_class.from_bsa(bsa_formal_relint)

    @abstract_method
    def closure(self, bsa_class):
        r"""
        Return the basic semialgebraic set that is the topological closure
        of ``self``.
        """

    @abstract_method
    def relint(self, bsa_class):
        r"""
        Return the basic semialgebraic set that is the topological relative interior
        of ``self``.
        """

    def intersection(self, *bsa_list, **kwds):
        r"""
        Return the basic semialgebraic set that is the intersection of ``self``
        with the basic semialgebraic sets in ``bsa_list``.

        By default, the intersection is represented by an instance of class
        ``BasicSemialgebraicSet_intersection``; use the argument ``bsa_class``
        to choose another class.  See ``_bsa_class`` for the allowed class nicknames.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P.<x,y,z> = QQ[]
            sage: bsa1 = BasicSemialgebraicSet_eq_lt_le_sets(le=[-x, -y])
            sage: bsa2 = BasicSemialgebraicSet_eq_lt_le_sets(eq=[x + y], le=[-z])
            sage: bsa3 = BasicSemialgebraicSet_eq_lt_le_sets(lt=[y + z])
            sage: bsa123 = bsa1.intersection(bsa2, bsa3)
            sage: sorted(bsa123.le_poly()), sorted(bsa123.eq_poly()), sorted(bsa123.lt_poly())
            ([-z, -y, -x], [x + y], [y + z])
        """
        from .basic_semialgebraic_intersection import BasicSemialgebraicSet_intersection
        bsa_class = _bsa_class(kwds.pop('bsa_class', 'intersection'))
        bsa_intersection = BasicSemialgebraicSet_intersection([self] + list(bsa_list), **kwds)
        return bsa_class.from_bsa(bsa_intersection)

    def add_space_dimensions_and_embed(self, space_dim_to_add=None, names=None):
        r"""
        Mutate ``self`` by injecting it into a higher dimensional space.
        """
        old_dim = self._ambient_dim
        self._ambient_dim += space_dim_to_add
        new_dim = self._ambient_dim
        if names is None:
            if space_dim_to_add is None:
                raise ValueError("need to provide space_dim_to_add or names")
            names = [ self._default_name(i) for i in range(old_dim, new_dim) ]
        else:
            names = list(names)
            if len(names) != space_dim_to_add:
                raise ValueError("space_dim_to_add and len(names) do not match")
        self._poly_ring = PolynomialRing(self.base_ring(), names=list(self.poly_ring().gens()) + names)
        #self._ambient_space = self._ambient_space.base_ring() ** new_dim

    def _linear_polynomial(self, form, constant=0):
        return sum(coeff * gen for coeff, gen in zip(form, self.poly_ring().gens())) + constant

    def is_empty(self):
        """
        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: from cutgeneratingfunctionology.spam.semialgebraic_mathematica import BasicSemialgebraicSet_mathematica
            sage: P.<x,y> = QQ[]
            sage: BasicSemialgebraicSet_eq_lt_le_sets(QQ, poly_ring=P).is_empty()
            False
            sage: BasicSemialgebraicSet_eq_lt_le_sets(lt=[x-x]).is_empty()
            True
            sage: BasicSemialgebraicSet_mathematica(lt=[x-x]).is_empty()   # optional - mathematica
            True
            sage: BasicSemialgebraicSet_eq_lt_le_sets(eq=[y-y], lt=[x-x-1], le=[x-x]).is_empty()
            Traceback (most recent call last):
            ...
            NotImplementedError...
            sage: BasicSemialgebraicSet_mathematica(eq=[y-y], lt=[x-x-1], le=[x-x]).is_empty()    # optional - mathematica
            False
            sage: BasicSemialgebraicSet_eq_lt_le_sets(lt=[x^2+y^2]).is_empty()
            Traceback (most recent call last):
            ...
            NotImplementedError...
            sage: BasicSemialgebraicSet_mathematica(lt=[x^2+y^2]).is_empty()  # optional - mathematica
            True
        """
        if not self.eq_poly() and not self.lt_poly() and not self.le_poly():
            return False
        for l in self.eq_poly():
            if l in self.base_ring() and l != 0:
                return True
        for l in self.lt_poly():
            if l in self.base_ring() and l >= 0:
                return True
        for l in self.le_poly():
            if l in self.base_ring() and l > 0:
                return True
        raise NotImplementedError

    def is_universe(self):
        """
        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: from cutgeneratingfunctionology.spam.semialgebraic_mathematica import BasicSemialgebraicSet_mathematica
            sage: P.<x,y> = QQ[]
            sage: BasicSemialgebraicSet_eq_lt_le_sets(QQ, 2).is_universe()
            True
            sage: BasicSemialgebraicSet_eq_lt_le_sets(lt=[x-x]).is_universe()
            Traceback (most recent call last):
            ...
            NotImplementedError...
            sage: BasicSemialgebraicSet_mathematica(lt=[x-x]).is_universe()   # optional - mathematica
            False
            sage: BasicSemialgebraicSet_eq_lt_le_sets(eq=[y-y], lt=[x-x-1], le=[x-x]).is_universe()
            True
            sage: BasicSemialgebraicSet_eq_lt_le_sets(le=[-x^2-y^2]).is_universe()
            Traceback (most recent call last):
            ...
            NotImplementedError...
            sage: BasicSemialgebraicSet_mathematica(le=[-x^2-y^2]).is_universe()  # optional - mathematica
            True
        """
        if not self.eq_poly() and not self.lt_poly() and not self.le_poly():
            return True
        for l in self.eq_poly():
            if l != 0:
                raise NotImplementedError
        for l in self.lt_poly():
            if not l in self.base_ring() or not l < 0:
                raise NotImplementedError
        for l in self.le_poly():
            if not l in self.base_ring() or not l <= 0:
                raise NotImplementedError
        return True

    def add_linear_constraint(self, lhs, cst, op):
        r"""
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
        self.add_polynomial_constraint(self._linear_polynomial(lhs, cst), op)

    def is_linear_constraint_valid(self, lhs, cst, op):
        r"""
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
        r"""
        ``lhs`` should be a polynomial.
        Add the constraint ``lhs``(x) ``op`` 0,
        where ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.
        """

    @abstract_method
    def is_polynomial_constraint_valid(self, lhs, op):
        r"""
        Check if the constraint ``lhs``(x) ``op`` 0 is satisfied
        for all points of ``self``, where ``lhs`` is a polynomial, and
        ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.

        Raise an error if the information provided does not suffice to decide
        the validity of the constraint.
        """

    def linear_function_upper_bound(self, form):
        r"""
        Find an upper bound for ``form`` (a vector) on ``self``.

        The default implementation just returns +oo.
        Subclasses should provide more useful bounds.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: bsa = BasicSemialgebraicSet_base(QQ, 2)
            sage: bsa.linear_function_upper_bound(vector(QQ, [1, 1]))
            +Infinity

        """
        return self.polynomial_function_upper_bound(self._linear_polynomial(form))

    def linear_function_lower_bound(self, form):
        r"""
        Find a lower bound for ``form`` (a vector) on ``self``.

        This implementation delegates to ``linear_function_upper_bound``.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: bsa = BasicSemialgebraicSet_base(QQ, 2)
            sage: bsa.linear_function_lower_bound(vector(QQ, [1, 1]))
            -Infinity

        """
        form = vector(form)
        return -(self.linear_function_upper_bound(-form))

    def polynomial_function_upper_bound(self, polynomial):
        r"""
        Find an upper bound for ``polynomial`` on ``self``.

        The default implementation just returns +oo for nonconstant polynomials.
        Subclasses should provide more useful bounds.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: bsa = BasicSemialgebraicSet_base(QQ, 2)
            sage: bsa.polynomial_function_upper_bound(0)
            0
            sage: bsa.polynomial_function_upper_bound(1)
            1
            sage: bsa.polynomial_function_upper_bound(bsa.poly_ring().gen(0))
            +Infinity
        """
        polynomial = self.poly_ring()(polynomial)
        if polynomial.degree() <= 0:
            return polynomial.constant_coefficient()
        else:
            return +Infinity

    def polynomial_function_lower_bound(self, polynomial):
        r"""
        Find a lower bound for ``polynomial`` on ``self``.

        The default implementation just returns -oo for nonconstant polynomials.
        Subclasses should provide more useful bounds.
        """
        return -(self.polynomial_function_upper_bound(-polynomial))

    @abstract_method
    def find_point(self):
        r"""
        Find a point in ``self``.
        """
        # default implementation could go through self.closure_polyhedron()

    @abstract_method
    def coordinate_projection(self, coordinates, bsa_class='projection'):
        r"""
        Compute the projection to ``coordinates`` (a list or tuple of indices or
        variables of ``self.poly_ring``).

        This is a semialgebraic set, but in general not a basic semialgebraic set.
        """

    def section(self, section_polynomial_map, bsa_class='section', **kwds):
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

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P_upstairs.<x,y> = QQ[]
            sage: P_downstairs.<t> = QQ[]
            sage: polynomial_map = [75/19*t, t]
            sage: bsa = BasicSemialgebraicSet_eq_lt_le_sets(lt=[21*x - 8, -x, 950*x^2 - 3700*x*y - 225*y^2 - 133*x])
            sage: bsa_section_1 = bsa.section(polynomial_map)
            sage: sorted(bsa_section_1.lt_poly())
            [-75/19*t, 1575/19*t - 8, -525/19*t^2 - 525*t]
            sage: bsa_section_2 = bsa.section(polynomial_map, bsa_class='veronese')
            sage: sorted(bsa_section_2.lt_poly())
            [-t, 1575*t - 152, -t^2 - 19*t]

            sage: P.<x,y,z> = QQ[]
            sage: upstairs_bsa = BasicSemialgebraicSet_eq_lt_le_sets(le=[x*y + y*z])
            sage: polynomial_map = [x, x, z]
            sage: downstairs_bsa_section = upstairs_bsa.section(polynomial_map)
            sage: downstairs_bsa_veronese = upstairs_bsa.section(polynomial_map, bsa_class='veronese')

        """
        bsa_section = BasicSemialgebraicSet_section(self, section_polynomial_map)
        bsa_class = _bsa_class(bsa_class)
        return bsa_class.from_bsa(bsa_section, **kwds)

    def plot(self, alpha=0.5, plot_points=300, slice_value=None, color='blue', fill_color='blue',
             constraints_color=None, constraints_fill_color=None,
             eq_constraints_kwds={'linestyle': 'dotted' },
             le_constraints_kwds={'linestyle': 'solid' },
             lt_constraints_kwds={'linestyle': 'dashed' },
             **kwds):
        r"""
        Plot the semialgebraic set or a slice (section) of it.

        - If slice_value is given, it is either a polynomial_map that defines a section, or a list of fixed parameter values with two of them being None. Plot the section.
        - plot_points controls the quality of the plotting.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: bsa = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(ambient_dim=2)
            sage: bsa.add_linear_constraint([1, 1], -1, operator.lt)
            sage: bsa.plot()                                         # not tested

         Plot the slice in (x,y)-space with z=4::

            sage: P.<x,y,z> = QQ[]
            sage: bsa = BasicSemialgebraicSet_eq_lt_le_sets(le=[-x, -y, -z, x+y-1, z-6])
            sage: bsa.plot(slice_value=[None, None, 4])              # not tested

        Plot the section with x=y::

            sage: Q.<u,v> = QQ[]
            sage: bsa.plot(slice_value=[u, u, v], xmin=-1, xmax=2, ymin=-1, ymax=8) # not tested

        Plot by showing all constraints::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P.<x,y> = QQ[]
            sage: bsa = BasicSemialgebraicSet_eq_lt_le_sets(le=[-x, -y, x+y-1, x-2], lt=[x^2 + y^2 - 5])
            sage: kwds = {'xmin': -3, 'xmax': 3, 'ymin': -3, 'ymax': 3}
            sage: bsa.plot(**kwds)   # not tested
            sage: bsa.plot(color=None, fill_color=None, constraints_color='black', constraints_fill_color='black', **kwds)   # not tested
            sage: bsa.plot(constraints_color='black', constraints_fill_color='black', **kwds)   # not tested

            sage: P.<x,y> = QQ[]
            sage: bsa = BasicSemialgebraicSet_eq_lt_le_sets(le=[-x, -y, x-2], eq=[x+y-1], lt=[x^2 + y^2 - 5])
            sage: bsa.plot(**kwds)   # not tested
            sage: bsa.plot(color='red', constraints_color='black', constraints_fill_color='black', eq_constraints_kwds={'linestyle': 'dotted', 'thickness': 3}, **kwds)   # not tested

        """
        ## Refactor SemialgebraicComplexComponent.plot and plot2dslice through this method.
        if (not slice_value) and (self.ambient_dim() != 2):
            raise NotImplementedError("Plotting with dimension not equal to 2 is not implemented. Try `slice_value` to plot a slice.")
        if slice_value:
            if slice_value.count(None) == 2:
                PM = QQ['x','y']; (x, y) = PM.gens()  # invalid syntax PM.<x,y>=QQ[]
                i = slice_value.index(None);
                j = slice_value.index(None, i+1);
                polynomial_map = [PM(pm) for pm in slice_value]
                polynomial_map[i] = x
                polynomial_map[j] = y
                section = self.section(polynomial_map)
            elif slice_value.count(None) == 0: # slice_value is a polynomial_map
                section = self.section(slice_value)
            else:
                raise NotImplementedError("Plotting with dimension not equal to 2 is not implemented.")
            return section.plot(alpha=alpha, plot_points=plot_points, slice_value=None,
                                color=color, fill_color=fill_color,
                                constraints_color=constraints_color,
                                constraints_fill_color=constraints_fill_color,
                                eq_constraints_kwds=eq_constraints_kwds,
                                le_constraints_kwds=le_constraints_kwds,
                                lt_constraints_kwds=lt_constraints_kwds, **kwds)
        g = Graphics()
        for bv in 'xmin', 'xmax', 'ymin', 'ymax':
            b = kwds.get(bv, None)
            if b is not None:
                getattr(g, bv)(b)
            else:
                kwds.pop(bv, None) # remove xmin=None from kwds, so that the xmin etc. below work.
        xmin = max(self.linear_function_lower_bound([1,0]), kwds.get('xmin', -Infinity))
        xmax = min(self.linear_function_upper_bound([1,0]), kwds.get('xmax', +Infinity))
        ymin = max(self.linear_function_lower_bound([0,1]), kwds.get('ymin', -Infinity))
        ymax = min(self.linear_function_upper_bound([0,1]), kwds.get('ymax', +Infinity))
        if (xmin > xmax) or (ymin > ymax):
            return g
        if xmin is -Infinity:
            xmin = 0
        if xmax is +Infinity:
            xmax = 1
        if ymin is -Infinity:
            ymin = 0
        if ymax is +Infinity:
            ymax = 1
        if (xmin > xmax) or (ymin > ymax):
            raise ValueError("Please provide bounds for plotting.")
        from sage.symbolic.ring import SR
        x, y = SR.var('x, y')
        x_range = (x, xmin-0.01, xmax+0.01)
        y_range = (y, ymin-0.01, ymax+0.01)
        eq_constraints = [ l(x, y) == 0 for l in self.eq_poly() ]
        lt_constraints = [ l(x, y) <  0 for l in self.lt_poly() ]
        le_constraints = [ l(x, y) <= 0 for l in self.le_poly() ]
        if constraints_fill_color:
            for c in chain(lt_constraints, le_constraints):
                if not isinstance(c, bool):
                    g += region_plot([c], x_range, y_range,
                                     alpha=0.1, plot_points=plot_points,
                                     incol=constraints_fill_color, bordercol=None, **kwds)
        if color or fill_color:
            all_constraints = eq_constraints + lt_constraints + le_constraints
            non_trivial_constraints = [ c for c in all_constraints if not isinstance(c, bool) ]
            if not any(c is False for c in all_constraints):
                g += region_plot(non_trivial_constraints, x_range, y_range,
                                 alpha=alpha, plot_points=plot_points,
                                 incol=fill_color, bordercol=color, **kwds)
        if constraints_color:
            for polys, plot_kwds in [(self.eq_poly(), eq_constraints_kwds),
                                     (self.lt_poly(), lt_constraints_kwds),
                                     (self.le_poly(), le_constraints_kwds)]:
                for l in polys:
                    c = l(x, y) == 0
                    if not isinstance(c, bool):
                        effective_kwds = copy(kwds)
                        effective_kwds['color'] = constraints_color
                        effective_kwds.update(plot_kwds)
                        g += implicit_plot(c, x_range, y_range, **effective_kwds)
        return g

class BasicSemialgebraicSet_polyhedral(BasicSemialgebraicSet_base):

    """
    An abstract class of polyhedral basic semialgebraic sets.

    """

    @abstract_method
    def add_linear_constraint(self, lhs, cst, op):
        r"""
        Add the constraint ``lhs`` * x + cst ``op`` 0,
        where ``lhs`` is a vector of length ``ambient_dim`` and
        ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.

        In subclasses of ``BasicSemialgebraicSet_polyhedral``,
        this method should be defined.
        """

    def add_polynomial_constraint(self, lhs, op):
        r"""
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
        lhs_vector = vector(lhs.monomial_coefficient(x) for x in self.poly_ring().gens())
        # univariate polynomials (Polynomial_rational_flint) does not define "coefficient", but has "monomial_coefficient".
        self.add_linear_constraint(lhs_vector, cst, op)

    def is_polynomial_constraint_valid(self, lhs, op):
        r"""
        Check if the constraint ``lhs``(x) ``op`` 0 is satisfied
        for all points of ``self``, where ``lhs`` is a polynomial, and
        ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.

        This implementation checks that ``lhs`` is a linear polynomial
        and then delegates to ``add_linear_constraint``.

        Raise an error if the information provided does not suffice to decide
        the validity of the constraint.
        """
        lhs = self.poly_ring()(lhs)   # convert if necessary
        if lhs.degree() > 1:
            raise NotImplementedError("{} is not a valid linear polynomial.".format(lhs))
        cst = lhs.constant_coefficient()
        lhs_vector = vector(lhs.monomial_coefficient(x) for x in self.poly_ring().gens())
        # univariate polynomials (Polynomial_rational_flint) does not define "coefficient", but has "monomial_coefficient".
        return self.is_linear_constraint_valid(lhs_vector, cst, op)

    def polynomial_function_upper_bound(self, polynomial):
        r"""
        Find an upper bound for ``polynomial`` on ``self``.

        This implementation checks that ``polynomial`` is a linear polynomial
        and then delegates to ``linear_function_upper_bound``.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(1)
            sage: f = 2 * P.poly_ring().gen() + 1
            sage: P.add_polynomial_constraint(f, operator.le)
            sage: P.polynomial_function_upper_bound(4*P.poly_ring().gen()+1)
            -1
        """
        polynomial = self.poly_ring()(polynomial)   # convert if necessary
        if polynomial.degree() > 1:
            raise NotImplementedError("{} is not a valid linear polynomial.".format(lhs))
        cst = polynomial.constant_coefficient()
        lhs_vector = vector(polynomial.monomial_coefficient(x) for x in self.poly_ring().gens())
        # univariate polynomials (Polynomial_rational_flint) does not define "coefficient", but has "monomial_coefficient".
        return self.linear_function_upper_bound(lhs_vector) + cst

## (1) In the first step, we implement the following class.  Everything is linear.
## Rewrite all direct uses of PPL in ParametricRealFieldElement, ParametricRealField
## using method calls to this class.

poly_is_included = Poly_Con_Relation.is_included()
#strictly_intersects = Poly_Con_Relation.strictly_intersects()
point_is_included = Poly_Gen_Relation.subsumes()
#con_saturates = Poly_Con_Relation.saturates()

from sage.arith.functions import lcm

class BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(BasicSemialgebraicSet_polyhedral):

    r"""
    A (possibly half-open) polyhedral basic semialgebraic set,
    represented by a PPL ``NNC_Polyhedron``

    """

    def __init__(self, ambient_dim=None, polyhedron=None, base_ring=None, poly_ring=None, **options):
        r"""
        Initialize a basic semialgebraic set as the universe in
        ``ambient_dim``, or, if ``polyhedron`` (an ``NNC_Polyhedron``,
        which after that belongs to this object) is provided, as
        that.
        
        TEST::
        
            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(2)
            sage: P.add_linear_constraint([0,1],0,operator.ge)
            sage: P.add_linear_constraint([1,0],0,operator.ge)
            sage: P.add_linear_constraint([2,3],-6,operator.lt)
            sage: P
            BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(Constraint_System {-2*x0-3*x1+6>0, x1>=0, x0>=0}, names=[x0, x1])
            sage: sorted(P.eq_poly())
            []
            sage: sorted(P.lt_poly())
            [2*x0 + 3*x1 - 6]
            sage: sorted(P.le_poly())
            [-x1, -x0]
        """
        if ambient_dim is None and polyhedron is not None:
            ambient_dim = polyhedron.space_dimension()
        if base_ring is None and poly_ring is None:
            base_ring = QQ
        poly_ring, base_ring, ambient_dim, names = self._poly_ring_from_options(
            ambient_dim=ambient_dim, base_ring=base_ring, poly_ring=poly_ring, **options)
        if base_ring is not QQ:
            raise ValueError("only base_ring=QQ is supported")
        super(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron, self).__init__(poly_ring=poly_ring)
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
        return self.__class__(polyhedron=NNC_Polyhedron(self._polyhedron), poly_ring=self.poly_ring())

    def _repr_(self):
        constraints = self._polyhedron.minimized_constraints()
        names = list(self.poly_ring().gens())
        return 'BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron({}, names={})'.format(
            constraints, names)

    def closure(self, bsa_class='formal_closure'):
        r"""
        Return the basic semialgebraic set that is the topological closure
        of ``self``.
        """
        # Because our description consists of minimized constraints, the closure is
        # just the formal closure.
        return self.formal_closure(bsa_class=bsa_class)

    def relint(self, bsa_class='formal_relint'):
        r"""
        Return the basic semialgebraic set that is the topological relative interior
        of ``self``.
        """
        # Because our description consists of minimized constraints, the relint is
        # just the formal relint.
        return self.formal_relint(bsa_class=bsa_class)

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
        r"""
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
        r"""
        Find a point in ``self``.
        
        EXAMPLES::
        
            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(2)
            sage: P.add_linear_constraint([0,1],0,operator.ge)
            sage: P.add_linear_constraint([1,0],0,operator.ge)
            sage: P.add_linear_constraint([2,3],-6,operator.lt)
            sage: P.find_point()
            (11/10, 11/15)
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
        r"""
        Mutate ``self`` by injecting it into a higher dimensional space.
        """
        super(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron, self).add_space_dimensions_and_embed(space_dim_to_add)
        self._polyhedron.add_space_dimensions_and_embed(space_dim_to_add)

    @staticmethod
    def _ppl_constraint(lhs, cst, op):
        r"""
        Make a PPL ``Constraint`` ``lhs`` * x + cst ``op`` 0,
        where ``lhs`` is be a vector of length ambient_dim.
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

    def linear_function_upper_bound(self, form):
        r"""
        Find an upper bound for ``form`` (a vector) on ``self``.
        This upper bound is the supremum.

        If ``self`` is empty, it returns -oo
        """
        def to_point(g):
            den = g.divisor()
            return vector(QQ, ( QQ(x)/den for x in g.coefficients() ))
        def to_vector(g):
            return vector(QQ, ( QQ(x) for x in g.coefficients() ))
        if self._polyhedron.is_empty():
            return -Infinity
        form = vector(form)
        for g in self._polyhedron.generators():
            if g.is_line():
                if to_vector(g) * form != 0:
                    return +Infinity
            if g.is_ray():
                if to_vector(g) * form > 0:
                    return +Infinity
        points = [ to_point(g) for g in self._polyhedron.generators()
                   if g.is_point() or g.is_closure_point() ]
        return max(p * form for p in points)

    def is_linear_constraint_valid(self, lhs, cst, op):
        r"""
        Whether the constraint ``lhs`` * x + cst ``op`` 0
        is satisfied for all points of ``self``,
        where ``lhs`` is be a vector of length ambient_dim.
        
        EXAMPLES::
        
            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(2)
            sage: P.add_linear_constraint([0,1],0,operator.ge)
            sage: P.add_linear_constraint([1,0],0,operator.ge)
            sage: P.add_linear_constraint([2,3],-6,operator.lt)
            sage: P.is_linear_constraint_valid([1,1],-3,operator.lt)
            True
            sage: P.is_linear_constraint_valid([0,1],0,operator.gt)
            False
        """
        lhs = vector(lhs)
        constraint = self._ppl_constraint(lhs, cst, op)
        return self._polyhedron.relation_with(constraint).implies(poly_is_included)

    def add_linear_constraint(self, lhs, cst, op):
        r"""
        Add the constraint ``lhs`` * x + cst ``op`` 0,
        where ``lhs`` is a vector of length ambient_dim, and
        ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``
        
        EXAMPLES::
        
            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(2)
            sage: P.add_linear_constraint([2,3],-6,operator.gt)
            sage: sorted(P.lt_poly())
            [-2*x0 - 3*x1 + 6]
        """
        lhs = vector(lhs)
        constraint = self._ppl_constraint(lhs, cst, op)
        self._polyhedron.add_constraint(constraint)

    def is_empty(self):
        """
        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(1)
            sage: S.add_linear_constraint([1], -1, operator.ge)
            sage: S.is_empty()
            False
            sage: S.add_linear_constraint([1], +1, operator.le)
            sage: S.is_empty()
            True
        """
        return self._polyhedron.is_empty()

    def is_universe(self):
        """
        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(1)
            sage: S.add_linear_constraint([0], 0, operator.eq)
            sage: S.is_universe()
            True
            sage: S.add_linear_constraint([1], +1, operator.le)
            sage: S.is_universe()
            False
        """
        return self._polyhedron.is_universe()

## class BasicSemialgebraicSet_polyhedral_ppl_MIP_Problem(BasicSemialgebraicSet_base):

##     """
##     A closed polyhedral basic semialgebraic set,
##     represented by a PPL ``MIP_Problem``.
##     """

### do we need the above, or should we just use the below, using solver='ppl'?


## (2) Then
class BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram(BasicSemialgebraicSet_polyhedral):

    r"""
    A closed polyhedral basic semialgebraic set,
    represented by a Sage ``MixedIntegerLinearProgram``
    with only continuous variables.
    """

    def __init__(self, base_ring, ambient_dim, poly_ring=None, names=None, solver=None):
        r"""
        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram(QQ, 3)

        """
        super(BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram, self).__init__(
            base_ring=base_ring, ambient_dim=ambient_dim, poly_ring=poly_ring, names=names)
        self._mip = MixedIntegerLinearProgram(solver=solver, maximization=True)

    def __copy__(self):
        raise NotImplementedError()

    def mip(self):
        return self._mip

    def mip_gens(self):
        r"""
        Return the components of the MIP variable corresponding to the
        space dimensions.
        """
        mip_var = self.mip().default_variable()
        for i in range(self.ambient_dim()):
            yield mip_var[i]

    def closure(self, bsa_class='mip'):
        r"""
        Return the basic semialgebraic set that is the topological closure
        of ``self``, which is ``self`` itself.
        """
        return self

    def is_linear_constraint_valid(self, lhs, cst, op):
        r"""
        Whether the constraint ``lhs`` * x + cst ``op`` 0
        is satisfied for all points of ``self``.

        Raise ValueError if ``lhs`` is not of the right dimension.

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
            # The default implementation checked
            # ``linear_function_upper_bound`` and
            # ``linear_function_lower_bound`` already.  These are the sup and
            # inf in our case, and because we are a closed polyhedral set, this
            # suffices to decide validity.
            return False

    def _mip_linear_function(self, form, cst=0):
        r"""
        Obtain a linear function object representing
        ``form`` * x + ``cst``.
        Raise ValueError if ``form`` is not of the right dimension.
        """
        if not len(form) ==  self.ambient_dim():
            raise ValueError("wrong dimension")
        return self.mip().sum(coeff * gen for coeff, gen in zip(form, self.mip_gens())) + cst

    def linear_function_upper_bound(self, form):
        r"""
        Find an upper bound for ``form`` (a vector) on ``self``.

        In this implementation, this is done by solving the LP,
        so this upper bound is the supremum.

        If ``self`` is empty, this returns -oo

        Raise ValueError if ``form`` is not of the right dimension.
        """
        mip = self.mip()
        objective = self._mip_linear_function(form)
        mip.set_objective(objective)
        try:
            return mip.solve()
        except MIPSolverException:
            # Here we distinguish infeasible and unbounded.
            mip_copy = copy(mip)
            mip_copy.set_objective(1)
            try:
                mip_copy.solve()
                return +Infinity
            except MIPSolverException:
                return -Infinity

    def add_linear_constraint(self, lhs, cst, op):
        r"""
        Add the constraint ``lhs`` * x + cst ``op`` 0,
        where ``lhs`` is a vector of length ``ambient_dim`` and
        ``op`` is one of ``operator.eq``, ``operator.le``, ``operator.ge``.

        Raise ValueError if ``lhs`` is not of the right dimension.

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

    def _repr_(self):
        return 'BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram with {}'.format(self._mip)

    def eq_poly(self):
        r"""
        Generate the polynomials `f` in equations `f(x) = 0`
        in the description of ``self``.

        Together, ``eq_poly`` and ``le_poly`` describe ``self``.
        (``lt_poly`` gives the empty list because ``self`` is closed.)

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram(QQ, 1, solver='ppl')
            sage: S.add_linear_constraint([1], -1, operator.eq)
            sage: sorted(S.eq_poly())
            [x0 - 1]
        """
        poly_ring = self.poly_ring()
        for lb, (indices, coeffs), ub in self.mip().constraints():
            if lb == ub:
                form = sum(coeff * poly_ring.gen(i) for i, coeff in zip(indices, coeffs))
                yield form - ub

    def le_poly(self):
        r"""
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
            [-x0, x0 - 1]
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
        r"""
        Return the empty list because ``self`` is closed.

        Together, ``eq_poly`` and ``le_poly`` describe ``self``.
        """
        return []

    def is_empty(self):
        """
        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram(QQ, 1, solver='ppl')
            sage: S.add_linear_constraint([1], -1, operator.ge)
            sage: S.is_empty()
            False
            sage: S.add_linear_constraint([1], +1, operator.le)
            sage: S.is_empty()
            True
        """
        try:
            self.mip().solve()
            return False
        except MIPSolverException:
            return True

    def is_universe(self):
        """
        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram(QQ, 1, solver='ppl')
            sage: S.add_linear_constraint([0], 0, operator.eq)
            sage: S.is_universe()
            True
            sage: S.add_linear_constraint([1], +1, operator.le)
            sage: S.is_universe()
            False
        """
        try:
            return super(BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram, self).is_universe()
        except NotImplementedError:
            return False

## (3) Then introduce the following class to simplify the code in parametric.sage
class BasicSemialgebraicSet_eq_lt_le_sets(BasicSemialgebraicSet_base):

    r"""
    A basic semialgebraic set, represented in a straightforward way
    as 3 finite sets of polynomial constraints `p(x) OP 0`.

    It does not remember the order of the polynomials provided at initialization.
    """

    def __init__(self, base_ring=None, ambient_dim=None, poly_ring=None, names=None, eq=[], lt=[], le=[]):
        r"""
        Construct a basic semialgebraic set.

        INPUT::

        - ``base_ring`` - the field in which the coefficients of the polynomials live.

        - ``ambient_dim`` - the dimension of the affine space.

        - ``poly_ring`` - a polynomial ring over ``base_ring`` in ``ambient_dim` variables.

        - ``eq`` - an iterable of finitely many polynomials `p(x)` in ``poly_ring`` that describe
          equations `p(x) = 0`.

        - ``lt`` - an iterable of finitely many polynomials `p(x)` in ``poly_ring`` that describe
          strict inequalities `p(x) < 0`.

        - ``le`` - an iterable of finitely many polynomials `p(x)` in ``poly_ring`` that describe
          inequalities `p(x) \leq 0`.

        If any polynomial is provided, ``poly_ring`` is inferred etc.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_eq_lt_le_sets(QQ, 3)
            sage: Q.<x0,x1,x2> = QQ[]
            sage: lhs = 27/113 * x0^2 + x1*x2 + 1/2
            sage: S.add_polynomial_constraint(lhs,operator.lt)
            sage: lhs = 27*x0 + 2*x1 + 1
            sage: S.add_polynomial_constraint(lhs,operator.gt)
            sage: lhs = 27*x1 + 2*x1*x0 - 1/2
            sage: S.add_polynomial_constraint(lhs,operator.eq)
            sage: lhs = x1^3 + x0
            sage: S.add_polynomial_constraint(lhs,operator.le)
            sage: S
            BasicSemialgebraicSet_eq_lt_le_sets(eq=[2*x0*x1 + 27*x1 - 1/2], lt=[-27*x0 - 2*x1 - 1, 27/113*x0^2 + x1*x2 + 1/2], le=[x1^3 + x0])
            sage: S.poly_ring()
            Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
            sage: S.eq_poly()
            {2*x0*x1 + 27*x1 - 1/2}
            sage: S.lt_poly()
            {-27*x0 - 2*x1 - 1, 27/113*x0^2 + x1*x2 + 1/2}
            sage: S.le_poly()
            {x1^3 + x0}
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
        super(BasicSemialgebraicSet_eq_lt_le_sets, self).__init__(base_ring=base_ring, ambient_dim=ambient_dim, poly_ring=poly_ring, names=names)
        self._eq = set(eq)
        self._lt = set(lt)
        self._le = set(le)

    def __copy__(self):
        r"""
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
        return self.__class__(base_ring=self.base_ring(), ambient_dim=self.ambient_dim(),
                              poly_ring=self.poly_ring(), eq=self._eq, lt=self._lt, le=self._le)

    # override the abstract methods

    def _repr_(self):
        return 'BasicSemialgebraicSet_eq_lt_le_sets(eq={}, lt={}, le={})'.format(sorted(self._eq), sorted(self._lt), sorted(self._le))

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

    def is_polynomial_constraint_valid(self, lhs, op):
        r"""
        Return True if the constraint ``lhs``(x) ``op`` 0 is a known valid
        inequality for ``self``. Return ``NotImplementedError`` otherwise.
        Input ``lhs`` is a polynomial, and ``op`` is one of ``operator.lt``,
        ``operator.gt``, ``operator.eq``,``operator.le``, ``operator.ge``.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_eq_lt_le_sets(QQ, 3)
            sage: Q.<x0,x1,x2> = QQ[]
            sage: lhs = 27/113 * x0^2 + x1*x2 + 1/2
            sage: S.add_polynomial_constraint(lhs,operator.lt)
            sage: S.is_polynomial_constraint_valid(lhs, operator.lt)
            True
            sage: S.is_polynomial_constraint_valid(lhs, operator.le)
            True
            sage: S.is_polynomial_constraint_valid(-lhs, operator.ge)
            True
            sage: S.is_polynomial_constraint_valid(lhs, operator.eq) #This equation is False but the code does not know.
            Traceback (most recent call last):
            ...
            NotImplementedError...
        """
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
        raise NotImplementedError

    def add_polynomial_constraint(self, lhs, op):
        r"""
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

    r"""
    Section of another ``BasicSemialgebraicSet``.

    See ``BasicSemialgebraicSet_base.section``.
    """

    def __init__(self, upstairs_bsa, polynomial_map, poly_ring=None, ambient_dim=None, base_ring=None, names=None):
        r"""
        EXAMPLES:

        McCormick::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: PY.<x,y,xy> = QQ[]
            sage: bsa = BasicSemialgebraicSet_eq_lt_le_sets(le=[-x, -y, -xy, x-2, y-3, xy-6])
            sage: PX.<u,v> = QQ[]
            sage: F = [u, v, u*v]
            sage: section = bsa.section(F); section
            BasicSemialgebraicSet_section(BasicSemialgebraicSet_eq_lt_le_sets(eq=[], lt=[], le=[-xy, xy - 6, -y, y - 3, -x, x - 2]), polynomial_map=[u, v, u*v])
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

        TESTS:

        Coercion to a common parent is applied to determine the effective ``poly_ring``::

            sage: PX_ZZ = ZZ['u', 'v']
            sage: u_zz, v_zz = PX_ZZ.gens()
            sage: F = [ZZ(1), v_zz, u]
            sage: section = bsa.section(F)
            sage: poly_ring = section.poly_ring(); poly_ring
            Multivariate Polynomial Ring in u, v over Rational Field
            sage: all(f.parent() is poly_ring for f in section.polynomial_map())
            True
        """
        if len(polynomial_map) != upstairs_bsa.ambient_dim():
            raise ValueError("polynomial_map must have the same length as the ambient dimension of the upstairs bsa")
        base_ring = upstairs_bsa.base_ring()
        if poly_ring is None:
            cm = get_coercion_model()
            poly_ring = cm.common_parent(*polynomial_map)
            if not is_PolynomialRing(poly_ring) and not is_MPolynomialRing(poly_ring):
                poly_ring = PolynomialRing(base_ring, [])
        polynomial_map = [ poly_ring(f) for f in polynomial_map ]
        if poly_ring and ambient_dim is None:
            ambient_dim = poly_ring.ngens()
        if not all(poly.parent().ngens() == ambient_dim for poly in polynomial_map):
            raise ValueError("elements in polynomial_map must come from a polynomial ring with the same number of variables as the ambient dimension")
        super(BasicSemialgebraicSet_section, self).__init__(base_ring=base_ring, ambient_dim=ambient_dim, names=names, poly_ring=poly_ring)
        self._upstairs_bsa = upstairs_bsa
        self._polynomial_map = polynomial_map

    def __copy__(self):
        r"""
        Make a copy of ``self``.

        """
        return self.__class__(upstairs_bsa=copy(self.upstairs()), polynomial_map=copy(self.polynomial_map()),
                              poly_ring=self.poly_ring(), ambient_dim=self.ambient_dim())

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
        r"""
        A list that maps the index of each generator in ``self.upstairs()``
        to a monomial.
        """
        return self._polynomial_map

    def plot_upstairs(self):
        r"""
        Plot the upstairs basic semialgebraic set and the variety::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: PY.<x,y> = QQ[]
            sage: upstairs = BasicSemialgebraicSet_eq_lt_le_sets(le=[-x-1, -y-1, x+y-3/2, x-1, y-1])
            sage: PX.<t> = QQ[]
            sage: F = [t, t^3 + 3/2 * t^2 - 5/2 * t + 5/4]
            sage: downstairs = upstairs.section(F)
            sage: downstairs.plot_upstairs()                  # not tested

        """
        from sage.plot.plot import parametric_plot
        gu = self.upstairs().plot()
        def range_arg(form, gen):
            a = self.linear_function_lower_bound(form)
            if a == -Infinity:
                a = -1
            b = self.linear_function_upper_bound(form)
            if b == +Infinity:
                b = 1
            return gen, a, b

        range_args = [ range_arg(form, gen)
                       for form, gen in zip(self.upstairs().linear_forms_space().basis(),
                                            self.poly_ring().gens()) ]
        return gu + parametric_plot(self.polynomial_map(), *range_args, color='yellow', thickness=2)

    def find_point(self):
        r"""
        Find a point in ``self``. FIXME: Lazy implementation. Raise NotImplementedError very often.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(3)
            sage: P.add_linear_constraint([1,0,0],0,operator.ge)
            sage: P.add_linear_constraint([0,1,0],0,operator.ge)
            sage: P.add_linear_constraint([0,0,1],0,operator.ge)
            sage: P.add_linear_constraint([1,0,0],-2,operator.le)
            sage: P.add_linear_constraint([0,1,0],-3,operator.le)
            sage: P.add_linear_constraint([0,0,1],-6,operator.le)
            sage: Q.<u,v> = QQ[]
            sage: F = [u, v, -3*u+4*v]
            sage: section = P.section(F)
            sage: pt = section.find_point()
            sage: pt
            (1, 3/2)
            sage: pt in section
            True
        """
        point_upstairs = self.upstairs().find_point()
        for f in self.polynomial_map():
            if f.degree() > 1:
                raise NotImplementedError
        A = matrix(QQ, [[f.monomial_coefficient(xi) for xi in self.poly_ring().gens()] for f in self.polynomial_map()])
        try:
            return A.solve_right(point_upstairs)
        except ValueError:
            raise NotImplementedError

## (4) Later... introduce a class that takes care of the monomial lifting etc.

class BasicSemialgebraicSet_veronese(BasicSemialgebraicSet_section):

    r"""
    A basic semialgebraic set that delegates to another semialgebraic set
    via a Veronese embedding (reformulation-linearization, RLT).

    Expand the polynomial in the standard monomial basis and replace each monomial by a new variable.
    Record monomials in polynomial_map and their corresponding variables in v_dict. Record the bounds on the monomials in bounds.
    The resulting linear expression in the extended space will be provided as inequality or equation in upstairs_bsa, which could, for example be represented by a PPL not-necessarily-closed polyhedron.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
        sage: upstairs_bsa_ppl = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(ambient_dim=0)
        sage: P.<f>=QQ[]
        sage: veronese = BasicSemialgebraicSet_veronese(upstairs_bsa_ppl, polynomial_map=[], v_dict=dict(), poly_ring=P)
        sage: lt_poly = [2*f - 2, f - 2, f^2 - f, -2*f, f - 1, -f - 1, -f, -2*f + 1]
        sage: for lhs in lt_poly: veronese.add_polynomial_constraint(lhs, operator.lt)
        sage: veronese.polynomial_map()
        [f, f^2]
        sage: sorted(veronese.eq_poly())
        []
        sage: sorted(veronese.le_poly())
        []
        sage: sorted(veronese.lt_poly())
        [-2*f + 1, f - 1, f^2 - f]
    """

    def __init__(self, upstairs_bsa=None, polynomial_map=None, v_dict=None,
                 upstairs_bsa_class=None, base_ring=None, poly_ring=None, ambient_dim=None, names=None):
        r"""
        EXAMPLES:

        Trivial initialization of the universe in dimension 3::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: upstairs_bsa_ppl = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(ambient_dim=0)
            sage: veronese = BasicSemialgebraicSet_veronese(upstairs_bsa_ppl, [], dict(), ambient_dim=3)
            sage: veronese.linear_forms_space()
            Vector space of dimension 3 over Rational Field

        Adding initial space dimensions::

            sage: names = ['u', 'v']
            sage: n = len(names)
            sage: P = PolynomialRing(QQ, names)
            sage: polynomial_map = list(P.gens())
            sage: v_dict = {P.gens()[i]:i for i in range(n)}
            sage: polyhedron = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(ambient_dim=n)
            sage: veronese = BasicSemialgebraicSet_veronese(polyhedron, polynomial_map, v_dict)
            sage: veronese.linear_forms_space()
            Vector space of dimension 2 over Rational Field

        Check error checking::

            sage: polynomial_map = [P.gen(0) - P.gen(1), 2 * P.gen(1)]
            sage: veronese = BasicSemialgebraicSet_veronese(polyhedron, polynomial_map, v_dict)
            Traceback (most recent call last):
            ...
            ValueError: all polynomials in polynomial_map must be monomials with coefficient 1

        Check polynomial_map and v_dict are set correctly in init::

            sage: P.<x,y,z>=QQ[]
            sage: veronese = BasicSemialgebraicSet_veronese(poly_ring=P, polynomial_map=list(P.gens()))
            sage: veronese.v_dict()
            {z: 2, y: 1, x: 0}
        """
        if poly_ring is None:
            if polynomial_map:
                cm = get_coercion_model()
                poly_ring = cm.common_parent(*polynomial_map)
                if not is_PolynomialRing(poly_ring) and not is_MPolynomialRing(poly_ring):
                    poly_ring = PolynomialRing(base_ring, [])
                polynomial_map = [ poly_ring(f) for f in polynomial_map ]
            elif base_ring is not None and ambient_dim is not None:
                poly_ring = PolynomialRing(base_ring, 'x', ambient_dim)
        if base_ring is None:
            if poly_ring is not None:
                base_ring = poly_ring.base_ring()
            elif upstairs_bsa is not None:
                base_ring = upstairs_bsa.base_ring()
        if upstairs_bsa is None:
            if upstairs_bsa_class is None:
                upstairs_bsa_class = 'ppl'
            upstairs_bsa_class = _bsa_class(upstairs_bsa_class)
            if polynomial_map is None:
                upstairs_base_ring = base_ring
                upstairs_ambient_dim = 0
                polynomial_map = []
                v_dict = {}
            else:
                upstairs_base_ring = poly_ring.base_ring()
                upstairs_ambient_dim = len(polynomial_map)
                if v_dict is None:
                    v_dict = {polynomial_map[i]:i for i in range(upstairs_ambient_dim)}
            upstairs_bsa = upstairs_bsa_class(base_ring=upstairs_base_ring, ambient_dim=upstairs_ambient_dim)
        if upstairs_bsa_class:
            if not isinstance(upstairs_bsa, upstairs_bsa_class):
                raise ValueError("upstairs_bsa is not an instance of upstairs_bsa_class")
        super(BasicSemialgebraicSet_veronese, self).__init__(
            upstairs_bsa, polynomial_map,
            poly_ring=poly_ring, base_ring=base_ring, ambient_dim=ambient_dim, names=names)
        if not all(len(f.monomials()) == 1 and f.lc() == 1 for f in self.polynomial_map()):
            raise ValueError("all polynomials in polynomial_map must be monomials with coefficient 1")
        self._v_dict = v_dict

    @classmethod
    def from_bsa(cls, bsa, poly_ring=None, **init_kwds):
        if poly_ring is None:
            poly_ring = bsa.poly_ring()
        return super(BasicSemialgebraicSet_veronese, cls).from_bsa(bsa, poly_ring=poly_ring, **init_kwds)

    def __copy__(self):
        r"""
        Make a copy of ``self``.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: bddbsa = BasicSemialgebraicSet_veronese(poly_ring=PolynomialRing(QQ, ['f','z']))
            sage: bddbsa.add_linear_constraint((0,-1), 0, operator.lt)
            sage: bddbsa.add_linear_constraint((-2,0), 1, operator.lt)
            sage: bddbsa.add_linear_constraint((1,6), -1, operator.lt)
            sage: bddbsa
            BasicSemialgebraicSet_veronese(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(Constraint_System {-6*x0-x1+1>0, x0>0, 2*x1-1>0}, names=[x0, x1]), polynomial_map=[z, f])
            sage: bddbsa_copy = copy(bddbsa); bddbsa_copy
            BasicSemialgebraicSet_veronese(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(Constraint_System {-6*x0-x1+1>0, x0>0, 2*x1-1>0}, names=[x0, x1]), polynomial_map=[z, f])
        """
        return self.__class__(upstairs_bsa=copy(self.upstairs()), polynomial_map=copy(self.polynomial_map()),
                              v_dict=copy(self.v_dict()), poly_ring=self.poly_ring(), ambient_dim=self.ambient_dim())

    def _repr_(self):
        return 'BasicSemialgebraicSet_veronese({}, polynomial_map={})'.format(self._upstairs_bsa, self._polynomial_map)

    def v_dict(self):
        r"""
        A dictionary that maps each monomial to the index of its corresponding generator
        in ``self.upstairs()``.
        """
        return self._v_dict

    def _to_upstairs_linear_constraint(self, lhs, allow_adding_upstairs_space_dimensions=True):
        lhs = self.poly_ring()(lhs)
        space_dim_to_add = 0
        upstairs_lhs_coeff = list(self.upstairs().linear_forms_space().zero())
        upstairs_lhs_cst = self.upstairs().base_ring().zero()
        for m in lhs.monomials():
            coeffm = lhs.monomial_coefficient(m)/1
            if m == 1:
                upstairs_lhs_cst = coeffm
            else:
                nv = self.v_dict().get(m, None)
                if nv is None:
                    if not allow_adding_upstairs_space_dimensions:
                        raise ValueError("this constraint requires adding upstairs space dimensions")
                    nv = len(self.polynomial_map())
                    self.v_dict()[m] = nv
                    self.polynomial_map().append(m)
                    space_dim_to_add += 1
                    upstairs_lhs_coeff.append(coeffm)
                else:
                    upstairs_lhs_coeff[nv] = coeffm
        if space_dim_to_add:
            self.upstairs().add_space_dimensions_and_embed(space_dim_to_add)
        return upstairs_lhs_coeff, upstairs_lhs_cst

    def add_polynomial_constraint(self, lhs, op):
        r"""
        ``lhs`` should be a polynomial.
        Add the constraint ``lhs``(x) ``op`` 0,
        where ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: upstairs_bsa_ppl = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(ambient_dim=0)
            sage: Q.<x0,x1,x2> = QQ[]
            sage: veronese = BasicSemialgebraicSet_veronese(upstairs_bsa_ppl, [], dict(), poly_ring=Q)
            sage: lhs = 27/113 * x0^2 + x1*x2 + 1/2
            sage: veronese.add_polynomial_constraint(lhs, operator.lt)
            sage: list(veronese.lt_poly())
            [54*x0^2 + 226*x1*x2 + 113]
            sage: veronese.polynomial_map()
            [x0^2, x1*x2]
            sage: veronese.v_dict()
            {x1*x2: 1, x0^2: 0}
            sage: lhs2 = x0 + 1/3*x1*x2
            sage: veronese.add_polynomial_constraint(lhs2, operator.lt)
            sage: sorted(veronese.lt_poly())
            [x1*x2 + 3*x0, 54*x0^2 + 226*x1*x2 + 113]
            sage: veronese.polynomial_map()
            [x0^2, x1*x2, x0]
            sage: veronese.v_dict()
            {x0: 2, x1*x2: 1, x0^2: 0}

            sage: P.<x,y,z> = QQ[]
            sage: bsa = BasicSemialgebraicSet_veronese(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(ambient_dim=0), [], dict(), poly_ring=P)
            sage: t = 27/113 * x^2 + y*z + 1/2
            sage: bsa.add_polynomial_constraint(t, operator.lt)
            sage: bsa.polynomial_map()
            [x^2, y*z]
            sage: bsa.v_dict()
            {y*z: 1, x^2: 0}
            sage: tt = x + 1/3 * y*z
            sage: bsa.add_polynomial_constraint(tt, operator.le)
            sage: bsa.polynomial_map()
            [x^2, y*z, x]
            sage: bsa.v_dict()
            {x: 2, y*z: 1, x^2: 0}
            sage: sorted(bsa.lt_poly()), sorted(bsa.le_poly())
            ([54*x^2 + 226*y*z + 113], [y*z + 3*x])
        """
        upstairs_lhs_coeff, upstairs_lhs_cst = self._to_upstairs_linear_constraint(lhs, True)
        self.upstairs().add_linear_constraint(upstairs_lhs_coeff, upstairs_lhs_cst, op)

    def polynomial_function_upper_bound(self, polynomial, allow_adding_upstairs_space_dimensions=False):
        r"""
        Find an upper bound for ``polynomial`` on ``self``.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: Q.<x0,x1,x2> = QQ[]
            sage: veronese = BasicSemialgebraicSet_veronese(poly_ring=Q)
            sage: lhs = 27/113 * x0^2 + x1*x2 + 1/2
            sage: veronese.add_polynomial_constraint(lhs, operator.lt)
            sage: veronese.polynomial_function_upper_bound(x1*x2)
            +Infinity
            sage: veronese.tighten_upstairs_by_mccormick(max_iter=1)
            sage: veronese.polynomial_function_upper_bound(x1*x2)
            -1/2
            sage: veronese._bounds
            [(0, +Infinity), (-Infinity, -1/2)]
            sage: veronese._bounds
            [(0, +Infinity), (-Infinity, -1/2)]
            sage: veronese.polynomial_function_upper_bound(x0*x0*x1*x2)
            +Infinity
            sage: veronese.polynomial_function_upper_bound(x0*x0*x1*x2, allow_adding_upstairs_space_dimensions=True)
            0
            sage: veronese._bounds   #allow_adding_upstairs_space_dimensions=True modified upstairs
            [(0, +Infinity), (-Infinity, -1/2), (-Infinity, 0)]
        """
        try:
            upstairs_coeff, upstairs_cst = self._to_upstairs_linear_constraint(polynomial, False)
        except ValueError:
            if allow_adding_upstairs_space_dimensions is True:
                upstairs_coeff, upstairs_cst = self._to_upstairs_linear_constraint(polynomial, True)
                self.tighten_upstairs_by_mccormick() # caution: modify upstairs
            else:
                return +Infinity
        return self.upstairs().linear_function_upper_bound(upstairs_coeff) + upstairs_cst

    def is_polynomial_constraint_valid(self, lhs, op, allow_adding_upstairs_space_dimensions=False):
        r"""
        Check if the constraint ``lhs``(x) ``op`` 0 is satisfied 
        for all points of ``self``, where ``lhs`` is a polynomial, and
        ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.

        Raise an error if the information provided does not suffice to decide
        the validity of the constraint.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: upstairs_bsa_ppl = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(ambient_dim=0)
            sage: Q.<x0,x1,x2> = QQ[]
            sage: veronese = BasicSemialgebraicSet_veronese(upstairs_bsa_ppl, [], dict(), poly_ring=Q)
            sage: lhs = 27/113 * x0^2 + x1*x2 + 1/2
            sage: veronese.add_polynomial_constraint(lhs, operator.lt)
            sage: veronese.is_polynomial_constraint_valid(lhs, operator.lt)
            True
            sage: lhs2 = x0 + 1/3*x1*x2
            sage: veronese.is_polynomial_constraint_valid(lhs2, operator.lt)
            Traceback (most recent call last):
            ...
            NotImplementedError...
        """
        try:
            upstairs_lhs_coeff, upstairs_lhs_cst = self._to_upstairs_linear_constraint(lhs, False)
        except ValueError:
            if allow_adding_upstairs_space_dimensions is True:
                upstairs_lhs_coeff, upstairs_lhs_cst = self._to_upstairs_linear_constraint(lhs, True)
                self.tighten_upstairs_by_mccormick() # caution: modify upstairs
            else:
                raise NotImplementedError
        if self.upstairs().is_linear_constraint_valid(upstairs_lhs_coeff, upstairs_lhs_cst, op):
            return True
        elif all(m.degree() < 2 for m in self._polynomial_map):
            return False
        else:
            raise NotImplementedError

    def _add_mccormick_bound(self, i, i1, i2, c1, c2, op):
        r"""
        If x[i] - c1*x[i1] - c2*x[i2] + c1*c2 ``op`` 0 is redundant in self.upstairs(), return False, otherwise return True and add this new constraint to self.upstairs().
        """
        if c1 in (-Infinity, +Infinity) or c2 in (-Infinity, +Infinity):
            # unbounded
            return False
        lhs = vector(self.upstairs().base_ring(), self.upstairs().ambient_dim())
        lhs[i] = 1; lhs[i1] -= c1; lhs[i2] -= c2; cst = c1 * c2;
        if self.upstairs().is_linear_constraint_valid(lhs, cst, op):
            return False
        self.upstairs().add_linear_constraint(lhs, cst, op)
        return True

    def _compute_bounds_of_the_ith_monomial(self, i):
        r"""
        Return the lower and upper bounds of the i-th variable in self.upstaris.
        """
        form = self.upstairs().linear_forms_space().basis()[i]
        lb = self.upstairs().linear_function_lower_bound(form)
        ub = self.upstairs().linear_function_upper_bound(form)
        return (lb, ub)

    def tighten_upstairs_by_mccormick(self, max_iter=1):
        r"""
        Recursively add McCormick inequalites on the monomials and do upward and downward bounds propagation until max_iter is attained or no more changes of bounds occur.
        Calling with max_iter=0 computes the bounds.
        Calling with max_iter=1 performs one round of upward bounds propagation and one round of downward bounds propagation.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: upstairs_bsa_mip = BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram(QQ, 0, solver='ppl')
            sage: P.<x,y> = QQ[]
            sage: veronese = BasicSemialgebraicSet_veronese(upstairs_bsa_mip, [], dict(), poly_ring=P)
            sage: veronese.add_polynomial_constraint(x, operator.ge)
            sage: veronese.add_polynomial_constraint(x-1, operator.le)
            sage: veronese.add_polynomial_constraint(y-1, operator.ge)
            sage: veronese.add_polynomial_constraint(y-2, operator.le)
            sage: veronese.add_polynomial_constraint(x*y+1000, operator.ge)
            sage: veronese.tighten_upstairs_by_mccormick(max_iter=0); veronese._bounds
            [(0, 1), (1, 2), (-1000, +Infinity)]
            sage: veronese.tighten_upstairs_by_mccormick(max_iter=1); veronese._bounds
            [(0, 1), (1, 2), (0, 2)]
            sage: veronese.add_polynomial_constraint(2*x*y^2-1, operator.ge)
            sage: veronese.polynomial_map()
            [x, y, x*y, x*y^2]
            sage: veronese.tighten_upstairs_by_mccormick(max_iter=0); veronese._bounds
            [(0, 1), (1, 2), (0, 2), (1/2, +Infinity)]
            sage: veronese.tighten_upstairs_by_mccormick(max_iter=1); veronese._bounds
            [(1/8, 1), (1, 2), (1/4, 2), (1/2, 4)]
            sage: veronese.add_polynomial_constraint(y^2-1, operator.le)
            sage: veronese.polynomial_map()
            [x, y, x*y, x*y^2, y^2]
            sage: veronese.tighten_upstairs_by_mccormick(max_iter=0); veronese._bounds
            [(1/8, 1), (1, 2), (1/4, 2), (1/2, 4), (0, 1)]
            sage: veronese.tighten_upstairs_by_mccormick(max_iter=1); veronese._bounds
            [(1/2, 1), (1, 1), (1/2, 1), (1/2, 1), (1, 1)]
        """
        for i in range(len(self._polynomial_map)):
            if self._polynomial_map[i].is_square():
                form = self.upstairs().linear_forms_space().basis()[i]
                self.upstairs().add_linear_constraint(form, QQ(0), operator.ge)
        self._bounds = [self._compute_bounds_of_the_ith_monomial(i) for i in range(len(self._polynomial_map))]
        bounds_propagation_iter = 0
        tightened = True # = bool(len(self.var_value) > 1)
        while bounds_propagation_iter < max_iter and tightened:
            tightened = False
            # upward bounds propagation
            for m in self._polynomial_map:
                if m.degree() < 2:
                    continue
                i = self._v_dict[m]
                for v1 in self._polynomial_map:
                    (v2, rem) = m.quo_rem(v1)
                    if rem != 0 or v1 > v2 or (not v2 in self._v_dict):
                        # don't want the recursive McCormicks to create a new monomial.
                        # v1, v2 symmetric
                        continue
                    i1 = self._v_dict[v1]
                    i2 = self._v_dict[v2]
                    lb1, ub1 = self._bounds[i1]
                    lb2, ub2 = self._bounds[i2]
                    if self._add_mccormick_bound(i, i1, i2, lb2, lb1, operator.ge):
                        tightened = True
                    if self._add_mccormick_bound(i, i1, i2, ub2, ub1, operator.ge):
                        tightened = True
                    if self._add_mccormick_bound(i, i1, i2, lb2, ub1, operator.le):
                        tightened = True
                    if self._add_mccormick_bound(i, i1, i2, ub2, lb1, operator.le):
                        tightened = True
                if tightened:
                    self._bounds[i] = self._compute_bounds_of_the_ith_monomial(i)
            if tightened:
                tightened = False
                # downward bounds propagation
                for i in range(self.upstairs().ambient_dim()):
                    (lb, ub) = self._bounds[i]
                    self._bounds[i] = self._compute_bounds_of_the_ith_monomial(i)
                    if (self._bounds[i][0] > lb + 0.001) or (self._bounds[i][1] < ub - 0.001):
                        tightened = True
            #if max_iter != 0 and bounds_propagation_iter >= max_iter:
            #    logging.warning("max number %s of bounds propagation iterations has attained." % max_iter)
            bounds_propagation_iter += 1
