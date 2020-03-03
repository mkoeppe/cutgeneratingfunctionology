"""
    sage: import warnings
    sage: warnings.filterwarnings('ignore', 'Matplotlib is building the font cache using fc-list. This may take a moment.')
"""

try:
    from ppl import Variable, Constraint, Linear_Expression, Constraint_System, NNC_Polyhedron, Poly_Con_Relation, Poly_Gen_Relation, Generator, MIP_Problem, point as ppl_point
except ImportError:
    # old Sage
    from sage.libs.ppl import Variable, Constraint, Linear_Expression, Constraint_System, NNC_Polyhedron, Poly_Con_Relation, Poly_Gen_Relation, Generator, MIP_Problem, point as ppl_point

from six.moves import zip
from six.moves import range

from sage.structure.sage_object import SageObject
import time

from cutgeneratingfunctionology.spam.basic_semialgebraic import *
from cutgeneratingfunctionology.spam.basic_semialgebraic_local import BasicSemialgebraicSet_local
from cutgeneratingfunctionology.spam.semialgebraic_mathematica import BasicSemialgebraicSet_mathematica, from_mathematica
from cutgeneratingfunctionology.spam.basic_semialgebraic_groebner_basis import BasicSemialgebraicSet_groebner_basis
from cutgeneratingfunctionology.spam.polyhedral_complex import PolyhedralComplex
from .parametric_family import Classcall, ParametricFamily_base, ParametricFamily

debug_new_factors = False
debug_cell_exceptions = False

def bigcellify_igp():
    """
    Load igp as follows to override primitives with big_cells primitives::

        sage: import cutgeneratingfunctionology.igp as igp; reload(igp); igp.bigcellify_igp(); from cutgeneratingfunctionology.igp import *   # not tested

    """
    import cutgeneratingfunctionology.igp as igp
    import cutgeneratingfunctionology.spam.big_cells as big_cells
    big_cells.bigcellify_module(igp)

###############################
# Parametric Real Number Field
###############################

from cutgeneratingfunctionology.spam.parametric_real_field_element import ParametricRealFieldElement, is_parametric_element

from sage.rings.ring import Field
import sage.rings.number_field.number_field_base as number_field_base
from sage.structure.coerce_maps import CallableConvertMap


class ParametricRealFieldFrozenError(ValueError):
    pass

class ParametricRealFieldInconsistencyError(ValueError):
    pass

class ParametricRealFieldRefinementError(ValueError):
    pass

from contextlib import contextmanager

class FactorUndetermined(Exception):
    pass

allow_refinement_default = True
big_cells_default = 'if_not_allow_refinement'
mutable_values_default = False

class ParametricRealField(Field):
    r"""
    A Metaprogramming trick for parameter space analysis.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: K.<f> = ParametricRealField([4/5])
        sage: h = gmic(f, field=K)
        sage: I_list = J_list = h.intervals()
        sage: K_list = I_list + [ (a+1, b+1) for (a, b) in h.intervals() ]
        sage: _ = [ verts(I, J, KK) for I in I_list for J in J_list for KK in K_list ]
        sage: K.get_eq()
        set()
        sage: R = f._sym.parent().ring()
        sage: sorted([ p for p in K.get_lt() if p in R ])   # filter out the rational function 1/(f^2 - f), which is normalized differently starting Sage 8.4b2
        [-2*f, -2*f + 1, -f - 1, -f, f - 2, f - 1, 2*f - 2]
        sage: K.get_eq_factor()
        set()
        sage: K.get_lt_factor()
        {-f, -f + 1/2, f - 1}

        sage: K.<f, lam> = ParametricRealField([4/5, 1/6])
        sage: h = gj_2_slope(f, lam, field=K)
        sage: K.get_lt()
        {(-1/2)/(-1/2*f^2*lam - 1/2*f^2 + f*lam + 1/2*f - 1/2*lam),
         (-f*lam - f + lam)/(-f + 1),
         -lam,
         lam - 1,
         -f,
         f - 1,
         -1/2*f*lam - 1/2*f + 1/2*lam,
         f*lam - lam}
        sage: K.get_lt_factor()
        {-lam, lam - 1, -f, f - 1, -f*lam - f + lam}

        sage: K.<f,alpha> = ParametricRealField([4/5, 3/10])
        sage: h=dg_2_step_mir(f, alpha, field=K, conditioncheck=False)
        sage: extremality_test(h)
        True

        sage: K.<f,a1,a2,a3> = ParametricRealField([4/5, 1, 3/10, 2/25])
        sage: h = kf_n_step_mir(f, (a1, a2, a3), conditioncheck=False)
        sage: extremality_test(h)
        True

        sage: K.<f> = ParametricRealField([1/5])
        sage: h = drlm_3_slope_limit(f, conditioncheck=False)
        sage: extremality_test(h)
        True
        sage: K.get_lt_factor()
        {-f, f - 1, f - 1/2, f - 1/3}

    Elements coerce to RDF, RR, float to enable plotting of functions::

        sage: K.<f> = ParametricRealField([4/5])
        sage: RDF(f)
        0.8
        sage: RR(f)
        0.800000000000000
        sage: float(f)
        0.8
        sage: 0.2 - f
        -0.600000000000000

    Plotting will show symbolic labels on the axes::

        sage: plot_with_colored_slopes(gmic(f))
        Graphics object...

    But they do not coerce or convert into any exact fields::

        sage: QQ(f)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert f~ to a rational
        sage: AA(f)
        Traceback (most recent call last):
        ...
        TypeError: Illegal initializer for algebraic number

    Test that values can also be initialized by vectors (internally we make it tuples)::

        sage: f = vector(QQ, [1, 2])
        sage: K.<f0, f1> = ParametricRealField(f)
        sage: f = vector(f0.parent(), [f0, f1]); f
        ((f0)~, (f1)~)
        sage: f.parent()
        Vector space of dimension 2 over ParametricRealField(names = ['f0', 'f1'], values = [1, 2])
        sage: f[0]*f[1] <= 4
        True

    Test-point free mode (limited functionality and MUCH slower because of many more polynoial
    evaluations via libsingular)::

        sage: K.<a, b> = ParametricRealField(None, mutable_values=True)
        sage: a <= 2
        Traceback (most recent call last):
        ...
        FactorUndetermined: a cannot be evaluated because the test point is not complete
        sage: K.assume_comparison(a.sym(), operator.le, 3)

    Partial test point mode::

        sage: K.<a, b> = ParametricRealField([None, 1], mutable_values=True)
        sage: a <= 2
        Traceback (most recent call last):
        ...
        FactorUndetermined: a cannot be evaluated because the test point is not complete
        sage: b <= 11
        True

    """
    Element = ParametricRealFieldElement

    def __init__(self, values=None, names=None, allow_coercion_to_float=True,
                 mutable_values=None, allow_refinement=None, big_cells=None,
                 base_ring=None, sym_ring=None, bsa=None):
        Field.__init__(self, self)

        if mutable_values is None:
            mutable_values = mutable_values_default
        if allow_refinement is None:
            allow_refinement = allow_refinement_default
        if big_cells is None:
            big_cells = big_cells_default

        if big_cells == 'if_not_allow_refinement':
            big_cells = not allow_refinement
        if mutable_values == 'if_big_cells':
            mutable_values = big_cells

        assert isinstance(allow_refinement, bool)
        assert isinstance(mutable_values, bool)
        assert isinstance(big_cells, bool)

        if not allow_refinement and not big_cells:
            raise ValueError("if allow_refinement=False, must have big_cells=True")

        self._allow_refinement = allow_refinement
        self._mutable_values = mutable_values
        self._big_cells = big_cells

        self._zero_element = ParametricRealFieldElement(self, 0)
        self._one_element =  ParametricRealFieldElement(self, 1)
        ## REFACTOR: Maybe replace this by an instance of BasicSemialgebraicSet_eq_lt_le_sets - but careful - this class right now assumes polynomials
        self._eq = set([])
        self._lt = set([])
        self._le = set([])

        if sym_ring is None:
            if bsa is not None:
                sym_ring = bsa.poly_ring()
        if base_ring is None:
            if sym_ring is not None:
                base_ring = sym_ring.base_ring()
            else:
                base_ring = QQ
        if names is None:
            if sym_ring is not None:
                names = sym_ring.gens()
            else:
                raise ValueError("must provide one of names, sym_ring, or bsa")
        if sym_ring is None:
            #sym_ring = PolynomialRing(base_ring, names, implementation='generic')
            sym_ring = PolynomialRing(base_ring, names)
        self._factor_bsa = BasicSemialgebraicSet_eq_lt_le_sets(poly_ring=sym_ring)
        self._sym_field = sym_ring.fraction_field()
        if values is None:
            values = [ None for n in names ]
        else:
            assert len(values) == len(names)
        vnames = self._sym_field.gens()
        self._gens = [ ParametricRealFieldElement(self, value, name) for (value, name) in zip(values, vnames) ]
        self._names = tuple(names)
        if mutable_values:
            self._values = list(values)
        else:
            self._values = tuple(values)

        self._frozen = False
        self._record = True

        if bsa is None:
            # do the computation of the polyhedron incrementally,
            # rather than first building a huge list and then in a second step processing it.
            # the upstairs polyhedron defined by all constraints in self._eq/lt_factor
            polyhedron = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(0)
            # monomial_list records the monomials that appear in self._eq/lt_factor.
            # v_dict is a dictionary that maps each monomial to the index of its corresponding Variable in polyhedron
            bsa = BasicSemialgebraicSet_veronese(polyhedron, polynomial_map=[], poly_ring=sym_ring, v_dict={})
        self._bsa = bsa

        self.allow_coercion_to_float = allow_coercion_to_float
        if allow_coercion_to_float:
            RDF.register_coercion(sage.structure.coerce_maps.CallableConvertMap(self, RDF, lambda x: RDF(x.val()), parent_as_first_arg=False))
            RR.register_coercion(sage.structure.coerce_maps.CallableConvertMap(self, RR, lambda x: RR(x.val()), parent_as_first_arg=False))
        logging.info("Initialized {}".format(self))

    def __copy__(self):
        logging.warning("copy(%s) is invoked" % self)
        Kcopy = self.__class__(self._values, self._names, allow_coercion_to_float=self.allow_coercion_to_float, mutable_values=self._mutable_values, allow_refinement=self._allow_refinement, big_cells=self._big_cells, base_ring=self._sym_field.base_ring())
        Kcopy._eq.update(self._eq)
        Kcopy._lt.update(self._lt)
        Kcopy._le.update(self._le)
        Kcopy._factor_bsa = copy(self._factor_bsa)
        Kcopy._bsa = copy(self._bsa)
        Kcopy._frozen = self._frozen
        Kcopy._record = self._record
        return Kcopy

    def ppl_polyhedron(self):
        return self._bsa.upstairs()._polyhedron

    def monomial_list(self):
        return self._bsa.polynomial_map()

    def v_dict(self):
        return self._bsa.v_dict()

    def _first_ngens(self, n):
        for i in range(n):
            yield self._gens[i]
    def ngens(self):
        return len(self._gens)
    def _an_element_impl(self):
        return ParametricRealFieldElement(self, 1)
    def _coerce_map_from_(self, S):
        """
        TESTS:

        Test that elements of different ``ParametricRealField``s have no coercion.

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: K.<a> = ParametricRealField([0])
            sage: L.<b> = ParametricRealField([1])
            sage: a + b
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s)...

        Test that real number field elements can be upgraded to ``ParametricRealFieldElement``s.
        Note that this requires setting up the ParametricRealField with a specific base ring, 
        because there is no common parent of QQ(x) and a RealNumberField``::

            sage: sqrt2, = nice_field_values([sqrt(2)])
            sage: K.<f> = ParametricRealField([0], base_ring=sqrt2.parent())
            sage: f + sqrt2
            (f + 1.414213562373095?)~

        This currently does not work for Sage's built-in embedded number field elements...
        """
        if isinstance(S, ParametricRealField) and self is not S:
            return None
        if  S is sage.interfaces.mathematica.MathematicaElement or isinstance(S, RealNumberField_absolute) or isinstance(S, RealNumberField_quadratic) or AA.has_coerce_map_from(S):
            # Does the test with MathematicaElement actually work?
            # We test whether S coerces into AA. This rules out inexact fields such as RDF.
            return True
    def __repr__(self):
        return 'ParametricRealField(names = {}, values = {})'.format(list(self._names), list(self._values))

    def freeze(self):
        self._frozen = True

    def unfreeze(self):
        self._frozen = False

    @contextmanager
    def frozen(self):
        was_frozen = self._frozen
        self._frozen = True
        try:
            yield True
        finally:
            self._frozen = was_frozen

    @contextmanager
    def unfrozen(self):
        was_frozen = self._frozen
        self._frozen = False
        try:
            yield True
        finally:
            self._frozen = was_frozen

    ## def wallcross_to_lt(self, comparison):
    ##     record_to_lt(self, comparison)
    ##     ..... compute new testpoint, error if fail.

    def is_point_consistent(self, new_values):
        """
        Check whether ``new_values`` satisfies all recorded assumptions.
        """
        if new_values is None:
            return True
        if any(x is None for x in new_values):
            if all(x is None for x in new_values):
                return True
            logging.warning("Consistency checking not implemented if some test point coordinates are None")
            return True
        ### Check that values satisfy all constraints.
        return new_values in self._bsa

    def change_test_point(self, new_values):
        if not self._mutable_values:
            raise ValueError("ParametricRealField is set up with mutable_values=False")
        new_values = list(new_values)
        if not self.is_point_consistent(new_values):
            raise ParametricRealFieldInconsistencyError("New test point {} does not satisfy the recorded constraints".format(new_values))
        self._values = new_values

    def change_values(self, **values):
        """
        Convenience interface for ``change_test_point``.
        """
        if not self._mutable_values:
            raise ValueError("ParametricRealField is set up with mutable_values=False")
        new_values = copy(self._values)
        for key, value in values.items():
            new_values[self._names.index(key)] = value
        self.change_test_point(new_values)

    def remove_test_point(self):
        """
        Switch ``self`` to test-point free mode.

        This requires a ParametricRealField set up with ``mutable_values=True``.

        Not many things are implemented in this mode::

        EXAMPLE::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: K.<a,b> = ParametricRealField([4, 1], big_cells=True, mutable_values=True, allow_refinement=False)
            sage: with K.changed_values():
            ....:     K.remove_test_point()
            ....:     with K.temporary_assumptions():
            ....:         K.assume_comparison(a.sym(), operator.le, 3)
            ....:         a <= 4
            Traceback (most recent call last):
            ...
            FactorUndetermined: a cannot be evaluated because the test point is not complete...
        """
        self._values = [ None for n in self._names ]

    @contextmanager
    def removed_test_point(self):
        """
        Context manager for temporarily switching to test-point free mode.

        This requires a ParametricRealField set up with ``mutable_values=True``.

        EXAMPLE::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: K.<a> = ParametricRealField([4], big_cells=True, mutable_values=True, allow_refinement=False)
            sage: with K.removed_test_point():
            ....:     with K.temporary_assumptions():
            ....:         K.assume_comparison(a.sym(), operator.le, 3)
            ....:         K.find_test_point()
            ....:         with K.frozen():
            ....:             a <= 4
            True
        """
        with self.changed_values():
            self.remove_test_point()
            yield True

    def find_test_point(self):
        """
        Sets a new test point that is consistent with the recorded constraints.

        This can fail with ``NotImplementedError'' or ``ParametricRealFieldInconsistencyError''.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: K.<x,y,z> = ParametricRealField([1, 2, 3], mutable_values=True, big_cells=True, allow_refinement=False)
            sage: assert 0 <= x <= y <= z <= 6
            sage: K.remove_test_point()
            sage: K.find_test_point()
            sage: x.val(), y.val(), z.val()
            (3/2, 3, 9/2)

            sage: K.<x,y,z> = ParametricRealField(mutable_values=True, big_cells=True, allow_refinement=False)
            sage: K.assume_comparison(0, operator.le, x.sym())
            sage: K.assume_comparison(x.sym(), operator.le, y.sym())
            sage: K.assume_comparison(y.sym(), operator.le, z.sym())
            sage: K.assume_comparison(z.sym(), operator.le, 4)
            sage: K.find_test_point()
            sage: x.val(), y.val(), z.val()
            (1, 2, 3)
        """
        p = self._bsa.find_point()
        self.change_test_point(p)

    @contextmanager
    def changed_values(self, **values):
        """
        Context manager for temporarily switching to another consistent test point.

        This requires setting up the field with ``mutable_values=True``.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: K.<f> = ParametricRealField([4/5], mutable_values=True)
            sage: 0 <= f <= 1
            True
            sage: with K.changed_values(f=1/5):
            ....:     assert f <= 1/2
            Traceback (most recent call last):
            ...
            ParametricRealFieldInconsistencyError: Old test point... does not satisfy the recorded constraints
        """
        save_values = self._values
        try:
            self.change_values(**values)
            yield True
        finally:
            self._values = save_values
            if not self.is_point_consistent(save_values):
                raise ParametricRealFieldInconsistencyError("Old test point {} does not satisfy the recorded constraints".format(save_values))

    @contextmanager
    def off_the_record(self):
        was_recording = self._record
        self._record = False
        try:
            yield True
        finally:
            self._record = was_recording

    @contextmanager
    def temporary_assumptions(self, case_id=None):
        if case_id is not None:
            logging.info("Entering case {}.".format(case_id))
        save_eq = self._eq
        self._eq = copy(save_eq)
        save_lt = self._lt
        self._lt = copy(save_lt)
        save_le = self._le
        self._le = copy(save_le)
        save_factor_bsa = self._factor_bsa
        self._factor_bsa = copy(save_factor_bsa)
        save_bsa = self._bsa
        self._bsa = copy(self._bsa)
        try:
            yield True
        finally:
            self._eq = save_eq
            self._lt = save_lt
            self._le = save_le
            self._factor_bsa = save_factor_bsa
            self._bsa = save_bsa
            if case_id is not None:
                logging.info("Finished case {}.".format(case_id))

    def get_eq(self):
        return self._eq
    def get_lt(self):
        return self._lt
    def get_le(self):
        return self._le

    def get_eq_factor(self):
        return self._factor_bsa.eq_poly()
    def get_lt_factor(self):
        return self._factor_bsa.lt_poly()
    def get_le_factor(self):
        return self._factor_bsa.le_poly()

    def _eval_factor(self, fac):
        """
        Evaluate ``fac`` on the test point.

        If there is no test point or the test point has some ``None`` coordinates
        that are needed for the evaluation, raise ``FactorUndetermined``.
        """
        base_ring = self._sym_field.base_ring()
        if fac in base_ring:
            return base_ring(fac)
        if self._values is not None and not all(x is None for x in self._values):
            try:
                return fac(self._values)
            except TypeError:             # 'None' components
                pass
        raise FactorUndetermined("{} cannot be evaluated because the test point is not complete".format(fac))

    def _factor_sign(self, fac):
        """
        Determine the sign of ``fac`` evaluated on the test point.

        If there is no test point or the test point has some ``None`` coordinates
        that are needed for the evaluation, raise ``FactorUndetermined``.
        """
        return sign(self._eval_factor(fac))

    def assume_comparison(self, lhs, op, rhs=0):
        r"""
        Record the assumption ``lhs op rhs``.

        If this assumption is not satisfied by the current test point,
        a ``ParametricRealFieldInconsistencyError`` is raised.

        ``lhs`` and ``rhs`` should be elements of ``self._sym_field``
        or be parametric elements (elements of ``self``).

        TESTS for consistency checks::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)
            sage: K.<f> = ParametricRealField([4/5])
            sage: K.assume_comparison(f.sym(), operator.gt, 1)
            Traceback (most recent call last):
            ...
            ParametricRealFieldInconsistencyError: New constraint...
            sage: K.assume_comparison(K(0).sym(), operator.gt, 1)
            Traceback (most recent call last):
            ...
            ParametricRealFieldInconsistencyError: New constant constraint...

        TESTS for consistency checks when testpoint values is None::

            sage: K.<a,b> = ParametricRealField([2, 1], big_cells=True, mutable_values=True, allow_refinement=False)
            sage: assert b>0
            sage: K.remove_test_point()
            sage: K.assume_comparison(b.sym(), operator.lt, 0)
            Traceback (most recent call last):
            ...
            ParametricRealFieldInconsistencyError...

        User code should not rely on whether assumptions regarding
        the nonvanishing of denominators are recorded::

            sage: K.<f> = ParametricRealField([4/5])
            sage: assert f < 1
            sage: K.freeze()
            sage: K.assume_comparison(1/f.sym(), operator.gt, 1)   # not tested

        Therefore, also the implicit cancellation in the field of rational
        functions does not change semantics when we transform lhs op rhs
        to lhs - rhs op 0::

            sage: K.<f> = ParametricRealField([4/5])
            sage: K.freeze()
            sage: K.assume_comparison(1/f.sym(), operator.eq, 1/f.sym())

        User code should call assume_comparison only for operands that do
        not have vanishing denominators.

            sage: K.<f> = ParametricRealField([4/5])
            sage: assert 0 < f < 1
            sage: K.freeze()
            sage: K.assume_comparison(1/f.sym(), operator.gt, 1)

        TEST that parametric elements are allowed::

            sage: K.<f> = ParametricRealField([4/5])
            sage: K.assume_comparison(f, operator.le, 1)
            sage: K.get_le_factor()
            {f - 1}

        TESTS for allow_refinement=True:

        Strict inequalities - even-multiplicity factors do matter::

            sage: K.<f> = ParametricRealField([4/5])
            sage: K.freeze()
            sage: f^2 > 0
            Traceback (most recent call last):
            ...
            ParametricRealFieldFrozenError...

        Weak inequalities::

            sage: K.<f> = ParametricRealField([4/5], big_cells=True)
            sage: f >= 4/5
            True
            sage: K.freeze()
            sage: f == 4/5
            Traceback (most recent call last):
            ...
            ParametricRealFieldFrozenError...

            sage: K.<f> = ParametricRealField([4/5], big_cells=True)
            sage: f <= 4/5
            True
            sage: K.freeze()
            sage: f == 4/5
            Traceback (most recent call last):
            ...
            ParametricRealFieldFrozenError...

        Inequations::

            sage: K.<f> = ParametricRealField([4/5], big_cells=True)
            sage: assert f != 0
            sage: K.get_lt_factor()
            {-f}
            sage: K.<f> = ParametricRealField([4/5], big_cells=False)
            sage: assert f != 0
            sage: K.get_lt_factor()
            {-f}

        TESTS for allow_refinement=False:

        Strict inequalities::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)
            sage: K.<f> = ParametricRealField([4/5], allow_refinement=False)
            sage: f * (1 - f) > 0
            Traceback (most recent call last):
            ...
            ParametricRealFieldRefinementError...
            sage: K.<f> = ParametricRealField([4/5], allow_refinement=False)
            sage: f > 0
            True
            sage: 1 - f > 0
            True

        Strict inequalities - even-multiplicity factors do matter::

            sage: K.<f> = ParametricRealField([4/5], allow_refinement=False)
            sage: f^2 > 0
            Traceback (most recent call last):
            ...
            ParametricRealFieldRefinementError...


        Equalities::

            sage: K.<f> = ParametricRealField([4/5], allow_refinement=False)
            sage: f * (f - 4/5) == 0
            Traceback (most recent call last):
            ...
            ParametricRealFieldRefinementError...
            sage: f == 4/5
            True
            sage: f * (f - 4/5) == 0
            True

            sage: K.<f> = ParametricRealField([4], allow_refinement=False)
            sage: f >= 3
            True
            sage: f * (f - 1) * (f - 2) * (f - 4) == 0
            True

        Weak inequalities::

            sage: K.<f> = ParametricRealField([4/5], allow_refinement=False)
            sage: f * (f - 4/5) >= 0
            Traceback (most recent call last):
            ...
            ParametricRealFieldRefinementError...
            sage: f == 4/5
            True
            sage: f * (f - 4/5) >= 0
            True

            sage: K.<x,y> = ParametricRealField([2, 3], allow_refinement=False)
            sage: x >= 2
            True
            sage: x^2 * (x - 1) * y >= 0
            True

        Weak inequalities - Even-multiplicity factors do matter::

            sage: K.<o> = ParametricRealField([0], allow_refinement=False)
            sage: o^2 <= 0
            True
            sage: K.freeze()
            sage: o == 0
            True

            sage: K.<x> = ParametricRealField([0], allow_refinement=False)
            sage: x^2 >= 0
            True
            sage: K.freeze()
            sage: x == 0
            Traceback (most recent call last):
            ...
            ParametricRealFieldFrozenError...

            sage: K.<f> = ParametricRealField([2], allow_refinement=False)
            sage: f^2 * (1 - f) <= 0    # disjunction 0 union [1, oo).
            Traceback (most recent call last):
            ...
            ParametricRealFieldRefinementError...

        Not-equal::

            sage: K.<x> = ParametricRealField([1], allow_refinement=False)
            sage: x != 0
            Traceback (most recent call last):
            ...
            ParametricRealFieldRefinementError...
            sage: x >= 0
            True
            sage: x != 0
            True

            sage: K.<x> = ParametricRealField([-1], allow_refinement=False)
            sage: x <= 0
            True
            sage: x != 0
            True

        Bug example::
            sage: K.<f, bkpt> = ParametricRealField([1236451991/5178221721, 3733953/12155000], big_cells=True, allow_refinement=True)
            sage: for lhs in [2*f - 1, bkpt - 1, (-bkpt)/(f + 1), f - 1, -bkpt, bkpt/(-f - 1), -f + 2*bkpt - 1, (f - bkpt + 1)/(-f - 1), f + bkpt - 1, -f, (-f + bkpt - 1)/(f + 1), f - bkpt]:
            ....:     assert(lhs <=0)
            sage: for lhs in [(-1)/(f^4 - 2*f^3*bkpt + f^2*bkpt^2 + 2*f^3 - 4*f^2*bkpt + 2*f*bkpt^2 + f^2 - 2*f*bkpt + bkpt^2), -f^2 - 2*f - 1, -f + 2*bkpt - 1, (-bkpt^2)/(f^6 - 2*f^5*bkpt + f^4*bkpt^2 + 2*f^5 - 4*f^4*bkpt + 2*f^3*bkpt^2 + f^4 - 2*f^3*bkpt + f^2*bkpt^2), -2*f + bkpt, -f^2 + 4*f*bkpt - 4*bkpt^2 - 2*f + 4*bkpt - 1, -4*f^2 + 8*f - 4, -f^2, -f, (-f^2 + 2*f*bkpt - bkpt^2 - 2*f + 2*bkpt - 1)/(f^2 + 2*f + 1), f - bkpt, -f^2 + 2*f*bkpt - bkpt^2, -f^2 + 2*f - 1, (-bkpt^2)/(f^2 + 2*f + 1), -bkpt^2 + 2*bkpt - 1]:
            ....:     assert(lhs < 0)
            sage: -f^2 - 2*f*bkpt - bkpt^2 + 2*f + 2*bkpt - 1 < 0
            True
            sage: sorted(K._bsa.eq_poly())  # bug was output [-1]
            []
            sage: sorted(K._bsa.lt_poly())  # bug was output []
            [-2*f + bkpt, f - bkpt, f + bkpt - 1]
            sage: sorted(K._bsa.le_poly())
            []
            sage: sorted(K._factor_bsa.lt_poly()) # bug was output [-f - bkpt + 1, -f + 2*bkpt - 1, -f - 1, -f, -2*f + bkpt, f - bkpt]
            [-2*f + bkpt, -f, -f - 1, -f + 2*bkpt - 1, f - bkpt, f + bkpt - 1]
            sage: sorted(K._factor_bsa.le_poly())
            [-bkpt, bkpt - 1, -f, -f + 2*bkpt - 1, f - bkpt, f + bkpt - 1, 2*f - 1]
            sage: sorted(K._factor_bsa.eq_poly())
            []
        """
        if not self._record:
            return
        if op == operator.gt:
            op = operator.lt
            lhs, rhs = rhs, lhs
        elif op == operator.ge:
            op = operator.le
            lhs, rhs = rhs, lhs
        elif op == operator.ne:
            op = operator.lt
            lhs = - (lhs - rhs) ** 2
            rhs = 0
        comparison = lhs - rhs
        # FIXME: The line above may cancel denominators.  If assume_comparison is called from _richcmp_, that's fine because
        # _div_ records denominators.
        ### Caching via self._eq, _lt, _le
        if is_parametric_element(comparison):
            comparison_sym = comparison.sym()
        else:
            comparison_sym = comparison
        if op == operator.eq:
            if (comparison_sym in self._eq) or (-comparison_sym in self._eq):
                return
        elif op == operator.lt:
            if comparison_sym in self._lt:
                return
        elif op == operator.le:
            if comparison_sym in self._le:
                return
        ### Test consistency with test point
        base_ring = self._sym_field.base_ring()
        if is_parametric_element(comparison):
            # Fast path: Use the precomputed val
            try:
                comparison_val = comparison.val()
            except FactorUndetermined:
                comparison_val = None
            comparison = comparison.sym()
        else:
            comparison = self._sym_field(comparison)
            try:
                comparison_val = self._eval_factor(comparison)
            except FactorUndetermined:
                comparison_val = None
        if comparison_val is not None:
            if not op(comparison_val, 0):
                if comparison in base_ring:
                    raise ParametricRealFieldInconsistencyError("New constant constraint {} {} {} is not satisfied".format(lhs, op, rhs))
                else:
                    raise ParametricRealFieldInconsistencyError("New constraint {} {}  {} is not satisfied by the test point".format(lhs, op, rhs))
        if comparison in base_ring:
            return
        if comparison.denominator() == 1 and comparison.numerator().degree() == 1:
            # Fast path for linear
            numerator = comparison.numerator()
            # We use two different normalizations for the univariate and the multivariate case,
            # just to match the behavior of factor() for the doctests.
            if self.ngens() == 1:
                unit = abs(comparison.numerator().lc())
            else:
                the_lcm = lcm([coeff.denominator() for coeff in numerator.coefficients()])
                numerator *= the_lcm
                unit = 1
            factors = Factorization([(numerator / unit, 1)], unit=unit)
        elif comparison.denominator() == 1:
            # Works around a sage bug in 9.1.beta1:
            numerator = comparison.numerator()
            factors = numerator.factor()
        else:
            factors = comparison.factor()
        if op in (operator.eq, operator.le):
            for (fac, d) in factors:
                if d > 0:
                    if self.is_factor_known(fac, operator.eq):
                    # Comparison is already known true, nothing to record.
                        return
        if len(factors) == 1 and factors[0][1] == 1 and comparison_val is not None:
            the_fac, d = factors[0]
            the_sign = sign(factors.unit() * comparison_val) 
            def factor_sign(fac):
                if fac == the_fac:
                    return the_sign
                else:
                    return self._factor_sign(fac)
        else:
            factor_sign = self._factor_sign
        unit_sign = base_ring(factors.unit()).sign()
        if op == operator.eq:
            eq_factors = []
            ne_factors = []
            unknown_factors = []
            # Record all factors for which testpoint gives zero.
            for (fac, d) in factors: #only one of -fac and fac is in the factorization.
                if d > 0:
                    if self.is_factor_known(fac, operator.lt) or self.is_factor_known(-fac, operator.lt):
                        pass
                    else:
                        try:
                            fac_sign = factor_sign(fac)
                            if fac_sign == 0:
                                eq_factors.append(fac) # cannot happen that -fac was in.
                            else:
                                ne_factors.append(fac)
                        except FactorUndetermined:
                            unknown_factors.append(fac)
            assert len(eq_factors + unknown_factors) > 0
            if self._allow_refinement:
                # Record all factors for which testpoint gives zero.
                for fac in eq_factors:
                    self.record_factor(fac, operator.eq)
            else:
                all_factors = eq_factors + ne_factors + unknown_factors
                if len(all_factors) > 1:
                    raise ParametricRealFieldRefinementError("{} == 0 has several new factors: {}".format(comparison, all_factors))
                self.record_factor(all_factors[0], op)
            logging.debug("New element in %s._eq: %s" % (repr(self), comparison))
            self._eq.add(comparison)
        elif op == operator.lt:
            lt_factors = []
            even_factors = []
            unknown_factors = []
            for (fac, d) in factors: #only one of -fac and fac is in the factorization.
                if self.is_factor_known(fac, operator.eq):
                    raise ParametricRealFieldInconsistencyError("New constraint {} {} {} is not satisfied".format(lhs, op, rhs))
                if d % 2 == 1:
                    if self.is_factor_known(fac, operator.lt):
                        unit_sign = -unit_sign
                    elif self.is_factor_known(-fac, operator.lt):
                        pass
                    else:
                        try:
                            fac_sign = factor_sign(fac)
                            if fac_sign == -1:
                                lt_factors.append(fac)
                                unit_sign = -unit_sign
                            elif fac_sign == +1:
                                lt_factors.append(-fac)
                            else:
                                # assert fac_sign == 0
                                raise ParametricRealFieldInconsistencyError("New constraint {} {} {} is not satisfied".format(lhs, op, rhs))
                        except FactorUndetermined:
                            unknown_factors.append(fac)
                else:
                    if not self.is_factor_known(fac, operator.lt) and not self.is_factor_known(-fac, operator.lt):
                        even_factors.append(fac)
            if (unit_sign == 1) and (not even_factors) and (not unknown_factors):
                raise ParametricRealFieldInconsistencyError("New constraint {} {} {} is not satisfied".format(lhs, op, rhs))
            if not self._allow_refinement:
                if len(lt_factors) + len(unknown_factors) > 1:
                    # or just record the product
                    raise ParametricRealFieldRefinementError("{} < 0 has several new factors: {}".format(comparison, lt_factors))
            for new_fac in lt_factors:
                self.record_factor(new_fac, op)
            if len(unknown_factors) > 1:
                raise NotImplementedError()
            if unknown_factors:
                if unit_sign > 0:
                    self.record_factor(unknown_factors[0], op)
                elif unit_sign < 0:
                    self.record_factor(-unknown_factors[0], op)
            for new_fac in even_factors:
                if not self._allow_refinement:
                    if self.is_factor_known(new_fac, operator.le):
                        self.record_factor(new_fac, operator.lt)
                    elif self.is_factor_known(-new_fac, operator.le):
                        self.record_factor(-new_fac, operator.lt)
                    else:
                        raise ParametricRealFieldRefinementError("{} < 0 has factor {} != 0".format(comparison, new_fac))
                else:
                    try:
                        fac_sign = factor_sign(new_fac)
                        if fac_sign == -1:
                            self.record_factor(new_fac, operator.lt)
                        elif fac_sign == 1:
                            self.record_factor(-new_fac, operator.lt)
                        else:
                            #assert fac_sign == 0
                            raise ParametricRealFieldInconsistencyError("New constraint {} {} {} is not satisfied".format(lhs, op, rhs))
                    except FactorUndetermined:
                        # break the cell ? raise NotImplementedError()?
                        self.record_factor(new_fac, operator.lt)
            logging.debug("New element in %s._lt: %s" % (repr(self), comparison))
            self._lt.add(comparison)
        elif op == operator.le:
            lt_factors = []
            eq_factors = []
            even_factors = []
            unknown_factors = []
            for (fac, d) in factors: #only one of -fac and fac is in the factorization.
                # self.is_factor_known(fac, operator.eq) or self.is_factor_known(-fac, operator.eq) were treated above already.
                if d % 2 == 1:
                    if self.is_factor_known(fac, operator.lt):
                        unit_sign = -unit_sign
                    elif self.is_factor_known(-fac, operator.lt):
                        pass
                    else:
                        try:
                            fac_sign = factor_sign(fac)
                            if fac_sign == 0:
                                eq_factors.append(fac) # cannot happen that -fac was in.
                            elif fac_sign == -1:
                                lt_factors.append(fac)
                                unit_sign = -unit_sign
                            else:
                                assert fac_sign == +1
                                lt_factors.append(-fac)
                        except FactorUndetermined:
                            unknown_factors.append(fac)
                else:
                    if not self.is_factor_known(fac, operator.lt) and not self.is_factor_known(-fac, operator.lt):
                        even_factors.append(fac)
            if (unit_sign == 1) and (not even_factors) and (not eq_factors) and (not unknown_factors):
                raise ParametricRealFieldInconsistencyError("New constraint {} {} {} is not satisfied".format(lhs, op, rhs))
            if not self._allow_refinement:
                if len(even_factors) + len(lt_factors) + len(eq_factors) + len(unknown_factors) > 1:
                    raise ParametricRealFieldRefinementError("{} <= 0 has several new factors: {}".format(comparison, even_factors+lt_factors+eq_factors))
                if even_factors:
                    if unit_sign > 0:
                        try:
                            assert factor_sign(even_factors[0]) == 0
                        except FactorUndetermined:
                            pass
                        self.record_factor(even_factors[0], operator.eq)
                if lt_factors:
                    self.record_factor(lt_factors[0], op)
                if eq_factors:
                    if unit_sign > 0:
                        self.record_factor(eq_factors[0], op)
                    elif unit_sign < 0:
                        self.record_factor(-eq_factors[0], op)
                if unknown_factors:
                    if unit_sign > 0:
                        self.record_factor(unknown_factors[0], op)
                    elif unit_sign < 0:
                        self.record_factor(-unknown_factors[0], op)
            else:
                if len(unknown_factors) > 1:
                    raise NotImplementedError()
                for new_fac in even_factors:
                    try:
                        if factor_sign(new_fac) == 0:
                            eq_factors.append(new_fac)
                    except FactorUndetermined:
                        eq_factors.append(new_fac) #??? or pass?
                        #This doesn't matter for now because under this branch  self._allow_refinement=True, and testpoint value is not removed, so  the exception FactorUndetermined should not happen.  In the future if the code in big_cells_impl changes and this exception happens, it is still safer to eq_factors.append(new_fac). This could result in smaller cell (which is allowed). With "pass" it may raise error in corner cases, such as when the removed testpoint value was the only eq factor without which there would be a sign contradiction.
                for new_fac in lt_factors:
                    self.record_factor(new_fac, operator.le)
                if not self._big_cells:
                    for new_fac in eq_factors:
                        self.record_factor(new_fac, operator.eq)
                else: #self._big_cells is True, self._allow_refinement is True.
                    undecided_eq = []
                    for new_fac in eq_factors:
                        if self.is_factor_known(new_fac, operator.le):
                            unit_sign = -unit_sign
                        elif not self.is_factor_known(-new_fac, operator.le):
                            undecided_eq.append(new_fac) #potentially record new_fac >=0, keep the sign
                    if not undecided_eq:
                        if unit_sign > 0:
                            self.record_factor(new_fac, operator.eq) #overwrite last le factor to eq
                    else:
                        if unit_sign < 0:
                            self.record_factor(-undecided_eq[0], operator.le)
                        else:
                            self.record_factor(undecided_eq[0], operator.le)
                        for new_fac in undecided_eq[1::]:
                            self.record_factor(-new_fac, operator.le)
            logging.debug("New element in %s._le: %s" % (repr(self), comparison))
            self._le.add(comparison)
        else:
            raise NotImplementedError("Not implemented operator: {}".format(op))

    def is_factor_known(self, fac, op):
        for bsa in (self._factor_bsa, self._bsa):
            try:
                return bsa.is_polynomial_constraint_valid(fac, op)
            except NotImplementedError:
                pass
        return False

    def record_factor(self, fac, op):
        if not self.is_factor_known(fac, op):
            if op == operator.lt:
                formatted_constraint = "%s < 0" % fac
            elif op == operator.eq:
                formatted_constraint = "%s == 0" % fac
            elif op == operator.le:
                formatted_constraint = "%s <= 0" % fac
            else:
                raise ValueError("{} is not a supported operator".format(op))
            if self._frozen:
                raise ParametricRealFieldFrozenError("Cannot prove that constraint is implied: {} ".format(formatted_constraint))
            self._bsa.add_polynomial_constraint(fac, op)
            self._factor_bsa.add_polynomial_constraint(fac, op)
            logging.info("New constraint: {}".format(formatted_constraint))
            if debug_new_factors:
                import pdb
                pdb.set_trace()

    def make_proof_cell(self, **opt):
        r"""
        Make a :class:`SemialgebraicComplexComponent` from a :class:`ParametricRealField`.
        
        In **opt, one can provide: region_type, function, find_region_type, default_var_bound, bddbsa, kwds_dict.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)

            sage: def foo(x,y):
            ....:     return (x+y < 2) and (y^2 < x)
            sage: K.<x,y> = ParametricRealField([1,1/2])
            sage: region_type = foo(*K.gens())
            sage: c = K.make_proof_cell(region_type=region_type)
            sage: sorted(c.bsa.lt_poly())
            [x + y - 2, y^2 - x]
            sage: c.plot()    #not tested

            sage: K.<f,bkpt> = ParametricRealField([1/5,3/7])
            sage: h = drlm_backward_3_slope(f, bkpt)
            sage: region_type = find_region_type_igp(K, h)
            sage: c1 = K.make_proof_cell(region_type=region_type)
            sage: c1.plot()  # not tested
            sage: sorted(c1.bsa.lt_poly())
            [-2*f + 3*bkpt - 1, f - 3*bkpt + 1, 2*f - bkpt]
            sage: c2 = K.make_proof_cell(region_type=region_type, function=h, find_region_type=None)
        """
        function = opt.pop('function', None)
        find_region_type = opt.pop('find_region_type', return_result)
        region_type = opt.pop('region_type', True)
        return SemialgebraicComplexComponent(self, region_type)

    def plot(self, *options, **kwds):
        return self.make_proof_cell().plot(*options, **kwds)


###############################
# TO REFACTOR using section.
###############################
def find_polynomial_map(eqs=[], poly_ring=None):
    """
    BAD FUCNTION! It is used in 'mathematica' approach for non-linear case. Can we avoid it?
    Return a polynomial map that eliminates linear variables in eqs, and a dictionary recording which equations were used to eliminate those linear variables.
    Assume that gaussian elimination has been performed by PPL.minimized_constraints() on the input list of equations eqs.
    It is only called in SemialgebraicComplex.add_new_component in the case polynomial_map is not provided but bddbsa has equations.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.WARN)
        sage: P.<a,b,c>=QQ[]
        sage: find_polynomial_map([a+2*b+1])
        [a, -1/2*a - 1/2, c]
        sage: find_polynomial_map([2*a+b, b+c-1/2])
        [a, -2*a, 2*a + 1/2]
        sage: find_polynomial_map([a*a-b, b*b-c*c])
        [a, a^2, c]
        sage: find_polynomial_map([a^2+4*a+1])
        [a, b, c]
        sage: P.<d>=QQ[]
        sage: find_polynomial_map([1-d^3])
        [d]
        sage: P.<f,a,b> = QQ[]
        sage: find_polynomial_map([-f*a + 3*a^2 + 3*f*b - 9*a*b, -11*b + 2, -11*f + 3, -11*a + 1])
        [3/11, 1/11, 2/11]
    """
    if (not eqs) and (not poly_ring):
        raise ValueError("Neither eqs nor poly_ring is not provided.")
    if not eqs:
        return list(poly_ring.gens())
    if not poly_ring:
        poly_ring = eqs[0].parent()
    # assert all(l.parent() == poly_ring for l in eqs)
    variables = list(poly_ring.gens())
    #independent_variables = set(variables)
    polynomial_map = copy(variables)
    n = len(variables)
    m = len(eqs)
    #FIXME: perhaps for i in range(m-1, -1, -1) and for j in range(n-1, -1, -1) is better, because ppl eliminates variables with larger indices in the inequalities, also because we want to keep f as variable.
    for i in range(m-1, -1, -1):
        for j in range(n-1, -1, -1):
            v = variables[j]
            if (eqs[i] // v).degree() == 0: # v is a linear variable in eqs[i].
                coef = eqs[i].monomial_coefficient(v) # type is rational
                #independent_variables.remove(v) # eliminate v
                v_mapped_to = v - eqs[i] / coef
                polynomial_map[j] = v_mapped_to
                for jj in range(n):
                    if j != jj:
                        polynomial_map[jj] = polynomial_map[jj].subs({v:v_mapped_to})
                for ii in range(m):
                    if i != ii:
                        eqs[ii] = eqs[ii].subs({v:v_mapped_to})
                break
    #P = PolynomialRing(poly_ring.base_ring(), list(independent_variables))
    #return [P(p) for p in polynomial_map]
    # didn't take sub-polynomialring to make the trick ineq orthogonal to eliminated variables work.
    return polynomial_map

######################################
# Functions with ParametricRealField K
######################################

from sage.misc.sageinspect import sage_getargspec, sage_getvariablename

def read_default_args(function, **opt_non_default):
    r"""
    Return the default values of arguments of the function.

    Override the default values if opt_non_default is given.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: read_default_args(gmic)
        {'conditioncheck': True, 'f': 4/5, 'field': None}
        sage: read_default_args(drlm_backward_3_slope, **{'bkpt': 1/5})
        {'bkpt': 1/5, 'conditioncheck': True, 'f': 1/12, 'field': None}
    """
    if function is None:
        return {}
    try:
        default_args = dict(function.default_values())
    except AttributeError:
        args, varargs, keywords, defaults = sage_getargspec(function)
        default_args = {}
        if defaults is not None:
            for i in range(len(defaults)):
                default_args[args[-i-1]]=defaults[-i-1]
    for (opt_name, opt_value) in opt_non_default.items():
        if opt_name in default_args:
            default_args[opt_name] = opt_value
    return default_args

###########################################
# Proof cells and proof complex:
###########################################
class SemialgebraicComplexComponent(SageObject):    # FIXME: Rename this to be more specific
    r"""
    A proof cell for parameter space analysis.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.WARN)
        sage: def foo(x,y):
        ....:     return (x+y < 2) and (y^2 < x)
        sage: K.<x,y> = ParametricRealField([1,1/2])
        sage: region_type = foo(*K.gens())
        sage: component = SemialgebraicComplexComponent(K, region_type)
        sage: list(component.bsa.lt_poly())
        [x + y - 2, y^2 - x]
        sage: component.plot(xmin=0, xmax=4, ymin=-3, ymax=3)        # not tested
        sage: new_points = component.find_neighbour_candidates(1/4, 'heuristic', goto_lower_dim=False, pos_poly=None)
        sage: new_pts = sorted(new_points[0].keys())
        sage: new_pts
        [(30093/42985, 54831/63571), (11/8, 7/8)]
        sage: len(new_pts)
        2
        sage: all((xx + yy - 2)*(yy^2 - xx) < 0 for (xx, yy) in new_pts)
        True
        sage: more_new_points = component.find_neighbour_candidates(1/4, 'mathematica', goto_lower_dim=True, pos_poly=None)  # optional - mathematica
        sage: sorted(more_new_points[0].keys())  # optional - mathematica
        [(17/8, 0), (2065/512, -33/16)]
        sage: sorted(more_new_points[1].keys())  # optional - mathematica
        [(0, 0), (2, 0)]

        sage: K.<x,y> = ParametricRealField([1,1/2])
        sage: region_type = foo(*K.gens())
        sage: x == 2*y
        True
        sage: sorted(K._bsa.eq_poly()), sorted(K._bsa.lt_poly()), sorted(K._bsa.le_poly())
        ([x - 2*y], [3*x - 4, y^2 - x], [])

    Without using BasicSemialgebraicSet_local, we got::

        sage: xx = K._bsa.poly_ring().gens()[0]; polynomial_map = [xx, 1/2*xx]
        sage: component = SemialgebraicComplexComponent(K, region_type, polynomial_map=polynomial_map)
        sage: sorted(component.bsa.lt_poly())
        [3*x - 4, y^2 - x]

    With using BasicSemialgebraicSet_local::

        sage: component = SemialgebraicComplexComponent(K, region_type)
        sage: sorted(component.bsa.lt_poly())
        [-x, 3*x - 4]

    In ProofCell region_type should alreay consider the polynomial_map::

        sage: K.<x,y> = ParametricRealField([1,1/2])
        sage: assert(x == 2*y)
        sage: region_type = foo(K.gens()[0], 1/2*K.gens()[0])
        sage: sorted(K._bsa.eq_poly()), sorted(K._bsa.lt_poly()), sorted(K._bsa.le_poly())
        ([x - 2*y], [-x, 3*x - 4], [])
        sage: xx = K._bsa.poly_ring().gens()[0]; polynomial_map = [xx, 1/2*xx]
        sage: component = SemialgebraicComplexComponent(K, region_type, polynomial_map=polynomial_map)
        sage: sorted(component.bsa.lt_poly())
        [-x, 3*x - 4]
    """

    def __init__(self, K, region_type, bddbsa=None, polynomial_map=None):
        self.var_name = list(K.variable_names()) #same as [ str(g) for g in K._bsa.poly_ring().gens() ]
        poly_ring = PolynomialRing(QQ, self.var_name)  #same as K._bsa.poly_ring()
        self.var_value = K._values # var_value was input of __init__. It doesn't seem necessary.
        self.region_type = region_type
        #Default bddbsa=None and polynomial_map=None are just for the convenience of doctests.
        #This doesn't happen when adding new components to complex.
        # # can use bddbsa of class ppl bsa later if poly_ring is implemented there. # Previously BasicSemialgebraicSet_eq_lt_le_sets class bddbsa was used.
        if polynomial_map is None: # for doctest purpose
            polynomial_map = list(poly_ring.gens())
        if bddbsa is None:   # for doctest purpose
            bddbsa = BasicSemialgebraicSet_veronese(poly_ring=poly_ring)
        self.bsa = copy(K._bsa)
        self.bddbsa = bddbsa
        # Polynomial map equalities and bddbsa were asserted in K in ProofCell.__init__. The inequalities in self.bsa and bddbsa should already been simplified through section with polynomial_map. This was done in _construct_field_and_test_point.
        # In lower dim proof cell or non-linear equations case, some equations of K._bsa are not presented in polynomial_map.
        eqs = list(K._bsa.eq_poly())
        if not all(l(polynomial_map) == 0 for l in eqs):
            polynomial_map = find_polynomial_map(eqs=eqs, poly_ring=poly_ring)
            #self.bsa = K._bsa.section(polynomial_map, bsa_class='veronese', poly_ring=poly_ring)  # this is a bigger_bsa
            self.bsa = BasicSemialgebraicSet_veronese.from_bsa(BasicSemialgebraicSet_local(K._bsa.section(polynomial_map, poly_ring=poly_ring), self.var_value)) # TODO:, polynomial_map=list(poly_ring.gens()))
            # WHY is this input polynomial_map sometimes not compatible with the variable elimination done in bddbsa? Because upstairs ppl bsa eliminates large x_i in the inequalities, and x_i doesn't necessarily correspond to the i-th variable in poly_ring. Since polynomial_map and v_dict were not given at the initialization of veronese, the variable first encounted in the constraints is considered as x0 by upstairs ppl bsa. # In old code, we fixed the order of upstairs variables by adding initial space dimensions. We don't do that in the current code. Instead, we take the section of bddbsa to eliminate the varibles in the equations. # Is the given bddbsa required to be veronese with upstairs being ppl_bsa? Convert it anyway. # It's the same as BasicSemialgebraicSet_veronese.from_bsa(bddbsa.section(self.polynomial_map), poly_ring=poly_ring)
            self.bddbsa = BasicSemialgebraicSet_veronese.from_bsa(BasicSemialgebraicSet_local(bddbsa.section(polynomial_map, poly_ring=poly_ring), self.var_value))
            # Taking section forgets the equations. Then add back the equations  # Finally self.bsa should be the same as K._bsa, but its inequalities don't have variables eliminated by polynomial map, so that heuristic wall crossing can be done later.
            for i in range(len(self.var_name)):
                if polynomial_map[i] != poly_ring.gens()[i]:
                    l = polynomial_map[i]-poly_ring.gens()[i]
                    (self.bsa).add_polynomial_constraint(l, operator.eq)
                    (self.bddbsa).add_polynomial_constraint(l, operator.eq)
        self.polynomial_map = polynomial_map
        self.neighbour_points = []  # just for plotting

    def __repr__(self):
        s = "SemialgebraicComplexComponent(var_value={}, region_type={}".format(self.var_value, self.region_type)
        if self.is_polyhedral():
            s += ", polyhedral"
        if self.bsa.eq_poly():
            s += ", {} eqs".format(len(list((self.bsa.eq_poly()))))
        if self.bsa.lt_poly():
            s += ", {} strict-ins".format(len(list(self.bsa.lt_poly())))
        if self.bsa.le_poly():
            s += ", {} nonstrict-ins".format(len(list(self.bsa.le_poly())))
        s += ")"
        return s

    def plot(self, alpha=0.5, plot_points=300, slice_value=None, default_var_bound=None, show_testpoints=True, goto_lower_dim=False, zorder=0, **kwds):
        r"""
        Plot the cell.

        - If slice_value is given, it is either a polynomial_map that defines a section, or a list of fixed parameter values with two of them being None. Plot the section. See examples in ``SemialgebraicComplex.plot()``.
        - show_testpoints controls whether to plot the testpoint in this cell.
        - plot_points controls the quality of the plotting.
        - goto_lower_dim controls whether strict and non-strict inequalies are plotted with different colors.
        """
        if not default_var_bound:
            default_var_bound = [None, None]
        xmin = kwds.pop('xmin', default_var_bound[0])
        xmax = kwds.pop('xmax', default_var_bound[1])
        ymin = kwds.pop('ymin', default_var_bound[0])
        ymax = kwds.pop('ymax', default_var_bound[1])
        color = kwds.pop('color', find_region_color(self.region_type))
        #bsa = self.bsa.intersection(self.bddbsa) # bsa_class is 'intersection'
        #not needed any more because bddbsa is recorded in self.bsa
        if not goto_lower_dim:
            g = (self.bsa).plot(alpha=alpha, plot_points=plot_points, slice_value=slice_value, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, color=color, fill_color=color, zorder=3*zorder, **kwds)
        elif not list(self.bsa.eq_poly()):
            # plot border with color white, then add non-strict inequalies.
            g = (self.bsa).plot(alpha=alpha, plot_points=plot_points, slice_value=slice_value, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, color='white', fill_color=color, zorder=3*zorder, **kwds)
            for l in self.bsa.le_poly():
                new_bsa = copy(self.bsa)
                new_bsa.add_polynomial_constraint(l, operator.eq)
                g += new_bsa.plot(alpha=alpha, plot_points=plot_points, slice_value=slice_value, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, color=color, fill_color=color, zorder=3*zorder+2, **kwds)
        else:
            # plot border with color, then add strict inequalies with color white, then add non-strict inequalies with color.
            g = (self.bsa).plot(alpha=alpha, plot_points=plot_points, slice_value=slice_value, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, color=color, fill_color=color, zorder=3*zorder, **kwds)
            for l in self.bsa.lt_poly():
                new_bsa = BasicSemialgebraicSet_eq_lt_le_sets(eq=list(self.bsa.eq_poly())+[l], lt=[ll for ll in self.bsa.lt_poly() if ll != l], le=list(self.bsa.le_poly()))
                g += new_bsa.plot(alpha=alpha, plot_points=plot_points, slice_value=slice_value, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, color='white', fill_color='white', zorder=3*zorder+1, **kwds)
            for l in self.bsa.le_poly():
                new_bsa = copy(self.bsa)
                new_bsa.add_polynomial_constraint(l, operator.eq)
                g += new_bsa.plot(alpha=alpha, plot_points=plot_points, slice_value=slice_value, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, color=color, fill_color=color, zorder=3*zorder+2, **kwds)
        if show_testpoints and not slice_value:
            pt = self.var_value
            if color == 'white':
                ptcolor = 'black'
            else:
                ptcolor = 'white'
            if (xmin <= pt[0] <= xmax) and (ymin <= pt[1] <= ymax):
                g += point(pt, color = ptcolor, size = 2, zorder=10)
        return g

    def find_neighbour_candidates(self, flip_ineq_step, wall_crossing_method='heuristic', goto_lower_dim=False, pos_poly=None):
        r"""
        Try flipping exactly one inequality at one time, to reach a new testpoint as neighbour candidate.

        - flip_ineq_step defines the step length
        - wall_crossing_method is 'heuristic' or 'mathematica' or 'heuristic_then_mathematica'
        - if goto_lower_dim=False, the cell is considered as its formal closure, so no recursion into test points in lower dimensional cells.
        - pos_poly is a polynomial. The return test point must satisfy pos_poly(new test point) > 0.

        OUTPUT new_points is a dictionary of dictionaries. The new_points[i] is a dictionay whose keys = candidate neighbour testpoints, values = (bddbsa whose eq_poly has i elements, polynomial_map, no_crossing_l) of the candidate neighbour cell that contains the candidate neighbour testpoint. bddbsa is recorded so that (1) complex.bddbsa is always respected; and (2) can recursively go into lower dimensional cells. polynomial_map is recorded and passed to the constructor of the neighbour cell. no_crossing is passed to the neighour cell for its find_neighbour_candidates method. We no longer update self.bsa by removing (obvious) redundant eq, lt, le constraints from its description at the end, even when 'mathematica' is used.
        """
        bsa_eq_poly = list(self.bsa.eq_poly())
        bsa_le_poly = list(self.bsa.le_poly())
        bsa_lt_poly = list(self.bsa.lt_poly())
        num_eq = len(bsa_eq_poly) #was len(list(self.bddbsa.eq_poly()))
        new_points = {}  #dictionary with key=num_eq, value=dictionay of pt: (bddbsa, polynomial_map).
        #bddbsa = copy(self.bddbsa)
        #for l in bsa_eq_poly: # should be already in bddbsa
        #    bddbsa.add_polynomial_constraint(l, operator.eq)
        # first flip l <= 0 to l > 0.
        for l in bsa_le_poly:
            bsa = copy(self.bddbsa) #veronese
            if pos_poly is not None:
                bsa.add_polynomial_constraint(pos_poly, operator.gt)
            # inequalites are orthogonal to the equations
            for ll in bsa_lt_poly:
                bsa.add_polynomial_constraint(ll, operator.lt)
            for ll in bsa_le_poly:
                if ll != l:
                    bsa.add_polynomial_constraint(ll, operator.le)
            # FIXME: should simply bsa, for example 4*f^2-4*f-1 <= 0 should be simplified to f = 1/2. More crucial for pt_on_wall. # now self.bddbsa and self.bsa were obtained through bsa_local, should be fine if bsa is now lower dimensional.
            try:
                if bsa.is_polynomial_constraint_valid(l, operator.le):
                    # can't flip l<=0, as bsa intersection with l > 0 is empty.
                    continue
            except NotImplementedError:
                pass
            pt_across_wall = None
            # find point in intersection(bsa, l>0, l<flip_ineq_step)
            if wall_crossing_method == 'heuristic' or wall_crossing_method is None or wall_crossing_method == 'heuristic_then_mathematica':
                # do not distinguish lt_poly and le_poly as we wish test_points are in general positions
                pt = find_point_flip_ineq_heuristic(self.var_value, l, list(bsa.lt_poly())+list(bsa.le_poly()), flip_ineq_step)
                if pt is not None:
                    # Find a new point, use polynomial map to recover the values of those eliminated variables.
                    pt_across_wall = tuple(p(pt) for p in self.polynomial_map)                 
            if wall_crossing_method == 'mathematica' or wall_crossing_method == 'heuristic_then_mathematica' and (pt_across_wall is None):
                bsa_mathematica = bsa.formal_relint(bsa_class='mathematica') # was BasicSemialgebraicSet_mathematica.from_bsa(bsa)
                bsa_mathematica.add_polynomial_constraint(l, operator.gt)
                bsa_mathematica.add_polynomial_constraint(l - flip_ineq_step, operator.lt) #needed? yes, otherwise often find (0,0)
                pt = bsa_mathematica.find_point()
                if pt is not None:
                    pt_across_wall = tuple(pt)
            if pt_across_wall is not None:
                if num_eq not in new_points:
                    new_points[num_eq] = OrderedDict()
                new_points[num_eq][pt_across_wall] = (self.bddbsa, self.polynomial_map, l)
        # now flip l < 0 to l >= 0
        for l in bsa_lt_poly:
            bsa = copy(self.bddbsa) #veronese
            if pos_poly is not None:
                bsa.add_polynomial_constraint(pos_poly, operator.gt)
            # inequalites are orthogonal to the equations
            for ll in bsa_le_poly:
                bsa.add_polynomial_constraint(ll, operator.le)
            for ll in bsa_lt_poly:
                if ll != l:
                    bsa.add_polynomial_constraint(ll, operator.lt)
            has_pt_across_wall = True
            has_pt_on_wall = True
            pt_across_wall = None
            pt_on_wall = None
            try:
                if bsa.is_polynomial_constraint_valid(l, operator.lt):
                    # can't flip l<0, as bsa intersection with l >= 0 is empty.
                    continue
                if bsa.is_polynomial_constraint_valid(l, operator.le):
                    # can only get l==0, but bsa intersection with l > 0 is empty.
                    has_pt_across_wall = False
            except NotImplementedError:
                pass
            # find point in intersection(bsa, l>0, l<flip_ineq_step)
            if has_pt_across_wall:
                if wall_crossing_method == 'heuristic' or wall_crossing_method is None or wall_crossing_method == 'heuristic_then_mathematica':
                    # do not distinguish lt_poly and le_poly as we wish test_points are in general positions
                    pt = find_point_flip_ineq_heuristic(self.var_value, l, list(bsa.lt_poly())+list(bsa.le_poly()), flip_ineq_step)
                    if pt is not None:
                        #find a new point, use polynomial map to recover the values of those eliminated variables.
                        pt_across_wall = tuple(p(pt) for p in self.polynomial_map)
                if wall_crossing_method == 'mathematica' or wall_crossing_method == 'heuristic_then_mathematica' and (pt_across_wall is None):
                    bsa_mathematica = bsa.formal_relint(bsa_class='mathematica') # was BasicSemialgebraicSet_mathematica.from_bsa(bsa)
                    bsa_mathematica.add_polynomial_constraint(l, operator.gt)
                    bsa_mathematica.add_polynomial_constraint(l - flip_ineq_step, operator.lt) #needed? yes
                    pt = bsa_mathematica.find_point()
                    if pt is not None:
                        pt_across_wall = tuple(pt)
                if pt_across_wall is not None:
                    if num_eq not in new_points:
                        new_points[num_eq] = OrderedDict()
                    new_points[num_eq][pt_across_wall] = (self.bddbsa, self.polynomial_map, l)
            # find point in intersection(bsa, l==0)
            if goto_lower_dim:
                if (wall_crossing_method == 'heuristic' or wall_crossing_method is None or wall_crossing_method == 'heuristic_then_mathematica') and (l.degree() == 1):  # FIXME: too restrictive. Consider rational map. Take section then find testpoint then map back.
                    current_point = vector(self.var_value)
                    # pt_on_wall satisfies one more equation l == 0 than the points in self.bsa
                    # eliminate v in strict_ineqs and nonstrict_ineqs
                    eq = l
                    polynomial_map = self.polynomial_map
                    bsa_section = bsa
                    # there was bug example: igp.big_cells_default=True; cpl r=regions[13].
                    while has_pt_on_wall:
                        lhs = vector(eq.monomial_coefficient(x) for x in l.parent().gens())
                        cst = eq.constant_coefficient()
                        if cst == 0:
                            v = eq.monomials()[-1]
                        else:
                            v = eq.monomials()[-2] # order??
                        v_mapped_to = v - eq / (eq.monomial_coefficient(v))  # eliminate v
                        polynomial_map = [p.subs({v: v_mapped_to}) for p in polynomial_map]
                        bsa_section = bsa_section.section(polynomial_map, bsa_class='veronese', poly_ring=l.parent())
                        # There was bug example. Fixed by asking bsa_section.upstairs()._polyhedron.is_empty().
                        # sage: lt = [-f + 2*z, 2*f^2 + 2*f*z - 6*z^2 - f + z]
                        # sage: le = [2*f + 3*z - 1, -3*f - 9*z + 2, -6*f*z - 18*z^2 + 2*f + 7*z - 1]
                        # sage: l = 3*f - 1
                        # There was bug example, CPl bigcell Cell 16
                        # sage: bsa = BasicSemialgebraicSet_eq_lt_le_sets(lt=[2*f - 1, -9*f + 2],le=[4*f^2 - 4*f + 1])
                        # sage: veronese = BasicSemialgebraicSet_veronese.from_bsa(bsa)
                        # sage: veronese.upstairs()._polyhedron.is_empty() #False, but veronese is empty.
                        # Should simply bsa, for example 4*f^2-4*f-1 <= 0 should be simplified to f = 1/2.
                        if bsa_section.upstairs()._polyhedron.is_empty():
                            has_pt_on_wall = False
                        else:
                            lts = []; les = []; eqs = []
                            for ll in list(bsa_section.lt_poly()):
                                factors = ll.factor()
                                if len(factors) == 1:
                                    (fac, d) = factors[0]
                                    if d % 2 == 0:
                                        if factors.unit() > 0:
                                            has_pt_on_wall = False
                                            break
                                        else:
                                            lts.append(-fac * fac)
                                    else:
                                        lts.append(factors.unit() * fac)
                                else:
                                    lts.append(ll)
                            if not has_pt_on_wall:
                                break
                            for ll in list(bsa_section.le_poly()):
                                factors = ll.factor()
                                if len(factors) == 1:
                                    (fac, d) = factors[0]
                                    if d % 2 == 0:
                                        if factors.unit() > 0:
                                            eqs.append(fac)
                                        else:
                                            pass # - fac ^ even <= 0 always true.
                                    else: # d % 1 == 0:
                                        les.append(factors.unit() * fac)
                                else:
                                    les.append(ll) # getting rid of even degrees may break cell.
                            for ll in list(bsa_section.eq_poly()):
                                factors = ll.factor()
                                eqs.append(factors.unit() * product([fac for (fac, d) in factors]))
                            bsa_section = BasicSemialgebraicSet_eq_lt_le_sets(lt=lts, le=les, eq=eqs, poly_ring=l.parent())
                            if eqs:
                                ###import pdb; pdb.set_trace()
                                eq = eqs[0]
                            else:
                                break
                    # if list(bsa_section.eq_poly()):
                    #     import pdb; pdb.set_trace() # should update polynomial_map. Do this later, in bddbsa for new point.
                    #strict_ineqs = [ll.subs({v: v_mapped_to}) for ll in bsa.lt_poly()]
                    #nonstrict_ineqs = [ll.subs({v: v_mapped_to}) for ll in bsa.le_poly() if ll.subs({v: v_mapped_to}) != 0]
                    # l is linear wall. try to find a point on the wall using heuristic gradient descent method.
                    if has_pt_on_wall:
                        current_point -= (l(*current_point)*vector(gradient(l))) / (vector(gradient(l))*vector(gradient(l)))
                        pt = adjust_pt_to_satisfy_ineqs(vector(RR(x) for x in current_point), None, list(bsa_section.lt_poly()), list(bsa_section.le_poly()), flip_ineq_step)
                        #pt = find_point_on_ineq_heuristic(self.var_value, l, strict_ineqs, nonstrict_ineqs, flip_ineq_step)
                        if pt is not None:
                            # restrict the bddbsa of the candidate neighbour cell.
                            # univariate polynomials (Polynomial_rational_flint) does not define "coefficient", but has "monomial_coefficient".
                            new_bddbsa = copy(self.bddbsa)
                            new_bddbsa.add_linear_constraint(lhs, cst, operator.eq)
                            eqs = list(new_bddbsa.eq_poly())
                            new_num_eq = len(eqs)
                            if new_num_eq != num_eq + 1:
                                #import pdb; pdb.set_trace()
                                polynomial_map = find_polynomial_map(eqs=eqs)
                            pt_on_wall = tuple(p(pt) for p in polynomial_map)
                if wall_crossing_method == 'mathematica' or wall_crossing_method == 'heuristic_then_mathematica' and (has_pt_on_wall is True) and (pt_on_wall is None):
                    bsa.add_polynomial_constraint(l, operator.eq) # add eq first then take formal_relint. Okay to change bsa here.
                    bsa_mathematica = bsa.formal_relint(bsa_class='mathematica') # was BasicSemialgebraicSet_mathematica.from_bsa(bsa)
                    pt = bsa_mathematica.find_point()
                    if pt is not None:
                        pt_on_wall = tuple(pt)
                    new_bddbsa = copy(self.bddbsa)
                    new_bddbsa.add_polynomial_constraint(l, operator.eq)
                    eqs = list(new_bddbsa.eq_poly())
                    new_num_eq = len(eqs)
                    # Although mathematica doesn't care about polynomial_map, it is good record the section map, so as to simplify K and to make len(list(new_component.bsa.eq_poly())) accurate. (for example, [lambda_1 - 1, 4*f - 3, 4*f*lambda_1 - 3] would be [lambda_1 - 1, 4*f - 3]. This example is generated by sage: complex = SemialgebraicComplex(gj_2_slope, ['f','lambda_1']); complex.bfs_completion(var_value=[3/5,1/6],wall_crossing_method='mathematica',goto_lower_dim=True)); c = list(complex.cells_containing_point((3/4,1)))[0]; sorted(c.bsa.eq_poly())
                    polynomial_map = find_polynomial_map(eqs=eqs)
                if pt_on_wall is not None:
                    if new_num_eq not in new_points:
                        new_points[new_num_eq] = OrderedDict()
                    new_points[new_num_eq][pt_on_wall] = (new_bddbsa, polynomial_map, None)
        return new_points

    def is_polyhedral(self):
        for l in list(self.bsa.eq_poly()) + list(self.bsa.lt_poly()) + list(self.bsa.le_poly()):
            if l.degree() > 1:
                return False
        return True

    def sage_polyhedron(self):
        if not self.is_polyhedral():
            raise NotImplementedError("The cell is not polyhedral. Construct the polyhedron in the lifted space.")
        ieqs = [ [-l.constant_coefficient()]+[-l.monomial_coefficient(m) for m in l.args()] for l in list(self.bsa.lt_poly()) + list(self.bsa.le_poly()) ]
        eqns = [ [-l.constant_coefficient()]+[-l.monomial_coefficient(m) for m in l.args()] for l in list(self.bsa.eq_poly())]
        return Polyhedron(ieqs=ieqs, eqns=eqns)

def _max_x_y(x, y):    # A global, for testing pickling
    return max(x, y)

class ProofCell(SemialgebraicComplexComponent, Classcall):

    r"""
    A proof cell for parameter space analysis.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import ProofCell, ParametricFamily, _max_x_y, result_symbolic_expression
        sage: import logging
        sage: logging.disable(logging.WARN)
        sage: C = ProofCell(ParametricFamily(_max_x_y), [1, 2], find_region_type=result_symbolic_expression)
        sage: sorted(C.bsa.lt_poly())
        [x - y]

    Proof cells remember their construction::

        sage: C._init_args
        (<class 'cutgeneratingfunctionology.igp.ProofCell'>,
         (ParametricFamily(_max_x_y, ...),
          (1, 2),
          <...result_symbolic_expression...>,
          ...),
         {})

    Proof cells can be pickled if ``family`` is a global variable or can be pickled otherwise
    (lambda functions cannot be pickled)::

        sage: p_C = dumps(C)
        sage: explain_pickle(p_C)      # not tested
        sage: C_copy = loads(p_C)
        sage: sorted(C_copy.bsa.lt_poly())
        [x - y]
        sage: C_copy._init_args == C._init_args
        True

    (We do not test for equality C_copy == C -- we have not even decided yet what the semantics 
    of equality of bsa is.)
    """

    @staticmethod
    def __classcall__(cls, family, var_value, find_region_type, bddbsa=None, polynomial_map=None):
        if not isinstance(family, ParametricFamily_base):
            family = ParametricFamily(family)
        var_value = tuple(var_value)
        return super(ProofCell, cls).__classcall__(cls, family, var_value, find_region_type, bddbsa, polynomial_map)

    @staticmethod
    def _construct_field_and_test_point(family, var_value, polynomial_map=None):
        r"""
        Construct a :class:`ParametricRealField` K using ``var_value`` and ``polynomial_map``.

        ``var_value`` is a list parallel to ``family.names()``.

        Construct a test_point of type dictionary, which maps each parameter of the family to the corresponding :class:`ParametricRealFieldElement` if this is a parameter of K, otherwise maps to the default argument value.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)
            sage: family=ParametricFamily(gmic); var_value=[1/2]; P.<f>=QQ[]; polynomial_map = [f]
            sage: K, pt, test_point = ProofCell._construct_field_and_test_point(family, var_value, polynomial_map)
            sage: K
            ParametricRealField(names = ['f'], values = [1/2])
            sage: pt
            [f~]
            sage: test_point
            {'conditioncheck': False,
                 'f': f~,
                 'field': ParametricRealField(names = ['f'], values = [1/2])}

        Bug example::

            sage: P.<f,z>=QQ[]
            sage: theta = [3*f/(2*f+2), 3*f/(2*f+2)]
            sage: family = cpl_n_group_function(3, True, theta=lambda x, y: tuple([t(x, *y) for t in theta]))
            sage: var_value = (111/1000, 111/1000)
            sage: K, pt, test_point = ProofCell._construct_field_and_test_point(family, var_value)
            sage: list(K._bsa.eq_poly())
            []
            sage: pt
            [f~, z~]
            sage: polynomial_map = [f, f]
            sage: K, pt, test_point = ProofCell._construct_field_and_test_point(family, var_value, polynomial_map)
            sage: list(K._bsa.eq_poly())
            [f - z]
            sage: pt
            [f~, f~]
        """
        var_name = family.names()
        K = ParametricRealField(var_value, var_name)
        if polynomial_map is None:
            pt = K.gens()
        else:
            pt = [p(*K.gens()) for p in polynomial_map]
            for i in range(len(K.gens())):
                assert K.gens()[i] == polynomial_map[i](*K.gens())  # record equations in K
        test_point = family.args_from_point(pt)
        args_set = family.default_values()
        if 'field' in args_set:
            test_point['field'] = K
        if 'conditioncheck' in args_set:
            test_point['conditioncheck'] = False
        if 'merge' in args_set:
            test_point['merge'] = False
        #if bddbsa is None:
        #    bddbsa =  BasicSemialgebraicSet_eq_lt_le_sets(poly_ring=PolynomialRing(QQ, var_name)) # class doesn't matter, just to record the constraints on K._bsa
        #The code seems to still work after commenting out the following.
        #K.add_initial_space_dim() #so that the parameters var_name are the first ones in the monomial list. Needed for ppl variable elimination, otherwise ineqs are not orthogonal to eliminated variables. FIXME: Get rid of this requirement. Also needed for Mccormicks?
        #assert (K.gens() in bddbsa) # record boundary constraints. # Move to __init__
        return K, pt, test_point

    def __init__(self, family, var_value, find_region_type, bddbsa, polynomial_map):
        """
        Bug example::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)
            sage: P.<f,z>=QQ[]
            sage: theta = [3*f/(2*f+2), 3*f/(2*f+2)]
            sage: family = cpl_n_group_function(3, True, theta=lambda x, y: tuple([t(x, *y) for t in theta]))
            sage: var_value = (111/1000, 111/1000)
            sage: find_region_type = find_region_type_igp_extreme_big_cells
            sage: bddbsa = BasicSemialgebraicSet_veronese(poly_ring=P)
            sage: bddbsa.add_polynomial_constraint(f - z, operator.eq)
            sage: bddbsa.add_polynomial_constraint(8*f - 1, operator.lt)
            sage: bddbsa.add_polynomial_constraint(-f, operator.lt)
            sage: polynomial_map = [f, f]
            sage: cell = ProofCell(family, var_value, find_region_type, bddbsa, polynomial_map)
            sage: sorted(cell.bsa.lt_poly()), sorted(cell.bsa.eq_poly())  # was ([-11*f + 1, 9*f - 1, -11*f^2 + f, 5*f^2 - f], [f - z]))
            ([-11*f + 1, 9*f - 1], [f - z])
        """
        K, pt, test_point = self._construct_field_and_test_point(family, var_value, polynomial_map)
        try:
            h = family(**test_point)
        except Exception:
            # Function is non-constructible at this random point.
            if debug_cell_exceptions:
                import pdb; pdb.post_mortem()
            h = None
        region_type = find_region_type(K, h)
        if bddbsa is not None:
            assert (pt in bddbsa) # record boundary constraints. #NEW: moved from ProofCell._construct_field_and_test_point
        super(ProofCell, self).__init__(K, region_type, bddbsa, polynomial_map)
        self.family = family

from collections import OrderedDict

class SemialgebraicComplex(SageObject):
    r"""
    A proof complex for parameter space analysis.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.WARN)

        sage: def vol(a,b):
        ....:     P = Polyhedron(ieqs=[(0,0,1),(0,1,0),(1,-1,0),(1,0,-1),(a,-1,0),(b,0,-1)])
        ....:     return P.volume()
        sage: complex = SemialgebraicComplex(vol, ['a','b'], find_region_type=result_symbolic_expression, default_var_bound=(-1,3))

    Breadth-first-search to complete the complex, starting at the point (a,b)=(2,1/2), using heuristic wall-crossing, considering full-dimensional cells only::

        sage: complex.bfs_completion(var_value=[2, 1/2], flip_ineq_step=1/1000, check_completion=False, wall_crossing_method='heuristic', goto_lower_dim=False)
        sage: sorted(c.region_type for c in complex.components)
        [(0,), (0,), (0,), (0,), (1,), (b,), (a,), (a*b,), (a*b,)]
        sage: complex.plot()                                  # not tested

    Instead of heuristic method, we can use Mathematica's ``FindInstance`` to look for uncovered points in wall-crossing::

        sage: complex = SemialgebraicComplex(vol, ['a','b'], find_region_type=result_symbolic_expression, default_var_bound=(-1,3))                         # optional - mathematica
        sage: complex.bfs_completion(var_value=[2, 1/2], flip_ineq_step=1/1000, check_completion=False, wall_crossing_method='mathematica', goto_lower_dim=True)        # optional - mathematica
        sage: len(complex.components)                         # optional - mathematica
        25
        sage: complex.plot(goto_lower_dim=True)               # not tested

    The entire parameter space is covered by cells::

        sage: complex.is_complete()                           # optional - mathematica
        True
        
    Example with non-linear wall::

        sage: complex = SemialgebraicComplex(lambda x,y: max(x,y^2), ['x','y'], find_region_type=result_symbolic_expression, default_var_bound=(-3,3))  # optional - mathematica
        sage: complex.bfs_completion(var_value=[1,1/2], check_completion=True, wall_crossing_method='mathematica', goto_lower_dim=True)                             # optional - mathematica
        sage: len(complex.components)                         # optional - mathematica
        3
        sage: complex.components[0].region_type               # optional - mathematica
        (x,)
        sage: complex.components[1].region_type               # optional - mathematica
        (y^2,)
        sage: complex.plot()                                  # not tested


    Analyse the extreme/minimal/valid regions of a function in the Gomory-Johnson infinite group problem.

    Use random shooting method to complete the complex. See more options in the method ``shoot_random_points``::

        sage: complex = SemialgebraicComplex(drlm_backward_3_slope, ['f','bkpt'])
        sage: complex.shoot_random_points(50)

    Use breadth-first-search to complete the complex. See more options in the method ``bfs_completion``::

        sage: complex = SemialgebraicComplex(drlm_backward_3_slope, ['f','bkpt'])
        sage: complex.bfs_completion(var_value=[4/5,1/2])
        sage: len(complex.components)
        17
        sage: complex.is_polyhedral()
        True
        sage: complex.is_face_to_face()
        False
        sage: complex2 = complex.subcomplex_of_cells_with_given_region_types({'is_extreme','not_extreme','not_minimal'})
        sage: complex2.is_face_to_face()
        True
        sage: extc = complex.subcomplex_of_cells_with_given_region_types()
        sage: len(extc.components)
        2
        sage: extc.is_polyhedral()
        True
        sage: extc.is_face_to_face()
        True
        sage: boundary = extc.guess_boundary()
        sage: sorted(boundary, key=str)
        [-f, -f + 4*bkpt - 1, f - bkpt]
        sage: extc.bddbsa = BasicSemialgebraicSet_eq_lt_le_sets(lt=boundary)
        sage: extc.is_complete(formal_closure=False) # optional - mathematica
        False
        sage: extc.is_complete(formal_closure=True) # optional - mathematica
        True
        sage: pc = extc.polyhedral_complex()
        sage: sorted( sorted(pol.vertices_list()) for pol in pc.cells_list() )
        [[[0, 0]],
         [[0, 0], [0, 1/4]],
         [[0, 0], [0, 1/4], [1/7, 2/7]],
         [[0, 0], [1/7, 2/7]],
         [[0, 0], [1/7, 2/7], [1/3, 1/3]],
         [[0, 0], [1/3, 1/3]],
         [[0, 1/4]],
         [[0, 1/4], [1/7, 2/7]],
         [[1/7, 2/7]],
         [[1/7, 2/7], [1/3, 1/3]],
         [[1/3, 1/3]]]

    Consider low-dimensional cells in the polyhedral case::

        sage: complex = SemialgebraicComplex(drlm_backward_3_slope, ['f','bkpt'])   #long time
        sage: complex.bfs_completion(var_value=[4/5,1/2],goto_lower_dim=True)       #long time
        sage: len(complex.components)                                               #long time
        43
        sage: (sum(c.sage_polyhedron().plot() for c in complex.components if len(list(c.bsa.eq_poly())) == 1) + sum(c.sage_polyhedron().plot(color='red') for c in complex.components if len(list(c.bsa.eq_poly())) == 2)).show()    # not tested

        sage: extc = complex.subcomplex_of_cells_with_given_region_types()          #long time
        sage: len(extc.components)                                                  #long time
        6
        sage: boundary = extc.guess_boundary()                                      #long time
        sage: sorted(boundary)                                                      #long time
        [-f, -f + 4*bkpt - 1, f - bkpt]
        sage: extc.bddbsa = BasicSemialgebraicSet_eq_lt_le_sets(lt=boundary)        #long time
        sage: extc.is_complete()                                                    #long time, optional - mathematica
        True
    """
    def __init__(self, family, var_name=None, find_region_type=None, default_var_bound=(-0.1,1.1), bddbsa=None, polynomial_map=None, kwds_dict={}, cell_class=None, **opt_non_default):
        r"""
        Construct a SemialgebraicComplex.

        - ``family``:  A class or constructor, defining a parametric family of computations
        - var_name: a subset of the parameters of ``family``
        - ``find_region_type``: maps the object constructed by ``family`` to a type of the parameter region;
            - set to ``None`` for functions in Gomory-Johnson model; 
            - often set to ``result_symbolic_expression`` or ``result_concrete_value`` for other functions;
        - default_var_bound: need if we use random shooting method to complete the complex; If given, the bound is also used in plotting

        The following argument defines the boundary of the ``SemialgebraicComplex``, so that bfs won't go beyond the region. They might be useful in the CPL3 examples.

        - ``bddbsa``: a BasicSemialgebraicSet that contains the points in the complex;
        - ``polynomial_map``: need when bddbsa is lower dimensional.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.WARN)
            sage: complex = SemialgebraicComplex(lambda x,y: max(x,y), ['x','y'], find_region_type=result_symbolic_expression, default_var_bound=(-10,10))
        """
        #self.num_components = 0
        self.components = []
        if isinstance(family, ParametricFamily_base) and var_name is None and not kwds_dict and not opt_non_default:
            pass
        else:
            default_values = copy(opt_non_default)
            default_values.update(kwds_dict)
            family = ParametricFamily(family, default_values=default_values, names=var_name)
        self.family = family
        var_name = family.names()
        self.d = len(var_name)
        self.graph = Graphics()
        self.num_plotted_components = 0
        self.points_to_test = OrderedDict() # a dictionary of dictionaries of the form {testpoint: (bddbsa, polynomial_map)}. The dictionary with key i corresponds to the cells that are expected to have i equations.
        # we record bddbsa for each testpoint, since we want to restrict to lower dimensional cells when  goto_lower_dim is set to True.
        self.tested_points = OrderedDict()
        if find_region_type is None:
            find_region_type = find_region_type_igp
        self.find_region_type = find_region_type
        self.default_var_bound = default_var_bound
        if bddbsa is None:   #HAS BUG: r = regions[31]; theta = thetas[16]; cpl_complex = cpl_fill_region_given_theta(r, theta); self.bddbsa = BasicSemialgebraicSet_veronese(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(Constraint_System {x0+6*x1-1==0, -x0+1>0, 2*x0-1>0}), polynomial_map=[f, z]); self.bddbsa.polynomial_map() is [f, z]; self.bddbsa.ambient_dim() is 1.
            poly_ring = PolynomialRing(QQ, var_name)
            self.bddbsa = BasicSemialgebraicSet_veronese(poly_ring=poly_ring)
        else:
            self.bddbsa = bddbsa
        if polynomial_map is None:
            eqs = list(self.bddbsa.eq_poly())
            if eqs:
                #Only useful for treating careless input with low dim bddbsa provided but not polynomial_map. # should not happen in cpl because we pass r.polynomal_map to cpl_complex.
                #import pdb; pdb.set_trace()
                self.polynomial_map = find_polynomial_map(eqs, poly_ring=PolynomialRing(QQ, var_name))
            else:
                self.polynomial_map = list(self.bddbsa.poly_ring().gens())
        else:
            self.polynomial_map = polynomial_map
        if cell_class is None:
            cell_class = ProofCell
        self._cell_class = cell_class

    def __repr__(self):
        return "SemialgebraicComplex with {} components".format(len(self.components))

    def generate_random_var_value(self, var_bounds=None):
        r"""
        Return a random point that satisfies var_bounds and is in self.bsa.

        - If var_bounds is not specified, self.default_var_bound is taken. 
        - var_bounds can be a list of 2-tuples whose length equals to the number of parameters, or lambda functions.
        - It is used in random shooting method for functions like ``dg_2_step_mir``, which involve floor/ceil operations. We try to plot one layer for each n = floor(...) and superimpose the layers at the end to get the whole picture.

        Notice that if self.bddbsa.eq_poly() is not empty, it never ends.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.WARN)
            sage: complex = SemialgebraicComplex(lambda x,y: max(x,y), ['x','y'], find_region_type=result_symbolic_expression, default_var_bound=(-10,10))
            sage: random_point = complex.generate_random_var_value()
            sage: all(complex.default_var_bound[0] < random_point[i] < complex.default_var_bound[1] for i in [0,1])
            True
            sage: var_bounds=[(0,1),((lambda x: x/2), (lambda x: x))]
            sage: random_point = complex.generate_random_var_value(var_bounds=var_bounds)
            sage: (0 < random_point[0] < 1) and (random_point[0]/2 < random_point[1] < random_point[0])
            True
        """
        while True:
            var_value = []
            for i in range(self.d):
                if not var_bounds:
                    x = QQ(uniform(self.default_var_bound[0], self.default_var_bound[1]))
                else:
                    if hasattr(var_bounds[i][0], '__call__'):
                        l =  var_bounds[i][0](*var_value)
                    else:
                        l = var_bounds[i][0]
                    if hasattr(var_bounds[i][1], '__call__'):
                        u =  var_bounds[i][1](*var_value)
                    else:
                        u = var_bounds[i][1]
                    if l > u:
                        x = None
                    else:
                        x = QQ(uniform(l, u))
                if x is None:
                    break
                var_value.append(x)
            # if random point is not in self.bddbsa, continue while.
            if (x is not None) and (var_value in self.bddbsa):
                return var_value

    def is_point_covered(self, var_value):
        r"""
        Return whether the given point var_value is contained in any cell of the complex.

        Inequalities are considered strict.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.WARN)
            sage: complex = SemialgebraicComplex(lambda x,y: max(x,y), ['x','y'], find_region_type=result_symbolic_expression, default_var_bound=(-10,10))
            sage: complex.add_new_component([1,2], bddbsa=None, flip_ineq_step=0, wall_crossing_method=None, goto_lower_dim=False) # the cell {(x,y): x<y}
            sage: complex.is_point_covered([2,3])
            True
            sage: complex.is_point_covered([3,2])
            False
            sage: complex.is_point_covered([2,2])
            False
            sage: complex.is_point_covered([sqrt(2), sqrt(3)])
            True
        """
        for c in self.cells_containing_point(var_value):
            return True
        return False

    def cells_containing_point(self, var_value):
        r"""
        Yield the cells of the complex that contain the given point var_value.

        Inequalities are considered strict.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.WARN)
            sage: complex = SemialgebraicComplex(lambda x,y: max(x,y), ['x','y'], find_region_type=result_symbolic_expression, default_var_bound=(-10,10))
            sage: complex.add_new_component([1,2], bddbsa=None, flip_ineq_step=0, wall_crossing_method=None, goto_lower_dim=False) # the cell {(x,y): x<y}
            sage: cell = next(complex.cells_containing_point([2,3]))
            sage: list(cell.var_value)
            [1, 2]
            sage: cell = next(complex.cells_containing_point([sqrt(2), sqrt(3)]))
            sage: list(cell.var_value)
            [1, 2]
        """
        for c in self.components:
            if var_value in c.bsa:
                yield c
        
    def find_uncovered_random_point(self, var_bounds=None, max_failings=10000):
        r"""
        Return a random point that satisfies the bounds and is uncovered by any cells in the complex.
        Return ``None`` if the number of attemps > max_failings.
 
        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.WARN)
            sage: complex = SemialgebraicComplex(lambda x,y: max(x,y), ['x','y'], find_region_type=result_symbolic_expression, default_var_bound=(-10,10))
            sage: complex.add_new_component([1,2], bddbsa=None, flip_ineq_step=0, wall_crossing_method=None, goto_lower_dim=False) # the cell {(x,y): x<y}
            sage: var_value = complex.find_uncovered_random_point(max_failings=100)
            sage: complex.is_point_covered(var_value)
            False
        """
        num_failings = 0
        while not max_failings or num_failings < max_failings:
            var_value = self.generate_random_var_value(var_bounds=var_bounds)
            # This point is not already covered.
            if self.is_point_covered(var_value):
                num_failings += 1
            else:
                return var_value
        logging.warning("The complex has %s cells. Cannot find one more uncovered point by shooting %s random points" % (len(self.components), max_failings))
        return None

    def find_uncovered_point_mathematica(self, formal_closure=False):
        r"""
        Call Mathematica ``FindInstance`` to get a point that satisfies self.bddbsa
        and is uncovered by any cells in the complex.

        The argument formal_closure whether inequalities are treated as <= 0 or as < 0.
        If such point does not exist, return ``None``.
 
        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.WARN)
            sage: complex = SemialgebraicComplex(lambda x,y: max(x,y), ['x','y'], find_region_type=result_symbolic_expression, default_var_bound=(-10,10)) # optional - mathematica
            sage: complex.add_new_component([1,2], bddbsa=None, flip_ineq_step=0, wall_crossing_method=None, goto_lower_dim=False) # the cell {(x,y): x<y} # optional - mathematica
            sage: var_value = complex.find_uncovered_point_mathematica(formal_closure=True) # optional - mathematica
            sage: var_value # optional - mathematica
            (0, -1)
            sage: complex.is_point_covered(var_value) # optional - mathematica
            False
        """
        bddbsa = BasicSemialgebraicSet_mathematica.from_bsa(self.bddbsa)
        constr = bddbsa.constraints_string()
        if not constr:
            constr = 'True'
        for c in self.components:
            if formal_closure:
                bsa = BasicSemialgebraicSet_mathematica.from_bsa(c.bsa.formal_closure())
            else:
                bsa = BasicSemialgebraicSet_mathematica.from_bsa(c.bsa)
            constr_c = bsa.constraints_string()
            if constr_c:
                constr += ' && !(' + constr_c + ')'
            else:
                # universal bsa covers the whole space.
                return None
        if constr == 'True':
            return tuple(self.bddbsa.ambient_space().zero())
        pt_math = mathematica.FindInstance(constr, bddbsa.variables_string(), 'Reals')
        if len(pt_math) == 0:
            return None
        pt = vector(from_mathematica(pt_math[1][i+1][2]) for i in range(self.d))
        return tuple(bddbsa.ambient_space(field=pt.parent().base_ring())(pt))

    def add_new_component(self, var_value, bddbsa=None, polynomial_map=None, pos_poly=None, flip_ineq_step=0, wall_crossing_method=None, goto_lower_dim=False, num_eq=None):
        r"""
        Compute one proof cell around var_value. Append this cell to the complex.

        - bddbsa defines the boundary constraints satisfied by the cell. By default, bddbsa=None, and complex.bddbsa is taken.
        - If flip_ineq_step = 0, do not search for neighbour testpoints. Used in ``shoot_random_points()``. wall_crossing_method = ``None`` in this case.
        - If flip_ineq_step > 0, search for neighbour testpoints using wall_crossing_method = 'mathematica' or 'heuristic' or 'heuristic_then_mathematica'. Used in bfs. If goto_lower_dim is ``False`` or (goto_lower_dim is ``True`` and wall_crossing method = 'heuristic' but wall is non-linear), then find new testpoints across the wall only.
        - pos_poly is a polynomial. var_value was found by crossing pos_poly < 0 from another cell. The new test points return by new_component.find_neighbour_candidates must satisfy pos_poly(new test point) > 0.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.WARN)
            sage: complex = SemialgebraicComplex(lambda x,y: max(x,y), ['x','y'], find_region_type=result_symbolic_expression, default_var_bound=(-10,10))
            sage: complex.add_new_component([1,2], bddbsa=None, flip_ineq_step=1/10, wall_crossing_method='heuristic', goto_lower_dim=True) # the cell {(x,y): x<y}
            sage: sorted(complex.components[0].bsa.lt_poly())
            [x - y]
            sage: complex.components[0].neighbour_points
            [(31/20, 29/20), (3/2, 3/2)]
            sage: list(complex.points_to_test[0].keys())
            [(31/20, 29/20)]
            sage: list(complex.points_to_test[1].keys())
            [(3/2, 3/2)]

            sage: complex = SemialgebraicComplex(lambda x,y: max(x,y^2), ['x','y'], find_region_type=result_symbolic_expression, default_var_bound=(-10,10))        # optional - mathematica
            sage: complex.add_new_component([1,1/2], bddbsa=None, flip_ineq_step=1/10, wall_crossing_method='mathematica', goto_lower_dim=False) # the cell {(x,y): x > y^2}      # optional - mathematica
            sage: list(complex.components[0].bsa.lt_poly()) # optional - mathematica
            [y^2 - x]
            sage: list(complex.points_to_test[0].keys()) # optional - mathematica
            [(19/20, 1)]
        """
        if bddbsa is None:
            bddbsa = self.bddbsa
        if polynomial_map is None:
            polynomial_map = self.polynomial_map
        if num_eq is None:
            num_eq = len(list(bddbsa.eq_poly()))
        if not num_eq in self.points_to_test:
            self.points_to_test[num_eq] = OrderedDict()
            if not num_eq in self.tested_points:
                self.tested_points[num_eq] = set([])
        new_component = self._cell_class(self.family, var_value, 
                                         find_region_type=self.find_region_type, bddbsa=bddbsa, polynomial_map=polynomial_map)

        new_num_eq = len(list(new_component.bsa.eq_poly()))
        if  new_num_eq > num_eq:
            logging.warning("The cell around %s defined by %s has more equations than boundary %s" %(new_component.var_value, new_component.bsa, bddbsa))
            #import pdb; pdb.set_trace()
            # bsa is lower dimensional as it has more equations than bddbsa, 
            # so we try to perturb the testpoint to obtain a
            # new testpoint in bddbsa that does not fall into a lower dimensional cell.
            # Heuristic code using gradient desecent.  #FIXME.
            for l in (set(new_component.bsa.eq_poly())- set(bddbsa.eq_poly())):
                ineqs = list(new_component.bddbsa.lt_poly())+list(new_component.bddbsa.le_poly())
                pts = [find_point_flip_ineq_heuristic(var_value, l, ineqs, 1/2017), find_point_flip_ineq_heuristic(var_value, -l, ineqs, 1/2017)]
                for pt in pts:
                    if pt is not None:
                        # Find a new point, use polynomial map to recover the values of those eliminated variables.
                        pert_value = tuple(p(pt) for p in polynomial_map)
                        if not (pert_value in self.tested_points[num_eq]) and not (pert_value in self.points_to_test[num_eq]) and (pert_value in bddbsa):
                            self.points_to_test[num_eq][pert_value] = (copy(bddbsa), copy(polynomial_map), pos_poly)
            if not goto_lower_dim:
                return

        if (flip_ineq_step != 0) and (new_component.region_type != 'stop'):
            # when using random shooting, don't generate neighbour points;
            new_points = new_component.find_neighbour_candidates(flip_ineq_step, wall_crossing_method, goto_lower_dim, pos_poly)
            new_component.neighbour_points = []
            for (n, new_points_n_eq) in new_points.items():
                for pt, (new_bddbsa, new_poly_map, no_crossing) in new_points_n_eq.items():
                    if not n in self.points_to_test:
                        self.points_to_test[n] = OrderedDict()
                        if not n in self.tested_points:
                            self.tested_points[n] = set([])
                    if not pt in self.tested_points[n]:
                        if new_component.polynomial_map != new_poly_map:
                            # neighbour cell has more equation. simplify its bddbsa.
                            poly_ring = new_poly_map[0].parent()
                            new_bddbsa = BasicSemialgebraicSet_veronese.from_bsa(BasicSemialgebraicSet_local(new_bddbsa.section(new_poly_map, poly_ring=poly_ring), pt))
                            for i in range(len(new_poly_map)):
                                l = new_poly_map[i]-poly_ring.gens()[i]
                                new_bddbsa.add_polynomial_constraint(l, operator.eq)
                        self.points_to_test[n][pt] = (new_bddbsa, new_poly_map, no_crossing)
                new_component.neighbour_points += list(new_points_n_eq.keys())
        self.components.append(new_component)

    def shoot_random_points(self, num, var_bounds=None, max_failings=1000):
        r"""
        Complete the complex by randomly shooting points.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.WARN)                   # not tested
            sage: complex = SemialgebraicComplex(lambda x,y: max(x,y), ['x','y'], find_region_type=result_symbolic_expression, default_var_bound=(-10,10))      # not tested
            sage: complex.shoot_random_points(100)                # not tested
            sage: complex.plot()                                  # not tested

            sage: complex = SemialgebraicComplex(lambda x,y: max(x,y^2), ['x','y'], find_region_type=result_symbolic_expression, default_var_bound=(-10,10))
            sage: complex.shoot_random_points(100)                # not tested
        """
        for i in range(num):
            var_value = self.find_uncovered_random_point(var_bounds=var_bounds, max_failings=max_failings)
            if var_value is None:
                return
            else:
                self.add_new_component(var_value, bddbsa=self.bddbsa, flip_ineq_step=0, goto_lower_dim=False)

    def plot(self, alpha=0.5, plot_points=300, slice_value=None, goto_lower_dim=False, restart=True, **kwds):
        r"""
        Plot the complex and store the graph.

        - If restart is ``False``, plot the newly added cells on top of the last graph; otherwise, start a new graph.
        - If slice_value is given, it is either a polynomial_map that defines a section, or a list of fixed parameter values with two of them being None. Plot the section. 
        - plot_points controls the quality of the plotting.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.WARN)
            sage: complex = SemialgebraicComplex(lambda x,y,z: min(x^2,y^2,z), ['x','y','z'], find_region_type=result_symbolic_expression, default_var_bound=(-10,10))    # not tested
            sage: complex.bfs_completion(goto_lower_dim=True)             # not tested
            sage: Q.<u,v> = QQ[]                                          # not tested

        Plot the slice in (x,y)-space with z=4 (the following two ways give the same thing)::

            sage: complex.plot(slice_value=[None, None, 4])          # not tested
            sage: complex.plot(slice_value=[u, v, 4])                # not tested

        Plot the slice in (y,z)-space with x=4 (the following two ways give the same thing)::

            sage: complex.plot(slice_value=[4, None, None])          # not tested
            sage: complex.plot(slice_value=[4, u, v])                # not tested

        Plot the slice in (x,y)-space with z=y::

            sage: complex.plot(slice_value=[u, v, v])                # not tested

        Plot the slice in (x,z)-space with y=x/2::

            sage: complex.plot(slice_value=[u, u/2, v])              # not tested
        """
        if restart:
            self.graph = Graphics()
            self.num_plotted_components = 0
        self.graph.set_aspect_ratio(1)
        if 'xmin' in kwds:
            self.graph.xmin(kwds['xmin'])
            xmin = kwds['xmin']
        else:
            xmin = self.default_var_bound[0]   # special treatement in the case goto_lower_dim which uses bsa.plot() instead of component.plot() because zorder is broken in region_plot/ContourPlot.
        if 'xmax' in kwds:
            self.graph.xmax(kwds['xmax'])
            xmax = kwds['xmax']
        else:
            xmax = self.default_var_bound[1]
        if 'ymin' in kwds:
            self.graph.ymin(kwds['ymin'])
            ymin = kwds['ymin']
        else:
            ymin = self.default_var_bound[0]
        if 'ymax' in kwds:
            self.graph.ymax(kwds['ymax'])
            ymax = kwds['ymax']
        else:
            ymax = self.default_var_bound[1]
        # # FIXME: zorder is broken in region_plot/ContourPlot.
        # for c in self.components[self.num_plotted_components::]:
        #     num_eq = len(list(c.bsa.eq_poly()))
        #     gc = c.plot(alpha=alpha, plot_points=plot_points, slice_value=slice_value, default_var_bound=self.default_var_bound, goto_lower_dim=goto_lower_dim, zorder=num_eq, **kwds)  
        #     if gc: # need this because (empty g + empty gc) forgets about xmin xmax ymin ymax.
        #         self.graph += gc
        # Workaround.
        # FIXME: bug example plot(goto_lower_dim=True) in doctest of SemialgebraicComplex
        components = self.components[self.num_plotted_components::]
        for c in components:
            if not list(c.bsa.eq_poly()):
                if not goto_lower_dim:
                    self.graph += c.plot(alpha=alpha, plot_points=plot_points, slice_value=slice_value, default_var_bound=self.default_var_bound, goto_lower_dim=False, zorder=0, **kwds)
                else:
                    color = find_region_color(c.region_type)
                    self.graph += (c.bsa).plot(alpha=alpha, plot_points=plot_points, slice_value=slice_value, color='white', fill_color=color, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
        if goto_lower_dim:
            for c in components:
                if not list(c.bsa.eq_poly()):
                    color = find_region_color(c.region_type)
                    for l in c.bsa.le_poly():
                        new_bsa = copy(c.bsa)
                        new_bsa.add_polynomial_constraint(l, operator.eq)
                        self.graph += new_bsa.plot(alpha=alpha, plot_points=plot_points, slice_value=slice_value, color=color, fill_color=color, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
        for c in components:
            if len(list(c.bsa.eq_poly()))==1:
                self.graph += c.plot(alpha=alpha, plot_points=plot_points, slice_value=slice_value, default_var_bound=self.default_var_bound, goto_lower_dim=False, zorder=0, **kwds)
        if goto_lower_dim:
            for c in components:
                if len(list(c.bsa.eq_poly()))==1:
                    color = find_region_color(c.region_type)
                    for l in c.bsa.lt_poly():
                        new_bsa = BasicSemialgebraicSet_eq_lt_le_sets(eq=list(c.bsa.eq_poly())+[l], lt=[ll for ll in c.bsa.lt_poly() if ll != l], le=list(c.bsa.le_poly()))
                        self.graph += new_bsa.plot(alpha=alpha, plot_points=plot_points, slice_value=slice_value, color='white', fill_color='white', xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
                    for l in c.bsa.le_poly():
                        new_bsa = copy(c.bsa)
                        new_bsa.add_polynomial_constraint(l, operator.eq)
                        self.graph += new_bsa.plot(alpha=alpha, plot_points=plot_points, slice_value=slice_value, color=color, fill_color=color, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
        for c in components:
            if len(list(c.bsa.eq_poly()))==2:
                ptcolor = find_region_color(c.region_type)
                self.graph += point(c.var_value, color = ptcolor, zorder=10)
        self.num_plotted_components = len(self.components)
        return self.graph

    def plot_bfs_tree(self, **kwds):
        r"""
        Plot the bfs tree provided that the complex is 2-dimensional.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.WARN)
            sage: complex = SemialgebraicComplex(drlm_backward_3_slope, ['f','bkpt'])
            sage: complex.bfs_completion(var_value=[4/5,1/2],flip_ineq_step=1/20)    #long time
            sage: g = complex.plot_bfs_tree()
        """
        g = self.plot(**kwds)
        for i in range(len(self.components)):
            c = self.components[i]
            var_value = c.var_value
            for pt in c.neighbour_points:
                g += line([var_value, pt],color='black',linestyle=":")
                g += point([pt], color='black', size=5)
            for j in range(i):
                if tuple(var_value) in self.components[j].neighbour_points:
                    g += line([var_value, self.components[j].var_value],color='red',linestyle="-")
                    break
            g += point([var_value], color='red', size=15, zorder=20)
        return g

    def bfs_completion(self, var_value=None, flip_ineq_step=1/100, check_completion=False, wall_crossing_method='heuristic', goto_lower_dim=False, max_failings=0):
        r"""
        Breadth-first-search to complete the complex.

        - var_value: the starting point. If not given, start with a random point.
        - flip_ineq_step: a small positive number that controls the step length in wall-crossing.
        - check_completion: When check_completion is ``True``, after bfs terminates, check whether the entire parameter space is covered by cells, using Mathematica's ``FindInstance`` (if max_failings=0) or random shooting (if max_failings>0). If an uncovered point has been found, restart the bfs from this point.
        - wall_crossing_method: 'mathematica' or 'heuristic' or 'heuristic_then_mathematica' 
              - wall_crossing_method='heuristic_then_mathematica': try heuristic method first. If it does not find a new testpoint, then use Mathematica.
        - goto_lower_dim: whether lower dimensional cells are considered. If goto_lower_dim is ``False`` or if goto_lower_dim is ``True`` and wall_crossing method is 'heuristic' but the wall is non-linear, then find new testpoint across the wall only.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.WARN)
            sage: complex = SemialgebraicComplex(lambda x,y: min(x^2+ 4*y^2, 4), ['x','y'], find_region_type=result_symbolic_expression, default_var_bound=(-3,3))    # optional - mathematica

        When using mathematica's FindInstance, the returned test point may not be rational. In this example, complex.components[1].var_value is (0, sqrt(799/800)), which is then cast to AA using the function from_mathematica()::

            sage: complex.bfs_completion(check_completion=True, wall_crossing_method='mathematica', goto_lower_dim=True)  # optional - mathematica
            sage: len(complex.components)                               # optional - mathematica
            3
            sage: complex.plot(goto_lower_dim=True)                     # not tested

            sage: complex = SemialgebraicComplex(lambda x,y: min(x+ 4*y, 4), ['x','y'], find_region_type=result_symbolic_expression, default_var_bound=(-5,5))
            sage: complex.bfs_completion(var_value=[1,1])
            sage: len(complex.components)
            2
            sage: complex.plot()                                        # not tested

        See more examples in the docstring of the class :class:`SemialgebraicComplex`.
        """
        if not any(dic for dic in self.points_to_test.values()) and not var_value:
            var_value = self.find_uncovered_random_point()
        num_eq = len(list(self.bddbsa.eq_poly()))
        if var_value:
            if num_eq not in self.points_to_test:
                self.points_to_test[num_eq] = OrderedDict()
                if num_eq not in self.tested_points:
                    self.tested_points[num_eq] = set([])
            if not tuple(var_value) in self.tested_points[num_eq]:
            # put given var_value to num_eq so that it pops out first in bfs.
                self.points_to_test[num_eq][tuple(var_value)] = (self.bddbsa, None, None)
        max_num_eq = max(self.points_to_test.keys())
        while num_eq <= max_num_eq:
            while num_eq in self.points_to_test and self.points_to_test[num_eq]:
                var_value, (bddbsa, polynomial_map, no_crossing) = self.points_to_test[num_eq].popitem(last=False)  # BFS
                self.tested_points[num_eq].add(var_value)
                var_value = list(var_value)
                if not self.is_point_covered(var_value):
                    self.add_new_component(var_value, bddbsa=bddbsa, polynomial_map=polynomial_map, pos_poly=no_crossing, flip_ineq_step=flip_ineq_step, wall_crossing_method=wall_crossing_method, goto_lower_dim=goto_lower_dim, num_eq=num_eq)
            max_num_eq = max(self.points_to_test.keys())
            num_eq += 1
        if check_completion:
            if max_failings == 0:
                #raise NotImplementedError()
                uncovered_pt = self.find_uncovered_point_mathematica(formal_closure=not goto_lower_dim)
            else: # assume that check_completion is an integer.
                uncovered_pt = self.find_uncovered_random_point(max_failings=max_failings)
            if uncovered_pt is not None:
                logging.warning("After bfs, the complex has uncovered point %s." % (uncovered_pt,))
                self.bfs_completion(var_value=uncovered_pt, \
                                    flip_ineq_step=flip_ineq_step, \
                                    check_completion=check_completion, \
                                    wall_crossing_method=wall_crossing_method, \
                                    goto_lower_dim=goto_lower_dim)

    def is_complete(self, formal_closure=False):
        r"""
        Return whether the entire parameter space satisfying self.bddbsa is covered by cells.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.WARN)
            sage: complex = SemialgebraicComplex(lambda x,y: min(x+ 4*y, 4), ['x','y'], find_region_type=result_symbolic_expression, default_var_bound=(-5,5))    # optional - mathematica
            sage: complex.bfs_completion(var_value=[1,1], goto_lower_dim=False)  # optional - mathematica
            sage: complex.is_complete(formal_closure=True)                       # optional - mathematica
            True
            sage: complex.is_complete()   # optional - mathematica
            False
        """
        if self.find_uncovered_point_mathematica(formal_closure=formal_closure) is None:
            return True
        else:
            return False

    def is_polyhedral(self):
        for c in self.components:
            if not c.is_polyhedral():
                return False
        return True

    def subcomplex_of_cells_with_given_region_types(self, given_region_types={'is_extreme'}):
        subcomplex = copy(self)
        subset_components = [c for c in self.components if bool(c.region_type in given_region_types)]
        subcomplex.components = subset_components
        subcomplex.graph = Graphics()
        subcomplex.num_plotted_components = 0
        return subcomplex

    def subcomplexes_of_cells_with_given_region_types(self):
        """
        Return a dictionary mapping region types to subcomplexes.
        """
        region_types = set(c.region_type for c in self.components)
        return { region_type:
                 self.subcomplex_of_cells_with_given_region_types({region_type})
                 for region_type in region_types }

    def guess_boundary(self):
        extlin = [l for c in self.components if not list(c.bsa.eq_poly()) for l in list(c.bsa.lt_poly()) + list(c.bsa.le_poly())]
        boundaries = set(extlin).difference(set([-l for l in extlin]))
        return list(boundaries)

    def is_face_to_face(self):
        sagepolyhedra = []
        for c in self.components:
            p1 = c.sage_polyhedron()
            for p2 in sagepolyhedra:
                p = p1.intersection(p2)
                d = p.dimension()
                if d == -1: # empty
                    continue
                is_a_face = False
                for f in p1.faces(d):
                    if p == f.as_polyhedron():
                        is_a_face = True
                        break
                if not is_a_face:
                    return False
                is_a_face = False
                for f in p2.faces(d):
                    if p == f.as_polyhedron():
                        is_a_face = True
                        break
                if not is_a_face:
                    return False
            sagepolyhedra.append(p1)
        return True

    def polyhedral_complex(self):
        r"""
        Assume that the cells in the :class:`SemialgebraicComplex` are polyhedral and face-to-face.
        Return the :class:`PolyhedralComplex`.
        """
        return PolyhedralComplex([c.sage_polyhedron() for c in self.components])

###########################################
# Helper functions for SemialgebraicComplex
###########################################
def gradient(ineq):
    r"""
    Return the gradient of the polynomial ineq.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: P.<x,y>=QQ[]
        sage: gradient(2*x^2*y+x^3+4*y+5)
        [3*x^2 + 4*x*y, 2*x^2 + 4]
        sage: P.<z>=QQ[]
        sage: gradient(z^3+3*z^2+1)
        [3*z^2 + 6*z]
    """
    if hasattr(ineq, 'gradient'):
       return ineq.gradient()
    else:
       return [ineq.derivative()]

####################################
# Find region type and region color
####################################

def find_region_type_igp(K, h, region_level='extreme', is_minimal=None):
    r"""
    Find the type of a igp function h in the :class:`ParametricRealField` K;
    (is it constructible? is it minimal? is it extreme?)
    Record the comparisons in K.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: K.<f> = ParametricRealField([4/5])
        sage: h = gmic(f, field=K)
        sage: find_region_type_igp(K, h)
        'is_extreme'
        sage: K.get_lt_factor()
        {-f, -f + 1/2, f - 1}

        sage: K.<f,bkpt>=ParametricRealField([1/7,3/7])
        sage: h = drlm_backward_3_slope(f, bkpt, field=K)
        sage: find_region_type_igp(K, h)
        'not_extreme'
        sage: K.get_lt_factor()
        {2*bkpt - 1,
         -f,
         -f + 2*bkpt - 1,
         f - 4*bkpt + 1,
         f - 3*bkpt + 1,
         f - bkpt,
         2*f - 1,
         2*f - bkpt}
    """
    ## Note: region_level = 'constructible' / 'minimal'/ 'extreme'. test cases see find_parameter_region()
    #if hasattr(K, '_big_cells') and K._big_cells:
    #    return find_region_type_igp_extreme_big_cells(K, h)
    if h is None:
        return 'not_constructible'
    if region_level == 'constructible':
        return 'is_constructible'
    if is_minimal is None:
        is_minimal = minimality_test(h, full_certificates=False)
    if is_minimal:
        if region_level == 'minimal':
            return 'is_minimal'
        generator = generate_perturbations(h, full_certificates=False)
        for perturbation in generator:
            return 'not_extreme'
        return 'is_extreme'
    else:
        return 'not_minimal'

def find_region_type_igp_extreme(K, h):
    if find_region_type_igp(K, h) == 'is_extreme':
        return 'is_extreme'
    else:
        return 'stop'

def coarse_regions_from_arrangement_of_bkpts(K, h):
    if h is None:
        return 'not_constructible'
    comparison = lambda a, b: True
    grid_vertices = unique_list(itertools.chain(generate_type_1_vertices(h, comparison), generate_type_2_vertices(h, comparison)))
    for (x, y, z, xeps, yeps, zeps) in grid_vertices:
        d = delta_pi_general(h, x, y, (xeps, yeps, zeps))
    return 'is_constructible'

def claimed_region_type(igp_function=gmic, condition_according_to_literature=True, **kwds):
    r"""
    Return if igp_function is 'not_constructible' or 'constructible' or 'extreme' for the values of parameters given by kwds (use default values if kwds is not provided).

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: claimed_region_type(igp_function=chen_4_slope)
        'extreme'
        sage: claimed_region_type(igp_function=chen_4_slope, condition_according_to_literature=True, f=1/2, s_pos=5, s_neg=-5, lam1=1/5, lam2=1/5)
        'constructible'
        sage: claimed_region_type(igp_function=chen_4_slope, condition_according_to_literature=True, f=7/10, s_pos=2, s_neg=-4, lam1=1/100, lam2=49/100)
        'extreme'
        sage: claimed_region_type(igp_function=chen_4_slope, condition_according_to_literature=False, f=7/10, s_pos=2, s_neg=-4, lam1=1/100, lam2=49/100)
        'constructible'

    The following examples show how claimed_region_type can be used in computing the complex of igp_function. The cells of the complex represent the regions where igp_function is 'extreme', 'constructible' or 'not_constructible' according to the claimed extremality conditions for igp_function::

        sage: complex = SemialgebraicComplex(claimed_region_type, ['f','bkpt'], find_region_type=return_result, igp_function=drlm_backward_3_slope)
        sage: complex.bfs_completion(var_value=[4/5,1/2])
        sage: [c.region_type for c in complex.cells_containing_point([1/12,1/6])]
        ['extreme']

        sage: complex = SemialgebraicComplex(claimed_region_type, ['lam1', 'lam2'], find_region_type=return_result, igp_function=chen_4_slope, condition_according_to_literature=True)
        sage: complex.bfs_completion(var_value=[3/10, 45/101])
        sage: [c.region_type for c in complex.cells_containing_point([1/100,49/100])]
        ['extreme']

        sage: complex = SemialgebraicComplex(claimed_region_type, ['lam1', 'lam2'], find_region_type=return_result, igp_function=chen_4_slope, condition_according_to_literature=False)
        sage: complex.bfs_completion(var_value=[3/10, 45/101])
        sage: [c.region_type for c in complex.cells_containing_point([1/100,49/100])]
        ['constructible']

        sage: complex = SemialgebraicComplex(claimed_region_type, ['lam1', 'lam2'], find_region_type=return_result, igp_function=chen_4_slope, condition_according_to_literature=False, kwds_dict={'f':1/2, 's_pos':5, 's_neg':-5})
        sage: complex.bfs_completion(var_value=[1/5, 1/4])
        sage: [c.region_type for c in complex.cells_containing_point([1/5,1/5])]
        ['extreme']

        sage: complex = SemialgebraicComplex(claimed_region_type,['f','alpha'],find_region_type=return_result, igp_function=dg_2_step_mir)
        sage: complex.bfs_completion(var_value=[4/5,3/10]) #not tested #infinite many cells.
    """
    if not isinstance(igp_function, ParametricFamily_base):
        igp_function = ParametricFamily(igp_function)
    if condition_according_to_literature:
        source = 'literature'
    else:
        source = 'best'
    return igp_function.parameter_attribute(according_to=source, **kwds)

def return_result(field, result):
    return result

def result_concrete_value(field, result):
    r"""
    Return the concrete values in result as a tuple. See also ``result_symbolic_expression()``.
 
    This function can provided to ``find_region_type`` when setting up a :class:`SemialgebraicComplex`. 
    In this way, one can compare result of type :class:`ParametricRealFieldElement` or list of :class:`ParametricRealFieldElement`
    with the previous elements in ``region_type_color_map`` which do not necessarily have the same parent.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.WARN)
        sage: def vol(a,b):
        ....:     P = Polyhedron(ieqs=[(1,0,-1),(0,0,1),(1,-1,0),(0,1,0),(1,-a,-b)])
        ....:     return P.volume()
        sage: K.<a,b>=ParametricRealField([2,3])
        sage: result = vol(a,b)
        sage: result_concrete_value(K, result)
        (1/12,)
    """
    concrete_value = tuple(elt.val() if hasattr(elt, 'val') else elt for elt in flatten([result]))
    return concrete_value

def result_symbolic_expression(field, result):
    r"""
    Return the symbolic expressions in result as a tuple
 
    This function can provided to ``find_region_type`` when setting up a :class:`SemialgebraicComplex`. 
    In this way, one can compare result of type :class:`ParametricRealFieldElement` or list of :class:`ParametricRealFieldElement`
    with the previous elements in ``region_type_color_map`` which do not necessarily have the same parent.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.WARN)
        sage: def vol(a,b):
        ....:     P = Polyhedron(ieqs=[(1,0,-1),(0,0,1),(1,-1,0),(0,1,0),(1,-a,-b)])
        ....:     return P.volume()
        sage: K1.<a,b>=ParametricRealField([2,3])
        sage: vol1 = vol(a,b)
        sage: sym_exp1 = result_symbolic_expression(K1, vol1)
        sage: sym_exp1
        (1/(2*a*b),)
        sage: K2.<a,b>=ParametricRealField([5,4])
        sage: vol2 = vol(a,b)
        sage: sym_exp2 = result_symbolic_expression(K2, vol2)
        sage: sym_exp2
        (1/(2*a*b),)
        sage: vol1 == vol2
        False
        sage: sym_exp1 == sym_exp2
        True    
    """
    symbolic_expression = tuple(elt._sym if hasattr(elt, '_sym') else elt for elt in flatten([result]))
    return symbolic_expression

def find_region_type_igp_extreme_big_cells(K, h):
    """
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.WARN)
        sage: complex = SemialgebraicComplex(drlm_backward_3_slope, ['f','bkpt'], find_region_type=find_region_type_igp_extreme_big_cells)
        sage: complex.bfs_completion(var_value=[4/5,1/2])
        sage: len(complex.components)
        11
        sage: complex.plot() # not tested
    """
    hcopy = copy(h)
    if h is None:
        return 'not_constructible'
    is_extreme = True
    with K.off_the_record():
        h = copy(hcopy)
        for x in h.values_at_end_points():
            if (x < 0) or (x > 1):
                is_extreme  = False
                break
    if not is_extreme:
        assert (x < 0) or (x > 1)
        return False
    f = find_f(h, no_error_if_not_minimal_anyway=True)
    if f is None:
        return False
    bkpt = h.end_points()
    with K.off_the_record():
        h = copy(hcopy)
        for (x, y, z, xeps, yeps, zeps) in generate_nonsubadditive_vertices(h, reduced=True):
            is_extreme = False
            break
    if not is_extreme:
        assert delta_pi_general(h, x, y, (xeps, yeps, zeps)) < 0
        return False
    with K.off_the_record():
        h = copy(hcopy)
        if h(f) != 1:
            is_extreme = False
    if not is_extreme:
        assert h(f) != 1
        return False
    with K.off_the_record():
        h = copy(hcopy)
        for (x, y, xeps, yeps) in generate_nonsymmetric_vertices(h, f):
            is_extreme = False
            break
    if not is_extreme:
        assert h(x)+h(y) != 1
        return False
    # function is minimal.
    with K.off_the_record():
        h = copy(hcopy)
        num_slope = number_of_slopes(h)
    if h.is_continuous() and num_slope == 2:
        # is_extreme  #record #shortcut?
        minimality_test(h, full_certificates=False)
        return True
    ### testing non-extremality with strategical uncovered interval.
    with K.off_the_record():
        h = copy(hcopy)
        ucs = generate_uncovered_components(h)
        f = find_f(h)
        for uncovered_pt in [f/2, (f+1)/2]:
            if any((i[0] == uncovered_pt or i[1] == uncovered_pt) for uc in ucs for i in uc if len(uc)==2):
                uncovered_pts = [uncovered_pt]
                is_extreme = False
                break
        ### not sure, maybe running extremality test gives bigger cells.
        ### the following trick performs good on kzh_3_slope_param_extreme_1, bad on drlm... when f/2 trick is turned off.
        if is_extreme and ucs:
            uc = min(ucs, key=len)
            uncovered_pts = [(i[0]+i[1])/2 for i in uc]
            #uncovered_pts = list(chain(*[[i[0], i[1]] for i in uc]))
            is_extreme = False
    if not is_extreme:
        h = copy(hcopy)
        bkpt1 = h.end_points()
        bkpt2 = h.end_points()[:-1] + [ x+1 for x in h.end_points() ]
        for uncovered_pt in uncovered_pts:
            for y in bkpt1:
                (delta_pi(h, uncovered_pt, y) > 0) or (delta_pi(h, uncovered_pt, y) == 0)
            for z in bkpt2:
                if uncovered_pt < z < 1+uncovered_pt:
                    (delta_pi(h, uncovered_pt, z-uncovered_pt) > 0) or (delta_pi(h, uncovered_pt, z-uncovered_pt) == 0)
            for x in bkpt1:
                (delta_pi(h, x, uncovered_pt-x) > 0) or (delta_pi(h, x, uncovered_pt-x) == 0)
        return False
    # maybe is_extreme
    h = copy(hcopy)
    is_minimal = minimality_test(h, full_certificates=False)
    assert is_minimal
    generator = generate_perturbations(h, full_certificates=False)
    for perb in generator:
        logging.warning("Function is not extreme but f/2 and (1+f)/2 are covered.")
        return False
    return True

region_type_color_map = [('not_constructible', 'lightgrey'), ('is_constructible', 'black'), ('not_minimal', 'orange'), ('is_minimal', 'darkgrey'),('not_extreme', 'green'), ('is_extreme', 'blue'), ('stop', 'grey'), (True, 'blue'), (False, 'red'), ('constructible', 'darkgrey'), ('extreme', 'red')]

def find_region_color(region_type):
    r"""
    Return the color of the region according to the global dictionary ``region_type_color_map``.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: find_region_color('is_extreme')
        'blue'
        sage: find_region_color(False)
        'red'
    """
    # if region_type is a color, return directly.
    if region_type in colors:
        return region_type
    # Assume that there will be at most 7 different types other than the ones like 'is_extreme' that were in region_type_color_map.
    global region_type_color_map
    n = len(region_type_color_map)
    i = 0
    while (i < n) and (region_type != region_type_color_map[i][0]):
        i += 1
    if i == n:
        region_color = color_of_ith_region_type(n - 9)  # the initial map includes 9 igp region types.
        region_type_color_map.append((region_type, region_color))
    else:
        region_color = region_type_color_map[i][1]
    return region_color

def color_of_ith_region_type(i):
    r"""
    Return a color in the rainbow.
    """
    j = (4 * i) % 7
    c = rainbow(7)[j]
    return c


##########################
# wall crossing heuristic
##########################

def find_point_flip_ineq_heuristic(current_var_value, ineq, ineqs, flip_ineq_step):
    r"""
    The current_var_value satisfies that l(current_var_value) <= 0 for l=ineq and for every l in ineqs, 
    where ineq is a polynomial and ineqs is a list of polynomials.

    Use heuristic method (gradient descent method with given small positive step length flip_ineq_step)
    to find a new_point (type is tuple) such that 
    ineq(new_point) > 0, l(new_point) < 0 for all l in ineqs
    Return new_point, or ``None`` if it fails to find one.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: P.<a,b>=QQ[]
        sage: find_point_flip_ineq_heuristic([1,1/2], a+b-2, [-a+b^2], 1/4)
        (11/8, 7/8)

    After walking towards ineq==0, ineqs < 0 is violated. Need to adjust::

        sage: x, y = find_point_flip_ineq_heuristic([1,9/10], a+b-2, [-a+b^2], 1/2)  # got (123/85, 179/170)? (253/170, 86/85)?
        sage: (x + y - 2 > 0) and (-x + y^2 < 0)
        True
        sage: x, y = find_point_flip_ineq_heuristic([11/40,1/2], a+b-2, [-a+b^2], 1/4) # got (39295901/31739294, 125037049/123564610), now get (45793/36987, 34307/33903)? (27439/22240, 22601/22240)?
        sage: (x + y - 2 > 0) and (-x + y^2 < 0)
        True

    Ineq is a redundant inequality, it's impossible to cross it without crossing any other ineqs.
    Thus, return ``None``::

        sage: pt = find_point_flip_ineq_heuristic([1,1/2], a+b-2, [-a+b^2, a-1], 1/4); pt is None
        True

    Bug example in cpl Cell 9 with test point (499/1250, 488072439572/4866126017667). Output was (333/1000, 167/1500) which violates ineqs[-1](pt) < 0. BUG! Notice that using sorted(ineqs) with the same current_var_value input gave a correct output. See def adjust_pt_to_satisfy_ineqs below::

        sage: P.<f,z>=QQ[]; ineq = -3*f + 1; flip_ineq_step = 1/1000; current_var_value = (499/1250, 488072439572/4866126017667)
        sage: ineqs = [2*f + 2*z - 1, f + 5*z - 1, -f - 6*z + 1, -2*f - 3*z + 1]
        sage: pt = find_point_flip_ineq_heuristic(current_var_value, ineq, ineqs, flip_ineq_step) #got pt (333/1000, 67/600)
        sage: all(l(pt) < 0 for l in ineqs) and ineq(pt)>0
        True

    Bug example from cpl cell 14 with test point (553/443147, 664/882853), theta 15 = ((-2*z)/(f - 1), 0), fixed. new bug: Bigcell r14, theta29, different order in ineqs::

        sage: P.<f,z>=QQ[]; current_var_value = (553/443147, 664/882853); ineq = f^2 + f*z - f + z; ineqs = [2*f + 6*z - 1, f - 2*z, -f + z]; flip_ineq_step = 1/1000
        sage: pt = find_point_flip_ineq_heuristic(current_var_value, ineq, ineqs, flip_ineq_step)  #(7927/71253, 39430/357637)? (17182/149045, 34610/302851)?
        sage: all(l(pt) < 0 for l in ineqs) and ineq(pt)>0
        True
        sage: current_var_value = (553/443147, 664/882853); ineq = f^2 + f*z - f + z; ineqs = [f - 2*z, -f + z, 2*f + 6*z - 1]
        sage: pt = find_point_flip_ineq_heuristic(current_var_value, ineq, ineqs, flip_ineq_step) # got (17182/149045, 34610/302851)
        sage: all(l(pt) < 0 for l in ineqs) and ineq(pt)>0
        True

    Bug example from cpl cell 17 with test point (1007/8000, 999/8000), theta 16 ((-2*z)/(f - 1), 0)::

        sage: P.<f,z>=QQ[]; current_var_value = (1007/8000, 999/8000); ineq = -f^2 - f*z + f - z; ineqs = [2*f + 5*z - 1, f - 2*z, -2*f - 6*z + 1, -f + z]; flip_ineq_step = 1/1000
        sage: pt = find_point_flip_ineq_heuristic(current_var_value, ineq, ineqs, flip_ineq_step) # got pt (13117/83401, 22445/196184)? (6824/43185, 11717/102636)
        sage: all(l(pt) < 0 for l in ineqs) and ineq(pt)>0
        True

    Bug example in cpl cell 15 with test point (664/882853, 553/443147) theta 6 (f/(f - 4*z + 1), 0). Return None::

        sage: P.<f,z>=QQ[]; current_var_value = (8608/76197, 11223/98473); ineq = -f^2 + 3*f*z + 4*z^2 - z; ineqs = [2*f + 6*z - 1, f - z, -2*f*z + 8*z^2 + f - 2*z, f*z - 4*z^2 - f + z]; flip_ineq_step = 1/1000
        sage: pt = find_point_flip_ineq_heuristic(current_var_value, ineq, ineqs, flip_ineq_step)
        sage: pt is None
        True

    Bug example in cpl bigcell r27 with test point (8242/139557, 18896/107013) theta9 ((-z)/(f + 2*z - 1), 0)::

        sage: P.<f,z>=QQ[]; current_var_value = (8242/139557, 18896/107013); ineq = -f^2 - 3*f*z - 2*z^2 + f; ineqs = [-f, f + 5*z - 1, f - z, -2*f - 5*z + 1]; flip_ineq_step = 1/1000
        sage: pt = find_point_flip_ineq_heuristic(current_var_value, ineq, ineqs, flip_ineq_step) # got pt (8269/72753, 16950/109541)
        sage: all(l(pt) < 0 for l in ineqs) and ineq(pt)>0
        True

    Bug example::

        sage: P.<f,z>=QQ[]; current_var_value = (111/1000, 111/1000); ineq = -11*f^2 + f; ineqs = [9*f- 1, 5*f^2 - f, -11*f + 1]; flip_ineq_step=1/1000 # Can't find a pt whose f < 1/11 because -11*f^2 + f and -11*f + 1 are contradiction.
        sage: pt = find_point_flip_ineq_heuristic(current_var_value, ineq, ineqs, flip_ineq_step); pt is None
        True

    Bug example, need (ineq_value <= 1e-10) rather than <=0. Without that the return was new_point = (7001/26000, 3001/26000) which was wrong as ineq(new_point) == 0::

        sage: P.<bkpt,f>=QQ[]; current_var_value = [7001/26000, 3001/26000]; ineq = -13*bkpt + 13*f + 2; ineqs = [26*bkpt - 11, -169*bkpt*f + 169*f^2 + 3, -26*bkpt + 7]; flip_ineq_step = 1/2017;
        sage: pt = find_point_flip_ineq_heuristic(current_var_value, ineq, ineqs, flip_ineq_step) # got pt (3738/13883, 7795/67523)
        sage: all(l(pt) < 0 for l in ineqs) and ineq(pt)>0
        True
    """
    # heuristic method.
    ineq_gradient = gradient(ineq)
    current_point = vector(RR(x) for x in current_var_value) # Real numbers, faster than QQ
    ineq_value = ineq(*current_point)
    try_before_fail =  min(ceil(2/flip_ineq_step), 2000)  # define maximum number of walks.
    while (ineq_value <= 1e-10) and (try_before_fail > 0):
        ineq_direction = vector(g(*current_point) for g in ineq_gradient)
        if ineq.degree() == 1:
            step_length = (-ineq(*current_point)+flip_ineq_step) / (ineq_direction * ineq_direction)
        else:
            step_length = flip_ineq_step / (ineq_direction * ineq_direction) # ineq_value increases by flip_ineq_step=0.001 roughly
            if step_length > 1:
                step_length = 1  # ensure that distance of move <= sqrt(flip_ineq_step) = 0.1 in each step
        current_point += step_length * ineq_direction
        new_ineq_value = ineq(*current_point)
        if new_ineq_value <= ineq_value:
            #print (False)
            return None
        ineq_value = new_ineq_value
        try_before_fail -= 1
        #print (current_point, RR(ineq_value))
    if ineq_value <= 0:
        return None
    new_point = adjust_pt_to_satisfy_ineqs(current_point, ineq, ineqs, [], flip_ineq_step)
    return new_point #type is tuple

def adjust_pt_to_satisfy_ineqs(current_point, ineq, strict_ineqs, nonstrict_ineqs, flip_ineq_step):
    r"""
    Walk from current_point (type=vector) in the direction perpendicular to 
    the gradient of ineq with small positive step length flip_ineq_step, 
    while maintaining the value of ineq(*current_point) which is >= 0.
    until get a new point such that l(new point)<0 for all l in strict_ineqs and l(new point)<=0 for all l in nonstrict_ineqs

    Return new_point, or ``None`` if it fails to find one.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: P.<a,b>=QQ[]
        sage: x, y = adjust_pt_to_satisfy_ineqs(vector([13/10,12/10]), a+b-2, [-a+b^2], [], 1/2) # was (123/85, 179/170), now (253/170, 86/85)
        sage: (x + y - 2 > 0) and (-x + y^2 < 0)
        True
        sage: x, y = adjust_pt_to_satisfy_ineqs(vector([71/80, 89/80]), a+b-2, [-a+b^2], [], 1/4) # was (171073319/163479120, 155884921/163479120), (22159/20640, 19121/20640)  # ineq(71/80, 89/80) == 0. now (446132279/359982240, 363827761/359982240).
        sage: (x + y - 2 >= 0) and (-x + y^2 < 0)
        True

    If impossible, return ``None``::

        sage: pt = adjust_pt_to_satisfy_ineqs(vector([11/8, 7/8]), a+b-2, [-a+b^2], [a-1], 1/4); pt is None
        True

    Bug example in cpl Cell 9 with test point (499/1250, 488072439572/4866126017667). Without converting input to QQ, output was (0.333000000000000, 0.111333333333333) with  -2*f - 3*z + 1 = -5.55111512312578e-17 , the QQ of the output=(333/1000, 167/1500) has -2*f - 3*z + 1 == 0. Revise the code to take QQ input current point:: 

        sage: P.<f,z>=QQ[]
        sage: current_point = vector([0.333000000000000, 0.100300000000000]); ineq = -3*f + 1; strict_ineqs = [2*f + 2*z - 1, f + 5*z - 1, -f - 6*z + 1, -2*f - 3*z + 1]; nonstrict_ineqs = []; flip_ineq_step = 1/1000
        sage: pt = adjust_pt_to_satisfy_ineqs(current_point, ineq, strict_ineqs, nonstrict_ineqs, flip_ineq_step) # got (333/1000, 67/600)
        sage: all(l(pt) < 0 for l in strict_ineqs) and all(l(pt) <= 0 for l in nonstrict_ineqs) and ineq(pt)>0
        True

    Bug example in cpl cell 3 and 25th theta=((f+5z-1)/(2f+12z-2), z/(2f+12z-2)). Computing with QQ results in huge denominator, never ends:
        sage: P.<f,z>=QQ[]
        sage: current_point = vector([RR(30136191997/49655508552), RR(3903863311/49655508552)])
        sage: strict_ineqs = [f - 1, -2*f + 1, f^2 - 2*f + 1] #used polynomial sub instead of section, had more ineqs: [1/5*f - 1/5, -2*f + 1, 8/25*f^2 - 6/25*f - 2/25, 2/25*f^2 - 4/25*f + 2/25] # empty because (f-1)^2<0.
        sage: nonstrict_ineqs = []
        sage: pt = adjust_pt_to_satisfy_ineqs(current_point, None, strict_ineqs, nonstrict_ineqs, flip_ineq_step=1/1000); pt is None
        True

    Bug example in cpl bigcell 16 with test point (12219/26000, 24/1625). Redo with QQ had infinite loop. Bug comes from find_neighbour_point where it calls bsa_section.upstairs()._polyhedron.is_empty(), which is not strong enough. If we could test bsa_section is empty (perhaps by tighten_upstairs_by_mccormick), then this example should not appear.
        sage: P.<f,z>=QQ[]; 
        sage: current_point = vector((71582788/143165577, 4673/377000)) # came from vector((RR(70727/150800), RR(4673/377000))), 
        sage: ineq=None; strict_ineqs=[2*f - 1, -9*f + 2]; nonstrict_ineqs=[4*f^2 - 4*f + 1]; flip_ineq_step=1/1000
        sage: pt = adjust_pt_to_satisfy_ineqs(current_point, None, strict_ineqs, nonstrict_ineqs, flip_ineq_step=1/1000); pt is None #long time
        True
    """
    #current_point is a vector
    if ineq is not None:
        ineq_gradient = gradient(ineq)
    if all(x.parent()==QQ for x in current_point):
        max_walks = min(ceil(2/flip_ineq_step), 20)
    else:
        max_walks = min(ceil(2/flip_ineq_step), 200) #1000? # define maximum number of walks.
    ineqs_and_strictness = [(l, True) for l in strict_ineqs] + [(l, False) for l in nonstrict_ineqs]
    for i in range(len(ineqs_and_strictness)):
        (l, strictness) = ineqs_and_strictness[i]
        l_gradient = gradient(l)
        l_value = l(*current_point)
        try_before_fail = max_walks
        while ((strictness and (l_value > -1e-10)) or (not strictness and (l_value > 0))) and (try_before_fail > 0):
            l_direction = vector(-g(*current_point) for g in l_gradient) #decrease l_value
            if ineq is not None: # and ineq(*current_point) < 2 * flip_ineq_step??
                ineq_direction = vector(g(*current_point) for g in ineq_gradient)
                if ineq(*current_point) > flip_ineq_step / 2: # want that ineq_value remains the same
                    s = (ineq_direction * l_direction) / (ineq_direction * ineq_direction)
                    projected_direction = l_direction - s * ineq_direction
                    step_length = (l_value+flip_ineq_step) / (projected_direction * l_direction)
                else: # want that l_value remains the same, ineq_value increases by roughly flip_ineq_step
                    s = (ineq_direction * l_direction) / (l_direction * l_direction)
                    projected_direction = ineq_direction - s * l_direction
                    step_length = flip_ineq_step / (projected_direction * ineq_direction)
            else:
                projected_direction = l_direction
                step_length = (l_value+flip_ineq_step) / (projected_direction * l_direction)
            if projected_direction == 0:
                #print ('projected_direction == 0')
                return None
            if step_length == Infinity:  # projected_direction almost 0
                #print (step_length, 'step_length == Infinity')
                return None
            for j in range(len(ineqs_and_strictness)):
                ll = ineqs_and_strictness[j][0]
                if (j == i) or ll(*current_point) > 0:
                    continue
                change_in_ll_value = vector(g(*current_point) for g in gradient(ll)) * projected_direction
                if change_in_ll_value > 0 and -ll(*current_point)*9/10/change_in_ll_value < step_length:
                    step_length = -ll(*current_point)*9/10/change_in_ll_value
                while (try_before_fail > 0) and (ll(*(current_point + step_length * projected_direction)) > 0):
                    step_length = step_length / 2
                    try_before_fail -= 1
            while (ineq is not None) and (try_before_fail > 0) and ineq(*(current_point + step_length * projected_direction)) < 0:
                step_length = step_length / 2
                try_before_fail -= 1
            if abs(step_length) < 1e-15:
                #print (step_length,  'step_length == 0')
                return None
            current_point += step_length * projected_direction
            l_value = l(*current_point)
            try_before_fail -= 1
            #print (try_before_fail, vector([RR(x) for x in current_point]), RR(step_length), RR(l_value))
        if (try_before_fail < 0) or (strictness and (l_value >= 0)) or (not strictness and (l_value > 0)):
            #print ("bad")
            return None
    for l in strict_ineqs:
        if l(*current_point) >= 0:
            return None
    for l in nonstrict_ineqs:
        if l(*current_point) > 0:
            return None
    if ineq is not None and ineq(*current_point) < 0:
        return None
    if all(x.parent()==QQ for x in current_point):
        return tuple(current_point)
    else:
        current_point = vector(QQ(x.n(30)) for x in current_point)
        if ineq is not None and ineq(*current_point) < 0:
            return None
        if all(l(*current_point) < 0 for l in strict_ineqs) and all(l(*current_point) <= 0 for l in nonstrict_ineqs):
            return tuple(current_point)
        # Redo it with QQ input.
        return adjust_pt_to_satisfy_ineqs(current_point, ineq, strict_ineqs, nonstrict_ineqs, flip_ineq_step)

################################################
#  Is the given function contained in a family?
################################################

def embed_function_into_family(given_function, parametric_family, check_completion=False, **opt):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.WARN)
        sage: embed_function_into_family(gmic(3/13), mlr_cpl3_a_2_slope)
        {'r0': 3/13, 'z1': 3/26}
        sage: given_function = drlm_backward_3_slope(1/12-1/100, 1/6+1/100)
        sage: embed_function_into_family(given_function, drlm_backward_3_slope)
        {'bkpt': 53/300, 'f': 11/150}
        sage: embed_function_into_family(gmic(4/13), mlr_cpl3_a_2_slope)
        {'r0': 4/13, 'z1': 3/26}
        sage: embed_function_into_family(gmic(4/13), drlm_backward_3_slope)
        {}
        sage: embed_function_into_family(gmic(4/13), drlm_backward_3_slope, check_completion=True)  # optional - mathematica
        {}

        sage: given_function = mlr_cpl3_b_3_slope(r0=3/26, z1=1/13)
        sage: parametric_family = drlm_backward_3_slope
        sage: param_values = embed_function_into_family(given_function, parametric_family)
        sage: param_values
        {'bkpt': 7/26, 'f': 3/26}
        sage: given_function == parametric_family(**param_values)
        True

        sage: given_function = automorphism(cpl3_8(f=3/5, z1=1/25, z2=1/125))
        sage: parametric_family = gj_forward_3_slope
        sage: embed_function_into_family(given_function, parametric_family, check_completion=True) # optional - mathematica, long time - 240s
        {'f': 2/5, 'lambda_1': 6/25, 'lambda_2': 2/75}
        # sage: embed_function_into_family(given_function, parametric_family, wall_crossing_method='mathematica') # not tested # optional - mathematica #long time  #TODO: does not terminate.
        # {'f': 2/5, 'lambda_1': 6/25, 'lambda_2': 2/75}
    """
    if 'var_name' in opt and 'var_value' in opt:
        var_name = opt.pop('var_name')
        var_value = opt.pop('var_value')
    else:
        default_args = read_default_args(parametric_family)
        var_name = []
        var_value = []
        for (name, value) in default_args.items():
            if not isinstance(value, bool) and not value is None:
                try:
                    RR(value)
                    var_name.append(name)
                    var_value.append(value)
                except:
                    pass
    def frt(K, h):
        if h is None:
            return False
        else:
            return h == given_function
    complex = SemialgebraicComplex(parametric_family, var_name, find_region_type=frt)
    # terminate complex.bfs_completion(var_value=var_value, goto_lower_dim=True) immediately when a cell has region_type == True.
    init_num_eq = len(list(complex.bddbsa.eq_poly()))
    complex.points_to_test[init_num_eq] = OrderedDict()
    complex.tested_points[init_num_eq] = set([])
    complex.points_to_test[init_num_eq][tuple(var_value)] = (None, None, None)
    flip_ineq_step = opt.pop('flip_ineq_step', 1/1000)
    goto_lower_dim = opt.pop('goto_lower_dim', True)
    is_complete = False
    while not is_complete:
        num_eq = init_num_eq
        max_num_eq = max(complex.points_to_test.keys())
        while num_eq <= max_num_eq:
            while num_eq in complex.points_to_test and complex.points_to_test[num_eq]:
                var_value, (bddbsa, polynomial_map, no_crossing) = complex.points_to_test[num_eq].popitem(last=False)
                complex.tested_points[num_eq].add(var_value)
                var_value = list(var_value)
                if not complex.is_point_covered(var_value):
                    complex.add_new_component(var_value, bddbsa=bddbsa, polynomial_map=polynomial_map, pos_poly=no_crossing, flip_ineq_step=flip_ineq_step, goto_lower_dim=goto_lower_dim, num_eq=num_eq, **opt)
                    if complex.components and complex.components[-1].region_type is True:
                        dic = {var_name[i]:var_value[i] for i in range(len(var_name))}
                        #import pdb; pdb.set_trace()
                        return dic
            max_num_eq = max(complex.points_to_test.keys())
            num_eq += 1
        if check_completion:
            var_value = complex.find_uncovered_point_mathematica(formal_closure=False)
            if not var_value:
                is_complete = True
            else:
                if not (tuple(var_value) in complex.tested_points[init_num_eq]):
                    complex.points_to_test[init_num_eq][tuple(var_value)] = (None, None, None)
        else:
            is_complete = True
    # plot_cpl_components(complex.components)
    return {}

"""
EXAMPLES::

    sage: from cutgeneratingfunctionology.igp import *
    sage: K.<f> = ParametricRealField([4/5])
    sage: h = gmic(f, field=K)
    sage: _ = extremality_test(h)
    sage: sorted(K.get_eq_factor()), sorted(K.get_lt_factor()), sorted(K.get_le_factor())
    ([], [-f, -f + 1/2, f - 1], [])
    sage: sorted(K._bsa.eq_poly()), sorted(K._bsa.lt_poly()), sorted(K._bsa.le_poly())
    ([], [-2*f + 1, f - 1], [])

    sage: K.<f> = ParametricRealField([1/5])
    sage: h = drlm_3_slope_limit(f, conditioncheck=False)
    sage: _ = extremality_test(h)
    sage: sorted(K._bsa.eq_poly()), sorted(K._bsa.lt_poly()), sorted(K._bsa.le_poly())
    ([], [-f, 3*f - 1], [])

    sage: K.<f, lam> = ParametricRealField([4/5, 1/6])
    sage: h = gj_2_slope(f, lam, field=K, conditioncheck=False)
    sage: sorted(K._bsa.lt_poly())
    [-lam, -f, f - 1, -f*lam - f + lam]
    sage: for lhs in [-lam, lam - 1, 2*lam - 1, -f, f - 1, -3*f*lam - f + 3*lam, -f*lam - 3*f + lam + 2, -f*lam - f + lam, f*lam - 3*f - lam + 2, f*lam - f - lam, 3*f*lam - f - 3*lam, 3*f*lam + f - 3*lam - 2]: assert lhs < 0
    sage: sorted(K._bsa.lt_poly())
    [-lam,
     2*lam - 1,
     f - 1,
     -3*f*lam - f + 3*lam,
     -f*lam - 3*f + lam + 2,
     f*lam - 3*f - lam + 2,
     3*f*lam - f - 3*lam]

    sage: K.<f,alpha> = ParametricRealField([4/5, 3/10])             # Bad example! parameter region = {given point}.
    sage: h=dg_2_step_mir(f, alpha, field=K, conditioncheck=False)
    sage: _ = extremality_test(h)
    sage: sorted(K._bsa.eq_poly()), sorted(K._bsa.lt_poly()), sorted(K._bsa.le_poly())
    ([10*alpha - 3, 5*f - 4], [], [])
"""
