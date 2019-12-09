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

    # def add_initial_space_dim(self):
    #     if self._bsa.polynomial_map():
    #         # the ParametricRealField already has monomials recorded. Not brand-new.
    #         return
    #     n = len(self._names)
    #     P = PolynomialRing(QQ, self._names)
    #     monomial_list = list(P.gens())
    #     v_dict = {P.gens()[i]:i for i in range(n)}
    #     polyhedron = BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(ambient_dim=n)
    #     self._bsa = BasicSemialgebraicSet_veronese(polyhedron, ambient_dim=n,
    #                                                polynomial_map=monomial_list, v_dict=v_dict)

    def plot(self, *options, **kwds):
        return self.make_proof_cell().plot(*options, **kwds)


###############################
# TO REFACTOR using section. DELETE
###############################
def find_polynomial_map(eqs=[], poly_ring=None):
    """
    HAHAHA, GOODBYE find_polynomial_map
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

# def substitute_llt_lle(llt, lle, var_map, var_name, var_value):
#     r"""
#     Return a list of strict inequalities and a list of non-strict inequalities
#     after substitution using var_map and simplification using the Reformulation-linearization trick.

#     Used in ``SemialgebraicComplexComponent.__init__``

#     EXAMPLES::

#         sage: from cutgeneratingfunctionology.igp import *
#         sage: logging.disable(logging.INFO)
#         sage: P.<x,y>=QQ[]
#         sage: var_map = {y: y, x: 75/19*y}
#         sage: llt = [21*x - 8, -x, 950*x^2 - 3700*x*y - 225*y^2 - 133*x]
#         sage: llt, lle = substitute_llt_lle(llt, [], var_map, ['x','y'], [38/100, 96/1000])
#         sage: sorted(llt)
#         [-y, 1575*y - 152]
#     """
#     # Shall we factorize ineq in lins after substitution using var_map?
#     # Yes, otherwise, got -525/19*lam2^2 - 525*lam2 in chen's 4 slope
#     # at (lam1, lam2) = [20015786415699825/52611587391975857, 5070665891977289/52611587391975857]
#     K = ParametricRealField(var_value, var_name)
#     K.add_initial_space_dim() # needed?
#     for l in llt:
#         ineq = l.parent()(l.subs(var_map))
#         assert ineq(K.gens())< 0 #always True
#     for l in lle:
#         ineq = l.parent()(l.subs(var_map))
#         assert ineq(K.gens())<= 0 #always True
#     bsa = BasicSemialgebraicSet_eq_lt_le_sets.from_bsa(K._bsa)
#     return list(bsa.lt_poly()), list(bsa.le_poly())

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

        sage: complex = SemialgebraicComplex(foo, ['x','y'], find_region_type=lambda r:r, default_var_bound=(-5,5))
        sage: K.<x,y> = ParametricRealField([1,1/2])
        sage: region_type = foo(*K.gens())
        sage: component = SemialgebraicComplexComponent(K, region_type)
        sage: list(component.bsa.lt_poly())
        [x + y - 2, y^2 - x]
        sage: component.plot()                                  # not tested
        sage: new_points = component.find_neighbour_candidates(1/4, 'heuristic', goto_lower_dim=False)
        sage: new_pts = list(new_points[0].keys()) # was [(11/8, 7/8), (19959383/28510088, 24590405/28510088)], now with .n(30) get [(11/8, 7/8), (30093/42985, 54831/63571)]
        sage: len(new_pts)
        2
        sage: all((xx + yy - 2)*(yy^2 - xx) < 0 for (xx, yy) in new_pts)
        True


    # component.find_walls_and_new_points(1/4, 'mathematica', goto_lower_dim=True)  # optional - mathematica
    # ([x + y - 2, y^2 - x],
    #  {(0, 0): [y^2 - x],
    #   (2, 0): [x + y - 2],
    #   (17/8, 0): [],
    #   (2065/512, -33/16): []})

    # FIXME: old code gave ([x - 2*y], [-y, 3*y - 2], []) by simplification through K.
    Test variable elimination::

        #sage: complex = SemialgebraicComplex(foo, ['x','y'], find_region_type=lambda r:r, default_var_bound=(-5,5))
        sage: K.<x,y> = ParametricRealField([1,1/2])
        sage: region_type = foo(*K.gens())
        sage: x == 2*y
        True
        sage: sorted(K._bsa.eq_poly()), sorted(K._bsa.lt_poly()), sorted(K._bsa.le_poly())
        ([x - 2*y], [3*x - 4, y^2 - x], [])
        sage: l = list(K._bsa.eq_poly())[0]; polynomial_map = [2*l.parent().gens()[1], l.parent().gens()[1]]
        sage: s1 = K._bsa.section(polynomial_map)
        sage: sorted(s1.eq_poly()), sorted(s1.lt_poly()), sorted(s1.le_poly())
        ([0], [6*y - 4, y^2 - 2*y], [])
        sage: s2 = K._bsa.section(polynomial_map, bsa_class='veronese')
        sage: sorted(s2.eq_poly()), sorted(s2.lt_poly()), sorted(s2.le_poly())
        ([], [3*y - 2, y^2 - 2*y], [])
        sage: s2.ambient_dim()  # because poly_ring has dim 2.
        2
        sage: x^2 == 8*y^3
        True
        sage: s3 = K._bsa.section(polynomial_map, bsa_class='veronese')
        sage: sorted(s3.eq_poly()), sorted(s3.lt_poly()), sorted(s3.le_poly())
        ([2*y^3 - y^2], [3*y - 2, y^3 - y], [])

        sage: sorted(K._bsa.eq_poly()), sorted(K._bsa.lt_poly()), sorted(K._bsa.le_poly())
        ([x - 2*y, 8*y^3 - x^2], [3*x - 4, y^2 - x], [])
        sage: component = SemialgebraicComplexComponent(K, region_type, polynomial_map=polynomial_map)
        sage: sorted(component.bsa.eq_poly()), sorted(component.bsa.lt_poly()), sorted(component.bsa.le_poly())
        ([-x + 2*y, 2*y^3 - y^2], [3*y - 2, y^3 - y], [])  
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
        self.polynomial_map = polynomial_map
        # Take input polynomial_map so that the stupid function find_polynomial_map is not needed any more.
        if bddbsa is None:   # for doctest purpose
            self.bddbsa = BasicSemialgebraicSet_veronese(poly_ring=poly_ring)
        else:
            # WHY is this input polynomial_map sometimes not compatible with the variable elimination done in bddbsa? Because upstairs ppl bsa eliminates large x_i in the inequalities, and x_i doesn't necessarily correspond to the i-th variable in poly_ring. Since polynomial_map and v_dict were not given at the initialization of veronese, the variable first encounted in the constraints is considered as x0 by upstairs ppl bsa.
            # In old code, we fixed the order of upstairs variables by adding initial space dimensions. We don't do that in the current code. Instead, we take the section of bddbsa to eliminate the varibles in the equations.
            # Is the given bddbsa required to be veronese with upstairs being ppl_bsa? Convert it anyway. # It's the same as BasicSemialgebraicSet_veronese.from_bsa(bddbsa.section(self.polynomial_map), poly_ring=poly_ring) # Taking section forgets the equations
            self.bddbsa = bddbsa.section(self.polynomial_map, bsa_class='veronese', poly_ring=poly_ring)
            for i in range(len(self.var_name)):
                if polynomial_map[i] != poly_ring.gens()[i]:
                    (self.bddbsa).add_polynomial_constraint(polynomial_map[i]-poly_ring.gens()[i], operator.eq)                   
        assert (K.gens() in self.bddbsa) # record boundary constraints. #NEW: moved from ProofCell._construct_field_and_test_point
        #need to apply self.polynomial_map to record bdd constraints with correct variables!
        #assert all(f(self.polynomial_map)(K.gens()) == 0 for f in self.bddbsa.eq_poly()) and all(f(self.polynomial_map)(K.gens()) <= 0 for f in self.bddbsa.le_poly()) and all(f(self.polynomial_map)(K.gens()) < 0 for f in self.bddbsa.lt_poly())

        # Use veronese to simplify the inequalities.
        # The equations are expected to cancel in veronese (otherwise the cell is lower dimensional).
        self.bsa = K._bsa.section(polynomial_map, bsa_class='veronese', poly_ring=poly_ring)  # this is a bigger_bsa
        # Then add back the equations from K._bsa
        # Finally self.bsa should be the same as K._bsa, but its inequalities don't have variables eliminated by polynomial map, so that heuristic wall crossing can be done later.
        #The following are the equations recorded by the old code.
        #for l in K._bsa.eq_poly():
        #    (self.bsa).add_polynomial_constraint(l, operator.eq)
        #Should be the same as adding poly_ring.gens()[i] == polynomial_map[i] for all i.
        for i in range(len(self.var_name)):
            if polynomial_map[i] != poly_ring.gens()[i]:
                (self.bsa).add_polynomial_constraint(polynomial_map[i]-poly_ring.gens()[i], operator.eq)
        #Was BasicSemialgebraicSet_eq_lt_le_sets class was used. self.bsa = BasicSemialgebraicSet_eq_lt_le_sets(poly_ring=poly_ring, eq=K._bsa.eq_poly(), lt=self.bsa.lt_poly(), le=bigger_bsa.le_poly())
        # FIXME: (ignore, as we don't need find_polynomial_map any more) Haven't the upstairs ppl bsa of K._bsa already eliminated the variables with larger indices in the inequaliteis? # Shall we make find_polynomial_map compatible with this echelon form by using for i in range(m, -1, -1), for j in...? 
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

    def plot(self, alpha=0.5, plot_points=300, slice_value=None, default_var_bound=None, show_testpoints=True, **kwds):
        r"""
        Plot the cell.

        - If slice_value is given, it is either a polynomial_map that defines a section, or a list of fixed parameter values with two of them being None. Plot the section. See examples in ``SemialgebraicComplex.plot()``.
        - show_testpoints controls whether to plot the testpoint in this cell.
        - plot_points controls the quality of the plotting.
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
        g = (self.bsa).plot(alpha=alpha, plot_points=plot_points, slice_value=slice_value, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, color=color, fill_color=color, **kwds)
        if show_testpoints and not slice_value:
            pt = self.var_value
            if color == 'white':
                ptcolor = 'black'
            else:
                ptcolor = 'white'
            if (xmin <= pt[0] <= xmax) and (ymin <= pt[1] <= ymax):
                g += point(pt, color = ptcolor, size = 2, zorder=10)
        return g

    def find_neighbour_candidates(self, flip_ineq_step, wall_crossing_method='heuristic', goto_lower_dim=False):
        r"""
        Try flipping exactly one inequality at one time, to reach a new testpoint as neighbour candidate.

        - flip_ineq_step defines the step length
        - wall_crossing_method is 'heuristic' or FIXME: 'mathematica' or 'heuristic_with_check'
        - if goto_lower_dim=False, the cell is considered as its formal closure, so no recursion into test points in lower dimensional cells.

        OUTPUT new_points is a dictionary of dictionaries. The new_points[i] is a dictionay whose keys = candidate neighbour testpoints, values = (bddbsa whose eq_poly has i elements, polynomial_map) of the candidate neighbour cell that contains the candidate neighbour testpoint. bddbsa is recorded so that (1) complex.bddbsa is always respected; and (2) can recursively go into lower dimensional cells. polynomial_map is recorded and passed to the constructor of the neighbour cell.
        """
        num_eq = len(list(self.bddbsa.eq_poly()))
        new_points = {}  #dictionary with key=num_eq, value=dictionay of pt: (bddbsa, polynomial_map).
        # first flip l <= 0 to l > 0.
        for l in list(self.bsa.le_poly()):
            bsa = copy(self.bddbsa) #veronese
            #for ll in self.bsa.eq_poly():
            # not needed because self.bsa and self.bddbsa have the same eq_poly at this moment.
            # inequalites are orthogonal to the equations
            for ll in list(self.bsa.lt_poly()):
                bsa.add_polynomial_constraint(ll, operator.lt)
            for ll in list(self.bsa.le_poly()):
                if ll != l:
                    bsa.add_polynomial_constraint(ll, operator.le)
            try:
                if bsa.is_polynomial_constraint_valid(l, operator.le):
                    # can't flip l<=0, as bsa intersection with l > 0 is empty.
                    continue
            except NotImplementedError:
                pass          
            # find point in intersection(bsa, l>0, l<flip_ineq_step)
            if wall_crossing_method == 'mathematica':
                raise NotImplementedError
            elif wall_crossing_method == 'heuristic_with_check':
                raise NotImplementedError
            else: #None is treated as 'heuristic'
                # do not distinguish lt_poly and le_poly as we wish test_points are in general positions
                pt = find_point_flip_ineq_heuristic(self.var_value, l, list(bsa.lt_poly())+list(bsa.le_poly()), flip_ineq_step)
                if pt is not None:
                    # Find a new point, use polynomial map to recover the values of those eliminated variables.
                    pt_across_wall = tuple(p(pt) for p in self.polynomial_map)
                    if num_eq not in new_points:
                        new_points[num_eq] = OrderedDict()
                    new_points[num_eq][pt_across_wall] = (self.bddbsa, self.polynomial_map)
        # now flip l < 0 to l >= 0
        for l in list(self.bsa.lt_poly()):
            bsa = copy(self.bddbsa) #veronese
            # adding self.bsa.eq_poly() is not needed because self.bsa and self.bddbsa have the same eq_poly at this moment.
            # inequalites are orthogonal to the equations
            for ll in list(self.bsa.le_poly()):
                bsa.add_polynomial_constraint(ll, operator.le)
            for ll in list(self.bsa.lt_poly()):
                if ll != l:
                    bsa.add_polynomial_constraint(ll, operator.lt)
            try:
                if bsa.is_polynomial_constraint_valid(l, operator.lt):
                    # can't flip l<0, as bsa intersection with l >= 0 is empty.
                    continue
            except NotImplementedError:
                pass          
            # find point in intersection(bsa, l>=0, l<flip_ineq_step)
            if wall_crossing_method == 'mathematica':
                raise NotImplementedError
            elif wall_crossing_method == 'heuristic_with_check':
                raise NotImplementedError
            else: #None is treated as 'heuristic'
                # do not distinguish lt_poly and le_poly as we wish test_points are in general positions
                pt = find_point_flip_ineq_heuristic(self.var_value, l, list(bsa.lt_poly())+list(bsa.le_poly()), flip_ineq_step)
                if pt is not None:
                    #find a new point, use polynomial map to recover the values of those eliminated variables.
                    pt_across_wall = tuple(p(pt) for p in self.polynomial_map)
                    if num_eq not in new_points:
                        new_points[num_eq] = OrderedDict()
                    new_points[num_eq][pt_across_wall] = (self.bddbsa, self.polynomial_map)
                if (goto_lower_dim is True) and (l.degree() == 1):
                    # l is linear wall. try to find a point on the wall using heuristic gradient descent method.
                    pt = find_point_on_ineq_heuristic(self.var_value, l, list(bsa.lt_poly()), list(bsa.le_poly()), flip_ineq_step)
                    if pt is not None:
                        pt_on_wall = tuple(p(pt) for p in self.polynomial_map)
                        # pt_on_wall satisfies one more equation l == 0 than the points in self.bsa,
                        # next, compute the polynomial map and restrict the bddbsa of the candidate neighbour cell.
                        # univariate polynomials (Polynomial_rational_flint) does not define "coefficient", but has "monomial_coefficient".
                        lhs = vector(l.monomial_coefficient(x) for x in l.parent().gens())
                        cst = l.constant_coefficient()
                        bddbsa = copy(self.bddbsa)
                        bddbsa.add_linear_constraint(lhs, cst, operator.eq)
                        if cst == 0:
                            v = l.monomials()[-1]
                        else:
                            v = l.monomials()[-2] # order??
                        v_mapped_to = v - l / (l.monomial_coefficient(v))  # eliminate v
                        polynomial_map = [p.subs({v: v_mapped_to}) for p in self.polynomial_map]
                        if (num_eq + 1) not in new_points:
                            new_points[num_eq+1] = OrderedDict()
                        new_points[num_eq+1][pt_on_wall] = (bddbsa, polynomial_map)
        return new_points

    # DELETE. TODO: complete Mathematica methods in find_neighbour_candidates
    # def find_walls_and_new_points(self, flip_ineq_step, wall_crossing_method, goto_lower_dim=False):
    #     r"""
    #     Try flipping exactly one inequality at one time, to reach a new testpoint in a neighbour cell.

    #     Discard the wall if it is impossible to do so, as the wall is redundant.

    #     - flip_ineq_step defines the step length
    #     - wall_crossing_method is 'heuristic' or 'mathematica' or 'heuristic_with_check'
    #     - goto_lower_dim tells whether it also finds new testpoints on the walls.
    #     """
    #     walls = []
    #     new_points = {}
    #     bddllt = list(self.bddbsa.lt_poly())
    #     bddlle = list(self.bddbsa.le_poly())
    #     if self.bsa.eq_poly():
    #         bddllt, bddlle = substitute_llt_lle(bddllt, bddlle, self.var_map, self.var_name, self.var_value)
    #         bddlin = bddllt+bddlle
    #     else:
    #         bddlin = list(bddllt) + list(bddlle)
    #     #FIXME: Mixing the strict and non-strict inequalities
    #     selflin = list(self.bsa.lt_poly()) + list(self.bsa.le_poly())
    #     # decide which inequalities among self.lin are walls (irredundant).
    #     for i in range(len(selflin)):
    #         ineq = selflin[i]
    #         ineqs = selflin[i+1::] + bddlin
    #         if ineq in ineqs:
    #             continue
    #         if wall_crossing_method == 'mathematica':
    #             ineqs = walls + ineqs
    #             condstr_others = write_mathematica_constraints(self.bsa.eq_poly(), ineqs)
    #             # maybe shouldn't put self.bsa.eq_poly() into FindInstance, but solve using var_map later.
    #             condstr_ineq = '0<'+str(ineq)+'<'+str(flip_ineq_step)
    #             pt_across_wall = find_instance_mathematica(condstr_others + condstr_ineq, self.var_name)
    #         else:
    #             if wall_crossing_method == 'heuristic_with_check':
    #                 ineqs = walls + ineqs
    #             else:
    #                 #less clever, more careful, for 'heuristic'
    #                 ineqs = selflin[:i] + ineqs
    #             pt = find_point_flip_ineq_heuristic(self.var_value, ineq, ineqs, flip_ineq_step)
    #             if pt is None:
    #                 if wall_crossing_method == 'heuristic_with_check':
    #                     condstr = write_mathematica_constraints(self.bsa.eq_poly(), ineqs) + \
    #                               '0<'+str(ineq)+'<'+str(flip_ineq_step) # Need the last inequality? see hyperbole ex.
    #                     pt_across_wall = find_instance_mathematica(condstr, self.var_name)
    #                 else:
    #                     pt_across_wall = None
    #             else:
    #                 if self.bsa.eq_poly():
    #                     pt_across_wall = tuple(self.var_map[v](pt) for v in ineq.args())
    #                 else:
    #                     pt_across_wall = pt
    #         if pt_across_wall is None:
    #             # ineq is not a wall
    #             continue
    #         walls.append(ineq)
    #         if not ((wall_crossing_method == 'heuristic_with_check') and (pt is None)):
    #             # Hope that the new pt_across_wall is in a general position, so that running the function on it does not generate new equations. Otherwise bfs got stuck in low dimension. Two solutions: take smaller flip_ineq_step, or do bfs with check_completion=True.
    #             new_points[pt_across_wall] = list(self.bsa.eq_poly())
    #         if goto_lower_dim is True:
    #             pt_on_wall = None
    #             if wall_crossing_method == 'mathematica':
    #                 # wall could be non-linear, contrasting the assumption above. 
    #                 condstr_ineq = str(ineq) + '==0'
    #                 pt_on_wall = find_instance_mathematica(condstr_others + condstr_ineq, self.var_name)
    #             elif ineq.degree() == 1:
    #                 # is linear wall. try to find a point on the wall using heuristic gradient descent method.
    #                 pt = find_point_on_ineq_heuristic(pt_across_wall, ineq, ineqs, [], flip_ineq_step)
    #                 if pt is None:
    #                     pt_on_wall = None
    #                 elif self.bsa.eq_poly():
    #                     pt_on_wall = tuple(self.var_map[v](pt) for v in ineq.args())
    #                 else:
    #                     pt_on_wall = pt
    #             if not pt_on_wall is None:
    #                 new_points[pt_on_wall] = list(self.bsa.eq_poly()) + [ineq]
    #     return walls, new_points

    # TODO
    # def discard_redundant_inequalities(self, flip_ineq_step=0):
    #     r"""
    #     Use Mathematica's ``FindInstance`` to discard redundant inequalities in the cell description.
    #     TO UPDATE with bsa
    #     """
    #     walls = []
    #     bddllt = self.bddbsa.lt_poly()
    #     bddlle = self.bddbsa.le_poly()
    #     if self.bsa.eq_poly():
    #         bddllt, bddlle = substitute_llt_lle(bddllt, bddlle, self.var_map, self.var_name, self.var_value)
    #         bddlin = bddllt+bddlle
    #     else:
    #         bddlin = list(bddllt) + list(bddlle)
    #     #FIXME: Mixing the strict and non-strict inequalities
    #     selflin = list(self.bsa.lt_poly()) + list(self.bsa.le_poly())
    #     # decide which inequalities among self.lin are walls (irredundant).
    #     for i in range(len(selflin)):
    #         ineq = selflin[i]
    #         ineqs = selflin[i+1::] + bddlin
    #         if ineq in ineqs:
    #             continue
    #         ineqs = walls + ineqs
    #         pt_across_wall = None
    #         for pt in self.neighbour_points:
    #             if (ineq(pt) > 0) and all(l(pt) < 0 for l in ineqs):
    #                 pt_across_wall = pt
    #                 break
    #         if pt_across_wall is None:
    #             condstr_others = write_mathematica_constraints(self.bsa.eq_poly(), ineqs)
    #             # maybe shouldn't put self.leq into FindInstance, but solve using var_map later.
    #             condstr_ineq = '0<'+str(ineq)
    #             if flip_ineq_step > 0:  # experimental. see hyperbole example where wall x=1/2 is not unique; discard it in this case.
    #                 condstr_ineq += '<'+str(flip_ineq_step)
    #             pt_across_wall = find_instance_mathematica(condstr_others + condstr_ineq, self.var_name)
    #         if not (pt_across_wall is None):
    #             walls.append(ineq)
    #     self._lt = set(walls)

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
    def _construct_field_and_test_point(family, var_value):
        r"""
        Construct a :class:`ParametricRealField` K using ``var_value``.

        ``var_value`` is a list parallel to ``family.names()``.

        Construct a test_point of type dictionary, which maps each parameter of the family to the corresponding :class:`ParametricRealFieldElement` if this is a parameter of K, otherwise maps to the default argument value.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: family=ParametricFamily(gmic); var_value=[1/2];
            sage: K, test_point = ProofCell._construct_field_and_test_point(family, var_value)
            sage: K
            ParametricRealField(names = ['f'], values = [1/2])
            sage: test_point
            {'conditioncheck': False,
                 'f': f~,
                 'field': ParametricRealField(names = ['f'], values = [1/2])}
        """
        var_name = family.names()
        K = ParametricRealField(var_value, var_name)
        test_point = family.args_from_point(K.gens())
        args_set = family.default_values()
        if 'field' in args_set:
            test_point['field'] = K
        if 'conditioncheck' in args_set:
            test_point['conditioncheck'] = False
        if 'merge' in args_set:
            test_point['merge'] = False
        #if bddbsa is None:
        #    bddbsa =  BasicSemialgebraicSet_eq_lt_le_sets(poly_ring=PolynomialRing(QQ, var_name)) # class doesn't matter, just to record the constraints on K._bsa
        #DELETE and check if the code still works.
        #K.add_initial_space_dim() #so that the parameters var_name are the first ones in the monomial list. Needed for ppl variable elimination, otherwise ineqs are not orthogonal to eliminated variables. FIXME: Get rid of this requirement. Also needed for Mccormicks?
        #assert (K.gens() in bddbsa) # record boundary constraints.
        return K, test_point

    def __init__(self, family, var_value, find_region_type, bddbsa, polynomial_map):
        K, test_point = self._construct_field_and_test_point(family, var_value)
        try:
            h = family(**test_point)
        except Exception:
            # Function is non-constructible at this random point.
            if debug_cell_exceptions:
                import pdb; pdb.post_mortem()
            h = None
        region_type = find_region_type(K, h)
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

    The entire parameter space is covered by cells::

        sage: complex.is_complete(strict=True)                # optional - mathematica
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
        sage: extc.is_complete(bddlin=boundary,strict=True) # optional - mathematica
        False
        sage: extc.is_complete(bddlin=boundary,strict=False) # optional - mathematica
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
        sage: extc.is_complete(bddlin=boundary,strict=True)                         #long time, optional - mathematica
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
        if isinstance(family, ParametricFamily_base):
            assert var_name is None
            assert not opt_non_default
            assert not kwds_dict
        else:
            default_values = copy(opt_non_default)
            default_values.update(kwds_dict)
            family = ParametricFamily(family, default_values=default_values, names=var_name)
        self.family = family
        var_name = family.names()
        self.d = len(var_name)
        self.graph = Graphics()
        self.num_plotted_components = 0
        self.points_to_test = [OrderedDict() for i in range(self.d + 1)] # a list of dictionaries of the form {testpoint: (bddbsa, polynomial_map)}. The i-th dictionary in the list corresponds to the cells that are expected to have i equations.
        # we record bddbsa for each testpoint, since we want to restrict to lower dimensional cells when  goto_lower_dim is set to True.
        if find_region_type is None:
            find_region_type = find_region_type_igp
        self.find_region_type = find_region_type
        self.default_var_bound = default_var_bound
        if bddbsa is None:   #HAS BUG: r = regions[31]; theta = thetas[16]; cpl_complex = cpl_fill_region_given_theta(r, theta); self.bddbsa = BasicSemialgebraicSet_veronese(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(Constraint_System {x0+6*x1-1==0, -x0+1>0, 2*x0-1>0}), polynomial_map=[f, z]); self.bddbsa.polynomial_map() is [f, z]; self.bddbsa.ambient_dim() is 1.
            self.bddbsa = BasicSemialgebraicSet_veronese(poly_ring=PolynomialRing(QQ, var_name))
        else:
            self.bddbsa = bddbsa
        if polynomial_map is None:
            eqs = list(self.bddbsa.eq_poly())
            if eqs:
                #Only useful for treating careless input with low dim bddbsa provided but not polynomial_map. # should not happen in cpl because we pass r.polynomal_map to cpl_complex.
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
            if (x is not None) and (x in self.bddbsa):
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
            # if any(dic for dic in self.points_to_test):
            #     var_value = list(self.points_to_test.popitem(last=False)[0]) # TO CHANGE  # last=False -> BFS.
            # else:
            var_value = self.generate_random_var_value(var_bounds=var_bounds)
            # This point is not already covered.
            if self.is_point_covered(var_value):
                num_failings += 1
            else:
                return var_value
        logging.warning("The complex has %s cells. Cannot find one more uncovered point by shooting %s random points" % (len(self.components), max_failings))
        return None

    # TODO
    # def find_uncovered_point_mathematica(self, strict=True, bddstrict=True):
    #     r"""
    #     Call Mathematica ``FindInstance`` to get a point that satisfies self.bddbsa
    #     and is uncovered by any cells in the complex.

    #     The argument strict controls whether inequalities are treated as <= 0 or as < 0.
    #     If such point does not exist, return ``None``.
 
    #     EXAMPLES::

    #         sage: from cutgeneratingfunctionology.igp import *
    #         sage: logging.disable(logging.WARN)
    #         sage: complex = SemialgebraicComplex(lambda x,y: max(x,y), ['x','y'], find_region_type=result_symbolic_expression, default_var_bound=(-10,10))            # optional - mathematica
    #         sage: complex.add_new_component([1,2], bddbsa=None, flip_ineq_step=0, wall_crossing_method=None, goto_lower_dim=False) # the cell {(x,y): x<y}                          # optional - mathematica
    #         sage: var_value = complex.find_uncovered_point_mathematica(strict=False)  # optional - mathematica
    #         sage: complex.is_point_covered(var_value)                   # optional - mathematica
    #         False
    #     """
    #     if not bddleq and not bddlin:
    #         #FIXME: mixing strict and non-strict inequalities
    #         condstr = write_mathematica_constraints(self.bddbsa.eq_poly(), (self.bddbsa.lt_poly()).union(bddbsa.le_poly()), strict=True) #why strict = strict doesn't work when goto_lower_dim=False?
    #     else:
    #         condstr = write_mathematica_constraints(bddleq, bddlin, strict=bddstrict)
    #     for c in self.components:
    #         condstr_c = write_mathematica_constraints(c.bsa.eq_poly(), (c.bsa.lt_poly()).union(c.bsa.le_poly()), strict=strict)
    #         if condstr_c:
    #             condstr += '!(' + condstr_c[:-4] + ') && '
    #         else:
    #             # c.lin == [] and c.leq == [], cell covers the whole space.
    #             return None
    #     if not condstr:
    #         return tuple([0]*len(self.var_name))
    #     return find_instance_mathematica(condstr[:-4], self.var_name)


    def add_new_component(self, var_value, bddbsa=None, polynomial_map=None, flip_ineq_step=0, wall_crossing_method=None, goto_lower_dim=False):
        r"""
        Compute one proof cell around var_value. Append this cell to the complex.

        - bddbsa defines the boundary constraints satisfied by the cell. By default, bddbsa=None, and complex.bddbsa is taken.
        - If flip_ineq_step = 0, do not search for neighbour testpoints. Used in ``shoot_random_points()``. wall_crossing_method = ``None`` in this case.
        - If flip_ineq_step > 0, search for neighbour testpoints using wall_crossing_method = 'mathematica' or 'heuristic' or 'heuristic_with_check'. Used in bfs. If goto_lower_dim is ``False`` or (goto_lower_dim is ``True`` and wall_crossing method = 'heuristic' but wall is non-linear), then find new testpoints across the wall only.

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
            sage: complex.components[0].bsa.lt_poly()                 # optional - mathematica
            [y^2 - x]
            sage: complex.points_to_test                              # optional - mathematica
            OrderedDict([((19/20, 1), (None, None))])
        """
        if bddbsa is None:
            bddbsa = self.bddbsa
        if polynomial_map is None:
            polynomial_map = self.polynomial_map
        new_component = self._cell_class(self.family, var_value, 
                                         find_region_type=self.find_region_type, bddbsa=bddbsa, polynomial_map=polynomial_map)
        num_eq = len(list(new_component.bddbsa.eq_poly()))
        if len(list(new_component.bsa.eq_poly())) > num_eq:
            # bsa is lower dimensional as it has more equations than bddbsa, 
            # so we try to perturb the testpoint to obtain a
            # new testpoint in bddbsa that does not fall into a lower dimensional cell.
            # Heuristic code using gradient desecent.
            for l in (set(new_component.bsa.eq_poly())- set(new_component.bddbsa.eq_poly())):
                pert_value = tuple(var_value[i] + gradient(l)[i](var_value) / 1000 for i in range(self.d))
                if not (pert_value in self.points_to_test[num_eq]) and (pert_value in new_component.bddbsa):
                    self.points_to_test[num_eq][pert_value] = (copy(new_component.bddbsa), copy(new_component.polynomial_map))
                pert_value = tuple(var_value[i] - gradient(l)[i](var_value) / 1000 for i in range(self.d))
                if not (pert_value in self.points_to_test[num_eq]) and (pert_value in new_component.bddbsa):
                    self.points_to_test[num_eq][pert_value] = (copy(new_component.bddbsa), copy(new_component.polynomial_map))
            logging.warning("The cell around %s defined by %s has more equations than boundary %s" %(new_component.var_value, new_component.bsa, new_component.bddbsa))
            #if not goto_lower_dim:  # quit anyway because bddbsa dim unknown.
            return
        elif len(list(new_component.bsa.eq_poly())) < num_eq:
            raise ValueError()
        elif (flip_ineq_step != 0) and (new_component.region_type != 'stop'):
            # when using random shooting, don't generate neighbour points;
            new_points = new_component.find_neighbour_candidates(flip_ineq_step, wall_crossing_method, goto_lower_dim)
            new_component.neighbour_points = []
            for (n, new_points_n_eq) in new_points.items():
                self.points_to_test[n].update(new_points_n_eq)
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

    def plot(self, alpha=0.5, plot_points=300, slice_value=None, restart=True, **kwds):
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
        if 'xmax' in kwds:
            self.graph.xmax(kwds['xmax'])
        if 'ymin' in kwds:
            self.graph.ymin(kwds['ymin'])
        if 'ymax' in kwds:
            self.graph.ymax(kwds['ymax'])
        for c in self.components[self.num_plotted_components::]:
            gc = c.plot(alpha=alpha, plot_points=plot_points, slice_value=slice_value, default_var_bound=self.default_var_bound, **kwds)
            if gc: # need this because (empty g + empty gc) forgets about xmin xmax ymin ymax.
                self.graph += gc
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
        - wall_crossing_method: 'mathematica' or 'heuristic' or 'heuristic_with_check' 
              - wall_crossing_method='heuristic_with_check': if heuristic wall-crossing does not find a new testpoint, then use Mathematica to check if the ineq is not a wall).
        - goto_lower_dim: whether lower dimensional cells are considered. If goto_lower_dim is ``False`` or if goto_lower_dim is ``True`` and wall_crossing method is 'heuristic' but the wall is non-linear, then find new testpoint across the wall only.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.WARN)
            sage: complex = SemialgebraicComplex(lambda x,y: min(x^2+ 4*y^2, 4), ['x','y'], find_region_type=result_symbolic_expression, default_var_bound=(-3,3))    # optional - mathematica
            sage: complex.bfs_completion(check_completion=True, wall_crossing_method='mathematica', goto_lower_dim=True)                                                          # optional - mathematica
            sage: len(complex.components)                               # optional - mathematica
            3
            sage: complex.plot()                                        # not tested

            sage: complex = SemialgebraicComplex(lambda x,y: min(x+ 4*y, 4), ['x','y'], find_region_type=result_symbolic_expression, default_var_bound=(-5,5))
            sage: complex.bfs_completion(var_value=[1,1])
            sage: len(complex.components)
            2
            sage: complex.plot()                                        # not tested

        See more examples in the docstring of the class :class:`SemialgebraicComplex`.
        """
        if not any(dic for dic in self.points_to_test) and not var_value:
            var_value = self.find_uncovered_random_point()
        if var_value:
            # put given var_value to num_eq=0 so that it pops out first in bfs.
            self.points_to_test[0][tuple(var_value)] = (self.bddbsa, None) 
        num_eq = 0
        while num_eq <= self.d:
            if self.points_to_test[num_eq]: # and len(self.components)<10:
                var_value, (bddbsa, polynomial_map) = self.points_to_test[num_eq].popitem(last=False)  # BFS
                var_value = list(var_value)
                if not self.is_point_covered(var_value):
                    self.add_new_component(var_value, bddbsa=bddbsa, polynomial_map=polynomial_map, flip_ineq_step=flip_ineq_step, wall_crossing_method=wall_crossing_method, goto_lower_dim=goto_lower_dim)
            else:
                num_eq += 1
        if check_completion:
            if max_failings == 0:
                raise NotImplementedError()
                #uncovered_pt = self.find_uncovered_point_mathematica(strict=goto_lower_dim)
            else: # assume that check_completion is an integer.
                uncovered_pt = self.find_uncovered_random_point(max_failings=max_failings)
            if uncovered_pt is not None:
                logging.warning("After bfs, the complex has uncovered point %s." % (uncovered_pt,))
                self.bfs_completion(var_value=uncovered_pt, \
                                    flip_ineq_step=flip_ineq_step, \
                                    check_completion=check_completion, \
                                    wall_crossing_method=wall_crossing_method, \
                                    goto_lower_dim=goto_lower_dim)

    # TODO
    # def is_complete(self, strict=False, bddstrict=True):
    #     r"""
    #     Return whether the entire parameter space satisfying self.bddbsa is covered by cells.

    #     EXAMPLES::

    #         sage: from cutgeneratingfunctionology.igp import *
    #         sage: logging.disable(logging.WARN)
    #         sage: complex = SemialgebraicComplex(lambda x,y: min(x+ 4*y, 4), ['x','y'], find_region_type=result_symbolic_expression, default_var_bound=(-5,5))    # optional - mathematica
    #         sage: complex.bfs_completion(var_value=[1,1])           # optional - mathematica
    #         sage: complex.is_complete()                             # optional - mathematica
    #         True
    #         sage: complex.is_complete(strict=True)                  # optional - mathematica
    #         False
    #     """
    #     if self.find_uncovered_point_mathematica(strict=strict, bddstrict=bddstrict) is None:
    #         return True
    #     else:
    #         return False

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

    # def discard_redundant_inequalities(self, flip_ineq_step=0):
    #     # TO UPDATE with bsa
    #     for c in self.components:
    #         c.discard_redundant_inequalities(flip_ineq_step=flip_ineq_step)

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

region_type_color_map = [('not_constructible', 'white'), ('is_constructible', 'black'), ('not_minimal', 'orange'), ('is_minimal', 'darkgrey'),('not_extreme', 'green'), ('is_extreme', 'blue'), ('stop', 'grey'), (True, 'blue'), (False, 'red'), ('constructible', 'darkgrey'), ('extreme', 'red')]

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
    The current_var_value satisfies that l(current_var_value) < 0 for l=ineq and for every l in ineqs, 
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

        sage: x, y = find_point_flip_ineq_heuristic([1,9/10], a+b-2, [-a+b^2], 1/2)  # got (123/85, 179/170)
        sage: float(x), float(y)    # abs tol 1e-10
        (1.4470588235294117, 1.0529411764705883)
        sage: x, y = find_point_flip_ineq_heuristic([11/40,1/2], a+b-2, [-a+b^2], 1/4) # got (39295901/31739294, 125037049/123564610), now get (45793/36987, 34307/33903)
        sage: (x + y - 2 > 0) and (-x + y^2 < 0)
        True

    Ineq is a redundant inequality, it's impossible to cross it without crossing any other ineqs.
    Thus, return ``None``::

        sage: find_point_flip_ineq_heuristic([1,1/2], a+b-2, [-a+b^2, a-1], 1/4)

    Bug example in cpl Cell 9 with test point (499/1250, 488072439572/4866126017667). Output was (333/1000, 167/1500) which violates ineqs[-1](pt) < 0. BUG! Notice that using sorted(ineqs) with the same current_var_value input gave a correct output. See def adjust_pt_to_satisfy_ineqs below::

        sage: P.<f,z>=QQ[]; ineq = -3*f + 1; flip_ineq_step = 1/1000; current_var_value = (499/1250, 488072439572/4866126017667)
        sage: ineqs = [2*f + 2*z - 1, f + 5*z - 1, -f - 6*z + 1, -2*f - 3*z + 1]
        sage: pt = find_point_flip_ineq_heuristic(current_var_value, ineq, ineqs, flip_ineq_step)
        sage: all(l(pt) < 0 for l in ineqs) and ineq(pt)>0
        True
    """
    # heuristic method.
    ineq_gradient = gradient(ineq)
    current_point = vector(RR(x) for x in current_var_value) # Real numbers, faster than QQ
    ineq_value = ineq(*current_point)
    try_before_fail = 10000 # define maximum number of walks.
    while (ineq_value <= 0) and (try_before_fail > 0):
        ineq_direction = vector(g(*current_point) for g in ineq_gradient)
        if ineq.degree() == 1:
            step_length = (-ineq(*current_point)+flip_ineq_step) / (ineq_direction * ineq_direction)
        else:
            step_length = flip_ineq_step / (ineq_direction * ineq_direction) # ineq_value increases by flip_ineq_step=0.01 roughly
            if step_length > 1:
                step_length = 1  # ensure that distance of move <= sqrt(flip_ineq_step) = 0.1 in each step
        current_point += step_length * ineq_direction
        new_ineq_value = ineq(*current_point)
        if new_ineq_value <= ineq_value:
            return None
        ineq_value = new_ineq_value
        try_before_fail -= 1
        #print current_point, RR(ineq_value)
    if ineq_value <= 0:
        return None
    new_point = adjust_pt_to_satisfy_ineqs(current_point, ineq, ineqs, [], flip_ineq_step)
    # if new_point is not None and ineq(*new_point) <= 0:
    #     #logging.info("Didn't add %s because it violates %s > 0" % (new_point, ineq))
    #     return None
    return new_point #type is tuple


def find_point_on_ineq_heuristic(current_var_value, ineq, strict_ineqs, nonstrict_ineqs, flip_ineq_step):
    r"""
    The current_var_value satisfies that l(current_var_value) < 0 for l=ineq,
    l(current_var_value) < 0 for all l in strict_ineqs, and
    l(current_var_value) <= 0 for all l in strict_ineqs, and
    where ineq is a polynomial and ineqs is a list of polynomials.

    Assume that ineq is linear.

    Use heuristic method (gradient descent method with given small positive step length flip_ineq_step)
    to find a new_point (type is tuple) such that
    ineq(new_point) == 0 and l(new_point) < 0 for all l in strict_ineqs, and l(new_point) <= 0 for all l in nonstrict_ineqs.
    Return new_point, or ``None`` if it fails to find one.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: P.<a,b>=QQ[]
        sage: find_point_on_ineq_heuristic([1,1/2], a+b-2, [-a+b^2], [], 1/4)
        (5/4, 3/4)
        sage: find_point_on_ineq_heuristic([11/40,1/2], a+b-2, [-a+b^2], [], 1/4)
        (171073319/163479120, 155884921/163479120)
    """
    ineq_gradient = gradient(ineq)
    current_point = vector(current_var_value)
    ineq_direction = vector(g(*current_point) for g in ineq_gradient)
    step_length = -ineq(*current_point) / (ineq_direction * ineq_direction)
    current_point += step_length * ineq_direction
    ineq_value = ineq(*current_point)
    if ineq_value != 0:
        return None
    new_point = adjust_pt_to_satisfy_ineqs(current_point, ineq, strict_ineqs, nonstrict_ineqs, flip_ineq_step)
    if new_point is not None and ineq(*new_point) != 0:
        #logging.info("Didn't add %s because it violates %s == 0" % (new_point, ineq))
        return None
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
        sage: ineq = a+b-2
        sage: adjust_pt_to_satisfy_ineqs(vector([13/10,12/10]), ineq, [-a+b^2], [], 1/2)
        (123/85, 179/170)
        sage: adjust_pt_to_satisfy_ineqs(vector([71/80, 89/80]), ineq, [-a+b^2], [], 1/4)
        (171073319/163479120, 155884921/163479120)

    If impossible, return ``None``::

        sage: adjust_pt_to_satisfy_ineqs(vector([11/8, 7/8]), ineq,[-a+b^2], [a-1], 1/4)

    Bug example in cpl Cell 9 with test point (499/1250, 488072439572/4866126017667). Without converting input to QQ, output was (0.333000000000000, 0.111333333333333) with  -2*f - 3*z + 1 = -5.55111512312578e-17 , the QQ of the output=(333/1000, 167/1500) has -2*f - 3*z + 1 == 0. Revise the code to take QQ input current point:: 

        sage: P.<f,z>=QQ[]
        sage: current_point = vector([0.333000000000000, 0.100300000000000]); ineq = -3*f + 1; strict_ineqs = [2*f + 2*z - 1, f + 5*z - 1, -f - 6*z + 1, -2*f - 3*z + 1]; nonstrict_ineqs = []; flip_ineq_step = 1/1000
        sage: pt = adjust_pt_to_satisfy_ineqs(current_point, ineq, strict_ineqs, nonstrict_ineqs, flip_ineq_step)
        sage: all(l(pt) < 0 for l in strict_ineqs) and all(l(pt) <= 0 for l in nonstrict_ineqs) and ineq(pt)>0
        True
    """
    #current_point is a vector
    ineq_gradient = gradient(ineq)
    max_walks = min(ceil(2/flip_ineq_step), 1000) # define maximum number of walks.
    ineqs_and_strictness = [(l, True) for l in strict_ineqs] + [(l, False) for l in nonstrict_ineqs]
    for (l, strictness) in ineqs_and_strictness:
        l_gradient = gradient(l)
        l_value = l(*current_point)
        try_before_fail = max_walks
        while ((strictness and (l_value >= 0)) or (not strictness and (l_value > 0))) and (try_before_fail > 0):
            l_direction = vector(-g(*current_point) for g in l_gradient) #decrease l_value
            ineq_direction = vector(g(*current_point) for g in ineq_gradient)
            s = (ineq_direction * l_direction) / (ineq_direction * ineq_direction)
            projected_direction = l_direction - s * ineq_direction # want that ineq_value remains the same
            if projected_direction == 0:
                return None
            if l.degree() == 1:
                step_length = (l_value+flip_ineq_step) / (projected_direction * l_direction)
            else:
                step_length = flip_ineq_step / (projected_direction * l_direction) # l_value decreases by 0.01 roughly
                # if step_length * norm(projected_direction) >= 1:  # move too far  # is 1 a good value here?? why this if?
                #     return None
            #flip_ineq_step = flip_ineq_step / 2
            current_point += step_length * projected_direction
            l_value = l(*current_point)
            try_before_fail -= 1
            #print current_point, RR(l_value)
        if (strictness and (l_value >= 0)) or (not strictness and (l_value > 0)) or (ineq(*current_point) < 0):
            return None
    for l in strict_ineqs:
        if l(*current_point) >= 0:
            return None
    for l in nonstrict_ineqs:
        if l(*current_point) > 0:
            return None
    if all(x.parent()==QQ for x in current_point):
        return tuple(current_point)
    else:
        current_point = vector(QQ(x.n(30)) for x in current_point)
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
        sage: embed_function_into_family(given_function, parametric_family, check_completion=True) # optional - mathematica
        {'f': 2/5, 'lambda_1': 6/25, 'lambda_2': 2/75}
    """
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
    complex.points_to_test[0][tuple(var_value)] = (None, None) # FIXME: num_eq = 0?? It's okay, num_eq = 0 let it pops out first.
    flip_ineq_step = opt.pop('flip_ineq_step', 1/1000)
    is_complete = False
    while not is_complete:
        num_eq = 0
        while num_eq <= complex.d:
            if complex.points_to_test[num_eq]: # and len(self.components)<10:
                var_value, (bddbsa, polynomial_map) = complex.points_to_test[num_eq].popitem(last=False)
                var_value = list(var_value)
                if not complex.is_point_covered(var_value):
                    complex.add_new_component(var_value, bddbsa=bddbsa, polynomial_map=polynomial_map, flip_ineq_step=flip_ineq_step, goto_lower_dim=True, **opt)
                    if complex.components and complex.components[-1].region_type is True:
                        dic = {var_name[i]:var_value[i] for i in range(len(var_name))}
                        return dic
            else:
                num_eq += 1
        if check_completion:
            var_value = complex.find_uncovered_point_mathematica(strict=True)
            if not var_value:
                is_complete = True
            else:
                complex.points_to_test[0][var_value]=(None, None)
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
