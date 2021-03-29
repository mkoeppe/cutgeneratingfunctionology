"""
Elements of parametric real fields.
"""

from __future__ import print_function, division, absolute_import

from sage.structure.element import FieldElement
from sage.structure.richcmp import richcmp, op_LT, op_LE, op_EQ, op_NE, op_GT, op_GE
from sage.rings.real_mpfr import RR
from sage.functions.other import ceil, floor
from sage.functions.generalized import sign
import operator

def richcmp_op_negation(op):
    if op == op_LT:
        return op_GE
    elif op == op_LE:
        return op_GT
    elif op == op_EQ:
        return op_NE
    elif op == op_NE:
        return op_EQ
    elif op == op_GT:
        return op_LE
    elif op == op_GE:
        return op_LT
    else:
        raise ValueError("{} is not a valid richcmp operator".format(op))

def format_richcmp_op(op):
    if op == op_LT:
        return '<'
    elif op == op_LE:
        return '<='
    elif op == op_EQ:
        return '=='
    elif op == op_NE:
        return '!='
    elif op == op_GT:
        return '>'
    elif op == op_GE:
        return '>='
    else:
        raise ValueError("{} is not a valid richcmp operator".format(op))

class ParametricRealFieldElement(FieldElement):
    r"""
    A :class:`ParametricRealFieldElement` stores a symbolic expression of the parameters in the problem and a concrete value, which is the evaluation of this expression on the given parameter tuple.

    When a comparison takes place on elements of the class, their concrete values are compared to compute the Boolean return value of the comparison. The constraint on the parameters that gives the same Boolean return value is recorded.
    """

    def __init__(self, parent, value, symbolic=None):
        FieldElement.__init__(self, parent) ## this is so that canonical_coercion works.
        if not parent._mutable_values and value is not None:
            ## Test coercing the value to RR, so that we do not try to build a ParametricRealFieldElement
            ## from something like a tuple or vector or list or variable of a polynomial ring
            ## or something else that does not make any sense.
            # FIXME: parent(value) caused SIGSEGV because of an infinite recursion
            from sage.structure.coerce import py_scalar_parent
            if hasattr(value, 'parent'):
                if not RR.has_coerce_map_from(value.parent()):
                    raise TypeError("Value is of wrong type")
            else:
                if not RR.has_coerce_map_from(py_scalar_parent(type(value))):
                    raise TypeError("Value is of wrong type")
            self._val = value
        if symbolic is None:
            self._sym = value # changed to not coerce into SR. -mkoeppe
        else:
            self._sym = symbolic

    def sym(self):
        return self._sym

    def val(self):
        try:
            return self._val
        except AttributeError:
            return self.parent()._eval_factor(self._sym)

    def _richcmp_(left, right, op):
        r"""
        Examples for traditional cmp semantics::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: K.<f> = ParametricRealField([0], big_cells=False, allow_refinement=True)
            sage: f >= 0
            True
            sage: K._eq
            {f}

        Examples for big_cells semantics::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: K.<f> = ParametricRealField([0], big_cells=True, allow_refinement=False)
            sage: f >= 0
            True
            sage: K._le
            {-f}
            sage: K.get_le_factor()
            {-f}
            sage: K.get_eq_factor()
            set()
            sage: K.get_lt_factor()
            set()

        Examples for big_cells=True, allow_refinement=True::

            sage: K.<f> = ParametricRealField([1], big_cells=True, allow_refinement=True)
            sage: f >= 0
            True
            sage: K._le
            {-f}
            sage: K.get_le_factor()
            {-f}
            sage: K.get_eq_factor()
            set()
            sage: K.get_lt_factor()
            set()

        """
        if not left.parent() is right.parent():
            # shouldn't really happen, within coercion
            raise TypeError("comparing elements from different fields")
        if left.parent()._big_cells:
            result = richcmp(left.val(), right.val(), op)
            if result:
                true_op = op
            else:
                true_op = richcmp_op_negation(op)
            if true_op == op_LT:
                # left.sym() - right.sym() may cancel denominators, but that is
                # OK because _div_ makes sure that denominators are nonzero.
                left.parent().assume_comparison(left - right, operator.lt)
            elif true_op == op_GT:
                left.parent().assume_comparison(left - right, operator.gt)
            elif true_op == op_EQ:
                left.parent().assume_comparison(right - left, operator.eq)
            elif true_op == op_LE:
                left.parent().assume_comparison(left - right, operator.le)
            elif true_op == op_GE:
                left.parent().assume_comparison(left - right, operator.ge)
            elif true_op == op_NE:
                left.parent().assume_comparison(right - left, operator.ne)
            else:
                raise ValueError("{} is not a valid richcmp operator".format(op))
            return result
        else:
            # Traditional cmp semantics.
            if (left.val() == right.val()):
                left.parent().assume_comparison(right, operator.eq, left)
            elif (left.val() < right.val()):
                left.parent().assume_comparison(left, operator.lt, right)
            else:
                left.parent().assume_comparison(right, operator.lt, left)
            return richcmp(left.val(), right.val(), op)

    def __abs__(self):
        """
        Examples for traditional cmp semantics::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: K.<f> = ParametricRealField([-1], big_cells=False, allow_refinement=True)
            sage: abs(f) + abs(-f)
            (-2*f)~
            sage: K._lt
            {f}

        Examples for big_cells semantics::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: K.<f> = ParametricRealField([-1], big_cells=True, allow_refinement=True)
            sage: abs(f) + abs(-f)
            (-2*f)~
            sage: K._lt
            set()
            sage: K._le
            {f}

        """
        preferred_sign = 1
        if self.parent()._big_cells:
            with self.parent().off_the_record():
                if self >= 0:
                    preferred_sign = 1
                else:
                    preferred_sign = -1
        if preferred_sign == 1:
            if self >= 0:
                return self
            else:
                return -self
        else:
            if self <= 0:
                return -self
            else:
                return self

    def sign(self):
        """
        Examples for traditional cmp semantics::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: K.<f, z> = ParametricRealField([-1, 0], big_cells=False, allow_refinement=True)
            sage: sign(f)
            -1
            sage: sign(z)
            0
            sage: K._lt
            {f}
            sage: K._le
            set()
            sage: K._eq
            {z}

        Test that the same result is obtained in big_cells, allow_refinement=True semantics::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: K.<f, z> = ParametricRealField([-1, 0], big_cells=True, allow_refinement=True)
            sage: sign(f)
            -1
            sage: sign(z)
            0
            sage: K._lt
            {f}
            sage: K._le
            set()
            sage: K._eq
            {z}

        Test that the same result is obtained for big_cells, allow_refinement=False semantics::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: K.<f, z> = ParametricRealField([-1, 0], big_cells=True, allow_refinement=False)
            sage: sign(f)
            -1
            sage: sign(z)
            0
            sage: K._lt
            {f}
            sage: K._le
            set()
            sage: K._eq
            {z}
        """
        if self.parent()._big_cells:
            preferred_sign = sign(self.val())    # off the record
        else:
            preferred_sign = 1
        if preferred_sign == 1:
            if self > 0:
                return 1
        elif preferred_sign == -1:
            if self < 0:
                return -1
        else:
            if self == 0:
                return 0
        if self > 0:
            return 1
        elif self < 0:
            return -1
        else:
            return 0

    def floor(self):
        result = floor(self.val())
        result <= self < result + 1
        return result

    def ceil(self):
        result = ceil(self.val())
        result - 1 < self <= result
        return result

    def __float__(self):
        if self.parent().allow_coercion_to_float:
            return float(self.val())
        else:
            raise ValueError("Conversion to float is not allowed")

    def _repr_(self):
        s = self._sym
        r = repr(s)
        if s in RR:
            return r
        if len(r) > 1:
            return '('+r+')~'
        else:
            return r+'~'

    def _latex_(self):
        from sage.misc.latex import latex
        return "%s" % (latex(self._sym))

    def _add_(self, other):
        if not isinstance(other, ParametricRealFieldElement):
            other = ParametricRealFieldElement(self.parent(), other)
        try:
            val = self._val + other._val
        except AttributeError:
            val = None
        return ParametricRealFieldElement(self.parent(), val, self._sym + other._sym)

    def _sub_(self, other):
        if not isinstance(other, ParametricRealFieldElement):
            other = ParametricRealFieldElement(self.parent(), other)
        try:
            val = self._val - other._val
        except AttributeError:
            val = None
        return ParametricRealFieldElement(self.parent(), val, self._sym - other._sym)

    def _neg_(self):
        try:
            val = -self._val
        except AttributeError:
            val = None
        return ParametricRealFieldElement(self.parent(), val, -self._sym)

    def _mul_(self, other):
        if not isinstance(other, ParametricRealFieldElement):
            try:
                other = ParametricRealFieldElement(self.parent(), other)
            except TypeError:
                # For example when other is a vector
                return other * self
        try:
            val = self._val * other._val
        except AttributeError:
            val = None
        return ParametricRealFieldElement(self.parent(), val, self._sym * other._sym)

    def _div_(self, other):
        r"""
        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: K.<f> = ParametricRealField([0])
            sage: 1 / f
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Division by 0 in ParametricRealField
            sage: f / f
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Division by 0 in ParametricRealField
        """
        if not isinstance(other, ParametricRealFieldElement):
            other = ParametricRealFieldElement(self.parent(), other)
        if other == 0:
            raise ZeroDivisionError("Division by 0 in ParametricRealField")
        try:
            val = self._val / other._val
        except AttributeError:
            val = None
        return ParametricRealFieldElement(self.parent(), val, self._sym / other._sym)

    def __hash__(self):
        r"""
        The hash function of these elements must be so that
        elements that would compare equal have the same hash value.

        The constant hash function would do the job.  Instead we use the
        hash of the .val() (because equality implies equality of _val).
        It is not correct to use the hash of the ._sym, or to compare
        the __repr__, because then the user could check for equality
        (for example, by testing the cardinality of a set, as in the
        tests below) without the equality being recorded in the field.

        The correctness of this implementation depends on the guarantee
        of equal _val elements having the same hash value.  If in
        doubt, make sure that the _val elements all come from the same
        field, by ``nice_field_values``.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: K.<f> = ParametricRealField([1])
            sage: s = {f, K(1)}
            sage: len(s)
            1
            sage: K._eq
            {f - 1}
            sage: K.<f> = ParametricRealField([1])
            sage: s = {f, K(2)}
            sage: len(s)
            2
        """
        return hash(self.val())

def is_parametric_element(x):
    # We avoid using isinstance here so that this is robust even if parametric.sage is reloaded.
    # For example, apparently in the test suite.
    return hasattr(x, '_sym')
