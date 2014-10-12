import sage.structure.element
from sage.structure.element import FieldElement

default_symbolic_field = None

class SymbolicRNFElement(FieldElement):

    def __init__(self, value, symbolic=None, parent=None):
        if parent is None:
            raise ValueError, "SymbolicRNFElement invoked with parent=None. That's asking for trouble"
            parent = default_symbolic_field
        FieldElement.__init__(self, parent) ## this is so that canonical_coercion works.
        self._val = value
        if symbolic is None:
            self._sym = value # changed to not coerce into SR. -mkoeppe
        else:
            self._sym = symbolic
        self._parent = parent ## this is so that .parent() works.

    def sym(self):
        return self._sym

    def val(self):
        return self._val

    def parent(self):
        return self._parent

    def __cmp__(left, right):
        result = cmp(left._val, right._val)
        # print  "__cmp__(%s, %s) = %s" %(left, right, result)
        if result == 0:
            left.parent().record_to_eq_list(left.sym() - right.sym())
        elif result == -1:
            left.parent().record_to_lt_list(left.sym() - right.sym())
        elif result == 1:
            left.parent().record_to_lt_list(right.sym() - left.sym())
        return result

    def __richcmp__(left, right, op):
        result = left._val.__richcmp__(right._val, op)
        # print  "__richcmp__(%s, %s, %s) = %s" %(left, right, op, result)
        return result

    def __abs__(self):
        if self.sign() >= 0:
            return self
        else:
            return -self

    def sign(self):
        parent = self._val.parent()
        return cmp(self._val, parent._zero_element)

    def __repr__(self):
        return repr(self._sym)

    def _latex_(self):
        return "%s" % (latex(self._sym))

    def __add__(self, other):
        if not isinstance(other, SymbolicRNFElement):
            other = SymbolicRNFElement(other, parent=self.parent())
        return SymbolicRNFElement(self._val + other._val, self._sym + other._sym, parent=self.parent())

    def __sub__(self, other):
        if not isinstance(other, SymbolicRNFElement):
            other = SymbolicRNFElement(other, parent=self.parent())
        return SymbolicRNFElement(self._val - other._val, self._sym - other._sym, parent=self.parent())

    def __neg__(self):
        return SymbolicRNFElement(-self._val, -self._sym, parent=self.parent())

    def __mul__(self, other):
        if not isinstance(other, SymbolicRNFElement):
            other = SymbolicRNFElement(other, parent=self.parent())
        return SymbolicRNFElement(self._val * other._val, self._sym * other._sym, parent=self.parent())

    def __div__(self, other):
        if not isinstance(other, SymbolicRNFElement):
            other = SymbolicRNFElement(other, parent=self.parent())
        return SymbolicRNFElement(self._val / other._val, self._sym / other._sym, parent=self.parent())


from sage.rings.ring import Field

import sage.rings.number_field.number_field_base as number_field_base

from sage.structure.coerce_maps import CallableConvertMap

from itertools import izip

class SymbolicRealNumberField(number_field_base.NumberField):
    """
    Parametric search:
    EXAMPLES::

        sage: K.<f> = SymbolicRealNumberField([4/5])
        sage: h = gmic(f, field=K)
        sage: generate_maximal_additive_faces(h)
        sage: K.get_eq_list()
        set([0])
        sage: K.get_lt_list()
        set([2*f - 2, f - 2, -f, f - 1, -1/f, -f - 1, -2*f, -2*f + 1, -1/(-f^2 + f), -1])

        sage: K.<f, lam> = SymbolicRealNumberField([4/5, 1/6])
        sage: h = gj_2_slope(f, lam, field=K)
    """

    def __init__(self, values=[], names=()):
        NumberField.__init__(self)
        self._element_class = SymbolicRNFElement
        self._zero_element = SymbolicRNFElement(0, parent=self)
        self._one_element =  SymbolicRNFElement(1, parent=self)
        self._eq = set([])
        #self._le = set([])
        self._lt = set([])
        #self._ne = set([])
        vnames = PolynomialRing(QQ, names).fraction_field().gens();
        self._gens = [ SymbolicRNFElement(value, name, parent=self) for (value, name) in izip(values, vnames) ]

    def _first_ngens(self, n):
        for i in range(n):
            yield self._gens[i]

    def _an_element_impl(self):
        return SymbolicRNFElement(1, parent=self)
    def _coerce_map_from_(self, S):
        # print "_coerce_map_from: self = %s, S = %s" % (self, S)
        # FIXME: need a coerce map from S to self
        # return CallableConvertMap(S, default_symbolic_field, SymbolicRNFElement, parent_as_first_arg=False)
        return CallableConvertMap(S, self, lambda s: SymbolicRNFElement(s, parent=self), parent_as_first_arg=False)
    def get_eq_list(self):
        return self._eq
    #def get_le_list(self):
    #    return self._le
    def get_lt_list(self):
        return self._lt
    #def get_ne_list(self):
    #    return self._ne
    def record_to_eq_list(self, comparison):
        self._eq.add(comparison)
    #def record_to_le_list(self, comparison):
    #    self._le.add(comparison)
    def record_to_lt_list(self, comparison):
        self._lt.add(comparison)
    #def record_to_ne_list(self, comparison):
    #    self._ne.add(comparison)

default_symbolic_field = SymbolicRealNumberField()
