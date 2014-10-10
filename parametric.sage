## fixed: BUG
# z = SymbolicRNFElement(0); zero = z.parent()._zero_element; z == zero
# returns False, should be True
# dir(zero) 
# raises AttributeError: 'NoneType' object has no attribute 'category' 

import sage.structure.element
from sage.structure.element import FieldElement

default_symbolic_field = None

class SymbolicRNFElement(FieldElement):

    def __init__(self, value, symbolic=None, parent=None):
        if parent is None:
            parent = default_symbolic_field
        FieldElement.__init__(self, parent) ## this is so that canonical_coercion works.
        self._val = value
        if symbolic is None:
            self._sym = SR(value)
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
        #if isinstance(left, SymbolicRNFElement):
        #    leftval = left._val
        #else:
        #    leftval = left
        #if isinstance(right, SymbolicRNFElement):
        #    rightval = right._val
        #else:
        #    rightval = right
        #return cmp(leftval, rightval)
        return cmp(left._val, right._val)

    def __richcmp__(left, right, op):
        #if isinstance(left, SymbolicRNFElement):
        #    leftval = left._val
        #else:
        #    leftval = left
        #if isinstance(right, SymbolicRNFElement):
        #    rightval = right._val
        #else:
        #    rightval = right
        #return leftval.__richcmp__(rightval, op)
        return left._val.__richcmp__(right._val, op)

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
            other = SymbolicRNFElement(other)
        return SymbolicRNFElement(self._val + other._val, self._sym + other._sym)

    def __sub__(self, other):
        if not isinstance(other, SymbolicRNFElement):
            other = SymbolicRNFElement(other)
        return SymbolicRNFElement(self._val - other._val, self._sym - other._sym)

    def __neg__(self):
        return SymbolicRNFElement(-self._val, -self._sym)

    def __mul__(self, other):
        if not isinstance(other, SymbolicRNFElement):
            other = SymbolicRNFElement(other)
        return SymbolicRNFElement(self._val * other._val, self._sym * other._sym)

    def __div__(self, other):
        if not isinstance(other, SymbolicRNFElement):
            other = SymbolicRNFElement(other)
        return SymbolicRNFElement(self._val / other._val, self._sym / other._sym)


from sage.rings.ring import Field

import sage.rings.number_field.number_field_base as number_field_base

from sage.structure.coerce_maps import CallableConvertMap

class SymbolicRealNumberField(number_field_base.NumberField):
    def __init__(self):
        NumberField.__init__(self)
        self._element_class = SymbolicRNFElement
        self._zero_element = SymbolicRNFElement(0, parent=self)
        self._one_element =  SymbolicRNFElement(1, parent=self)
        self._eq = []
        self._le = []
        self._lt = []
        self._ne = []
    def _an_element_impl(self):
        return SymbolicRNFElement(1)
    def _coerce_map_from_(self, S):
        return CallableConvertMap(S, default_symbolic_field, SymbolicRNFElement, parent_as_first_arg=False)
    def get_eq_list(self):
        return self._eq
    def get_le_list(self):
        return self._le
    def get_lt_list(self):
        return self._lt
    def get_ne_list(self):
        return self._ne   
    def record_to_eq_list(self, comparison):
        self._eq.append(comparison)
    def record_to_le_list(self, comparison):
        self._le.append(comparison)
    def record_to_lt_list(self, comparison):
        self._lt.append(comparison)
    def record_to_ne_list(self, comparison):
        self._ne.append(comparison)
    def initialize_comparison_list(self):
        self._eq = []
        self._le = []
        self._lt = []
        self._ne = []

default_symbolic_field = SymbolicRealNumberField()

