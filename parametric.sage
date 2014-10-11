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

class SymbolicRealNumberField(number_field_base.NumberField):
    def __init__(self):
        NumberField.__init__(self)
        self._element_class = SymbolicRNFElement
        self._zero_element = SymbolicRNFElement(0, parent=self)
        self._one_element =  SymbolicRNFElement(1, parent=self)
        self._eq = set([])
        self._le = set([])
        self._lt = set([])
        self._ne = set([])
    def _an_element_impl(self):
        return SymbolicRNFElement(1, parent=self)
    def _coerce_map_from_(self, S):
        print "_coerce_map_from: self = %s, S = %s" % (self, S)
        # FIXME: need a coerce map from S to self
        # return CallableConvertMap(S, default_symbolic_field, SymbolicRNFElement, parent_as_first_arg=False)
        return CallableConvertMap(S, self, lambda s: SymbolicRNFElement(s, parent=self), parent_as_first_arg=False)
    def get_eq_list(self):
        return self._eq
    def get_le_list(self):
        return self._le
    def get_lt_list(self):
        return self._lt
    def get_ne_list(self):
        return self._ne
    def record_to_eq_list(self, comparison):
        self._eq.add(comparison)
    def record_to_le_list(self, comparison):
        self._le.add(comparison)
    def record_to_lt_list(self, comparison):
        self._lt.add(comparison)
    def record_to_ne_list(self, comparison):
        self._ne.add(comparison)

default_symbolic_field = SymbolicRealNumberField()

def parameters_give_same_2d_diagram(fn, f=None):
    """
    Record symbolic comparisons in default_symbolic_field._eq /._le /._lt /._ne
    This version only support continuous functions.

    EXAMPLES::

        sage: f = SymbolicRNFElement(4/5, var('f'))
        sage: h = gmic(f, field=default_symbolic_field)
        sage: default_symbolic_field.initialize_comparison_list()
        sage: parameters_give_same_2d_diagram(h)
        sage: default_symbolic_field.get_eq_list()
        [0, 0, 0, 0, -f/(f - 1) + 1/(f - 1) + 1, -f/(f - 1) + 1/(f - 1) + 1, -f/(f - 1) + 1/(f - 1) + 1, -f/(f - 1) + 1/(f - 1) + 1, 0, 0, 0]
        sage: default_symbolic_field.get_le_list()
        [0, -1, -1, 0, f/(f - 1) - 1/(f - 1) - 1, -f/(f - 1) + 1/(f - 1)]
        sage: default_symbolic_field.get_lt_list()
        [f - 1, -2*f + 1, (f - 1)/f - f/(f - 1) + 1/(f - 1), (2*f - 1)/f - 2]
        sage: default_symbolic_field.get_ne_list()
        []
        sage: lambda_1 = SymbolicRNFElement(1/6, var('lambda_1'))
        sage: h = gj_2_slope(f, lambda_1, field=default_symbolic_field)
        sage: default_symbolic_field.initialize_comparison_list()
        sage: parameters_give_same_2d_diagram(h)
    """
    bkpt = fn.end_points()
    # work on continuous functions
    values = fn.values_at_end_points()

    # breakpoints are orderd from left to right, the first one is 0, the last one is 1
    default_symbolic_field.record_to_eq_list(bkpt[0].sym())
    for i in range(1,len(bkpt)-1):
        default_symbolic_field.record_to_lt_list(bkpt[i].sym() - bkpt[i+1].sym())
    default_symbolic_field.record_to_eq_list(bkpt[-1].sym() - 1)

    # diagonal hyperplanes' position
    for i in range(2, len(bkpt)):
        x_hor_intersection = bkpt[1: i] # sorted list
        x_ver_intersection = [bkpt[i] - bkpt[i - j] for j in range(1, i)] # sorted list
        # x_intersection records the x-coordinates of vertices on diagonal x+y = bkpt[i]
        x_intersection = list(merge(x_hor_intersection, x_ver_intersection))
        for j in range(0, len(x_intersection) - 1):
            if x_intersection[j] < x_intersection[j + 1]:
                default_symbolic_field.record_to_lt_list(x_intersection[j].sym() - x_intersection[j + 1].sym())
            else: # must have x_intersection[j] == x_intersection[j+1]
                default_symbolic_field.record_to_eq_list(x_intersection[j].sym() - x_intersection[j + 1].sym())
    for i in range(1, len(bkpt)-2):
        x_hor_intersection = bkpt[i + 1: len(bkpt) - 1]
        x_ver_intersection = [1 + bkpt[i] - bkpt[len(bkpt) - j] for j in range(i + 1, len(bkpt) - 1)]
        # x_intersection records the x-coordinates of vertices on diagonal x+y = bkpt[i] + 1
        x_intersection = list(merge(x_hor_intersection, x_ver_intersection))
        for j in range(0, len(x_intersection) - 1):
            if x_intersection[j] < x_intersection[j + 1]:
                default_symbolic_field.record_to_lt_list(x_intersection[j].sym() - x_intersection[j + 1].sym())
            else: # must have x_intersection[j] == x_intersection[j+1]
                default_symbolic_field.record_to_eq_list(x_intersection[j].sym() - x_intersection[j + 1].sym())

    # values_at_end_points stay in the range [0,1]
    for i in range(len(bkpt)):
        default_symbolic_field.record_to_le_list(- values[i].sym()) # fn(x) >= 0
        default_symbolic_field.record_to_le_list(values[i].sym() - 1) # fn(x) <= 1

    if f is None:
        f = find_f(fn)
    default_symbolic_field.record_to_eq_list(values[0].sym()) # fn(0) == 0
    # fn(f) == 1 is included in symmetric conditions
    # default_symbolic_field.record_to_eq_list(fn(f).sym() - 1)

    # symmetric conditions
    for i in range(len(bkpt)):
        x = bkpt[i]
        if x == f:
            continue
        if x < f:
            y = f - x
        else:
            y = f + 1 - x
        default_symbolic_field.record_to_eq_list(values[i].sym() + fn(y).sym() - 1)

    # additivity / subadditivity conditions
    additive_vertices_set = generate_additive_vertices(fn)
    for (x, y, z, xeps, yeps, zeps) in additive_vertices_set:
        default_symbolic_field.record_to_eq_list(fn(x).sym() + fn(y).sym() - fn(fractional(z)).sym())
    strict_subadditive_vertices_set = set(itertools.chain( \
                                            generate_type_1_vertices(fn, operator.gt),\
                                            generate_type_2_vertices(fn, operator.gt)))
    for (x, y, z, xeps, yeps, zeps) in strict_subadditive_vertices_set:
        default_symbolic_field.record_to_lt_list(fn(fractional(z)).sym() - fn(x).sym() - fn(y).sym())
    # suppose that fn is extreme for given parameter values, so nonsubadditive_vertiecs_set is empty.
