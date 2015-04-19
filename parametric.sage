# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

import sage.structure.element
from sage.structure.element import FieldElement

default_symbolic_field = None

class SymbolicRNFElement(FieldElement):

    def __init__(self, value, symbolic=None, parent=None):
        if parent is None:
            raise ValueError, "SymbolicRNFElement invoked with parent=None. That's asking for trouble"
            parent = default_symbolic_field
        FieldElement.__init__(self, parent) ## this is so that canonical_coercion works.
        ## Test coercing the value to RR, so that we do not try to build a SymbolicRNFElement
        ## from something like a tuple or something else that does not make any sense.
        RR(value)
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

    def floor(self):
        result = floor(self._val)
        result <= self < result + 1
        return result

    def ceil(self):
        result = ceil(self._val)
        result - 1 < self <= result
        return result

    def __float__(self):
        return float(self._val)

    def __repr__(self):
        r = repr(self._sym)
        if len(r) > 1:
            return '('+r+')~'
        else:
            return r+'~'

    def _latex_(self):
        return "%s" % (latex(self._sym))

    def __add__(self, other):
        if not isinstance(other, SymbolicRNFElement):
            other = SymbolicRNFElement(other, parent=self.parent())
        return SymbolicRNFElement(self._val + other._val, self._sym + other._sym, parent=self.parent())
    def _add_(self, other):
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
    def _mul_(self, other):
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

#class SymbolicRealNumberField(number_field_base.NumberField):
class SymbolicRealNumberField(Field):
    """
    Parametric search:
    EXAMPLES::

        sage: K.<f> = SymbolicRealNumberField([4/5])
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = gmic(f, field=K)
        sage: _ = generate_maximal_additive_faces(h);
        sage: list(K.get_eq_list())
        [0]
        sage: list(K.get_eq_poly())
        []
        sage: list(K.get_lt_list())
        [2*f - 2, f - 2, -f, f - 1, -1/f, -f - 1, -2*f, -2*f + 1, -1/(-f^2 + f), -1]
        sage: list(K.get_lt_poly())
        [2*f - 2, f - 2, f^2 - f, -f, f - 1, -f - 1, -2*f, -2*f + 1]
        sage: list(K.get_lt_factor())
        [f - 2, -f + 1/2, f - 1, -f - 1, -f]

        sage: K.<f, lam> = SymbolicRealNumberField([4/5, 1/6])
        sage: h = gj_2_slope(f, lam, field=K)
        sage: list(K.get_eq_list())
        [0]
        sage: list(K.get_eq_poly())
        []
        sage: list(K.get_lt_list())
        [-1/2*f*lam - 1/2*f + 1/2*lam, lam - 1, f - 1, -lam, (-f*lam - f + lam)/(-f + 1), f*lam - lam, (-1/2)/(-1/2*f^2*lam - 1/2*f^2 + f*lam + 1/2*f - 1/2*lam), -f]
        sage: list(K.get_lt_poly())
        [f - 1, -f*lam - f + lam, -1/2*f*lam - 1/2*f + 1/2*lam, 1/2*f^2*lam + 1/2*f^2 - f*lam - 1/2*f + 1/2*lam, -lam, f*lam - lam, -f, lam - 1]
        sage: list(K.get_lt_factor())
        [-lam, f - 1, -f*lam - f + lam, -f, lam - 1]

        sage: K.<f,alpha> = SymbolicRealNumberField([4/5, 3/10])
        sage: h=dg_2_step_mir(f, alpha, field=K, conditioncheck=False)
        sage: extremality_test(h)
        True

        sage: K.<f,a1,a2,a3> = SymbolicRealNumberField([4/5, 1, 3/10, 2/25])
        sage: h = kf_n_step_mir(f, (a1, a2, a3), conditioncheck=False)
        sage: extremality_test(h)
        True

        sage: K.<f> = SymbolicRealNumberField([1/5])
        sage: h = drlm_3_slope_limit(f, conditioncheck=False)
        sage: extremality_test(h)
        True
        sage: list(K.get_lt_factor())
        [f - 1, f - 1/2, f - 1/3, f - 2, -f - 1, -f]

    """

    def __init__(self, values=[], names=()):
        NumberField.__init__(self)
        self._element_class = SymbolicRNFElement
        self._zero_element = SymbolicRNFElement(0, parent=self)
        self._one_element =  SymbolicRNFElement(1, parent=self)
        self._eq = set([])
        self._lt = set([])
        self._eq_poly = set([])
        self._lt_poly = set([])
        self._eq_factor = set([])
        self._lt_factor = set([])
        vnames = PolynomialRing(QQ, names).fraction_field().gens();
        self._gens = [ SymbolicRNFElement(value, name, parent=self) for (value, name) in izip(values, vnames) ]
        self._names = names
        self._values = values

    def __copy__(self):
        Kcopy = self.__class__(self._values, self._names)
        Kcopy._eq.update(self._eq)
        Kcopy._lt.update(self._lt)
        Kcopy._eq_poly.update(self._eq_poly)
        Kcopy._lt_poly.update(self._lt_poly)
        Kcopy._eq_factor.update(self._eq_factor)
        Kcopy._lt_factor.update(self._lt_factor)
        return Kcopy

    def _first_ngens(self, n):
        for i in range(n):
            yield self._gens[i]

    def _an_element_impl(self):
        return SymbolicRNFElement(1, parent=self)
    def _coerce_map_from_(self, S):
        # print "_coerce_map_from: self = %s, S = %s" % (self, S)
        return CallableConvertMap(S, self, lambda s: SymbolicRNFElement(s, parent=self), parent_as_first_arg=False)
    def __repr__(self):
        return 'SymbolicRNF%s' %repr(self.gens())
    def __call__(self, elt):
        if parent(elt) == self:
            return elt
        try: 
            QQ_elt = QQ(elt)
            return SymbolicRNFElement(QQ_elt, parent=self)
        except:
            raise ValueError, "SymbolicRealNumberField called with element", elt

    def _coerce_impl(self, x):
        return self(x)
    def get_eq_list(self):
        return self._eq
    def get_lt_list(self):
        return self._lt
    def get_eq_poly(self):
        return self._eq_poly
    def get_lt_poly(self):
        return self._lt_poly
    def get_eq_factor(self):
        return self._eq_factor
    def get_lt_factor(self):
        return self._lt_factor
    def record_to_eq_list(self, comparison):
        if not comparison in self._eq:
            logging.info("New element in %s._eq: %s" % (repr(self), comparison))
            self._eq.add(comparison)
            self.record_poly(comparison.numerator())
            self.record_poly(comparison.denominator())
    def record_to_lt_list(self, comparison):
        if not comparison in self._lt:
            logging.info("New element in %s._lt: %s" % (repr(self), comparison))
            self._lt.add(comparison)
            self.record_poly(comparison.numerator())
            self.record_poly(comparison.denominator())
    def record_poly(self, poly):
        if not poly in QQ and poly.degree() > 0:
            v = poly(self._values)
            if v == 0:
                self.record_to_eq_poly(poly)
            elif v < 0:
                self.record_to_lt_poly(poly)
            else:
                self.record_to_lt_poly(-poly)
    def record_to_eq_poly(self, poly):
        if not poly in self._eq_poly:
            #logging.info("New element in %s._eq_poly: %s" % (repr(self), poly))
            self._eq_poly.add(poly)
            for (fac, d) in poly.factor():
                # record the factor if it's zero
                if fac(self._values) == 0 and not fac in self._eq_factor:
                    #logging.info("New element in %s._eq_factor: %s" % (repr(self), fac))
                    self._eq_factor.add(fac)
    def record_to_lt_poly(self, poly):
        if not poly in self._lt_poly:
            #logging.info("New element in %s._lt_poly: %s" % (repr(self), poly))
            self._lt_poly.add(poly)
            for (fac, d) in poly.factor():
                # record the factor if it's raised to an odd power.
                if d % 2 == 1:
                    if fac(self._values) < 0:
                        new_fac = fac
                    else:
                        new_fac = -fac
                    if not new_fac in self._lt_factor:
                        #logging.info("New element in %s._lt_factor: %s" % (repr(self), new_fac))
                        self._lt_factor.add(new_fac)

default_symbolic_field = SymbolicRealNumberField()

from sage.libs.ppl import Variable, Constraint, Linear_Expression, Constraint_System, NNC_Polyhedron

def polynomial_to_linexpr(t, monomial_list, v_dict):
    # coefficients in ppl constraint must be integers.
    lcd = lcm([x.denominator() for x in t.coefficients()])
    linexpr = Linear_Expression(0)
    #print 't=%s' % t
    if len(t.args()) <= 1:
        # sage.rings.polynomial.polynomial_rational_flint object has no attribute 'monomials'
        for (k, c) in t.dict().items():
            #print '(k, c) = (%s, %s)' % (k, c);
            m = (t.args()[0])^k
            if m in v_dict.keys():
                v = v_dict[m]
            elif k == 0:
                # constant term, don't construct a new Variable for it.
                v = 1
            else:
                nv = len(monomial_list)
                v = Variable(nv)
                v_dict[m] = v
                monomial_list.append(m)
            #print '(m, v) = (%s, %s)' % (m, v);
            linexpr += (lcd * c) * v
    else:
        for m in t.monomials():
            if m in v_dict.keys():
                v = v_dict[m]
            elif m == 1:
                v = 1
            else:
                nv = len(monomial_list)
                v = Variable(nv)
                v_dict[m] = v
                monomial_list.append(m)
            coeffv = t.monomial_coefficient(m)
            linexpr += (lcd * coeffv) * v
    #print 'linexpr=%s' % linexpr
    return linexpr

def cs_of_eq_lt_poly(eq_poly, lt_poly):
    monomial_list = []
    v_dict ={}
    cs = Constraint_System()
    for t in eq_poly:
        linexpr = polynomial_to_linexpr(t, monomial_list, v_dict)
        cs.insert( linexpr == 0 )
    for t in lt_poly:
        linexpr = polynomial_to_linexpr(t, monomial_list, v_dict)
        cs.insert( linexpr < 0 )
    return cs, monomial_list, v_dict

def simplify_eq_lt_poly_via_ppl(eq_poly, lt_poly):
    """
    Given polymonial equality and inequality lists.
    Treat each monomial as a new variable.
    This gives a linear inequality system.
    Remove redundant inequalities using PPL.

    EXAMPLES::

        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: K.<f> = SymbolicRealNumberField([4/5])
        sage: h = gmic(f, field=K)
        sage: _ = extremality_test(h)
        sage: eq_poly =list(K.get_eq_poly())
        sage: lt_poly =list(K.get_lt_poly())
        sage: (eq_poly, lt_poly)
        ([], [2*f - 2, f - 2, f^2 - f, -2*f, f - 1, -f - 1, -f, -2*f + 1])
        sage: simplify_eq_lt_poly_via_ppl(eq_poly, lt_poly)
        ([], [f - 1, -2*f + 1, f^2 - f])

        sage: eq_factor =list(K.get_eq_factor())
        sage: lt_factor =list(K.get_lt_factor())
        sage: (eq_factor, lt_factor)
        ([], [f - 2, -f + 1/2, f - 1, -f - 1, -f])
        sage: simplify_eq_lt_poly_via_ppl(eq_factor, lt_factor)
        ([], [f - 1, -2*f + 1])

        sage: K.<f, lam> = SymbolicRealNumberField([4/5, 1/6])
        sage: h = gj_2_slope(f, lam, field=K, conditioncheck=False)
        sage: simplify_eq_lt_poly_via_ppl(list(K.get_eq_poly()), list(K.get_lt_poly()))
        ([], [f - 1, f*lam - lam, f^2*lam + f^2 - 2*f*lam - f + lam, -f*lam - f + lam])
        sage: simplify_eq_lt_poly_via_ppl(list(K.get_eq_factor()), list(K.get_lt_factor()))
        ([], [-f*lam - f + lam, f - 1, -lam])

        sage: _ = extremality_test(h)
        sage: simplify_eq_lt_poly_via_ppl(list(K.get_eq_poly()), list(K.get_lt_poly()))
        ([], [f^2*lam + f^2 - 2*f*lam - f + lam, -2*f*lam + f + 2*lam - 1, -f*lam - 3*f + lam + 2, -lam, lam - 1, f*lam - lam])
        sage: simplify_eq_lt_poly_via_ppl(list(K.get_eq_factor()), list(K.get_lt_factor()))
        ([], [f*lam - 3*f - lam + 2, -f*lam - 3*f + lam + 2, 3*f*lam - f - 3*lam, -3*f*lam - f + 3*lam, -lam, f - 1, 2*lam - 1])

        sage: K.<f,alpha> = SymbolicRealNumberField([4/5, 3/10])             # Bad example! parameter region = {given point}.
        sage: h=dg_2_step_mir(f, alpha, field=K, conditioncheck=False)
        sage: _ = extremality_test(h)
        sage: simplify_eq_lt_poly_via_ppl(list(K.get_eq_poly()), list(K.get_lt_poly()))
        ([-10*alpha + 3, -5*f + 4], [5*f^2 - 10*f*alpha - 1])
        sage: simplify_eq_lt_poly_via_ppl(list(K.get_eq_factor()), list(K.get_lt_factor()))
        ([-10*alpha + 3, -5*f + 4], [])

        sage: K.<f> = SymbolicRealNumberField([1/5])
        sage: h = drlm_3_slope_limit(f, conditioncheck=False)
        sage: _ = extremality_test(h)
        sage: simplify_eq_lt_poly_via_ppl(list(K.get_eq_poly()), list(K.get_lt_poly()))
        ([], [3*f - 1, -f^2 - f, -f])
        sage: simplify_eq_lt_poly_via_ppl(list(K.get_eq_factor()), list(K.get_lt_factor()))
        ([], [3*f - 1, -f])

    """
    mineq = []
    minlt = []
    cs, monomial_list, v_dict = cs_of_eq_lt_poly(eq_poly, lt_poly)
    #print 'cs = %s' % cs
    p = NNC_Polyhedron(cs)
    mincs = p.minimized_constraints()
    #print 'mincs = %s' % mincs
    #print 'monomial_list = %s' % monomial_list
    #print 'v_dict = %s' % v_dict
    for c in mincs:
        coeff = c.coefficients()
        # constraint is written with '>', while lt_poly records '<' relation
        t = sum([-x*y for x, y in itertools.izip(coeff, monomial_list)]) - c.inhomogeneous_term()
        if c.is_equality():
            mineq.append(t)
        else:
            minlt.append(t)
    # note that polynomials in mineq and minlt can have leading coefficient != 1
    return mineq, minlt

