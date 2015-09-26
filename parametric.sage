# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

import sage.structure.element
from sage.structure.element import FieldElement

from sage.libs.ppl import Variable, Constraint, Linear_Expression, Constraint_System, NNC_Polyhedron, Poly_Con_Relation
poly_is_included = Poly_Con_Relation.is_included()

###############################
# Symbolic Real Number Field
###############################

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
        if result == 0:
            left.parent().record_to_eq_list(left.sym() - right.sym())
        elif result == -1:
            left.parent().record_to_lt_list(left.sym() - right.sym())
        elif result == 1:
            left.parent().record_to_lt_list(right.sym() - left.sym())
        return result

    def __richcmp__(left, right, op):
        result = left._val.__richcmp__(right._val, op)
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

class SymbolicRealNumberField(Field):
    """
    Parametric search:
    EXAMPLES::

        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: K.<f> = SymbolicRealNumberField([4/5])
        sage: h = gmic(f, field=K)
        sage: _ = generate_maximal_additive_faces(h);
        sage: list(K.get_eq_list()), list(K.get_lt_list())
        ([], [2*f - 2, f - 2, -f, f - 1, -1/f, -f - 1, -2*f, -2*f + 1, -1/(-f^2 + f)])
        sage: list(K.get_eq_poly()), list(K.get_lt_poly())
        ([], [2*f - 2, f - 2, f^2 - f, -f, f - 1, -f - 1, -2*f, -2*f + 1])
        sage: list(K.get_eq_factor()), list(K.get_lt_factor())
        ([], [-f + 1/2, f - 1, -f])

        sage: K.<f, lam> = SymbolicRealNumberField([4/5, 1/6])
        sage: h = gj_2_slope(f, lam, field=K)
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
        [f - 1/2, f - 1/3, f - 1, -f]
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

        # do the computation of the polyhedron incrementally,
        # rather than first building a huge list and then in a second step processing it.
        # the polyhedron defined by all constraints in self._eq/lt_factor
        self.polyhedron = NNC_Polyhedron(0, 'universe')
        # records the monomials that appear in self._eq/lt_factor
        self.monomial_list = []
        # a dictionary that maps each monomial to the index of its corresponding Variable in self.polyhedron
        self.v_dict = {}
        # record QQ_linearly_independent of pairs. needs simplification
        self._independent_pairs = set([])
        self._dependency= []
        self._independency=[]
        self._zero_kernel=set([])

    def __copy__(self):
        logging.warn("copy(%s) is invoked" % self)
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
    def ngens(self):
        return len(self._gens)
    def _an_element_impl(self):
        return SymbolicRNFElement(1, parent=self)
    def _coerce_map_from_(self, S):
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
        if not comparison.is_zero() and not comparison in QQ and not comparison in self._eq:
            logging.debug("New element in %s._eq: %s" % (repr(self), comparison))
            self._eq.add(comparison)
            self.record_poly(comparison.numerator())
            self.record_poly(comparison.denominator())
    def record_to_lt_list(self, comparison):
        if not comparison in QQ and not comparison in self._lt:
            logging.debug("New element in %s._lt: %s" % (repr(self), comparison))
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
            self._eq_poly.add(poly)
            for (fac, d) in poly.factor():
                # record the factor if it's zero
                if fac(self._values) == 0 and not fac in self._eq_factor:
                    self.record_factor(fac, operator.eq)

    def record_to_lt_poly(self, poly):
        if not poly in self._lt_poly:
            self._lt_poly.add(poly)
            for (fac, d) in poly.factor():
                # record the factor if it's raised to an odd power.
                if d % 2 == 1:
                    if fac(self._values) < 0:
                        new_fac = fac
                    else:
                        new_fac = -fac
                    if not new_fac in self._lt_factor:
                        self.record_factor(new_fac, operator.lt)

    def record_factor(self, fac, op):
        #print "add %s, %s to %s" % (fac, op, self.polyhedron.constraints())
        space_dim_old = len(self.monomial_list)
        linexpr = polynomial_to_linexpr(fac, self.monomial_list, self.v_dict)
        space_dim_to_add = len(self.monomial_list) - space_dim_old
        if op == operator.lt:
            constraint_to_add = (linexpr < 0)
        else:
            constraint_to_add = (linexpr == 0)
        #print "constraint_to_add = %s" % constraint_to_add
        if space_dim_to_add:
            self.polyhedron.add_space_dimensions_and_embed(space_dim_to_add)
            add_new_element = True
        else:
            add_new_element = not self.polyhedron.relation_with(constraint_to_add).implies(poly_is_included)
        if add_new_element:
            self.polyhedron.add_constraint(constraint_to_add)
            #print " add new constraint, %s" %self.polyhedron.constraints()
            if op == operator.lt:
                logging.info("New constraint: %s < 0" % fac)
                self._lt_factor.add(fac)
            else:
                logging.info("New constraint: %s == 0" % fac)
                self._eq_factor.add(fac)

    def record_independence_of_pair(self, numbers, is_independent):
        if len(numbers) != 2:
            raise NotImplementedError, "%s has more than two elements. Not implemented." % numbers
        t1 = affine_linear_form_of_symbolicrnfelement(numbers[0])
        t2 = affine_linear_form_of_symbolicrnfelement(numbers[1])
        vector_space = VectorSpace(QQ,len(t1))
        t1 = vector_space(t1)
        t2 = vector_space(t2)
        pair_space = vector_space.subspace([t1, t2])
        if pair_space.dimension() <= 1:
            if is_independent:
                raise ValueError, "Contradiction: (%s, %s) are not linearly independent in Q." % (t1, t2)
        else:
            if is_independent:
                self._independent_pairs.add(pair_space)
                self._zero_kernel.add(vector_space.subspace([t1]).gen(0))
                self._zero_kernel.add(vector_space.subspace([t2]).gen(0))
            else:
                self._dependency = update_dependency(self._dependency, pair_space)

    def construct_independency(self):
        self._independency, self._zero_kernel = construct_independency(self._independent_pairs, self._dependency, self._zero_kernel)

    def get_reduced_independent_pairs(self):
        if not self._independency:
            self._independency, self._zero_kernel = construct_independency(\
                self._independent_pairs, self._dependency, self._zero_kernel)
        reduced_independent_pairs = get_independent_pairs_from_independency(self._independency)
        return reduced_independent_pairs

default_symbolic_field = SymbolicRealNumberField()


###############################
# Simplify polynomials
###############################

def polynomial_to_linexpr(t, monomial_list, v_dict):
    """
    sage: P.<x,y,z> = QQ[]
    sage: monomial_list = []; v_dict = {};
    sage: t = 27/113 * x^2 + y*z + 1/2
    sage: polynomial_to_linexpr(t, monomial_list, v_dict)
    54*x0+226*x1+113
    sage: monomial_list
    [x^2, y*z]
    sage: v_dict
    {y*z: x1, x^2: x0}

    sage: tt = x + 1/3 * y*z
    sage: polynomial_to_linexpr(tt, monomial_list, v_dict)
    x1+3*x2
    sage: monomial_list
    [x^2, y*z, x]
    sage: v_dict
    {x: x2, y*z: x1, x^2: x0}
    """
    # coefficients in ppl constraint must be integers.
    lcd = lcm([x.denominator() for x in t.coefficients()])
    linexpr = Linear_Expression(0)
    if len(t.args()) <= 1:
        # sage.rings.polynomial.polynomial_rational_flint object has no attribute 'monomials'
        for (k, c) in t.dict().items():
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
    return linexpr

def cs_of_eq_lt_poly(eq_poly, lt_poly):
    """
    sage: P.<f>=QQ[]
    sage: eq_poly =[]; lt_poly = [2*f - 2, f - 2, f^2 - f, -2*f, f - 1, -f - 1, -f, -2*f + 1]
    sage: cs, monomial_list, v_dict = cs_of_eq_lt_poly(eq_poly, lt_poly)
    sage: cs
    Constraint_System {-x0+1>0, -x0+2>0, x0-x1>0, x0>0, -x0+1>0, x0+1>0, x0>0, 2*x0-1>0}
    sage: monomial_list
    [f, f^2]
    sage: v_dict
    {f: x0, f^2: x1}
    """
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
        ([], [-f + 1/2, f - 1, -f])
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
        ([],
         [-3*f*lam - f + 3*lam,
          3*f*lam - f - 3*lam,
          -f*lam - 3*f + lam + 2,
          f*lam - 3*f - lam + 2,
          f - 1,
          2*lam - 1,
          -lam])

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
    cs, monomial_list, v_dict = cs_of_eq_lt_poly(eq_poly, lt_poly)
    p = NNC_Polyhedron(cs)
    return read_leq_lin_from_polyhedron(p, monomial_list, v_dict)


def read_leq_lin_from_polyhedron(p, monomial_list, v_dict):
    """
    sage: P.<f>=QQ[]
    sage: eq_poly =[]; lt_poly = [2*f - 2, f - 2, f^2 - f, -2*f, f - 1, -f - 1, -f, -2*f + 1]
    sage: cs, monomial_list, v_dict = cs_of_eq_lt_poly(eq_poly, lt_poly)
    sage: p = NNC_Polyhedron(cs)
    sage: read_leq_lin_from_polyhedron(p, monomial_list, v_dict)
    ([], [f - 1, -2*f + 1, f^2 - f])
    """
    mineq = []
    minlt = []
    mincs = p.minimized_constraints()
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

def read_simplified_leq_lin(K, level="factor"):
    """
    sage: K.<f> = SymbolicRealNumberField([4/5])
    sage: h = gmic(f, field=K)
    sage: _ = extremality_test(h)
    sage: read_simplified_leq_lin(K)
    ([], [f - 1, -2*f + 1])

    sage: K.<f> = SymbolicRealNumberField([1/5])
    sage: h = drlm_3_slope_limit(f, conditioncheck=False)
    sage: _ = extremality_test(h)
    sage: read_simplified_leq_lin(K)
    ([], [3*f - 1, -f])
    """
    if level == "factor":
        #leq, lin = simplify_eq_lt_poly_via_ppl(K.get_eq_factor(), K.get_lt_factor())
        # Since we update K.polyhedron incrementally,
        # just read leq and lin from its minimized constraint system.
        leq, lin = read_leq_lin_from_polyhedron(K.polyhedron, K.monomial_list, K.v_dict)
    elif level == "poly":
        leq, lin = simplify_eq_lt_poly_via_ppl(K.get_eq_poly(), K.get_lt_poly())
    else:
        leq = list(K.get_eq_list())
        lin = list(K.get_lt_list())
    if leq:
        logging.warn("equation list %s is not empty!" % leq)
    return leq, lin


######################################
# Extreme functions with the magic K
######################################

from sage.misc.sageinspect import sage_getargspec, sage_getvariablename

def read_default_args(function, **opt_non_default):
    """
    sage: read_default_args(gmic)
    {'conditioncheck': True, 'f': 4/5, 'field': None}
    sage: read_default_args(drlm_backward_3_slope, **{'bkpt': 1/5})
    {'bkpt': 1/5, 'conditioncheck': True, 'f': 1/12, 'field': None}
    """
    args, varargs, keywords, defaults = sage_getargspec(function)
    default_args = {}
    for i in range(len(defaults)):
        default_args[args[-i-1]]=defaults[-i-1]
    for (opt_name, opt_value) in opt_non_default.items():
        if opt_name in default_args:
            default_args[opt_name] = opt_value
    return default_args

def construct_field_and_test_point(function, var_name, var_value, default_args):
    """
    sage: function=gmic; var_name=['f']; var_value=[1/2];
    sage: default_args = read_default_args(function)
    sage: K, test_point = construct_field_and_test_point(function, var_name, var_value, default_args)
    sage: K
    SymbolicRNF[f~]
    sage: test_point
    {'conditioncheck': False, 'f': f~, 'field': SymbolicRNF[f~]}
    """
    K = SymbolicRealNumberField(var_value, var_name)
    test_point = copy(default_args)
    for i in range(len(var_name)):
        test_point[var_name[i]] = K.gens()[i]
    test_point['field'] = K
    test_point['conditioncheck'] = False
    return K, test_point

def simplified_extremality_test(function):
    """
    function has rational bkpts; function is known to be minimal.
    """
    f = find_f(function, no_error_if_not_minimal_anyway=True)
    covered_intervals = generate_covered_intervals(function)
    uncovered_intervals = generate_uncovered_intervals(function)
    if uncovered_intervals:
        logging.info("Function has uncovered intervals, thus is NOT extreme.")
        return False
    else:
        components = covered_intervals
    field = function(0).parent().fraction_field()
    symbolic = generate_symbolic(function, components, field=field)
    equation_matrix = generate_additivity_equations(function, symbolic, field, f=f)
    slope_jump_vects = equation_matrix.right_kernel_matrix()
    sol_dim = slope_jump_vects.nrows()
    if sol_dim > 0:
        logging.info("Finite dimensional test: Solution space has dimension %s" % sol_dim)
        logging.info("Thus the function is NOT extreme.")
        return False
    else:
        logging.info("The function is extreme.")
        return True

##############################
# Find regions
#############################

def find_region_according_to_literature(function, var_name, var_value, region_level, default_args=None, \
                                        level="factor", non_algrebraic_warn=True, **opt_non_default):
    """
    sage: find_region_according_to_literature(gmic, ['f'], [4/5], 'constructible')
    ([], [f - 1, -f])
    sage: find_region_according_to_literature(drlm_backward_3_slope, ['f','bkpt'], [1/12, 2/12], 'extreme')
    ([], [f - bkpt, -f + 4*bkpt - 1, -f])
    """
    ## Note: region_level = 'constructible' / 'extreme'
    if default_args is None:
        default_args = read_default_args(function, **opt_non_default)
    if not isinstance(function, ExtremeFunctionsFactory):
        K, test_point = construct_field_and_test_point(function, var_name, var_value, default_args)
        if region_level == 'extreme':
            test_point['conditioncheck'] = True
        h = function(**test_point)
        leq, lin =  read_simplified_leq_lin(K, level=level)
        if leq:
            logging.warn("Bad default arguments!")
        return leq, lin
    else:
        if non_algrebraic_warn:
            logging.warn("This function involves non-algebraic operations. Parameter region according to literature was plotted separately.")
        test_point = copy(default_args)
        del test_point['field']
        del test_point['conditioncheck']
        for v in var_name:
            del test_point[v]
        x, y = var('x, y')
        if region_level == 'constructible':
            if len(var_name) == 2:
                return lambda x, y: function.claimed_parameter_attributes(**dict({var_name[0]: x, var_name[1]: y}, **test_point)) != 'not_constructible'
            else:
                return lambda x, y: function.claimed_parameter_attributes(**dict({var_name[0]: x}, **test_point )) != 'not_constructible'
        elif region_level == 'extreme':
            if len(var_name) == 2:
                return lambda x, y: function.claimed_parameter_attributes(**dict({var_name[0]: x, var_name[1]: y}, **test_point)) == 'extreme'
            else:
                return lambda x, y: function.claimed_parameter_attributes(**dict({var_name[0]: x}, **test_point )) == 'extreme'
        else:
            raise ValueError, "Bad argument region_level = %s" % region_level

def find_region_around_given_point(K, h, level="factor", region_level='extreme', is_minimal=None, use_simplified_extremality_test=True):
    ## Note: region_level = 'constructible' / 'minimal'/ 'extreme'. test cases see find_parameter_region()
    if region_level == 'constructible':
        leq, lin = read_simplified_leq_lin(K, level=level)
        return 'is_constructible', leq, lin
    if is_minimal is None:
        is_minimal = minimality_test(h)
    if is_minimal:
        if region_level == 'minimal':
            leq, lin = read_simplified_leq_lin(K, level=level)
            return 'is_minimal', leq, lin
        if use_simplified_extremality_test:
            is_extreme = simplified_extremality_test(h)
        else:
            is_extreme = extremality_test(h)
        leq, lin = read_simplified_leq_lin(K, level=level)
        if is_extreme:
            return 'is_extreme', leq, lin
        else:
            return 'not_extreme', leq, lin
    else:
        leq, lin = read_simplified_leq_lin(K, level=level)
        return 'not_minimal', leq, lin

def find_parameter_region(function=drlm_backward_3_slope, var_name=['f'], var_value=None, use_simplified_extremality_test=True,\
                        level="factor", region_level='extreme', **opt_non_default):
    """
    sage: find_parameter_region(drlm_backward_3_slope, ['f','bkpt'], [1/12+1/30, 2/12], region_level='constructible')
    ('is_constructible', [], [f - bkpt, -f + 2*bkpt - 1, -f])
    sage: find_parameter_region(drlm_backward_3_slope, ['f','bkpt'], [1/12+1/30, 2/12], region_level='minimal')
    ('is_minimal', [], [f - bkpt, -f + 3*bkpt - 1, -2*f + bkpt])
    sage: find_parameter_region(drlm_backward_3_slope, ['f','bkpt'], [1/12+1/30, 2/12], region_level='extreme')
    ('is_extreme', [], [f - bkpt, -f + 4*bkpt - 1, -2*f + bkpt, f^2 - f*bkpt - 1])
    """
    # similar to find_region_around_given_point(), but start from constructing K and test point.
    default_args = read_default_args(function, **opt_non_default)
    if not var_value:
        var_value = [default_args[v] for v in var_name]
    K, test_point = construct_field_and_test_point(function, var_name, var_value, default_args)
    try:
        h = function(**test_point)
        region_type, leq, lin = find_region_around_given_point(K, h, level=level, region_level=region_level, \
                        is_minimal=None, use_simplified_extremality_test=use_simplified_extremality_test)
    except:
        region_type, leq, lin = 'not_constructible', [], []
    return region_type, leq, lin

##############################
# Plot regions
#############################

from sage.plot.contour_plot import equify
from sage.plot.contour_plot import ContourPlot
from sage.plot.primitive import GraphicPrimitive
from sage.misc.decorators import options, suboptions
from sage.plot.colors import rgbcolor, get_cmap
from sage.misc.misc import xsrange
import operator

@options(plot_points=100, incol='blue', outcol=None, bordercol=None, borderstyle=None, borderwidth=None,frame=False,axes=True, legend_label=None, aspect_ratio=1, alpha=1)
def region_plot_patch(f, xrange, yrange, plot_points, incol, outcol, bordercol, borderstyle, borderwidth, alpha, **options):
    from sage.plot.all import Graphics
    from sage.plot.misc import setup_for_eval_on_grid
    from sage.symbolic.expression import is_Expression
    import numpy

    if not isinstance(f, (list, tuple)):
        f = [f]

    feqs = [equify(g) for g in f if is_Expression(g) and g.operator() is operator.eq and not equify(g).is_zero()]
    f = [equify(g) for g in f if not (is_Expression(g) and g.operator() is operator.eq)]
    neqs = len(feqs)
    if neqs > 1:
        logging.warn("There are at least 2 equations; If the region is degenerated to points, plotting might show nothing.")
        feqs = [sum([fn**2 for fn in feqs])]
        neqs = 1
    if neqs and not bordercol:
        bordercol = incol
    if not f:
        return implicit_plot(feqs[0], xrange, yrange, plot_points=plot_points, fill=False, \
                             linewidth=borderwidth, linestyle=borderstyle, color=bordercol, **options)
    f_all, ranges = setup_for_eval_on_grid(feqs + f, [xrange, yrange], plot_points)
    xrange,yrange=[r[:2] for r in ranges]

    xy_data_arrays = numpy.asarray([[[func(x, y) for x in xsrange(*ranges[0], include_endpoint=True)]
                                     for y in xsrange(*ranges[1], include_endpoint=True)]
                                    for func in f_all[neqs::]],dtype=float)
    xy_data_array=numpy.abs(xy_data_arrays.prod(axis=0))
    # Now we need to set entries to negative iff all
    # functions were negative at that point.
    neg_indices = (xy_data_arrays<0).all(axis=0)
    xy_data_array[neg_indices]=-xy_data_array[neg_indices]

    from matplotlib.colors import ListedColormap
    incol = rgbcolor(incol)
    if outcol:
        outcol = rgbcolor(outcol)
        cmap = ListedColormap([incol, outcol])
        cmap.set_over(outcol, alpha=alpha)
    else:
        outcol = rgbcolor('white')
        cmap = ListedColormap([incol, outcol])
        cmap.set_over(outcol, alpha=0)
    cmap.set_under(incol, alpha=alpha)

    g = Graphics()

    # Reset aspect_ratio to 'automatic' in case scale is 'semilog[xy]'.
    # Otherwise matplotlib complains.
    scale = options.get('scale', None)
    if isinstance(scale, (list, tuple)):
        scale = scale[0]
    if scale == 'semilogy' or scale == 'semilogx':
        options['aspect_ratio'] = 'automatic'

    g._set_extra_kwds(Graphics._extract_kwds_for_show(options, ignore=['xmin', 'xmax']))

    if neqs == 0:
        g.add_primitive(ContourPlot(xy_data_array, xrange,yrange,
                                dict(contours=[-1e-20, 0, 1e-20], cmap=cmap, fill=True, **options)))
    else:
        mask = numpy.asarray([[elt > 0 for elt in rows] for rows in xy_data_array], dtype=bool)
        xy_data_array = numpy.asarray([[f_all[0](x, y) for x in xsrange(*ranges[0], include_endpoint=True)]
                                        for y in xsrange(*ranges[1], include_endpoint=True)], dtype=float)
        xy_data_array[mask] = None
    if bordercol or borderstyle or borderwidth:
        cmap = [rgbcolor(bordercol)] if bordercol else ['black']
        linestyles = [borderstyle] if borderstyle else None
        linewidths = [borderwidth] if borderwidth else None
        g.add_primitive(ContourPlot(xy_data_array, xrange, yrange,
                                    dict(linestyles=linestyles, linewidths=linewidths,
                                         contours=[0], cmap=[bordercol], fill=False, **options)))

    return g

def plot_region_given_dim_and_description(d, region_description, color, alpha=0.5, xmin=-0.1, xmax=1.1, ymin=-0.1, ymax=1.1, plot_points=1000):
    x, y = var('x, y')
    if type(region_description) is tuple:
        leq, lin = region_description
        if d == 2:
            return region_plot_patch([ lhs(x, y) == 0 for lhs in leq ] + [ lhs(x, y) < 0 for lhs in lin ], \
                    (x, xmin, xmax), (y, ymin, ymax), incol=color, alpha=alpha, plot_points=plot_points, bordercol=color)
        else:
            return region_plot_patch([ lhs(x) == 0 for lhs in leq ] + [ lhs(x) < 0 for lhs in lin ] + [y >= -0.01, y <= 0.01], \
                    (x, xmin, xmax), (y, -0.1, 0.3), incol=color, alpha=alpha, plot_points=plot_points, bordercol=color, ticks=[None,[]])
    else:
        # region_description is a lambda function,
        # such as find_region_according_to_literature(function from factory classes)
        if d == 2:
            return region_plot_patch([region_description], (x, xmin, xmax), (y, ymin, ymax), \
                    incol=color, alpha=alpha, plot_points=plot_points, bordercol=color)
        else:
            return region_plot_patch([region_description] + [y >= -0.01, y <= 0.01], (x, xmin, xmax), (y, -0.1, 0.3), \
                    incol=color, alpha=alpha, plot_points=plot_points, bordercol=color, ticks=[None,[]])

def plot_region_according_to_literature(function, var_name, default_args=None, level="factor", \
                                        alpha=0.5, xmin=-0.1, xmax=1.1, ymin=-0.1, ymax=1.1, plot_points=1000, **opt_non_default):
    """
    sage: plot_region_according_to_literature(drlm_backward_3_slope, ['f','bkpt'], plot_points=1000) # not tested
    sage: plot_region_according_to_literature(gj_forward_3_slope, ['lambda_1', 'lambda_2'], f = 2/3, plot_points=500)  # not tested
    sage: plot_region_according_to_literature(dg_2_step_mir, ['f','alpha'], plot_points=500)  # not tested
    """
    if default_args is None:
        default_args = read_default_args(function, **opt_non_default)
    var_value = [default_args[v] for v in var_name]
    g = Graphics()
    d = len(var_name)
    region_constructible = find_region_according_to_literature(function, var_name, var_value, 'constructible', default_args, level=level)
    g += plot_region_given_dim_and_description(d, region_constructible, "orange", \
                    alpha=alpha, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, plot_points=plot_points)
    g += line([(0,0),(0,0.01)], color = "orange", legend_label="constructible", zorder=-10)
    region_extreme = find_region_according_to_literature(function, var_name, var_value, 'extreme', default_args, level=level, non_algrebraic_warn=False)
    g += plot_region_given_dim_and_description(d, region_extreme, "red", \
                    alpha=alpha, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, plot_points=plot_points)
    g += line([(0,0),(0,0.01)], color ="red", legend_label="claimed_extreme", zorder=-10)
    return g

def plot_covered_regions(regions, tested_points, d, alpha=0.5, xmin=-0.1, xmax=1.1, ymin=-0.1, ymax=1.1, plot_points=1000, show_points=2):
    covered_type_color = {'not_minimal': 'orange', 'not_extreme': 'green', 'is_extreme': 'blue'}
    g = Graphics()
    for (covered_type, leq, lin) in regions:
        g += plot_region_given_dim_and_description(d, (leq, lin), color=covered_type_color[covered_type], \
                            alpha=alpha, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, plot_points=plot_points)
    # show_points: 0 -- don't show points, 1 -- show white points; 2 -- show all points
    if show_points:
        for (p_val, p_col) in tested_points:
            if show_points == 1 and p_col == 'black':
                # don't plot black points out of the constructible region.
                continue
            if d == 2:
                g += point([p_val], color = p_col, size = 2, zorder=10)
            else:
                g += point([(p_val[0], 0)], color = p_col, size = 2, zorder=10)
    return g

def plot_parameter_region(function=drlm_backward_3_slope, var_name=['f'], var_value=None, use_simplified_extremality_test=True,\
                        alpha=0.5, xmin=-0.1, xmax=1.1, ymin=-0.1, ymax=1.1, level="factor", plot_points=1000, **opt_non_default):

    """
    sage: logging.disable(logging.INFO) # not tested
    sage: g, fname, fparam = plot_parameter_region(drlm_backward_3_slope, ['f','bkpt'], [1/12-1/30, 2/12]) # not tested
    sage: g.save(fname+".pdf", title=fname+fparam, legend_loc=7) # not tested

    sage: plot_parameter_region(drlm_backward_3_slope, ['f','bkpt'], [1/12+1/30, 2/12], plot_points=500) # not tested
    sage: plot_parameter_region(drlm_backward_3_slope, ['f','bkpt'], [1/10-1/100, 3/10], plot_points=500) # not tested

    sage: plot_parameter_region(gj_2_slope, ['f','lambda_1'], plot_points=100) # not tested
    sage: plot_parameter_region(gj_2_slope, ['f','lambda_1'], [4/5, 4/6], plot_points=100) # not tested

    sage: plot_parameter_region(gj_forward_3_slope, ['f','lambda_1'], [9/10, 3/10]) # not tested
    sage: plot_parameter_region(gj_forward_3_slope, ['lambda_1','lambda_2']) # not tested
    sage: plot_parameter_region(gj_forward_3_slope, ['lambda_1','lambda_2'], f=3/4, lambda_1=2/3, lambda_2=3/5) # not tested
    sage: plot_parameter_region(gj_forward_3_slope, ['lambda_1','lambda_2'], [3/5, 7/5], f=1/2, lambda_1=4/9, lambda_2=1/3) # not tested
    sage: plot_parameter_region(gj_forward_3_slope, ['lambda_1','lambda_2'], [9/10, 3/7], f=1/2, lambda_1=4/9, lambda_2=1/3) # not tested

    sage: plot_parameter_region(chen_4_slope, ['lam1', 'lam2'], [1/4, 1/4], plot_points=100) # not tested

    sage: plot_parameter_region(dg_2_step_mir, ['f','alpha'], [3/5,2/5-1/23],  plot_points=500) # not tested
    sage: plot_parameter_region(dg_2_step_mir_limit, ['f'], [8/13], d=2, plot_points=500) # not tested

    sage: plot_parameter_region(kf_n_step_mir, ['f'], [1/13], plot_points=100) # not tested
    """
    d = len(var_name)
    if d >= 3:
        #TODO
        raise NotImplementedError, "More than three parameters. Not implemented."
    region_type_color = {'is_minimal': 'green', 'not_minimal': 'darkslategrey', 'is_extreme': 'blue', 'not_extreme': 'cyan'}
    default_args = read_default_args(function, **opt_non_default)
    g = plot_region_according_to_literature(function, var_name, default_args, level=level, \
                                alpha=alpha, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, plot_points=plot_points)
    if not var_value:
        var_value = [default_args[v] for v in var_name]
    K, test_point = construct_field_and_test_point(function, var_name, var_value, default_args)
    try:
        h = function(**test_point)
        color_p = "white"
        region_type, leq, lin = find_region_around_given_point(K, h, level=level, region_level='minimal', \
                        is_minimal=None, use_simplified_extremality_test=use_simplified_extremality_test)
        g += plot_region_given_dim_and_description(d, (leq, lin), color=region_type_color[region_type], alpha=alpha, \
                        xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, plot_points=plot_points)
        g += line([(0,0),(0,0.01)], color=region_type_color[region_type], legend_label=region_type, zorder=-10)
        if region_type == 'is_minimal':
            region_type, leq, lin = find_region_around_given_point(K, h, level=level, region_level='extreme', \
                            is_minimal=True, use_simplified_extremality_test=use_simplified_extremality_test)
            g += plot_region_given_dim_and_description(d, (leq, lin), color=region_type_color[region_type], alpha=alpha, \
                            xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, plot_points=plot_points)
            g += line([(0,0),(0,0.01)], color=region_type_color[region_type], legend_label=region_type, zorder=-10)
    except:
        logging.warn("%s = %s is out of the constructible region." % (var_name, var_value))
        color_p = "black"
    if d == 2:
        p = point([var_value], color = color_p, size = 10, zorder=10)
        fun_param = "(%s=%s, %s=%s)" % (var_name[0], var_value[0], var_name[1], var_value[1])
    else:
        p = point([(var_value[0], 0)], color = color_p, size = 10, zorder=10)
        fun_param = "(%s=%s)" % (var_name[0], var_value[0])
    g += p

    fun_names = sage_getvariablename(function)
    i = 0
    while fun_names[i] == 'function':
        i += 1
    fun_name = fun_names[i]

    g.show(title=fun_name+fun_param) #show_legend=False, axes_labels=var_name)
    return g, fun_name, fun_param

###########################
# Randomly pick some points.
# Can their parameter regions cover the extreme region according to literature?
###########################

def cover_parameter_region(function=drlm_backward_3_slope, var_name=['f'], random_points=10, xmin=-0.1, xmax=1.1, ymin=-0.1, ymax=1.1, \
                            plot_points=0, show_points=2, regions=None, tested_points=None, **opt_non_default):
    """
    sage: regions, tested_points, last_n_points_was_covered = cover_parameter_region(gmic, ['f'], random_points=10, plot_points=100) # not tested
    sage: _ = cover_parameter_region(drlm_backward_3_slope, ['f','bkpt'], random_points=200, plot_points=500) # not tested
    sage: _ = cover_parameter_region(gj_2_slope, ['f', 'lambda_1'], 300, plot_points=100, show_points=1) # not tested

    sage: _ = cover_parameter_region(gj_forward_3_slope, ['lambda_1', 'lambda_2'], random_points=500, plot_points=500) # not tested
    sage: _ = cover_parameter_region(gj_forward_3_slope, ['f', 'lambda_1'], random_points=100, plot_points=500) # not tested

    sage: regions, tested_points, last_n_points_was_covered = cover_parameter_region(chen_4_slope, ['lam1', 'lam2'], 20, plot_points=0) # not tested
    sage: _ = cover_parameter_region(chen_4_slope, ['lam1', 'lam2'], 20, plot_points=500, show_points=1, regions=regions, tested_points=tested_points) # not tested
    """
    d = len(var_name)
    if d >= 3:
        #TODO
        raise NotImplementedError, "More than three parameters. Not implemented."
    default_args = read_default_args(function, **opt_non_default)
    if regions is None:
        regions = []
    if tested_points is None:
        tested_points = []
    last_n_points_were_covered = 0
    while random_points > 0:
        x, y = QQ(uniform(xmin, xmax)), QQ(uniform(ymin, ymax))
        if d == 2:
            var_value = [x, y]
        else:
            var_value = [x]
        # NOTE: When the function is not_constructible on this point,
        # it's costly to go through all the previously found regions and then return "point_is_covered=False"
        # So, first test if it's constructible
        K, test_point = construct_field_and_test_point(function, var_name, var_value, default_args)
        try:
            h = function(**test_point)
            # Function is contructible at this random point.
            tested_points.append((var_value, 'white'))
            random_points -= 1
            # Test if this point is already covered.
            point_is_covered = False
            for (region_type, leq, lin) in regions:
                if (d == 2) and all(l(x, y) == 0 for l in leq) and all(l(x, y) < 0 for l in lin) or \
                   (d == 1) and all(l(x) == 0 for l in leq) and all(l(x) < 0 for l in lin):
                    point_is_covered = True
                    break
            if point_is_covered:
                last_n_points_were_covered += 1
            else:
                region_type, leq, lin =  find_region_around_given_point(K, h, level="factor", region_level='extreme', \
                                                                is_minimal=None,use_simplified_extremality_test=True)
                regions.append((region_type, leq, lin))
                last_n_points_were_covered = 0
        except:
            # Function is not constructible at this random point.
            tested_points.append((var_value, 'black'))
    if plot_points > 0:
        g = plot_covered_regions(regions, tested_points, d=d, alpha=0.5, \
            xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, plot_points=plot_points, show_points=show_points)
        fun_names = sage_getvariablename(function)
        i = 0
        while fun_names[i] == 'function':
            i += 1
        fun_name = fun_names[i]
        g.show(title=fun_name+'%s' % var_name)
    return regions, tested_points, last_n_points_were_covered

##############################
# linearly independent in Q
##############################
def affine_linear_form_of_symbolicrnfelement(number):
    sym_frac_poly = number.sym()
    if sym_frac_poly.denominator() != 1:
        raise NotImplementedError, "%s is a fraction of polynomials. Not implemented." % sym_frac_poly
    poly = sym_frac_poly.numerator()
    coef = [poly.monomial_coefficient(v) for v in poly.parent().gens()]
    return coef + [poly.constant_coefficient()]

def update_dependency(dependency, parallel_space):
    """
    Update the dependency list given a new parallel relation f // g.

    Inputs:

        - `dependency` is a list [E1, E2, ..., En],
        where Ei is a subspace over Q of affine linear forms in \Q^{k+1},
        such that for any l1, l2 in Ei, it is known that l1 // l2.
        (For any li in Ei, lj in Ej, it is unknown whether l1 // l2.)

        - `parallel_space` is the vector space over Q generated by f and g, f // g.


    EXAMPLES::

        sage: V = VectorSpace(QQ,5)
        sage: E1 = V.subspace([V([1,0,0,0,0]), V([0,0,1,0,0])])
        sage: E2 = V.subspace([V([0,1,0,0,0]), V([0,0,0,1,0])])
        sage: pair1 = V.subspace([V([1,1,0,0,0]), V([0,0,0,0,1])])
        sage: dependency = update_dependency([E1, E2], pair1)
        sage: dependency
        [Vector space of degree 5 and dimension 2 over Rational Field
         Basis matrix:
         [1 0 0 0 0]
         [0 0 1 0 0], Vector space of degree 5 and dimension 2 over Rational Field
         Basis matrix:
         [0 1 0 0 0]
         [0 0 0 1 0], Vector space of degree 5 and dimension 2 over Rational Field
         Basis matrix:
         [1 1 0 0 0]
         [0 0 0 0 1]]
        sage: pair2 = V.subspace([V([1,0,0,0,0]), V([0,1,0,0,0])])
        sage: update_dependency(dependency, pair2)
        [Vector space of degree 5 and dimension 5 over Rational Field
        Basis matrix:
         [1 0 0 0 0]
         [0 1 0 0 0]
         [0 0 1 0 0]
         [0 0 0 1 0]
         [0 0 0 0 1]]
        sage: update_dependency([E1, E2], pair2)
        [Vector space of degree 5 and dimension 4 over Rational Field
         Basis matrix:
         [1 0 0 0 0]
         [0 1 0 0 0]
         [0 0 1 0 0]
         [0 0 0 1 0]]

    Update dependency sequentially:

        sage: V = VectorSpace(QQ,7)
        sage: m = matrix.identity(QQ, 7)
        sage: e = [V(m[i]) for i in range(7)]
        sage: dep_pairs = [V.subspace([e[0]+e[3], 3*e[0]]), \
        ....: V.subspace([e[6], 2*e[2]-3*e[1]]), \
        ....: V.subspace([2*e[6], e[4]-3*e[1]])]
        sage: dependency = []
        sage: for dep_pair in dep_pairs:
        ....:     dependency = update_dependency(dependency, dep_pair)
        sage: dependency
        [Vector space of degree 7 and dimension 2 over Rational Field
         Basis matrix:
         [1 0 0 0 0 0 0]
         [0 0 0 1 0 0 0], Vector space of degree 7 and dimension 3 over Rational Field
         Basis matrix:
         [   0    1    0    0 -1/3    0    0]
         [   0    0    1    0 -1/2    0    0]
         [   0    0    0    0    0    0    1]]
    """
    parallel = parallel_space
    new_dependency = []
    s = None
    while len(new_dependency) < len(dependency):
        if not s is None:
            dependency = new_dependency
            new_dependency = []
        for s in dependency:
            if s.intersection(parallel).dimension() > 0:
                parallel += s
            else:
                new_dependency.append(s)
    new_dependency.append(parallel)
    return new_dependency

def update_independency(independency, simultaneously_independency):
    """
    `independency` is a list [S1, S2, ...]
    Si is a set of pairs (X1, Y1), (X2, Y2) such that Xj _|_ Yj simultaneously.
    Xj, Yj are either some Ek (parallel subspace) or an affine linear form.
    Ei _|_ Ej means for any f in Ei, g in Ej such that ker f = ker g = {0},
    we have f _|_ g.
    """
    sim_ind = simultaneously_independency
    new_independency = []
    for s in independency:
        if s.intersection(sim_ind):
            sim_ind.update(s)
        else:
            new_independency.append(s)
    new_independency.append(sim_ind)
    return new_independency

def construct_independency(independent_pairs, dependency, seen_linear_forms):
    """
    EXAMPLES:

        sage: V = VectorSpace(QQ,7)
        sage: m = matrix.identity(QQ, 7)
        sage: e = [V(m[i]) for i in range(7)]

    Imagine e[0] =1, e[1] = sqrt(2), e[2] = sqrt(3), e[3] = 2,
            e[4] = 2 * sqrt(3), e[5] = sqrt(5), e[6] = sqrt(3)-3*sqrt(2).

    Define 3 dep_pairs:

        sage: dep_pairs = []
        sage: dep_pairs = [V.subspace([e[0], e[3]]), \
        ....:             V.subspace([e[6], 2*e[2]-3*e[1]]), \
        ....:             V.subspace([e[2], e[4]]) ]
        sage: dependency = []
        sage: for dep_pair in dep_pairs: \
        ....:     dependency = update_dependency(dependency, dep_pair)
        sage: dependency
        [Vector space of degree 7 and dimension 2 over Rational Field
         Basis matrix:
         [1 0 0 0 0 0 0]
         [0 0 0 1 0 0 0], Vector space of degree 7 and dimension 2 over Rational Field
         Basis matrix:
         [   0    1 -2/3    0    0    0    0]
         [   0    0    0    0    0    0    1], Vector space of degree 7 and dimension 2 over Rational Field
         Basis matrix:
         [0 0 1 0 0 0 0]
         [0 0 0 0 1 0 0]]

    Define 19 ind_pairs:

        sage: ind_pairs = [V.subspace([e[0], e[1]]), \
        ....:              V.subspace([e[0], e[2]]), \
        ....:              V.subspace([e[0], e[4]]), \
        ....:              V.subspace([e[0], e[5]]), \
        ....:              V.subspace([e[0], e[6]]), \
        ....:              V.subspace([e[1], e[2]]), \
        ....:              V.subspace([e[1], e[3]]), \
        ....:              V.subspace([e[1], e[4]]), \
        ....:              V.subspace([e[1], e[5]]), \
        ....:              V.subspace([e[1], e[6]]), \
        ....:              V.subspace([e[2], e[3]]), \
        ....:              V.subspace([e[2], e[5]]), \
        ....:              V.subspace([e[2], e[6]]), \
        ....:              V.subspace([e[3], e[4]]), \
        ....:              V.subspace([e[3], e[5]]), \
        ....:              V.subspace([e[3], e[6]]), \
        ....:              V.subspace([e[4], e[5]]), \
        ....:              V.subspace([e[4], e[6]]), \
        ....:              V.subspace([e[5], e[6]])   ]

    Define `seen_linear_forms` as the set of the generators of `ind_pairs`:

        sage: seen_linear_forms = set(e)


        sage: independency, zero_kernel = \
        ....:     construct_independency(ind_pairs, dependency, seen_linear_forms)

        sage: len(independency)
        8

    Only need 8 pairs to represent the above independency (which had 19 pairs):

        sage: get_independent_pairs_from_independency(independency)
        [((1, 0, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0, 0)),
         ((0, 1, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1, 0)),
         ((1, 0, 0, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0, 0)),
         ((1, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1, 0)),
         ((1, 0, 0, 0, 0, 0, 0), (0, 1, -2/3, 0, 0, 0, 0)),
         ((0, 0, 1, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1, 0)),
         ((0, 1, 0, 0, 0, 0, 0), (0, 0, 1, 0, 0, 0, 0)),
         ((0, 1, -2/3, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1, 0))]

     With the additional conditions that the following forms have kernel = {0}.
        sage: zero_kernel
        {(0, 0, 0, 0, 0, 0, 1),
         (0, 0, 0, 0, 1, 0, 0),
         (0, 0, 0, 1, 0, 0, 0),
         (0, 0, 1, 0, 0, 0, 0),
         (0, 1, -2/3, 0, 0, 0, 0),
         (1, 0, 0, 0, 0, 0, 0)}
    """
    independency = []
    zero_kernel = set([]) # zero_kernel = seen_linear_forms?
    for t in seen_linear_forms:
        to_add = True
        for d in dependency:
            if t in d:
                to_add = False
                break
        if to_add:
            vector_space = VectorSpace(QQ,len(t))
            dependency += [vector_space.subspace([vector_space(t)])]
    for ind_pair_1 in independent_pairs:
        for ind_pair_2 in independent_pairs:
            if ind_pair_1 < ind_pair_2:
                i = ind_pair_1.intersection(ind_pair_2)
                if i.dimension() > 0:
                    to_add = True
                    for d in dependency:
                        if i.gen(0) in d:
                            to_add = False
                            break
                    if to_add:
                        dependency += [i]
    for ind_pair in independent_pairs:
        simultaneously_independency, intersections = construct_simultaneously_independency(ind_pair, dependency)
        independency = update_independency(independency, simultaneously_independency)
        zero_kernel.update(intersections)
    return independency, zero_kernel

def construct_simultaneously_independency(ind_pair, dependency):
    """
    EXAMPLE_1:

        sage: V = VectorSpace(QQ,5)
        sage: E1 = V.subspace([V([1,0,0,0,0]), V([0,0,1,0,0])])
        sage: E2 = V.subspace([V([0,1,0,0,0]), V([0,0,0,1,0])])
        sage: E3 = V.subspace([V([1,1,0,0,0]), V([0,0,0,0,1])])
        sage: dependency = [E1, E2, E3]
        sage: ind_pair = V.subspace([V([1,0,0,0,0]), V([0,1,0,0,0])])
        sage: sim_ind, intersections = \
        ....:     construct_simultaneously_independency(ind_pair, dependency)
        sage: sim_ind == set([(E1, E2), (E1, E3), (E2, E3)])
        True
        sage: intersections
        {(0, 1, 0, 0, 0), (1, 0, 0, 0, 0), (1, 1, 0, 0, 0)}

    EXAMPLE_2:

    Imagine e[0] =1, e[1] = sqrt(2), e[2] = sqrt(3), e[3] = 2,
            e[4] = 2 * sqrt(3), e[5] = sqrt(5), e[6] = sqrt(3)-3*sqrt(2).

    Use the same dep_pairs as in the last example of update_dependency()
    to get the input value for dependency.
    Define 4 pairs of independent relations as follows.

        sage: V = VectorSpace(QQ,7)
        sage: dependency = [ \
        ....:     V.subspace([V([1,0,0,0,   0,0,0]), \
        ....:                 V([0,0,0,1,   0,0,0])]), \
        ....:     V.subspace([V([0,1,0,0,-1/3,0,0]),\
        ....:                 V([0,0,1,0,-1/2,0,0]),\
        ....:                 V([0,0,0,0,   0,0,1])]) ]

    Define ind_pair whose generators are already in dependency:

        sage: ind_pair = V.subspace([V([1,0,0,0,0,0,0]), V([0,0,0,0,0,0,1])])
        sage: sim_ind, intersections = \
        ....:     construct_simultaneously_independency(ind_pair, dependency)
        sage: sim_ind
        {(Vector space of degree 7 and dimension 2 over Rational Field
          Basis matrix:
          [1 0 0 0 0 0 0]
          [0 0 0 1 0 0 0], Vector space of degree 7 and dimension 3 over Rational Field
          Basis matrix:
          [   0    1    0    0 -1/3    0    0]
          [   0    0    1    0 -1/2    0    0]
          [   0    0    0    0    0    0    1])}
        sage: intersections
        {(0, 0, 0, 0, 0, 0, 1), (1, 0, 0, 0, 0, 0, 0)}

    Define other ind_pairs whose generators are not all in dependency: 

        sage: ind_pair =  V.subspace([V([0,0,0,0,0,1,0]), V([0,0,0,0,0,0,1])])
        sage: construct_simultaneously_independency(ind_pair, dependency)
        ({(Vector space of degree 7 and dimension 3 over Rational Field
           Basis matrix:
           [   0    1    0    0 -1/3    0    0]
           [   0    0    1    0 -1/2    0    0]
           [   0    0    0    0    0    0    1], (0, 0, 0, 0, 0, 1, 0))},
         {(0, 0, 0, 0, 0, 0, 1)})

        sage: ind_pair = V.subspace([V([0,0,0,0,5,0,0]), V([1,0,0,0,0,0,0])])
        sage: construct_simultaneously_independency(ind_pair, dependency)
        ({(Vector space of degree 7 and dimension 2 over Rational Field
           Basis matrix:
           [1 0 0 0 0 0 0]
           [0 0 0 1 0 0 0], (0, 0, 0, 0, 1, 0, 0))}, {(1, 0, 0, 0, 0, 0, 0)})

        sage: ind_pair = V.subspace([V([0,0,0,0,0,1,0]), V([0,9,8,0,0,0,0])])
        sage: construct_simultaneously_independency(ind_pair, dependency)
        ({((0, 1, 8/9, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1, 0))}, set())

    Try again with generators of ind_pairs being added to 'dependency'.
    (This is what construct_independency() does first.) 

        sage: ind_pair = V.subspace([V([0,0,0,0,5,0,0]), V([1,0,0,0,0,0,0])])
        sage: construct_simultaneously_independency(ind_pair, \
        ....:     dependency + [V.subspace([V([0,0,0,0,5,0,0])])])
        ({(Vector space of degree 7 and dimension 2 over Rational Field
           Basis matrix:
           [1 0 0 0 0 0 0]
           [0 0 0 1 0 0 0],
           Vector space of degree 7 and dimension 1 over Rational Field
           Basis matrix:
           [0 0 0 0 1 0 0])},
         {(1, 0, 0, 0, 0, 0, 0)})
    """
    intersecting = []
    intersections = []
    for s in dependency:
        i = s.intersection(ind_pair)
        if i.dimension() > 0:
            intersecting.append(s)
            if s.dimension() > 1: # valid condition?
                intersections.append(i.gen(0))
    # NOTE: if construct_independency() was called, generators of ind_pair must be in `dependency`.
    # Then len(intersecting) > 1, the first two cases won't happen.
    if len(intersecting) == 0:
        return set([(ind_pair.gen(0), ind_pair.gen(1))]), set(intersections)
    elif len(intersecting) == 1:
        dependent_space = intersecting[0]
        # find orth_vec in (dependent_space + ind_pair) that is orthogonal to dependent_space
        orth_vec = find_orthogonal_vector(dependent_space, dependent_space + ind_pair)
        return set([(dependent_space, orth_vec)]), set(intersections)
    else:
        sim_ind_list = [(intersecting[i], intersecting[j]) for i in range(len(intersecting)) for j in range(i+1, len(intersecting))]
        return set(sim_ind_list), set(intersections)

def find_orthogonal_vector(s1, s2):
    s3 = s1.complement()
    orth_vec = s2.intersection(s3).gen(0)
    return orth_vec

def get_independent_pairs_from_independency(independency):
    independent_pairs = []
    for sim_ind in independency:
        s1, s2 = list(sim_ind)[0]
        t1 = s1.gen(0)
        t2 = s2.gen(0)
        vector_space = VectorSpace(QQ,len(t1))
        pair_space = vector_space.subspace([t1, t2])
        independent_pairs.append((pair_space.gen(0), pair_space.gen(1)))
    return independent_pairs
