# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

import sage.structure.element
from sage.structure.element import FieldElement

from sage.libs.ppl import Variable, Constraint, Linear_Expression, Constraint_System, NNC_Polyhedron, Poly_Con_Relation, Poly_Gen_Relation, Generator, MIP_Problem
poly_is_included = Poly_Con_Relation.is_included()
#strictly_intersects = Poly_Con_Relation.strictly_intersects()
point_is_included = Poly_Gen_Relation.subsumes()
#con_saturates = Poly_Con_Relation.saturates()

from sage.structure.sage_object import SageObject

###############################
# Parametric Real Number Field
###############################

default_parametric_field = None

class ParametricRealFieldElement(FieldElement):

    def __init__(self, value, symbolic=None, parent=None):
        if parent is None:
            raise ValueError, "ParametricRealFieldElement invoked with parent=None. That's asking for trouble"
            parent = default_parametric_field
        FieldElement.__init__(self, parent) ## this is so that canonical_coercion works.
        ## Test coercing the value to RR, so that we do not try to build a ParametricRealFieldElement
        ## from something like a tuple or vector or list or variable of a polynomial ring
        ## or something else that does not make any sense.
        if not isinstance(value, sage.interfaces.mathematica.MathematicaElement):
            RR(value)
        self._val = value
        if symbolic is None:
            self._sym = value # changed to not coerce into SR. -mkoeppe
        else:
            self._sym = symbolic

    def sym(self):
        return self._sym

    def val(self):
        return self._val

    def _cmp_(left, right):
        if not left.parent() is right.parent():
            raise TypeError, "comparing elements from different fields"
        result = cmp(left._val, right._val)
        if result == 0:
            left.parent().record_to_eq_list(left.sym() - right.sym())
        elif result == -1:
            left.parent().record_to_lt_list(left.sym() - right.sym())
        elif result == 1:
            left.parent().record_to_lt_list(right.sym() - left.sym())
        return result

    __cmp__ = _cmp_

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
        return "%s" % (latex(self._sym))

    def _add_(self, other):
        if not isinstance(other, ParametricRealFieldElement):
            other = ParametricRealFieldElement(other, parent=self.parent())
        return ParametricRealFieldElement(self._val + other._val, self._sym + other._sym, parent=self.parent())

    def _sub_(self, other):
        if not isinstance(other, ParametricRealFieldElement):
            other = ParametricRealFieldElement(other, parent=self.parent())
        return ParametricRealFieldElement(self._val - other._val, self._sym - other._sym, parent=self.parent())

    def __neg__(self):
        return ParametricRealFieldElement(-self._val, -self._sym, parent=self.parent())

    def _mul_(self, other):
        if not isinstance(other, ParametricRealFieldElement):
            try:
                other = ParametricRealFieldElement(other, parent=self.parent())
            except TypeError:
                # For example when other is a vector
                return other * self
        return ParametricRealFieldElement(self._val * other._val, self._sym * other._sym, parent=self.parent())

    def _div_(self, other):
        if not isinstance(other, ParametricRealFieldElement):
            other = ParametricRealFieldElement(other, parent=self.parent())
        return ParametricRealFieldElement(self._val / other._val, self._sym / other._sym, parent=self.parent())

    def __hash__(self):
        """
        The hash function of these elements must be so that 
        elements that would compare equal have the same hash value. 

        The constant hash function would do the job.  Instead we use the
        hash of the ._val (because equality implies equality of _val).
        It is not correct to use the hash of the ._sym, or to compare
        the __repr__, because then the user could check for equality
        (for example, by testing the cardinality of a set, as in the
        tests below) without the equality being recorded in the field.

        The correctness of this implementation depends on the guarantee
        of the that equal _val elements have the same hash value.  If in
        doubt, make sure that the _val elements all come from the same
        field, by `nice_field_values`.

        EXAMPLES::

            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: K.<f> = ParametricRealField([1])
            sage: s = {f, K(1)}
            sage: len(s)
            1
            sage: s
            {1}
            sage: K.<f> = ParametricRealField([1])
            sage: s = {f, K(2)}
            sage: len(s)
            2
        """
        return hash(self._val)

from sage.rings.ring import Field
import sage.rings.number_field.number_field_base as number_field_base
from sage.structure.coerce_maps import CallableConvertMap
from itertools import izip

class ParametricRealField(Field):
    """
    Parametric search:
    EXAMPLES::

        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: K.<f> = ParametricRealField([4/5])
        sage: h = gmic(f, field=K)
        sage: _ = generate_maximal_additive_faces(h);
        sage: K.get_eq_list()
        set()
        sage: K.get_lt_list()
        {-1/(-f^2 + f), -1/f, -2*f, -2*f + 1, -f - 1, -f, f - 2, f - 1, 2*f - 2}
        sage: K.get_eq_poly()
        set()
        sage: K.get_lt_poly()
        {-2*f, -2*f + 1, -f - 1, -f, f - 2, f - 1, 2*f - 2, f^2 - f}
        sage: K.get_eq_factor()
        set()
        sage: K.get_lt_factor()
        {-f, -f + 1/2, f - 1}

        sage: K.<f, lam> = ParametricRealField([4/5, 1/6])
        sage: h = gj_2_slope(f, lam, field=K)
        sage: K.get_lt_list()
        {(-1/2)/(-1/2*f^2*lam - 1/2*f^2 + f*lam + 1/2*f - 1/2*lam),
         (-f*lam - f + lam)/(-f + 1),
         -lam,
         lam - 1,
         -f,
         f - 1,
         -1/2*f*lam - 1/2*f + 1/2*lam,
         f*lam - lam}
        sage: K.get_lt_poly()
        {-lam,
         lam - 1,
         -f,
         f - 1,
         -f*lam - f + lam,
         -1/2*f*lam - 1/2*f + 1/2*lam,
         f*lam - lam,
         1/2*f^2*lam + 1/2*f^2 - f*lam - 1/2*f + 1/2*lam}
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
    """
    Element = ParametricRealFieldElement

    def __init__(self, values=[], names=()):
        Field.__init__(self, self)
        self._zero_element = ParametricRealFieldElement(0, parent=self)
        self._one_element =  ParametricRealFieldElement(1, parent=self)
        self._eq = set([])
        self._lt = set([])
        self._eq_poly = set([])
        self._lt_poly = set([])
        self._eq_factor = set([])
        self._lt_factor = set([])
        vnames = PolynomialRing(QQ, names).fraction_field().gens();
        self._gens = [ ParametricRealFieldElement(value, name, parent=self) for (value, name) in izip(values, vnames) ]
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
        return ParametricRealFieldElement(1, parent=self)
    def _coerce_map_from_(self, S):
        if isinstance(S, ParametricRealField) and self is not S:
            return None
        return CallableConvertMap(S, self, lambda s: ParametricRealFieldElement(s, parent=self), parent_as_first_arg=False)
    def __repr__(self):
        return 'ParametricRealField%s' %repr(self.gens())
    def _element_constructor_(self, elt):
        if parent(elt) == self:
            return elt
        try: 
            QQ_elt = QQ(elt)
            return ParametricRealFieldElement(QQ_elt, parent=self)
            #return self.element_class(self, QQ_elt)
        except:
            raise TypeError, "ParametricRealField called with element %s" % elt

    def _coerce_impl(self, x):
        return self._element_constructor_(x)
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
        if not comparison in QQ and not comparison in self._eq:
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


default_parametric_field = ParametricRealField()

###############################
# Simplify polynomials
###############################

def polynomial_to_linexpr(t, monomial_list, v_dict):
    """
    Reformulation–linearization: Expand the polynomial in the standard monomial basis and replace each monomial by a new variable. Record monomials in monomial_list and their corresponding variables in v_dict. The resulting linear expression in the extended space will be provided as inequality or equation in the linear system that describes a PPL not-necessarily-closed polyhedron.

    EXAMPLES::

        sage: P.<x,y,z> = QQ[]
        sage: monomial_list = []; v_dict = {};
        sage: t = 27/113 * x^2 + y*z + 1/2
        sage: polynomial_to_linexpr(t, monomial_list, v_dict)
        54*x0+226*x1+113
        sage: monomial_list
        [x^2, y*z]
        sage: v_dict
        {y*z: 1, x^2: 0}

        sage: tt = x + 1/3 * y*z
        sage: polynomial_to_linexpr(tt, monomial_list, v_dict)
        x1+3*x2
        sage: monomial_list
        [x^2, y*z, x]
        sage: v_dict
        {x: 2, y*z: 1, x^2: 0}
    """
    # coefficients in ppl constraint must be integers.
    lcd = lcm([x.denominator() for x in t.coefficients()])
    linexpr = Linear_Expression(0)
    if len(t.args()) <= 1:
        # sage.rings.polynomial.polynomial_rational_flint object has no attribute 'monomials'
        for (k, c) in t.dict().items():
            m = (t.args()[0])^k
            if m in v_dict.keys():
                v = Variable(v_dict[m])
            elif k == 0:
                # constant term, don't construct a new Variable for it.
                v = 1
            else:
                nv = len(monomial_list)
                v = Variable(nv)
                v_dict[m] = nv
                monomial_list.append(m)
            linexpr += (lcd * c) * v
    else:
        for m in t.monomials():
            if m in v_dict.keys():
                v = Variable(v_dict[m])
            elif m == 1:
                v = 1
            else:
                nv = len(monomial_list)
                v = Variable(nv)
                v_dict[m] = nv
                monomial_list.append(m)
            coeffv = t.monomial_coefficient(m)
            linexpr += (lcd * coeffv) * v
    return linexpr

def cs_of_eq_lt_poly(eq_poly, lt_poly):
    """
    Reformulation–linearization: Expand the polynomials in the standard monomial basis and replace each monomial by a new variable. Construct a linear constraint system in the extended space, which describes a PPL not-necessarily-closed polyhedron. Record monomials in monomial_list and their corresponding variables in v_dict.

    EXAMPLES::

        sage: P.<f>=QQ[]
        sage: eq_poly =[]; lt_poly = [2*f - 2, f - 2, f^2 - f, -2*f, f - 1, -f - 1, -f, -2*f + 1]
        sage: cs, monomial_list, v_dict = cs_of_eq_lt_poly(eq_poly, lt_poly)
        sage: cs
        Constraint_System {-x0+1>0, -x0+2>0, x0-x1>0, x0>0, -x0+1>0, x0+1>0, x0>0, 2*x0-1>0}
        sage: monomial_list
        [f, f^2]
        sage: v_dict
        {f: 0, f^2: 1}
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
        sage: K.<f> = ParametricRealField([4/5])
        sage: h = gmic(f, field=K)
        sage: _ = extremality_test(h)
        sage: eq_poly = K.get_eq_poly()
        sage: lt_poly = K.get_lt_poly()
        sage: (eq_poly, lt_poly)
        (set(), {-2*f, -2*f + 1, -f - 1, -f, f - 2, f - 1, 2*f - 2, f^2 - f})
        sage: simplify_eq_lt_poly_via_ppl(eq_poly, lt_poly)
        ([], [f - 1, -2*f + 1, f^2 - f])

        sage: eq_factor = K.get_eq_factor()
        sage: lt_factor = K.get_lt_factor()
        sage: (eq_factor, lt_factor)
        (set(), {-f, -f + 1/2, f - 1})
        sage: simplify_eq_lt_poly_via_ppl(eq_factor, lt_factor)
        ([], [f - 1, -2*f + 1])

        sage: K.<f, lam> = ParametricRealField([4/5, 1/6])
        sage: h = gj_2_slope(f, lam, field=K, conditioncheck=False)
        sage: leq, lin = simplify_eq_lt_poly_via_ppl(K.get_eq_poly(), K.get_lt_poly())
        sage: set(lin)
        {-lam, f - 1, -f*lam - f + lam, f*lam - lam, f^2*lam + f^2 - 2*f*lam - f + lam}
        sage: leq, lin = simplify_eq_lt_poly_via_ppl(list(K.get_eq_factor()), list(K.get_lt_factor()))
        sage: set(lin)
        {-lam, -f, f - 1, -f*lam - f + lam}

        sage: _ = extremality_test(h)
        sage: leq, lin = simplify_eq_lt_poly_via_ppl(K.get_eq_poly(), K.get_lt_poly())
        sage: set(lin)
        {-lam,
         lam - 1,
         -2*f*lam + f + 2*lam - 1,
         -f*lam - 3*f + lam + 2,
         f*lam - lam,
         f^2*lam + f^2 - 2*f*lam - f + lam}
        sage: leq, lin = simplify_eq_lt_poly_via_ppl(K.get_eq_factor(), K.get_lt_factor())
        sage: set(lin)
        {-lam,
         2*lam - 1,
         f - 1,
         -3*f*lam - f + 3*lam,
         -f*lam - 3*f + lam + 2,
         f*lam - 3*f - lam + 2,
         3*f*lam - f - 3*lam}

        sage: K.<f,alpha> = ParametricRealField([4/5, 3/10])             # Bad example! parameter region = {given point}.
        sage: h=dg_2_step_mir(f, alpha, field=K, conditioncheck=False)
        sage: _ = extremality_test(h)
        sage: leq, lin = simplify_eq_lt_poly_via_ppl(K.get_eq_poly(), K.get_lt_poly())
        sage: set(leq), set(lin)
        ({-10*alpha + 3, -5*f + 4}, {5*f^2 - 10*f*alpha - 1})

        sage: leq, lin = simplify_eq_lt_poly_via_ppl(K.get_eq_factor(), K.get_lt_factor())
        sage: set(leq), set(lin)
        ({-10*alpha + 3, -5*f + 4}, set())

        sage: K.<f> = ParametricRealField([1/5])
        sage: h = drlm_3_slope_limit(f, conditioncheck=False)
        sage: _ = extremality_test(h)
        sage: leq, lin = simplify_eq_lt_poly_via_ppl(K.get_eq_poly(), K.get_lt_poly())
        sage: set(leq), set(lin)
        (set(), {-f, 3*f - 1, -f^2 - f})
        sage: leq, lin = simplify_eq_lt_poly_via_ppl(list(K.get_eq_factor()), list(K.get_lt_factor()))
        sage: set(leq), set(lin)
        (set(), {-f, 3*f - 1})
    """
    cs, monomial_list, v_dict = cs_of_eq_lt_poly(eq_poly, lt_poly)
    p = NNC_Polyhedron(cs)
    return read_leq_lin_from_polyhedron(p, monomial_list, v_dict)


def read_leq_lin_from_polyhedron(p, monomial_list, v_dict, tightened_mip=None):
    """
    Given a PPL polyhedron p, map the minimal constraints system of p back to polynomial equations and inequalities in the orginal space. If a constraint in the minimal constraint system of p is not tight for the tightened_mip, then this constraint is discarded.

    EXAMPLES::

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
        if tightened_mip is not None and is_not_a_downstairs_wall(c, tightened_mip):
            # McCormick trash, don't put in minlt.
            continue
        coeff = c.coefficients()
        # observe: coeffients in a constraint of NNC_Polyhedron could have gcd != 1.
        # take care of this.
        gcd_c = gcd(gcd(coeff), c.inhomogeneous_term())
        # constraint is written with '>', while lt_poly records '<' relation
        t = sum([-(x/gcd_c)*y for x, y in itertools.izip(coeff, monomial_list)]) - c.inhomogeneous_term()/gcd_c
        if c.is_equality():
            mineq.append(t)
        else:
            minlt.append(t)
    # note that polynomials in mineq and minlt can have leading coefficient != 1
    return mineq, minlt

def read_simplified_leq_lin(K, level="factor"):
    """
    Use the reformulation–linearization techinque to remove redundant inequalties and equations recorded in ParametricRealField K.

    EXAMPLES::

        sage: K.<f> = ParametricRealField([4/5])
        sage: h = gmic(f, field=K)
        sage: _ = extremality_test(h)
        sage: read_simplified_leq_lin(K)
        ([], [f - 1, -2*f + 1])

        sage: K.<f> = ParametricRealField([1/5])
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

def find_variable_mapping(leqs, lins):
    """
    Assume that leqs, lins are the results from read_leq_lin_from_polyhedron(),
so that gaussian elimination has been performed by PPL on the list of equations. If an equation has a linear variable that does not appear in other equations, then eliminate this variable.
    FIXME: This function is bad; it assumes too many things.

    EXAMPLES::

        sage: logging.disable(logging.WARN)
        sage: P.<a,b,c>=QQ[]
        sage: leqs = [a+2*b+1]; lins=[a+b+c, b+2*c-3]
        sage: find_variable_mapping(leqs, lins)
        {c: c, b: b, a: -2*b - 1}
        sage: leqs = [2*a+b, b+c-1/2]; lins=[a+b+c]
        sage: find_variable_mapping(leqs, lins)
        {c: -b + 1/2, b: b, a: -1/2*b}
        sage: leqs = [a*a-b, b*b-c*c]; lins=[a+b+c]
        sage: find_variable_mapping(leqs, lins)
        {c: c, b: b, a: a}
        sage: P.<d>=QQ[]
        sage: leqs = [1-d^3]; lins = []
        sage: find_variable_mapping(leqs, lins)
        {d: d}
    """
    if leqs:
        variables = leqs[0].args()
    elif lins:
        variables = lins[0].args()
    else:
        raise ValueError, "constraints are not provided."
    var_map = {}
    for v in variables:
        var_map[v] = v
    if not leqs:
        return var_map
    n = len(leqs)
    if len(variables) == 1: # workaround for single variable 'Polynomial_rational_flint'
        v = variables[0]
        for i in range(n):
            if leqs[i].degree() == 1:
                var_map[v] = -leqs[i].list()[0]/leqs[i].list()[1]
                return var_map
        logging.warn("Can't solve for %s in the system %s == 0, %s < 0. Heurist wall crossing may fail." % (v, leqs, lins))
        return var_map
    for i in range(n):
        found_pivot = False
        for v in variables:
            if leqs[i].coefficient(v).degree() == 0: # v is a linear variable in leqs[i].
                coef = leqs[i].monomial_coefficient(v) # type is rational
                if all((v not in leqs[j].variables()) for j in range(n) if j != i):
                    found_pivot = True
                    var_map[v] = v - leqs[i] / coef # eliminate v
                    break
        if not found_pivot:
            logging.warn("Can't find linear variable in %s == 0 to eliminate in the system %s == 0, %s < 0. Heurist wall crossing may fail." % (leqs[i], leqs, lins))
    return var_map

######################################
# Functions with ParametricRealField K
######################################

from sage.misc.sageinspect import sage_getargspec, sage_getvariablename

def read_default_args(function, **opt_non_default):
    """
    Return the default values of arguments of the function.
    Override the default values if opt_non_default is given.

    EXAMPLES::

        sage: read_default_args(gmic)
        {'conditioncheck': True, 'f': 4/5, 'field': None}
        sage: read_default_args(drlm_backward_3_slope, **{'bkpt': 1/5})
        {'bkpt': 1/5, 'conditioncheck': True, 'f': 1/12, 'field': None}
    """
    args, varargs, keywords, defaults = sage_getargspec(function)
    default_args = {}
    if defaults is not None:
        for i in range(len(defaults)):
            default_args[args[-i-1]]=defaults[-i-1]
    for (opt_name, opt_value) in opt_non_default.items():
        if opt_name in default_args:
            default_args[opt_name] = opt_value
    return default_args

def construct_field_and_test_point(function, var_name, var_value, default_args):
    """
    Construct a ParametricRealField K using var_name and var_value.
    var_name and var_value are two parallel lists.
    Construct a test_point of type dictionary, which maps each parameter of the function to the corresponding ParametricRealFieldElement if this is a parameter of K, otherwise maps to the default argument value.

    EXAMPLES::

        sage: function=gmic; var_name=['f']; var_value=[1/2];
        sage: default_args = read_default_args(function)
        sage: K, test_point = construct_field_and_test_point(function, var_name, var_value, default_args)
        sage: K
        ParametricRealField[f~]
        sage: test_point
        {'conditioncheck': False, 'f': f~, 'field': ParametricRealField[f~]}
    """
    K = ParametricRealField(var_value, var_name)
    test_point = copy(default_args)
    for i in range(len(var_name)):
        test_point[var_name[i]] = K.gens()[i]
    args_set = set(sage_getargspec(function)[0])
    if 'field' in args_set:
        test_point['field'] = K
    if 'conditioncheck' in args_set:
        test_point['conditioncheck'] = False
    return K, test_point

def simplified_extremality_test(function):
    """
    A simplified version of the extremality test, which assumes that
    the given function is minimal valid and has rational bkpts.
    Return True or False, without computing the perturbations.
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

###########################
# Super class that stores the computation so far
###########################
class SemialgebraicComplexComponent(SageObject):

    def __init__(self, parent, K, var_value, region_type):
        self.parent = parent
        self.var_value = var_value
        self.region_type = region_type

        #note that self.polyhedron and K.polyhedron,
        # self.parent.monomial_list and K.monomial_list,
        # self.parent.v_dict and K.v_dict change simultaneously while lifting.
        self.polyhedron = K.polyhedron
        self.bounds, tightened_mip = self.bounds_propagation(self.parent.max_iter)
        if self.parent.max_iter == 0:
            tightened_mip = None

        leqs, lins = read_leq_lin_from_polyhedron(self.polyhedron, \
                                                          self.parent.monomial_list, self.parent.v_dict, tightened_mip)
        self.var_map = find_variable_mapping(leqs, lins)
        self.leq = leqs
        self.lin = []
        for l in lins:
            ineq = l.subs(self.var_map)
            if ineq.degree() > 0:
                self.lin.append(ineq)

    def bounds_propagation(self, max_iter):
        tightened_mip = construct_mip_of_nnc_polyhedron(self.polyhedron)
        # Compute LP bounds first
        bounds = [find_bounds_of_variable(tightened_mip, i) for i in range(len(self.parent.monomial_list))]

        bounds_propagation_iter = 0
        # TODO: single parameter polynomial_rational_flint object needs special treatment. ignore bounds propagation in this case for now.  #tightened = True
        tightened = bool(len(self.var_value) > 1)

        while bounds_propagation_iter < max_iter and tightened:
            bounds_propagation_iter += 1
            tightened = False
            # upward bounds propagation
            for i in range(len(self.var_value), len(self.parent.monomial_list)):
                m = self.parent.monomial_list[i] # m has degre >= 2
                if update_mccormicks_for_monomial(m, tightened_mip, self.polyhedron, \
                                                  self.parent.monomial_list, self.parent.v_dict, bounds):
                    tightened = True
            if tightened:
                tightened = False
                # downward bounds propagation
                for i in range(len(self.parent.monomial_list)):
                    (lb, ub) = bounds[i]
                    bounds[i] = find_bounds_of_variable(tightened_mip, i)
                    # Not sure about i < len(self.var_value) condition,
                    # but without it, bounds_propagation_iter >= max_iter is often attained.
                    #if (i < len(self.var_value)) and \
                    #   ((lb < bounds[i][0]) or (lb is None) and (bounds[i][0] is not None) or \
                    #    (bounds[i][1] < ub) or (ub is None) and (bounds[i][1] is not None)):
                    if (bounds[i][0] is not None) and ((lb is None) or (bounds[i][0] - lb > 0.001)) or \
                       (bounds[i][1] is not None) and ((ub is None) or (ub - bounds[i][1] > 0.001)):
                        tightened = True
            if max_iter != 0 and bounds_propagation_iter >= max_iter:
                logging.warn("max number %s of bounds propagation iterations has attained." % max_iter)
        #print bounds_propagation_iter
        return bounds, tightened_mip

    def plot(self, alpha=0.5, plot_points=300, slice_value=None, show_testpoints=True):
        g = Graphics()
        if not slice_value:
            d = len(self.var_value)
            if d > 2:
                raise NotImplementedError, "Plotting region with dimension > 2 is not implemented. Provide `slice_value` to plot a slice of the region."
            leqs = self.leq
            lins = self.lin + self.parent.bddlin
            var_bounds = [bounds_for_plotting(self.bounds[i], self.parent.default_var_bound) for i in range(d)]
        else:
            # assert (len(slice_value) == len(self.var_value))
            d = 0
            var_bounds = []
            for (i, z) in enumerate(slice_value):
                if z is None:
                    d += 1
                    var_bounds.append(bounds_for_plotting(self.bounds[i], self.parent.default_var_bound))
                elif not is_value_in_interval(z, self.bounds[i]):
                    # empty slice
                    return g
            if d == 1:
                P.<unknown_x> = QQ[]
            elif d == 2:
                P.<unknown_x, unknown_y> = QQ[]
            else:
                raise NotImplementedError, "Plotting region with dimension > 2 is not implemented. Provide `slice_value` to plot a slice of the region."
            unknowns=iter(P.gens())
            parametric_point = [unknowns.next() if z is None else z for z in slice_value]
            leqs = []
            for leq in self.leq:
                l = leq(*parametric_point)
                if l in QQ:
                    if l != 0:
                        return g
                else:
                   leqs.append(l)
            lins = []
            for lin in self.lin + self.parent.bddlin:
                l = lin(*parametric_point)
                if l in QQ:
                    if l >= 0:
                        return g
                else:
                    lins.append(l)

        innercolor = find_region_color(self.region_type)
        bordercolor = innercolor
        if innercolor == 'white':
            ptcolor = 'black'
        else:
            ptcolor = 'white'

        x, y = var('x, y')

        if d == 2:
            if leqs or lins:
                g += region_plot([l(x, y) == 0 for l in leqs] + [l(x, y) < 0 for l in lins], \
                                 (x, var_bounds[0][0], var_bounds[0][1]), (y, var_bounds[1][0], var_bounds[1][1]), \
                                 incol=innercolor, alpha=alpha, plot_points=plot_points, bordercol=bordercolor)
            if show_testpoints and not slice_value:
                if (var_bounds[0][0] <= self.var_value[0] <= var_bounds[0][1]) and \
                   (var_bounds[1][0] <= self.var_value[1] <= var_bounds[1][1]):
                    g += point(self.var_value, color = ptcolor, size = 2, zorder=10)
        else:
            if leqs or lins:
                g += region_plot([l(x) == 0 for l in leqs] + [l(x) < 0 for l in lins] + [y >= -0.01, y <= 0.01], \
                                 (x, var_bounds[0][0], var_bounds[0][1]), (y, -0.1, 0.3), \
                                 incol=innercolor, alpha=alpha, plot_points=plot_points, bordercol=bordercolor, ticks=[None,[]])
            if show_testpoints and not slice_value:
                if (var_bounds[0][0] <= self.var_value[0] <= var_bounds[0][1]):
                    g += point([self.var_value[0], 0], color = ptcolor, size = 2, zorder=10)
        return g

    def find_walls_and_new_points(self, flip_ineq_step, wall_crossing_method, goto_lower_dim=True):
        #self.lin self.leq, self.parent.bddlin
        walls = []
        new_points = {}
        # decide which inequalities among self.lin are walls (irredundant).
        for i in range(len(self.lin)):
            ineq = self.lin[i]
            ineqs = walls + self.lin[i+1::] + self.parent.bddlin
            if wall_crossing_method == 'mathematica':
                condstr_others = write_mathematica_constraints(self.leq, ineqs)
                # maybe shouldn't put self.leq into FindInstance, but solve using var_map later.
                condstr_ineq = '0<'+str(ineq)+'<'+str(flip_ineq_step)
                pt_across_wall = find_instance_mathematica(condstr_others + condstr_ineq, self.parent.var_name)
            else:
                pt = find_point_flip_ineq_heuristic(self.var_value, ineq, ineqs, flip_ineq_step)
                if pt is None:
                    pt_across_wall = None
                else:
                    pt_across_wall = tuple(self.var_map[v](pt) for v in ineq.args())
            if pt_across_wall is None:
                # ineq is not an wall
                continue
            walls.append(ineq)
            new_points[pt_across_wall] = copy(self.leq) # assume all linear equations.
            # or maybe new_points[pt_across_wall] = copy(self.parent.bddleq)?
            if goto_lower_dim is True:
                pt_on_wall = None
                if wall_crossing_method == 'mathematica':
                    # wall could be non-linear, contrasting the assumption above. 
                    condstr_ineq = str(ineq) + '==0'
                    pt_on_wall = find_instance_mathematica(condstr_others + condstr_ineq, self.parent.var_name)
                elif ineq.degree() == 1:
                    # is linear wall. try to find a point on the wall using heuristic gradient descent method.
                    pt = find_point_on_ineq_heuristic(pt_across_wall, ineq, ineqs, flip_ineq_step)
                    if pt is None:
                        pt_on_wall = None
                    else:
                        pt_on_wall = tuple(self.var_map[v](pt) for v in ineq.args())
                if not pt_on_wall is None:
                    new_points[pt_on_wall] = copy(self.leq) + [ineq]
        return walls, new_points

class SemialgebraicComplex(SageObject):
    """
    EXAMPLES::

        sage: logging.disable(logging.WARN)
        sage: def vol(a,b):
        ....:     P = Polyhedron(ieqs=[(1,0,-1),(0,0,1),(1,-1,0),(0,1,0),(1,-a,-b)])
        ....:     return P.volume()
        sage: K.<a,b>=ParametricRealField([2,3])

    """

    def __init__(self, function, var_name, max_iter=8, find_region_type=None, default_var_bound=(-0.1,1.1), bddleq=[], bddlin=[], **opt_non_default):
        #self.num_components = 0
        self.components = []

        self.function = function
        self.d = len(var_name)
        self.var_name = var_name
        self.default_args = read_default_args(function, **opt_non_default)

        self.monomial_list = []
        self.v_dict = {}
        K = ParametricRealField([0]*self.d, var_name)
        for i in range(self.d):
            v = K.gens()[i].sym().numerator()
            self.monomial_list.append(v)
            self.v_dict[v] = i
        self.graph = Graphics()
        self.num_plotted_components = 0
        self.points_to_test = {} # a dictionary of the form {testpoint: bddleq}
        self.max_iter = max_iter
        if find_region_type is None:
            self.find_region_type = find_region_type_igp
        else:
            self.find_region_type = find_region_type
        self.default_var_bound = default_var_bound
        self.bddleq = bddleq
        self.bddlin = bddlin

    def generate_random_var_value(self, var_bounds=None):
        # var_bounds could be defined as lambda functions, see the testcase dg_2_step_mir in param_graphics.sage.
        # FIXME: if self.bddleq is not empty, it never ends.
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
                    x = QQ(uniform(l, u))
                var_value.append(x)
            # if random point doesn't satisfy self.bddleq or self.bddlin, continue while.
            if point_satisfies_bddleq_bddlin(var_value, self.bddleq, self.bddlin, strict=False):
                return var_value

    def is_point_covered(self, var_value):
        if all(x in QQ for x in var_value):
            #FIXME: is going through ppl the best way?
            monomial_value = [m(var_value) for m in self.monomial_list]
            # coefficients in ppl point must be integers.
            lcm_monomial_value = lcm([x.denominator() for x in monomial_value])
            #print [x * lcm_monomial_value for x in monomial_value]
            pt = Generator.point(Linear_Expression([x * lcm_monomial_value for x in monomial_value], 0), lcm_monomial_value)
            for c in self.components:
                # Check if the random_point is contained in the box.
                if c.region_type == 'not_constructible' and c.leq == [] and c.lin == []:
                    continue
                if is_point_in_box(monomial_value, c.bounds):
                    # Check if all eqns/ineqs are satisfied.
                    if c.polyhedron.relation_with(pt).implies(point_is_included):
                        return True
        else:
            for c in self.components:
                if point_satisfies_bddleq_bddlin(var_value, c.leq, c.lin, strict=True):
                    return True
        return False
        
    def find_uncovered_random_point(self, var_bounds=None, max_failings=1000):
        num_failings = 0
        while not max_failings or num_failings < max_failings:
            if self.points_to_test:
                var_value = list(self.points_to_test.popitem()[0])
            else:
                var_value = self.generate_random_var_value(var_bounds=var_bounds)
            # This point is not already covered.
            if self.is_point_covered(var_value):
                num_failings += 1
            else:
                return var_value
        logging.warn("The graph has %s components. Cannot find one more uncovered point by shooting %s random points" % (len(self.components), max_failings))
        return False

    def find_uncovered_point_mathematica(self, strict=True):
        condstr = write_mathematica_constraints(self.bddleq, self.bddlin, strict=True) #why strict = strict doesn't work when goto_lower_dim=False?
        for c in self.components:
            condstr_c = write_mathematica_constraints(c.leq, c.lin, strict=strict)
            if condstr_c:
                condstr += '!(' + condstr_c[:-4] + ') && '
        if not condstr:
            return tuple([0]*len(self.var_name))
        return find_instance_mathematica(condstr[:-4], self.var_name)

    def add_new_component(self, var_value, bddleq=[], flip_ineq_step=0, wall_crossing_method=None, goto_lower_dim=True):
        # Remark: the sign of flip_ineq_step indicates how to search for neighbour testpoints:
        # if flip_ineq_step = 0, don't search for neighbour testpoints. Used in shoot_random_points().
        # if flip_ineq_step < 0, we assume that the walls of the cell are linear eqn/ineq over original parameters.(So, gradient is constant; easy to find a new testpoint on the wall and another testpoint (-flip_ineq_step away) across the wall.) Used in bfs.
        # if flip_ineq_step > 0, we don't assume the walls are linear. Apply generate_one_point_by_flipping_inequality() with flip_ineq_step to find new testpoints across the wall only. Used in bfs.

        unlifted_space_dim =  len(self.monomial_list)
        K, test_point = construct_field_and_test_point(self.function, self.var_name, var_value, self.default_args)
        K.monomial_list = self.monomial_list # change simultaneously while lifting
        K.v_dict = self.v_dict # change simultaneously while lifting
        K.polyhedron.add_space_dimensions_and_embed(len(K.monomial_list))
        for l in bddleq:
            # need to put these equations in K, so call comparaison.
            if not l(*K.gens()) == 0:
                logging.warn("Test point %s doesn't satisfy %s == 0." % (var_value, l))
                return
        try:
            h = self.function(**test_point)
        except: # Dangerous!!
            # Function is non-contructible at this random point.
            h = None
        region_type = self.find_region_type(K, h)
        new_component = SemialgebraicComplexComponent(self, K, var_value, region_type)
        #if see new monomial, lift polyhedrons of the previously computed components.
        dim_to_add = len(self.monomial_list) - unlifted_space_dim
        if dim_to_add > 0:
            for c in self.components:
                c.polyhedron.add_space_dimensions_and_embed(dim_to_add)
        if flip_ineq_step != 0:
            # when using random shooting, don't generate neighbour points; don't remove redundant walls.
            walls, new_points = new_component.find_walls_and_new_points(flip_ineq_step, wall_crossing_method, goto_lower_dim)
            new_component.lin = walls
            self.points_to_test.update(new_points)
        self.components.append(new_component)

    def shoot_random_points(self, num, var_bounds=None, max_failings=1000):
        for i in range(num):
            var_value = self.find_uncovered_random_point(var_bounds=var_bounds, max_failings=max_failings)
            if var_value is False:
                return
            else:
                self.add_new_component(var_value, bddleq=[], flip_ineq_step=0, goto_lower_dim=False)

    def plot(self, alpha=0.5, plot_points=300, slice_value=None, restart=False):
        if restart:
            self.graph = Graphics()
            self.num_plotted_components = 0
        for c in self.components[self.num_plotted_components::]:
            self.graph += c.plot(alpha=alpha, plot_points=plot_points, slice_value=slice_value)
        self.num_plotted_components = len(self.components)
        return self.graph

    def bfs_completion(self, var_value=None, flip_ineq_step=1/100, check_completion=False, wall_crossing_method='heuristic', goto_lower_dim=False):
        if not self.components and not self.points_to_test and not var_value:
            var_value = self.find_uncovered_random_point()
        if var_value and not tuple(var_value) in self.points_to_test:
            self.points_to_test[tuple(var_value)] = copy(self.bddleq)
        while self.points_to_test:
            var_value, bddleq = self.points_to_test.popitem()
            var_value = list(var_value)
            if not self.is_point_covered(var_value):
                self.add_new_component(var_value, flip_ineq_step=flip_ineq_step, wall_crossing_method=wall_crossing_method, goto_lower_dim=goto_lower_dim)
        if check_completion:
            uncovered_pt = self.find_uncovered_point_mathematica(strict=goto_lower_dim)
            if uncovered_pt is not None:
                logging.warn("After bfs, the complex has uncovered point %s." % (uncovered_pt,))
                self.bfs_completion(var_value=uncovered_pt, \
                                    flip_ineq_step=flip_ineq_step, \
                                    check_completion=check_completion)

def gradient(ineq):
    # need this function since when K has only one variable,
    # got AttributeError: 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint' object has no attribute 'gradient'
    if hasattr(ineq, 'gradient'):
       return ineq.gradient()
    else:
       return [ineq.derivative()]

def point_satisfies_bddleq_bddlin(var_value, bddleq, bddlin, strict=True):
    # for functions involving ceil/floor, might be a good to devide region of search, then glue components together.
    for l in bddleq:
        if not l(var_value) == 0:
            return False
    for l in bddlin:
        if l(var_value) > 0 or (strict and l(var_value)==0):
            return False
    return True

def is_value_in_interval(v, (lb, ub)):
    return ((lb is None) or (lb <= v)) and ((ub is None) or (v <= ub))

def is_point_in_box(monomial_value, bounds):
    # note: monomial_value can have length bigger than len(bounds),
    # just ignore the tailing monomial values which come from lifting.
    return all(is_value_in_interval(monomial_value[i], bounds[i]) for i in range(len(bounds)))

def bounds_for_plotting((lb, ub), default_var_bound):
    if not lb is None:
        l = lb - 0.01
    else:
        l = default_var_bound[0]
    if not ub is None:
        u = ub + 0.01
    else:
        u = default_var_bound[1]
    return (l, u)

def construct_mip_of_nnc_polyhedron(nncp):
    min_cs = nncp.minimized_constraints()
    cs = Constraint_System()
    for c in min_cs:
        if c.is_equality():
            cs.insert(Linear_Expression(c.coefficients(), c.inhomogeneous_term()) == 0)
        else:
            cs.insert(Linear_Expression(c.coefficients(), c.inhomogeneous_term()) >= 0)
    mip = MIP_Problem(nncp.space_dimension())
    mip.add_constraints(cs)
    mip.set_optimization_mode('minimization')
    return mip

def find_bounds_of_variable(mip, i):
    linexpr = Linear_Expression(Variable(i))
    lb = find_lower_bound_of_linexpr(mip, linexpr)
    ub = find_upper_bound_of_linexpr(mip, linexpr)
    return (lb, ub)

def find_lower_bound_of_linexpr(mip, linexpr):
    # assume mip.set_optimization_mode('minimization')
    mip.set_objective_function(linexpr)
    try:
        lb = mip.optimal_value()
    except:
        # unbounded
        lb = None
    return lb

def find_upper_bound_of_linexpr(mip, linexpr):
    # assume mip.set_optimization_mode('minimization')
    mip.set_objective_function(-linexpr)
    try:
        ub = -mip.optimal_value()
    except:
        # unbounded
        ub = None
    return ub

def is_not_a_downstairs_wall(c, mip):
    linexpr = Linear_Expression(c.coefficients(), c.inhomogeneous_term())
    # we know lb exists and is >= 0
    lb = find_lower_bound_of_linexpr(mip, linexpr)
    return bool(lb >  0)

def add_mccormick_bound(mip, x3, x1, x2, c1, c2, is_lowerbound):
    """
    Try adding x3 > c1*x1 + c2*x2 - c1*c2 (if is_lowerbound, else x3 < c1*x1 + c2*x2 - c1*c2 to mip. Return True this new constraint is not redundnant.
    """
    if c1 is None or c2 is None:
        # unbounded
        return False
    d = c1.denominator() * c2.denominator() #linear expression needs integer coefficients
    linexpr = Linear_Expression(d*x3 - d*c1*x1 - d*c2*x2 + d*c1*c2)
    if is_lowerbound:
        lb = find_lower_bound_of_linexpr(mip, linexpr)
        if (lb is not None) and (lb >= 0):
            return False
        else:
            mip.add_constraint(linexpr >= 0)
            return True
    else:
        ub = find_upper_bound_of_linexpr(mip, linexpr)
        if (ub is not None) and (ub <=0):
            return False
        else:
            mip.add_constraint(linexpr <= 0)
            return True

def update_mccormicks_for_monomial(m, tightened_mip, original_polyhedron, monomial_list, v_dict, bounds):
    # the argument original_polyhedron is not needed if we assume that
    # recursive McCormicks does not create new monomials.
    # Expect that the monomials in monomial_list have non-decreasing degrees.
    if m.degree() < 2:
        return False
    i = v_dict[m]
    v = Variable(i)
    tightened = False
    for v1 in m.variables():
        i_1 = v_dict[v1]
        v_1 = Variable(i_1)
        lb_1, ub_1 = bounds[i_1]

        v2 = (m/v1).numerator()
        if v2 in v_dict.keys():
            i_2 = v_dict[v2]
        else:
            logging.warn("new monomial %s is needed during recursive McCormick" % v2)
            i_2 = len(monomial_list)
            v_dict[v2]= i_2
            monomial_list.append(v2)
            original_polyhedron.add_space_dimensions_and_embed(1)
            tightened_mip.add_space_dimensions_and_embed(1)
            bounds.append((None, None))
            if update_mccormicks_for_monomial(v2, tightened_mip, original_polyhedron, \
                                              monomial_list, v_dict, bounds):
                tightened = True
        v_2 = Variable(i_2)
        lb_2, ub_2 = bounds[i_2]
        if add_mccormick_bound(tightened_mip, v, v_1, v_2, lb_2, lb_1, True):
            tightened = True
        if add_mccormick_bound(tightened_mip, v, v_1, v_2, ub_2, ub_1, True):
            tightened = True
        if add_mccormick_bound(tightened_mip, v, v_1, v_2, lb_2, ub_1, False):
            tightened = True
        if add_mccormick_bound(tightened_mip, v, v_1, v_2, ub_2, lb_1, False):
            tightened = True
        if m.degree() == 2:
            break
    if tightened:
        bounds[i] = find_bounds_of_variable(tightened_mip, i)
    return tightened

####################################
# Find region type and region color
####################################

def find_region_type_igp(K, h, region_level='extreme', is_minimal=None, use_simplified_extremality_test=True):
    """
    Find the type of a igp function h in the ParametricRealField K;
    (is it constructible? is it minimal? is it extreme?)
    Record the comparisons in K.

    EXAMPLES::

        sage: logging.disable(logging.INFO)
        sage: K.<f> = ParametricRealField([4/5])
        sage: h = gmic(f, field=K)
        sage: find_region_type_igp(K, h)
        'is_extreme'
        sage: K._lt_factor
        {-f, -f + 1/2, f - 1}

        sage: K.<f,bkpt>=ParametricRealField([1/7,3/7])
        sage: h = drlm_backward_3_slope(f, bkpt, field=K)
        sage: find_region_type_igp(K, h)
        'not_extreme'
        sage: K._lt_factor
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
    if h is None:
        return 'not_constructible'
    if region_level == 'constructible':
        return 'is_constructible'
    if is_minimal is None:
        is_minimal = minimality_test(h, stop_if_fail=True)
    if is_minimal:
        if region_level == 'minimal':
            return 'is_minimal'
        if use_simplified_extremality_test:
            is_extreme = simplified_extremality_test(h)
        else:
            is_extreme = extremality_test(h)
        if is_extreme:
            return 'is_extreme'
        else:
            return 'not_extreme'
    else:
        return 'not_minimal'


def result_concrete_value(field, result):
    """
    Return the concrete values in result as a tuple. See also result_symbolic_expression()
 
    This function can provided to find_region_type when setting up SemialgebraicComplex. 
    In this way, one can compare result of type ParametricRealFieldElement or list of ParametricRealFieldElements
    with the previous elements in region_type_color_map which do not necessairy have the same parent.

    EXAMPLES::

        sage: logging.disable(logging.WARN)
        sage: def vol(a,b):
        ....:     P = Polyhedron(ieqs=[(1,0,-1),(0,0,1),(1,-1,0),(0,1,0),(1,-a,-b)])
        ....:     return P.volume()
        sage: K.<a,b>=ParametricRealField([2,3])
        sage: result = vol(a,b)
        sage: result_concrete_value(K, result)
        (1/12,)
    """
    concrete_value = tuple(elt._val if hasattr(elt, '_val') else elt for elt in flatten([result]))
    return concrete_value

def result_symbolic_expression(field, result):
    """
    Return the symbolic expressions in result as a tuple
 
    This function can provided to find_region_type when setting up SemialgebraicComplex. 
    In this way, one can compare result of type ParametricRealFieldElement or list of ParametricRealFieldElements
    with the previous elements in region_type_color_map which do not necessairy have the same parent.

    EXAMPLES::

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

region_type_color_map = [('not_constructible', 'white'), ('not_minimal', 'orange'), ('not_extreme', 'green'), ('is_extreme', 'blue'), (True, 'blue'), (False, 'red')]

def find_region_color(region_type):
    """
    Return the color of the region according to the global dictionary region_type_color_map.

    EXAMPLES::

        sage: find_region_color('is_extreme')
        'blue'
        sage: find_region_color(False)
        'red'
        sage: find_region_color(4)
        '#ff0000'
        sage: find_region_color(10)
        '#0091ff'
    """
    # at most 7 different types other than the ones like 'is_extreme' that were in region_type_color_map.
    global region_type_color_map
    n = len(region_type_color_map)
    i = 0
    while (i < n) and (region_type != region_type_color_map[i][0]):
        i += 1
    if i == n:
        region_color = color_of_ith_region_type(n - 6)  # the initial map includes 4 igp region types and True and False.
        region_type_color_map.append((region_type, region_color))
    else:
        region_color = region_type_color_map[i][1]
    return region_color

def color_of_ith_region_type(i):
    """
    Return a color in the rainbow.
    """
    j = (4 * i) % 7
    c = rainbow(7)[j]
    return c


#######################
# interface mathematica
#######################

def write_mathematica_constraints(eqs, ineqs, strict=True):
    """
    Write polynomial constraints in the mathematica format. 
    Notice that the string ends with ' && '; in practice, often take condstr[:-4]

    EXAMPLES::

        sage: P.<x,y,z>=QQ[]
        sage: eqs = [z]
        sage: ineqs = [-x, x-1, -y, y-1]
        sage: write_mathematica_constraints(eqs, ineqs, strict=True)
        'z == 0 && -x < 0 && x - 1 < 0 && y - 1 < 0 && -y < 0 && '
        sage: write_mathematica_constraints(eqs, ineqs, strict=False)
        'z == 0 && -x <= 0 && x - 1 <= 0 && y - 1 <= 0 && -y <= 0 && '
    """
    condstr = ''
    for l in set(eqs):
        condstr += str(l) + ' == 0 && '
    for l in set(ineqs):
        if strict:
            condstr += str(l) + ' < 0 && '
        else:
            condstr += str(l) + ' <= 0 && '
    return condstr

def write_mathematica_variables(var_name):
    """
    Write the variables in the mathematica format.

    EXAMPLES::

        sage: var_name = ['x','y','z']
        sage: write_mathematica_variables(var_name)
        '{x, y, z}'
    """
    varstr = var_name[0]
    for v in var_name[1::]:
        varstr = varstr + ', ' + v
    return '{' + varstr + '}'

def find_instance_mathematica(condstr, var_name):
    """
    Call the Mathematica's FindInstance to get a point that satisfies the given conditions.

    EXAMPLES::

        sage: condstr = 'z == 0 && -x < 0 && x - 1 < 0 && y - 1 < 0 && -y < 0'
        sage: var_name = ['x','y','z']
        sage: find_instance_mathematica(condstr, var_name)
        (1/2, 1/2, 0)
    """
    varstr =  write_mathematica_variables(var_name)
    pt_math = mathematica.FindInstance(condstr, varstr)
    if len(pt_math) == 0:
        return None
    n = len(var_name)
    pt = []
    for i in range(n):
        try:
            pt_i = QQ(pt_math[1][i+1][2])
        except TypeError:
            pt_i = pt_math[1][i+1][2]
        pt.append(pt_i)
    return tuple(pt)

##########################
# wall crossing heuristic
##########################

def find_point_flip_ineq_heuristic(current_var_value, ineq, ineqs, flip_ineq_step):
    # heuristic method.
    ineq_gradient = gradient(ineq)
    current_point = vector([RR(x) for x in current_var_value]) # Real numbers, faster than QQ
    ineq_value = ineq(*current_point)
    while ineq_value <= 0:
        ineq_direction = vector([g(*current_point) for g in ineq_gradient])
        if ineq.degree() == 1:
            step_length = (-ineq(*current_point)+flip_ineq_step) / (ineq_direction * ineq_direction)
        else:
            step_length = flip_ineq_step / (ineq_direction * ineq_direction) # ineq_value increases by flip_ineq_step=0.01 roughly
            if step_length > 1:
                step_length = 1  # ensure that distance of move <= sqrt(flip_ineq_step) = 0.1 in each step
        current_point += step_length * ineq_direction
        ineq_value = ineq(*current_point)
        #print current_point, RR(ineq_value)
    new_point = adjust_pt_to_satisfy_ineqs(current_point, ineq_gradient, ineqs, flip_ineq_step)
    return new_point #type is tuple


def find_point_on_ineq_heuristic(current_var_value, ineq, ineqs, flip_ineq_step):
    ineq_gradient = gradient(ineq)
    current_point = vector(current_var_value)
    ineq_direction = vector([g(*current_point) for g in ineq_gradient])
    step_length = -ineq(*current_point) / (ineq_direction * ineq_direction)
    current_point += step_length * ineq_direction
    ineq_value = ineq(*current_point)
    new_point = adjust_pt_to_satisfy_ineqs(current_point, ineq_gradient, ineqs, flip_ineq_step)
    return new_point #type is tuple

def adjust_pt_to_satisfy_ineqs(current_point, ineq_gradient, ineqs, flip_ineq_step):
    #current_point is a vector
    for l in ineqs:
        l_gradient = gradient(l)
        l_value = l(*current_point)
        while l_value >= 0:
            l_direction = vector([-g(*current_point) for g in l_gradient]) #decrease l_value
            ineq_direction = vector([g(*current_point) for g in ineq_gradient])
            s = (ineq_direction * l_direction) / (ineq_direction * ineq_direction)
            projected_direction = l_direction - s * ineq_direction # want that ineq_value remains the same
            if projected_direction == 0:
                return None
            if l.degree() == 1:
                step_length = (l_value+flip_ineq_step) / (projected_direction * l_direction)
            else:
                step_length = flip_ineq_step / (projected_direction * l_direction) # l_value decreases by 0.01 roughly
                if step_length * norm(projected_direction) >= 1:  # move too far  # is 1 a good value here?? why this if?
                    return None
            current_point += step_length * projected_direction
            l_value = l(*current_point)
            #print current_point, RR(l_value)
    for l in ineqs:
        if l(*current_point) >= 0:
            return None
    return tuple(QQ(x) for x in current_point)
