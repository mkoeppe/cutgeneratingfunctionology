# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

import sage.structure.element
from sage.structure.element import FieldElement


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
            logging.info("New element in %s._eq: %s" % (repr(self), comparison))
            self._eq.add(comparison)
            self.record_poly(comparison.numerator())
            self.record_poly(comparison.denominator())
    def record_to_lt_list(self, comparison):
        if not comparison in QQ and not comparison in self._lt:
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
            self._eq_poly.add(poly)
            for (fac, d) in poly.factor():
                # record the factor if it's zero
                if fac(self._values) == 0 and not fac in self._eq_factor:
                    self._eq_factor.add(fac)
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
                        self._lt_factor.add(new_fac)

default_symbolic_field = SymbolicRealNumberField()


###############################
# Simplify polynomials
###############################

from sage.libs.ppl import Variable, Constraint, Linear_Expression, Constraint_System, NNC_Polyhedron

def polynomial_to_linexpr(t, monomial_list, v_dict):
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
    p = NNC_Polyhedron(cs)
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
    if level == "factor":
        leq, lin = simplify_eq_lt_poly_via_ppl(K.get_eq_factor(), K.get_lt_factor())
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
                return lambda x, y: function.check_conditions(**dict({var_name[0]: x, var_name[1]: y}, **test_point)) != 'not_constructible'
            else:
                return lambda x, y: function.check_conditions(**dict({var_name[0]: x}, **test_point )) != 'not_constructible'
        elif region_level == 'extreme':
            if len(var_name) == 2:
                return lambda x, y: function.check_conditions(**dict({var_name[0]: x, var_name[1]: y}, **test_point)) == 'extreme'
            else:
                return lambda x, y: function.check_conditions(**dict({var_name[0]: x}, **test_point )) == 'extreme'
        else:
            raise ValueError, "Bad argument region_level = %s" % region_level

def find_region_around_given_point(K, h, level="factor", region_level='extreme', is_minimal=None, use_simplified_extremality_test=True):
    ## Note: region_level = 'constructible' / 'minimal'/ 'extreme'
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
    sage: plot_region_according_to_literature(drlm_backward_3_slope, ['f','bkpt'], plot_points=1000)
    sage: plot_region_according_to_literature(gj_forward_3_slope, ['lambda_1', 'lambda_2'], f = 2/3, plot_points=500)
    sage: plot_region_according_to_literature(dg_2_step_mir, ['f','alpha'], plot_points=500)
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
    sage: logging.disable(logging.INFO)
    sage: g, fname, fparam = plot_parameter_region(drlm_backward_3_slope, ['f','bkpt'], [1/12-1/30, 2/12])
    sage: g.save(fname+".pdf", title=fname+fparam, legend_loc=7)

    sage: plot_parameter_region(drlm_backward_3_slope, ['f','bkpt'], [1/12+1/30, 2/12], plot_points=500)
    sage: plot_parameter_region(drlm_backward_3_slope, ['f','bkpt'], [1/10-1/100, 3/10], plot_points=500)

    sage: plot_parameter_region(gj_2_slope, ['f','lambda_1'], plot_points=100)
    sage: plot_parameter_region(gj_2_slope, ['f','lambda_1'], [4/5, 4/6], plot_points=100)

    sage: plot_parameter_region(gj_forward_3_slope, ['f','lambda_1'], [9/10, 3/10])
    sage: plot_parameter_region(gj_forward_3_slope, ['lambda_1','lambda_2'])
    sage: plot_parameter_region(gj_forward_3_slope, ['lambda_1','lambda_2'], f=3/4, lambda_1=2/3, lambda_2=3/5)
    sage: plot_parameter_region(gj_forward_3_slope, ['lambda_1','lambda_2'], [3/5, 7/5], f=1/2, lambda_1=4/9, lambda_2=1/3)
    sage: plot_parameter_region(gj_forward_3_slope, ['lambda_1','lambda_2'], [9/10, 3/7], f=1/2, lambda_1=4/9, lambda_2=1/3)

    sage: plot_parameter_region(chen_4_slope, ['lam1', 'lam2'], [1/4, 1/4], plot_points=100)

    sage: plot_parameter_region(dg_2_step_mir, ['f','alpha'], [3/5,2/5-1/23],  plot_points=500)
    sage: plot_parameter_region(dg_2_step_mir_limit, ['f'], [8/13], d=2, plot_points=500)

    sage: plot_parameter_region(kf_n_step_mir, ['f'], [1/13], plot_points=100)
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
    sage: regions, tested_points, last_n_points_was_covered = cover_parameter_region(gmic, ['f'], random_points=10, plot_points=100)
    sage: _ = cover_parameter_region(drlm_backward_3_slope, ['f','bkpt'], random_points=200, plot_points=500)
    sage: _ = cover_parameter_region(gj_2_slope, ['f', 'lambda_1'], 300, plot_points=100, show_points=1)

    sage: _ = cover_parameter_region(gj_forward_3_slope, ['lambda_1', 'lambda_2'], random_points=500, plot_points=500)
    sage: _ = cover_parameter_region(gj_forward_3_slope, ['f', 'lambda_1'], random_points=100, plot_points=500)

    sage: regions, tested_points, last_n_points_was_covered = cover_parameter_region(chen_4_slope, ['lam1', 'lam2'], 20, plot_points=0)
    sage: _ = cover_parameter_region(chen_4_slope, ['lam1', 'lam2'], 20, plot_points=500, show_points=1, regions=regions, tested_points=tested_points)
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
