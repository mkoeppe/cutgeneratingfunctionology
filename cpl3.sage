# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *
from sage.interfaces.qepcad import _qepcad_atoms
from sage.misc.parser import Parser, Tokenizer

def cpl3_group_function(r0=1/6, z1=1/12, o1=1/5, o2=0, merge=True):
    """
    sage: h = cpl3_group_function()
    sage: extremality_test(h)
    True
    """
    if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if not (bool(0 <= o1) & bool(0 <= o2) & bool(o1+o2 <= 1/2)):
        raise ValueError, "Bad thetas parameters. function value outside [0,1]."

    if z1 < (1-r0)/4:
        bkpts = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
        phi_values = [0, 0, o1, o1+o2, 1-(o1+o2), 1-o1, 1]
    else:
        bkpts = [0, r0, r0+z1, 1-z1, 1]
        phi_values = [0, 0, o1, 1-o1, 1]
    values = [(bkpts[i] - phi_values[i])/r0 for i in range(len(bkpts))]
    #return piecewise_function_from_breakpoints_and_values(bkpts, values, field=field)
    # don't want the comparaison values[i] != values[i+1]
    # in piecewise_function_from_breakpoints_slopes_and_values and in FastPiecewise(merge=True)
    slopes = [(values[i+1]-values[i])/(bkpts[i+1]-bkpts[i]) for i in range(len(bkpts)-1)]
    intercepts = [ values[i] - slopes[i]*bkpts[i] for i in range(len(slopes)) ]
    list_of_pairs = [ [(bkpts[i], bkpts[i+1]), fast_linear_function(slopes[i], intercepts[i])] \
                      for i in range(len(bkpts)-1) ]
    return FastPiecewise(list_of_pairs, periodic_extension=True, merge=merge)

class Cpl3Complex(SageObject):

    def __init__(self, var_name, theta=None, max_iter=8, bddleq=[], bddlin=[]):
        self.components = []
        self.d = len(var_name)
        self.var_name = var_name
        self.monomial_list = []
        self.v_dict = {}
        K = SymbolicRealNumberField([0]*self.d, var_name)
        for i in range(self.d):
            v = K.gens()[i].sym().numerator()
            self.monomial_list.append(v)
            self.v_dict[v] = i
        self.graph = Graphics()
        self.num_plotted_components = 0
        self.points_to_test = set()
        self.bddleq_of_testpoint={}
        self.theta = theta
        self.max_iter = max_iter
        self.bddleq = bddleq
        self.bddlin = bddlin


    def is_point_covered(self, var_value):
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
        return False

    def add_new_component(self, var_value, bddleq, flip_ineq_step=-1/100):
        # Remark: the sign of flip_ineq_step indicates how to search for neighbour testpoints:
        # if flip_ineq_step = 0, don't search for neighbour testpoints. Used in shoot_random_points().
        # if flip_ineq_step < 0, we assume that the walls of the cell are linear eqn/ineq over original parameters. (So, gradient is constant; easy to find a new testpoint on the wall and another testpoint (-flip_ineq_step away) across the wall.) Used in bfs.
        # if flip_ineq_step > 0, we don't assume the walls are linear. Apply generate_one_point_by_flipping_inequality() with flip_ineq_step to find new testpoints across the wall only. Used in bfs.
        unlifted_space_dim =  len(self.monomial_list)
        K = SymbolicRealNumberField(var_value, self.var_name)
        K.monomial_list = self.monomial_list # change simultaneously while lifting
        K.v_dict = self.v_dict # change simultaneously while lifting
        K.polyhedron.add_space_dimensions_and_embed(len(K.monomial_list))
        r0, z1 = K.gens()[0], K.gens()[1]
        for l in self.bddlin:
            if not l(*K.gens()) < 0:
                return
        for l in bddleq:
            if not l(*K.gens()) == 0:
                return
        (r0mapping, z1mapping) = mapping_r0_z1(bddleq)
        r0m = r0mapping(r0, z1) * K.one(); z1m = z1mapping(r0, z1) * K.one()
        if self.theta is None:
            try:
                o1 = o2 = z1 / (1 - r0)
                h = cpl3_group_function(r0, z1, o1, o2, merge=False)
                subadditivity_test(h) # always True
                region_type = 'is_constructible'
            except:
                region_type = 'not_constructible'
        else:
            try:
                o1 = self.theta[0](r0m, z1m)
                o2 = self.theta[1](r0m, z1m)
                h = cpl3_group_function(r0m, z1m, o1, o2)
            except:
                h = None
            region_type =  find_region_type_around_given_point(K, h)
        new_component = SemialgebraicComplexComponent(self, K, var_value, region_type)
        if new_component.leq and max([l.degree() for l in new_component.leq]) > 1:
            logging.warn("The cell around testpoint %s defined by %s leqs and %s lins has higher order equations." % (new_component.var_value, new_component.leq, new_component.lin)) # Might be problem to variable elemination.
        # TODO: call qepcad to get a minimal discription.
        # ignore equations for now to keep number of variables small.
        if new_component.lin and max([l.degree() for l in new_component.lin]) > 1:
            new_component.lin = remove_redundancy_using_qepcad(new_component.lin)
        # put non-linear eq/ineq to somewhere else.
        if flip_ineq_step < 0:
            new_component.nleq = [l for l in new_component.leq if l.degree()>1]
            new_component.leq = [l for l in new_component.leq if l.degree()<=1]
            new_component.nlin = [l for l in new_component.lin if l.degree()>1]
            new_component.lin = [l for l in new_component.lin if l.degree()<=1]
        #if see new monomial, lift polyhedrons of the previously computed components.
        dim_to_add = len(self.monomial_list) - unlifted_space_dim
        if dim_to_add > 0:
            for c in self.components:
                c.polyhedron.add_space_dimensions_and_embed(dim_to_add)
        self.components.append(new_component)
        neighbour_points = new_component.generate_neighbour_points(flip_ineq_step, self.bddlin)
        if (flip_ineq_step > 0):
            (self.points_to_test).update(neighbour_points)
        elif (flip_ineq_step < 0):
            for (new_point, new_bddleq) in neighbour_points:
                (self.points_to_test).add(new_point)
                (self.bddleq_of_testpoint)[new_point] = bddleq + new_bddleq

    def plot(self, alpha=0.5, plot_points=300, slice_value=None, restart=False):
        if restart:
            self.graph = Graphics()
            self.num_plotted_components = 0
        for c in self.components[self.num_plotted_components::]:
            self.graph += c.plot(alpha=alpha, plot_points=plot_points, slice_value=slice_value)
        self.num_plotted_components = len(self.components)
        return self.graph

    def bfs_completion(self, var_value=None, flip_ineq_step=-1/100):
        # See remark about flip_ineq_step in def add_new_component().
        if var_value:
            (self.points_to_test).add(tuple(var_value))
            self.bddleq_of_testpoint[tuple(var_value)] = copy(self.bddleq)
        while self.points_to_test:
            var_value = list(self.points_to_test.pop())
            bddleq = self.bddleq_of_testpoint[tuple(var_value)]
            if not self.is_point_covered(var_value):
                self.add_new_component(var_value, bddleq, flip_ineq_step=flip_ineq_step)

def remove_redundancy_using_qepcad(lins):
    """
    sage: PR2.<r0,z1>=QQ[]
    sage: lins = [r0 + 8*z1 - 1, -z1, -r0*z1 - z1, r0^2 + r0*z1 - r0 + z1, r0*z1 - r0 - z1, -2*r0 + 1, r0^2*z1 - 2*r0*z1 - z1]
    sage: remove_redundancy_using_qepcad(lins)
    [r0 + 8*z1 - 1, -2*r0 + 1, -z1]
    sage: lins = [r0 + 4*z1 - 1, -r0 - 8*z1 + 1, -2*r0 + 1, 4*r0*z1 - 3*r0 + 4*z1 - 1, -4*r0*z1 + 4*z1 - 1]
    sage: remove_redundancy_using_qepcad(lins)
    [-r0 - 8*z1 + 1, -2*r0 + 1, r0 + 4*z1 - 1]
    """
    # TODO: check on the qepcad option measure-zero-error
    # which is supposed to bring huge reduction in time and space.
    # lins are lists of Multivariate Polynomial. (we ignored leqs)
    qf = qepcad_formula
    qe = qepcad([qf.atomic(l, operator.lt) for l in lins], memcells=50000000)
    # note: type(qe) = <class 'sage.interfaces.interface.AsciiArtString'>
    # note: with interact=True, call qe.finish() and qe.answer()
    #qe_interact = qepcad([qf.atomic(l, operator.lt) for l in lins], memcells=50000000, interact=True)
    #qe_interact.finish()
    #qe = qe_interact.answer()
    if ("\\/" in qe):
        logging.warn("The qepcad output of %s is %s, which has \\/." % (lins, qe))
    # see how many inequalities qepcad managed to reduce.
    #if (qe.count("/\\")+1 != len(lins)):
    #    print len(lins), lins
    #    print qe.count("/\\")+1, qe
    # not complete....to parse qe
    # assume there is no disjuction.
    atomset = _qepcad_atoms(qe)
    p = Parser(make_var=var)
    # assume rhs is 0, operator is either < or >
    new_lins = []
    for a in atomset:
        q = p.p_eqn(Tokenizer(a))
        if q.operator() == operator.lt:
            qexpr = q.lhs()-q.rhs()
        elif q.operator() == operator.gt:
            qexpr = q.rhs()-q.lhs()
        else:
            raise ValueError, "The min description of %s is %s. It has atom %s which is not a strict inequality." % (lins, qe, a)
        qpoly = QQ['r0, z1'](qexpr)
        new_lins.append(qpoly)
    del(qe) # useful?
    return new_lins

def bddlin_cpl():
    """
    boundary of (r0, z1) in cpl3: r0>0, z1>0, r0+4*z1<1.
    """
    K.<r0,z1>=QQ[]
    return [-r0, -z1, r0+4*z1-1]

def regions_r0_z1_from_arrangement_of_bkpts(max_iter=8):
    """
    Got regions[0:30]: 2-dim; regions[30:73]: 1-dim; regions[73:87]: 0-dim.

    sage: logging.disable(logging.INFO)  # not tested
    sage: regions = regions_r0_z1_from_arrangement_of_bkpts() # not tested
    sage: len(regions) #not tested
    87
    """
    arr_complex=Cpl3Complex(['r0','z1'], theta=None, bddlin=bddlin_cpl(), max_iter=max_iter)
    arr_complex.bfs_completion(var_value=[6/10,4/100])
    regions = arr_complex.components
    regions.sort(key=lambda r: len(r.leq))
    return regions

def mapping_r0_z1(leq):
    """
    sage: PR2.<r0,z1>=QQ[]
    sage: mapping_r0_z1([])
    (r0, z1)
    sage: mapping_r0_z1([-2*r0 + 1])
    (1/2, z1)
    sage: mapping_r0_z1([-r0 - 6*z1 + 1])
    (r0, -1/6*r0 + 1/6)
    sage: mapping_r0_z1([-12*z1 + 1, -2*r0 + 1])
    (1/2, 1/12)
    """
    PR2.<r0,z1>=QQ[]
    if not leq:
        return (r0, z1)
    elif len(leq) == 1:
        l = leq[0]
        c_r0 = l.monomial_coefficient(r0)
        c_z1 = l.monomial_coefficient(z1)
        if c_z1 != 0:
            return (r0, z1 - l / c_z1)
        elif c_r0 != 0:
            return (r0 - l / c_r0, z1)
        else:
            raise ValueError
    elif len(leq) == 2:
        l1 = leq[0]; l2 = leq[1]
        (a11, a12, b1) = (l1.monomial_coefficient(r0), l1.monomial_coefficient(z1), -l1.constant_coefficient())
        (a21, a22, b2) = (l2.monomial_coefficient(r0), l2.monomial_coefficient(z1), -l2.constant_coefficient())
        determinant = a11 * a22 - a12 * a21
        if determinant == 0:
            if not a11 == a12 == 0:
                return mapping_r0_z1([l1])
            else:
                return mapping_r0_z1([l2])
        r0_val = (b1 * a22 - a12 * b2) / determinant
        z1_val = (a11 * b2 - b1 * a21) / determinant
        return (r0_val * PR2.one(), z1_val * PR2.one())

def symbolic_subbadditivity_constraints_of_cpl3_given_region(r):
    """
    sage: regions = regions_r0_z1_from_arrangement_of_bkpts()
    sage: r = regions[0]
    sage: symbolic_subbadditivity_constraints_of_cpl3_given_region(r)
    [1/r0,
     (-o1 + 1)/r0,
     (-o1 - o2 + 1)/r0,
     (o1 + o2)/r0,
     o1/r0,
     (-2*o1 + 1)/r0,
     (-2*o1 - o2 + 1)/r0,
     o2/r0,
     (-2*o1 - 2*o2 + 1)/r0,
     (o1 - o2)/r0,
     (r0*o1 + r0*o2 + 2*z1 - o1 - o2)/(r0^2 + 4*r0*z1 - r0),
     (r0*o1 + 2*z1*o1 - 2*z1*o2 + z1 - o1)/(r0^2 + 4*r0*z1 - r0)]
    """
    if r.region_type == 'not_constructible':
        return []
    K.<r0,z1,o1,o2>=SymbolicRealNumberField([r.var_value[0], r.var_value[1], 0, 0])
    r0s = r0.sym(); z1s = z1.sym(); o1s = o1.sym(); o2s = o2.sym()
    (r0mapping, z1mapping) = mapping_r0_z1(r.leq)
    r0m = r0mapping(r0s, z1s); z1m = z1mapping(r0s, z1s)

    h = cpl3_group_function(r0, z1, o1, o2, merge=False)
    bkpts = h.end_points()
    bkpts2 = bkpts[:-1] + [ x+1 for x in bkpts]
    constraints = [] # can't use set becasue of the difference between '==' and 'is'
    for x in bkpts:
        for y in bkpts:
            if 0 < x <= y < 1:
                c = delta_pi(h, x, y).sym()(r0m, z1m, o1s, o2s)
                if (c not in QQ) and all(c != cc for cc in constraints):
                    constraints.append(c)
    for x in bkpts:
        for z in bkpts2:
            if 0 < x < z < 1+x:
                c = delta_pi(h, x, z-x).sym()(r0m, z1m, o1s, o2s)
                if (c not in QQ) and all(c != cc for cc in constraints):
                    constraints.append(c)
    return constraints

def coefficients_of_linear_constraint_in_thetas(c):
    """
    sage: PR4.<r0,z1,o1,o2>=QQ[]
    sage: c = (r0*o1 + r0*o2 + 2*z1 - o1 - o2)/(r0^2 + 4*r0*z1 - r0)
    sage: coefficients_of_linear_constraint_in_thetas(c)
    ((r0 - 1)/(r0^2 + 4*r0*z1 - r0),
     (r0 - 1)/(r0^2 + 4*r0*z1 - r0),
     (-2*z1)/(r0^2 + 4*r0*z1 - r0))
    sage: c = (r0*o1 + 2*z1*o1 - 2*z1*o2 + z1 - o1)/(r0^2 + 4*r0*z1 - r0)
    sage: coefficients_of_linear_constraint_in_thetas(c)
    ((r0 + 2*z1 - 1)/(r0^2 + 4*r0*z1 - r0),
     (-2*z1)/(r0^2 + 4*r0*z1 - r0),
     (-z1)/(r0^2 + 4*r0*z1 - r0))
    """
    # assume c.parent() is Fraction Field of Multivariate Polynomial Ring in r0, z1, o1, o2 over Rational Field
    # assume c is linear wrt o1 and o2.
    PR4.<r0,z1,o1,o2>=QQ[]
    cn = c.numerator(); cd = c.denominator()
    c_o1 = cn.coefficient(o1) / cd
    c_o2 = cn.coefficient(o2) / cd
    c_rhs = - cn.mod([o1,o2]) / cd
    PR2.<r0,z1>=QQ[]
    return (c_o1(r0, z1, 0, 0), c_o2(r0, z1, 0, 0), c_rhs(r0, z1, 0, 0))

def symbolic_subadditivity_constraint_coefficients_of_thetas_given_region(r):
    """
    sage: regions = regions_r0_z1_from_arrangement_of_bkpts()
    sage: r = regions[0]
    sage: coefficients_of_thetas = symbolic_subadditivity_constraint_coefficients_of_thetas_given_region(r)
    sage: coefficients_of_thetas[-1]
    ((r0 + 2*z1 - 1)/(r0^2 + 4*r0*z1 - r0),
     (-2*z1)/(r0^2 + 4*r0*z1 - r0),
     (-z1)/(r0^2 + 4*r0*z1 - r0))
    """
    constraints = symbolic_subbadditivity_constraints_of_cpl3_given_region(r)
    coefficients_of_thetas = []
    for c in constraints:
        (c_o1, c_o2, c_rhs) = coefficients_of_linear_constraint_in_thetas(c)
        if not (c_o1 == c_o2 == 0):
            coefficients_of_thetas.append((c_o1, c_o2, c_rhs))
    return coefficients_of_thetas

def get_pretty_fraction_polynomial(fp):
    """
    sage: PR2.<r0,z1>=QQ[]
    sage: fp = (-2*z1)/(-2*r0 + 4*z1 - 6)
    sage: get_pretty_fraction_polynomial(fp)
    z1/(r0 - 2*z1 + 3)
    sage: fp = (r0/4+z1/3)/(1/12)
    sage: get_pretty_fraction_polynomial(fp)
    3*r0 + 4*z1
    """
    pn = fp.numerator()
    pd = fp.denominator()
    if type(pn)==sage.rings.rational.Rational or type(pn)==sage.rings.integer.Integer:
        gcd_pn = abs(pn)
    else:
        gcd_pn = gcd(pn.coefficients())
    if type(pd)==sage.rings.rational.Rational or type(pd)==sage.rings.integer.Integer:
        gcd_pd = abs(pd)
    else:
        gcd_pd = gcd(pd.coefficients())
    gcd_n_d = gcd([gcd_pn, gcd_pd])
    if type(pd)==sage.rings.rational.Rational or type(pd)==sage.rings.integer.Integer:
        if pd < 0:
            gcd_n_d = - gcd_n_d
    elif pd.lc() < 0:
        gcd_n_d = - gcd_n_d
    pretty_fp = (pn/gcd_n_d)/(pd/gcd_n_d)
    assert pretty_fp == fp
    return pretty_fp

def generate_thetas_of_region(r):
    """
    sage: regions = regions_r0_z1_from_arrangement_of_bkpts()
    sage: r = regions[0]
    sage: thetas = generate_thetas_of_region(r)
    sage: len(thetas)
    24
    sage: thetas[-1]
    ((-z1)/(r0 - 1), (-z1)/(r0 - 1))
    """
    thetas = []  # can't use set becasue of the difference between '==' and 'is'
    coefficients_of_thetas = symbolic_subadditivity_constraint_coefficients_of_thetas_given_region(r)
    n = len(coefficients_of_thetas)
    for i in range(n):
        for j in range(i+1, n):
            a11, a12, b1 = coefficients_of_thetas[i]
            a21, a22, b2 = coefficients_of_thetas[j]
            determinant = a11 * a22 - a12 * a21
            if determinant == 0:
                continue
            theta1 = (b1 * a22 - a12 * b2) / determinant
            theta2 = (a11 * b2 - b1 * a21) / determinant
            theta = (get_pretty_fraction_polynomial(theta1), get_pretty_fraction_polynomial(theta2))
            if all(theta != t for t in thetas):
                thetas.append(theta)
    return thetas

def fill_region_given_theta(r, theta, max_iter=8):
    """
    sage: regions = regions_r0_z1_from_arrangement_of_bkpts()
    sage: r = regions[0]
    sage: thetas = generate_thetas_of_region(r)
    sage: theta = thetas[-1]
    sage: theta
    ((-z1)/(r0 - 1), (-z1)/(r0 - 1))
    sage: cpl_complex = fill_region_given_theta(r, theta)
    sage: len(cpl_complex.components)
    1
    sage: cpl_complex.components[0].region_type
    'is_extreme'
    """
    cpl_complex = Cpl3Complex(['r0','z1'], theta=theta, bddleq=copy(r.leq), bddlin=copy(r.lin), max_iter=max_iter)
    cpl_complex.bfs_completion(var_value=tuple(r.var_value), flip_ineq_step=-1/100)
    return cpl_complex


def cpl_regions_with_thetas_and_components(keep_extreme_only=False):
    """
    sage: regions = cpl_regions_with_thetas_and_components()
    output warning:
    0 [3/5, 1/25]
    ...
    9 [243/520, 6/325]
    WARNING: 2016-02-17 07:44:25,836 max number 8 of bounds propagation iterations has attained.
    ...
    24 [319/2000, 169/1000]
    WARNING: 2016-02-17 07:46:44,986 crossed non-linear wall 2*r0^2 + 9*r0*z1 + 7*z1^2 - 3*r0 - 6*z1 + 1 < 0 while flipping -4*r0 - 14*z1 + 3 < 0
    ...
    86 [1/8, 1/8]
    """
    regions = regions_r0_z1_from_arrangement_of_bkpts()
    for i in range(len(regions)):
        r = regions[i]
        print i, r.var_value
        r.thetas = {}
        thetas_of_r = generate_thetas_of_region(r)
        for theta in thetas_of_r:
            cpl_complex = fill_region_given_theta(r, theta)
            if keep_extreme_only:
                components = [c for c in cpl_complex.components if c.region_type=='is_extreme']
            else:
                components = cpl_complex.components
            if components:
                r.thetas[theta] = components
    return regions

def cpl_thetas_and_regions_extreme(regions):
    """
    sage: regions = cpl_regions_with_thetas_and_components()
    sage: thetas_and_regions = cpl_thetas_and_regions_extreme(regions)
    sage: len(thetas_and_regions)
    9
    """
    thetas_and_regions = {}
    for r in regions:
        (r0m, z1m) = mapping_r0_z1(r.leq)
        for (theta, components) in (r.thetas).items():
            extreme_regions = [c for c in components if c.region_type=='is_extreme']
            if not extreme_regions:
                continue
            new_theta = True
            for t in thetas_and_regions.keys():
                try:
                    if theta == (t[0](r0m, z1m), t[1](r0m, z1m)):
                        thetas_and_regions[t] += extreme_regions
                        new_theta = False
                except ZeroDivisionError:
                    continue
            if new_theta:
                thetas_and_regions[theta] = extreme_regions
    return thetas_and_regions

def cpl_regions_fix_theta(regions, theta):
    components = []
    for r in regions:
        (r0m, z1m) = mapping_r0_z1(r.leq)
        try:
            tt =  (theta[0](r0m, z1m), theta[1](r0m, z1m))
        except ZeroDivisionError:
            continue
        for (t, ct) in (r.thetas).items():
            if t == tt:
                components += ct
    return components

def plot_cpl_components(components, show_testpoints=False):
    """
    sage: regions = regions_r0_z1_from_arrangement_of_bkpts()
    sage: g = plot_cpl_components(regions)
    sage: g.show(xmin=0, xmax=1, ymin=0, ymax=1/4) #not tested
    """
    covered_type_color = {'is_constructible': 'black', 'not_constructible': 'white', 'not_minimal': 'orange', 'not_extreme': 'green', 'is_extreme': 'blue'}
    g = Graphics()
    for c in components:
        if not c.leq:
            g += c.plot(show_testpoints=show_testpoints)
    for c in components:
        if len(c.leq)==1:
            g += c.plot(show_testpoints=show_testpoints)
    for c in components:
        if len(c.leq)==2:
            ptcolor = covered_type_color[c.region_type]
            g += point(c.var_value, color = ptcolor, zorder=10)
    return g

def save_cpl_extreme_theta_regions(thetas_and_regions):
    """
    To plot only blue regions
    sage: regions = cpl_regions_with_thetas_and_components()
    sage: thetas_and_regions = cpl_thetas_and_regions_extreme(regions)
    sage: save_cpl_extreme_theta_regions(thetas_and_regions)
    sage: len(thetas_and_regions)
    9

    To plot colorful regions corresponding to each extreme thetas.
    sage: thetas_and_components = {}
    sage: for theta in thetas_and_regions.keys():
    ....:     components = cpl_regions_fix_theta(regions, theta)
    ....:     thetas_and_components[theta]=components
    sage: save_cpl_extreme_theta_regions(thetas_and_components)
    """
    k = 0
    for (theta, components) in thetas_and_regions.items():
        k += 1
        print (k, theta)
        if len(str(theta)) < 100:
            title = "extreme point %s:  theta = %s" % (k, theta)
        else:
            title = "extreme point %s:\ntheta = (%s,\n %s)" % (k, theta[0], theta[1])
        g = plot_cpl_components(components)
        g += text(title, (0.5, 1/4), color='black')
        g.save("cpl_theta_%s.png" %k, xmin=0, xmax=1, ymin=0, ymax=1/4)
