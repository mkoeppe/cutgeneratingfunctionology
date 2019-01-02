from six.moves import range
def cpl_n_group_function(n, cpleq=False, merge=True):
    return CPLFunctionsFactory(n, cpleq, merge)

class CPLFunctionsFactory:
    """
    A Factory of CPL functions.

    EXAMPLES::

        sage: logging.disable(logging.INFO)
        sage: cpl3 = cpl_n_group_function(3)
        sage: h = cpl3(f=1/6, z=(1/12, 1/12), o=(1/5, 0))
        sage: extremality_test(h)
        True

    Omit the argument o, then it takes the mapping self._theta.
    With the default value of self._theta, it creates a gmic function::

        sage: g = cpl3(1/6, (1/12, 1/12))
        sage: g == gmic(1/6)

    Change self._theta::

        sage: cpl3._theta = lambda f, z: (1/5, 0)
        sage: hh = cpl3(1/6, (1/12, 1/12))
        sage: h == hh
        True
    """
    def __init__(self, n, cpleq=False, merge=True):
        self._n = n
        self._cpleq = cpleq # in cpl3eq, z1=z2
        if not cpleq:
            self._theta = lambda f, z: tuple([z[i]/(1-f) for i in range(self._n-1)])
        else:
            self._theta = lambda f, z: tuple([z[0]/(1-f) for i in range(self._n-1)])
        self._merge = merge # need  merge=False in cpl_regions_from_arrangement_of_bkpts
    def __call__(self, f, z, o=None, field=None):
        if o is None:
            o = self._theta(f, z)
        if not self._cpleq:
            zz = z
        else:
            zz = [z[0] for i in range(self._n-1)]
        if not (bool(0 < f < 1) & bool(all(0 < zi for zi in zz)) & bool (sum(zz) <= (1-f)/2)):
            raise ValueError("Bad parameters. Unable to construct the function.")
        if not (bool(0 <= oi for oi in o) & bool(sum(o) <= 1/2)):
             raise ValueError("Bad thetas parameters. function value outside [0,1].")
        if sum(zz) < (1-f)/2:
            m = self._n
        else:
            m = self._n - 1 
        bkpts = [0] + [f + sum(zz[0:i]) for i in range(m)] + \
                [1 - sum(zz[0:i]) for i in range(m - 1, -1, -1)]
        phi_values = [0] + [sum(o[0:i]) for i in range(m)] + \
                     [1 - sum(o[0:i]) for i in range(m - 1, -1, -1)]
        values = [(bkpts[i] - phi_values[i])/f for i in range(len(bkpts))]
        return piecewise_function_from_breakpoints_and_values(bkpts, values, field=field, merge=self._merge)

def cpl_regions_from_arrangement_of_bkpts(n=3, cpleq=True, max_iter=0, flip_ineq_step=1/1000, check_completion=False, wall_crossing_method='heuristic', goto_lower_dim=True, max_failings=0):
    """
    Got regions[0:30]: 2-dim; regions[30:78]: 1-dim; regions[78:96]: 0-dim.

    EXAMPLES::

        sage: logging.disable(logging.INFO)                                  # not tested
        sage: regions = cpl_regions_from_arrangement_of_bkpts(3, cpleq=True) # not tested
        sage: len(regions)                                                   # not tested
        96
        sage: regions = cpl_regions_from_arrangement_of_bkpts(n=3, cpleq=True, wall_crossing_method='mathematica', flip_ineq_step=1/1000)                          # not tested
        sage: len(regions)                                                   # not tested
        96
    """
    cpln = cpl_n_group_function(n, cpleq, merge=False)
    if cpleq:
        var_name = ['f', 'z']
        var_value = [3/5, 3/25/n] #arbitary values
    else:
        var_name = ['f']+['z%s' % i for i in range(1, n)]
        var_value = [3/5]+[1/(5^i) for i in range(2, n+1)]  #arbitary values
    def frt(K, h):
        return find_region_type_igp(K, h, region_level='minimal', is_minimal=None)
    arr_complex=SemialgebraicComplex(cpln, var_name, max_iter=max_iter, find_region_type=frt)
    if flip_ineq_step != 0:
        arr_complex.bfs_completion(var_value, \
                                   flip_ineq_step=flip_ineq_step, \
                                   check_completion=check_completion, \
                                   wall_crossing_method=wall_crossing_method, \
                                   goto_lower_dim=goto_lower_dim, \
                                   max_failings = max_failings)
    else:
        arr_complex.shoot_random_points(max_failings)
    regions = [c for c in arr_complex.components if c.region_type != 'not_constructible']
    regions.sort(key=lambda r: len(r.leq))
    return regions

def plot_cpl_components(components, show_testpoints=False):
    """
    EXAMPLES::

        sage: regions = cpl_regions_from_arrangement_of_bkpts(3, cpleq=True) # not tested
        sage: g = plot_cpl_components(regions)                               # not tested
        sage: g.show(xmin=0, xmax=1, ymin=0, ymax=1/4)                       # not tested
    """
    g = Graphics()
    for c in components:
        if not c.leq:
            g += c.plot(show_testpoints=show_testpoints)
    for c in components:
        if len(c.leq)==1:
            g += c.plot(show_testpoints=show_testpoints)
    for c in components:
        if len(c.leq)==2:
            ptcolor = find_region_color(c.region_type)
            g += point(c.var_value, color = ptcolor, zorder=10)
    return g

def symbolic_subbadditivity_constraints_of_cpl_given_region(r):
    """
    EXAMPLES::

        sage: regions = cpl_regions_from_arrangement_of_bkpts(3, cpleq=True) # not tested
        sage: r = regions[0]                                                 # not tested
        sage: symbolic_subbadditivity_constraints_of_cpl_given_region(r)     # not tested
        [1/f,
         (-o1 + 1)/f,
         (-o1 - o2 + 1)/f,
         (o1 + o2)/f,
         o1/f,
         (-2*o1 + 1)/f,
         (-2*o1 - o2 + 1)/f,
         o2/f,
         (-2*o1 - 2*o2 + 1)/f,
         (o1 - o2)/f,
         (f*o1 + f*o2 + 2*z - o1 - o2)/(f^2 + 4*f*z - f),
         (f*o1 + 2*z*o1 - 2*z*o2 + z - o1)/(f^2 + 4*f*z - f)]
    """
    if r.region_type == 'not_constructible':
        return []
    cpln = r.parent.function
    n = cpln._n
    var_value = r.var_value
    var_name = r.parent.var_name
    o_str = ['o%s' % i for i in range(1, n)]
    K = ParametricRealField(var_value+[0]*(n-1), var_name+o_str)
    f_K = K.gens()[0]
    z_K = K.gens()[1:len(var_name)]
    o_K = K.gens()[len(var_name)::]
    fz = [v._sym for v in [f_K]+z_K]
    P = PolynomialRing(QQ, var_name)
    var_sym = [r.var_map[v](*fz) for v in P.gens()] + [v._sym for v in o_K]
    h = cpln(f_K, z_K, o_K)
    bkpts = h.end_points()
    bkpts2 = bkpts[:-1] + [ x+1 for x in bkpts]
    constraints = [] # can't use set becasue of the difference between '==' and 'is'
    for x in bkpts:
        for y in bkpts:
            if 0 < x <= y < 1:
                c = delta_pi(h, x, y).sym()(var_sym)
                if (c not in QQ) and all(c != cc for cc in constraints):
                    constraints.append(c)
    for x in bkpts:
        for z in bkpts2:
            if 0 < x < z < 1+x:
                c = delta_pi(h, x, z-x).sym()(var_sym)
                if (c not in QQ) and all(c != cc for cc in constraints):
                    constraints.append(c)
    return constraints

def coefficients_of_theta_and_rhs_in_constraint(c, var_name, n):
    """
    EXAMPLES::

        sage: n = 3
        sage: var_name = ['f','z1','z2']
        sage: P.<f,z1,z2,o1,o2>=QQ[]
        sage: c = (f*o1 + f*o2 + z1 + z2 - o1 - o2)/(f^2 + 2*f*z1 + 2*f*z2 - f)
        sage: coefficients_of_theta_and_rhs_in_constraint(c, var_name, n)
        (((f - 1)/(f^2 + 2*f*z1 + 2*f*z2 - f), (f - 1)/(f^2 + 2*f*z1 + 2*f*z2 - f)),
         (-z1 - z2)/(f^2 + 2*f*z1 + 2*f*z2 - f))
        sage: c = (f*o1 + 2*z2*o1 - 2*z1*o2 + z1 - o1)/(f^2 + 2*f*z1 + 2*f*z2 - f)
        sage: coefficients_of_theta_and_rhs_in_constraint(c, var_name, n)
        (((f + 2*z2 - 1)/(f^2 + 2*f*z1 + 2*f*z2 - f),
          (-2*z1)/(f^2 + 2*f*z1 + 2*f*z2 - f)),
         (-z1)/(f^2 + 2*f*z1 + 2*f*z2 - f))

        sage: var_name = ['f','z']
        sage: P.<f,z,o1,o2>=QQ[]
        sage: c = (f*o1 + 2*z*o1 - 2*z*o2 + z - o1)/(f^2 + 4*f*z - f)
        sage: coefficients_of_theta_and_rhs_in_constraint(c, var_name, n)
        (((f + 2*z - 1)/(f^2 + 4*f*z - f), (-2*z)/(f^2 + 4*f*z - f)),
         (-z)/(f^2 + 4*f*z - f))
    """
    # assume c.parent() is Fraction Field of Multivariate Polynomial Ring in var_names and (n-1) theta variables over Rational Field.
    # assume c is linear over the (n-1) theta variables.
    P = PolynomialRing(QQ, var_name)
    P1 = P.one(); P0 = P.zero()
    cn = c.numerator()
    var_sym = list(P.gens()) + [P0] * (n - 1)
    cd = c.denominator()(*var_sym)
    rhs = - cn(*var_sym)
    coeffs = []
    for i in range(n-1):
        var_sym[len(var_name)+i] = P1
        coeffs.append(cn(*var_sym) + rhs)
        var_sym[len(var_name)+i] = P0
    return tuple([get_pretty_fraction_polynomial(coeff / cd) for coeff in coeffs]), get_pretty_fraction_polynomial(rhs / cd)
    #return tuple([coeff / cd for coeff in coeffs]), rhs / cd

def get_pretty_fraction_polynomial(fp):
    """
    EXAMPLES::

        sage: PR2.<f,z1,z2>=QQ[]
        sage: fp = (-2*z1)/(-2*f + 4*z1 - 6)
        sage: get_pretty_fraction_polynomial(fp)
        z1/(f - 2*z1 + 3)
        sage: fp = (r0/4+z1/3+z2/3)/(1/12)
        sage: get_pretty_fraction_polynomial(fp)
        3*f + 4*z1 + 4*z2
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
    EXAMPLES::

        sage: logging.disable(logging.INFO)                                  # not tested
        sage: regions = cpl_regions_from_arrangement_of_bkpts(3, cpleq=True) # not tested
        sage: r = regions[0]                                                 # not tested
        sage: thetas = generate_thetas_of_region(r)                          # not tested
        sage: len(thetas)                                                    # not tested
        24
        sage: thetas[-1]                                                     # not tested
        ((-z)/(f - 1), (-z)/(f - 1))
    """
    constraints = symbolic_subbadditivity_constraints_of_cpl_given_region(r)
    m = len(constraints)
    n = r.parent.function._n
    var_name = r.parent.var_name
    coefficients_of_thetas = []
    constants = []
    for c in constraints:
        coeffs, rhs = coefficients_of_theta_and_rhs_in_constraint(c, var_name, n)
        coefficients_of_thetas.append(coeffs)
        constants.append(rhs)
    coeff_matrix = matrix(coefficients_of_thetas) # dim = m * (n-1)
    rows_indicies_subsets = Set(list(range(m))).subsets(n-1)
    thetas = []  # can't use set becasue of the difference between '==' and 'is'
    for rows_indicies in rows_indicies_subsets:
        indices = list(rows_indicies)
        rows = coeff_matrix.matrix_from_rows(indices)
        if rows.rank() != (n - 1):
            continue
        rhs = vector([constants[i] for i in indices])
        sol = rows * BackslashOperator() * rhs
        theta = tuple(get_pretty_fraction_polynomial(fp) for fp in sol)
        if all(theta != t for t in thetas):
            thetas.append(theta)
    return thetas

def cpl_fill_region_given_theta(r, theta, max_iter=0, flip_ineq_step=1/1000, check_completion=False, wall_crossing_method='heuristic', goto_lower_dim=True):
    """
    EXAMPLES::

        sage: logging.disable(logging.INFO)                                  # not tested
        sage: regions = cpl_regions_from_arrangement_of_bkpts(3, cpleq=True) # not tested
        sage: r = regions[0]                                                 # not tested
        sage: thetas = generate_thetas_of_region(r)                          # not tested
        sage: theta = thetas[-1]                                             # not tested
        sage: theta                                                          # not tested
        ((-z)/(f - 1), (-z)/(f - 1))
        sage: cpl_complex = cpl_fill_region_given_theta(r, theta)
        sage: len(cpl_complex.components)
        1
        sage: cpl_complex.components[0].region_type
        'is_extreme'
    """
    cpln = cpl_n_group_function(r.parent.function._n, r.parent.function._cpleq)
    cpln._theta = lambda f, z: tuple([t(f, *z) for t in theta])
    var_value = r.var_value
    var_name = r.parent.var_name
    cpl_complex=SemialgebraicComplex(cpln, var_name, max_iter=max_iter, bddleq=copy(r.leq), bddlin=copy(r.lin))
    if flip_ineq_step != 0:
        cpl_complex.bfs_completion(var_value, flip_ineq_step, check_completion, wall_crossing_method, goto_lower_dim)
    else:
        cpl_complex.shoot_random_points(1000)
    return cpl_complex

def cpl_regions_with_thetas_and_components(n=3, cpleq=True, keep_extreme_only=False,\
                                           max_iter=0, flip_ineq_step=1/1000, \
                                           check_completion=False, \
                                           wall_crossing_method='heuristic', \
                                           goto_lower_dim=True, \
                                           regions = None):
    """
    Divide the space into cells where the arrangement of breakpoints of `\\pi` is combinatorially the same.

    For each region, find theta solutions, then subdivide into smaller components by running bfs with parametric field.

    EXAMPLES::

        sage: regions = cpl_regions_with_thetas_and_components(3, True, False, 0, 1/1000, False, 'mathematica', True)   # not tested  # 25-30 mins
        sage: regions = cpl_regions_with_thetas_and_components(3, True, False, 0, 1/1000, False, 'heuristic', True)     # not tested  # 20 mins, 15 mins
    """
    if regions is None:
        regions = cpl_regions_from_arrangement_of_bkpts(n, cpleq, max_iter, flip_ineq_step, False, wall_crossing_method, goto_lower_dim) # Remark: check_completion=False for arr_complex.
    for i in range(len(regions)):
        r = regions[i]
        logging.warn("Cell %s with test point %s." %(i, r.var_value))
        r.thetas = {}
        thetas_of_r = generate_thetas_of_region(r)
        for theta in thetas_of_r:
            cpl_complex = cpl_fill_region_given_theta(r, theta, max_iter, flip_ineq_step, check_completion, wall_crossing_method, goto_lower_dim)
            if keep_extreme_only:
                components = [c for c in cpl_complex.components if c.region_type=='is_extreme']
            else:
                components = copy(cpl_complex.components)
            if components:
                r.thetas[theta] = components
    return regions


def cpl_thetas_and_regions_extreme(regions):
    """
    Gather the blue components that correspond to the same expression of theta together.

    EXAMPLES::

        sage: regions = cpl_regions_with_thetas_and_components()             # not tested
        sage: thetas_and_regions = cpl_thetas_and_regions_extreme(regions)   # not tested
        sage: len(thetas_and_regions)                                        # not tested
        9
    """
    thetas_and_regions = {}
    for r in regions:
        for (theta, components) in (r.thetas).items():
            extreme_regions = [c for c in components if c.region_type=='is_extreme']
            if not extreme_regions:
                continue
            new_theta = True
            for t in thetas_and_regions.keys():
                try:
                    if theta == tuple([ti.subs(r.var_map) for ti in t]):
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
        try:
            tt =  tuple([ti.subs(r.var_map) for ti in theta])
        except ZeroDivisionError:
            continue
        for (t, ct) in (r.thetas).items():
            if t == tt:
                components += ct
    return components

def cpl_thetas_and_regions(regions, thetas_and_regions):
    """
    Gather colorful components that correspond to the same expression of theta together.

    EXAMPLES::

        sage: regions = cpl_regions_with_thetas_and_components()             # not tested
        sage: thetas_and_regions = cpl_thetas_and_regions_extreme(regions)   # not tested
        sage: thetas_and_components = cpl_thetas_and_regions(regions, thetas_and_regions) # not tested
    """
    thetas_and_components = {}
    for theta in thetas_and_regions.keys():
        components = cpl_regions_fix_theta(regions, theta)
        thetas_and_components[theta]=components
    return thetas_and_components

def save_cpl_extreme_theta_regions(thetas_and_regions, name="cpl_theta"):
    """
    To plot only blue regions.
    Get diagrams "cpl_ext_theta_i" that show only blue regions.

    EXAMPLES::

        sage: regions = cpl_regions_with_thetas_and_components()             # not tested
        sage: thetas_and_regions = cpl_thetas_and_regions_extreme(regions)   # not tested
        sage: save_cpl_extreme_theta_regions(thetas_and_regions, name="cpl_ext_theta")  # not tested

    To plot colorful regions corresponding to each extreme thetas.
    Get diagrams "cpl_theta_i", that show a bit more::

        sage: thetas_and_components = cpl_thetas_and_regions(regions, thetas_and_regions) # not tested
        sage: save_cpl_extreme_theta_regions(thetas_and_components, name="cpl_theta") # not tested 
        (1, (z/(15*z - 2), (6*z - 1)/(15*z - 2)))
        (2, ((-2*z)/(f - 1), 0))
        (3, (1/4, 1/4))
        (4, ((f + z)/(f + 1), z/(f + 1)))
        (5, ((-z)/(f + 2*z - 1), 0))
        (6, (2*z/(12*z - 1), (z - 1/6)/(2*z - 1/6)))
        (7, (1/2, 0))
        (8, ((f + 2*z)/(2*f + 2), (f + 2*z)/(2*f + 2)))
        (9, ((-z)/(f - 1), (-z)/(f - 1)))

    If wall_crossing_mehtod='mathematica' in ``cpl_regions_with_thetas_and_components()``, then the output is::

        (1, (z/(15*z - 2), (6*z - 1)/(15*z - 2)))
        (2, ((-2*z)/(f - 1), 0))
        (3, (1/4, 1/4))
        (4, ((f + z)/(f + 1), z/(f + 1)))
        (5, ((-z)/(f + 2*z - 1), 0))
        (6, (2*z/(12*z - 1), (z - 1/6)/(2*z - 1/6)))
        (7, ((f + 2*z)/(2*f + 2), (f + 2*z)/(2*f + 2)))
        (8, (1/2, 0))
        (9, ((-z)/(f - 1), (-z)/(f - 1)))
    """
    k = 0
    for (theta, components) in thetas_and_regions.items():
        k += 1
        print((k, theta))
        if len(str(theta)) < 100:
            title = "extreme point %s:  theta = %s" % (k, theta)
        else:
            n_theta = len(theta)
            title = "extreme point %s:\n" % k
            for i in range(n_theta):
                title += "theta_%s = %s \n" % (i+1, theta[i])
        g = plot_cpl_components(components)
        g += text(title, (0.5, 1/4), color='black')
        g.save(name+"_%s.png" %k, xmin=0, xmax=1, ymin=0, ymax=1/4)

def collect_and_save_cpl_extreme_theta_regions(regions, name="cpl_theta"):
    thetas_and_regions = cpl_thetas_and_regions_extreme(regions)
    thetas_and_components = cpl_thetas_and_regions(regions, thetas_and_regions)
    save_cpl_extreme_theta_regions(thetas_and_components, name=name)
