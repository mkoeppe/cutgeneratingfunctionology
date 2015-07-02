# attach kslope_ppl_mip.py

vertex_enumeration_dim_threshold = 10

def pattern_setup_lp(l, more_ini_additive=False, objcoef=None):
    r"""
    Set up a MixedIntegerLinearProgram() with respect to the prescribled painting (pattern=0, sym=True).

    Returns:
        - fn: function values at 1/q\Z in terms of s_0, s_1, ... , s_{l+1}
    Global variables:
        - pattern_lp: MixedIntegerLinearProgram()
        - var_slope: slope variables s_0, s_1, ... , s_{l+1}
        - var_delta: auxiliary variables deltafn;
                     controle additive/subadd by pattern_lp.set_max(var_delta[deltafn], 0) or pattern_lp.set_max(var_delta[deltafn], None)      
        - deltafn_dic: a dictionary, key = delta, value = a list [(i,j) such that \delta\pi[i,j]=key]

    EXAMPLES::

        sage: fn = pattern_setup_lp(1, more_ini_additive=False, objcoef=None)
        sage: fn[28]
        2*x_0 + 12*x_1 + 14*x_2
        sage: pattern_lp
        sage: var_slope
        MIPVariable of dimension 1.
        sage: var_delta
        MIPVariable of dimension 1.
        Mixed Integer Program  ( minimization, 95 variables, 121 constraints )
        sage: len(deltafn_dic.items())
        92
        sage: deltafn_dic.items()[0]
        ((2, 10, 13), [(25, 29), (25, 33)])
    """
    global pattern_lp, var_slope, var_delta, deltafn_dic
    q, f = pattern_q_and_f(l, 0)
    pattern_lp = MixedIntegerLinearProgram(maximization=True, solver = "GLPK")
    var_slope = pattern_lp.new_variable(real=True, nonnegative=False)
    s = [var_slope[0], var_slope[1]]
    for k in range(1, l + 1):
        s += [var_slope[k], var_slope[k + 1]] * 3
    s += [var_slope[l + 1], var_slope[l + 1], var_slope[l + 1]]
    for k in range(l, 0, -1):
        s += [var_slope[k], var_slope[k + 1], var_slope[k]]
    if objcoef is None:
        objcoef = [16^i for i in range(l+2)]
    objfun = sum([objcoef[i] * var_slope[i] for i in range(l+2)])
    s = s + [s[0]] + s[-1::-1]
    fn = [0] * (q + 1)
    for k in range(f):
        fn[k + 1] = fn[k] + s[k]
    for k in range(f):
        fn[q - k] = fn[k]
    vertices_color = pattern_vertices_color(l, 0, more_ini_additive=more_ini_additive)
    pattern_lp.add_constraint(fn[f] == 1)
    for i in range(1, f):
        pattern_lp.add_constraint(fn[i], max=1, min=0)
    var_delta = pattern_lp.new_variable(real=True, nonnegative=True)
    deltafn_dic = {}
    for x in range(1, f + 1):
        for y in range(x, q - x + 1):
            # cannot use mutable key: deltafn = fn[x] + fn[y] - fn[x + y]
            deltafn = tuple_of_deltafn(l+2, fn[x] + fn[y] - fn[x + y])
            if deltafn in deltafn_dic.keys():
                deltafn_dic[deltafn].append((x,y))
            else:
                deltafn_dic[deltafn] = [(x,y)]
                pattern_lp.add_constraint(fn[x] + fn[y] - fn[x + y] - var_delta[deltafn], min=0, max=0)
    for deltafn, xypairs in deltafn_dic.items():
        is_additive = False
        for x, y in xypairs:
            if vertices_color[x, y] == 0:
                is_additive = True
                break
        if is_additive:
            pattern_lp.set_max(var_delta[deltafn], 0)
        else:
            pattern_lp.set_max(var_delta[deltafn], None)
    pattern_lp.set_objective(objfun)
    pattern_lp.solver_parameter("primal_v_dual", "GLP_PRIMAL")
    return fn

def pattern_positive_zero_undecided_deltafn(vertices_color):
    """
    According to the prescribled painting, compute the lists positive_deltafn, zero_deltafn and undecided_deltafn.
    zero_deltafn does not include the implied zero deltafn.
    pattern_lp, vertices_color might be modified.

    EXAMPLES::

        sage: l = 1; q, f = pattern_q_and_f(l, 0);
        sage: vertices_color = pattern_vertices_color(l, pattern=0, more_ini_additive=False)
        sage: fn = pattern_setup_lp(l)
        sage: positive_deltafn, zero_deltafn, undecided_deltafn = pattern_positive_zero_undecided_deltafn(vertices_color)
        ...
        sage: positive_deltafn
        [(1, -1, 0), (1, 0, 1), (1, 0, -1), (1, 1, 0), (0, 1, -1), (1, 0, 0)]
        sage: zero_deltafn
        [(0, 0, 0)]
        sage: undecided_deltafn
        [(1, 7, 12), (1, 2, 5), (1, 2, -3), (1, 1, -2), (1, 8, 13), (1, 11, 14)]
    """
    positive_deltafn = []
    zero_deltafn = []
    undecided_deltafn = []
    global pattern_lp, var_delta, deltafn_dic

    # look for implied additivities.
    #pattern_lp.solver_parameter("obj_upper_limit", 0.01)
    pattern_lp.solver_parameter("primal_v_dual", "GLP_PRIMAL")
    pattern_lp.get_backend().set_sense(1) # maximization
    for deltafn, xypairs in deltafn_dic.items():
        v = var_delta[deltafn]
        if pattern_lp.get_max(v) == 0:
            zero_deltafn.append(deltafn)
        else:
            pattern_lp.set_objective(v)
            if glpk_simplex_exact_solve(pattern_lp) == 0: # maximization
                pattern_lp.set_max(v, 0)
        if pattern_lp.get_max(v) == 0:
            for x, y in xypairs:
                vertices_color[x, y] = 0
    #pattern_lp.solver_parameter("obj_upper_limit", None)
    pattern_lp.get_backend().set_sense(-1) # minimization
    for deltafn, xypairs in deltafn_dic.items():
        v = var_delta[deltafn]
        # For deltafn that is not already 0,
        # if (x,y) in xypairs forms an edge with initial additive points, then deltafn belongs to positive_deltafn,
        # otherwise deltafn belongs to undecided
        if pattern_lp.get_max(v)!= 0:
            is_undecided = True
            for x, y in xypairs:
                if pattern_will_form_edge(x, y, vertices_color):
                    is_undecided = False
                    positive_deltafn.append(deltafn)
                    break
            if is_undecided:
                pattern_lp.set_objective(v)
                if glpk_simplex_exact_solve(pattern_lp) == 0: # minimization
                    undecided_deltafn.append(deltafn)
    # FIXME: sort undecided_deltafn randomly ??
    return positive_deltafn, zero_deltafn, undecided_deltafn

def pattern_will_form_edge(x, y, vertices_color):
    r"""
    
    """
    # assume 1<= x <=f, x <= y <= q-x
    for (i, j) in [(-1,0), (-1, 1), (0, 1), (1, 0), (1, 1), (0, -1)]:
        if vertices_color[x-1, y] == 0:
            return True
    return False

def tuple_of_deltafn(n, linfun):
    t = [0]*n
    for i,j in linfun.dict().items():
        if i != -1:
            t[i] = int(j)
    d = gcd(t)
    if d == 0:
        return tuple(t)
    return tuple([x/d for x in t])

def pattern_backtrack_polytope(l, k_slopes):
    """
    sage: igp.vertex_enumeration_dim_threshold = 1
    sage: pattern_backtrack_polytope(1, 6)
    sage: igp.vertex_enumeration_dim_threshold = 3
    sage: pattern_backtrack_polytope(4, 10)
    """
    global deltafn_dic

    q, f = pattern_q_and_f(l, 0)
    vertices_color = pattern_vertices_color(l, pattern=0, more_ini_additive=False)
    fn = pattern_setup_lp(l)
    positive_deltafn, zero_deltafn, undecided_deltafn = pattern_positive_zero_undecided_deltafn(vertices_color)
    exp_dim_ = l+1
    gen = pattern_branch_on_deltafn(positive_deltafn, zero_deltafn, undecided_deltafn, vertices_color, exp_dim_)

    # Create the ppl polytope (variables = var_slope[0]...var_slope[l+1]).
    # First put in the trivial constraints for minimal.
    cs_trivial = Constraint_System()
    cs_trivial.insert(convert_linfun_to_linexp(fn[f]) == 1)
    for i in range(1, f):
        cs_trivial.insert(convert_linfun_to_linexp(fn[i]) >= 0)
        cs_trivial.insert(convert_linfun_to_linexp(fn[i]) <= 1)
    for deltafn, xypairs in deltafn_dic.items():
        x, y = xypairs[0]
        cs_trivial.insert(convert_linfun_to_linexp(fn[x] + fn[y] - fn[x + y]) >= 0)
    # Then Add more constraints corresponding to deltafn = 0.
    for zero_delta_list, exp_dim in gen:
        #print zero_delta_list
        polytope = C_Polyhedron(cs_trivial)
        for zero_delta in zero_delta_list:
            x, y = deltafn_dic[zero_delta][0]
            polytope.add_constraint(convert_linfun_to_linexp(fn[x] + fn[y] - fn[x + y]) == 0)
        print polytope
        print "exp_dim = %s" % exp_dim
        #Call pattern_extreme() to enumerate the vertices and record >=k slope functions.
        max_num_slopes = pattern_extreme(l, k_slopes, pattern=0, show_plots=False,
                    test_extremality=False, polytope = polytope,
                    more_ini_additive=False, count_components=False, use_sha1=True)
        #If >=k slope vertex function is found, stop backtracking.
        if max_num_slopes > 0:
            return max_num_slopes

def convert_linfun_to_linexp(linfun):
    """
    convert MILP's Linear_Function to PPL's Linear_Eexpression
    """
    return sum([ Variable(i)*j for i,j in linfun.dict().items() if i != -1])

def pattern_branch_on_deltafn(positive_deltafn, zero_deltafn, undecided_deltafn, vertices_color, exp_dim):
    global pattern_lp, var_delta, deltafn_dic

    if exp_dim <= vertex_enumeration_dim_threshold:
        yield zero_deltafn, exp_dim
    elif undecided_deltafn:
        branch_delta = undecided_deltafn[0]
        # try setting branch_delta = 0
        # check: no green edge caused by branch_delta = 0
        is_feasible = True
        for x, y in deltafn_dic[branch_delta]:
            vertices_color[x, y] = 0
            if pattern_will_form_edge(x, y, vertices_color):
                is_feasible = False
                break
        # check: pattern_lp is feasible
        if is_feasible:
            pattern_lp.set_max(var_delta[branch_delta], 0)
            pattern_lp.solver_parameter("simplex_or_intopt", "simplex_only")
            pattern_lp.solver_parameter("primal_v_dual", "GLP_DUAL")
            try:
                pattern_lp.solve(objective_only=True) # maximization or minimation doesn't matter
            except MIPSolverException:
                # If infeasible, stop recursion
                is_feasible = False
        # check: none of the positive_deltafn becomes 0
        if is_feasible:
            pattern_lp.solver_parameter("primal_v_dual", "GLP_PRIMAL")
            #pattern_lp.solver_parameter("obj_upper_limit", 0.01)
            pattern_lp.get_backend().set_sense(1) # maximization
            for pos_delta in positive_deltafn:
                pattern_lp.set_objective(var_delta[pos_delta])
                if glpk_simplex_exact_solve(pattern_lp) == 0: # maximization
                    is_feasible = False
                    break
        # NOT NEEDED  check: none of the zero_deltafn becomes >0, since otherwise problem is infeasible
        #if is_feasible:
        #    pattern_lp.solver_parameter("obj_upper_limit", None)
        #    pattern_lp.get_backend().set_sense(-1) # minimization
        #    for zer_delta in zero_deltafn:
        #        pattern_lp.set_objective(var_delta[zer_delta])
        #        if glpk_simplex_exact_solve(pattern_lp) != 0: # minimization
        #            is_feasible = False
        #            break

        # check: no green edge caused by implied_zero_deltafn
        implied_zero_deltafn = []
        still_undecided_deltafn = []
        if is_feasible:
            for deltafn in undecided_deltafn[1::]:
                pattern_lp.set_objective(var_delta[deltafn])
                pattern_lp.get_backend().set_sense(1) # maximization
                if glpk_simplex_exact_solve(pattern_lp) == 0: #maximation
                    implied_zero_deltafn.append(deltafn)
                else:
                    pattern_lp.get_backend().set_sense(-1) # minimization
                    if glpk_simplex_exact_solve(pattern_lp) == 0: # minimization
                        still_undecided_deltafn.append(deltafn)
            for deltafn in implied_zero_deltafn:
                for x, y in deltafn_dic[deltafn]:
                    vertices_color[x, y] = 0
                    pattern_lp.set_max(var_delta[deltafn], 0)
                    if pattern_will_form_edge(x, y, vertices_color):
                        is_feasible = False
                        break
                if not is_feasible:
                    break
        if is_feasible:
            for result in pattern_branch_on_deltafn(positive_deltafn, zero_deltafn+[branch_delta], \
                                            still_undecided_deltafn, vertices_color, exp_dim-1):
                yield result
        for deltafn in [branch_delta]+implied_zero_deltafn:
            pattern_lp.set_max(var_delta[deltafn], None)
            for x, y in deltafn_dic[deltafn]:
                vertices_color[x, y] = 1
        for result in  pattern_branch_on_deltafn(positive_deltafn+[branch_delta], zero_deltafn, \
                                            undecided_deltafn[1::], vertices_color, exp_dim):
            yield result

def glpk_simplex_exact_solve(lp):
    lp.solver_parameter("simplex_or_intopt", "simplex_only")
    try:
        optval = lp.solve()
    except MIPSolverException:
        return 'Infeasible'
    lp.solver_parameter("simplex_or_intopt", "exact_simplex_only")
    optval = lp.solve()
    return optval

def pattern_glpk_lp(l, more_ini_additive=False, exact_arithmetic=True, simplex_first=True, reconstruct_rational=False, objcoef=None):
    # pattern=0, guess obj function, solve lp
    fn = pattern_setup_lp(l, more_ini_additive=more_ini_additive, objcoef=objcoef)
    if exact_arithmetic and not simplex_first:
        lp.solver_parameter("simplex_or_intopt", "exact_simplex_only")
    else:
        lp.solver_parameter("simplex_or_intopt", "simplex_only")
    lp.solver_parameter("primal_v_dual", "GLP_PRIMAL")
    try:
        optval = lp.solve()
    except MIPSolverException:
        return 'NA', 'NA', 'NA', 'NA'
    if exact_arithmetic and simplex_first:
        lp.solver_parameter("simplex_or_intopt", "exact_simplex_only")
        optval = lp.solve()
    if reconstruct_rational:
        b = lp.get_backend() 
        optsol = exact_optsol(b)
        k_slope = len(set(optsol) | set([-i for i in optsol]))
        optval = 0
        for i in range(b.ncols()):
            optval += QQ(b.objective_coefficient(i)) * optsol[i]
    else:
        optsol = lp.get_values(var_slope)
        slopes = optsol.values()
        k_slope = len(set(slopes) | set([-i for i in slopes]))
    v = [0]
    for i in range(1, q):
        d = fn[i].dict()
        d.pop(-1, None)
        value = 0
        for j in d:
            value += optsol[j] * QQ(d[j])
        v.append(value)
    v.append(0)
    #bkpt = [i/q for i in range(q+1)]
    #h = piecewise_function_from_breakpoints_and_values(bkpt, v)
    #return optval, optsol, k_slope, h
    return optval, optsol, k_slope, v

def exact_optsol(b):
    #sage_input(b)
    ncol = b.ncols()
    nrow = b.nrows()
    A = matrix(QQ, ncol + nrow, ncol + nrow, sparse = True)
    for i in range(nrow):
        r = b.row(i)
        for (j, c) in itertools.izip(r[0], r[1]):
            A[i, j] = QQ(c)
        A[i, ncol + i] = -1
    n = nrow
    Y = zero_vector(QQ, ncol + nrow)
    for i in range(ncol):
        status =  b.get_col_stat(i)
        if status > 1:
            A[n, i] = 1
            if status == 2:
                Y[n] = b.col_bounds(i)[0]
            else:
                Y[n] = b.col_bounds(i)[1]
            n += 1
    for i in range(nrow):
        status =  b.get_row_stat(i)
        if status > 1:
            A[n, ncol + i] = 1
            if status == 2:
                Y[n] = b.row_bounds(i)[0]
            else:
                Y[n] = b.row_bounds(i)[1]
            n += 1

    #filename = open(dir_math+"profiler/solveAXisY", "w")
    #print >> filename, "A =",
    #print >> filename, sage_input(A)
    #print >> filename
    #print >> filename, "Y =",
    #print >> filename, sage_input(Y)
    #filename.close()

    #X = A \ Y
    X = A.solve_right(Y, check=False)
    return X[0:ncol]    
            

def pattern_ppl_lp(l, more_ini_additive=False, objcoef=None): #TOO SLOW
    pattern = 0
    q, f = pattern_q_and_f(l, pattern)
    vertices_color = pattern_vertices_color(l, pattern, more_ini_additive=more_ini_additive)
    s = pattern_s(l, pattern)
    fn = pattern_fn(l, pattern)
    #objfun = Linear_Expression(0)
    #objfun = sum([16^i * s[i] for i in range(len(s))])
    if objcoef is None:
        objcoef = [16^i for i in range(l+2)]
    objfun = sum([objcoef[i] * var_slope[i] for i in range(l+2)])
    cs = Constraint_System()
    cs.insert(fn[0] == 0)
    cs.insert(fn[f] == 1)
    for i in range(1, f):
        cs.insert(fn[i] >= 0)
        cs.insert(fn[i] <= 1)
    for x in range(1, f + 1):
        for y in range(x, q - x + 1):
            if vertices_color[x, y] == 0:
                cs.insert(fn[x] + fn[y] == fn[x + y])
            else:
                cs.insert(fn[x] + fn[y] >= fn[x + y])
    lp = MIP_Problem(l + 2, cs, objfun)
    optval = lp.optimal_value()
    optsol = lp.optimizing_point()
    s_coe = optsol.coefficients()
    k_slope = len(set(s_coe) | set([-i for i in s_coe]))
    v_div = optsol.divisor()   
    v = [sum(vf * vs for vf, vs in zip(fn[i].coefficients(), s_coe))/v_div for i in range(q+1)]
    #bkpt = [i/q for i in range(q+1)]
    #h = piecewise_function_from_breakpoints_and_values(bkpt, v)
    #return optval, optsol, k_slope, h, v_div
    return optval, optsol, k_slope, v, v_div

def pattern_glpk_test(l_list, more_ini_additive=False, exact_arithmetic=True, simplex_first=True, reconstruct_rational=False):
    slopes = []
    for l in l_list:
        start_cpu_t = time.clock();
        optval, optsol, k, v_glpk = pattern_glpk_lp(l, more_ini_additive=more_ini_additive, \
                                                    exact_arithmetic=exact_arithmetic, simplex_first = simplex_first, \
                                                    reconstruct_rational=reconstruct_rational);
        cpu_t = time.clock();
        print l, k, cpu_t - start_cpu_t
        slopes.append(k)
    return slopes

def pattern_ppl_test(l_list, more_ini_additive=False):
    slopes = []
    for l in l_list:
        start_cpu_t = time.clock();
        try:
            optval, optsol, k, v_ppl, v_div = pattern_ppl_lp(l, more_ini_additive=more_ini_additive);
            cpu_t = time.clock();
            print l, k, cpu_t - start_cpu_t
            slopes.append(k)
        except ValueError:
            cpu_t = time.clock();
            print l, "NA", cpu_t - start_cpu_t
            slopes.append(-1)
    return slopes
    
def piecewise_function_from_values_and_evenlyspaced_bkpts(v):
    q = len(v) - 1
    bkpt = [i/q for i in range(q + 1)]
    h = piecewise_function_from_breakpoints_and_values(bkpt, v)
    return h

    
