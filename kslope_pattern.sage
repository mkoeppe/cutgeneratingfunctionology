# attach kslope_ppl_mip.py

def pattern_glpk_lp(l, more_ini_additive=True, exact_arithmetic=True, reconstruct_rational=True):
    # pattern=0, guess obj function, solve lp
    q, f = pattern_q_and_f(l, 0)
    lp = MixedIntegerLinearProgram(maximization=True, solver = "GLPK")
    var_slope = lp.new_variable(real=True, nonnegative=False)
    s = [var_slope[0], var_slope[1]]
    for k in range(1, l + 1):
        s += [var_slope[k], var_slope[k + 1]] * 3
    s += [var_slope[l + 1], var_slope[l + 1], var_slope[l + 1]]
    for k in range(l, 0, -1):
        s += [var_slope[k], var_slope[k + 1], var_slope[k]]
    objfun = sum([16^i * s[i] for i in range(len(s))])
    s = s + [s[0]] + s[-1::-1]
    fn = [0] * (q + 1)
    for k in range(f):
        fn[k + 1] = fn[k] + s[k]
    for k in range(f):
        fn[q - k] = fn[k]
    vertices_color = pattern_vertices_color(l, 0, more_ini_additive=more_ini_additive)
    lp.add_constraint(fn[f] == 1)
    for i in range(1, f):
        lp.add_constraint(fn[i], max=1, min=0)
    for x in range(1, f + 1):
        for y in range(x, q - x + 1):
            if vertices_color[x, y] == 0:
                lp.add_constraint(fn[x] + fn[y] == fn[x + y])
            else:
                lp.add_constraint(fn[x] + fn[y] >= fn[x + y]) 
    lp.set_objective(objfun)
    #lp.show()
    if exact_arithmetic:
        lp.solver_parameter("simplex_or_intopt", "exact_simplex_only")
    else:
        lp.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
    lp.solver_parameter("primal_v_dual", "GLP_PRIMAL")  ### ???
    try:
        optval = lp.solve()
    except MIPSolverException:
        return 'NA', 'NA', 'NA', 'NA'
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
            

def pattern_ppl_lp(l, more_ini_additive=True): #TOO SLOW
    pattern = 0
    q, f = pattern_q_and_f(l, pattern)
    vertices_color = pattern_vertices_color(l, pattern, more_ini_additive=more_ini_additive)
    s = pattern_s(l, pattern)
    fn = pattern_fn(l, pattern)
    #objfun = Linear_Expression(0)
    objfun = sum([16^i * s[i] for i in range(len(s))])
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

def pattern_glpk_test(l_list, more_ini_additive=True, exact_arithmetic=True, reconstruct_rational=True):
    slopes = []
    for l in l_list:
        start_cpu_t = time.clock();
        optval, optsol, k, v_glpk = pattern_glpk_lp(l, more_ini_additive=more_ini_additive, \
                                                    exact_arithmetic=exact_arithmetic, \
                                                    reconstruct_rational=reconstruct_rational);
        cpu_t = time.clock();
        print l, k, cpu_t - start_cpu_t
        slopes.append(k)
    return slopes

def pattern_ppl_test(l_list, more_ini_additive=True):
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

    
