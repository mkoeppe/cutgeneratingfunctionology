from six.moves import range
from six.moves import zip
vertex_enumeration_dim_threshold = 10

import sys
sys.setrecursionlimit(50000)

try:
    from ppl import Variable, Constraint, Linear_Expression, Constraint_System, NNC_Polyhedron, Poly_Con_Relation, Poly_Gen_Relation, Generator, MIP_Problem, point as ppl_point
except ImportError:
    # old Sage
    from sage.libs.ppl import Variable, Constraint, Linear_Expression, Constraint_System, NNC_Polyhedron, Poly_Con_Relation, Poly_Gen_Relation, Generator, MIP_Problem, point as ppl_point

def pattern_setup_lp(l, more_ini_additive=False, objcoef=None, use_auxiliary_delta=True):
    r"""
    Set up a ``MixedIntegerLinearProgram()`` with respect to the prescribled painting,
    where pattern=0, sym=True and the size parameter is `l`.

    Returns:
        - fn: a list of LinearFunction in terms of `s_0, s_1, \dots , s_{l+1}`, corresponding to the function values at `(1/q)Z`.
    Global variables:
        - pattern_lp: ``MixedIntegerLinearProgram()``
        - var_slope: slope variables`s_0, s_1, \dots , s_{l+1}`

    When use_auxiliary_delta=True:
        - var_delta: auxiliary variables deltafn;
                     controle additive/subadd by ``pattern_lp.set_max(var_delta[deltafn], 0)``
                     or ``pattern_lp.set_max(var_delta[deltafn], None)``
        - deltafn_dic: a dictionary, key = delta, value = a list [(i,j) such that `\delta\pi`[i,j]=key]

    EXAMPLES::

        sage: import cutgeneratingfunctionology.igp as igp
        sage: fn = igp.pattern_setup_lp(1, more_ini_additive=False, objcoef=None)
        sage: fn[28]
        2*x_0 + 12*x_1 + 14*x_2
        sage: igp.pattern_lp.number_of_variables()
        95
        sage: igp.var_slope
        MIPVariable ...
        sage: igp.var_delta
        MIPVariable ...
        sage: len(igp.deltafn_dic.items())
        92
        sage: igp.deltafn_dic[(2, 10, 13)]
        [(25, 29), (25, 33)]
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
    if use_auxiliary_delta:
        var_delta = pattern_lp.new_variable(real=True, nonnegative=True)
        deltafn_dic = {}
        for x in range(1, f + 1):
            for y in range(x, q - x + 1):
                # cannot use mutable key: deltafn = fn[x] + fn[y] - fn[x + y]
                deltafn = tuple_of_deltafn(l+2, fn[x] + fn[y] - fn[x + y])
                if deltafn in list(deltafn_dic.keys()):
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
    else:
        for x in range(1, f + 1):
            for y in range(x, q - x + 1):
                if vertices_color[x, y] == 0:
                    pattern_lp.add_constraint(fn[x] + fn[y] == fn[x + y])
                else:
                    pattern_lp.add_constraint(fn[x] + fn[y] >= fn[x + y])
    pattern_lp.set_objective(objfun)
    pattern_lp.solver_parameter("primal_v_dual", "GLP_PRIMAL")
    return fn

def pattern_positive_zero_undecided_deltafn(vertices_color):
    r"""
    According to the prescribed painting, compute the lists positive_deltafn, zero_deltafn and undecided_deltafn.

    zero_deltafn does not include the implied zero deltafn.
    pattern_lp, vertices_color might be modified.

    See also ``pattern_positive_zero_undecided_deltafn_trivial()``.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: l = 1; q, f = pattern_q_and_f(l, 0);
        sage: vertices_color = pattern_vertices_color(l, pattern=0, more_ini_additive=False)
        sage: print("glp_exact noise follows in old sage versions"); fn = pattern_setup_lp(l)
        glp_exact...
        sage: print("glp_exact noise follows in old sage versions"); positive_deltafn, zero_deltafn, undecided_deltafn = pattern_positive_zero_undecided_deltafn(vertices_color)
        glp_exact...
        sage: sorted(positive_deltafn)
        [(0, 1, -1), (1, -1, 0), (1, 0, -1), (1, 0, 0), (1, 0, 1), (1, 1, 0)]
        sage: zero_deltafn
        [(0, 0, 0)]
        sage: sorted(undecided_deltafn)
        [(1, 1, -2), (1, 2, -3), (1, 2, 5), (1, 7, 12), (1, 8, 13), (1, 11, 14)]
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

def pattern_positive_zero_undecided_deltafn_trivial(vertices_color):
    r"""
    Compute the list zero_deltafn of deltafn that corresponds to green vertices on the prescribled painting.

    The other deltafn are put in the list undecided_deltafn.
    positive_deltafn is an empty list.
    Modifiy vertices_color accordingly.

    Unlike ``pattern_positive_zero_undecided_deltafn()``,
    this function does not call ``glpk_simplex_exact_solve()`` and ``pattern_will_form_edge()``
    to compute implied green and explicit white vertices.
    Therefore, when ``vertex_enumeration_dim_threshold`` is large enough for given `q` so that
    no backtracking is needed in ``pattern_backtrack_polytope()``, it will not waste time on ``exact_solve``.
    In this case, ``pattern_backtrack_polytope()`` should perform as good as ``pattern_glpk_lp()``.
    (Hopefully even better, since same deltafn value is accounted only once in constraints.)

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: l = 1; q, f = pattern_q_and_f(l, 0);
        sage: vertices_color = pattern_vertices_color(l, pattern=0, more_ini_additive=False)
        sage: fn = pattern_setup_lp(l)
        sage: positive_deltafn, zero_deltafn, undecided_deltafn = pattern_positive_zero_undecided_deltafn_trivial(vertices_color)
        sage: positive_deltafn
        []
        sage: zero_deltafn
        [(0, 0, 0)]
        sage: len(undecided_deltafn)
        91
    """
    positive_deltafn = []
    zero_deltafn = []
    undecided_deltafn = []
    global pattern_lp, var_delta, deltafn_dic

    for deltafn, xypairs in deltafn_dic.items():
        v = var_delta[deltafn]
        if pattern_lp.get_max(v) == 0:
            zero_deltafn.append(deltafn)
            for x, y in xypairs:
                vertices_color[x, y] = 0
        else:
            undecided_deltafn.append(deltafn)
    # FIXME: sort undecided_deltafn randomly ??
    return positive_deltafn, zero_deltafn, undecided_deltafn

def pattern_will_form_edge(x, y, vertices_color):
    r"""
    Return whether paiting (x, y) green on vertices_color will create a green edge. Assume x <= y.
    
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: l = 1; q, f = pattern_q_and_f(l, 0);
        sage: vertices_color = pattern_vertices_color(l, pattern=0, more_ini_additive=False)
        sage: pattern_will_form_edge(1, 1, vertices_color)
        True
        sage: pattern_will_form_edge(2, q, vertices_color)
        True
        sage: pattern_will_form_edge(2, q-1, vertices_color)
        False
    """
    # assume 1<= x <=f, x <= y <= q-x
    for (i, j) in [(-1,0), (-1, 1), (0, 1), (1, 0), (1, 1), (0, -1)]:
        if vertices_color[x-1, y] == 0:
            return True
    return False

def tuple_of_deltafn(n, linfun):
    r"""
    Convert a Linear Function (with variable index < n) to a n-tuple,
    dividing by the gcd of all elements.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: fn = pattern_setup_lp(1, more_ini_additive=False, objcoef=None)
        sage: linfun = fn[5]; linfun
        x_0 + 3*x_1 + x_2
        sage: tuple_of_deltafn(3, linfun)
        (1, 3, 1)
        sage: linfun = fn[18]+fn[19]-fn[37]; linfun
        2*x_0 + 8*x_1 + 6*x_2
        sage: tuple_of_deltafn(3, linfun)
        (1, 4, 3)
    """
    t = [0]*n
    for i,j in linfun.dict().items():
        if i != -1:
            t[i] = int(j)
    d = gcd(t)
    if d == 0:
        return tuple(t)
    return tuple([x/d for x in t])

def pattern_backtrack_polytope(l, k_slopes):
    r"""
    (MAIN FUNCTION) Look for >= k_slopes vertex-functoins on prescribed painting
    according to pattern=0 and size parameter `l`.

    Backtrack to find some green points such that they do not create new green edge.
    Start the vertex-enumeration when exp_dim of the polytope is at most ``vertex_enumeration_dim_threshold``.

    Output the >= k_slopes vertex-functions (write to file), and max_num_slopes.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: import cutgeneratingfunctionology.igp as igp
        sage: vertex_enumeration_dim_threshold = 1
        sage: pattern_backtrack_polytope(1, 6) # long time
        #####  ...
        6
        sage: pattern_backtrack_polytope(1, 8) # long time
    """
    global deltafn_dic

    q, f = pattern_q_and_f(l, 0)
    vertices_color = pattern_vertices_color(l, pattern=0, more_ini_additive=False)
    fn = pattern_setup_lp(l)
    positive_deltafn, zero_deltafn, undecided_deltafn = pattern_positive_zero_undecided_deltafn_trivial(vertices_color)
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
        #print polytope
        #print "exp_dim = %s" % exp_dim
        #Call pattern_extreme() to enumerate the vertices and record >=k slope functions.
        max_num_slopes = pattern_extreme(l, k_slopes, pattern=0, show_plots=False,
                    test_extremality=False, polytope = polytope, exp_dim = exp_dim,
                    more_ini_additive=False, count_components=False, use_sha1=True)
        #If >=k slope vertex function is found, stop backtracking.
        if max_num_slopes > 0:
            return max_num_slopes

def convert_linfun_to_linexp(linfun):
    r"""
    convert MILP's ``Linear_Function`` to PPL's ``Linear_Expression``.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: fn = pattern_setup_lp(1, more_ini_additive=False, objcoef=None)
        sage: linfun = fn[5]; linfun
        x_0 + 3*x_1 + x_2
        sage: type(linfun)
        <type 'sage.numerical.linear_functions.LinearFunction'>
        sage: linexp = convert_linfun_to_linexp(linfun); linexp
        x0+3*x1+x2
        sage: type(linexp)
        <type '...Linear_Expression'>
    """
    return sum([ Variable(i)*j for i,j in linfun.dict().items() if i != -1])

def pattern_branch_on_deltafn(positive_deltafn, zero_deltafn, undecided_deltafn, vertices_color, exp_dim):
    r"""
    Backtracking subroutine of ``pattern_backtrack_polytope()``.
    """
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
                pattern_lp.set_max(var_delta[deltafn], 0)
                for x, y in deltafn_dic[deltafn]:
                    vertices_color[x, y] = 0
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
    r"""
    Solve lp by ``glp_simplex`` + ``glp_exact``

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: lp = MixedIntegerLinearProgram(solver = 'GLPK', maximization = False)
        sage: x, y = lp[0], lp[1]
        sage: lp.add_constraint(-2*x + y <= 1)
        sage: lp.add_constraint(x - y <= 1)
        sage: lp.add_constraint(x + y >= 2)
        sage: lp.set_objective(x + y)
        sage: print("glp_exact noise follows in old sage versions"); glpk_simplex_exact_solve(lp)
        glp_exact...
        2.0
    """
    lp.solver_parameter("simplex_or_intopt", "simplex_only")
    try:
        optval = lp.solve()
    except MIPSolverException:
        return 'Infeasible'
    lp.solver_parameter("simplex_or_intopt", "exact_simplex_only")
    optval = lp.solve()
    return optval

def pattern_glpk_lp(l, more_ini_additive=False, exact_arithmetic=True, simplex_first=True, reconstruct_rational=False, objcoef=None):
    r"""
    Find an extreme point of the polytope corresponding to the prescribed painting (pattern=0 and size parameter `l`),
    by solving the LP (given objective funtion's coefficient) using ``glp_simplex`` (+ ``glp_exact``, with or without retrieving rational solution).

    We hope to obtain many-slope vertex functions by choosing objcoef nicely.

    Returns:
        - optval: optimal value
        - optsol: optimal solution `s0, s1, \dots , s_{l+1}`
        - k_slope: number of slopes of this vertex function
        - v: vertex function's values at (1/q)Z.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: print("glp_exact noise follows in old sage versions"); optval, optsol, k_slope, v = pattern_glpk_lp(1, objcoef=(12,55,0))
        glp_exact...
        sage: k_slope
        6
        sage: optsol
        {0: 0.23404255319148934, 1: 0.14893617021276595, 2: -0.10638297872340424}

        sage: print("glp_exact noise follows in old sage versions"); optval, optsol, k_slope, v = pattern_glpk_lp(1, reconstruct_rational=True,objcoef=(12,55,0))
        glp_exact...
        sage: optsol
        (11/47, 7/47, -5/47)
    """
    # pattern=0, guess obj function, solve lp
    fn = pattern_setup_lp(l, more_ini_additive=more_ini_additive, objcoef=objcoef, use_auxiliary_delta=False)
    q = len(fn) - 1
    if exact_arithmetic and not simplex_first:
        pattern_lp.solver_parameter("simplex_or_intopt", "exact_simplex_only")
    else:
        pattern_lp.solver_parameter("simplex_or_intopt", "simplex_only")
    pattern_lp.solver_parameter("primal_v_dual", "GLP_PRIMAL")
    try:
        optval = pattern_lp.solve()
    except MIPSolverException:
        return 'NA', 'NA', 'NA', 'NA'
    if exact_arithmetic and simplex_first:
        pattern_lp.solver_parameter("simplex_or_intopt", "exact_simplex_only")
        optval = pattern_lp.solve()
    if reconstruct_rational:
        b = pattern_lp.get_backend()
        optsol = exact_optsol(b)
        k_slope = len(set(optsol) | set([-i for i in optsol]))
        optval = 0
        for i in range(b.ncols()):
            optval += QQ(b.objective_coefficient(i)) * optsol[i]
    else:
        optsol = pattern_lp.get_values(var_slope)
        slopes = list(optsol.values())
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
    r"""
    Reconstruct exact rational basic solution. (solver = ``glp_simplex``)

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: lp = MixedIntegerLinearProgram(solver = 'GLPK', maximization = False)
        sage: x, y = lp[0], lp[1]
        sage: lp.add_constraint(-2*x + y <= 1)
        sage: lp.add_constraint(x - y <= 1)
        sage: lp.add_constraint(x + y >= 2)
        sage: lp.set_objective(x + y)
        sage: lp.solver_parameter("simplex_or_intopt", "simplex_only")
        sage: lp.solve()
        2.0
        sage: lp.get_values(x)
        1.5
        sage: lp.get_values(y)
        0.5
        sage: b = lp.get_backend()
        sage: exact_optsol(b)
        (3/2, 1/2)
    """
    #sage_input(b)
    ncol = b.ncols()
    nrow = b.nrows()
    A = matrix(QQ, ncol + nrow, ncol + nrow, sparse = True)
    for i in range(nrow):
        r = b.row(i)
        for (j, c) in zip(r[0], r[1]):
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
    r"""
    Similar to ``pattern_glpk_lp()``. Uses PPL's ``MIP_Problem`` class for solving lp.

    Exact arithmetics, but too slow.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: optval, optsol, k_slope, v, v_div = pattern_ppl_lp(1, objcoef=(12,55,0))
        sage: k_slope
        6
        sage: optsol
        point(11/47, 7/47, -5/47)
        sage: QQ(optval)
        11
        sage: QQ(v_div)
        47
    """
    pattern = 0
    q, f = pattern_q_and_f(l, pattern)
    vertices_color = pattern_vertices_color(l, pattern, more_ini_additive=more_ini_additive)
    s = pattern_s(l, pattern)
    fn = pattern_fn(l, pattern)
    #objfun = Linear_Expression(0)
    #objfun = sum([16^i * s[i] for i in range(len(s))])
    if objcoef is None:
        objcoef = [16^i for i in range(l+2)]
    objfun = sum([objcoef[i] * Variable(i) for i in range(l+2)])
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
    r"""
    Test the running times for ``pattern_glpk_lp()`` with ``glp_simplex`` + ``glp_exact`` + reconstruction.
    See table in :trac:`18735` comment 7.

    Print l, number of slopes for the optimal solution, running time.
    Return a list of numbers of slopes.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: print("Ignore deprecation warnings for time.clock"); pattern_glpk_test(range(1,4), more_ini_additive=False, exact_arithmetic=False, simplex_first=False, reconstruct_rational=True)
        Ignore...
        1 2 ...
        [2, 2, 2]
    """
    slopes = []
    for l in l_list:
        start_cpu_t = time.clock();
        optval, optsol, k, v_glpk = pattern_glpk_lp(l, more_ini_additive=more_ini_additive, \
                                                    exact_arithmetic=exact_arithmetic, simplex_first = simplex_first, \
                                                    reconstruct_rational=reconstruct_rational);
        cpu_t = time.clock();
        print(l, k, cpu_t - start_cpu_t)
        slopes.append(k)
    return slopes

def pattern_ppl_test(l_list, more_ini_additive=False):
    r"""
    Test the performance of pattern_ppl_lp(). Similar to pattern_glpk_test().
    """
    slopes = []
    for l in l_list:
        start_cpu_t = time.clock();
        try:
            optval, optsol, k, v_ppl, v_div = pattern_ppl_lp(l, more_ini_additive=more_ini_additive);
            cpu_t = time.clock();
            print(l, k, cpu_t - start_cpu_t)
            slopes.append(k)
        except ValueError:
            cpu_t = time.clock();
            print(l, "NA", cpu_t - start_cpu_t)
            slopes.append(-1)
    return slopes

##### The following uses PPL, Variable, Linear_Expression...

def pattern_q_and_f(l, pattern):
    r"""
    Returns the values of q and f for the prescribed pattern.

    pattern == 0 corresponds to the pattern described in the paper.
    In the paper, r is what is called l here.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: pattern_q_and_f(1, pattern=0)
        (58, 29)
    """
    if pattern <= 3:
        f = 18 * l + 11
    elif pattern == 4 or pattern == 5:
        f = 18 * l + 5
    elif pattern == 6:
        f = 18 * l + 13
    elif pattern == 7:
        f = 18 * l + 15
    q = 2 * f
    return q, f

def pattern_additive_vertices(l, pattern):
    r"""
    Return a list of pairs (i, j) of integers i <= j
    such that (i, j)/q is a prescribed additivity.

    These additivities form the additive triangles.

    pattern == 0 corresponds to the pattern described in the paper.
    Remaining patterns are undocumented.
    In the paper, r is what is called l here.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: additive_vertices = pattern_additive_vertices(1, pattern=0)
        sage: len(additive_vertices)
        31
    """
    q, f = pattern_q_and_f(l, pattern)
    f2 = f // 2
    f3 = f // 3
    changed_vertices = []
    # impose some initial green vertices
    for k in range(f3 + 2, f2, 3):
        changed_vertices += [(k, k), (k, k + 1), (k, k + 2), (k - 1, k - 1), (k - 1, k), (k - 1, k + 2), (k + 1, k + 2)]
        if pattern <= 1 or pattern == 4 or pattern == 5 and k < f2 - 2 or pattern == 6 or pattern == 7:
            changed_vertices +=  [(k, k + 3)]
        if pattern == 3:
            changed_vertices += [(k - 1, k + 1), (k + 1, k + 1)]
    changed_vertices += [(1, f2), (f2, f2)]
    if pattern <= 3 or pattern == 6 or pattern == 7:
        changed_vertices += [(f2 - 1, f2 - 1), (f2 - 1, f2), (2, f2 - 1), (2, f2), (3, f2 - 1)]
    if pattern == 3:
        changed_vertices += [(f2 - 1, f2 + 1), (1, f2 - 1), (1, f2 + 1)]
    if pattern == 0:
        changed_vertices += [(f3, f3), (f3, f3 + 2)]
    if pattern == 6:
        changed_vertices += [(f3, f3), (f3, f3 + 1)]
    if pattern == 7:
        changed_vertices += [(f3 - 1, f3 - 1), (f3 - 1, f3), (f3 - 1, f3 + 1), (f3 - 1, f3 + 2), (f3, f3 + 1)]
    for k in range(1, l + 1):
        if pattern <= 3 or pattern == 6 or pattern == 7:
            i = 6 * k - 1
            j = f2 - 3 * k + 1
        elif pattern == 4 or pattern == 5:
            i = 6 * k - 3
            j = f2 - 3 * k + 2 
        changed_vertices += [(i - 1, j), (i - 1, j + 1), (i, j - 1), (i, j + 1), \
                         (i + 1, j - 2), (i + 1, j - 1), (i + 1, j), (i + 1, j + 1), \
                         (i + 2, j - 1), (i + 3, j - 2), (i + 3, j - 1), (i + 4, j - 2)]
        if pattern <= 1 or pattern == 4 or pattern == 5 and k > 1 or pattern == 6 or pattern == 7:
            changed_vertices += [(i - 1, j - 1), (i - 1, j + 2)]
        if pattern == 3:
            changed_vertices += [(i, j), (i + 2, j), (i + 2, j - 2)]
    return changed_vertices

def pattern_more_additive_vertices(l, pattern):
    r"""
    Extra additive points, to be added to ``pattern_additive_vertices`` to reduce the dimension.
    Experimental code.
    """
    q, f = pattern_q_and_f(l, pattern)
    f2 = f // 2
    f3 = f // 3
    more_additivity = [(f3 + 2, q - 4), (f3 + 4, q - 10), (f2 - 2 , f2 - 2), (f3 + 4, q - 12), \
                       (f2 - 8, f2 - 8), (f2 - 17, f2 - 17), (f2 - 23, f2 - 23), \
                       (f3 + 8, q - 22), (f3 + 7, q - 19), (f2 - 14, f2 - 14)]
    ## components merge: l>=3, s1=s2; l>=5, s3=s4; l>=10, s6=s7; l>=13, s8=s9; l>=16, s12=s13.
    l_seuil = [0, 1, 3, 4, 5, 10, 13, 15, 15, 16]
    changed_vertices = [more_additivity[i] for i in range(len(l_seuil)) if l >= l_seuil[i]]
    return changed_vertices

def pattern_vertices_color(l, pattern=0, more_ini_additive=False):
    r"""
    Returns a (q+1)*(q+1) 0-1 array that describes the color of vertices in the prescribed painting
    according to `pattern` and size parameter `l`.
    """
    q, f = pattern_q_and_f(l, pattern)
    vertices_color = initial_vertices_color(q, f)
    changed_vertices = pattern_additive_vertices(l, pattern)
    if pattern == 0 and more_ini_additive:
        changed_vertices += pattern_more_additive_vertices(l, pattern)
    # impose their symmetric vertices too
    for (i, j) in changed_vertices:
        vertices_color[i, j] = vertices_color[q - j, q - i] = 0
    return vertices_color

def pattern_s(l, pattern):
    r"""
    Suppose the prescribed painting (according to pattern and size parameter `l`) has c connected components,
    with slope values q*Variable(0),..., q*Variable(c-1).

    Returns a list s of length f.
    s[i] = fn[i+1]-fn[i] is in {Variable(0),..., Variable(c-1)}, for i=0,...,f-1.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: s = pattern_s(1, pattern=0)
        sage: s[0]
        x0
        sage: len(s)
        29
    """
    if pattern <= 1 or pattern == 6 or pattern == 7:
        s = [Variable(0), Variable(1)]
        for k in range(1, l + 1):
            s += [Variable(k), Variable(k + 1)] * 3
        if pattern == 0:
            s += [Variable(l + 1), Variable(l + 1), Variable(l + 1)]
        elif pattern == 1:
            s += [Variable(l + 1), Variable(l + 2), Variable(l + 1)]
        elif pattern == 6:
            s += [Variable(l + 1), Variable(l + 2), Variable(l + 2), Variable(l + 1)]
        elif pattern == 7:
            s += [Variable(l + 1), Variable(l + 2), Variable(l + 1), Variable(l + 2), Variable(l + 1)]
        for k in range(l, 0, -1):
            s += [Variable(k), Variable(k + 1), Variable(k)]
    elif pattern == 2:
        s = [Variable(0), Variable(1), Variable(1)]
        for k in range(1, l + 1):
            s += [Variable(2 * k), Variable(2 * k + 1), Variable(2 * k), Variable(2 * k + 1), Variable(2 * k), Variable(2 * k)]
        s += [Variable(2 * l + 2)]
        for k in range(l, 0, -1):
            s += [Variable(2 * k), Variable(2 * k + 1), Variable(2 * k)]
        s += [Variable(1)]
    elif pattern == 3:
        s = [Variable(0)] * 3
        for k in range(1, l + 1):
            s += [Variable(k)] * 6
        #s += [Variable(l + 1)]
        s += [-Variable(0)]
        for k in range(l, 0, -1):
            s += [Variable(k)] * 3
        s += [Variable(0)]
    elif pattern == 4:
        s = []
        for k in range(l):
            s += [Variable(k), Variable(k + 1)] * 3
        for k in range(l, 0, -1):
            s += [Variable(k), Variable(k + 1), Variable(k)]
        s += [Variable(0), Variable(1)]
    elif pattern == 5: # inverse
        s = [Variable(l + 2)]
        for k in range(l, 0, -1):
            s += [Variable(k), Variable(k + 1), Variable(k), Variable(k + 1), Variable(k), Variable(k)]
        s += [Variable(0)]
        for k in range(1, l+1):
            s += [Variable(k), Variable(k + 1), Variable(k)]
    return s + [s[0]] + s[-1::-1]

def pattern_fn(l, pattern):
    r"""
    Suppose the prescribed painting (according to pattern and size parameter `l`) has c connected components,
    with slope values q*Variable(0),..., q*Variable(c-1).

    Returns a list fn of length (q+1).
    fn[i] is a Linear_Expression in terms of Variable(0),..., Variable(c-1), for i=0,...,q.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: fn = pattern_fn(1, pattern=0)
        sage: fn[28]
        2*x0+12*x1+14*x2
        sage: len(fn)
        59
        sage: fn[58]
        0
    """
    s = pattern_s(l, pattern)
    f = len(s)
    q = 2 * f
    fn = [Linear_Expression(0)] * (q + 1)
    for k in range(f):
        fn[k + 1] = fn[k] + s[k]
    for k in range(f):
        fn[q - k] = fn[k]
    return fn

def pattern_polytope(vertices_color, fn):
    r"""
    Set up a PPL polytope with respect to the vertices_color.

    fn is a list of ``Linear_Expressions`` in terms of Variable(0),..., Variable(c-1), of length (q+1).
    The polytope has variables {Variable(0),..., Variable(c-1)}.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: vertices_color = pattern_vertices_color(1, pattern=0, more_ini_additive=False)
        sage: fn = pattern_fn(1, pattern=0)
        sage: pattern_polytope(vertices_color, fn)
        A 2-dimensional polyhedron in QQ^3 defined as the convex hull of 6 points
    """
    q = len(fn) - 1
    f = q // 2
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
    polytope = C_Polyhedron(cs)
    return polytope

def pattern_extreme(l, k_slopes, pattern=0, show_plots=False,
                    test_extremality=False, polytope = None, exp_dim = None,
                    more_ini_additive=True, count_components=False, use_sha1=True):
    r"""
    Computes the functions (with >= k_slopes slopes) corresponding to the extreme points of the polytope.

    If polytope is ``None``, set up a polytope corresponding to subadditivities and prescribed additivities,
    according to pattern and size parameter `l` (called `r` in the paper).

    If test_extremality is ``True`` (default is ``False``), check extremality.

    Prints lines:
    NUM-SLOPES  SV (input to ``pattern0_sym_fn``)  DENOMINATOR

    Creates files in the directory 'sym_mode_2d_diagrams/patterns_PATTERN/'
    (this directory is created if it does not exist).

    Returns the max number of slopes found, regardless of extremality test.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: pattern_extreme(3, 4, 0, test_extremality=True, count_components=True)
        #####  ...
        q = 130; f = 65; num_components = 8; num_slopes = 8; divisor = 205; extreme = True
        l = 3; sv = (49, 9, 9, 3, -19)
        h = pattern0_sym_fn(l, sv)
        #####  ...
        q = 130; f = 65; num_components = 6; num_slopes = 6; divisor = 33; extreme = True
        l = 3; sv = (9, 1, 1, 1, -3)
        h = pattern0_sym_fn(l, sv)
        8
    """
    q, f = pattern_q_and_f(l, pattern)
    fn = pattern_fn(l, pattern)
    if polytope is None:
        vertices_color = pattern_vertices_color(l, pattern, more_ini_additive=more_ini_additive)
        polytope = pattern_polytope(vertices_color, fn)
        exp_dim = -1 # no preprocessing, ppl vertex enumeration
    if exp_dim is None:
        # no preprocessing, ppl vertex enumeration
        exp_dim = -1
    v_set = set([])
    vv = []
    nn = []
    destdir = output_dir+"sym_mode_2d_diagrams/"+"patterns_%s/" % pattern
    mkdir_p(destdir)
    logging.disable(logging.INFO)
    #print polytope
    for v in vertex_enumeration(polytope, exp_dim=exp_dim, vetime=False):
        #v.coefficients() is numerator of component's slope value
        v_n = [sum(p*q for p, q in zip(fn[i].coefficients(), v.coefficients())) for i in range(q+1)]
        num = len(set([v_n[i+1] - v_n[i] for i in range(q)]))
        if num >= k_slopes:
            if not tuple(v_n) in v_set:
                v_set.add(tuple(v_n))
                vv.append(v_n)
                nn.append(num)
                #print num, v.coefficients(), v.divisor()
                # Note: Imposing the invariance under x -> 1-x may
                # result in a non-extreme all covered vertex-function.
                h_is_extreme = 'untested'
                h = None
                if test_extremality or count_components or show_plots or use_sha1:
                    # all covered, continuous => extreme iff system of equations has unique solution
                    h = h_from_vertex_values(v_n)
                if test_extremality:
                    h_is_extreme = simple_finite_dimensional_extremality_test(h, show_plots=False, f=None, oversampling=None, order=q, show_all_perturbations=False)
                num_components = 'not computed'
                if count_components:
                    num_components = len(generate_covered_intervals(h))
                if use_sha1:
                    id = h.sha1()
                else:
                    id = len(vv)
                sage_name = "%sq%s_%s.sage" %(num, q, id)
                info = "q = %s; f = %s; num_components = %r; num_slopes = %s; divisor = %s; extreme = %r\n" % (q, f, num_components, num, v.divisor(), h_is_extreme)
                info += "l = %s; " % l
                info += "sv = %s\n" % (tuple(ZZ(x) for x in v.coefficients()),)
                if pattern != 0:
                    info += "v_n = %s\n" % (v_n,)
                    info += "h = h_from_vertex_values(v_n)\n"
                else:
                    info += "h = pattern%s_sym_fn(l, sv)\n" % pattern  # this function only exists for pattern=0
                print("##### ", destdir + sage_name)
                print(info, end=' ')
                with open(destdir + sage_name, "w") as sage_file:
                    print(info, file=sage_file)
                if show_plots: # and h_is_extreme:
                    name = "%sq%s_%s.png" %(num, q, id)
                    g = plot_2d_diagram(h, colorful=True)
                    figsize = 10 * l
                    g.save(destdir + name, figsize = figsize, show_legend=False)
                    print("# Plot saved in", destdir + name)
    logging.disable(logging.NOTSET)
    #return vv, nn
    if len(nn) > 0:
        return max(nn)
    else:
        return 0

def plot_pattern(l, vertices_color=None, more_ini_additive=True, show_plots=True):
    r"""
    Plots prescribled painting (pattern=0) according to size parameter `l` and more_ini_additive.
    If vertices_color is given, use it.
    """
    f = 18 * l + 11
    q = 2 * f
    s=[0, 1]
    for k in range(1, l + 1):
        s += [k, k + 1] * 3
    s += [l + 1] * 3
    for k in range(l, 0, -1):
        s += [k, k + 1, k]
    sk = s + [s[0]] + s[-1::-1]
    nc = max(sk)+1
    sc = [nc - i - 1 for i in sk] + [2*nc - i - 1 for i in sk]
    colors = rainbow(2*nc)
    if vertices_color is None:
        vertices_color = pattern_vertices_color(l, pattern=0, more_ini_additive=more_ini_additive)
    for i in range(q):
        for j in range (i):
            vertices_color[i,j]=vertices_color[j,i]
    for i in range(q):
        vertices_color[q][i]=vertices_color[i][q]=0
    vertices_color[q][q]=0
    g = Graphics()

    vertical_edge = set([])
    horizontal_edge = set([])
    pts = set([])
    for i in range(q):
        for j in range(q):
            #if vertices_color[i,j]==0:
            #    g += points([(i,j)])
            if vertices_color[i+1,j]==0 and vertices_color[i,j+1]==0:
                if vertices_color[i,j]==0:
                    g += polygon([(i,j),(i+1,j),(i,j+1)], color=colors[sc[i]], fill=True)
                    pts.add((i,j));
                    vertical_edge.add((i,j))
                    horizontal_edge.add((i,j))
                if vertices_color[i+1, j+1]==0:
                    g += polygon([(i+1,j+1),(i+1,j),(i,j+1)], color=colors[sc[i]], fill=True)
                    pts.add((i+1, j+1))
                    vertical_edge.add((i+1,j))
                    horizontal_edge.add((i,j+1))
                if vertices_color[i,j] == 1 and vertices_color[i+1, j+1] == 1:
                    g += line([(i+1,j),(i,j+1)], color='springgreen',linestyle='-')
                pts.add((i+1,j));
                pts.add((i,j+1));
            if vertices_color[i,j]==0 and vertices_color[i,j+1]==0 and not ((i,j) in vertical_edge):
                g += line([(i,j),(i,j+1)], color='springgreen',linestyle='-')
                vertical_edge.add((i,j))
                pts.add((i,j))
                pts.add((i,j+1))
            if vertices_color[i,j]==0 and vertices_color[i+1,j]==0 and not ((i,j) in horizontal_edge):
                g += line([(i,j),(i+1,j)], color='springgreen',linestyle='-')
                horizontal_edge.add((i,j))
                pts.add((i,j))
                pts.add((i+1,j))
            if vertices_color[i,j]==0 and not ((i,j) in pts):
                g += points([(i,j)], color='springgreen')
                pts.add((i,j))
    if show_plots:
        g.show(gridlines=[list(range(q+1)),list(range(q+1))],gridlinesstyle=dict(color="grey", linestyle=":"),figsize=10*l)
    else:
        return g

def check_pattern_polytope_has_full_dim(l):
    #if l == 0:
    #    print "l == 0"
    #    return True
    vertices_color = pattern_vertices_color(l, pattern=0, more_ini_additive=False)
    fn = pattern_fn(l, pattern=0)
    polytope = pattern_polytope(vertices_color, fn)
    linearity =  [3,12]+[18]*(l-1)+[14]
    #sumvariables = 100*l*sum([Variable(i) for i in range(l+2)])
    #p = sage.libs.ppl.point(sumvariables, 100*l*(18*l+11))
    #if not polytope.relation_with(p).implies(sage.libs.ppl.Poly_Gen_Relation.subsumes()):
    #    print "gmic(1/2) is not in polytope"
    #    return False
    #for i in range(l+1):
    #    sumvariables += 14*Variable(i) - linearity[i]*Variable(l+1)
    #    p = sage.libs.ppl.point(sumvariables, 100*l*(18*l+11))
    #    if not polytope.relation_with(p).implies(sage.libs.ppl.Poly_Gen_Relation.subsumes()):
    #        print "perturbation %s is not in polytope" % i
    #        return False
    for i in range(1,l+3):
        d = sum(linearity[0:i])
        p = ppl_point(sum([Variable(j) for j in range(i)]), d)
        if not polytope.relation_with(p).implies(Poly_Gen_Relation.subsumes()):
            print("vector [1]*%s+[0]*%s is not in the polytope" % (i, l+2-i))
            return False
    return True
