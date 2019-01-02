def test_posdef(a, b):
    return is_positive_definite(matrix([[a, b], [b, 1/4]]))

def is_positive_definite(self):
    # copyed from "sage/matrix/matrix2.pyx":12383, modified L12579.
    imaginary = True
    a = self.base_ring().an_element()
    try:
        a.conjugate()
    except AttributeError:
        imaginary = False
    if imaginary:
        if not self.is_hermitian():
            return False
    else:
        if not self.is_symmetric():
            return False
    try:
        if imaginary:
            _, d = self._indefinite_factorization('hermitian', check=False)
        else:
            _, d = self._indefinite_factorization('symmetric', check=False)
    except ValueError as e:
        # a zero singular leading principal submatrix is one
        # indicator that the matrix is not positive definite
        if str(e).find('leading principal submatrix is singular') != -1:
            return False
        else:
            raise ValueError(e)
    # Now have diagonal entries (hopefully real) and so can
    # test with a generator (which will short-circuit)
    # positive definite iff all entries of d are positive
    try:
        posdef = all( x > 0 for x in d )  # change RR(x)>0 to x>0, x is always real.
    except TypeError:
        universe = Sequence(d).universe()
        msg = "cannot convert computations from {0} into real numbers"
        raise TypeError(msg.format(universe))
    return posdef

complex = SemialgebraicComplex(test_posdef, ['a', 'b'], find_region_type=return_result, default_var_bound=(-2,2))

complex.shoot_random_points(10)
#complex.bfs_completion(var_value=[-1,1], check_completion=True, goto_lower_dim=False)  # 3 components
#complex.bfs_completion(var_value=[-1,1], check_completion=False, goto_lower_dim=True)  # 4 components
#complex.bfs_completion(var_value=[-1,1], check_completion=True, goto_lower_dim=True)   # 5 components
#complex.bfs_completion(var_value=[-1,1], wall_crossing_method='mathematica', check_completion=True, goto_lower_dim=False)  # 3 components
#complex.bfs_completion(var_value=[-1,1], wall_crossing_method='mathematica', check_completion=False, goto_lower_dim=True)  # 5 components
#complex.bfs_completion(var_value=[-1,1], wall_crossing_method='mathematica', check_completion=True, goto_lower_dim=True)   # 5 components
plot(complex)

def opt_sol_to_varying_obj_func(c0, c1):
    global A, b
    K = c0.parent()
    p = MixedIntegerLinearProgram(maximization=True, base_ring=K)
    v = p.new_variable(real=True, nonnegative=True)
    for i in range(len(A)):
        p.add_constraint(v[0]*A[i][0]+v[1]*A[i][1] <= b[i])
    p.set_objective(v[0]*c0+v[1]*c1)
    opt_val = p.solve()
    opt_sol = p.get_values(v).values()
    return opt_sol


A = [[1,0],[-1,0],[0,1],[0,-1]];
b = [2, -1, 2, -1]
complex = SemialgebraicComplex(opt_sol_to_varying_obj_func, ['c0', 'c1'], find_region_type=find_region_type, default_var_bound=(-4,4))
complex.bfs_completion(var_value=[2, -3], flip_ineq_step=1/2, wall_crossing_method='mathematica')
plot(complex)

#copmlex.bfs_completion(var_value=[2, 3])
#sage: complex.bfs_completion(var_value=[2, 3],goto_lower_dim=True)
#sage: complex.bfs_completion(var_value=[2, 3], wall_crossing_method='mathematica')



#sage: A = [[1,-3],[-3,1],[1,1]]
#sage: b = [0,0,4]


A = [[0,-1],[-4,1],[-1,2],[4,1]]
b = [0,0,5,16]
complex = SemialgebraicComplex(opt_sol_to_varying_obj_func, ['c0', 'c1'], find_region_type=find_region_type, default_var_bound=(-10,10))
complex.bfs_completion(var_value=[-2, -3], flip_ineq_step=1/3, wall_crossing_method='mathematica')
plot(complex)



###################

def subdetmip(a, delta=None, subsystem=False, base_ring=QQ, solver=None, integer=True):
    "If delta is None, optimize for it"
    m = MixedIntegerLinearProgram(maximization=False, base_ring=base_ring, solver=solver)
    x = m.new_variable(integer=integer)
    n = len(a)
    # create variables:
    [ x[i] for i, ai in enumerate(a) ]
    if delta is None:
        delta = m.new_variable()[0]
        m.add_constraint(sum(x[i] for i, ai in enumerate(a)) >= 1)
        m.set_objective(delta)
    if subsystem:
        m.add_constraint(-delta <= x[n-1] <= delta)
        for i, ai in enumerate(a):
            if i < n-1:
                m.add_constraint(-delta <= x[n-1] * ai - x[i] * a[n-1] <= delta)
    else:
        for i, ai in enumerate(a):
            m.add_constraint(-delta <= x[i] <= delta)
        for i, ai in enumerate(a):
            for j, aj in enumerate(a):
                if i < j:
                    # matrix([[ai, aj], [x[i], x[j]]]).det() gives "unable to find a common ring for all elements"
                    m.add_constraint(-delta <= x[j] * ai - x[i] * aj <= delta)
    return m

volume_engine='lrs'

def subdetpoly(a, delta, subsystem=False, base_ring=QQ, integer=False):
    return subdetmip(a, delta, subsystem=subsystem, base_ring=base_ring, integer=integer).polyhedron()

def volsubdetpoly(a1=1, a2=1, a3=1):
    if not bool(1 <= a1 <= a2 <= a3):
        raise ValueError ("Assumption 1 <= a1 <= a2 <= a3 is not satisfied")
    K = a1.parent()
    P = subdetpoly([a1, a2, a3], 1, base_ring=K, integer=False)
    return P.volume()


sage: complex = SemialgebraicComplex(volsubdetpoly, ['a1','a2','a3'], max_iter=0, find_region_type=result_symbolic_expression, default_var_bound=(0,10))

K.<a1,a2,a3>=ParametricRealField([1,2,3])
volsubdetpoly(a1, a2, a3)
