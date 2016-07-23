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

complex = SemialgebraicComplex(test_posdef, ['a', 'b'], find_region_type=find_region_type, default_var_bound=(-2,2))

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
    #p = MixedIntegerLinearProgram(maximization=True, solver='InteractiveLP', base_ring=K)
    v = p.new_variable(real=True, nonnegative=True)
    for i in range(len(A)):
        p.add_constraint(v[0]*A[i][0]+v[1]*A[i][1] <= b[i])
        # a = A[i]
        # linexp = 0
        # for j in range(len(a)):
        #     if a[j] != 0:
        #         linexp += v[j]*a[j]
        # p.add_constraint(linexp <= b[i])
    p.set_objective(v[0]*c0+v[1]*c1)
    opt_val = p.solve()
    opt_sol = p.get_values(v).values()
    return tuple(opt_sol)

# def find_opt_sol_color(K, opt_sol):
#     global A, b
#     p = MixedIntegerLinearProgram(maximization=True, solver = "InteractiveLP")
#     v = p.new_variable(real=True, nonnegative=True)
#     for i in range(len(A)):
#         p.add_constraint(v[0]*A[i][0]+v[1]*A[i][1] <= b[i])
#     poly = p.polyhedron()
#     vertices = poly.vertices_list()
#     n = len(vertices)
#     i = vertices.index(opt_sol) #list(opt_sol)
#     vertices_colors = rainbow(n)
#     return vertices_colors[i]


A = [[1,0],[-1,0],[0,1],[0,-1]];
b = [2, -1, 2, -1]
#complex = SemialgebraicComplex(opt_sol_to_varying_obj_func, ['c0', 'c1'], find_region_type=find_opt_sol_color, default_var_bound=(-4,4))
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
#complex = SemialgebraicComplex(opt_sol_to_varying_obj_func, ['c0', 'c1'], find_region_type=find_opt_sol_color, default_var_bound=(-10,10))
complex = SemialgebraicComplex(opt_sol_to_varying_obj_func, ['c0', 'c1'], find_region_type=find_region_type, default_var_bound=(-10,10))
complex.bfs_completion(var_value=[-2, -3], flip_ineq_step=1/3, wall_crossing_method='mathematica')
plot(complex)




