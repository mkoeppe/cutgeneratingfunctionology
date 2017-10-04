# Theorem 2.2. in Finiteness and an explicit description in dimension three
M1 =Polyhedron(vertices=[[0, 0, 0] ,[2, 0, 0],[0, 3, 0],[0, 0, 6]])
M2 =Polyhedron(vertices=[[0, 0, 0],[2, 0, 0],[0, 4, 0],[0, 0, 4]])
M3 =Polyhedron(vertices=[[0, 0, 0],[3, 0, 0],[0, 3, 0], [0, 0, 3]])
M4 = Polyhedron(vertices=[[0, 0, 0], [1, 0, 0], [2, 4, 0], [3, 0, 4]])
M5 = Polyhedron(vertices=[[0, 0, 0], [1, 0, 0], [2, 5, 0], [3, 0, 5]])
M6 = Polyhedron(vertices=[[0, 0, 0], [3, 0, 0], [1, 3, 0], [2, 0, 3]])
M7 = Polyhedron(vertices=[[0, 0, 0], [4, 0, 0], [1, 2, 0], [2, 0, 4]])

# plot(M) + sum([point(int_pt, color='red') for int_pt in bounding_box(M).integral_points()]) 

M8 = Polyhedron(vertices=[[2, 0, 0], [-2, 0, 0], [0, 2, 0], [0, -2, 0], [1, 1, 2]]) # pyramid
M9 = Polyhedron(vertices=[[-1, 0, 0], [0, -1, 0], [2, 0, 0], [0, 2, 0], [1, 1, 3]]) # pyramid
M10 = Polyhedron(vertices=[[1, 0, 0], [0, 1, 0], [-1, -1, 0], [2, 2, 3], [1, 3, 3], [0, 1, 3]]) # prism
#M11 = Polyhedron(vertices=[[1, 0, 0], [-1, 0, 0], [2, 0, 0], [2, 0, 2], [0, 0, 2], [3, 0, 2]]) # prism
M11 = Polyhedron(vertices=[[-1, 0, 0], [0, 0, 2], [0, 2, 0], [1, 0, 0], [1, 2, 2], [2, 0, 2]])
M12 = Polyhedron(vertices=[[-1, 1, 0], [0, 0, 0], [0, 2, 0], [0, 2, 2], [1, 1, 0], [1, 1, 2], [1, 3, 2], [2, 2, 2]])

triangle = Polyhedron(vertices=[[-3/13, 21/13], [1 - 4/10, 3], [3/2, 3/4]])
not_square = Polyhedron(vertices=[[-1/5, 1/2], [2/5, 2], [13/10, 1/2], [2/5, -2/5]])

def volume_is_affine_of_f(B, tests=5):
    for test_i in range(tests):
        f_list = []
        volume_list = []   
        for i in range(B.ambient_dim() + 1):
            # needs B.ambient_dim() + 1 points to construct the hyperplane
            f = pick_f(B)
            f_list.append(f)
            volume_list.append([find_volme_of_B_with_f(B, f)])
            print i
        pts_list = [f + volume for f, volume in zip(f_list, volume_list)]
        print pts_list
        eqns = Polyhedron(vertices=pts_list).equations_list()
        # take the equations from the polyhedron to construct the hyperplane
        print Polyhedron(vertices=pts_list).equations()
        h = Polyhedron(eqns=eqns)
        new_f = pick_f(B)
        new_f_volume = new_f + [find_volme_of_B_with_f(B, new_f)]
        if not h.contains(new_f_volume):
            print "not affine!"
            return False
    print "can't find a new f to show vol(f) is not affine of f"
    return None

volume_is_affine_of_f(triangle)
volume_is_affine_of_f(not_square)
print "test M8\n\n"
M8_result = volume_is_affine_of_f(M8)
print "test M9\n\n"
M9_result = volume_is_affine_of_f(M9)





