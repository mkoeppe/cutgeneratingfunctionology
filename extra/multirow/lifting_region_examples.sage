triangle1 = Polyhedron(vertices=[[-3/13, 21/13], [1 - 4/10, 3], [3/2, 3/4]])
triangle2 = Polyhedron(vertices=[[0, 0], [0, 2], [2, 0]])

diamond = Polyhedron(vertices=[[-1/2, 1/2], [1/2, 3/2], [1/2, -1/2], [3/2, 1/2]])
not_square = Polyhedron(vertices=[[-1/5, 1/2], [2/5, 2], [13/10, 1/2], [2/5, -2/5]])
not_maximal_1 = Polyhedron(vertices=[[-1/5, 1/2], [2/5, 2], [13/10, 1/2], [2/5, -2/5]])
not_maximal_2 = Polyhedron(vertices=[[-1/5, 1/2], [2/5, 2], [13/10, 1/2]])
not_maximal_3 = Polyhedron(vertices=[[-1/5, 1/2], [2/5, 2], [13/10, 1/2], [1, 0]])
not_maximal_4 = Polyhedron(vertices=[[0, 0], [0, 1], [1, 0]])
point = Polyhedron(vertices=[[0, 0]])

not_maximal = Polyhedron(vertices=[[0, 1], [-3, 0], [3, 0]])
maximal = Polyhedron(vertices=[[-3/5, 6/5], [-3, 0], [3, 0]])
not_maximal_1 = Polyhedron(vertices=[[0, 1], [-2, 0], [2, 0]])

lifting_graphics(not_maximal)
lifting_graphics(not_maximal_1)
lifting_graphics(maximal)

# pyramid whose base has 10 lattice points.
B = Polyhedron(vertices=((-37/24, 1/2, 0),(1/8, -7/6, 0),(1/8, 13/6, 0),(83/24, 1/2, 0),(13/24, 1/2, 2)))
