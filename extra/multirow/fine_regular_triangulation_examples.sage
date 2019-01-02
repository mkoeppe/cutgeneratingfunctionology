##########################################
# Tests For Fine Regular Triangulation
##########################################

# Examples from literature

# Example 2.2.5 in Jesus Triangulation
outside_circle = Polyhedron(vertices=[[4, 0, 0], [0, 4, 0], [0, 0, 4]])
inside_circle = Polyhedron(vertices=[[2, 1, 1], [1, 2, 1], [1, 1, 2]])
pa = PolyhedralArrangement([outside_circle, inside_circle])
pa.plot()
trigs = pa.fine_regular_triangulation()
pa.plot(trigs)


# Example 2.2.9 in Jesus Triangulation
square = Polyhedron(vertices=[[0, 0], [3, 0], [0, 3], [3, 3]])
point = Polyhedron(vertices=[[1, 1]])
pa = PolyhedralArrangement([square, point])
pa.plot()
trigs = pa.fine_regular_triangulation()
pa.plot(trigs)


# cube and a square (dissection)
square = Polyhedron(vertices=[[0, 0, 1/2], [1, 0, 1/2], [0, 1, 1/2], [1, 1, 1/2]])
cube = Polyhedron(vertices=[[1, 0, 0], [1, 1, 0], [2, 1, 0], [2, 0, 0], [1, 0, 1], [1, 1, 1], [2, 1, 1], [2, 0, 1]])
pa = PolyhedralArrangement([square, cube])
pa.plot()
trigs = pa.fine_regular_triangulation()
pa.plot(trigs)

# self constructed examples

# letter_C
top = Polyhedron(vertices=[[0, 0, 1], [0, 1, 1], [1, 0, 1], [1, 1, 1]])
botton = Polyhedron(vertices=[[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]])
left = Polyhedron(vertices=[[0, 0, 0], [1, 0, 0], [1, 0, 1], [0, 0, 1]])
letter_C = PolyhedralArrangement({top, botton, left})
letter_C.plot()
trigs = letter_C.fine_regular_triangulation()
letter_C.plot(trigs)

# box
box = Polyhedron(vertices=[[0, 0, 1], [0, 1, 1], [1, 0, 1], [1, 1, 1], [0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]])
pa = PolyhedralArrangement({box})
pa.plot()
trigs = pa.fine_regular_triangulation()
pa.plot(trigs)

# square with a cross
square = Polyhedron(vertices=[[0, 0], [0, 1], [1, 0], [1, 1]])
horizontal_line = Polyhedron(vertices=[[0, 1/2], [1, 1/2]])
vertical_line = Polyhedron(vertices=[[1/2, 0], [1/2, 1]])
pa = PolyhedralArrangement({square, horizontal_line, vertical_line})
pa.plot()
trigs = pa.fine_regular_triangulation()
pa.plot(trigs)
