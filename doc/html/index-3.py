from cutgeneratingfunctionology.multirow import *
f = (1/2, 1/2)
B = Polyhedron(vertices=[(0,0),(0,2),(2,0)])
v, g1, g2 = volume_of_lifting_region(B, f, return_plots=True)
sphinx_plot(graphics_array([g1, g2]), figsize=(8, 3))