from cutgeneratingfunctionology.igp import *
h = equiv7_example_xyz_2()
g = plot_completion_diagram(h)
g += point((2.5, 1), color='white', figsize=(8,4))
g.set_legend_options(title="Moves closure of equiv7_example_xyz_2()")
sphinx_plot(g)