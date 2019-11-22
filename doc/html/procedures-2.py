from cutgeneratingfunctionology.igp import *

def procedure_graph(fn, g):
    G1 = plot_with_colored_slopes(fn, show_legend=False, **only_f_ticks_keywords(fn))
    G2 = plot_with_colored_slopes(g, show_legend=False, **only_f_ticks_keywords(g))
    sphinx_plot(graphics_array([G1, G2]), figsize=(8, 1.5))

h = restrict_to_finite_group(gmic(f=QQ('4/5')))
procedure_graph(h, automorphism(h, 2))