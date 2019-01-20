from cutgeneratingfunctionology.igp import *

def procedure_graph(fn, g):
    G1 = plot_with_colored_slopes(fn, show_legend=False, **only_f_ticks_keywords(fn))
    G2 = plot_with_colored_slopes(g, show_legend=False, **only_f_ticks_keywords(g))
    sphinx_plot(graphics_array([G1, G2]), figsize=(8, 1.5))

h = multiplicative_homomorphism(gj_forward_3_slope(), -1)
procedure_graph(h, projected_sequential_merge(h))