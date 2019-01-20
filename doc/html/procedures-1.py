from cutgeneratingfunctionology.igp import *

def procedure_graph(procedure_name, fn, g=None):
    G1 = plot_with_colored_slopes(fn, show_legend=False, **only_f_ticks_keywords(fn))
    if g is None:
        proc = eval(procedure_name)
        g = proc(fn)
    G2 = plot_with_colored_slopes(g, show_legend=False, **only_f_ticks_keywords(g))
    sphinx_plot(graphics_array([G1, G2]), figsize=(8, 1.5))
    #sphinx_plot(G1); sphinx_plot(G2)

procedure_graph('multiplicative_homomorphism', gmic(), multiplicative_homomorphism(gmic(), 3))