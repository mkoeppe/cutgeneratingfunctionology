from cutgeneratingfunctionology.igp import *
h = discontinuous_facets_paper_example_psi_prime()
g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
sphinx_plot(g)