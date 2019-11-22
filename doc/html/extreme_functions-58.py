from cutgeneratingfunctionology.igp import *
h = kzh_extreme_and_weak_facet_but_not_facet()
g = h.plot(show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h.pi))
sphinx_plot(g)