from cutgeneratingfunctionology.igp import *
h = gj_2_slope_repeat()
g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
sphinx_plot(g)