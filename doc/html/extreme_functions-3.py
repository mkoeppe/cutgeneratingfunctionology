from cutgeneratingfunctionology.igp import *
h = bccz_counterexample()
h = psi_n_in_bccz_counterexample_construction(e=generate_example_e_for_psi_n(n=7))
g = plot(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=1, **only_f_ticks_keywords(h))
sphinx_plot(g)