from cutgeneratingfunctionology.igp import *
h = bcds_discontinuous_everywhere()
g = h.plot(show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, ticks=[[1],[1]], tick_formatter=[["$1$", "$1$"]])
sphinx_plot(g)