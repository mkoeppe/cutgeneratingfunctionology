import igp
from igp import *

try:
    destdir = algo_paper_output_dir  # defined in config.sage
except Exception:
    #destdir = "survey_graphics/algo_paper_graphics/"
    destdir = "/Users/yzh/Dropbox/basu-hildebrand-koeppe-papers-for-yuan/algo-paper/graphics-for-algo-paper/"

ftype = ".png"

logging.disable(logging.INFO)

igp.plot = plot_no_legend
igp.plot_kwds_hook = plot_kwds_hook_no_legend

# Two independently configurable style options for moves diagrams.
igp.show_translations_and_reflections_separately = False
igp.show_translations_and_reflections_by_color = True

igp.show_covered_components_as_rectangles = True
igp.show_moves_with_discontinuity_markers = True   ## do we want this?

igp.show_plots_figsize = 10
