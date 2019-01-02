import cutgeneratingfunctionology.igp as igp
from cutgeneratingfunctionology.igp import *

try:
    destdir = algo_paper_output_dir  # defined in config.sage
except Exception:
    #destdir = "survey_graphics/algo_paper_graphics/"
    destdir = "/Users/yzh/Dropbox/basu-hildebrand-koeppe-papers-for-yuan/algo-paper/graphics-for-algo-paper/"

ftype = ".png"
#ftype = ".pdf"

logging.disable(logging.INFO)

igp.plot = plot_no_legend

paper_plot_kwds = { 'dpi': 200, 'fontsize': 2 * 10 };

def plot_kwds_hook_paper(kwds):
    plot_kwds_hook_no_legend(kwds)
    kwds.update(paper_plot_kwds)

igp.plot_kwds_hook = plot_kwds_hook_paper

# Two independently configurable style options for moves diagrams.
igp.show_translations_and_reflections_separately = False
#igp.show_translations_and_reflections_by_color = True     # is the default

#igp.show_covered_components_as_rectangles = True          # is the default
#igp.show_moves_with_discontinuity_markers = True          # is the default; do we want this?

igp.show_plots_figsize = 10

igp.additive_color = 'forestgreen'

# override function to get darker colors suitable for print
igp.rainbow = dark_rainbow_general
