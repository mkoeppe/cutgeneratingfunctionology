from cutgeneratingfunctionology.igp import *

load("survey_graphics/tab_functions.sage")

h = hildebrand_discont_3_slope_1()
h, K = param_piecewise(h)
try:
    logging.info("Running facet_test")
    facet_test(h, known_extreme=True)
except NotImplementedError as e:
    logging.info("Exception: {} (this is normal)".format(e))
    pass

with open('/Users/mkoeppe/w/papers/basu-hildebrand-koeppe-papers/algo-paper/tab_hildebrand_discont_3_slope_1_finite_system.tex', 'w') as f:
    f.write('%% Automatically generated\n')
    caption = r'Homogeneous linear system for $\bar\psi$'
    extra_caption = None
    f.write(tabulate_finite_system(h, label='tab:hildebrand_discont_3_slope_1_finite_system',
                                   caption=caption, extra_caption=extra_caption))
