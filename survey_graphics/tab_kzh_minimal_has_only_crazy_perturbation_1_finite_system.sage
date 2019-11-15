import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *

load("survey_graphics/tab_functions.sage")

parametric=True

h = kzh_minimal_has_only_crazy_perturbation_1(parametric=parametric)

try:
    logging.info("Running facet_test")
    facet_test(h, known_extreme=True)
except NotImplementedError as e:
    logging.info("Exception: {} (this is normal)".format(e))
    pass

# Set names of the slope variables
assert len(h._facet_covered_components) == 2
if open_interval(h.end_points()[0], h.end_points()[1]) in h._facet_covered_components[0]:
    h._slope_names = 'c_3', 'c_1'
elif open_interval(h.end_points()[0], h.end_points()[1]) in h._facet_covered_components[1]:
    h._slope_names = 'c_1', 'c_3'
else:
    raise AssertionError

with open('/Users/mkoeppe/w/papers/basu-hildebrand-koeppe-papers/algo-paper/tab_kzh_minimal_has_only_crazy_perturbation_finite_system.tex', 'w') as f:
    f.write(r'%% Automatically generated.' + '\n')
    caption = r'Homogeneous linear system for determining the restriction of $\bar\pi$ to the non-special intervals in the proof of \autoref{lemma:discontinuous_examples_2}'
    format = r'*3{>{\tiny$}c<{$}@{}}@{\quad}*{%s}{>{\tiny$}c<{$}@{}}@{\quad}*{2}{>{\tiny$}c<{$}@{}}@{\quad}*2{>{\tiny$}c<{$}@{}}' % (len(h._facet_symbolic.basis) - 2 - 2)
    extra_caption = r'. The variables are the values of $\bar\pi$ at breakpoints ($\bullet$) and at midpoints of intervals between breakpoints ($-$), and the slopes $\bar c_1$, $\bar c_3$ of~$\bar\pi$ on the non-special intervals.  Matrix coefficients are abbreviated as ${+} = 1$, ${-} = -1$, and ${\cdot} = 0$.'
    f.write(tabulate_finite_system(h, label='tab:kzh_minimal_has_only_crazy_perturbation_finite_system', format=format, caption=caption, extra_caption=extra_caption))
