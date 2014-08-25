destdir = "/Users/mkoeppe/Dropbox/basu-hildebrand-koeppe-papers-for-yuan/survey/"

load("functions.sage")
load("extreme_functions_in_literature.sage")
load("compendium_procedures.sage")
load("survey_examples.sage")

for name in [ 'gmic', 'gj_2_slope', 'gj_2_slope_repeat', 'two_step_mir', 'n_step_mir', 'forward_3_slope', 'backward_3_slope', 'gj_2_slope_limit', 'gj_2_slope_limit_1_1', 'bhk_irrational', 'bccz_counterexample', 'three_slope_limit', 'projected_sequential_merge' ]:
    h = eval(name)()
    g = None
    try:
        g = plot_covered_intervals(h)
    except AttributeError:
        g = plot(h)
    g.save(destdir + "%s.pdf" % name)

#plot_2d_complex(gj_2_slope()).save(destdir + "%s-2d_complex.pdf" % "gj_2_slope")
for name in [ 'bhk_irrational', 'forward_3_slope', 'not_minimal_2' ]:
    h = eval(name)()
    plot_2d_diagram(h).save(destdir + "%s-2d_diagram.pdf" % name)

for name in [ 'not_extreme_1', 'bhk_irrational_extreme_limit_to_rational_nonextreme', 'drlm_not_extreme_1', 'drlm_gj_2_slope_extreme_limit_to_nonextreme' ]:
    h = eval(name)()
    extremality_test(h, show_plots=destdir + "%s-%%s.pdf" % name)

for name in [ 'bhk_irrational_extreme_limit_to_rational_nonextreme' ]:
    for n in [1, 2]:
        h = eval(name)(n)
        extremality_test(h, show_plots=destdir + "%s_%s-%%s.pdf" % (name, n))

for name in [ 'drlm_gj_2_slope_extreme_limit_to_nonextreme' ]:
    for s in [3, 50]:
        h = eval(name)(s)
        extremality_test(h, show_plots=destdir + "%s_%s-%%s.pdf" % (name, s))

