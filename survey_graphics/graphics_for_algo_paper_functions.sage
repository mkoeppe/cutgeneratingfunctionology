load("survey_graphics/graphics_for_algo_paper_init.sage")

## intro: gmic

igp.show_plots_figsize = 5
std_ticks_keywords = igp.ticks_keywords
igp.ticks_keywords = only_f_ticks_keywords

def gmic45auto():
    h = gmic(f=4/5)
    hr = restrict_to_finite_group(h)
    ha = automorphism(hr, 2)
    return interpolate_to_infinite_group(ha)

def gmic35():
    return gmic(f=3/5)

for name in ['gmic35', 'gmic45auto']:
    h = eval(name)()
    g = plot_with_colored_slopes(h)
    g.save(destdir + "{}_orig.png".format(name),
           figsize=igp.show_plots_figsize, aspect_ratio=0.3)
    for over in [1, 3]:
        g = plot_with_colored_slopes(restrict_to_finite_group(h, oversampling=over))
        g.save(destdir + "{}_restrict_{}.png".format(name, over),
               figsize=igp.show_plots_figsize, aspect_ratio=0.3)

igp.ticks_keywords = std_ticks_keywords

## Some standard functions with happy moves diagrams!!

def perturb_ticks_keywords(function):
    """
    Compute ``plot`` keywords for displaying the ticks.
    """
    xticks = function.end_points()
    try:
        f = find_f(function, no_error_if_not_minimal_anyway=True)
        if f is not None and not f in xticks:
            xticks.append(f)
    except ValueError:
        pass
    xticks = uniq(xticks)
    xtick_formatter = [ "$%s$" % latex(x) for x in xticks ]
    yticks = []
    ytick_formatter = [ "$%s$" % latex(y) for y in yticks ]
    return {'ticks': [xticks, yticks], \

            'gridlines': True, \
            'tick_formatter': [xtick_formatter, ytick_formatter]}


igp.show_plots_figsize = 10

names = 'equiv7_example_1', 'minimal_no_covered_interval', 'equiv7_minimal_2_covered_2_uncovered' # essential
names += 'kzh_7_slope_1', 'example7slopecoarse2', 'not_extreme_1', 'drlm_not_extreme_1', 'bhk_irrational_extreme_limit_to_rational_nonextreme', 'hildebrand_5_slope_22_1', 'rlm_dpl1_extreme_3a', 'hildebrand_discont_3_slope_1', 'zhou_two_sided_discontinuous_cannot_assume_any_continuity'  # additional, we will select later

for name in names:
    h = eval(name)()
    extremality_test(h, show_all_perturbations=True, show_plots=destdir + "%s-%%s.pdf" % name)
    show_plot(h._completion.plot(), destdir + "%s-%%s.pdf" % name, tag='completion-final')
    g = plot_covered_intervals(h, **ticks_keywords(h))
    g.save(destdir + "{}-only-function.png".format(name),
           figsize=igp.show_plots_figsize, aspect_ratio=0.3)
    for index, perturb in enumerate(h._perturbations):
        rescaled = rescale_to_amplitude(perturb, 4/10)
        g = plot(rescaled, color='magenta', **perturb_ticks_keywords(h))
        g.save(destdir + "{}-only-perturbation-{}.png".format(name, index+1),
               figsize=igp.show_plots_figsize, ymin=-0.5, ymax=0.5,
               aspect_ratio=0.3)
