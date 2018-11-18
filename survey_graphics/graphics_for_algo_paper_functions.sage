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

show_kwds = copy(paper_plot_kwds)
show_kwds['figsize'] = igp.show_plots_figsize
show_kwds['aspect_ratio'] = 0.3


for name in ['gmic35', 'gmic45auto']:
    h = eval(name)()
    g = plot_with_colored_slopes(h)
    g.save(destdir + "{}_orig".format(name) + ftype, **show_kwds)

    for over in [1, 3]:
        g = plot_with_colored_slopes(restrict_to_finite_group(h, oversampling=over))
        g.save(destdir + "{}_restrict_{}".format(name, over) + ftype,
               **show_kwds)

igp.ticks_keywords = std_ticks_keywords

## Some standard functions with happy moves diagrams!!

def perturb_ticks_keywords(function, x_labels=True):
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
    if x_labels:
        xtick_formatter = [ "$%s$" % latex(x) for x in xticks ]
    else:
        xtick_formatter = [ "" for x in xticks ]
    yticks = []
    ytick_formatter = [ "$%s$" % latex(y) for y in yticks ]
    return {'ticks': [xticks, yticks], \

            'gridlines': True, \
            'tick_formatter': [xtick_formatter, ytick_formatter]}

def plot_it(name):
    h = eval(name)()
    extremality_test(h, show_all_perturbations=True, show_plots=destdir + "%s-%%s" % name + ftype)
    show_plot(h._completion.plot(), destdir + "%s-%%s" % name + ftype, tag='completion-final')
    g = plot_covered_intervals(h, **ticks_keywords(h))
    g.save(destdir + "{}-only-function".format(name) + ftype,
           figsize=igp.show_plots_figsize, aspect_ratio=0.3, **paper_plot_kwds)
    c = h._completion
    zero_perturbation = zero_perturbation_partial_function(c.covered_components,
                                                           c.generate_zero_perturbation_points())
    for index, perturb in enumerate(h._perturbations):
        rescaled = rescale_to_amplitude(perturb, 1/10)
        for x_labels in [False, True, 'zeropert']:
            perturb_kwds = copy(paper_plot_kwds)
            perturb_kwds.update(perturb_ticks_keywords(zero_perturbation, x_labels=x_labels))
            g = plot(rescaled, color='magenta', thickness=2, **perturb_kwds)
            if x_labels is 'zeropert':
                g += plot(zero_perturbation, color='magenta', thickness=3, **perturb_kwds)
                fname = "{}-only-perturbation-zeropert-{}"
            elif x_labels:
                fname = "{}-only-perturbation-{}"
            else:
                fname = "{}-only-perturbation-nolabel-{}"
            g.save(destdir + fname.format(name, index+1) + ftype,
                   figsize=igp.show_plots_figsize, ymin=-0.11, ymax=0.11,
                   aspect_ratio=0.3, **paper_plot_kwds)


# figures for 1/2 linewidth
igp.show_plots_figsize = 5
paper_plot_kwds['fontsize'] = 15
for name in 'equiv7_example_1', 'minimal_no_covered_interval':
    plot_it(name)

# figure for 0.8 linewidth
igp.show_plots_figsize = 8
paper_plot_kwds['fontsize'] = 20
for name in 'equiv7_example_xyz_2', :
    plot_it(name)

# figures for full linewidth
igp.show_plots_figsize = 10
paper_plot_kwds['fontsize'] = 10   ## FIXME: For bigger fontsize, needs manual work to select subset of ticks
for name in 'equiv7_minimal_2_covered_2_uncovered', :
    plot_it(name)

    paper_plot_kwds['fontsize'] = 20
for name in 'equiv7_example_post_3', :
    plot_it(name)

#names += 'equiv7_example_3',  # Robert's proposed example
## names += 'equiv7_example_xyz_1', 'equiv7_example_xyz_3' # XYZ's new examples
## names += 'equiv7_example_post_1', 'equiv7_example_post_2', 'equiv7_example_post_4',
## names += 'kzh_7_slope_1', 'example7slopecoarse2', 'not_extreme_1', 'drlm_not_extreme_1', 'bhk_irrational_extreme_limit_to_rational_nonextreme', 'hildebrand_5_slope_22_1', 'rlm_dpl1_extreme_3a', 'hildebrand_discont_3_slope_1', 'zhou_two_sided_discontinuous_cannot_assume_any_continuity'  # additional, we will select later
