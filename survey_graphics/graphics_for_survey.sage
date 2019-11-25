import cutgeneratingfunctionology.igp as igp
from cutgeneratingfunctionology.igp import *

destdir = "survey_graphics/survey_graphics/"

emitted_names = set()

def emit_tex_sage_command(name):
    if name not in emitted_names:
        sage_commands.write(r'\pgfkeyssetvalue{/sagefunc/' + name + r'}{\href{\githubsearchurl?q=\%22def+' + name.replace('\\', r'\\') + '(\\%22}{\\sage{' + name.replace(r'_', r'\underscore{}') + '}}}%)' + '\n')
        emitted_names.add(name)

compendium_figsize = 2.6

survey_figsize = 4

orig_ticks_keywords = igp.ticks_keywords
orig_show_plots_figsize = igp.show_plots_figsize
orig_rainbow = sage.plot.colors.rainbow
orig_plot = sage.plot.plot.plot
orig_plot_kwds_hook = igp.plot_kwds_hook
orig_plot_rescaled_perturbation = igp.plot_rescaled_perturbation
orig_proj_plot_colors = igp.proj_plot_colors
orig_check_perturbation_plot_three_perturbations = igp.check_perturbation_plot_three_perturbations

def compendium_ticks_keywords(function, y_ticks_for_breakpoints=False):
    return {'ticks':[[], []], 
            'thickness': 3, 
            #'frame': True, 
            'axes': False, 
            #'fig_tight': False
    }

survey_ticks_keywords = only_f_ticks_keywords

def c7_ticks_keywords(function, y_ticks_for_breakpoints=False):
    xticks = [i/7 for i in range(7+1)]
    xtick_formatter = [ "$%s$" % latex(x) for x in xticks ]
    return {'ticks': [xticks, [1]],
            'tick_formatter': [xtick_formatter, ['$1$']]}

def c37_ticks_keywords(function, y_ticks_for_breakpoints=False):
    xticks = [i/37 for i in range(37+1)]
    xtick_formatter = [ "$%s$" % latex(x) for x in xticks ]
    return {'ticks': [xticks, [1]],
            'tick_formatter': [xtick_formatter, ['$1$']]}

def dark_rainbow(num):
    return ['darkblue', 'darkgreen', 'firebrick', 'darkcyan', 'darkmagenta', 'darkgrey', 'royalblue'][:num]

from sage.plot.colors import float_to_html, hsv_to_rgb

def dark_rainbow_general(n, format='hex'):
    """
    Adapted from sage function rainbow.

    Returns a list of colors sampled at equal intervals over the
    spectrum, from Hue-Saturation-Value (HSV) coordinates (0, 1, 1) to
    (1, 1, 1).  This range is red at the extremes, but it covers
    orange, yellow, green, cyan, blue, violet, and many other hues in
    between.  This function is particularly useful for representing
    vertex partitions on graphs.

    INPUT:

    - ``n`` - a number; the length of the list

    - ``format`` - a string (default: 'hex'); the output format for
      each color in the list; the other choice is 'rgbtuple'

    OUTPUT:

    - a list of strings or RGB 3-tuples of floats in the interval
      [0.0, 1.0]

    EXAMPLES::

        sage: from sage.plot.colors import rainbow
        sage: rainbow(7)
        ['#ff0000', '#ffda00', '#48ff00', '#00ff91', '#0091ff', '#4800ff', '#ff00da']
        sage: rainbow(7, 'rgbtuple')
        [(1.0, 0.0, 0.0), (1.0, 0.8571428571428571, 0.0), (0.2857142857142858, 1.0, 0.0), (0.0, 1.0, 0.5714285714285712), (0.0, 0.5714285714285716, 1.0), (0.2857142857142856, 0.0, 1.0), (1.0, 0.0, 0.8571428571428577)]

    AUTHORS:

    - Robert L. Miller

    - Karl-Dieter Crisman (directly use :func:`hsv_to_rgb` for hues)

    """
    from sage.rings.integer import Integer
    n = Integer(n) # In case n is a Python int and i/n below would give 0!
    R = []

    for i in range(n):
        R.append(tuple(map(float, hsv_to_rgb(i / n, 1, 0.4))))

    if format == 'rgbtuple':
        return R
    elif format == 'hex':
        for j in range(len(R)):
            R[j] = float_to_html(*R[j])
        return R

def plot_something(h):
    g = None
    try:
        g = plot_covered_intervals(h) #, with_legend=False)
    except AttributeError:
        g = plot(h, color='black', **igp.ticks_keywords(h))
    except ValueError:
        g = plot(h, color='black', **igp.ticks_keywords(h))
    return g

def procedure_graph(procedure_name, fn, g=None):
    print("##", procedure_name)
    emit_tex_sage_command(procedure_name)
    plot_something(fn).save(destdir + "%s-from.pdf" % procedure_name, figsize=compendium_figsize)
    if g is None:
        proc = eval(procedure_name)
        g = proc(fn)
    plot_something(g).save(destdir + "%s-to.pdf" % procedure_name, figsize=compendium_figsize)

with open(destdir + "sage-commands.tex", "w") as sage_commands:

    try:

        # override function to get darker colors suitable for print
        igp.rainbow = dark_rainbow_general
        # override
        igp.plot = plot_no_legend
        igp.plot_kwds_hook = plot_kwds_hook_no_legend

        # Graph
        
        print("## graphics_for_survey_poset")
        load('survey_graphics/graphics_for_survey_poset.sage')

        # override function!
        igp.ticks_keywords = compendium_ticks_keywords
        ## Compendium procedures table
        procedure_graph('automorphism', gmic())
        procedure_graph('multiplicative_homomorphism', gmic(), multiplicative_homomorphism(gmic(), 3))
        procedure_graph('projected_sequential_merge', multiplicative_homomorphism(gj_forward_3_slope(),-1))
        h = gj_2_slope(f=3/5, lambda_1=1/2)
        procedure_graph('restrict_to_finite_group', h)
        procedure_graph('restrict_to_finite_group_3', h, restrict_to_finite_group(h, oversampling=3))
        procedure_graph('interpolate_to_infinite_group', restrict_to_finite_group(gmic()))
        procedure_graph('two_slope_fill_in', restrict_to_finite_group(gmic()))

        ## Compendium tables
        for name in [ 'll_strong_fractional', 'hildebrand_2_sided_discont_1_slope_1', 'hildebrand_2_sided_discont_2_slope_1', 'hildebrand_discont_3_slope_1', 'dr_projected_sequential_merge_3_slope', 'chen_4_slope', 'gmic', 'gj_2_slope', 'gj_2_slope_repeat', 'dg_2_step_mir', 'kf_n_step_mir', 'gj_forward_3_slope', 'drlm_backward_3_slope', 'drlm_2_slope_limit', 'drlm_2_slope_limit_1_1', 'bhk_irrational', 'bccz_counterexample', 'drlm_3_slope_limit', 'dg_2_step_mir_limit', 'rlm_dpl1_extreme_3a', 'hildebrand_5_slope_22_1', 'hildebrand_5_slope_24_1', 'hildebrand_5_slope_28_1', 'kzh_7_slope_1', 'kzh_28_slope_1', 'bcdsp_arbitrary_slope', 'bcds_discontinuous_everywhere' ]:
            emit_tex_sage_command(name)
            print("##", name)
            h = eval(name)()
            g = plot_something(h)
            g.save(destdir + "%s.pdf" % name, figsize=compendium_figsize)

        for f in ['extremality_test', 'plot_2d_diagram', 'generate_example_e_for_psi_n', 'chen_3_slope_not_extreme', 'psi_n_in_bccz_counterexample_construction', 'gomory_fractional']:
            emit_tex_sage_command(f)

        ## Other figures.

        igp.show_plots_figsize = survey_figsize
        igp.ticks_keywords = survey_ticks_keywords
        igp.check_perturbation_plot_three_perturbations = False

        #plot_2d_complex(gj_2_slope()).save(destdir + "%s-2d_complex.pdf" % "gj_2_slope")
        for name in [ 'not_extreme_1', 'bhk_irrational_extreme_limit_to_rational_nonextreme', 'zhou_two_sided_discontinuous_cannot_assume_any_continuity']:
            emit_tex_sage_command(name)
            h = eval(name)()
            extremality_test(h, show_all_perturbations=True,
                             show_plots=destdir + "%s-%%s.pdf" % name)

        for name in [ 'bhk_irrational_extreme_limit_to_rational_nonextreme' ]:
            emit_tex_sage_command(name)
            for n in [1, 2]:
                h = eval(name)(n)
                extremality_test(h, show_plots=destdir + "%s_%s-%%s.pdf" % (name, n))

        igp.plot_rescaled_perturbation = dont_plot_rescaled_perturbation

        for name in [ 'drlm_gj_2_slope_extreme_limit_to_nonextreme' ]:
            emit_tex_sage_command(name)
            h = eval(name)()
            extremality_test(h, show_plots=destdir + "%s-%%s.pdf" % name)

        for name in [ 'drlm_gj_2_slope_extreme_limit_to_nonextreme' ]:
            emit_tex_sage_command(name)
            for s in [3, 50]:
                h = eval(name)(s)
                extremality_test(h, show_plots=destdir + "%s_%s-%%s.pdf" % (name, s))

        # Bccz figure
        print("## BCCZ")
        load("survey_graphics/graphics_for_survey_bccz.sage")

        # Plot or re-plot some 2d diagrams with a different style
        igp.proj_plot_colors = ['grey', 'grey', 'grey']
        igp.ticks_keywords = survey_ticks_keywords

        name = 'gmic'
        print("##", name)
        h = eval(name)(2/3)
        g = plot_2d_diagram(h)
        f = find_f(h)
        g += text("$F_1$", (f/3, f/3), axis_coords=False, vertical_alignment='center',
                  horizontal_alignment='center', color='black', fontsize=15)
        g += text("$F_2$", (1 - (1-f)/3, 1 - (1-f)/3), axis_coords=False, vertical_alignment='center',
                  horizontal_alignment='center', color='black', fontsize=15)
        g.save(destdir + "%s-2d_diagram.pdf" % name, figsize = 6)

        igp.ticks_keywords = no_labels_ticks_keywords
        
        for name in [ 'bhk_irrational', 'gj_forward_3_slope', 'not_minimal_2', 'not_extreme_1' ]:
            print("##", name)
            emit_tex_sage_command(name)
            h = eval(name)()
            plot_2d_diagram(h, True).save(destdir + "%s-2d_diagram.pdf" % name, figsize=6) # figsize??

        igp.ticks_keywords = c7_ticks_keywords

        for name in [ 'drlm_not_extreme_1' ]:
            print("##", name)
            emit_tex_sage_command(name)
            h = eval(name)()
            extremality_test(h, show_plots=destdir + "%s-%%s.pdf" % name)
            for oversampling in [1, 4]:
                hr = restrict_to_finite_group(h, oversampling=oversampling)
                plot_something(hr).save(destdir + "%s-restricted-%s.pdf" % (name, oversampling), figsize=survey_figsize)
                extremality_test(hr, show_plots=destdir + "%s-restricted-%s-%%s.pdf" % (name, oversampling))


        #### This one is too complicated to be pretty in the survey's style of figure. Omit.
        igp.ticks_keywords = c37_ticks_keywords

        for name in [ 'kzh_2q_example_1' ]:
            emit_tex_sage_command(name)
            ## h = eval(name)()
            ## extremality_test(h, show_plots=destdir + "%s-%%s.pdf" % name)
            ## for oversampling in [1, 2, 3, 4]:
            ##     hr = restrict_to_finite_group(h, oversampling=oversampling)
            ##     plot_something(hr).save(destdir + "%s-restricted-%s.pdf" % (name, oversampling), figsize=survey_figsize)
            ##     extremality_test(hr, show_plots=destdir + "%s-restricted-%s-%%s.pdf" % (name, oversampling))


    finally:
        igp.plot_rescaled_perturbation = orig_plot_rescaled_perturbation
        igp.show_plots_figsize = orig_show_plots_figsize
        igp.rainbow = orig_rainbow
        igp.ticks_keywords = orig_ticks_keywords
        igp.plot = orig_plot
        igp.plot_kwds_hook = orig_plot_kwds_hook
        igp.proj_plot_colors = orig_proj_plot_colors
        igp.check_perturbation_plot_three_perturbations = orig_check_perturbation_plot_three_perturbations

os.system("cd %s && (pdflatex -synctex=1 -src-specials -interaction=nonstopmode survey_graphics; pdflatex -synctex=1 -src-specials -interaction=nonstopmode compendium_graphics)" % (destdir,))
