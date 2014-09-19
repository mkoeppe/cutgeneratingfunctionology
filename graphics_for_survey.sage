import igp
from igp import *

destdir = "/Users/mkoeppe/Dropbox/basu-hildebrand-koeppe-papers-for-yuan/survey/"

emitted_names = set()

def emit_tex_sage_command(name):
    if name not in emitted_names:
        print >> sage_commands, '\\pgfkeyssetvalue{/sagefunc/' + name + '}{\\href{\\githubsearchurl?q=\\%22def+' + name.replace('\\', '\\\\') + '(\\%22}{\\sage{' + name.replace('_', '\\underscore{}') + '}}}%)' 
        emitted_names.add(name)

compendium_figsize = 2.6

survey_figsize = 4

orig_ticks_keywords = igp.ticks_keywords
orig_show_plots_figsize = igp.show_plots_figsize
orig_rainbow = sage.plot.colors.rainbow
orig_plot = sage.plot.plot.plot
orig_plot_rescaled_perturbation = igp.plot_rescaled_perturbation
orig_proj_plot_colors = igp.proj_plot_colors

def compendium_ticks_keywords(function, y_ticks_for_breakpoints=False):
    return {'ticks':[[], []], 
            'thickness': 3, 
            #'frame': True, 
            'axes': False, 
            #'fig_tight': False
    }

def survey_ticks_keywords(function, y_ticks_for_breakpoints=False):
    return {'ticks':[[0,1], [0,1]], 
            'gridlines': True,
            #'thickness': 3, 
            #'frame': True, 
            #'axes': False, 
            #'fig_tight': False
    }

def dark_rainbow(num):
    return ['darkblue', 'darkgreen', 'firebrick', 'darkcyan', 'darkmagenta'][:num]

def plot_no_legend(f, *args, **kwds):
    if 'legend_label' in kwds:
        del kwds['legend_label']
    if 'legend_title' in kwds:
        del kwds['legend_title']
    return orig_plot(f, *args, **kwds)

def dont_plot_rescaled_perturbation(*args, **kwds):
    return Graphics()

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
    emit_tex_sage_command(procedure_name)
    plot_something(fn).save(destdir + "%s-from.pdf" % procedure_name, figsize=compendium_figsize)
    if g is None:
        proc = eval(procedure_name)
        g = proc(fn)
    plot_something(g).save(destdir + "%s-to.pdf" % procedure_name, figsize=compendium_figsize)

with open(destdir + "sage-commands.tex", "w") as sage_commands:

    try:

        # override function!
        igp.ticks_keywords = compendium_ticks_keywords
        # override function to get darker colors suitable for print
        igp.rainbow = dark_rainbow
        # override
        igp.plot = plot_no_legend

        ## Compendium procedures table
        procedure_graph('automorphism', gmic())
        procedure_graph('multiplicative_homomorphism', gmic(), multiplicative_homomorphism(gmic(), 3))
        procedure_graph('projected_sequential_merge', multiplicative_homomorphism(gj_forward_3_slope(),-1))
        procedure_graph('restrict_to_finite_group', drlm_not_extreme_1())
        procedure_graph('restrict_to_finite_group_3', drlm_not_extreme_1(), restrict_to_finite_group(drlm_not_extreme_1(), oversampling=3))
        procedure_graph('interpolate_to_infinite_group', restrict_to_finite_group(gmic()))

        ## Compendium tables
        for name in [ 'll_strong_fractional', 'hildebrand_2_sided_discont_1_slope_1', 'hildebrand_2_sided_discont_2_slope_1', 'hildebrand_discont_3_slope_1', 'dr_projected_sequential_merge_3_slope', 'chen_4_slope', 'gmic', 'gj_2_slope', 'gj_2_slope_repeat', 'dg_2_step_mir', 'kf_n_step_mir', 'gj_forward_3_slope', 'drlm_backward_3_slope', 'drlm_2_slope_limit', 'drlm_2_slope_limit_1_1', 'bhk_irrational', 'bccz_counterexample', 'drlm_3_slope_limit', 'dg_2_step_mir_limit', 'rlm_dpl1_extreme_3a', 'hildebrand_5_slope_22_1', 'hildebrand_5_slope_24_1', 'hildebrand_5_slope_28_1' ]:
            emit_tex_sage_command(name)
            h = eval(name)()
            g = plot_something(h)
            g.save(destdir + "%s.pdf" % name, figsize=compendium_figsize)

        for f in ['extremality_test', 'plot_2d_diagram', 'generate_example_e_for_psi_n', 'chen_3_slope_not_extreme', 'psi_n_in_bccz_counterexample_construction', 'gomory_fractional']:
            emit_tex_sage_command(f)

        ## Other figures.

        igp.show_plots_figsize = survey_figsize
        igp.ticks_keywords = survey_ticks_keywords

        #plot_2d_complex(gj_2_slope()).save(destdir + "%s-2d_complex.pdf" % "gj_2_slope")
        for name in [ 'not_extreme_1', 'bhk_irrational_extreme_limit_to_rational_nonextreme', 'drlm_not_extreme_1' ]:
            emit_tex_sage_command(name)
            h = eval(name)()
            extremality_test(h, show_plots=destdir + "%s-%%s.pdf" % name)

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

        # Plot or re-plot some 2d diagrams with a different style
        igp.proj_plot_colors = ['grey', 'grey', 'grey']
        
        for name in [ 'bhk_irrational', 'gj_forward_3_slope', 'not_minimal_2', 'not_extreme_1' ]:
            emit_tex_sage_command(name)
            h = eval(name)()
            # FIXME: Get rid of legend!
            plot_2d_diagram(h, True).save(destdir + "%s-2d_diagram.pdf" % name, figsize=6) # figsize??

    finally:
        igp.plot_rescaled_perturbation = orig_plot_rescaled_perturbation
        igp.show_plots_figsize = orig_show_plots_figsize
        igp.rainbow = orig_rainbow
        igp.ticks_keywords = orig_ticks_keywords
        igp.plot = orig_plot
        igp.proj_plot_colors = orig_proj_plot_colors

os.system("cd %s && (pdflatex -synctex=1 -src-specials -interaction=nonstopmode igp-survey || pdflatex -synctex=1 -src-specials -interaction=nonstopmode igp-survey)" % (destdir,))
