destdir = "/Users/mkoeppe/Dropbox/basu-hildebrand-koeppe-papers-for-yuan/survey/"

load("functions.sage")
load("extreme_functions_in_literature.sage")
load("compendium_procedures.sage")
load("survey_examples.sage")

def emit_tex_sage_command(name):
    print >> sage_commands, '\\pgfkeyssetvalue{/sagefunc/' + name + '}{\\href{\\githubsearchurl?q=\\%22def+' + name.replace('\\', '\\\\') + '(\\%22}{\\sage{' + name.replace('_', '\\underscore{}') + '}}}%)' 

def plot_something(h):
    g = None
    try:
        g = plot_covered_intervals(h)
    except AttributeError:
        g = plot(h)
    except ValueError:
        g = plot(h)
    return g

def procedure_graph(procedure_name, fn, g=None):
    emit_tex_sage_command(procedure_name)
    plot_something(fn).save(destdir + "%s-from.pdf" % procedure_name)
    if g is None:
        proc = eval(procedure_name)
        g = proc(fn)
    plot_something(g).save(destdir + "%s-to.pdf" % procedure_name)

with open(destdir + "sage-commands.tex", "w") as sage_commands:
    ## Compendium procedures table
    procedure_graph('multiplicative_homomorphism', gmic(), multiplicative_homomorphism(gmic(), 3))
    procedure_graph('projected_sequential_merge', multiplicative_homomorphism(gj_forward_3_slope(),-1))
    procedure_graph('restrict_to_finite_group', gmic())
    procedure_graph('interpolate_to_infinite_group', restrict_to_finite_group(gmic()))

    ## Compendium tables
    for name in [ 'dr_projected_sequential_merge_3_slope', 'chen_4_slope', 'gmic', 'gj_2_slope', 'gj_2_slope_repeat', 'two_step_mir', 'n_step_mir', 'gj_forward_3_slope', 'drlm_backward_3_slope', 'drlm_2_slope_limit', 'drlm_2_slope_limit_1_1', 'bhk_irrational', 'bccz_counterexample', 'drlm_3_slope_limit', 'dg_2_step_mir_limit', 'hildebrand_5_slope_22_1', 'hildebrand_5_slope_24_1', 'hildebrand_5_slope_28_1' ]:
        emit_tex_sage_command(name)
        h = eval(name)()
        g = plot_something(h)
        g.save(destdir + "%s.pdf" % name)

    for f in ['extremality_test', 'plot_2d_diagram', 'generate_example_e_for_psi_n', 'chen_3_slope_not_extreme', 'rlm_dpl1_fig3_lowerleft', 'psi_n_in_bccz_counterexample_construction']:
        emit_tex_sage_command(f)

    ## Other figures.

    #plot_2d_complex(gj_2_slope()).save(destdir + "%s-2d_complex.pdf" % "gj_2_slope")
    for name in [ 'bhk_irrational', 'gj_forward_3_slope', 'not_minimal_2' ]:
        emit_tex_sage_command(name)
        h = eval(name)()
        plot_2d_diagram(h).save(destdir + "%s-2d_diagram.pdf" % name)

    for name in [ 'not_extreme_1', 'bhk_irrational_extreme_limit_to_rational_nonextreme', 'drlm_not_extreme_1', 'drlm_gj_2_slope_extreme_limit_to_nonextreme' ]:
        emit_tex_sage_command(name)
        h = eval(name)()
        extremality_test(h, show_plots=destdir + "%s-%%s.pdf" % name)

    for name in [ 'bhk_irrational_extreme_limit_to_rational_nonextreme' ]:
        emit_tex_sage_command(name)
        for n in [1, 2]:
            h = eval(name)(n)
            extremality_test(h, show_plots=destdir + "%s_%s-%%s.pdf" % (name, n))

    for name in [ 'drlm_gj_2_slope_extreme_limit_to_nonextreme' ]:
        emit_tex_sage_command(name)
        for s in [3, 50]:
            h = eval(name)(s)
            extremality_test(h, show_plots=destdir + "%s_%s-%%s.pdf" % (name, s))

