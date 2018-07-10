import igp
from igp import *

destdir = "survey_graphics/extreme_notes_graphics/"

emitted_names = set()

def emit_tex_sage_command(name):
    if name not in emitted_names:
        print >> sage_commands, '\\pgfkeyssetvalue{/sagefunc/' + name + '}{\\href{\\githubsearchurl?q=\\%22def+' + name.replace('\\', '\\\\') + '(\\%22}{\\sage{' + name.replace('_', '\\underscore{}') + '}}}%)' 
        emitted_names.add(name)

compendium_figsize = 2.6

survey_figsize = 4

diagram_figsize = 8

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

def survey_ticks_keywords(function, y_ticks_for_breakpoints=False):
    if y_ticks_for_breakpoints:
        return {'ticks':[[0,find_f(function),1], [0,find_f(function), 1]], 
                'tick_formatter': [['$0$', '$f$', '$1$'], ['$0$', '$f$', '$1$']],
                'gridlines': True,
                #'thickness': 3, 
                #'frame': True, 
                #'axes': False, 
                #'fig_tight': False
                }
    else:
        return {'ticks':[[0,find_f(function),1], [0,1]], 
                'tick_formatter': [['$0$', '$f$', '$1$'], ['$0$', '$1$']],
                'gridlines': True,
        }
def rlm_dpl1_extreme_3a_ticks_keywords(function, y_ticks_for_breakpoints=False):
    f = find_f(function)
    if y_ticks_for_breakpoints:
        return {'ticks':[[0, f, (1+f)/2, 1], [0, f, (1+f)/2, 1]], 
                'tick_formatter': [['$0$', '$f$', '$b$', '$1$'], ['$0$', '$f$', '$b$', '$1$']],
                'gridlines': True,
                'fontsize': 15,
                #'thickness': 3, 
                #'frame': True, 
                #'axes': False, 
                #'fig_tight': False
                }
    else:
        return {'ticks':[[0, f, (1+f)/2, 1], [0,1/2, 1]], 
                'tick_formatter': [['$0$', '$f$', '$b$', '$1$'], [ "$%s$" % latex(x) for x in [0, 1/2, 1]]],
                'gridlines': True,
                'fontsize': 15,
                }

def chen_3_slope_not_extreme_ticks_keywords(function, y_ticks_for_breakpoints=False):
    if y_ticks_for_breakpoints:
        return {'ticks':[[0,find_f(function), 1], [0,find_f(function), 1]], 
                'tick_formatter': [['$0$', '$f$', '$1$'], ['$0$', '$f$', '$1$']],
                'gridlines': True,
                'fontsize': 15,
                #'thickness': 3, 
                #'frame': True, 
                #'axes': False, 
                #'fig_tight': False
                }
    else:
        xticks = function.end_points()
        yticks = [0, 1/3, 2/3, 1]
        ytick_formatter = [ "$%s$" % latex(y) for y in yticks ]
        return {'ticks':[xticks, yticks], 
                'tick_formatter': [['$0$', '', '', '', '', '$f$', '', '','','','$1$'], \
                                   ytick_formatter ],
                'gridlines':[[0,find_f(function), 1], [0, 1]],
                'fontsize': 15
               }
def drlm_backward_3_slope_ticks_keywords(function, y_ticks_for_breakpoints=False):
    f = find_f(function)
    b = function.end_points()[2]
    if y_ticks_for_breakpoints:
        xticks = [0, f, b, 2*b, 1+f-2*b, 1+f-b, 1]
        xtick_formatter = ['$0$', '$f$', '$b$', '$2b$', '$1+f-2b$', '', '$1$']
        ytick_formatter = ['$0$', '$f$', '$b$', '$2b$', '$1+f-2b$', '$1+f-b$', '$1$']
        return {'ticks':[xticks, xticks], 
                'tick_formatter': [xtick_formatter, ytick_formatter],
                'gridlines': [[0, f, 1], [0,f,1]],
                'fontsize': 15,
                #'thickness': 3, 
                #'frame': True, 
                #'axes': False, 
                #'fig_tight': False
                }
    else:
        xticks = [0, f, b, (1+f)/4, 1+f-b,  1]
        xtick_formatter = [ '$0$', '$f$', '$b$', "$%s$" % LatexExpr(r"\frac{1+f}{4}"),"", '$1$' ] 
        return {'ticks':[xticks, [0,1]], 
                'tick_formatter': [xtick_formatter, ['$0$', '$1$']],
                'gridlines': [[0,f,1], [0,1]],
                'fontsize': 15,
               }
               
def c7_ticks_keywords(function, y_ticks_for_breakpoints=False):
    xticks = [i/7 for i in range(7+1)]
    xtick_formatter = [ "$%s$" % latex(x) for x in xticks ]
    return {'ticks': [xticks, [1]],
            'tick_formatter': [xtick_formatter, ['$1$']]}

def latex_formatter_or_empty(x, labels_list = [0, 1]):
    if x in labels_list:
        return "$%s$" % latex(x)
    else:
        return ""

def latex_formatter_or_f_or_empty(x, f, labels_list = [0, 1]):
    if x == f:
        return "$f$"
    elif x in labels_list:
        return "$%s$" % latex(x)
    else:
        return ""

def no_labels_ticks_keywords(function, y_ticks_for_breakpoints=False):
    xticks = function.end_points()
    xtick_formatter = [ latex_formatter_or_f_or_empty(x, find_f(function)) for x in xticks ]
    #xtick_formatter = 'latex'  # would not show rationals as fractions
    ytick_formatter = None
    if y_ticks_for_breakpoints:
        yticks = xticks
        ytick_formatter = xtick_formatter
    else:
        #yticks = 1/5
        yticks = uniq([ y for limits in function.limits_at_end_points() for y in limits if y is not None ])
        ytick_formatter = [ latex_formatter_or_empty(y) ]
    return {'ticks': [xticks, yticks],
            'gridlines': True,
            'tick_formatter': [xtick_formatter, ytick_formatter]}

def dark_rainbow(num):
    return ['darkblue', 'darkgreen', 'firebrick', 'darkcyan', 'darkmagenta'][:num]

def plot_something(h):
    g = None
    try:
        g = plot_covered_intervals(h) #, with_legend=False)
    except AttributeError:
        g = plot(h, color='black', **igp.ticks_keywords(h))
    except ValueError:
        g = plot(h, color='black', **igp.ticks_keywords(h))
    return g


#with open(destdir + "sage-commands.tex", "w") as sage_commands:

try:

        # override function to get darker colors suitable for print
        #igp.rainbow = dark_rainbow
        # override
        igp.plot = plot_no_legend
        igp.plot_kwds_hook = plot_kwds_hook_no_legend

        ## Other figures.

        igp.show_plots_figsize = survey_figsize
        igp.ticks_keywords = rlm_dpl1_extreme_3a_ticks_keywords
        name = 'rlm_dpl1_extreme_3a'
        h = eval(name)()
        g = plot_covered_intervals(h)
        g.save(destdir + "%s-covered_intervals.pdf" % name, figsize = survey_figsize)
        
        g = plot_2d_diagram(h)
        f = find_f(h)
        b = (1+f)/2
        g += text("$F_1$", (f/3, f/3), axis_coords=False, vertical_alignment='center', \
                  horizontal_alignment='center', color='black', fontsize=15)
        g += text("$F_2$", ((1-b)/4 + b, (1-b)/4 + b), axis_coords=False, vertical_alignment='center', \
                  horizontal_alignment='center', color='black', fontsize=15)
        g += text("$F_3$", (f + (b-2*f)/3, f + (b-2*f)/3), axis_coords=False, vertical_alignment='center', \
                  horizontal_alignment='center', color='black', fontsize=15)
        g.save(destdir + "%s-2d_diagram.pdf" % name, figsize = diagram_figsize)   
        
        name = 'chen_3_slope_not_extreme'
        h = eval(name)()
        igp.ticks_keywords = chen_3_slope_not_extreme_ticks_keywords
        igp.check_perturbation_plot_three_perturbations = False
        extremality_test(h, show_plots=destdir + "%s-%%s.pdf" % name)
        g = plot_covered_intervals(h)
        bkpt = h.end_points()
        xticks = ['$0$', '$A$', '$B$', '$BB$', '$AA$', '$f$', '$CC$', '$DD$','$D$','$C$','$1$']
        for i in range(len(xticks)):
            if xticks[i] != '$0$' and  xticks[i] != '$f$' and xticks[i] != '$1$':
                g += text(xticks[i], (bkpt[i], 0), axis_coords=True, vertical_alignment='top', \
                  horizontal_alignment='center', color='black', fontsize=8)   
        g.save(destdir + "%s-covered_intervals.pdf" % name, figsize = survey_figsize)
        #plot_2d_diagram(h).save(destdir + "%s-2d_diagram.pdf" % name, figsize = diagram_figsize)

        igp.ticks_keywords = drlm_backward_3_slope_ticks_keywords
        name = 'drlm_backward_3_slope'
        h = eval(name)()
        f = find_f(h)
        b = h.end_points()[2]
        #c = (1+f-b)/2
        
        g = plot_covered_intervals(h)
        g += text("$%s$" % LatexExpr(r"1+f-b"), (1+f-b, 0), axis_coords=True, vertical_alignment='top', \
                  horizontal_alignment='right', color='black', fontsize=15)
        g.save(destdir + "%s-covered_intervals.pdf" % name, figsize = survey_figsize)

        g = plot_2d_diagram(h)
        g += text("$F_1$", (0, 0), axis_coords=False, vertical_alignment='bottom', \
                  horizontal_alignment='left', color='black', fontsize=15)
        g += text("$F_2$", (b, 1), axis_coords=False, vertical_alignment='top', \
                  horizontal_alignment='right', color='black', fontsize=15)
        g += text("$F_3$", (2*b, 2*b), axis_coords=False, vertical_alignment='bottom', \
                  horizontal_alignment='left', color='black', fontsize=15)
        g.save(destdir + "%s-2d_diagram.pdf" % name, figsize = diagram_figsize)
        
        igp.ticks_keywords = survey_ticks_keywords
            
        name= 'kf_n_step_mir'
        f = 2/3 
        q = 4
        eta = 1/1000
        n = 2
        e = generate_example_e_for_psi_n(f=f, n=n, q=q, eta=eta)
        a1 = (f + e[0])/2
        a2 = (a1 - e[0] + e[1])/2
        h = psi_n_in_bccz_counterexample_construction(f=f, e=e)
        xticks = [a2, a1, f, 1]
        xtick_formatter = ["$a_2$", "$a_1$", "$f$", "$a_0 = 1$"]
        xticks += h.end_points()
        xtick_formatter += [ "" for x in h.end_points() ]
        yticks = [1]
        ytick_formatter = ["$1$"]
        g = plot(h, ticks=[xticks, yticks],
             tick_formatter=[xtick_formatter, ytick_formatter],
             gridlines=True, color='black')
        g += text("$\\epsilon_2$", (a2 - e[1]/2, 0), axis_coords=False, fontsize=15, \
                              vertical_alignment='bottom', horizontal_alignment='center', color='black')
        g += text("$\\epsilon_1$", (a1 - e[0]/2, 0), axis_coords=False, fontsize=15,\
                              vertical_alignment='bottom', horizontal_alignment='center', color='black')
        g += text("$\\epsilon_2$", (a1 + a2 - e[1]/2, 0), axis_coords=False, fontsize=15,\
                              vertical_alignment='bottom', horizontal_alignment='center', color='black')
        g.save(destdir + "%s-lemma5_1.pdf" % name, figsize = 5, fontsize=15)
    
        plot(bccz_counterexample(), color='black').save(destdir + "%s-function.pdf" % "bccz_counterexample")

finally:
        igp.plot_rescaled_perturbation = orig_plot_rescaled_perturbation
        igp.show_plots_figsize = orig_show_plots_figsize
        igp.rainbow = orig_rainbow
        igp.ticks_keywords = orig_ticks_keywords
        igp.plot = orig_plot
        igp.plot_kwds_hook = orig_plot_kwds_hook
        igp.proj_plot_colors = orig_proj_plot_colors
        igp.check_perturbation_plot_three_perturbations = orig_check_perturbation_plot_three_perturbations

os.system("cd %s && (pdflatex -synctex=1 -src-specials -interaction=nonstopmode extreme_notes_graphics)" % (destdir,)) 
