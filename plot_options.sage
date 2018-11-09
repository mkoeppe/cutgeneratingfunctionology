
def plot_no_legend(f, *args, **kwds):
    # really should rather use plot_kwds_hook everywhere in functions.sage
    plot_kwds_hook_no_legend(kwds)
    return sage.plot.plot.plot(f, *args, **kwds)

def plot_kwds_hook_no_legend(kwds):
    if 'legend_label' in kwds:
        del kwds['legend_label']
    if 'legend_loc' in kwds:
        del kwds['legend_loc']
    if 'legend_title' in kwds:
        del kwds['legend_title']

def dont_plot_rescaled_perturbation(*args, **kwds):
    return Graphics()

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
        ytick_formatter = [ latex_formatter_or_empty(y) for y in yticks]
    return {'ticks': [xticks, yticks],
            'gridlines': True,
            'tick_formatter': [xtick_formatter, ytick_formatter]}

def only_f_ticks_keywords(function, y_ticks_for_breakpoints=False):
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
