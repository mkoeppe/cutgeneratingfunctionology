
from six.moves import map
from six.moves import range

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
        yticks = sorted(set([ y for limits in function.limits_at_end_points() for y in limits if y is not None ]))
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
    from sage.plot.colors import float_to_html, hsv_to_rgb

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
