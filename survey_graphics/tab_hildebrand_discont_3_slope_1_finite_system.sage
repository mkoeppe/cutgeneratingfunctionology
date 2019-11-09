from cutgeneratingfunctionology.igp import *
h = hildebrand_discont_3_slope_1()
h, K = param_piecewise(h)
try:
    logging.info("Running facet_test")
    facet_test(h, known_extreme=True)
except NotImplementedError as e:
    logging.info("Exception: {} (this is normal)".format(e))
    pass

def format_variable(variable):
    if variable[0] == 'slope of component':
        return r'\bar c_{%s}' % (1 + variable[1], )
    elif variable[0] == 'function value at':
        return r'\bar\pi(%s)' % (latex(variable[1]), )
    else:
        return latex(variable)

def format_label(label):
    return latex(label)

def format_coefficient(coeff):
    if coeff == 0:
        return r' \cdot '
    else:
        return latex(coeff)

def begin_longtable(columns, num_columns=None, format=None, caption=None, extra_caption=None, label=None):
    if num_columns is None:
        num_columns = len(columns)
    if format is None:
        format = r'*{%s}C' % num_columns
    s = []
    s += [r'\def\arraystretch{1.17}']
    s += [r'\begin{longtable}{%s}' % format]
    if caption is not None:
        if extra_caption is None:
            extra_caption = ''
        s += [r'\caption{' + caption + extra_caption + r'}']
    if label is not None:
        s += [r'\label{%s}\\' % label]
    head  = [r'  \toprule']
    head += ['  ' + ' & '.join(columns) + r'\\']
    head += [r'  \midrule']
    s += head
    s += [r'\endfirsthead']
    s += [r'\caption{' + caption + r' (ctd.)}\\']
    s += head
    s += [r'\endhead']
    s += [r'  \bottomrule']
    s += [r'\endfoot']
    return '\n'.join(s)

def end_longtable():
    s = []
    s += [r'\end{longtable}']
    return '\n'.join(s)


def tabulate_finite_system(h):

    s = []
    caption = r'Finite system'
    extra_caption = None
    label = 'tab:hildebrand_discont_3_slope_1_finite_system'
    s += [begin_longtable(columns = [''] + [ format_variable(v) for v in h._facet_symbolic.basis ],
                          caption=caption, extra_caption=extra_caption, label=label)]

    for label, row in zip(h._facet_used_vertices, h._facet_equation_matrix):
        s += [r' & '.join([format_label(label)] + [ format_coefficient(coeff) for coeff in row ]) + r'\\']

    s += [end_longtable()]
    return '\n'.join(s)

with open('/Users/mkoeppe/w/papers/basu-hildebrand-koeppe-papers/algo-paper/tab_hildebrand_discont_3_slope_1_finite_system.tex', 'w') as f:
    f.write('%% Automatically generated\n')
    f.write(tabulate_finite_system(h))
