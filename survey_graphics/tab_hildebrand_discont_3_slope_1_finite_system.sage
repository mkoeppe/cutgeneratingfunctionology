from cutgeneratingfunctionology.igp import *

load("survey_graphics/tab_functions.sage")

h = hildebrand_discont_3_slope_1()
h, K = param_piecewise(h)
try:
    logging.info("Running facet_test")
    facet_test(h, known_extreme=True)
except NotImplementedError as e:
    logging.info("Exception: {} (this is normal)".format(e))
    pass

from cutgeneratingfunctionology.spam.real_set import RealSet
all_triples = list(generate_triples_with_projections_intersecting(h, RealSet([0, 1])))
veps_to_face_dict = {}
for T in all_triples:
    hF = Face(T)
    if is_additive_face_sans_limits(h, hF):
        for v in hF.vertices:
            for eps in generate_containing_eps_triple(v, hF.minimal_triple, with_limits=False):
                veps = (v[0], v[1], v[0]+v[1], eps[0], eps[1], eps[2])
                veps_to_face_dict[veps] = hF

def format_variable(variable):
    if variable[0] == 'slope of component':
        return r'\bar c_{%s}' % (1 + variable[1], )
    elif variable[0] == 'function value at':
        return r'\bar\pi(%s)' % (latex(variable[1]), )
    else:
        return latex(variable)

def format_vector(v):
    assert len(v) == 2
    return r'\ColVec{%s}{%s}' % (latex(v[0]), latex(v[1]))

def format_label(label):
    try:
        x, y, z, eps_x, eps_y, eps_z = label
        return format_face_triple_extra_fancy(veps_to_face_dict[label].minimal_triple, label) + [format_vector((x, y))]
    except Exception as e:
        logging.warn("format_label: {}".format(e))
        return [r'\multicolumn{4}{l}{%s}' % latex(label)]

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
    ###
    head += [r'  ' + ' & '.join([r'\multicolumn{4}{c}{Eqn.~$\Delta\bar\pi_F(x, y)=0$, $F=F(I, J, K)$}',
                                 r'\multicolumn{%s}{c}{Coefficients}' % (num_columns-4)]) + r'\\']
    head += [r'  \cmidrule{1-4}\cmidrule{5-%s}' % num_columns]
    ###
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
    caption = r'Homogeneous linear system'
    extra_caption = None
    label = 'tab:hildebrand_discont_3_slope_1_finite_system'
    s += [begin_longtable(columns = ['I', 'J', 'K', '(x, y)'] + [ format_variable(v) for v in h._facet_symbolic.basis ],
                          caption=caption, extra_caption=extra_caption, label=label)]

    for label, row in zip(h._facet_used_vertices, h._facet_equation_matrix):
        s += [r'  \relax' + ' & '.join(format_label(label) + [ format_coefficient(coeff) for coeff in row ]) + r'\\']

    s += [end_longtable()]
    return '\n'.join(s)

with open('/Users/mkoeppe/w/papers/basu-hildebrand-koeppe-papers/algo-paper/tab_hildebrand_discont_3_slope_1_finite_system.tex', 'w') as f:
    f.write('%% Automatically generated\n')
    f.write(tabulate_finite_system(h))
