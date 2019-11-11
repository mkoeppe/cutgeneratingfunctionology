### Uses global variable h

def compactop(l):
    # reduce spacing
    return l.replace(r'-', r'\compactop-').replace(r'+', r'\compactop+')

def format_interval(interval):
    return compactop(latex(realset_from_interval(interval)))

def boldmath(s):
    return r'\text{\boldmath$' + s + '$}'

def format_interval_extra_fancy(nominal_I, minimal_I, number=None):
    if len(nominal_I) < 2:
        s = format_interval(nominal_I)
        if nominal_I[0] == number:
            s = boldmath(s)
        return s
    minimal_I = list(interval_to_endpoints(minimal_I))
    left = ""
    if nominal_I[0] < minimal_I[0]:
        left += r'\langle '
    elif h.limit(fractional(nominal_I[0]), 0) == h.limit(fractional(nominal_I[0]), 1):
        left += r'['
    else:
        left += r'('
    left += compactop(latex(nominal_I[0]))
    if nominal_I[0] == number:
        left = boldmath(left)
    right = compactop(latex(nominal_I[1]))
    if minimal_I[1] < nominal_I[1]:
        right += r'\rangle '
    elif h.limit(fractional(nominal_I[1]), 0) == h.limit(fractional(nominal_I[1]), -1):
        right += r']'
    else:
        right += r')'
    if nominal_I[1] == number:
        right = boldmath(right)
    s = left + r', ' + right
    ## if not h.special_intervals.is_disjoint_from(realset_from_interval(interval_mod_1(minimal_I))):
    ##     s = specialinterval(s)
    return s

def format_face_triple_extra_fancy(triple, vector=None):
    face = Face(triple)
    if vector is None:
        xyz = [None, None, None]
    else:
        xyz = [vector[0], vector[1], vector[0] + vector[1]]
    return [ format_interval_extra_fancy(nominal_I, minimal_I, number)
             for nominal_I, minimal_I, number in zip(triple, face.minimal_triple, xyz) ]

def cmidrules(num_columns_list):
    s = ''
    n = 1
    for num_columns in num_columns_list[:-1]:
        s += r'\cmidrule(r){%s-%s}' % (n, n + num_columns - 1)
        n += num_columns
    num_columns = num_columns_list[-1]
    s += r'\cmidrule{%s-%s}' % (n, n + num_columns - 1)
    return s

def begin_longtable(columns=None, num_columns=None, column_head=None, format=None, caption=None,
                    extra_caption=None, label=None):
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
    if column_head is None:
        column_head = ['  ' + ' & '.join(columns) + r'\\']
    head += column_head
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

############## tabulate_finite_system ###############

def format_variable(variable, compact_head=False):
    if variable[0] == 'slope of component':
        slope_index = variable[1]
        if hasattr(h, '_slope_names'):
            l = r'\bar ' + h._slope_names[slope_index]
        else:
            l = r'\bar c_{%s}}' % (1 + slope_index, )
        return r'\multicolumn{1}{C}{%s}' % l
    elif variable[0] == 'function value at':
        x = variable[1]
        if compact_head:
            if x in h.end_points():
                l = None
                if x == h.a0:
                    l = r'a_0\!'
                elif x == h.a1:
                    l = r'a_1\!'
                elif x == h.a2:
                    l = r'a_2\!'
                else:
                    lx = latex(x)
                    if len(lx) == 1:
                        l = lx
                if l:
                    return r'\stackrel{%s\vphantom{{}_0}}{\bullet}' % l
                return r'\bullet'
            else:
                return r'-'
        else:
            return r'\bar\pi(%s)' % (latex(x), )
    else:
        return latex(variable)

def format_vector(v):
    assert len(v) == 2
    return r'\ColVec{%s}{%s}' % (latex(v[0]), latex(v[1]))

def veps_to_face_dict(h):
    if not hasattr(h, '_veps_to_face_dict'):
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
        h._veps_to_face_dict = veps_to_face_dict
    return h._veps_to_face_dict

def format_label_by_faces(label):
    try:
        x, y, z, eps_x, eps_y, eps_z = label
        return format_face_triple_extra_fancy(veps_to_face_dict(h)[label].minimal_triple, label) + [format_vector((x, y))]
    except Exception as e:
        logging.warn("format_label: {}".format(e))
        return [r'\multicolumn{4}{l}{%s}' % latex(label)]
format_label_by_faces.columns = ['I', 'J', 'K', '(x, y)']
format_label_by_faces.head = r'\multicolumn{4}{c}{Eqn.~$\Delta\bar\pi_F(x, y)=0$, $F=F(I, J, K)$}'

def format_label_by_veps(label):
    try:
        x, y, z, eps_x, eps_y, eps_z = label
        def format_eps(eps):
            return (r'{}^{\hphantom{+}}', r'{}^+', r'{}^-')[eps]
        return [ compactop(latex(t)) + format_eps(eps) for t, eps in [[x, eps_x], [y, eps_y], [z, eps_z]] ]
    except Exception as e:
        logging.warn("format_label: {}".format(e))
        return [r'\multicolumn{3}{l}{%s}' % latex(label)]
format_label_by_veps.columns = [ r'\multicolumn{1}{C}{%s}' % s for s in ['x', 'y', 'x+y'] ]
format_label_by_veps.head = r'\multicolumn{3}{c}{Eqn.~$\Delta\bar\pi_F(x, y)=0$}'


def format_coefficient(coeff):
    if coeff == 0:
        return r' \cdot '
    elif coeff == 1:
        return r' + '
    elif coeff == -1:
        return r' - '
    else:
        l = latex(coeff)
        if r'\frac' in l:
            l = r'\frac12(' + latex(2*coeff) + r')'
        return compactop(l)

def tabulate_finite_system(h, label=None, format=None, caption=None, extra_caption=None):
    s = []
    #format_label = format_label_by_faces
    format_label = format_label_by_veps
    basis = h._facet_symbolic.basis
    num_value_columns = len([ v for v in basis if v[0] == 'function value at' ])
    num_slope_columns = len([ v for v in basis if v[0] == 'slope of component' ])
    compact_head = num_value_columns >= 10
    columns = format_label.columns + [ format_variable(v, compact_head=compact_head) for v in basis ]
    num_columns = len(columns)
    column_head = [r'  ' + ' & '.join([format_label.head,
                                       r'\multicolumn{%s}{c}{Coefficients of $\bar\pi(x_i)$ and $\bar\pi(\frac{x_i+x_{i+1}}2)$}' % num_value_columns,
                                       r'\multicolumn{%s}{c}{Coefficients of slopes}' % num_slope_columns]) + r'\\']
    column_head += [cmidrules([len(format_label.columns), num_value_columns, num_slope_columns])]
    column_head += ['  ' + ' & '.join(columns) + r'\\']
    s += [begin_longtable(column_head=column_head, num_columns=num_columns,
                          format=format,
                          caption=caption, extra_caption=extra_caption, label=label)]

    v_M = zip(h._facet_used_vertices, h._facet_equation_matrix)
    # quasi-echelonize
    v_M = sorted(v_M, key=lambda label_row: label_row[1].support())

    for label, row in v_M:
        if abs(min(row)) > max(row):
            row = -row
        s += [r'  \relax ' + ' & '.join(format_label(label) + [ format_coefficient(coeff) for coeff in row ]) + r'\\']

    s += [end_longtable()]
    return '\n'.join(s)
