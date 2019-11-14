### Uses global variable h

def compactop(l):
    # reduce spacing
    return l.replace(r'-', r'\compactop-').replace(r'+', r'\compactop+')

def format_interval(interval):
    return compactop(latex(realset_from_interval(interval)))

def specialinterval(s):
    #return r'\text{\boldmath$' + s + r'$}'
    return r'\specialinterval{' + s + r'}'

def format_interval_markup_special(interval):
    s = format_interval(interval)
    if not h.special_intervals.is_disjoint_from(realset_from_interval(interval_mod_1(interval))):
        s = specialinterval(s)
    return s

def format_face_triple(triple):
    return [ format_interval_markup_special(I) for I in triple ]

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
    if hasattr(h, "special_intervals"):
        if minimal_I[0] < minimal_I[1] and not h.special_intervals.is_disjoint_from(realset_from_interval(interval_mod_1(minimal_I))):
            s = specialinterval(s)
    return s

def format_face_triple_extra_fancy(triple, vector=None):
    face = Face(triple)
    if vector is None:
        xyz = [None, None, None]
    else:
        xyz = [vector[0], vector[1], vector[0] + vector[1]]
    return [ format_interval_extra_fancy(nominal_I, minimal_I, number)
             for nominal_I, minimal_I, number in zip(triple, face.minimal_triple, xyz) ]
format_face_triple_extra_fancy.columns = ['I', 'J', 'K']
format_face_triple_extra_fancy.head = r'\multicolumn{3}{c}{Face $F = F(I, J, K)$}'

def cmidrules(num_columns_list, start=1, truncate_right=False):
    s = ''
    n = start
    for num_columns in num_columns_list[:-1]:
        s += r'\cmidrule(r){%s-%s}' % (n, n + num_columns - 1)
        n += num_columns
    num_columns = num_columns_list[-1]
    if truncate_right:
        s += r'\cmidrule(r){%s-%s}' % (n, n + num_columns - 1)
    else:
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
                l = lx = latex(x)
                # [-3] is the last one, just after 'f'.
                if x in (h.a0, h.a1, h.a2, h.end_points()[1], h.end_points()[-3]) or len(l) == 1:
                    if '_' in l:
                        if len(l) > 3:
                            l += r'\!\!\!'
                        else:
                            l += r'\!'
                    else:
                        l += r'\vphantom{{}_0}'
                    s = r'\stackrel{%s}{\bullet}' % l
                    if lx == 'u':
                        s += r'\rlap{\ \normalsize$\wr$}'
                    return s
                return r'\bullet'
            else:
                return r'-'
        else:
            return r'\bar\pi(%s)' % (latex(x), )
    else:
        return latex(variable)

def format_vector(v):
    assert len(v) == 2
    return [r'\ColVec{%s}{%s}' % (latex(v[0]), latex(v[1]))]
format_vector.columns = [ r'' ]

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
        return format_face_triple_extra_fancy(veps_to_face_dict(h)[label].minimal_triple, label) + format_vector((x, y))
    except Exception as e:
        logging.warn("format_label: {}".format(e))
        return [r'\multicolumn{4}{l}{%s}' % latex(label)]
format_label_by_faces.columns = ['I', 'J', 'K', '(u, v)']
format_label_by_faces.head = r'\multicolumn{4}{c}{Eqn.~$\Delta\bar\pi_F(u, v)=0$, $F=F(I, J, K)$}'

def format_label_by_veps(label):
    try:
        x, y, z, eps_x, eps_y, eps_z = label
        def format_eps(eps):
            return (r'{}^{\hphantom{+}}', r'{}^+', r'{}^-')[eps]
        return [ compactop(latex(t)) + format_eps(eps) for t, eps in [[x, eps_x], [y, eps_y], [z, eps_z]] ]
    except Exception as e:
        logging.warn("format_label: {}".format(e))
        return [r'\multicolumn{3}{l}{%s}' % latex(label)]
format_label_by_veps.columns = [ r'\multicolumn{1}{C}{%s}' % s for s in ['u', 'v', 'u+v'] ]
format_label_by_veps.head = r'\multicolumn{3}{c}{Eqn.~$\Delta\bar\pi_F(u, v)=0$}'


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




def tabulate_additive_faces(faces, dimension=None, show_used=False, show_slope=True, max_vertices=None, coordinate_format='C', **longtable_kwds):

    if max_vertices is None:
        if dimension == 1:
            max_vertices = 2
        else:
            max_vertices = 3
    s = []
    format_vertex = format_vertex_by_veps
    if show_used:
        vertex_heads = [r'\multicolumn{%s}{c}{selected vertex}' % len(format_vertex.columns),
                        r'\multicolumn{%s}{c}{other vertices of $F$}' % ((max_vertices - 1) * len(format_vertex.columns))]
    else:
        vertex_heads = [r'\multicolumn{%s}{c}{vertices of $F$}' % (max_vertices * len(format_vertex.columns))]
    column_head = [r'  ' + ' & '.join([format_face_triple_extra_fancy.head] + ([r''] if show_slope else []) + vertex_heads) + r'\\']
    if show_used:
        column_head += [cmidrules([len(format_face_triple_extra_fancy.columns)], truncate_right=True)
                        + cmidrules([len(format_vertex.columns), (max_vertices - 1) * len(format_vertex.columns)],
                                    start=len(format_face_triple_extra_fancy.columns) + 1 + (1 if show_slope else 0))]
    else:
        column_head += [cmidrules([len(format_face_triple_extra_fancy.columns)], truncate_right=True)
                    + cmidrules([max_vertices * len(format_vertex.columns) ], start=len(format_face_triple_extra_fancy.columns)+2)]
    columns = format_face_triple_extra_fancy.columns + ( [r'\text{slope}'] if show_slope else [] ) + max_vertices * format_vertex.columns
    column_head += ['  ' + ' & '.join(columns) + r'\\']
    format = r'*{%s}C' % ((1 if show_slope else 0) + len(format_face_triple_extra_fancy.columns)) + max_vertices * (('*{%s}{%s}' % (len(format_vertex.columns), coordinate_format)) + r'@{\qquad}')
    s += [begin_longtable(column_head=column_head,
                          format=format,
                          num_columns=len(columns),
                          **longtable_kwds)]

    for triple in sorted(min(list(triple),
                             [triple[1], triple[0], triple[2]])
                         for face, triple in faces.items()):
        F = Face(triple)
        if dimension is None or F.dimension() == dimension:
            s += ['{}  ' + ' & '.join(#format_face_triple(triple)
                                          format_face_triple_extra_fancy(triple)
                                    #+ [format_symbolic_slack(delta_pi_of_face_symbolic(h, x, y, F).sym())]
                                    + (format_component_slope(F) if show_slope else [])
                                    + format_vertices(F, format_vertex=format_vertex, first_is_used=show_used)
                                    )
                  + r'\\']
    s += [end_longtable()]
    return '\n'.join(s)


if 'faces_used_notice_dict' not in globals():
    faces_used_notice_dict = {}

def format_used_notice(face, slacks):
    if slacks[0] == h.s:
        #return r'(\ref{claim:n_F}) tight'
        return r'(tight)'
    if face.minimal_triple[2] in ([h._f], [h._f + 1]):
        return r'(symm.)'
    return faces_used_notice_dict.get(face, None) or faces_used_notice_dict.get(x_y_swapped_face(face), r'')

def tabulate_delta_pi(faces, dimension=None, show_used=False, max_vertices=None, **longtable_kwds):
    s = []
    if max_vertices is None:
        if dimension is None or dimension == 2:
            # eps describing 2d cones of a vertex
            eps_list = [ eps for eps in nonzero_eps if all(e != 0 for e in eps) ]
            max_vertices = len(eps_list)
        elif dimension == 1:
            max_vertices = 3   # yes, it's 2 but it formats more nicely with 3
        elif dimension == 0:
            max_vertices = 1
        else:
            raise ValueError("bad dimension")

    format = (r'*{%s}C' % 4) + r'{>{\tiny$}c<{$}}' + (r'*{%s}C' % max_vertices) + '@{}r'
    s += [begin_longtable(columns=['I', 'J', 'K',
                                   r'n_F',
                                   r'\multicolumn{1}{C}{\Delta\pi_{F}(x,y), \ (x,y) \in F = F(I, J, K)}',
                                   # [r'\Delta\pi_F(v_F^\bgroup{}{}{}\egroup)'.format(*[print_sign(e) for e in eps]) for eps in eps_list]    ### attributing vertices to their basis....... would need more work
                                   r'\multicolumn{%s}{L}{\Delta\pi_F(u,v), \ (u,v)\in\verts(F)}' % max_vertices,
                                   r''  # comment
                                   ],
                          num_columns=5+max_vertices + 1,
                          format=format,
                          **longtable_kwds)]

    for triple in sorted(min(list(triple),
                             [triple[1], triple[0], triple[2]])
                         for face, triple in faces.items()):
        F = Face(triple)
        if dimension is None or F.dimension() == dimension:
            slacks = sorted(delta_pi_of_face(h, v[0], v[1], F) for v in F.vertices)
            # Put a {} in front because the square brackets coming from formatting the intervals
            # would otherwise be taken as a latex optional argument for the preceding command!
            s += ['{}  ' + ' & '.join(#format_face_triple(triple)
                                      format_face_triple_extra_fancy(triple)
                                    + [ latex(number_of_projections_intersecting(F, h.special_intervals)) ]
                                    + [ format_face_symbolic_slack(F, slacks) ]
                                    + [ format_slack(slack) for slack in slacks ]
                                    + [ "" for i in range(max_vertices - len(slacks)) ]
                                    + [ format_used_notice(F, slacks) ]
                                      )
                    + r'\\']

    s += [end_longtable()]
    return '\n'.join(s)

def format_slack(slack):
    if slack == 0:
        return "0"
    fs = "{:.3f}".format(float(slack))
    if slack == h.s:
        return r'\tightslack{' + fs + r'}'
    return fs

def format_component_slope(F):
    for I in F.minimal_triple:
        if len(I) == 2:
            return [latex(h.which_function((I[0] + I[1])/2)._slope)]
    return [""]

def format_vector(v):
    assert len(v) == 2
    return r'\ColVec{%s}{%s}' % (latex(v[0]), latex(v[1]))

def format_vertex_as_vector(veps):
    return format_vector([veps[0], veps[1]])

def format_vertex_by_veps(veps):
    return format_label_by_veps(veps)
format_vertex_by_veps.columns = [ r'\multicolumn{1}{C}{%s}' % s for s in [r'u', r'v'] ] + [r'\multicolumn{1}{C@{\qquad}}{%s}' % r'u+v']

def x_y_swapped_veps(veps):
    x, y, z, xeps, yeps, zeps = veps
    return y, x, z, yeps, xeps, zeps

def format_vertices(F, format_vertex=format_vertex_as_vector, first_is_used=False):
    vertices_used = set(itertools.chain(h._vertices_used,
                                        (x_y_swapped_veps(veps) for veps in h._vertices_used)))
    def format_it(v):
        is_used = False
        for xeps, yeps, zeps in generate_containing_eps_triple(v, F.minimal_triple):
            veps = (v[0], v[1], v[0]+v[1], xeps, yeps, zeps)
            if not first_is_used:
                break
            if veps in vertices_used:
                is_used = True
                break
        return (0 if is_used else 1), (v[0], xeps, v[1], yeps, v[0]+v[1], zeps), format_vertex(veps)
    formatted_list = sorted( format_it(v) for v in F.vertices )    # lex key
    if first_is_used:
        assert formatted_list[0][0] == 0
        if len(formatted_list) > 1:
            assert formatted_list[1][0] != 0
    return list(itertools.chain.from_iterable(f for _, __, f in formatted_list))

def format_symbolic_slack(slack):
    l = latex(SR(slack).collect(SR.var('x')).collect(SR.var('y')).subs(sqrt(1/2)==1/2*sqrt(2)))
    return l.replace(r'} - c', r'} \compactop- c').replace(r'} + c', r'} \compactop+ c')

def format_face_symbolic_slack(F, slacks):
    if not any(slacks):
        return "0"
    else:
        return format_symbolic_slack(delta_pi_of_face_symbolic(h, x, y, F).sym())
