import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *

if 'h' not in globals():
    h = kzh_minimal_has_only_crazy_perturbation_1(parametric=True, extra_names=['x', 'y'], extra_values=[0, 0])
if 'all_special_faces' not in globals():
    all_special_faces = list(generate_triples_with_projections_intersecting(h, h.special_intervals, break_symmetry=True))

K = h.end_points()[1].parent()
x, y = K.gens()[-2:]

if 'faces' not in globals():
    faces = { Face(triple): triple for triple in all_special_faces }
    faces = sorted(faces.items(), key=lambda fi: fi[0])

def format_slack(slack):
    if slack == 0:
        return "0"
    fs = "{:.3f}".format(float(slack))
    if slack == h.s:
        return r'\tightslack{' + fs + r'}'
    return fs

def format_interval(interval):
    l = latex(realset_from_interval(interval))
    return l.replace(r'-', r'\compactop-').replace(r'+', r'\compactop+')  # reduce spacing

def format_symbolic_slack(slack):
    l = latex(SR(slack).collect(SR.var('x')).collect(SR.var('y')).subs(sqrt(1/2)==1/2*sqrt(2)))
    return l.replace(r'- c', r'\compactop- c').replace(r'+ c', r'\compactop+ c')

def format_interval_markup_special(interval):
    if not h.special_intervals.is_disjoint_from(realset_from_interval(interval_mod_1(interval))):
        #return r'\text{\boldmath$' + format_interval(interval) + r'$}'
        return r'\specialinterval{' + format_interval(interval) + r'}'
    else:
        return format_interval(interval)

latex.add_package_to_preamble_if_available("booktabs")

def tabulate_delta_pi(faces, dimension):
    s = []
    if dimension == 2:
        # eps describing 2d cones of a vertex
        eps_list = [ eps for eps in nonzero_eps if all(e != 0 for e in eps) ]
        max_vertices = len(eps_list)
    elif dimension == 1:
        max_vertices = 2
    elif dimension == 0:
        max_vertices = 1
    else:
        raise ValueError("bad dimension")

    num_columns = 5 + max_vertices

    s += [r'\def\arraystretch{1.17}']
    s += [r'\begin{longtable}{*{%s}C}' % num_columns]
    main_caption = r'Subadditivity slacks $\Delta\pi_F$ for $\dim F=%s$ and $n_F>0$' % dimension
    extra_caption = r'. An asterisk marks the special intervals.'
    s += [r'\caption{' + main_caption + extra_caption + r'}']
    #  of the piecewise linear function $\pi$ = \sage{kzh\_minimal\_has\_only\_crazy\_perturbation\_1}()
    s += [r'\label{tab:kzh_minimal_has_only_crazy_perturbation_1_delta_pi_dim_%s}\\' % dimension]

    head  = [r'  \toprule']
    ## head += ['  ' + ' & '.join([r'\multicolumn{3}{C}{F = F(I, J, K)}']) + r'\\']
    ## head += [r'  \cmidrule{1-3}']
    head += ['  ' + ' & '.join(['I', 'J', 'K',
                                r'n_F',
                                r'\Delta\pi_{F}(x,y), \ (x,y) \in F = F(I, J, K)',
                             # [r'\Delta\pi_F(v_F^\bgroup{}{}{}\egroup)'.format(*[print_sign(e) for e in eps]) for eps in eps_list]    ### attributing vertices to their basis....... would need more work
                             r'\multicolumn{%s}{C}{\Delta\pi_F(u,v), \ (u,v)\in\verts(F)}' % max_vertices
                            ])
          + r'\\']
    head += [r'  \midrule']

    s += head
    s += [r'\endfirsthead']

    s += [r'\caption{' + main_caption + r' (ctd.)}\\']
    s += head
    s += [r'\endhead']

    s += [r'  \bottomrule']
    s += [r'\endfoot']

    #faces = faces[:40]

    for face, triple in faces:

        F = Face(triple)
        if F.dimension() == dimension:
            slacks = sorted(delta_pi_of_face(h, v[0], v[1], F) for v in F.vertices)
            # Put a {} in front because the square brackets coming from formatting the intervals
            # would otherwise be taken as a latex optional argument for the preceding command!
            s += ['{}  ' + ' & '.join([ format_interval_markup_special(I) for I in triple ]
                                    + [ latex(number_of_projections_intersecting(F, h.special_intervals)) ]
                                    + [format_symbolic_slack(delta_pi_of_face_symbolic(h, x, y, F).sym())]
                                    + [ format_slack(slack) for slack in slacks ])
                    + r'\\']

    s += [r'\end{longtable}']
    return '\n'.join(s)

with open('/Users/mkoeppe/w/papers/basu-hildebrand-koeppe-papers/algo-paper/tab_kzh_minimal_has_only_crazy_perturbation_delta_pi.tex', 'w') as f:
    f.write(r'%% Automatically generated.' + '\n')
    f.write(r'%% Preamble: \usepackage{longtable,booktabs,array}\newcolumntype{C}{>{$}c<{$}}' + '\n')
    f.write(r'\providecommand\compactop[1]{\kern1pt{#1}\kern0.5pt\relax}' + '\n')
    #f.write(r'\providecommand\specialinterval[1]{{#1}\rlap{*}}' + '\n')
    f.write(r'\providecommand\specialinterval[1]{\hphantom*{#1}\text{*}}' + '\n')
    f.write(r'\providecommand\tightslack[1]{\llap{$\triangleright$\,}{#1}}' + '\n')
    #f.write(tabulate_delta_pi(faces, dimension=0))    # empty!
    f.write(tabulate_delta_pi(faces, dimension=1))   # Comes first because used earlier.
    f.write(tabulate_delta_pi(faces, dimension=2))
