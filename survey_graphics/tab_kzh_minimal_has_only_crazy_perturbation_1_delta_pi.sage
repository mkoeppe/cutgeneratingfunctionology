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

## def smallest_function_interval_containing(interval, function):
##     from cutgeneratingfunctionology.spam.real_set import RealSet   # param-safe version of RealSet
##     if len(interval) == 1:
##         if interval[0] in function.end_points():
##             return interval
##     rs = realset_from_interval(interval)
##     for I in function.intervals():
##         if I[0] < I[1] and rs.is_included_in(RealSet.closed(I[0], I[1])):
##             return I
##     raise ValueError("interval not contained in a breakpoint interval of function")

## def breakpoint_interval_triple(F, function):
##     """
##     The opposite of ``minimal_triple``.
##     """
##     triple = [ smallest_function_interval_containing(I, function) for I in F.minimal_triple ]
##     if F != Face(triple):
##         raise ValueError("not a face defined by breakpoint intervals")
##     return triple

def format_slack(slack):
    if slack == h.s:
        return r"\llap\bgroup$s = \bgroup\egroup$\egroup{:.4f}".format(float(slack))
    elif slack == 0:
        return "0"
    else:
        return "{:.4f}".format(float(slack))

def format_interval(interval):
    l = latex(realset_from_interval(interval))
    return l.replace(r'-', r'\,{-}\,').replace(r'+', r'\,{+}\,')  # reduce spacing

def format_interval_markup_special(interval):
    if not h.special_intervals.is_disjoint_from(realset_from_interval(interval_mod_1(interval))):
        #return r'\text{\boldmath$' + format_interval(interval) + r'$}'
        return format_interval(interval) + r'\rlap{*}'
    else:
        return format_interval(interval)

latex.add_package_to_preamble_if_available("booktabs")

def tabulate_delta_pi(faces, dimension):
    s = []
    if dimension == 2:
        # eps describing 2d cones of a vertex
        eps_list = [ eps for eps in nonzero_eps if all(e != 0 for e in eps) ]
    elif dimension == 1:
        eps_list = [ eps for eps in nonzero_eps if not all(e != 0 for e in eps) ]
    elif dimension == 0:
        eps_list = [[0, 0, 0]]
    else:
        raise ValueError("bad dimension")

    num_columns = 4 + len(eps_list)

    s += [r'\def\arraystretch{1.17}']
    s += [r'\begin{longtable}{*{%s}C}' % num_columns]
    s += [r'\caption{Subadditivity slacks $\Delta\pi_F(x,y)$, $x, y\in\verts(F)$ for faces $F$ with $\dim F=%s$ and $n_F>0$}' % dimension]
    #  of the piecewise linear function $\pi$ = \sage{kzh\_minimal\_has\_only\_crazy\_perturbation\_1}()
    s += [r'\label{tab:kzh_minimal_has_only_crazy_perturbation_1_delta_pi_dim_%s}\\' % dimension]

    ##  \toprule
    ##   I & J & K & \Delta\pi_{F}, \ F = F(I, J, K)} &
    ##   \multicolumn{6}{c}{$\Delta\pi_F(x,y),\ (x,y)\in\verts(F)$}
    ## \\
    ##  \midrule
    ## \endhead
    ##   \bottomrule
    ## \endfoot

    s += [r'  \toprule']
    ## s += ['  ' + ' & '.join([r'\multicolumn{3}{C}{F = F(I, J, K)}']) + r'\\']
    ## s += [r'  \cmidrule{1-3}']
    s += ['  ' + ' & '.join(['I', 'J', 'K',
                             r'\Delta\pi_{F}, \ F = F(I, J, K)',
                             # [r'\Delta\pi_F(v_F^\bgroup{}{}{}\egroup)'.format(*[print_sign(e) for e in eps]) for eps in eps_list]    ### attributing vertices to their basis....... would need more work
                             r'\multicolumn{%s}{C}{\Delta\pi_F(x,y),\ (x,y)\in\verts(F)}' % len(eps_list)
                            ])
          + r'\\']
    s += [r'  \midrule']
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
                                    + [latex(delta_pi_of_face(h, x, y, F).sym())]
                                    + [ format_slack(slack) for slack in slacks ])
                    + r'\\']

    s += [r'\end{longtable}']
    return '\n'.join(s)

with open('/Users/mkoeppe/w/papers/basu-hildebrand-koeppe-papers/algo-paper/tab_kzh_minimal_has_only_crazy_perturbation_delta_pi.tex', 'w') as f:
    f.write(r'%% Automatically generated.' + '\n')
    f.write(r'%% Preamble: \usepackage{longtable,booktabs,array}\newcolumntype{C}{>{$}c<{$}}' + '\n')
    f.write(tabulate_delta_pi(faces, dimension=2))
    f.write(tabulate_delta_pi(faces, dimension=1))
    #f.write(tabulate_delta_pi(faces, dimension=0))
