import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *

if 'h' not in globals():
    h = kzh_minimal_has_only_crazy_perturbation_1(parametric=True, extra_names=['x', 'y'], extra_values=[0, 0])

hh = kzh_minimal_has_only_crazy_perturbation_1()

for fn in (hh, ):
    fn._faces_used = []      # for instrumented code
    fn._vertices_used = []   # for instrumented code
    try:
        facet_test(fn, known_extreme=True)
    except NotImplementedError:
        pass

logging.warn("Checking and converting the faces used in the published proof")
## We go to some length here to be able to reuse part of the proof published in Equi VI (Appendix).
version="v0.9-274-g1342a5da"    # Version that reproduces published proof in Equi VI.
faces_used_published = load("/Users/mkoeppe/w/papers/basu-hildebrand-koeppe-papers/algo-paper/auto_proof_crazy_example-{}_faces_used.sobj".format(version))
p = hh(0).parent()
faces_used_published = set( Face([ [ p(x) for x in I ] for I in T ]) for T in faces_used_published )
# Unfortunately, our code cannot reproduce the same published faces for the covering test -- algorithm changed apparently!
###  assert faces_used_published == set(hh._faces_used)   # Not true
# But we make sure that the published faces are all additive-sans-limits etc.
assert faces_used_published.issubset(set(generate_additive_faces_sans_limits(hh)))
covered_components_published = generate_covered_components_strategically(hh, additive_faces=faces_used_published)
assert len(covered_components_published) == 2
uncovered_intervals = uncovered_intervals_from_covered_components(covered_components_published)
assert uncovered_intervals == [open_interval(hh.ucl, hh.ucr), open_interval(hh._f - hh.ucr, hh._f - hh.ucl)]
# Convert to faces of h.
faces_used = {}
from cutgeneratingfunctionology.spam.real_set import RealSet   # param-safe version of RealSet
all_triples = list(generate_triples_with_projections_intersecting(h, RealSet([0, 1])))
for T in all_triples:
    hF = Face(T)
    if hF in faces_used_published:
        faces_used[hF] = T
        faces_used_published.remove(hF)
assert not faces_used_published

h._faces_used = []
h._vertices_used = []   # for instrumented code
h._facet_covered_components = generate_covered_components_strategically(fn, additive_faces=faces_used)

logging.warn("Facet test for the parametric function")
try:
    facet_test(h, known_extreme=True)
except NotImplementedError:
    pass
vertices_used = set(h._vertices_used)

faces_of_vertices_used = set()
for


if 'all_special_faces' not in globals():
    all_special_faces = list(generate_triples_with_projections_intersecting(h, h.special_intervals, break_symmetry=True))

K = h.end_points()[1].parent()
x, y = K.gens()[-2:]

if 'faces' not in globals():
    faces = { Face(triple): triple for triple in all_special_faces }

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

def format_vertices(F, show_used=False):
    def is_used(v):
        return any((v[0], v[1], v[0]+v[1], xeps, yeps, zeps) in vertices_used
                   for xeps, yeps, zeps in generate_containing_eps_triple((v[0], v[1]), F.minimal_triple))
    def format_vertex(v):
        s = format_vector(v)
        if show_used and is_used(v):
            s = r'\llap{$\rightarrow$}' + s
        return s
    return [format_vertex(v) for v in F.vertices]

def format_interval(interval):
    l = latex(realset_from_interval(interval))
    return l.replace(r'-', r'\compactop-').replace(r'+', r'\compactop+')  # reduce spacing

def format_symbolic_slack(slack):
    l = latex(SR(slack).collect(SR.var('x')).collect(SR.var('y')).subs(sqrt(1/2)==1/2*sqrt(2)))
    return l.replace(r'} - c', r'} \compactop- c').replace(r'} + c', r'} \compactop+ c')

def format_face_symbolic_slack(F, slacks):
    if not any(slacks):
        return "0"
    else:
        return format_symbolic_slack(delta_pi_of_face_symbolic(h, x, y, F).sym())

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

def format_interval_extra_fancy(nominal_I, minimal_I):
    if len(nominal_I) < 2:
        return format_interval(nominal_I)
    s = ""
    minimal_I = list(interval_to_endpoints(minimal_I))
    if nominal_I[0] < minimal_I[0]:
        s += r'\langle '
    elif h.limit(fractional(nominal_I[0]), 0) == h.limit(fractional(nominal_I[0]), 1):
        s += r'['
    else:
        s += r'('
    s += latex(nominal_I[0]).replace(r'-', r'\compactop-').replace(r'+', r'\compactop+')  # reduce spacing
    s += ", "
    s += latex(nominal_I[1]).replace(r'-', r'\compactop-').replace(r'+', r'\compactop+')  # reduce spacing
    if minimal_I[1] < nominal_I[1]:
        s += r'\rangle '
    elif h.limit(fractional(nominal_I[1]), 0) == h.limit(fractional(nominal_I[1]), -1):
        s += r']'
    else:
        s += r')'
    if not h.special_intervals.is_disjoint_from(realset_from_interval(interval_mod_1(minimal_I))):
        s = specialinterval(s)
    return s

def format_face_triple_extra_fancy(triple):
    face = Face(triple)
    return [ format_interval_extra_fancy(nominal_I, minimal_I)
             for nominal_I, minimal_I in zip(triple, face.minimal_triple) ]

latex.add_package_to_preamble_if_available("booktabs")


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
    #  of the piecewise linear function $\pi$ = \sage{kzh\_minimal\_has\_only\_crazy\_perturbation\_1}()
    if label is not None:
        s += [r'\label{%s}\\' % label]

    head  = [r'  \toprule']
    ## head += ['  ' + ' & '.join([r'\multicolumn{3}{C}{F = F(I, J, K)}']) + r'\\']
    ## head += [r'  \cmidrule{1-3}']
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

def tabulate_additive_faces(faces, dimension=None, **longtable_kwds):
    faces = sorted(faces.items(), key=lambda fi: fi[0])
    max_vertices = 6
    s = []
    s += [begin_longtable(columns=['I', 'J', 'K',
                                   #r'\Delta\pi_{F}(x,y), \ (x,y) \in F = F(I, J, K)',
                                   r'\text{slope}',
                                   r'\multicolumn{%s}{C}{\verts(F)}' % max_vertices
                                   ],
                          num_columns = 4 + max_vertices,
                          **longtable_kwds)]
    for F, triple in faces:
        if dimension is None or F.dimension() == dimension:
            s += ['{}  ' + ' & '.join(#format_face_triple(triple)
                                          format_face_triple_extra_fancy(triple)
                                    #+ [format_symbolic_slack(delta_pi_of_face_symbolic(h, x, y, F).sym())]
                                    + format_component_slope(F)
                                    + format_vertices(F)
                                    )
                  + r'\\']
    s += [end_longtable()]
    return '\n'.join(s)

def tabulate_delta_pi(faces, dimension=None, **longtable_kwds):
    faces = sorted(faces.items(), key=lambda fi: fi[0])
    s = []
    if dimension is None or dimension == 2:
        # eps describing 2d cones of a vertex
        eps_list = [ eps for eps in nonzero_eps if all(e != 0 for e in eps) ]
        max_vertices = len(eps_list)
    elif dimension == 1:
        max_vertices = 2
    elif dimension == 0:
        max_vertices = 1
    else:
        raise ValueError("bad dimension")

    s += [begin_longtable(columns=['I', 'J', 'K',
                                   r'n_F',
                                   r'\Delta\pi_{F}(x,y), \ (x,y) \in F = F(I, J, K)',
                                   # [r'\Delta\pi_F(v_F^\bgroup{}{}{}\egroup)'.format(*[print_sign(e) for e in eps]) for eps in eps_list]    ### attributing vertices to their basis....... would need more work
                                   r'\multicolumn{%s}{C}{\Delta\pi_F(u,v), \ (u,v)\in\verts(F)}' % max_vertices],
                          num_columns = 5 + max_vertices,
                          **longtable_kwds)]

    for face, triple in faces:

        F = Face(triple)
        if dimension is None or F.dimension() == dimension:
            slacks = sorted(delta_pi_of_face(h, v[0], v[1], F) for v in F.vertices)
            # Put a {} in front because the square brackets coming from formatting the intervals
            # would otherwise be taken as a latex optional argument for the preceding command!
            s += ['{}  ' + ' & '.join(#format_face_triple(triple)
                                      format_face_triple_extra_fancy(triple)
                                    + [ latex(number_of_projections_intersecting(F, h.special_intervals)) ]
                                    + [ format_face_symbolic_slack(F, slacks) ]
                                    + [ format_slack(slack) for slack in slacks ])
                    + r'\\']

    s += [end_longtable()]
    return '\n'.join(s)

logging.warn("Writing tables")
with open('/Users/mkoeppe/w/papers/basu-hildebrand-koeppe-papers/algo-paper/tab_kzh_minimal_has_only_crazy_perturbation_delta_pi.tex', 'w') as f:
    f.write(r'%% Automatically generated.' + '\n')
    f.write(r'%% Preamble: \usepackage{longtable,booktabs,array}\newcolumntype{C}{>{$}c<{$}}' + '\n')
    f.write(r'\providecommand\compactop[1]{\kern0.2pt{#1}\kern0.2pt\relax}' + '\n')
    #f.write(r'\providecommand\specialinterval[1]{{#1}\rlap{*}}' + '\n')
    f.write(r'\providecommand\specialinterval[1]{\hphantom*{#1}\text{*}}' + '\n')
    f.write(r'\providecommand\tightslack[1]{\llap{$\triangleright$\,}{#1}}' + '\n')
    for dimension in (2, 1):   # 2 comes first (directly covering)
        caption = r'%s-dimensional faces $F$ with additivity on $%s$ for proving piecewise linearity outside of the special intervals' % ("One" if dimension == 1 else "Two",
                                                                                                                                          r'\relint(F)' if dimension == 1 else r'\intr(F)')
        label = 'tab:kzh_minimal_has_only_crazy_perturbation_1_faces_used_dim_%s' % dimension
        f.write(tabulate_additive_faces(faces_used, dimension=dimension,
                                        caption=caption, label=label))
    # Used vertices

    for dimension in (1, 2):      # 1 comes first because used earlier. 0 is empty.
        caption = r'Subadditivity slacks $\Delta\pi_F$ for $\dim F=%s$ and $n_F>0$' % dimension
        extra_caption = r'. An asterisk marks the special intervals.'
        label = r'tab:kzh_minimal_has_only_crazy_perturbation_1_delta_pi_dim_%s' % dimension
        f.write(tabulate_delta_pi(faces, dimension=dimension, caption=caption, extra_caption=extra_caption,
                                  label=label))
