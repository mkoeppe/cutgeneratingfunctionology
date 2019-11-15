import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *

load("survey_graphics/tab_functions.sage")

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
faces_of_vertices_used = {}
for T in all_triples:
    hF = Face(T)
    for v in hF.vertices:
        for eps in generate_containing_eps_triple(v, hF.minimal_triple, with_limits=False):
            veps = (v[0], v[1], v[0]+v[1], eps[0], eps[1], eps[2])
            if veps in vertices_used:
                assert is_additive_face_sans_limits(h, hF)
                faces_of_vertices_used[hF] = T
                vertices_used.remove(veps)
assert not vertices_used

if 'all_special_faces' not in globals():
    all_special_faces = list(generate_triples_with_projections_intersecting(h, h.special_intervals, break_symmetry=True, halfopen_cover_only=True))

K = h.end_points()[1].parent()
x, y = K.gens()[-2:]

if 'faces' not in globals():
    faces = { Face(triple): triple for triple in all_special_faces }


latex.add_package_to_preamble_if_available("booktabs")

for a in h.a0, h.a1, h.a2:
    faces_used_notice_dict[Face(([a], [h.ucl, h.ucr], [h._f - h.ucr, h._f - h.ucl]))] = r'(\ref{li:bar-pi-constant-on-cosets-on-special})'
faces_used_notice_dict[Face(([h.ucl, h.ucr], [h.ucl, h.ucr], [h.ucl + h.ucr]))] = r'(\ref{li:bar-pi-additivities-on-special-at-z=l+u--pre})'

def write_header(f):
    f.write(r'%% Automatically generated.' + '\n')
    f.write(r'%% Preamble: \usepackage{longtable,booktabs,array}\newcolumntype{C}{>{$}c<{$}}' + '\n')
    f.write(r'\providecommand\compactop[1]{\kern0.2pt{#1}\kern0.2pt\relax}' + '\n')
    #f.write(r'\providecommand\specialinterval[1]{{#1}\rlap{*}}' + '\n')
    f.write(r'\providecommand\specialinterval[1]{\hphantom*{#1}\text{*}}' + '\n')
    f.write(r'\providecommand\tightslack[1]{\llap{$\triangleright$\,}{#1}}' + '\n')

logging.warn("Writing tables")

def write_tables():
    intervals_explainer = r'All intervals $I$, $J$, $K$ are closed and elements of the complex~$\P$; notation $\langle a, b\rangle$: endpoints are not reached by the projection of the face; $(a, b)$: function $\pi$ is one-sided discontinuous at the endpoints from within the interval; $[a, b]$: function $\pi$ is one-sided continuous at the endpoints from within the interval.'

    for dimension in (2, 1):
        with open('/Users/mkoeppe/w/papers/basu-hildebrand-koeppe-papers/algo-paper/tab_kzh_minimal_has_only_crazy_perturbation_delta_pi_covering_dim{}.tex'.format(dimension), 'w') as f:
            write_header(f)
            caption = r'%s-dimensional faces $F$ with additivity on $%s$ for proving piecewise linearity outside of the special intervals' % ("One" if dimension == 1 else "Two",
                                                                                                                                              r'\relint(F)' if dimension == 1 else r'\intr(F)')
            extra_caption = r'. ' + intervals_explainer
            label = 'tab:kzh_minimal_has_only_crazy_perturbation_1_faces_used_dim_%s' % dimension
            f.write(tabulate_additive_faces(faces_used, dimension=dimension,
                                            caption=caption, extra_caption=extra_caption, label=label))


    with open('/Users/mkoeppe/w/papers/basu-hildebrand-koeppe-papers/algo-paper/tab_kzh_minimal_has_only_crazy_perturbation_delta_pi_rank.tex', 'w') as f:
        write_header(f)
        # Faces of vertices used for rank.
        caption = r'Faces $F$ with additivity on $\relint(F)$, one vertex of each providing an equation $\Delta\bar\pi_F(u, v)=0$'
        extra_caption = r', to form a full-rank homogeneous linear system in the proof of \autoref{lemma:discontinuous_examples_2}. ' + intervals_explainer
        label = 'tab:kzh_minimal_has_only_crazy_perturbation_1_faces_of_vertices_used'
        f.write(tabulate_additive_faces(faces_of_vertices_used, show_used=True, show_slope=False, max_vertices=3,
                                        coordinate_format=r'>{\tiny$}c<{$}',
                                        caption=caption, extra_caption=extra_caption, label=label))

    with open('/Users/mkoeppe/w/papers/basu-hildebrand-koeppe-papers/algo-paper/tab_kzh_minimal_has_only_crazy_perturbation_delta_pi.tex', 'w') as f:
        write_header(f)
        #
        for dimension in (1, 2):      # 1 comes first because used earlier. 0 is empty.
            caption = r'Subadditivity slacks $\Delta\pi_F$ for $\dim F=%s$ and $n_F>0$' % dimension
            extra_caption = r'. ' + intervals_explainer + ' An asterisk marks the special intervals.'
            label = r'tab:kzh_minimal_has_only_crazy_perturbation_1_delta_pi_dim_%s' % dimension
            f.write(tabulate_delta_pi(faces, dimension=dimension, caption=caption, extra_caption=extra_caption,
                                      label=label))

write_tables()
