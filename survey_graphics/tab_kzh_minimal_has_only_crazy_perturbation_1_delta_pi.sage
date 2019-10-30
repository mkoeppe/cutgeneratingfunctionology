import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *


if 'h' not in globals():
    h = kzh_minimal_has_only_crazy_perturbation_1(parametric=True)
if 'all_special_faces' not in globals():
    all_special_faces = [ F for F in generate_all_faces(h)
                          if number_of_projections_intersecting(F, h.special_intervals) != 0 ] ## LONG


faces = [ F for F in all_special_faces if F.is_2D() ]

# eps describing 2d cones of a vertex
eps_list = [ eps for eps in nonzero_eps if all(e != 0 for e in eps) ]

latex.add_package_to_preamble_if_available("booktabs")
s = []
num_columns = 4 + len(eps_list)
s += [r'\begin{array}{*%sc}' % num_columns]
s += [r'  \toprule']
s += ['  ' + ' & '.join([r'\multicolumn{3}{c}{F = F(I, J, K)}'])]
s += [r'  \cmidrule']
s += ['  ' + ' & '.join(['I', 'J', 'K',
                         r'\Delta\pi_F']
                         + [r'\Delta\pi_F(v_F^\bgroup{}{}{}\egroup)'.format(*[print_sign(e) for e in eps])
                            for eps in eps_list])
      + r'\\']
s += [r'  \midrule']

var('x, y')

def format_slack(slack):
    if slack == h.s:
        return "s = {:.4}".format(slack)
    elif slack == 0:
        return "0"
    else:
        return "{:.4}".format(slack)

for F in faces:

    s += ['  ' + ' & '.join([ latex(I) for I in F.minimal_triple ]
                            + [latex(delta_pi_of_face(h, x, y, F))]
                            + [ latex(delta_pi_of_face(h, v[0], v[1], F)) for v in F.vertices ]) ] # wrong order

s += [r'  \bottomrule']
s += [r'\end{array}']
'\n'.join(s)

with open('/Users/mkoeppe/w/papers/basu-hildebrand-koeppe-papers/algo-paper/tab_kzh_minimal_has_only_crazy_perturbation_delta_pi.tex', 'w') as f:
    f.write('%% Automatically generated.\n')
    f.write(s)
