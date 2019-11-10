import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *

parametric=True

h = kzh_minimal_has_only_crazy_perturbation_1(parametric=parametric)
components = generate_covered_intervals(h)[:2]
symbolic = generate_symbolic(h, components)
M, vs = generate_additivity_equations(h, symbolic, reduce_system=True, return_vertices=True, undefined_ok=True)
#

def latex_system(M):
    return latex(M).replace(r' 0 ', r' \cdot ').replace(r' -1 ', r' - ').replace(r' 1 ', r' + ')



with open('/Users/mkoeppe/w/papers/basu-hildebrand-koeppe-papers/algo-paper/tab_kzh_minimal_has_only_crazy_perturbation_finite_system.tex', 'w') as f:
    f.write(r'%% Automatically generated.' + '\n')
    f.write(latex_system(M))


# Reproduce the computation in the Appendix of Equivariant VI.
h6 = kzh_minimal_has_only_crazy_perturbation_1(parametric=parametric)
components6 = generate_covered_intervals(h6)
igp.generate_symbolic_two_sided_discontinuous_basis_functions = ('slopes', 'jumps')   # Restore classic behavior
symbolic6 = generate_symbolic(h6, components6)
M6, vs6 = generate_additivity_equations(h6, symbolic6, reduce_system=True, return_vertices=True, undefined_ok=True)

with open('/Users/mkoeppe/w/papers/basu-hildebrand-koeppe-papers/algo-paper/tab_kzh_minimal_has_only_crazy_perturbation_finite_system_equi6.tex', 'w') as f:
    f.write(r'%% Automatically generated.' + '\n')
    f.write(latex_system(M6))
