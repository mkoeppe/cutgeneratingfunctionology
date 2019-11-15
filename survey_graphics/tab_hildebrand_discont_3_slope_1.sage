from cutgeneratingfunctionology.igp import *
h = hildebrand_discont_3_slope_1()
with open('/Users/mkoeppe/w/papers/basu-hildebrand-koeppe-papers/algo-paper/tab_hildebrand_discont_3_slope_1.tex', 'w') as f:
    f.write('%% Automatically generated\n')
    f.write(h._latex_(table=True, name=r'\psi'))
