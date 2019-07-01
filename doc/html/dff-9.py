from cutgeneratingfunctionology.dff import *
import cutgeneratingfunctionology.igp as igp
K = ParametricRealField([QQ('1/5'), 2], names=['delta', 's'])
delta, s = K.gens()
phi = phi_s_delta(delta, s)
igp.plot_kwds_hook = plot_kwds_hook_no_legend
g = plot_with_colored_slopes(phi, thickness=2, figsize=(8, 2.5))
sphinx_plot(g, xmin=-1, xmax=2, ymin=-1.5, ymax=2.5, aspect_ratio=0.3)