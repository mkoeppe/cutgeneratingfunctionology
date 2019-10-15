import cutgeneratingfunctionology.igp as igp
from cutgeneratingfunctionology.dff import *

try:
    destdir = dff_output_dir
except Exception:
    destdir = "survey_graphics/dff_graphics/"


#ftype = ".png"
ftype = ".pdf"

## Set up style - no "legend"

paper_plot_kwds = { 'dpi': 200, 'fontsize': 10 };

def plot_kwds_hook_paper(kwds):
    plot_kwds_hook_no_legend(kwds)
    kwds.update(paper_plot_kwds)

igp.plot_kwds_hook = plot_kwds_hook_paper

## Where are the following used?

f=phi_simple(3/2)
g=plot(f)
g.save(destdir + "f_simple1" + ftype)

f=phi_simple(7/3)
g=plot(f)
g.save(destdir + "f_simple2" + ftype)

f=phi_simple(13/4)
g=plot(f)
g.save(destdir + "f_simple3" + ftype)

f=phi_bj_1(3/2)
g=plot(f)
g.save(destdir + "f_bj1" + ftype)

f=phi_bj_1(7/3)
g=plot(f)
g.save(destdir + "f_bj2" + ftype)

f=phi_bj_1(13/4)
g=plot(f)
g.save(destdir + "f_bj3" + ftype)

f=phi_ccm_1(3/2)
g=plot(f)
g.save(destdir + "f_ccm1" + ftype)

f=phi_ccm_1(7/3)
g=plot(f)
g.save(destdir + "f_ccm2" + ftype)

f=phi_ccm_1(13/4)
g=plot(f)
g.save(destdir + "f_ccm3" + ftype)

f=phi_fs_1(2)
g=plot(f)
g.save(destdir + "f_fs1" + ftype)

f=phi_fs_1(5)
g=plot(f)
g.save(destdir + "f_fs2" + ftype)

f=phi_fs_1(11)
g=plot(f)
g.save(destdir + "f_fs3" + ftype)

f=phi_vb_2(2)
g=plot(f)
g.save(destdir + "f_vb1" + ftype)

f=phi_vb_2(5)
g=plot(f)
g.save(destdir + "f_vb2" + ftype)

f=phi_vb_2(11)
g=plot(f)
g.save(destdir + "f_vb3" + ftype)

f=phi_ll_1(3/2,3)
g=plot(f)
g.save(destdir + "f_ll11" + ftype)

f=phi_ll_1(7/3,3)
g=plot(f)
g.save(destdir + "f_ll12" + ftype)

f=phi_ll_1(13/4,3)
g=plot(f)
g.save(destdir + "f_ll13" + ftype)

f=phi_ll_2(3/2,3)
g=plot(f)
g.save(destdir + "f_ll21" + ftype)

f=phi_ll_2(7/3,3)
g=plot(f)
g.save(destdir + "f_ll22" + ftype)

f=phi_ll_2(13/4,3)
g=plot(f)
g.save(destdir + "f_ll23" + ftype)

f=phi_dg_1(3/2,3)
g=plot(f)
g.save(destdir + "f_dg1" + ftype)

f=phi_dg_1(7/3,3)
g=plot(f)
g.save(destdir + "f_dg2" + ftype)

f=phi_dg_1(13/4,3)
g=plot(f)
g.save(destdir + "f_dg3" + ftype)


f=phi_bj_1(7/3)
g=plot_2d_diagram_dff(f)
g.save(destdir + "maximal_function_example" + ftype)

f=gmic(1/5)
g=plot_2d_diagram_dff(f)
g.save(destdir + "nonmaximal_function_example" + ftype)

f=phi_simple(7/3)
g=plot_2d_diagram_with_cones_dff(f)
g.save(destdir + "discontinuous_function_example" + ftype)

f=w_2slope_3covered()
g=plot_2d_diagram_dff_no_lable(f, colorful=True)
g.save(destdir + "continuous_2slope_3covered" + ftype)

# Figure 1 in dff-paper-long-version
igp.show_plots_figsize = 5
paper_plot_kwds['fontsize'] = 15
f=phi_bj_1(5/2)
g=plot_2d_diagram_dff_no_lable(f,show_projections=False)
g.save(destdir + "maximality_test_example" + ftype, figsize=igp.show_plots_figsize)


# Missing scripts for graphics used in dff-paper-long-version
## Figure 2: Graphs of $\phi_{BJ,1}$ \cite[Example 3.1]{alves-clautiaux-valerio-rietz-2016:dual-feasible-book} for $C=3/2$ (left) and $C=7/3$ (right). }
## 2slope_dff_example1.png
def fig2_ticks_keywords(function, y_ticks_for_breakpoints=False):
    xticks = yticks = [-2, -1, 0, 1, 2]
    xtick_formatter = [ latex_tick_formatter_x(x) for x in xticks ]
    ytick_formatter = xtick_formatter
    return { 'ticks': [xticks, yticks],
             'tick_formatter': [xtick_formatter, ytick_formatter],
             'gridlines': True }
igp.ticks_keywords = fig2_ticks_keywords
phi = phi_bj_1_gdff(3/2)
g = plot_with_colored_slopes(phi, thickness=2, figsize=(4, 4))
g.save(destdir + "2slope_dff_example1" + ftype, xmin=-2, xmax=2, ymin=-2, ymax=2, aspect_ratio=1)
## 2slope_dff_example2.png
phi = phi_bj_1_gdff(7/3, periods=5)
g = plot_with_colored_slopes(phi, thickness=2, figsize=(4, 4))
g.save(destdir + "2slope_dff_example2" + ftype, xmin=-2, xmax=2, ymin=-2, ymax=2, aspect_ratio=1)
igp.ticks_keywords = default_igp_ticks_keywords

## Figure 3
## w_2slope_3covered_nonextreme-2d_diagram.png
phi = w_2slope_3covered_nonextreme()
paper_plot_kwds['fontsize'] = 24
g = plot_2d_diagram_dff_no_lable(phi, colorful=True)
g.save(destdir + "w_2slope_3covered_nonextreme-2d_diagram" + ftype, figsize=8)

## Figure 4
## w_2slope_3covered_nonextreme-perturbations.png
gen = generate_perturbations_finite_dimensional_continuous_dff(phi)
perturbation = next(gen)
igp.find_epsilon_interval = find_epsilon_interval_continuous_dff   # just a hack so that plot_perturbation_diagram works for us
igp.check_perturbation_plot_three_perturbations = False
igp.check_perturbation_plot_rescaled_perturbation = True
paper_plot_kwds['fontsize'] = 12
g = plot_perturbation_diagram(phi, perturbation)
g.save(destdir + "w_2slope_3covered_nonextreme-perturbations" + ftype, aspect_ratio=0.3, figsize=(8, 2.5), show_legend=False, **paper_plot_kwds)

## Figure 5
##phi_s_delta.png
## The graph of $\phi_{s,\delta}$ for $s=\frac{1}{5}$ and $\delta=2$  <--- Caption was wrong
K = ParametricRealField([QQ('1/5'), 2], names=['delta', 's'])
delta, s = K.gens()
phi = phi_s_delta(delta, s)
paper_plot_kwds['fontsize'] = 12
g = plot_with_colored_slopes(phi, thickness=2, figsize=(8, 2.5))
g.save(destdir + "phi_s=2_delta=1_5" + ftype, xmin=-1, xmax=2, ymin=-1.5, ymax=2.5, aspect_ratio=0.3)


## Not used, no script, but in Dropbox: 2d_diagram_phi_s=2_delta=1_5.png
