import cutgeneratingfunctionology.igp as igp
from cutgeneratingfunctionology.igp import *

try:
    destdir = dff_output_dir
except Exception:
    destdir = "survey_graphics/dff_graphics/"


## Set up style - no "legend"

def plot_kwds_hook_paper(kwds):
    plot_kwds_hook_no_legend(kwds)
    kwds.update(paper_plot_kwds)

igp.plot_kwds_hook = plot_kwds_hook_paper

## Where are the following used?

f=phi_simple(3/2)
g=plot(f)
g.save(destdir + "f_simple1.png")

f=phi_simple(7/3)
g=plot(f)
g.save(destdir + "f_simple2.png")

f=phi_simple(13/4)
g=plot(f)
g.save(destdir + "f_simple3.png")

f=phi_bj_1(3/2)
g=plot(f)
g.save(destdir + "f_bj1.png")

f=phi_bj_1(7/3)
g=plot(f)
g.save(destdir + "f_bj2.png")

f=phi_bj_1(13/4)
g=plot(f)
g.save(destdir + "f_bj3.png")

f=phi_ccm_1(3/2)
g=plot(f)
g.save(destdir + "f_ccm1.png")

f=phi_ccm_1(7/3)
g=plot(f)
g.save(destdir + "f_ccm2.png")

f=phi_ccm_1(13/4)
g=plot(f)
g.save(destdir + "f_ccm3.png")

f=phi_fs_1(2)
g=plot(f)
g.save(destdir + "f_fs1.png")

f=phi_fs_1(5)
g=plot(f)
g.save(destdir + "f_fs2.png")

f=phi_fs_1(11)
g=plot(f)
g.save(destdir + "f_fs3.png")

f=phi_vb_2(2)
g=plot(f)
g.save(destdir + "f_vb1.png")

f=phi_vb_2(5)
g=plot(f)
g.save(destdir + "f_vb2.png")

f=phi_vb_2(11)
g=plot(f)
g.save(destdir + "f_vb3.png")

f=phi_ll_1(3/2,3)
g=plot(f)
g.save(destdir + "f_ll11.png")

f=phi_ll_1(7/3,3)
g=plot(f)
g.save(destdir + "f_ll12.png")

f=phi_ll_1(13/4,3)
g=plot(f)
g.save(destdir + "f_ll13.png")

f=phi_ll_2(3/2,3)
g=plot(f)
g.save(destdir + "f_ll21.png")

f=phi_ll_2(7/3,3)
g=plot(f)
g.save(destdir + "f_ll22.png")

f=phi_ll_2(13/4,3)
g=plot(f)
g.save(destdir + "f_ll23.png")

f=phi_dg_1(3/2,3)
g=plot(f)
g.save(destdir + "f_dg1.png")

f=phi_dg_1(7/3,3)
g=plot(f)
g.save(destdir + "f_dg2.png")

f=phi_dg_1(13/4,3)
g=plot(f)
g.save(destdir + "f_dg3.png")


f=phi_bj_1(7/3)
g=plot_2d_diagram_dff(f)
g.save(destdir + "maximal_function_example.png")

f=gmic(1/5)
g=plot_2d_diagram_dff(f)
g.save(destdir + "nonmaximal_function_example.png")

f=phi_simple(7/3)
g=plot_2d_diagram_with_cones_dff(f)
g.save(destdir + "discontinuous_function_example.png")

f=phi_bj_1(5/2)
g=plot_2d_diagram_dff_no_lable(f,show_projections=False)
g.save(destdir + "maximality_test_example.png")   # used in dff-paper-long-version

f=w_2slope_3covered()
g=plot_2d_diagram_dff_no_lable(f, colorful=True)
g.save(destdir + "continuous_2slope_3covered.png")


# Missing scripts for graphics used in dff-paper-long-version
## Figure 2: Graphs of $\phi_{BJ,1}$ \cite[Example 3.1]{alves-clautiaux-valerio-rietz-2016:dual-feasible-book} for $C=3/2$ (left) and $C=7/3$ (right). }
## 2slope_dff_example1.png
## 2slope_dff_example2.png

## Figure 3
## w_2slope_3covered_nonextreme-2d_diagram.png

## Figure 4
## w_2slope_3covered_nonextreme-perturbations.png

## Figure 5
##phi_s_delta.png
## The graph of $\phi_{s,\delta}$ for $s=\frac{1}{5}$ and $\delta=2$

## Not used, no script, but in Dropbox: 2d_diagram_phi_s=2_delta=1_5.png
