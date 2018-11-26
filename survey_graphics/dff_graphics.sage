f=phi_simple(3/2)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_simple1.png")

f=phi_simple(7/3)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_simple2.png")

f=phi_simple(13/4)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_simple3.png")

f=phi_bj_1(3/2)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_bj1.png")

f=phi_bj_1(7/3)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_bj2.png")

f=phi_bj_1(13/4)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_bj3.png")

f=phi_ccm_1(3/2)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_ccm1.png")

f=phi_ccm_1(7/3)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_ccm2.png")

f=phi_ccm_1(13/4)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_ccm3.png")

f=phi_fs_1(2)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_fs1.png")

f=phi_fs_1(5)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_fs2.png")

f=phi_fs_1(11)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_fs3.png")

f=phi_vb_2(2)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_vb1.png")

f=phi_vb_2(5)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_vb2.png")

f=phi_vb_2(11)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_vb3.png")

f=phi_ll_1(3/2,3)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_ll11.png")

f=phi_ll_1(7/3,3)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_ll12.png")

f=phi_ll_1(13/4,3)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_ll13.png")

f=phi_ll_2(3/2,3)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_ll21.png")

f=phi_ll_2(7/3,3)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_ll22.png")

f=phi_ll_2(13/4,3)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_ll23.png")

f=phi_dg_1(3/2,3)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_dg1.png")

f=phi_dg_1(7/3,3)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_dg2.png")

f=phi_dg_1(13/4,3)
g=plot(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/f_dg3.png")


f=phi_bj_1(7/3)
g=plot_2d_diagram_dff(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/maximal_function_example.png")

f=gmic(1/5)
g=plot_2d_diagram_dff(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/nonmaximal_function_example.png")

f=phi_simple(7/3)
g=plot_2d_diagram_with_cones_dff(f)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/discontinuous_function_example.png")

f=phi_bj_1(5/2)
g=plot_2d_diagram_dff_no_lable(f,show_projections=False)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/maximality_test_example.png")

f=w_2slope_3covered()
g=plot_2d_diagram_dff_no_lable(f, colorful=True)
g.save("/Users/arsenalcrown/Dropbox/jiawei-bhk-papers/jiawei/dff-paper/dff_paper_graphics/continuous_2slope_3covered.png",figsize=25)


