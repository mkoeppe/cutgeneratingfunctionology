import igp
from igp import *

destdir = "survey_graphics/algo_paper_graphics/"
ftype = ".png"

logging.disable(logging.INFO)


## ######## not used ###
## sage: g = line([(0,0), (0.2,0)],color='black')
## sage: g += plot_limit_cone_of_vertex(0, 0, epstriple_to_cone((1, 1, 1)), color =  "mediumspringgreen", r = 0.03)
## sage: g += plot_limit_cone_of_vertex(0.2, 0, epstriple_to_cone((-1, 1, -1)), color =  "mediumspringgreen", r = 0.03)
## sage: g.show(axes=False)
## sage: g.save("additive_edge_disct_case_1.png", axes=False, figsize=2)

## sage: g = line([(0,0), (0.2,0)],color='black')
## sage: g += plot_limit_cone_of_vertex(0, 0, epstriple_to_cone((1, 0, 1)), color =  "mediumspringgreen", r = 0.03)
## sage: g += plot_limit_cone_of_vertex(0.2, 0, epstriple_to_cone((-1, 0, -1)), color =  "mediumspringgreen", r = 0.03)
## sage: g.show(axes=False, figsize=2)
## sage: g.save("additive_edge_disct_case_2.png", axes=False, figsize=2)

## sage: g = line([(0,0), (0.2,0)],color='black')
## sage: g += plot_limit_cone_of_vertex(0, 0, epstriple_to_cone((1, -1, 1)), color =  "mediumspringgreen", r = 0.03)
## sage: g += plot_limit_cone_of_vertex(0.2, 0, epstriple_to_cone((-1, -1, -1)), color =  "mediumspringgreen", r = 0.03)
## sage: g.show(axes=False, figsize=2)
## sage: g.save("additive_edge_disct_case_3.png", axes=False, figsize=2)




######  algo-paper ####
sage: h1 = equiv7_example_1()
sage: h = minimal_no_covered_interval()

### copy some plot formatting from graphics_for-survey ###...
sage: extremality_test(h1, show_plots="one-sided-discontinuous-%s.png", show_all
....: _perturbations=True)
sage: extremality_test(h, show_plots="two-sided-discontinuous-%s.png", show_all_
....: perturbations=True)

## # crazy perturbation exists for h.
## sage: t1, t2 = nice_field_values([1, sqrt(2)])
## sage: rnf = t1.parent().fraction_field()

## # sage: h = FastPiecewise([singleton_piece(0, 0),open_piece((0,t1/2),(t1/2, t1/2)), singleton_piece(t1/2, t1), open_piece((t1/2,t1/2),(t1,t1/2)), singleton_piece(t1,0)])
## # sage: h1 = FastPiecewise([singleton_piece(0, 0),open_piece((0, t1/2), (t1/2, t1/2)),closed_piece((t1/2, t1),(t1,0))])

## sage: generators = [t1, t2]
## sage: pwl = piecewise_function_from_breakpoints_and_slopes([rnf(0),rnf(1)],[rnf(0)])
## sage: crazy_piece_1 = CrazyPiece((rnf(0), rnf(1/4)), generators, [(rnf(1/8), rnf(1)), (rnf(t2/8), rnf(-1))])
## sage: crazy_piece_2 = CrazyPiece((rnf(1/4), rnf(1/2)), generators, [(rnf(1/2-t2/8), rnf(1)), (rnf(1/2-1/8), rnf(-1))])
## sage: cp = PiecewiseCrazyFunction(pwl, [crazy_piece_1, crazy_piece_2])


## sage: find_epsilon_for_crazy_perturbation(h, cp)
## INFO: 2017-04-19 16:03:35,927 Rational case.
## 0.166666666666667

## sage: find_epsilon_for_crazy_perturbation(h1, cp)
## 0

sage: plot_background = polygon2d([[0,0], [0,1], [1,1], [1,0]], fill=False, color='grey')
sage: t1 = 2/15 #1/10
sage: t2 = 3/15 #3/10
sage: l1 = 4/15
sage: u1 = 3/4-t1
sage: l2 = 1/3
sage: u2 = 14/15-t2
sage: m1 = FunctionalDirectedMove([open_interval(l1, u1)], (1, t1))
sage: m2 = FunctionalDirectedMove([open_interval(l2, u2)], (1, t2))
sage: l = l1; ll = l2; uu = 3/4; u = 14/15
sage: c = DirectedMoveCompositionCompletion([m1, m2], show_plots=True, plot_background=plot_background)
sage: c.add_backward_moves()
sage: g = c.plot()+ line([(l,0),(l,u)], color='black', linestyle=':') +line([(ll,0),(ll,uu)], color='black', linestyle=':') + line([(uu,0),(uu,uu)], color='black', linestyle=':') + line([(u,0),(u,u)], color='black', linestyle=':') + line([(0,t1),(l,l+t1)], color='black', linestyle=':') + line([(0,t2),(ll,ll+t2)], color='black', linestyle=':') + line([(l,l),(u,l)], color='black', linestyle=':') + line([(ll,ll),(uu,ll)], color='black', linestyle=':')+line([(l,u),(u,u)], color='black', linestyle=':') +line([(ll,uu),(uu,uu)], color='black', linestyle=':') + text("$\\tau_1$", (0.5,0.58), axis_coords=False) +text("$\\tau_2$", (0.62,0.88), axis_coords=False) + text("$\\tau_1^{-1}$", (0.58,0.5), axis_coords=False) +text("$\\tau_2^{-1}$", (0.88, 0.62), axis_coords=False)
sage: tkx = [0,l,ll,uu,u,1]
sage: tky = [0,t1,t2,1]
sage: tkxf = ["$0$","$l$","$l'$","$u'$","$u$","$1$"]
sage: tkyf = ["$0$","$t_1$","$t_2$","$1$"]
sage: g.show(ticks=[tkx,tky], tick_formatter=[tkxf, tkyf], figsize=5)
sage: g.save("strip-lemma.png",ticks=[tkx,tky], tick_formatter=[tkxf, tkyf], figsize=5)
sage: c.complete()
sage: g1 = c.plot()+line([(l,0),(l,u)], color='black', linestyle=':') + line([(u,0),(u,u)], color='black', linestyle=':') + line([(0,l),(u,l)], color='black', linestyle=':') +line([(0,u),(u,u)], color='black', linestyle=':')
sage: g1.show(ticks=[[0,l,u,1],[0,l,u,1]],tick_formatter=[["$0$","$l$","$u$","$1$"], ["$0$","$l$","$u$","$1$"]],figsize=5)
sage: g1.save("strip-lemma-dependent.png",ticks=[[0,l,u,1],[0,l,u,1]],tick_formatter=[["$0$","$l$","$u$","$1$"], ["$0$","$l$","$u$","$1$"]],figsize=5)
sage: g2 = g1 + polygon(((l,l), (l,u), (u,u), (u,l)), color="cyan")
sage: g2.show(ticks=[[0,l,u,1],[0,l,u,1]],tick_formatter=[["$0$","$l$","$u$","$1$"], ["$0$","$l$","$u$","$1$"]],figsize=5)
sage: g2.save("strip-lemma-independent.png",ticks=[[0,l,u,1],[0,l,u,1]],tick_formatter=[["$0$","$l$","$u$","$1$"], ["$0$","$l$","$u$","$1$"]],figsize=5)


sage: plot_background = polygon2d([[0,0], [0,1], [1,1], [1,0]], fill=False, color='grey')
sage: a = 3/10; b=45/100; t1=1/4; t2=1-1/10; r=7/6
sage: F1 = Face([[a,b],[t1],[a+t1,b+t1]])
sage: F1v = Face([[t1],[a,b],[a+t1,b+t1]])
sage: F2 = Face([[a,b],[t2],[a+t2,b+t2]])
sage: F2v = Face([[t2],[a,b],[a+t2,b+t2]])
sage: F3 = Face([[a,b],[r-b,r-a],[r]])
sage: F3r = Face([[r-b,r-a],[a,b],[r]])
sage: IJK_kwds = [{'alpha': 0.350000000000000, 'color': 'grey', 'zorder': -10},{'alpha': 0.350000000000000, 'color': 'grey', 'zorder': -10},{'alpha': 0.350000000000000, 'color': 'grey', 'zorder': -10}]

sage: g = plot_background + F1.plot()+F1v.plot()+F2.plot()+F2v.plot()+F3.plot()+F3r.plot()
sage: g += line([(a, 0), (a, 1)], linestyle=':', color='grey')
sage: g += line([(b, 0), (b, 1)], linestyle=':', color='grey')
sage: I, J, K = F1.minimal_triple
sage: g += polygon([(I[0], 1), (I[1], 1), (I[1], 1 + proj_plot_width), (I[0], 1 + proj_plot_width)], **IJK_kwds[0])
sage: g += polygon([(K[0], 0), (K[1], 0), (K[1] + proj_plot_width, -proj_plot_width), (K[0] + proj_plot_width, -proj_plot_width)], **IJK_kwds[2])
sage: g += line([(t1, a), (K[0], 0)], linestyle=':', color='grey')
sage: g += line([(t1, b), (K[1], 0)], linestyle=':', color='grey')
sage: I, J, K = F2.minimal_triple
sage: g += polygon([(I[0], 1), (I[1], 1), (I[1], 1 + proj_plot_width), (I[0], 1 + proj_plot_width)], **IJK_kwds[0])
sage: g += polygon([(1, K[0]-1), (1, K[1]-1), (1 + proj_plot_width, K[1] - 1 - proj_plot_width), (1 + proj_plot_width, K[0] - 1 - proj_plot_width)], **IJK_kwds[2])
sage: g += line([(1, K[0]-1), (a, t2)], linestyle=':', color='grey')
sage: g += line([(1, K[1]-1), (b, t2)], linestyle=':', color='grey')
sage: g += line([(0, t1), (a, t1)], linestyle=':', color='grey')
sage: g += line([(0, t2), (a, t2)], linestyle=':', color='grey')
sage: g += line([(0, 0), (1, 1)], linestyle='-.', color='grey')
sage: g += text("$F_1$", ((a+b)/2,t1-0.02), axis_coords=False, vertical_alignment='top',horizontal_alignment='center', color='black')
sage: g += text("$F_2$", ((a+b)/2,t2+0.02), axis_coords=False, vertical_alignment='bottom',horizontal_alignment='center', color='black')
sage: g += text("$F_1'$", (t1-0.02, (a+b)/2), axis_coords=False, vertical_alignment='center',horizontal_alignment='right', color='black')
sage: g += text("$F_2'$", (t2+0.02, (a+b)/2), axis_coords=False, vertical_alignment='top',horizontal_alignment='left', color='black')
sage: g += text("$p_1(F_1)=p_1(F_2)=p_1(F_3)$", ((a+b)/2,1+0.03), axis_coords=False, vertical_alignment='bottom',horizontal_alignment='center', color='black')
sage: g += text("$p_3(F_1)$", ((a+b)/2+t1,0-0.02), axis_coords=False, vertical_alignment='top',horizontal_alignment='center', color='black')
sage: g += text("$p_3(F_2)$", (1+0.03,(a+b)/2+t2-1), axis_coords=False, vertical_alignment='center',horizontal_alignment='left', color='black')
sage: g += line([(b, r-b), (r-b,b)], linestyle=':', color='grey')+ line([(r-a, a),(r, 0)], linestyle=':', color='grey')
sage: g += text("$F_3$", ((a+b)/2,r-(a+b)/2-0.02), axis_coords=False, vertical_alignment='top',horizontal_alignment='center', color='black')
sage: g += text("$F_3'$", (r-(a+b)/2,(a+b)/2-0.02), axis_coords=False, vertical_alignment='top',horizontal_alignment='right', color='black')
sage: tkx = [0,a,b,1,r]
sage: tky = [0,t1,t2,1]
sage: tkxf = ["$0$","$a$","$b$","$1$","$r$"]
sage: tkyf = ["$0$","$t_1$","$t_2$","$1$"]
sage: g += text("", (0,t2-1), axis_coords=False, vertical_alignment='center',horizontal_alignment='left', color='black')
sage: g.show(ticks=[tkx,tky], tick_formatter=[tkxf, tkyf], figsize=5)
sage: g.save("moves-diagram-edges.png",ticks=[tkx,tky], tick_formatter=[tkxf, tkyf], figsize=5)


sage: fdms = [face.functional_directed_move() for face in [F1,F1v,F2,F2v,F3]]
sage: c = DirectedMoveCompositionCompletion(fdms, plot_background=plot_background)
sage: c.add_backward_moves()
sage: gg = c.plot() + line([(0, 0), (1, 1)], linestyle='-.', color='grey') + line([(a, 0), (a, 1)], linestyle=':', color='grey') + line([(b, 0), (b, 1)], linestyle=':', color='grey') + line([(0, a+t1), (a, a+t1)], linestyle=':', color='grey') + line([(0,b+t1), (b, b+t1)], linestyle=':', color='grey') + line([(0, a+t2-1), (a, a+t2-1)], linestyle=':', color='grey') + line([(0,b+t2-1), (b, b+t2-1)], linestyle=':', color='grey') + line([(0, t1), (a, a+t1)], linestyle=':', color='grey') + line([(0, t2-1), (a, a+t2-1)], linestyle=':', color='grey')
sage: gg += text("$\\tau_1^+$", ((a+b)/2,(a+b)/2+t1+0.02), axis_coords=False, vertical_alignment='center',horizontal_alignment='right', color='black')
sage: gg += text("$\\tau_2^+$", ((a+b)/2,(a+b)/2+t2-1-0.02), axis_coords=False, vertical_alignment='top',horizontal_alignment='center', color='black')
sage: gg += text("$\\tau_1^-$", ((a+b)/2+t1+0.02, (a+b)/2), axis_coords=False, vertical_alignment='top',horizontal_alignment='left', color='black')
sage: gg += text("$\\tau_2^-$", ((a+b)/2+t2-1, (a+b)/2+0.02), axis_coords=False, vertical_alignment='center',horizontal_alignment='right', color='black')
sage: gg += text("$\\rho_r|_{(a,b)}$", ((a+b)/2,r-(a+b)/2+0.02), axis_coords=False, vertical_alignment='bottom',horizontal_alignment='left', color='black')
sage: gg += text("$(\\rho_r|_{(a,b)})^{-1}$", (r-(a+b)/2-0.02,(a+b)/2), axis_coords=False, vertical_alignment='bottom',horizontal_alignment='left', color='black')
sage: gg += line([(b, r-b), (r, 0)], linestyle=':', color='grey', zorder=-10)
sage: mtkx = [0,a,b,1,r]
sage: mtky = [t2-1,0,a+t2-1,t1, b+t2-1,a+t1,b+t1,1]
sage: mtkxf = ["$0$","$a$","$b$","$1$","$r$"]
sage: mtkyf = ["$t_2-1$","$0$","$a+t_2-1$","$t_1$","$b+t_2-1$","$a+t_1$","$b+t_1$","$1$"]
sage: gg += text("", ((a+b)/2,1+0.03), axis_coords=False, vertical_alignment='bottom',horizontal_alignment='center', color='black')
sage: gg.show(ticks=[mtkx,mtky], tick_formatter=[mtkxf, mtkyf], figsize=5)
sage: gg.save("moves-diagram.png",ticks=[mtkx,mtky], tick_formatter=[mtkxf, mtkyf], figsize=5)

## ########## mip 2017 slides ###########
## igprainbow=igp.rainbow
## def dark_rainbow(num):
##     return ['darkblue', 'darkgreen', 'firebrick', 'darkcyan', 'darkmagenta', 'darkgrey', 'royalblue'][:num]
## igp.rainbow=dark_rainbow

## h = drlm_not_extreme_1()
## extremality_test(h)
## hb = h+ h._perturbations[0]*7/20
## hr = h -h._perturbations[0]*7/20
## g = plot(hb, color='blue')+plot(hr,color='red')+plot(h,color='black',thickness=2)
## xx =[i/7 for i in range(8)]; yy = [0,1]
## ticks=[xx,yy];tick_formatter=[["$%s$" % latex(x) for x in xx],["$%s$" % latex(x) for x in yy]]
## g.fontsize(10)
## g.show(ticks=ticks, tick_formatter=tick_formatter,aspect_ratio=1/6, show_legend=False, figsize=4)
## g.save("drlm_not_extreme_pert720.pdf",ticks=ticks, tick_formatter=tick_formatter,aspect_ratio=1/6, show_legend=False, figsize=4)


## h = automorphism(drlm_backward_3_slope(1/7, 1/7+196/2017))
## xx = h.end_points()
## ticks=[xx,yy];tick_formatter=[["$%s$" % latex(x) for x in xx],["$%s$" % latex(x) for x in yy]]
## g = plot(h,color='black',thickness=2)
## g.fontsize(10)
## g.show(ticks=ticks, tick_formatter=tick_formatter,aspect_ratio=1/6, show_legend=False, figsize=4)
## g.save("drlm_ext_big_grid.pdf",ticks=ticks, tick_formatter=tick_formatter,aspect_ratio=1/6, show_legend=False, figsize=4)


## bhk = bhk_irrational(f=4/5, d1=3/5, d2=1/10, a0=15/100, delta=(1/80, sqrt(2)/80))
## xx=[0,1]; yy = []
## ticks=[xx,yy];tick_formatter=[["$%s$" % latex(x) for x in xx],["$%s$" % latex(x) for x in yy]]
## gbhk = plot_with_colored_slopes(bhk)
## xs = bhk.end_points()
## a0 = xs[2]; f = xs[-2]; l = xs[8]; u = xs[9]; t1 = xs[4]-xs[2]; t2 = xs[6]-xs[2]
## g = gbhk+line([(l,0),(u,0)], thickness=2, color='magenta')+line([(f-l,0),(f-u,0)], thickness=2, color='magenta')
## g = g +  polygon2d([(l, 0.1), (u, 0.1),(u, 0.9), (l, 0.9)], color="yellow", alpha=0.2, ticks=[[],[]], zorder=-10) +  polygon2d([(f-l, 0.1), (f-u, 0.1),(f-u, 0.9), (f-l, 0.9)], color="yellow", alpha=0.2, ticks=[[],[]], zorder=-10)
## tx1 = text("$t_1=\\frac{1}{80}$", (a0,1/8), horizontal_alignment='left', color='black',fontsize=10)
## tx2 = text("$t_2 = \\frac{\\sqrt{2}}{80}$", (xs[4],1/2),horizontal_alignment='left', color='black',fontsize=10)
## g = g + tx1 + tx2
## d2 = t2-t1; ll = 31/100
## ar1 = arrow(path=[[(ll, 3/8),(ll+t1/2, 1/2), (ll+t1, 3/8)]],head=1, color='orange', width=1,arrowsize=1,zorder=2) + arrow(path=[[(ll+t1, 3/8),(ll+t1+t1/2, 1/2), (ll+t1+t1, 3/8)]],head=1, color = 'orange', width=1,arrowsize=1,zorder=2) + arrow(path=[[(ll+2*t1, 3/8),(ll+2*t1-d2/2, 1/4), (ll+2*t1 -d2, 3/8)]],head=1, color= 'orange', width=1,arrowsize=1,zorder=2)
## tx3 = text("densely covered", ((l+u)/2,0.1),color='magenta',fontsize=10,vertical_alignment='bottom')
## g = g+ar1+tx3
## g.save("bhk_densely_covered.pdf",ticks=ticks, tick_formatter=tick_formatter,aspect_ratio=1/8, show_legend=False)
