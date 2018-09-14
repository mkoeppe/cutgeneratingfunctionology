import igp
from igp import *

#destdir = "survey_graphics/algo_paper_graphics/"
destdir = "/Users/mkoeppe/w/papers/basu-hildebrand-koeppe-papers/algo-paper/graphics-for-algo-paper/"
ftype = ".png"

logging.disable(logging.INFO)

igp.plot = plot_no_legend
igp.plot_kwds_hook = plot_kwds_hook_no_legend

# Two independently configurable style options for moves diagrams.
igp.show_translations_and_reflections_separately = False
igp.show_translations_and_reflections_by_color = True

igp.show_covered_components_as_rectangles = True
igp.show_moves_with_discontinuity_markers = True   ## do we want this?

## ######## not used ###
## g = line([(0,0), (0.2,0)],color='black')
## g += plot_limit_cone_of_vertex(0, 0, epstriple_to_cone((1, 1, 1)), color =  "mediumspringgreen", r = 0.03)
## g += plot_limit_cone_of_vertex(0.2, 0, epstriple_to_cone((-1, 1, -1)), color =  "mediumspringgreen", r = 0.03)
## g.show(axes=False)
## g.save("additive_edge_disct_case_1.png", axes=False, figsize=2)

## g = line([(0,0), (0.2,0)],color='black')
## g += plot_limit_cone_of_vertex(0, 0, epstriple_to_cone((1, 0, 1)), color =  "mediumspringgreen", r = 0.03)
## g += plot_limit_cone_of_vertex(0.2, 0, epstriple_to_cone((-1, 0, -1)), color =  "mediumspringgreen", r = 0.03)
## g.show(axes=False, figsize=2)
## g.save("additive_edge_disct_case_2.png", axes=False, figsize=2)

## g = line([(0,0), (0.2,0)],color='black')
## g += plot_limit_cone_of_vertex(0, 0, epstriple_to_cone((1, -1, 1)), color =  "mediumspringgreen", r = 0.03)
## g += plot_limit_cone_of_vertex(0.2, 0, epstriple_to_cone((-1, -1, -1)), color =  "mediumspringgreen", r = 0.03)
## g.show(axes=False, figsize=2)
## g.save("additive_edge_disct_case_3.png", axes=False, figsize=2)




######  algo-paper ####

## Some standard functions with happy moves diagrams!!
names = 'equiv7_example_1', 'minimal_no_covered_interval' # essential
names += 'kzh_7_slope_1', 'example7slopecoarse2', 'not_extreme_1', 'drlm_not_extreme_1', 'bhk_irrational_extreme_limit_to_rational_nonextreme', 'hildebrand_5_slope_22_1', 'rlm_dpl1_extreme_3a', 'hildebrand_discont_3_slope_1', 'zhou_two_sided_discontinuous_cannot_assume_any_continuity'  # additional, we will select later

for name in names:
    h = eval(name)()
    extremality_test(h, show_all_perturbations=True, show_plots=destdir + "%s-%%s.pdf" % name)
    show_plot(h._completion.plot(), destdir + "%s-%%s.pdf" % name, tag='completion-final')

## # crazy perturbation exists for h.
## t1, t2 = nice_field_values([1, sqrt(2)])
## rnf = t1.parent().fraction_field()

## # h = FastPiecewise([singleton_piece(0, 0),open_piece((0,t1/2),(t1/2, t1/2)), singleton_piece(t1/2, t1), open_piece((t1/2,t1/2),(t1,t1/2)), singleton_piece(t1,0)])
## # h1 = FastPiecewise([singleton_piece(0, 0),open_piece((0, t1/2), (t1/2, t1/2)),closed_piece((t1/2, t1),(t1,0))])

## generators = [t1, t2]
## pwl = piecewise_function_from_breakpoints_and_slopes([rnf(0),rnf(1)],[rnf(0)])
## crazy_piece_1 = CrazyPiece((rnf(0), rnf(1/4)), generators, [(rnf(1/8), rnf(1)), (rnf(t2/8), rnf(-1))])
## crazy_piece_2 = CrazyPiece((rnf(1/4), rnf(1/2)), generators, [(rnf(1/2-t2/8), rnf(1)), (rnf(1/2-1/8), rnf(-1))])
## cp = PiecewiseCrazyFunction(pwl, [crazy_piece_1, crazy_piece_2])


## find_epsilon_for_crazy_perturbation(h, cp)
## INFO: 2017-04-19 16:03:35,927 Rational case.
## 0.166666666666667

## find_epsilon_for_crazy_perturbation(h1, cp)
## 0

############################################################
### Split lemma diagrams
############################################################
save_show_translations_and_reflections_separately = igp.show_translations_and_reflections_separately
igp.show_translations_and_reflections_separately = False # override -- show translations only!

plot_background = polygon2d([[0,0], [0,1], [1,1], [1,0]], fill=False, color='grey')
t1 = 2/15 #1/10
t2 = 3/15 #3/10
l1 = 4/15
u1 = 3/4-t1
l2 = 1/3
u2 = 14/15-t2
m1 = FunctionalDirectedMove([open_interval(l1, u1)], (1, t1))
m2 = FunctionalDirectedMove([open_interval(l2, u2)], (1, t2))
l = l1; ll = l2; uu = 3/4; u = 14/15
c = DirectedMoveCompositionCompletion([m1, m2], show_plots=False, plot_background=plot_background)
c.add_backward_moves()
g = c.plot()+ line([(l,0),(l,u)], color='black', linestyle=':') +line([(ll,0),(ll,uu)], color='black', linestyle=':') + line([(uu,0),(uu,uu)], color='black', linestyle=':') + line([(u,0),(u,u)], color='black', linestyle=':') + line([(0,t1),(l,l+t1)], color='black', linestyle=':') + line([(0,t2),(ll,ll+t2)], color='black', linestyle=':') + line([(l,l),(u,l)], color='black', linestyle=':') + line([(ll,ll),(uu,ll)], color='black', linestyle=':')+line([(l,u),(u,u)], color='black', linestyle=':') +line([(ll,uu),(uu,uu)], color='black', linestyle=':') + text("$\\tau_1$", (0.5,0.58), axis_coords=False) +text("$\\tau_2$", (0.62,0.88), axis_coords=False) + text("$\\tau_1^{-1}$", (0.58,0.5), axis_coords=False) +text("$\\tau_2^{-1}$", (0.88, 0.62), axis_coords=False)
tkx = [0,l,ll,uu,u,1]
tky = [0,t1,t2,1]
tkxf = ["$0$","$l$","$l'$","$u'$","$u$","$1$"]
tkyf = ["$0$","$t_1$","$t_2$","$1$"]
#g.show(ticks=[tkx,tky], tick_formatter=[tkxf, tkyf], figsize=5)
g.save(destdir+"strip-lemma.png",ticks=[tkx,tky], tick_formatter=[tkxf, tkyf], figsize=5)
c.complete()
g1 = c.plot()+line([(l,0),(l,u)], color='black', linestyle=':') + line([(u,0),(u,u)], color='black', linestyle=':') + line([(0,l),(u,l)], color='black', linestyle=':') +line([(0,u),(u,u)], color='black', linestyle=':')
#g1.show(ticks=[[0,l,u,1],[0,l,u,1]],tick_formatter=[["$0$","$l$","$u$","$1$"], ["$0$","$l$","$u$","$1$"]],figsize=5)
g1.save(destdir+"strip-lemma-dependent.png",ticks=[[0,l,u,1],[0,l,u,1]],tick_formatter=[["$0$","$l$","$u$","$1$"], ["$0$","$l$","$u$","$1$"]],figsize=5)
g2 = g1 + polygon(((l,l), (l,u), (u,u), (u,l)), color="cyan")
#g2.show(ticks=[[0,l,u,1],[0,l,u,1]],tick_formatter=[["$0$","$l$","$u$","$1$"], ["$0$","$l$","$u$","$1$"]],figsize=5)
g2.save(destdir+"strip-lemma-independent.png",ticks=[[0,l,u,1],[0,l,u,1]],tick_formatter=[["$0$","$l$","$u$","$1$"], ["$0$","$l$","$u$","$1$"]],figsize=5)
igp.show_translations_and_reflections_separately = save_show_translations_and_reflections_separately

############################################################
## Other moves diagrams:
############################################################

# 2d diagram
plot_background = polygon2d([[0,0], [0,1], [1,1], [1,0]], fill=False, color='grey')
#a = 3/10; b=45/100; t1=1/4; t2=1-1/10; r=7/6
K.<r,a,b,t1,t2> = ParametricRealField([7/6, 3/10, 45/100, 1/4, 1-1/10])  # r comes first, good term order for labels!
F1 = Face([[a,b],[t1],[a+t1,b+t1]])
F1v = Face([[t1],[a,b],[a+t1,b+t1]])
F2 = Face([[a,b],[t2],[a+t2,b+t2]])
F2v = Face([[t2],[a,b],[a+t2,b+t2]])
F3 = Face([[a,b],[r-b,r-a],[r]])
F3r = Face([[r-b,r-a],[a,b],[r]])
IJK_kwds = [{'alpha': 0.350000000000000, 'color': 'grey', 'zorder': -10},{'alpha': 0.350000000000000, 'color': 'grey', 'zorder': -10},{'alpha': 0.350000000000000, 'color': 'grey', 'zorder': -10}]

g = plot_background + F1.plot()+F1v.plot()+F2.plot()+F2v.plot()+F3.plot()+F3r.plot()
g += line([(a, 0), (a, 1)], linestyle=':', color='grey')
g += line([(b, 0), (b, 1)], linestyle=':', color='grey')
I, J, K = F1.minimal_triple
g += polygon([(I[0], 1), (I[1], 1), (I[1], 1 + proj_plot_width), (I[0], 1 + proj_plot_width)], **IJK_kwds[0])
g += polygon([(K[0], 0), (K[1], 0), (K[1] + proj_plot_width, -proj_plot_width), (K[0] + proj_plot_width, -proj_plot_width)], **IJK_kwds[2])
g += line([(t1, a), (K[0], 0)], linestyle=':', color='grey')
g += line([(t1, b), (K[1], 0)], linestyle=':', color='grey')
I, J, K = F2.minimal_triple
g += polygon([(I[0], 1), (I[1], 1), (I[1], 1 + proj_plot_width), (I[0], 1 + proj_plot_width)], **IJK_kwds[0])
g += polygon([(1, K[0]-1), (1, K[1]-1), (1 + proj_plot_width, K[1] - 1 - proj_plot_width), (1 + proj_plot_width, K[0] - 1 - proj_plot_width)], **IJK_kwds[2])
g += line([(1, K[0]-1), (a, t2)], linestyle=':', color='grey')
g += line([(1, K[1]-1), (b, t2)], linestyle=':', color='grey')
g += line([(0, t1), (a, t1)], linestyle=':', color='grey')
g += line([(0, t2), (a, t2)], linestyle=':', color='grey')
g += line([(0, 0), (1, 1)], linestyle='-.', color='grey')
g += text("$F_1$", ((a+b)/2,t1-0.02), axis_coords=False, vertical_alignment='top',horizontal_alignment='center', color='black')
g += text("$F_2$", ((a+b)/2,t2+0.02), axis_coords=False, vertical_alignment='bottom',horizontal_alignment='center', color='black')
g += text("$F_1'$", (t1-0.02, (a+b)/2), axis_coords=False, vertical_alignment='center',horizontal_alignment='right', color='black')
g += text("$F_2'$", (t2+0.02, (a+b)/2), axis_coords=False, vertical_alignment='top',horizontal_alignment='left', color='black')
g += text("$p_1(F_1)=p_1(F_2)=p_1(F_3)$", ((a+b)/2,1+0.03), axis_coords=False, vertical_alignment='bottom',horizontal_alignment='center', color='black')
g += text("$p_3(F_1)$", ((a+b)/2+t1,0-0.02), axis_coords=False, vertical_alignment='top',horizontal_alignment='center', color='black')
g += text("$p_3(F_2)$", (1+0.03,(a+b)/2+t2-1), axis_coords=False, vertical_alignment='center',horizontal_alignment='left', color='black')
g += line([(b, r-b), (r-b,b)], linestyle=':', color='grey')+ line([(r-a, a),(r, 0)], linestyle=':', color='grey')
g += text("$F_3$", ((a+b)/2,r-(a+b)/2-0.02), axis_coords=False, vertical_alignment='top',horizontal_alignment='center', color='black')
g += text("$F_3'$", (r-(a+b)/2,(a+b)/2-0.02), axis_coords=False, vertical_alignment='top',horizontal_alignment='right', color='black')
tkx = [0,a,b,1,r]
tky = [0,t1,t2,1]
tkxf = ["$0$","$a$","$b$","$1$","$r$"]
tkyf = ["$0$","$t_1$","$t_2$","$1$"]
g += text("", (0,t2-1), axis_coords=False, vertical_alignment='center',horizontal_alignment='left', color='black')
#g.show(ticks=[tkx,tky], tick_formatter=[tkxf, tkyf], figsize=5)
g.save(destdir+"moves-diagram-edges.png",ticks=[tkx,tky], tick_formatter=[tkxf, tkyf], figsize=5)

# moves diagram
plot_background = polygon2d([[0,0], [0,1], [1,1], [1,0]], fill=False, color='grey')
fdms = [face.functional_directed_move() for face in [F1,F1v,F2,F2v,F3]]
gg = plot_background
# perhaps split up into extra_translation_graphics, extra_reflection_graphics...
gg += line([(0, 0), (1, 1)], linestyle='-.', color='grey') + line([(a, 0), (a, 1)], linestyle=':', color='grey') + line([(b, 0), (b, 1)], linestyle=':', color='grey') + line([(0, a+t1), (a, a+t1)], linestyle=':', color='grey') + line([(0,b+t1), (b, b+t1)], linestyle=':', color='grey') + line([(0, a+t2-1), (a, a+t2-1)], linestyle=':', color='grey') + line([(0,b+t2-1), (b, b+t2-1)], linestyle=':', color='grey') + line([(0, t1), (a, a+t1)], linestyle=':', color='grey') + line([(0, t2-1), (a, a+t2-1)], linestyle=':', color='grey')
tau1plus_label = text("$\\tau_1^+$", ((a+b)/2,(a+b)/2+t1+0.02), axis_coords=False, vertical_alignment='center',horizontal_alignment='right', color='black')
gg += tau1plus_label
tau2plus_label = text("$\\tau_2^+$", ((a+b)/2,(a+b)/2+t2-1-0.02), axis_coords=False, vertical_alignment='top',horizontal_alignment='center', color='black')
gg += tau2plus_label
tau1minus_label = text("$\\tau_1^-$", ((a+b)/2+t1+0.02, (a+b)/2), axis_coords=False, vertical_alignment='top',horizontal_alignment='left', color='black')
gg += tau1minus_label
tau2minus_label = text("$\\tau_2^-$", ((a+b)/2+t2-1, (a+b)/2+0.02), axis_coords=False, vertical_alignment='center',horizontal_alignment='right', color='black')
gg += tau2minus_label
rhor_ab_label = text("$\\rho_r|_{(a,b)}$", ((a+b)/2,r-(a+b)/2+0.02), axis_coords=False, vertical_alignment='bottom',horizontal_alignment='left', color='black')
gg += rhor_ab_label
rhor_ab_inverse_label = text("$(\\rho_r|_{(a,b)})^{-1}$", (r-(a+b)/2-0.02,(a+b)/2), axis_coords=False, vertical_alignment='bottom',horizontal_alignment='left', color='black')
gg += rhor_ab_inverse_label
gg += line([(b, r-b), (r, 0)], linestyle=':', color='grey', zorder=-10)
c = DirectedMoveCompositionCompletion(fdms, plot_background=gg)
c.add_backward_moves()
mtkx = [0,a,b,1,r]
mtky = [t2-1,0,a+t2-1,t1, b+t2-1,a+t1,b+t1,1]
mtkxf = ["$0$","$a$","$b$","$1$","$r$"]
mtkyf = ["$t_2-1$","$0$","$a+t_2-1$","$t_1$","$b+t_2-1$","$a+t_1$","$b+t_1$","$1$"]
ticks_graphics = text("", ((a+b)/2,1+0.03), axis_coords=False, vertical_alignment='bottom',horizontal_alignment='center', color='black',ticks=[mtkx,mtky], tick_formatter=[mtkxf, mtkyf])
gg = c.plot(extra_graphics=ticks_graphics)
#gg.show(figsize=5)
gg.save(destdir+"moves-diagram.png", figsize=5)

# separate moves
tau1 = F1.functional_directed_move()
tau2 = F2.functional_directed_move()
rhoab = FunctionalDirectedMove([open_interval(a, b)], F3.functional_directed_move().directed_move)

def save_move_plot(move, name):
    g = Graphics()
    g += plot_background
    if not move:
        move = FastPiecewise([])
    g += move.plot(**ticks_keywords(move, extra_xticks=[1], extra_yticks=[1]))
    g.save(destdir + name + ".png", figsize=3)

save_move_plot(tau1, "move_tau1+")
save_move_plot(~tau1, "move_tau1-")
save_move_plot(tau2, "move_tau2+")
save_move_plot(~tau2, "move_tau2-")
save_move_plot(rhoab, "move_rhoab+")
save_move_plot(~rhoab, "move_rhoab-")

# some compositions
save_move_plot(tau1 * (~tau2), "move_tau1+_o_tau2-")
save_move_plot(rhoab * (~tau2), "move_rhoab+_o_tau2-")
save_move_plot(tau1 * (~rhoab), "move_tau1+_o_rhoab-")

# restrictions of identity
save_move_plot(tau1 * (~tau1), "move_tau1+_o_tau1-")
save_move_plot((~tau1) * tau1, "move_tau1-_o_tau1+")

############################################################
## Moves and covered intervals
############################################################

plot_background = polygon2d([[0,0], [0,1], [1,1], [1,0]], fill=False, color='grey')

t1 = 2/15
tau1 = FunctionalDirectedMove([open_interval(6/15, 11/15)], (1, t1))
moves = [tau1]
comp = [ open_interval(7/15, 10/15) ]
fname = destdir + 'reduces_moves_by_components_ex1' + "-%s.png"
c = DirectedMoveCompositionCompletion(moves, [comp],
                                      show_plots=fname, plot_background=plot_background)
#c.add_backward_moves()
show_plot(plot_background + plot_covered_components_as_rectangles([comp]) + sum(m.plot() for m in moves), fname, tag='completion-unreduced')
show_plot(c.plot(), fname, tag='completion-initial')


t1 = 2/15
tau1a = FunctionalDirectedMove([open_interval(6/15, 7/15)], (1, t1))
tau1b = FunctionalDirectedMove([open_interval(8/15, 11/15)], (1, t1))
comp = [ open_interval(7/15, 10/15) ]
fname = destdir + 'extend_components_by_moves_ex1' + "-%s.png"
c = DirectedMoveCompositionCompletion([tau1a, tau1b], [comp],
                                      show_plots=fname, plot_background=plot_background)
c.complete()
show_plot(c.plot(), fname, tag='completion-final')

t1 = 2/15
tau1a = FunctionalDirectedMove([open_interval(4/15, 7/15)], (1, t1))
tau1b = FunctionalDirectedMove([open_interval(8/15, 11/15)], (1, t1))
comp = [ open_interval(7/15, 10/15) ]
fname = destdir + 'extend_components_by_moves_ex1a' + "-%s.png"
c = DirectedMoveCompositionCompletion([tau1a, tau1b], [comp],
                                      show_plots=fname, plot_background=plot_background)
c.complete()
show_plot(c.plot(), fname, tag='completion-final')

t1 = 5/15
tau1a = FunctionalDirectedMove([open_interval(2/15, 3/15)], (1, t1))
tau1b = FunctionalDirectedMove([open_interval(4/15, 5/15)], (1, t1))
comp = [ open_interval(3/15, 4/15), (7/15, 10/15) ]
fname = destdir + 'extend_components_by_moves_ex2' + "-%s.png"
c = DirectedMoveCompositionCompletion([tau1a, tau1b], [comp],
                                      show_plots=fname, plot_background=plot_background)
c.complete()
show_plot(c.plot(), fname, tag='completion-final')

t1 = 6/15
tau1a = FunctionalDirectedMove([open_interval(0/15, 1/15)], (1, t1))
tau1b = FunctionalDirectedMove([open_interval(4/15, 5/15)], (1, t1))
comp = [ open_interval(1/15, 4/15), (7/15, 10/15) ]
fname = destdir + 'extend_components_by_moves_ex3' + "-%s.png"
c = DirectedMoveCompositionCompletion([tau1a, tau1b], [comp],
                                      show_plots=fname,
                                      plot_background=plot_background)
c.complete()
show_plot(c.plot(), fname, tag='completion-final')

t1 = 5/15
tau1 = FunctionalDirectedMove([open_interval(2/15, 5/15)], (1, t1))
comp = [ open_interval(3/15, 4/15) ]
fname = destdir + 'extend_components_by_moves_ex4' + "-%s.png"
c = DirectedMoveCompositionCompletion([tau1], [comp],
                                      show_plots=fname, plot_background=plot_background)
c.complete()
show_plot(c.plot(), fname, tag='completion-final')





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



## Face sampling

def sampled_interval(I, num_samples=5):
    a, b = interval_to_endpoints(I)
    delta = (b - a) / num_samples
    return [ a + i * delta for i in range(num_samples+1) ]

def sampled_face(F, num_samples=5):
    I, J, K = F.minimal_triple
    return [ Face([[x], J, K]) for x in sampled_interval(I, num_samples) ] \
           + [ Face([I, [y], K]) for y in sampled_interval(J, num_samples) ] \
           + [ Face([I, J, [z]]) for z in sampled_interval(K, num_samples) ]

def symmetric_sampled_faces(F, Fprime, num_samples=4):
    E_list = [ E for E in sampled_face(F, num_samples) if E.is_1D() ]
    E_list += [ E for E in sampled_face(Fprime, num_samples) if E.is_1D() ]
    return E_list

def ticks_keywords_for_faces(faces):
    L = []
    for F in faces:
        L += flatten(F.minimal_triple)
    L = unique_list(sorted(L))
    tick_formatter = [ "$%s$" % latex(x) for x in L ]
    return { 'ticks': (L, L), 'tick_formatter': (tick_formatter, tick_formatter) }

def plot_sampled_stuff(F_list, E_list, name):
    fname = destdir + name + "-%s.png"
    fname_sampled = destdir + name + "_sampled-%s.png"

    background = polygon(((0,0), (0,1), (1,1), (1,0)), color='grey', fill=False, aspect_ratio=1, zorder=-2, **ticks_keywords_for_faces(F_list))
    g = background + plot_projections_of_faces(additive_faces=F_list)
    g += sum(E.plot(edge_thickness=1) for E in E_list)
    g.save(fname_sampled % "2d_diagram", xmin=0, xmax=1, ymin=0, ymax=1, aspect_ratio=1,
       **ticks_keywords_for_faces(F_list))

    background = polygon(((0,0), (0,1), (1,1), (1,0)), color='grey', fill=False, aspect_ratio=1, zorder=-2, **ticks_keywords_for_faces(F_list))
    # FIXME: The ticks keywords don't seem to apply to the translation moves diagram....
    completion = DirectedMoveCompositionCompletion(fdms=[ E.functional_directed_move() for E in E_list ], show_plots=fname_sampled, plot_background=background, show_zero_perturbation=False)
    #show(completion.plot())
    completion.complete()
    show_plot(completion.plot(), fname_sampled, tag='completion-final')

    # unsampled
    completion = DirectedMoveCompositionCompletion(covered_components=[ F.covered_component() for F in F_list ], show_plots=fname, plot_background=background, show_zero_perturbation=False)
    completion.complete()
    show_plot(completion.plot(), fname, tag='completion-final')


###############
I, J, K = [2/19, 3/19], [6/19, 7/19], [8/19, 9/19]
F = Face([I, J, K]) # lower triangle
Fprime = Face([J, I, K]) # swapped
F_list = [F, Fprime]
E_list = symmetric_sampled_faces(F, Fprime)
plot_sampled_stuff(F_list, E_list, 'triangle')

#################
I2, J2, K2 = [4/19, 5/19], [11/19, 12/19], [16/19, 17/19]
F2 = Face([I2, J2, K2]) # upper triangle
F2prime = Face([J2, I2, K2]) # swapped
F2_list = F_list + [F2, F2prime]  # lower and upper together (separate components)
E2_list = E_list + symmetric_sampled_faces(F2, F2prime)
plot_sampled_stuff(F2_list, E2_list, 'two_triangles')

#################
I, J, K = [3/19, 7/19], [1/19, 2/19], [4/19, 6/19]
F = Face([I, J, K]) # lower triangle
Fprime = Face([J, I, K]) # swapped
F_list = [F, Fprime]
E_list = symmetric_sampled_faces(F, Fprime)
plot_sampled_stuff(F_list, E_list, 'quadrilateral_overlapping_projections')
