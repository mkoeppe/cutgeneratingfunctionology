load("survey_graphics/graphics_for_algo_paper_init.sage")

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
load("survey_graphics/graphics_for_algo_paper_functions.sage")

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
load("survey_graphics/graphics_for_algo_paper_strip_lemma.sage")

############################################################
## Other moves diagrams:
############################################################
load("survey_graphics/graphics_for_algo_paper_moves.sage")

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

load("survey_graphics/graphics_for_algo_paper_face_sampling.sage")
