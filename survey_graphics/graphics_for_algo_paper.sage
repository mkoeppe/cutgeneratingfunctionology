######
h = zhou_two_sided_discontinuous_cannot_assume_any_continuity()
g = plot(h, color='black', thickness=2, **ticks_keywords(h))
g += plot(per, color='magenta',alpha=0)
g.save("two_sided_discontinuous.png", figsize=6)
extremality_test(h)
perturbation = h._perturbations[0]
per = rescale_to_amplitude(perturbation, 1/9)
g = plot(h - 4/9 * perturbation, color='red')
g += plot(h + 4/9 * perturbation, color='blue')
g += plot(h, color='black', thickness=2, **ticks_keywords(h))
g += plot(per, color='magenta')
g.save("two_sided_discontinuous_perturbations.png", figsize=6)



######## not used ###
sage: g = line([(0,0), (0.2,0)],color='black')
sage: g += plot_limit_cone_of_vertex(0, 0, epstriple_to_cone((1, 1, 1)), color =  "mediumspringgreen", r = 0.03)
sage: g += plot_limit_cone_of_vertex(0.2, 0, epstriple_to_cone((-1, 1, -1)), color =  "mediumspringgreen", r = 0.03)
sage: g.show(axes=False)
sage: g.save("additive_edge_disct_case_1.png", axes=False, figsize=2)

sage: g = line([(0,0), (0.2,0)],color='black')
sage: g += plot_limit_cone_of_vertex(0, 0, epstriple_to_cone((1, 0, 1)), color =  "mediumspringgreen", r = 0.03)
sage: g += plot_limit_cone_of_vertex(0.2, 0, epstriple_to_cone((-1, 0, -1)), color =  "mediumspringgreen", r = 0.03)
sage: g.show(axes=False, figsize=2)
sage: g.save("additive_edge_disct_case_2.png", axes=False, figsize=2)

sage: g = line([(0,0), (0.2,0)],color='black')
sage: g += plot_limit_cone_of_vertex(0, 0, epstriple_to_cone((1, -1, 1)), color =  "mediumspringgreen", r = 0.03)
sage: g += plot_limit_cone_of_vertex(0.2, 0, epstriple_to_cone((-1, -1, -1)), color =  "mediumspringgreen", r = 0.03)
sage: g.show(axes=False, figsize=2)
sage: g.save("additive_edge_disct_case_3.png", axes=False, figsize=2)


######### lifting example ######
# random_pwl
sage: h = FastPiecewise([[(QQ(0), 1/8), FastLinearFunction(QQ(8), QQ(0))], [(1/8, 1/4), FastLinearFunction(-16/3, 5/3)], [(1/4, 3/8), FastLinearFunction(-4/3, 2/3)], [(3/8, 1/2), FastLinearFunction(4/3, -1/3)], [(1/2, 5/8), FastLinearFunction(8/3, -QQ(1))], [(5/8, 3/4), FastLinearFunction(4/3, -1/6)], [(3/4, 7/8), FastLinearFunction(-4/3, 11/6)], [(7/8, QQ(1)), FastLinearFunction(-16/3, 16/3)]])
sage: plot_covered_intervals(h).save("lift_not_min_random.png",show_legend=False, figsize=4)

#sage: g = plot_2d_diagram(h, known_minimal=True)+plot_2d_diagram_with_cones(h)
sage: g = plot_2d_diagram(h)
sage: g.show(show_legend=False)
sage: g.save("lift_not_minimal_2d_diagram.png",show_legend=False)

sage: extremality_test(h, phase_1=True)

sage: g = plot(h,color='black', linestyle='-')+plot(h._perturbations[0]*4/3,color='magenta')+plot(h+h._perturbations[0]*4/3, color='blue', thickness=2)
sage: g.save("lift_phase_1.png")

sage: hh = lift(h) #hmin
sage: finite_dimensional_extremality_test(hh, show_all_perturbations=True)

sage: g = plot(hh,color='black')+plot(hh._perturbations[0]*4/3,color='magenta')+plot(hh._perturbations[1]*2, color='cyan')
sage: g.save("lift_min_has_two_finite_dim_pert.png")

sage: plot_2d_diagram(hh, colorful=True).save('lift_min_2d_diagram.png',show_legend=False)

sage: finite_dimensional_extremality_test(hh)
sage: g = plot(hh,color='black')+plot(hh._perturbations[0]*4/3,color='magenta')+plot(hh+hh._perturbations[0]*4/3,color='blue')
sage: g.save('lift_until_extreme_1.png')
sage: hhh = lift(hh)

sage: finite_dimensional_extremality_test(hhh)
sage: g = plot(hhh,color='black')+plot(hhh._perturbations[0]*4/3,color='magenta')+plot(hhh+hhh._perturbations[0]*4/3,color='blue')
sage: g.save('lift_until_extreme_2.png')


sage: finite_dimensional_extremality_test(hhhh)
True
sage: extremality_test(hhhh)
sage: g = plot(hhhh,color='black')+plot(hhhh._perturbations[0]/3,color='magenta')+plot(hhhh+hhhh._perturbations[0]/3,color='blue')
sage: g.save('lift_until_extreme_3.png')


sage: hhhhh = lift(hhhh)
sage: hhhhh == hhhh+hhhh._perturbations[0]/3
True
sage: extremality_test(hhhhh)
True


sage: lift_until_extreme(hh, show_plots=True, finite_dimensional_extremality_test=True)
INFO: 2016-04-20 15:52:22,315 Finite dimensional test: Solution space has dimension 2
INFO: 2016-04-20 15:52:22,324 Finding epsilon interval for perturbation... done.  Interval is [-2, 4/3]
INFO: 2016-04-20 15:52:22,331 Finding epsilon interval for perturbation... done.  Interval is [-8/3, 2]
sage: perturbs = hh._perturbations
sage: P = perturbation_polyhedron(hh, perturbs)
sage: P = perturbation_polyhedron(hh, perturbs)
sage: P.plot()
Launched png viewer for Graphics object consisting of 6 graphics primitives
sage: P.vertices()
(A vertex at (4/3, -4),
 A vertex at (-4/5, 12/5),
 A vertex at (4/3, 4/3),
 A vertex at (-20/9, -4/9))


sage: igprainbow = igp.rainbow
sage: def new_rainbow(num):
....:     return ['red','cyan','blue','orange','purple'][:num]
....: 
sage: igp.rainbow=new_rainbow

sage: verts = P.vertices()
sage: verts[0]
A vertex at (4/3, -4)
Launched png viewer for Graphics object consisting of 51 graphics primitives
sage: h1 = hh + perturbation_corresponding_to_vertex(perturbs, verts[0])
sage: plot_covered_intervals(h1).save('lift_pert_poly_1.png', figsize=4, show_legend=False)
#sage: plot(h1) # extreme
sage: verts[3]
A vertex at (-20/9, -4/9)
sage: h2 = hh + perturbation_corresponding_to_vertex(perturbs, verts[3])
sage: plot_covered_intervals(h2).save('lift_pert_poly_2.png', figsize=4, show_legend=False)
# sage: plot(h2) # extreme
sage: verts[1]
A vertex at (-4/5, 12/5)
sage: h3 = hh + perturbation_corresponding_to_vertex(perturbs, verts[1])
sage: plot_covered_intervals(h3).save('lift_pert_poly_3.png', figsize=4, show_legend=False)
#sage: plot(h3) #not extreme
sage: verts[2]
A vertex at (4/3, 4/3)
sage: h4 = hh + perturbation_corresponding_to_vertex(perturbs, verts[2])
sage: plot_covered_intervals(h4).save('lift_pert_poly_4.png', figsize=4, show_legend=False)
#sage: plot(h4) # not extreme

sage: P.plot(color='springgreen')
g = polygon([(4/3, -4), (-20/9, -4/9), (-4/5, 12/5),(4/3, 4/3)], color="orange", edgecolor='black', alpha=0.5) + point([(4/3, -4), (-20/9, -4/9)], color='blue', size=50) + points([(-4/5, 12/5),(4/3, 4/3)], color='red', size=50)
g.save('lift_pert_space_polyhedron.png')

#polygon([(-2, -8/3), (-2, 2), (4/3, 2),(4/3, -8/3)], color="blue", fill=None)


########## facet paper ###########

sage: h = hildebrand_discont_3_slope_1()
sage: hh = discontinuous_interpolation([0,1/2,5/8,7/8],[0,1,3/4,1/4],[0,1/2,3/4,1/4],[1/2,1,3/4,1/4])
#sage: h3 = restrict_to_finite_group(h, oversampling=3)
#sage: extremality_test(h3, True)
sage: g = plot_2d_diagram_additive_domain_sans_limits(h)
sage: g.show(show_legend=False)
sage: g.save('2d_simple_E_pi_extreme_not_facet.png', show_legend=False)
Launched png viewer for Graphics object consisting of 98 graphics primitives
sage: gg = plot_2d_diagram_additive_domain_sans_limits(hh)
sage: gg.show(show_legend=False)
sage: gg.save('2d_simple_E_pi_prime_extreme_not_facet.png', show_legend=False)

sage: fn = kzh_minimal_has_only_crazy_perturbation_1()
sage: gen = set(itertools.chain(generate_type_1_vertices(fn, operator.gt), \
                                generate_type_2_vertices(fn, operator.gt)))
sage: dp = [delta_pi_general(fn, x, y, (xeps, yeps, zeps)) for (x, y, z, xeps, yeps, zeps) in gen]
sage: igp.show_RNFElement_by_embedding=False
sage: min(dp)
19/23998


sage: h = kzh_minimal_has_only_crazy_perturbation_1()
sage: x = h.end_points()
sage: l=x[17]; u = x[18]; ll = x[19]; uu=x[20]; f=x[37]
sage: g = plot_2d_complex(h)
sage: g += polygon([(0, l), (0, u), (1, u), (1, l)], color='yellow',fill=True, zorder=-5)
sage: g += polygon([(0, ll), (0, uu), (1, uu), (1, ll)], color='yellow',fill=True, zorder=-5)
sage: g += polygon([(ll,0), (uu,0), (uu,1), (ll,1)], color='yellow',fill=True, zorder=-5)
sage: g += polygon([(l,0), (u,0), (u,1), (l,1)], color='yellow',fill=True, zorder=-5)
sage: g += polygon([(l,0), (u,0), (0,u), (0,l)], color='yellow',fill=True, zorder=-5)
sage: g += polygon([(ll,0), (uu,0), (0,uu), (0,ll)], color='yellow',fill=True, zorder=-5)
sage: g += polygon([(l,1), (u,1), (1,u), (1,l)], color='yellow',fill=True, zorder=-5)
sage: g += polygon([(ll,1), (uu,1), (1,uu), (1,ll)], color='yellow',fill=True, zorder=-5)
sage: g += polygon([(l,1), (u,1), (u,1-w)], color='red',fill=True, zorder=-5) 
sage: g += polygon([(ll,1), (uu,1), (uu,1-w)], color='red',fill=True, zorder=-5)
sage: g += polygon([(0,l), (0,u), (w,l)], color='red',fill=True, zorder=-5)
sage: g += polygon([(0,ll), (0,uu), (w,ll)], color='red',fill=True, zorder=-5)
sage: g += polygon([(1,l), (1,u), (1-w,u)], color='red',fill=True, zorder=-5) 
sage: g += polygon([(1,ll), (1,uu), (1-w,uu)], color='red',fill=True, zorder=-5)
sage: g += polygon([(l,0), (u,0), (l,w)], color='red',fill=True, zorder=-5)
sage: g += polygon([(ll,0), (uu,0), (ll,w)], color='red',fill=True, zorder=-5)
sage: g += polygon([(l,l), (l,u), (u,u), (u,l)], color='red',fill=True, zorder=-5)
sage: g += polygon([(ll,l), (ll,u), (uu,u), (uu,l)], color='red',fill=True, zorder=-5)
sage: g += polygon([(l,ll), (l,uu), (u,uu), (u,ll)], color='red',fill=True, zorder=-5)
sage: g += polygon([(ll,ll), (ll,uu), (uu,uu), (uu,ll)], color='red',fill=True, zorder=-5)
sage: g += polygon([(l,ll-l), (l,uu-l), (u,uu-u), (u,ll-u)], color='red',fill=True, zorder=-5)
sage: g += polygon([(ll-l,l), (uu-l,l), (uu-u,u), (ll-u,u)], color='red',fill=True, zorder=-5)
sage: g += polygon([(ll,1+l-ll), (ll,1+u-ll), (uu,1+u-uu), (uu,1+l-uu)], color='red',fill=True, zorder=-5)
sage: g += polygon([(1+l-ll,ll), (1+u-ll,ll), (1+u-uu,uu), (1+l-uu,uu)], color='red',fill=True, zorder=-5)
sage: tk=[0,l,u,ll,uu,f,1];
sage: tkf=[0,"$l$","$u$","$f-l$","$f-u$","$f$",1]
sage: g.show(ticks=[tk,tk], tick_formatter=[tkf, tkf], show_legend=False)
sage: g.save('2d_crazy_nf.png',ticks=[tk,tk], tick_formatter=[tkf, tkf], show_legend=False)

######  algo-paper ####
sage: h1 = equiv7_example_1()
sage: h = minimal_no_covered_interval()

### copy some plot formatting from graphics_for-survey ###...
sage: extremality_test(h1, show_plots="one-sided-discontinuous-%s.png", show_all
....: _perturbations=True)
sage: extremality_test(h, show_plots="two-sided-discontinuous-%s.png", show_all_
....: perturbations=True)

# crazy perturbation exists for h.
sage: t1, t2 = nice_field_values([1, sqrt(2)])
sage: rnf = t1.parent().fraction_field()

# sage: h = FastPiecewise([singleton_piece(0, 0),open_piece((0,t1/2),(t1/2, t1/2)), singleton_piece(t1/2, t1), open_piece((t1/2,t1/2),(t1,t1/2)), singleton_piece(t1,0)])
# sage: h1 = FastPiecewise([singleton_piece(0, 0),open_piece((0, t1/2), (t1/2, t1/2)),closed_piece((t1/2, t1),(t1,0))])

sage: generators = [t1, t2]
sage: pwl = piecewise_function_from_breakpoints_and_slopes([rnf(0),rnf(1)],[rnf(0)])
sage: crazy_piece_1 = CrazyPiece((rnf(0), rnf(1/4)), generators, [(rnf(1/8), rnf(1)), (rnf(t2/8), rnf(-1))])
sage: crazy_piece_2 = CrazyPiece((rnf(1/4), rnf(1/2)), generators, [(rnf(1/2-t2/8), rnf(1)), (rnf(1/2-1/8), rnf(-1))])
sage: cp = PiecewiseCrazyFunction(pwl, [crazy_piece_1, crazy_piece_2])


sage: find_epsilon_for_crazy_perturbation(h, cp)
INFO: 2017-04-19 16:03:35,927 Rational case.
0.166666666666667

sage: find_epsilon_for_crazy_perturbation(h1, cp)
0

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

########## mip 2017 slides ###########
igprainbow=igp.rainbow
def dark_rainbow(num):
    return ['darkblue', 'darkgreen', 'firebrick', 'darkcyan', 'darkmagenta', 'darkgrey', 'royalblue'][:num]
igp.rainbow=dark_rainbow

h = drlm_not_extreme_1()
extremality_test(h)
hb = h+ h._perturbations[0]*7/20
hr = h -h._perturbations[0]*7/20
g = plot(hb, color='blue')+plot(hr,color='red')+plot(h,color='black',thickness=2)
xx =[i/7 for i in range(8)]; yy = [0,1]
ticks=[xx,yy];tick_formatter=[["$%s$" % latex(x) for x in xx],["$%s$" % latex(x) for x in yy]]
g.fontsize(10)
g.show(ticks=ticks, tick_formatter=tick_formatter,aspect_ratio=1/6, show_legend=False, figsize=4)
g.save("drlm_not_extreme_pert720.pdf",ticks=ticks, tick_formatter=tick_formatter,aspect_ratio=1/6, show_legend=False, figsize=4)


h = automorphism(drlm_backward_3_slope(1/7, 1/7+196/2017))
xx = h.end_points()
ticks=[xx,yy];tick_formatter=[["$%s$" % latex(x) for x in xx],["$%s$" % latex(x) for x in yy]]
g = plot(h,color='black',thickness=2)
g.fontsize(10)
g.show(ticks=ticks, tick_formatter=tick_formatter,aspect_ratio=1/6, show_legend=False, figsize=4)
g.save("drlm_ext_big_grid.pdf",ticks=ticks, tick_formatter=tick_formatter,aspect_ratio=1/6, show_legend=False, figsize=4)


bhk = bhk_irrational(f=4/5, d1=3/5, d2=1/10, a0=15/100, delta=(1/80, sqrt(2)/80))
xx=[0,1]; yy = []
ticks=[xx,yy];tick_formatter=[["$%s$" % latex(x) for x in xx],["$%s$" % latex(x) for x in yy]]
gbhk = plot_with_colored_slopes(bhk)
xs = bhk.end_points()
a0 = xs[2]; f = xs[-2]; l = xs[8]; u = xs[9]; t1 = xs[4]-xs[2]; t2 = xs[6]-xs[2]
g = gbhk+line([(l,0),(u,0)], thickness=2, color='magenta')+line([(f-l,0),(f-u,0)], thickness=2, color='magenta')
g = g +  polygon2d([(l, 0.1), (u, 0.1),(u, 0.9), (l, 0.9)], color="yellow", alpha=0.2, ticks=[[],[]], zorder=-10) +  polygon2d([(f-l, 0.1), (f-u, 0.1),(f-u, 0.9), (f-l, 0.9)], color="yellow", alpha=0.2, ticks=[[],[]], zorder=-10)
tx1 = text("$t_1=\\frac{1}{80}$", (a0,1/8), horizontal_alignment='left', color='black',fontsize=10)
tx2 = text("$t_2 = \\frac{\\sqrt{2}}{80}$", (xs[4],1/2),horizontal_alignment='left', color='black',fontsize=10)
g = g + tx1 + tx2
d2 = t2-t1; ll = 31/100
ar1 = arrow(path=[[(ll, 3/8),(ll+t1/2, 1/2), (ll+t1, 3/8)]],head=1, color='orange', width=1,arrowsize=1,zorder=2) + arrow(path=[[(ll+t1, 3/8),(ll+t1+t1/2, 1/2), (ll+t1+t1, 3/8)]],head=1, color = 'orange', width=1,arrowsize=1,zorder=2) + arrow(path=[[(ll+2*t1, 3/8),(ll+2*t1-d2/2, 1/4), (ll+2*t1 -d2, 3/8)]],head=1, color= 'orange', width=1,arrowsize=1,zorder=2)
tx3 = text("densely covered", ((l+u)/2,0.1),color='magenta',fontsize=10,vertical_alignment='bottom')
g = g+ar1+tx3
g.save("bhk_densely_covered.pdf",ticks=ticks, tick_formatter=tick_formatter,aspect_ratio=1/8, show_legend=False)


def dark_rainbow(num):
    return ['darkblue', 'firebrick', 'darkgreen', 'darkcyan', 'darkmagenta', 'darkgrey', 'royalblue'][:num]
igp.rainbow=dark_rainbow
h = kzh_minimal_has_only_crazy_perturbation_1()
g = plot_with_colored_slopes(h)
bkpts = h.end_points()
t1 = bkpts[10]-bkpts[6]
t2 = bkpts[13]-bkpts[6]
f = bkpts[37]
ucl = bkpts[17]
ucr = bkpts[18]
generators = [t1, t2]
pwl = piecewise_function_from_breakpoints_and_slopes([0,1],[0])
crazy_piece_1 = CrazyPiece((ucl, ucr), generators, [(ucl, 1), (ucr, -1)])
crazy_piece_2 = CrazyPiece((f-ucr, f-ucl), generators, [(f-ucr, 1), (f-ucl, -1)])
cp = PiecewiseCrazyFunction(pwl, [crazy_piece_1, crazy_piece_2])
tx = text("microperiodic", (f/2,0),color='magenta',fontsize=10,vertical_alignment='bottom')+text("perturbation", (f/2,0),color='magenta',fontsize=10,vertical_alignment='top')
g = g + (cp/10).plot() + tx + line([(0,0),(0,1.05)],color='white',zorder=-30)
xx=[0,ucl,ucr,f-ucr,f-ucl,f,1]
ticks=[xx,yy];tick_formatter=[["","","","","","",""],["$%s$" % latex(x) for x in yy]]
xxlabel = ["$0$","$l$","$u$","$f-u$","$f-l$","$f$","$1$"]
g.fontsize(10)
for i in range(1,6):
    g += text(xxlabel[i], (xx[i],-0.15),color='black',fontsize=10,vertical_alignment='top')
g.show(ticks=ticks, tick_formatter=tick_formatter,aspect_ratio=1/8, show_legend=False)g.show(ticks=ticks, tick_formatter=tick_formatter,aspect_ratio=1/8, show_legend=False)
g.save("has_only_crazy_perturbation.pdf",ticks=ticks, tick_formatter=tick_formatter,aspect_ratio=1/8, show_legend=False)
