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
