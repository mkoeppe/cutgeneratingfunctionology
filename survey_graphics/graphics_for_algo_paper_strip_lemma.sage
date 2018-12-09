load("survey_graphics/graphics_for_algo_paper_init.sage")

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
c = DirectedMoveCompositionCompletion([m1, m2], show_plots=False, plot_background=plot_background,pts_of_discontinuity=[], show_zero_perturbation=False)
c.add_backward_moves()
g = c.plot()+ line([(l,0),(l,u)], color='black', linestyle=':') +line([(ll,0),(ll,uu)], color='black', linestyle=':') + line([(uu,0),(uu,uu)], color='black', linestyle=':') + line([(u,0),(u,u)], color='black', linestyle=':') + line([(0,t1),(l,l+t1)], color='black', linestyle=':') + line([(0,t2),(ll,ll+t2)], color='black', linestyle=':') + line([(l,l),(u,l)], color='black', linestyle=':') + line([(ll,ll),(uu,ll)], color='black', linestyle=':')+line([(l,u),(u,u)], color='black', linestyle=':') +line([(ll,uu),(uu,uu)], color='black', linestyle=':') + text("$\\tau_1$", (0.5,0.58), axis_coords=False) +text("$\\tau_2$", (0.62,0.88), axis_coords=False) + text("$\\tau_1^{-1}$", (0.58,0.5), axis_coords=False) +text("$\\tau_2^{-1}$", (0.88, 0.62), axis_coords=False)
tkx = [0,l,ll,uu,u,1]
tky = [0,t1,t2,1]
tkxf = ["$0$","$l$","$l'$","$u'$","$u$","$1$"]
tkyf = ["$0$","$t_1$","$t_2$","$1$"]
#g.show(ticks=[tkx,tky], tick_formatter=[tkxf, tkyf], figsize=5)
g.save(destdir+"strip-lemma" + ftype,ticks=[tkx,tky], tick_formatter=[tkxf, tkyf], figsize=5)
c.complete()
g1 = c.plot()+line([(l,0),(l,u)], color='black', linestyle=':') + line([(u,0),(u,u)], color='black', linestyle=':') + line([(0,l),(u,l)], color='black', linestyle=':') +line([(0,u),(u,u)], color='black', linestyle=':')
#g1.show(ticks=[[0,l,u,1],[0,l,u,1]],tick_formatter=[["$0$","$l$","$u$","$1$"], ["$0$","$l$","$u$","$1$"]],figsize=5)
g1.save(destdir+"strip-lemma-dependent" + ftype,ticks=[[0,l,u,1],[0,l,u,1]],tick_formatter=[["$0$","$l$","$u$","$1$"], ["$0$","$l$","$u$","$1$"]],figsize=5)

t2e = 3/15 + 1/100000*sqrt(2)
m2e = FunctionalDirectedMove([open_interval(l2, u2)], (1, t2e))
c = DirectedMoveCompositionCompletion([m1, m2e], show_plots=False, plot_background=plot_background,pts_of_discontinuity=[], show_zero_perturbation=False)
c.complete()
g2 = c.plot() +line([(l,0),(l,u)], color='black', linestyle=':') + line([(u,0),(u,u)], color='black', linestyle=':') + line([(0,l),(u,l)], color='black', linestyle=':') +line([(0,u),(u,u)], color='black', linestyle=':')
g2.save(destdir+"strip-lemma-independent" + ftype,ticks=[[0,l,u,1],[0,l,u,1]],tick_formatter=[["$0$","$l$","$u$","$1$"], ["$0$","$l$","$u$","$1$"]],figsize=5)
igp.show_translations_and_reflections_separately = save_show_translations_and_reflections_separately


## (reflection, reflection) gives no strip lemma

## # l1 = 3/15
## # u1 = 7/15
## # l2 = 5/15
## # u2 = 9/15
## l1 = l2 = 0
## r1 = 17/15
## r2 = 18/15 + 1/100
## u1 = r1
## u2 = r2
## m1 = FunctionalDirectedMove([open_interval(l1, u1)], (-1, r1))
## m2 = FunctionalDirectedMove([open_interval(l2, u2)], (-1, r2))
## c = DirectedMoveCompositionCompletion([m1, m2], show_plots=destdir+'strip-lemma-r-r-dependent-%s.png', plot_background=plot_background,pts_of_discontinuity=[], show_zero_perturbation=False)
## c.complete()
