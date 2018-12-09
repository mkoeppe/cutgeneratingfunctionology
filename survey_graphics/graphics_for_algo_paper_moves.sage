############################################################
## Other moves diagrams:
############################################################

load("survey_graphics/graphics_for_algo_paper_init.sage")

igp.show_plots_figsize = 3
paper_plot_kwds['fontsize'] = 10   # restore to sage default

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
g.save(destdir+"moves-diagram-edges" + ftype,ticks=[tkx,tky], tick_formatter=[tkxf, tkyf], figsize=5)

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
c = DirectedMoveCompositionCompletion(fdms, plot_background=gg, pts_of_discontinuity=[])
c.add_backward_moves()
mtkx = [0,a,b,1,r]
mtky = [t2-1,0,a+t2-1,t1, b+t2-1,a+t1,b+t1,1]
mtkxf = ["$0$","$a$","$b$","$1$","$r$"]
mtkyf = ["$t_2-1$","$0$","$a+t_2-1$","$t_1$","$b+t_2-1$","$a+t_1$","$b+t_1$","$1$"]
ticks_graphics = text("", ((a+b)/2,1+0.03), axis_coords=False, vertical_alignment='bottom',horizontal_alignment='center', color='black',ticks=[mtkx,mtky], tick_formatter=[mtkxf, mtkyf])
gg = c.plot(extra_graphics=ticks_graphics)
#gg.show(figsize=5)
gg.save(destdir+"moves-diagram" + ftype, figsize=5)

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
    g.save(destdir + name + ftype, figsize=3)

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

plot_background = polygon2d([[0,0], [0,1], [1,1], [1,0]], fill=False, color='grey', ticks = [[], []])

## reduce_moves_by_components
#############################

def show_reduced_moves_by_components(moves, comp, name, pts_of_discontinuity=[]):
    fname = destdir + name + "-%s" + ftype
    show_plot(plot_background + plot_covered_components_as_rectangles(comp) + sum(m.plot() for m in moves) + plot_points_of_discontinuity_at_borders(pts_of_discontinuity), fname, tag='completion-unreduced')
    c = DirectedMoveCompositionCompletion(moves, comp,
                                          show_plots=fname, plot_background=plot_background,
                                          pts_of_discontinuity=pts_of_discontinuity,
                                          show_zero_perturbation=False)
    c.extend_components_by_continuity()
    c.reduce_moves_by_components()
    #c.add_backward_moves()
    ## FIXME: This should be an option as well.
    show_plot(c.plot(), fname, tag='completion-initial')
    ## if complete:
    ##     c.complete()
    ## # FIXME: Add an option to DirectedMoveCompositionCompletion.  Optional because within extremality test, don't want to show the last one
    ## # because it will be shown again with perturbation stuff.
    ## show_plot(c.plot(), fname, tag='completion-final')

t1 = 2/15
comp = [ [ open_interval(7/15, 12/15) ] ]

# old - is now a no-op in "long normal form" model
tau1 = FunctionalDirectedMove([open_interval(6/15, 11/15)], (1, t1))
moves = [tau1]
name = 'reduces_moves_by_components_ex1'   ######### TYPO -- change reduces to reduce...
# This one is no-op with algo-paper branch
show_reduced_moves_by_components(moves, comp, name)

# (1a) move pokes in a bit from one side. in "long normal form" gets extended to far boundary
tau1a = FunctionalDirectedMove([open_interval(6/15, 8/15)], (1, t1))
moves = [tau1a]
name = 'reduce_moves_by_components_ex1a'
show_reduced_moves_by_components(moves, comp, name)

# (1b) move pokes in a bit from two sides. In "long normal form" gets extended, joined
tau1a = FunctionalDirectedMove([open_interval(6/15, 8/15)], (1, t1))
tau1b = FunctionalDirectedMove([open_interval(9/15, 11/15)], (1, t1))
moves = [tau1a, tau1b]
name = 'reduce_moves_by_components_ex1b'
show_reduced_moves_by_components(moves, comp, name)

# (1c) move completely within square (gets removed)
tau1a = FunctionalDirectedMove([open_interval(7/15, 8/15)], (1, t1))
moves = [tau1a]
name = 'reduce_moves_by_components_ex1c'
show_reduced_moves_by_components(moves, comp, name)

# (2a) adjacent move at discontinuity, no-op.
tau1a = FunctionalDirectedMove([open_interval(6/15, 7/15)], (1, t1))
moves = [tau1a]
name = 'reduce_moves_by_components_ex2a'
show_reduced_moves_by_components(moves, comp, name, pts_of_discontinuity=[7/15])

# (2b) adjacent move at continuity, extend.
name = 'reduce_moves_by_components_ex2b'
show_reduced_moves_by_components(moves, comp, name)

# (2c) two squares joined at corner by a move
comp = [ [ open_interval(7/15, 12/15) ], [ open_interval(5/15, 7/15) ] ]
tau1 = FunctionalDirectedMove([open_interval(5/15, 11/15)], (1, 0))
moves = [tau1]
name = 'reduce_moves_by_components_ex2c'
show_reduced_moves_by_components(moves, comp, name, pts_of_discontinuity=[7/15])

# (2d) likewise with continuity
name = 'reduce_moves_by_components_ex2d'
show_reduced_moves_by_components(moves, comp, name)

# (e) .... Crucial drawback of the proposed half-open model: Graphs do not have
# the full information stored in moves.
#
# ----> 1) Need to decorate ----](-----
# for adjacent, not to be joined moves
#
# ----> 2) A "short normal form" is NOT obtained by 'taking out boxes from lines'.
#          Consider closed secant in an open box:  Reducing it would lead to singletons
#          at the boundary.

## extend_A operation
#####################

# FIXME: Would be nice to be able to bring to normal form with calling complete().

# merging moves. version without continuity
name = 'extend_moves_ex1a'
tau1a = FunctionalDirectedMove([open_interval(6/15, 7/15)], (1, t1))
tau1b = FunctionalDirectedMove([open_interval(7/15, 11/15)], (1, t1))
moves = [tau1a, tau1b]
show_reduced_moves_by_components(moves, [], name, pts_of_discontinuity=[7/15])

# likewise with continuity
name = 'extend_moves_ex1b'
show_reduced_moves_by_components(moves, [], name, pts_of_discontinuity=[])

# (3) merging rectangles at continuity. version without continuity.
name = 'merge_rectangles_ex3a'
moves = []
comp = [ [open_interval(6/15, 7/15), open_interval(7/15, 10/15) ] ]
show_reduced_moves_by_components(moves, comp, name, pts_of_discontinuity=[7/15])

# likewise with continuity
name = 'merge_rectangles_ex3b'
show_reduced_moves_by_components(moves, comp, name, pts_of_discontinuity=[])


## extend_components_by_moves
#############################

def show_extend_components_by_moves(moves, comp, name, pts_of_discontinuity=[]):
    fname = destdir + name + "-%s" + ftype
    c = DirectedMoveCompositionCompletion(moves, comp,
                                          show_plots=fname, plot_background=plot_background,
                                          pts_of_discontinuity=pts_of_discontinuity,
                                          show_zero_perturbation=False)
    c.complete()
    show_plot(c.plot(), fname, tag='completion-final')


t1 = 2/15
tau1a = FunctionalDirectedMove([open_interval(6/15, 7/15)], (1, t1))
tau1b = FunctionalDirectedMove([open_interval(8/15, 11/15)], (1, t1))
comp = [ [ open_interval(7/15, 10/15) ] ]
name = 'extend_components_by_moves_ex1'
show_extend_components_by_moves([tau1a, tau1b], comp, name)

t1 = 2/15
tau1a = FunctionalDirectedMove([open_interval(4/15, 7/15)], (1, t1))
tau1b = FunctionalDirectedMove([open_interval(8/15, 11/15)], (1, t1))
comp = [ [ open_interval(7/15, 10/15) ] ]
name = 'extend_components_by_moves_ex1a'
show_extend_components_by_moves([tau1a, tau1b], comp, name)

t1 = 5/15
tau1a = FunctionalDirectedMove([open_interval(2/15, 3/15)], (1, t1))
tau1b = FunctionalDirectedMove([open_interval(4/15, 5/15)], (1, t1))
comp = [ [ open_interval(3/15, 4/15), (7/15, 10/15) ] ]
name = 'extend_components_by_moves_ex2'
show_extend_components_by_moves([tau1a, tau1b], comp, name)

t1 = 6/15
tau1a = FunctionalDirectedMove([open_interval(0/15, 1/15)], (1, t1))
tau1b = FunctionalDirectedMove([open_interval(4/15, 5/15)], (1, t1))
comp = [ [ open_interval(1/15, 4/15), (7/15, 10/15) ] ]
name = 'extend_components_by_moves_ex3'
show_extend_components_by_moves([tau1a, tau1b], comp, name)

t1 = 5/15
tau1 = FunctionalDirectedMove([open_interval(2/15, 5/15)], (1, t1))
comp = [ [ open_interval(3/15, 4/15) ] ]
name = 'extend_components_by_moves_ex4'
show_extend_components_by_moves([tau1], comp, name)

# adjacent move at discontinuity
t1 = 2/15
tau1a = FunctionalDirectedMove([open_interval(6/15, 7/15)], (1, t1))
comp = [ [ open_interval(7/15, 12/15) ] ]
name = 'extend_components_by_moves_ex5'
show_extend_components_by_moves([tau1a], comp, name, pts_of_discontinuity=[7/15])
