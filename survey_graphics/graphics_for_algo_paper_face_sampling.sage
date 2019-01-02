## Face sampling

load("survey_graphics/graphics_for_algo_paper_init.sage")


igp.show_plots_figsize = 10
paper_plot_kwds['fontsize'] = 10   # restore to sage default


def sampled_interval(I, num_samples=5):
    a, b = interval_to_endpoints(I)
    delta = (b - a) / num_samples
    return [ a + i * delta for i in range(num_samples+1) ]

def sampled_face(F, num_samples=5, translations=True, reflections=True):
    I, J, K = F.minimal_triple
    result = []
    if translations:
        result = [ Face([[x], J, K]) for x in sampled_interval(I, num_samples) ] \
          + [ Face([I, [y], K]) for y in sampled_interval(J, num_samples) ]
    if reflections:
        result += [ Face([I, J, [z]]) for z in sampled_interval(K, num_samples) ]
    return result

def symmetric_sampled_faces(F, Fprime, num_samples=4, **kwds):
    E_list = [ E for E in sampled_face(F, num_samples, **kwds) if E.is_1D() ]
    E_list += [ E for E in sampled_face(Fprime, num_samples, **kwds) if E.is_1D() ]
    return E_list

def ticks_keywords_for_faces(faces):
    L = []
    for F in faces:
        L += flatten(F.minimal_triple)
    L = unique_list(sorted(L))
    tick_formatter = [ "$%s$" % latex(x) for x in L ]
    return { 'ticks': (L, L), 'tick_formatter': (tick_formatter, tick_formatter) }

def plot_sampled_stuff(F_list, E_list, name):
    fname = destdir + name + "-%s" + ftype
    fname_sampled = destdir + name + "_sampled-%s" + ftype

    background = polygon(((0,0), (0,1), (1,1), (1,0)), color='grey', fill=False, aspect_ratio=1, zorder=-2, **ticks_keywords_for_faces(F_list))
    g = background + plot_projections_of_faces(additive_faces=F_list)
    g += sum(E.plot(edge_thickness=1) for E in E_list)
    g.save(fname_sampled % "2d_diagram", xmin=0, xmax=1, ymin=0, ymax=1, aspect_ratio=1, dpi=200,
       **ticks_keywords_for_faces(F_list))

    background = polygon(((0,0), (0,1), (1,1), (1,0)), color='grey', fill=False, aspect_ratio=1, zorder=-2, **ticks_keywords_for_faces(F_list))
    # FIXME: The ticks keywords don't seem to apply to the translation moves diagram....
    completion = DirectedMoveCompositionCompletion(fdms=[ E.functional_directed_move() for E in E_list ], show_plots=fname_sampled, plot_background=background, pts_of_discontinuity=[], show_zero_perturbation=False)
    #show(completion.plot())
    completion.complete()
    show_plot(completion.plot(), fname_sampled, tag='completion-final')

    # unsampled
    g = background + plot_projections_of_faces(additive_faces=F_list)
    g += sum(F.plot(edge_thickness=1) for F in F_list)
    g.save(fname % "2d_diagram", xmin=0, xmax=1, ymin=0, ymax=1, aspect_ratio=1, dpi=200,
       **ticks_keywords_for_faces(F_list))

    completion = DirectedMoveCompositionCompletion(covered_components=[ F.covered_component() for F in F_list ], show_plots=fname, plot_background=background, pts_of_discontinuity=[], show_zero_perturbation=False)
    completion.complete()
    show_plot(completion.plot(), fname, tag='completion-final')


###############
I, J, K = [2/19, 3/19], [6/19, 7/19], [8/19, 9/19]
F = Face([I, J, K]) # lower triangle
Fprime = Face([J, I, K]) # swapped
F_list = [F, Fprime]
E_list = symmetric_sampled_faces(F, Fprime)
plot_sampled_stuff(F_list, E_list, 'triangle')

# reflections-only version
E_list = symmetric_sampled_faces(F, Fprime, translations=False)
plot_sampled_stuff(F_list, E_list, 'triangle_reflections')

# translations-only version
E_list = symmetric_sampled_faces(F, Fprime, reflections=False)
plot_sampled_stuff(F_list, E_list, 'triangle_translations')


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


# #############

xmax = 0.8
ymax = 0.8

igp.show_plots_figsize = 6

show_kwds = copy(paper_plot_kwds)
show_kwds.update({'ticks': [[], []], 'axes': False, 'frame': True})
show_kwds.update({'xmin': 0, 'xmax': xmax, 'ymin': 0, 'ymax': ymax})

background = Graphics()
## background = polygon(((0,0), (0,ymax), (xmax,ymax), (xmax,0)), color='grey',
##                      fill=False, aspect_ratio=1, zorder=-2, **show_kwds)

num_sample = 7
t = [35/100 + (50-35)/100/(num_sample-1) * i for i in range(num_sample)] # delta = 0.025
open_domains = []
for i in range(num_sample):
    x, y = var('x,y')
    eq1 = 40000*(x-0.2)^2-20000*(x-0.2)*(y-0.7)+10000*(y-0.7)^2-400*(x-0.2)+1600*(y-0.7)+1.0 == 0.0
    eq2 = (x-18/100)^2+(y-6/10)^2==(5/200)^2
    eq = y == x + t[i]
    solns = solve([eq1,eq], x, y)
    a1 = (QQ(solns[0][0].rhs().n(30)))
    b2 = (QQ(solns[1][0].rhs().n(30)))
    solns = solve([eq2,eq], x, y)
    if solns[0][0].rhs().is_real():
        b1 = (QQ(solns[0][0].rhs().n(30)))
        a2 = (QQ(solns[1][0].rhs().n(30)))
        open_domains.append([open_interval(a1,b1), open_interval(a2,b2)])
    else:
        open_domains.append([open_interval(a1,b2)])
fdms_tau = [FunctionalDirectedMove(open_domains[i],(1,t[i])) for i in range(num_sample)]
name = 'open_sets_of_tau'
fname = destdir + name + "-%s" + ftype
fname_sampled = destdir + name + "_sampled-%s" + ftype
x, y = var('x,y')
g = background
g += implicit_plot(40000*(x-0.2)^2-20000*(x-0.2)*(y-0.7)+10000*(y-0.7)^2-400*(x-0.2)+1600*(y-0.7)+1.0, (x, 0,1), (y,0,1),linewidth=1, linestyle=':', **show_kwds)+implicit_plot((x-18/100)^2+(y-6/10)^2==(5/200)^2, (x,0,xmax), (y,0,ymax),linewidth=1, linestyle=':', **show_kwds)
for fdm in fdms_tau:
    g += fdm.plot()
show_plot(g, fname_sampled, tag='', **show_kwds)
completion_tau = DirectedMoveCompositionCompletion(fdms=fdms_tau,show_plots=fname_sampled, plot_background=background, pts_of_discontinuity=[], show_zero_perturbation=False, show_kwds=show_kwds)
completion_tau.complete()
show_plot(completion_tau.plot(), fname_sampled, tag='completion-final', **show_kwds)

num_sample = 11 ;
r = [65/100 + (90-65)/100/(num_sample-1) * i for i in range(num_sample)]
open_domains = []
for i in range(num_sample):
    x, y = var('x,y')
    eq1 = 40000*(x-0.2)^2-20000*(x-0.2)*(y-0.7)+10000*(y-0.7)^2-400*(x-0.2)+1600*(y-0.7)+1.0 == 0.0
    eq2 = (x-18/100)^2+(y-6/10)^2==(5/200)^2
    eq = x + y == r[i]
    solns = solve([eq1,eq], x, y)
    a1 = (QQ(solns[0][0].rhs().n(30)))
    b2 = (QQ(solns[1][0].rhs().n(30)))
    solns = solve([eq2,eq], x, y)
    if solns[0][0].rhs().is_real():
        b1 = (QQ(solns[0][0].rhs().n(30)))
        a2 = (QQ(solns[1][0].rhs().n(30)))
        open_domains.append([open_interval(a1,b1), open_interval(a2,b2)])
    else:
        open_domains.append([open_interval(a1,b2)])
fdms_rho = [FunctionalDirectedMove(open_domains[i],(-1,r[i])) for i in range(num_sample)]
name = 'open_sets_of_rho'
fname = destdir + name + "-%s" + ftype
fname_sampled = destdir + name + "_sampled-%s" + ftype
x, y = var('x,y')
g = background
g += implicit_plot(40000*(x-0.2)^2-20000*(x-0.2)*(y-0.7)+10000*(y-0.7)^2-400*(x-0.2)+1600*(y-0.7)+1.0, (x, 0,xmax), (y,0,ymax),linewidth=1, linestyle=':',fillcolor='red', color='red')+implicit_plot((x-18/100)^2+(y-6/10)^2==(5/200)^2, (x,0,xmax), (y,0,ymax),linewidth=1, linestyle=':',fillcolor='red', color='red', **show_kwds)
for fdm in fdms_rho:
    g += fdm.plot()
show_plot(g, fname_sampled, tag='', **show_kwds)
completion_rho = DirectedMoveCompositionCompletion(fdms=fdms_rho,show_plots=fname_sampled, plot_background=background, pts_of_discontinuity=[], show_zero_perturbation=False, show_kwds=show_kwds)
completion_rho.complete()
show_plot(completion_rho.plot(), fname_sampled, tag='completion-final', **show_kwds)

num_sample = 11 ;
dx = 15/100; dy = -15/100  #dx = 16/100; dy = -16/100 better?
r2 = [18/100+6/10+dx+dy + 26/100/(num_sample-1) * (i-5) for i in range(num_sample)]
open_domains = []
for i in range(num_sample):
    x, y = var('x,y')
    eq1 = 40000*(x-dx-0.2)^2-20000*(x-dx-0.2)*(y-dy-0.7)+10000*(y-dy-0.7)^2-400*(x-dx-0.2)+1600*(y-dy-0.7)+1.0 == 0.0
    eq2 = (x-dx-18/100)^2+(y-dy-6/10)^2==(5/200)^2
    eq = x + y == r2[i]
    solns = solve([eq1,eq], x, y)
    a1 = (QQ(solns[0][0].rhs().n(30)))
    b2 = (QQ(solns[1][0].rhs().n(30)))
    if a1 > b2:
        a1, b2 = b2, a1
    solns = solve([eq2,eq], x, y)
    if solns[0][0].rhs().is_real():
        b1 = (QQ(solns[0][0].rhs().n(30)))
        a2 = (QQ(solns[1][0].rhs().n(30)))
        if b1 > a2:
            b1, a2 = a2, b1 
        open_domains.append([open_interval(a1,b1), open_interval(a2,b2)])
    else:
        open_domains.append([open_interval(a1,b2)])
fdms_rho2 = [FunctionalDirectedMove(open_domains[i],(-1,r2[i])) for i in range(num_sample)]
name = 'open_sets_of_rho2'
fname = destdir + name + "-%s" + ftype
fname_sampled = destdir + name + "_sampled-%s" + ftype
x, y = var('x,y')
g = background
g = implicit_plot(40000*(x-dx-0.2)^2-20000*(x-dx-0.2)*(y-dy-0.7)+10000*(y-dy-0.7)^2-400*(x-dx-0.2)+1600*(y-dy-0.7)+1.0, (x, 0,xmax), (y,0,ymax),linewidth=1, linestyle=':',fillcolor='red', color='red')+implicit_plot((x-dx-18/100)^2+(y-dy-6/10)^2==(5/200)^2, (x,0,xmax), (y,0,ymax),linewidth=1, linestyle=':',fillcolor='red', color='red', **show_kwds)
for fdm in fdms_rho2:
    g += fdm.plot()
show_plot(g, fname_sampled, tag='', **show_kwds)
completion_rho2 = DirectedMoveCompositionCompletion(fdms=fdms_rho2,show_plots=fname_sampled, plot_background=background, pts_of_discontinuity=[], show_zero_perturbation=False, show_kwds=show_kwds)
completion_rho2.complete()
show_plot(completion_rho2.plot(), fname_sampled, tag='completion-final', **show_kwds)

name = 'open_sets_of_moves'
fname = destdir + name + "-%s" + ftype
fname_sampled = destdir + name + "_sampled-%s" + ftype
x, y = var('x,y')
g = background
g += implicit_plot(40000*(x-0.2)^2-20000*(x-0.2)*(y-0.7)+10000*(y-0.7)^2-400*(x-0.2)+1600*(y-0.7)+1.0, (x, 0,xmax), (y,0,ymax),linewidth=1, linestyle=':',fillcolor='purple', color='purple')+implicit_plot((x-18/100)^2+(y-6/10)^2==(5/200)^2, (x,0,1), (y,0,1),linewidth=1, linestyle=':',fillcolor='purple', color='purple', **show_kwds)
fdms = fdms_tau + fdms_rho
for fdm in fdms:
    g += fdm.plot()
show_plot(g, fname_sampled, tag='', **show_kwds)
completion = DirectedMoveCompositionCompletion(fdms=fdms,show_plots=fname_sampled, plot_background=background, pts_of_discontinuity=[], show_zero_perturbation=False, show_kwds=show_kwds)
completion.complete()
show_plot(completion.plot(), fname_sampled, tag='completion-final', **show_kwds)

### Limits imply components
xmax = 1
ymax = 1
show_kwds = copy(paper_plot_kwds)
show_kwds.update({'ticks': [[], []], 'axes': False, 'frame': False})

name = 'lim_tau_dense'
fname = destdir + name + "-%s" + ftype
t = [0, 1/64, 1/32, 1/16]
fdms = [FunctionalDirectedMove([open_interval(1/8+ti,3/8)],(1,1/2-ti)) for ti in t]
background = polygon(((0,0), (0,ymax), (xmax,ymax), (xmax,0)), color='grey', fill=False, aspect_ratio=1, zorder=-2)
completion = DirectedMoveCompositionCompletion(fdms=fdms,show_plots=fname, plot_background=background, pts_of_discontinuity=[], show_zero_perturbation=False, show_kwds=show_kwds)
completion.complete()
show_plot(completion.plot(), fname, tag='completion-final', **show_kwds)

name = 'lim_rho_dense'
fname = destdir + name + "-%s" + ftype
fdms = [FunctionalDirectedMove([open_interval(1/8+ti,3/8)],(-1,1+ti)) for ti in t]
completion = DirectedMoveCompositionCompletion(fdms=fdms,show_plots=fname, plot_background=background, pts_of_discontinuity=[], show_zero_perturbation=False, show_kwds=show_kwds)
completion.complete()
show_plot(completion.plot(), fname, tag='completion-final', **show_kwds)

