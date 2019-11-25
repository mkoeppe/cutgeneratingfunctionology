#############################
## Paper: New computer-based search strategies for extreme functions
##        of the Gomory--Johnson infinite group problem
## in:    basu-hildebrand-koeppe-papers/yuan/mip-notes/mip-notes.tex
## Bib:   koeppe-zhou:extreme-search
## Math. Prog. Comp. DOI 10.1007/s12532-016-0115-9
## ArXiv: 1506.00017
#############################

import cutgeneratingfunctionology.igp as igp
from cutgeneratingfunctionology.igp import *

destdir = "survey_graphics/mip_notes_graphics/"

def dark_rainbow(num):
    return ['darkblue', 'darkgreen', 'firebrick', 'darkcyan', 'darkmagenta', 'darkgrey', 'royalblue'][:num]

orig_rainbow = sage.plot.colors.rainbow

def plotfig1():
    igp.rainbow = dark_rainbow
    h1 = gj_2_slope(3/5, 1/2)
    h1r=restrict_to_finite_group(h1,oversampling=1)
    g1_inf = plot_with_colored_slopes(h1)
    g1_fi= plot_covered_intervals(h1r,uncovered_color='red')
    g1_inf.save(destdir+"gj2s_interpolation.pdf", figsize=3,show_legend=False)
    g1_fi.save(destdir+"gj2s_restriction.pdf", figsize=3,show_legend=False)
    igp.rainbow = orig_rainbow

def plotfig2():
    h = gj_forward_3_slope(3/4,2/3,2/2);
    #h = piecewise_function_from_breakpoints_and_values([0,1/8,2/8,3/8,4/8,5/8,6/8,7/8,1],[0, 5/6, 2/6, 3/6, 4/6, 1/6, 6/6, 3/6, 0])
    #g = plot_with_colored_slopes(h)
    g = plot(h)
    txfh = [r'$\frac{%s}{8}$' % x for x in range(9)];
    tyfh = [r'$\frac{%s}{6}$' % x for x in range(7)];

    g += text(r"$\pi_0$", (0/8,0/6), axis_coords=False, vertical_alignment='bottom',horizontal_alignment='left')
    g += text(r"$\pi_1$", (1/8, 5/6), axis_coords=False, vertical_alignment='bottom',horizontal_alignment='center')
    g += text(r"$\pi_2$", (2/8, 2/6), axis_coords=False, vertical_alignment='top',horizontal_alignment='center')
    g += text(r"$\pi_3$", (3/8, 3/6), axis_coords=False, vertical_alignment='top',horizontal_alignment='left')
    g += text(r"$\pi_4$", (4/8, 4/6), axis_coords=False, vertical_alignment='bottom',horizontal_alignment='center')
    g += text(r"$\pi_5$", (5/8, 1/6), axis_coords=False, vertical_alignment='top',horizontal_alignment='center')
    g += text(r"$\pi_6$", (6/8, 6/6), axis_coords=False, vertical_alignment='bottom',horizontal_alignment='center')
    g += text(r"$\pi_7$", (7/8,3/6), axis_coords=False, vertical_alignment='bottom',horizontal_alignment='left')
    g += text(r"$\pi_8$", (8/8,0/6), axis_coords=False, vertical_alignment='bottom',horizontal_alignment='left') 

    s2 = arrow((1/8,5/6),(1/8,2/6), width=1,arrowsize=2,color='black')+line([(1/8,2/6),(2/8,2/6)],color='black')+ text("$s_2$", (1/8,7/12), axis_coords=False, vertical_alignment='center',horizontal_alignment='left', color='black')
    s4 = arrow((3/8,3/6),(3/8,4/6), width=1,arrowsize=2,color='black')+line([(3/8,4/6),(4/8,4/6)],color='black')+ text("$s_4$", (3/8,7/12), axis_coords=False, vertical_alignment='center',horizontal_alignment='right', color='black')
    s5 = arrow((4/8,4/6),(4/8,1/6), width=1,arrowsize=2,color='black')+line([(4/8,1/6),(5/8,1/6)],color='black')+ text("$s_5$", (4/8,5/12), axis_coords=False, vertical_alignment='center',horizontal_alignment='right', color='black')
    s6 = arrow((5/8,1/6),(5/8,6/6), width=1,arrowsize=2,color='black')+line([(5/8,6/6),(6/8,6/6)],color='black')+ text("$s_6$", (5/8,7/12), axis_coords=False, vertical_alignment='center',horizontal_alignment='right', color='black')

    (g+s2+s4+s5+s6).save(destdir+"q_v_grid.pdf",figsize=4, show_legend=False, ticks=[[x/8 for x in range(9)], [x/6 for x in range(7)]], tick_formatter=[txfh, tyfh],  gridlines=True)

def plotfig3():
    h28=kzh_28_slope_1()
    g28=plot_with_colored_slopes(h28)
    g28.save(destdir+"kzh_28_slope.pdf",aspect_ratio=0.3, show_legend=False, gridlines=False, ticks=[[0,1/2,1],[0,1]], tick_formatter=[["$%s$" % latex(x) for x in [0,1/2,1]],["$%s$" % latex(x) for x in [0,1]]])

def plotfig4():
    igp.rainbow = dark_rainbow
    h2q=kzh_2q_example_1()
    tx2q=h2q.end_points()
    txf2q = ["$%s$" % latex(x) for x in tx2q]
    gh2q=plot_covered_intervals(h2q)
    per2q=1/15*piecewise_function_from_breakpoints_and_values([0,30/111,31/111,32/111,33/111,42/111,43/111,44/111,45/111,1],[0,0,1,-1,0,0,1,-1,0,0])
    gper2q=plot(per2q, color='magenta')
    ghper2q = gh2q+gper2q
    g_int1 = polygon2d([(10/37, 0.1), (11/37, 0.1),(11/37, 0.9), (10/37, 0.9)], color="yellow", alpha=0.2, ticks=[[],[]], zorder=-10)
    g_int2 = polygon2d([(14/37, 0.1), (15/37, 0.1),(15/37, 0.9), (14/37, 0.9)], color="yellow", alpha=0.2, ticks=[[],[]], zorder=-10)
    g_arrow_trans = arrow(path=[[(10/37, 0.85),(12/37, 0.85), (14/37, 0.85)]],head = 1, color = 'orange', width=1,arrowsize=2,zorder=-5 )
    g_arrow_reflect = arrow(path = [[(21/74, 101/244), (25/74, 1),(29/74, 143/244)]],head = 2, color = 'orange', width=1,arrowsize=2,zorder=-5 )
    g_mid = line([(25/74, 0), (25/74, 1)],color='orange', linestyle='-.',zorder=-5)
    g_f = line([(25/37, 0), (25/37, 1)],color='black', linestyle=':',zorder=-5)+line([(0, 1), (1, 1)],color='black', linestyle=':',zorder=-5)+line([(1, 0), (1, 1)],color='black', linestyle=':',zorder=-5)
    g2q = ghper2q+g_arrow_reflect+g_mid+g_arrow_trans+g_int1+g_int2+g_f
    g2q.save(destdir+"kzh_2q_move.pdf", aspect_ratio = 0.2, show_legend =False, gridlines=False, ticks=[tx2q, [0,1]], tick_formatter=[txf2q, ["$%s$" % latex(x) for x in [0,1]]])
    igp.rainbow = orig_rainbow

def plotfig7():
     h = not_extreme_1()
     g = plot_2d_diagram(h, colorful = True)
     # not_extreme_1-2d_diagram.pdf is in survey_graphics
     # Don't save again.
     # g.save(destdir+"not_extreme_1-2d_diagram.pdf", figsize=6, show_legend =False, ticks=[[x/10 for x in range(10)] ,[0,1]], tick_formatter=[[r"$\frac{%s}{10}$" % x for x in range(10)], ["$%s$" % latex(x) for x in [0,1]]])

def plotfig9():     
     h = kzh_5_slope_fulldim_1()
     g = plot_covered_intervals(h)
     g.save(destdir+"kzh_5_slope_fulldim.pdf",figsize=3, show_legend=False, ticks=[[0,25/37,1],[0,1]], tick_formatter=[["$%s$" % latex(x) for x in [0,25/37,1]], ["$%s$" % latex(x) for x in [0,1]]],gridlines=True)
     g2d = plot_2d_diagram(h, show_function=False,colorful=True)
     g2d.save(destdir+"kzh_5_slope_fulldim-2d_diagram.pdf",figsize=6, show_legend =False, ticks=[[0,25/37,1],[0,1]], tick_formatter=[["$%s$" % latex(x) for x in [0,25/37,1]], ["$%s$" % latex(x) for x in [0,1]]])

def plotfig11():
    q=58;
    h1 = kzh_6_slope_1()
    g1= plot_2d_diagram(h1, show_function=False, show_projections=False, known_minimal=True, f=None, colorful=True)
    g1.save(destdir+"pattern_s6_q58.pdf", figsize=10, show_legend=False, ticks=[[0.001, (q-4)/6/q, (q-2)/4/q, 1/2, 1],[0,(q-4)/6/q, (q-2)/4/q, 1/2, 1]], tick_formatter=[["$0$",r"$\frac{q-4}{6q}$", r"$\frac{q-2}{4q}$", r"$\frac{1}{2}$","$1$"], ["$0$",r"$\frac{q-4}{6q}$",r"$\frac{q-2}{4q}$", r"$\frac{1}{2}$","$1$"]], fontsize=20, gridlines=False)


    q=166;
    h4=kzh_10_slope_1()
    g4= plot_2d_diagram(h4,show_function=False, show_projections=False, known_minimal=True, f=None, colorful=True)

    c=polygon([(21/q, 31/q),(31/q, 31/q),(31/q, 21/q)], color='black', fill=False, zorder=5)+polygon([(29/q, 29/q),(29/q, 34/q),(34/q, 34/q),(34/q, 29/q)], color='black', fill=False, zorder=5)+polygon([(32/q, 32/q),(32/q, 37/q),(37/q, 37/q),(37/q, 32/q)], color='black', fill=False, zorder=5)+polygon([(35/q, 35/q),(35/q, 40/q),(40/q, 40/q),(40/q, 35/q)], color='black', fill=False, zorder=5)+polygon([(15/q, 34/q),(20/q, 34/q),(25/q, 29/q),(20/q, 29/q)], color='black', fill=False, zorder=5)+polygon([(9/q, 37/q),(14/q, 37/q),(19/q, 32/q),(14/q, 32/q)], color='black', fill=False, zorder=5)+polygon([(3/q, 40/q),(8/q, 40/q),(13/q, 35/q),(8/q, 35/q)], color='black', fill=False, zorder=5)+polygon([(34/q, 15/q),(34/q, 20/q),(29/q, 25/q),(29/q, 20/q)], color='black', fill=False, zorder=5)+polygon([(37/q, 9/q),(37/q, 14/q),(32/q, 19/q),(32/q, 14/q)], color='black', fill=False, zorder=5)+polygon([(40/q, 3/q),(40/q, 8/q),(35/q, 13/q),(35/q, 8/q)], color='black', fill=False, zorder=5);

    c2=polygon([(21/q, 31/q),(31/q, 31/q),(31/q, 21/q)], color='black', fill=False, thickness=2, zorder=5)+polygon([(29/q, 29/q),(29/q, 34/q),(34/q, 34/q),(34/q, 29/q)], color='black', fill=False, thickness=2,zorder=5)+polygon([(32/q, 32/q),(32/q, 37/q),(37/q, 37/q),(37/q, 32/q)], color='black', fill=False, thickness=2, zorder=5)+polygon([(35/q, 35/q),(35/q, 40/q),(40/q, 40/q),(40/q, 35/q)], color='black', fill=False, thickness=2, zorder=5)+polygon([(15/q, 34/q),(20/q, 34/q),(25/q, 29/q),(20/q, 29/q)], color='black', fill=False, thickness=2, zorder=5)+polygon([(9/q, 37/q),(14/q, 37/q),(19/q, 32/q),(14/q, 32/q)], color='black', fill=False, thickness=2, zorder=5)+polygon([(3/q, 40/q),(8/q, 40/q),(13/q, 35/q),(8/q, 35/q)], color='black', fill=False, thickness=2, zorder=5)+polygon([(34/q, 15/q),(34/q, 20/q),(29/q, 25/q),(29/q, 20/q)], color='black', fill=False, thickness=2, zorder=5)+polygon([(37/q, 9/q),(37/q, 14/q),(32/q, 19/q),(32/q, 14/q)], color='black', fill=False, thickness=2, zorder=5)+polygon([(40/q, 3/q),(40/q, 8/q),(35/q, 13/q),(35/q, 8/q)], color='black', fill=False, thickness=2, zorder=5);

    (g4+c2).save(destdir+"pattern_s10_q166.pdf", figsize=10, show_legend=False, ticks=[[0.001, (q-4)/6/q, (q-2)/4/q, 1/2, 1],[0,(q-4)/6/q, (q-2)/4/q, 1/2, 1]], tick_formatter=[["$0$",r"$\frac{q-4}{6q}$", r"$\frac{q-2}{4q}$", r"$\frac{1}{2}$","$1$"], ["$0$",r"$\frac{q-4}{6q}$",r"$\frac{q-2}{4q}$", r"$\frac{1}{2}$","$1$"]], fontsize=20, gridlines=False, xmin=0, xmax=0.25, ymin=0, ymax=0.25)


### MAIN ###
plotfig1()
plotfig2()
plotfig3()
plotfig4()
plotfig7()
plotfig9()
#plotfig10()
plotfig11()

### THE TEX FILE ###
os.system("cd %s && (pdflatex -synctex=1 -src-specials -interaction=nonstopmode mip_notes_graphics)" % (destdir,))

