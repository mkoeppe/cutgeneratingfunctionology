import igp
from igp import *

destdir = "/Users/mkoeppe/Dropbox/basu-hildebrand-koeppe-papers-for-yuan/algo-paper/slides/"

orig_plot = sage.plot.plot.plot

igp.plot = plot_no_legend
igp.plot_kwds_hook = plot_kwds_hook_no_legend

bb = point((1.1, 0.5), color='white', alpha=0)

for name in ['chen_3_slope_not_extreme', 'not_extreme_1', 'drlm_not_extreme_1', 'hildebrand_discont_3_slope_1']:
    h = eval(name)()

    liftname = ''
    rounds = 0
    while True:

        Gh = (plot_2d_diagram_additive_domain_sans_limits(h, show_function=False) + plot_function_at_borders(h, color='black', thickness=3) + bb)
        Gh.save(destdir + '{}{}-2d_diagram_sans_limits.pdf'.format(name, liftname))

        if extremality_test(h):
            break
        #hl = lift_until_extreme(h)
        hl = lift(h)

        Ghl = (plot_2d_diagram_additive_domain_sans_limits(hl, show_function=False, alpha=0.35) + plot_function_at_borders(hl, color='red'))

        (Gh + Ghl).save(destdir + '{}{}-lifted-2d_diagram_sans_limits.pdf'.format(name, liftname))

        #(Gh + Ghl).show()
        h = hl
        rounds += 1
        liftname = '-lift{}'.format(rounds)


name = 'hildebrand_discont_3_slope_1'
h = eval(name)()
hl = discontinuous_interpolation([0,1/2,5/8,7/8],[0,1,3/4,1/4],[0,1/2,3/4,1/4],[1/2,1,3/4,1/4])
liftname = ''

Gh = (plot_2d_diagram_additive_domain_sans_limits(h, show_function=False) + plot_function_at_borders(h, color='black', thickness=3) + bb)
Ghl = (plot_2d_diagram_additive_domain_sans_limits(hl, show_function=False, alpha=0.35) + plot_function_at_borders(hl, color='red') + bb)
(Gh + Ghl).save(destdir + '{}{}-lifted-2d_diagram_sans_limits.pdf'.format(name, liftname))
Ghl.save(destdir + '{}{}-lift1-2d_diagram_sans_limits.pdf'.format(name, liftname))

Gh += plot_2d_diagram_with_cones(h, show_function=False)
Gh.save(destdir + '{}{}-2d_diagram_plus_limits.pdf'.format(name, liftname))
Ghl += plot_2d_diagram_with_cones(hl, show_function=False)
(Gh + Ghl).save(destdir + '{}{}-lifted-2d_diagram_plus_limits.pdf'.format(name, liftname))
Ghl.save(destdir + '{}{}-lift1-2d_diagram_plus_limits.pdf'.format(name, liftname))


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
