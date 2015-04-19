from sage.plot.contour_plot import equify
from sage.plot.contour_plot import ContourPlot
from sage.plot.primitive import GraphicPrimitive
from sage.misc.decorators import options, suboptions
from sage.plot.colors import rgbcolor, get_cmap
from sage.misc.misc import xsrange
import operator

@options(plot_points=100, incol='blue', outcol='white', bordercol=None, borderstyle=None, borderwidth=None,frame=False,axes=True, legend_label=None, aspect_ratio=1)
def region_plot_bis(f, xrange, yrange, plot_points, incol, outcol, bordercol, borderstyle, borderwidth, inalpha=0.5, outalpha=0, **options):
    from sage.plot.all import Graphics
    from sage.plot.misc import setup_for_eval_on_grid
    import numpy

    if not isinstance(f, (list, tuple)):
        f = [f]

    f = [equify(g) for g in f]

    g, ranges = setup_for_eval_on_grid(f, [xrange, yrange], plot_points)
    xrange,yrange=[r[:2] for r in ranges]

    xy_data_arrays = numpy.asarray([[[func(x, y) for x in xsrange(*ranges[0], include_endpoint=True)]
                                     for y in xsrange(*ranges[1], include_endpoint=True)]
                                    for func in g],dtype=float)
    xy_data_array=numpy.abs(xy_data_arrays.prod(axis=0))
    # Now we need to set entries to negative iff all
    # functions were negative at that point.
    neg_indices = (xy_data_arrays<0).all(axis=0)
    xy_data_array[neg_indices]=-xy_data_array[neg_indices]

    from matplotlib.colors import ListedColormap
    incol = rgbcolor(incol)
    outcol = rgbcolor(outcol)
    cmap = ListedColormap([incol, outcol])
    cmap.set_over(outcol, alpha=outalpha)
    cmap.set_under(incol, alpha=inalpha)

    g = Graphics()

    # Reset aspect_ratio to 'automatic' in case scale is 'semilog[xy]'.
    # Otherwise matplotlib complains.
    scale = options.get('scale', None)
    if isinstance(scale, (list, tuple)):
        scale = scale[0]
    if scale == 'semilogy' or scale == 'semilogx':
        options['aspect_ratio'] = 'automatic'

    g._set_extra_kwds(Graphics._extract_kwds_for_show(options, ignore=['xmin', 'xmax']))
    g.add_primitive(ContourPlot(xy_data_array, xrange,yrange,
                                dict(contours=[-1e-20, 0, 1e-20], cmap=cmap, fill=True, **options)))

    if bordercol or borderstyle or borderwidth:
        cmap = [rgbcolor(bordercol)] if bordercol else ['black']
        linestyles = [borderstyle] if borderstyle else None
        linewidths = [borderwidth] if borderwidth else None
        g.add_primitive(ContourPlot(xy_data_array, xrange, yrange,
                                    dict(linestyles=linestyles, linewidths=linewidths,
                                         contours=[0], cmap=[bordercol], fill=False, **options)))

    return g

def linesplot(K, fill=False, color="blue", linewidth=1, legend_label=None, xmin=-0.1, xmax=1.1, ymin=-0.1, ymax=1.1):
    x,y = var('x,y')
    leq, lin = simplify_eq_lt_poly_via_ppl(K.get_eq_factor(), K.get_lt_factor())
    if leq:
        print "WARNING: equation list is not empty!"
    if fill:
        g = region_plot([ lhs(x, y) < 0 for lhs in lin ], (x, xmin, xmax), (y, ymin, ymax), incol=color, plot_points=1000)
        return g
    g = Graphics()
    for l in lin:
        g += implicit_plot(l(x, y) == 0, (x, xmin, xmax), (y, ymin, ymax), color=color, linewidth=linewidth)
    g += line([(0,0),(0,1)], color = color, legend_label=legend_label, zorder=-10)
    return g

def regionplot(K, color="blue", alpha=0.5, legend_label=None, xmin=-0.1, xmax=1.1, ymin=-0.1, ymax=1.1, plot_points=1000):
    x,y = var('x,y')
    leq, lin = simplify_eq_lt_poly_via_ppl(K.get_eq_factor(), K.get_lt_factor())
    if leq:
        print "WARNING: equation list is not empty!"
    g = region_plot_bis([ lhs(x, y) < 0 for lhs in lin ], (x, xmin, xmax), (y, ymin, ymax), incol=color, inalpha=alpha, plot_points=plot_points, bordercol=color)
    g += line([(0,0),(0,1)], color = color, legend_label=legend_label, alpha = alpha, zorder=-10)
    return g
    
def regionplot_drlm_backward_3_slope(f_val=1/12-1/1000, b_val=2/12, plot_points=1000):
    K.<f, b> = SymbolicRealNumberField([f_val, b_val])
    h = drlm_backward_3_slope(f, b, field=K, conditioncheck=False)
    reg_cf = regionplot(K, color="orange", alpha=0.5, legend_label="construction", plot_points=plot_points)

    minimality_test(h)
    reg_cfm = regionplot(K, color="green", alpha=0.5, legend_label="min_test", plot_points=plot_points)

    extremality_test(h)
    reg_cfe = regionplot(K, color="blue", alpha=0.5, legend_label="ext_test", plot_points=plot_points)

    K.<f, b> = SymbolicRealNumberField([f_val, b_val])
    h = drlm_backward_3_slope(f, b, field=K, conditioncheck=True)
    reg_ct  = regionplot(K, color="red", alpha=0.5, legend_label="conditioncheck", plot_points=plot_points)

    p = point([K._values],color = "white", size = 10, zorder=5)

    g = reg_cf+reg_ct+reg_cfm+reg_cfe+p
    return g
    
#sage: logging.disable(logging.INFO)
#sage: g1 = regionplot_drlm_backward_3_slope(f_val=1/12-1/30, b_val=2/12)
#sage: g1.save("drlm_backward_3_slope_1.pdf", title="drlm_backward_3_slope(1/12-1/30, 2/12)",legend_loc=7)
#sage: g2= regionplot_drlm_backward_3_slope(f_val=1/12+1/30, b_val=2/12)
#sage: g2.save("drlm_backward_3_slope_2.pdf", title="drlm_backward_3_slope(1/12+1/30, 2/12)",legend_loc=7)
