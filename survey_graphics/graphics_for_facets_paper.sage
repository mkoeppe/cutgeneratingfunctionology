load("survey_graphics/graphics_for_algo_paper_init.sage")

try:
    destdir = facets_paper_output_dir  # defined in config.sage
except Exception:
    #destdir = "survey_graphics/algo_paper_graphics/"
    destdir = "/Users/yzh/Dropbox/basu-hildebrand-koeppe-papers-for-yuan/algo-paper/graphics-for-facets-paper/"

#ftype = ".pdf"

########## facet paper ###########

# figures for 1/2 linewidth
igp.show_plots_figsize = 8
paper_plot_kwds['fontsize'] = 20


bb = point((1.1, 0.5), color='white', alpha=0)

name = 'hildebrand_discont_3_slope_1'
h = eval(name)()
hl = discontinuous_facets_paper_example_psi_prime()
liftname = ''

Gh = (plot_2d_diagram_additive_domain_sans_limits(h, show_function=False) + plot_function_at_borders(h, color='black', thickness=2) + bb)
Gh.save(destdir + '{}{}-2d_diagram_sans_limits'.format(name, liftname) + ftype, figsize=igp.show_plots_figsize, **paper_plot_kwds)
Ghl = (plot_2d_diagram_additive_domain_sans_limits(hl, show_function=False) + plot_function_at_borders(hl, color='red', thickness=2) + bb)
## (Gh + Ghl).save(destdir + '{}{}-lifted-2d_diagram_sans_limits'.format(name, liftname) + ftype, figsize=igp.show_plots_figsize, **paper_plot_kwds)
Ghl.save(destdir + '{}{}-lift1-2d_diagram_sans_limits'.format(name, liftname) + ftype, figsize=igp.show_plots_figsize, **paper_plot_kwds)

Gh += plot_2d_diagram_with_cones(h, show_function=False)
Gh.save(destdir + '{}{}-2d_diagram_plus_limits'.format(name, liftname) + ftype, figsize=igp.show_plots_figsize, **paper_plot_kwds)
Ghl += plot_2d_diagram_with_cones(hl, show_function=False)
## (Gh + Ghl).save(destdir + '{}{}-lifted-2d_diagram_plus_limits'.format(name, liftname) + ftype, figsize=igp.show_plots_figsize, **paper_plot_kwds)
Ghl.save(destdir + '{}{}-lift1-2d_diagram_plus_limits'.format(name, liftname) + ftype, figsize=igp.show_plots_figsize, **paper_plot_kwds)


####

# figures for full linewidth
igp.show_plots_figsize = 12
paper_plot_kwds['fontsize'] = 16

h = kzh_minimal_has_only_crazy_perturbation_1()
x = h.end_points()
l=x[17]; u = x[18]; ll = x[19]; uu=x[20]; f=x[37]
assert l + uu == f and u + ll == f
w = u - l
color1 = 'lightblue'
color2 = 'mediumslateblue'
color3 = 'darkblue'
tk=[0,l,u,ll,uu,f,1];
tkfx=[0,"$l$","$u$","$f-u$","","$f$","$1$"]  # omit f-l (doesn't fit)
tkfy=[0,"$l$","$u$","$f-u$","$f-l$","$f$",""]  # omit 1 (doesn't fit)
#g.show(ticks=[tk,tk], tick_formatter=[tkf, tkf], show_legend=False)

all_special_faces = set(generate_faces_with_projections_intersecting(h, h.special_intervals, break_symmetry=False))
for F in all_special_faces:
    if F.is_2D():
        if number_of_projections_intersecting(F, h.special_intervals) == 1:
            g += F.plot(zorder=-5, fill_color=color1)
        else:
            g += F.plot(zorder=-5, fill_color=color2)
    else:
        g += F.plot(zorder=0, rgbcolor=color3, edge_thickness=1)

g.save(destdir + '2d_crazy_nf_with_func' + ftype,ticks=[tk,tk], tick_formatter=[tkfx, tkfy], show_legend=False, figsize=igp.show_plots_figsize, **paper_plot_kwds)
