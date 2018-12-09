load("survey_graphics/graphics_for_algo_paper_init.sage")

try:
    destdir = facets_paper_output_dir  # defined in config.sage
except Exception:
    #destdir = "survey_graphics/algo_paper_graphics/"
    destdir = "/Users/yzh/Dropbox/basu-hildebrand-koeppe-papers-for-yuan/algo-paper/graphics-for-facets-paper/"



########## facet paper ###########

### old files
## h = hildebrand_discont_3_slope_1()
## hh = discontinuous_interpolation([0,1/2,5/8,7/8],[0,1,3/4,1/4],[0,1/2,3/4,1/4],[1/2,1,3/4,1/4])
## #h3 = restrict_to_finite_group(h, oversampling=3)
## #extremality_test(h3, True)
## g = plot_2d_diagram_additive_domain_sans_limits(h)
## #g.show(show_legend=False)
## g.save(destdir + '2d_simple_E_pi_extreme_not_facet' + ftype, show_legend=False)
## #Launched png viewer for Graphics object consisting of 98 graphics primitives
## gg = plot_2d_diagram_additive_domain_sans_limits(hh)
## #gg.show(show_legend=False)
## gg.save(destdir + '2d_simple_E_pi_prime_extreme_not_facet' + ftype, show_legend=False)


# figures for 1/2 linewidth
igp.show_plots_figsize = 8
paper_plot_kwds['fontsize'] = 20


bb = point((1.1, 0.5), color='white', alpha=0)

name = 'hildebrand_discont_3_slope_1'
h = eval(name)()
hl = discontinuous_interpolation([0,1/2,5/8,7/8],[0,1,3/4,1/4],[0,1/2,3/4,1/4],[1/2,1,3/4,1/4])
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

### The following has been replaced by the more precise test kzh_minimal_has_only_crazy_perturbation_1_check_subadditivity_slacks()
## fn = kzh_minimal_has_only_crazy_perturbation_1()
## gen = set(itertools.chain(generate_type_1_vertices(fn, operator.gt), \
##                                 generate_type_2_vertices(fn, operator.gt)))
## dp = [delta_pi_general(fn, x, y, (xeps, yeps, zeps)) for (x, y, z, xeps, yeps, zeps) in gen]
## igp.show_RNFElement_by_embedding=False
## assert min(dp) == 19/23998


# figures for full linewidth
igp.show_plots_figsize = 12
paper_plot_kwds['fontsize'] = 16

h = kzh_minimal_has_only_crazy_perturbation_1()
x = h.end_points()
l=x[17]; u = x[18]; ll = x[19]; uu=x[20]; f=x[37]
assert l + uu == f and u + ll == f
w = u - l
color1 = 'lightblue'
color2 = 'blue'
g = plot_2d_complex(h)
g += plot_function_at_borders(h, covered_components=generate_covered_components(h), thickness=2)
g += polygon([(0, l), (0, u), (1, u), (1, l)], color=color1,fill=True, zorder=-5)
g += polygon([(0, ll), (0, uu), (1, uu), (1, ll)], color=color1,fill=True, zorder=-5)
g += polygon([(ll,0), (uu,0), (uu,1), (ll,1)], color=color1,fill=True, zorder=-5)
g += polygon([(l,0), (u,0), (u,1), (l,1)], color=color1,fill=True, zorder=-5)
g += polygon([(l,0), (u,0), (0,u), (0,l)], color=color1,fill=True, zorder=-5)
g += polygon([(ll,0), (uu,0), (0,uu), (0,ll)], color=color1,fill=True, zorder=-5)
g += polygon([(l,1), (u,1), (1,u), (1,l)], color=color1,fill=True, zorder=-5)
g += polygon([(ll,1), (uu,1), (1,uu), (1,ll)], color=color1,fill=True, zorder=-5)
g += polygon([(l,1), (u,1), (u,1-w)], color=color2,fill=True, zorder=-5)
g += polygon([(ll,1), (uu,1), (uu,1-w)], color=color2,fill=True, zorder=-5)
g += polygon([(0,l), (0,u), (w,l)], color=color2,fill=True, zorder=-5)
g += polygon([(0,ll), (0,uu), (w,ll)], color=color2,fill=True, zorder=-5)
g += polygon([(1,l), (1,u), (1-w,u)], color=color2,fill=True, zorder=-5)
g += polygon([(1,ll), (1,uu), (1-w,uu)], color=color2,fill=True, zorder=-5)
g += polygon([(l,0), (u,0), (l,w)], color=color2,fill=True, zorder=-5)
g += polygon([(ll,0), (uu,0), (ll,w)], color=color2,fill=True, zorder=-5)
g += polygon([(l,l), (l,u), (u,u), (u,l)], color=color2,fill=True, zorder=-5)
g += polygon([(ll,l), (ll,u), (uu,u), (uu,l)], color=color2,fill=True, zorder=-5)
g += polygon([(l,ll), (l,uu), (u,uu), (u,ll)], color=color2,fill=True, zorder=-5)
g += polygon([(ll,ll), (ll,uu), (uu,uu), (uu,ll)], color=color2,fill=True, zorder=-5)
g += polygon([(l,ll-l), (l,uu-l), (u,uu-u), (u,ll-u)], color=color2,fill=True, zorder=-5)
g += polygon([(ll-l,l), (uu-l,l), (uu-u,u), (ll-u,u)], color=color2,fill=True, zorder=-5)
g += polygon([(ll,1+l-ll), (ll,1+u-ll), (uu,1+u-uu), (uu,1+l-uu)], color=color2,fill=True, zorder=-5)
g += polygon([(1+l-ll,ll), (1+u-ll,ll), (1+u-uu,uu), (1+l-uu,uu)], color=color2,fill=True, zorder=-5)
tk=[0,l,u,ll,uu,f,1];
tkfx=[0,"$l$","$u$","$f-u$","","$f$","$1$"]  # omit f-l (doesn't fit)
tkfy=[0,"$l$","$u$","$f-u$","$f-l$","$f$",""]  # omit 1 (doesn't fit)
#g.show(ticks=[tk,tk], tick_formatter=[tkf, tkf], show_legend=False)
g.save(destdir + '2d_crazy_nf_with_func' + ftype,ticks=[tk,tk], tick_formatter=[tkfx, tkfy], show_legend=False, figsize=igp.show_plots_figsize, **paper_plot_kwds)
