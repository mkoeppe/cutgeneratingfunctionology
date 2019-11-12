load("survey_graphics/graphics_for_algo_paper_init.sage")

try:
    destdir = facets_paper_output_dir  # defined in config.sage
except Exception:
    #destdir = "survey_graphics/algo_paper_graphics/"
    destdir = "/Users/yzh/Dropbox/basu-hildebrand-koeppe-papers-for-yuan/algo-paper/graphics-for-facets-paper/"

########## facet paper ###########

def save_graphics(g, fname, **kwds):
    fname = destdir + fname + ftype
    logging.warn('Saving "{}"'.format(fname))
    all_kwds = copy(kwds)
    all_kwds.update(save_kwds)
    g.save(fname, **all_kwds)

igp.plot_limit_cone_style = 'arrows'

# figures for 1/2 linewidth
igp.show_plots_figsize = 8
paper_plot_kwds['fontsize'] = 20

ftype = ".pdf"

save_kwds = {'figsize': igp.show_plots_figsize}

bb = point((1.1, 0.5), color='white', alpha=0)

name = 'hildebrand_discont_3_slope_1'
h = eval(name)()
hl = discontinuous_facets_paper_example_psi_prime(merge=False)
liftname = ''

Gh = (plot_2d_diagram_additive_domain_sans_limits(h, show_function=False) + plot_function_at_borders(h, color='black', thickness=2) + bb)
save_graphics(Gh, '{}{}-2d_diagram_sans_limits'.format(name, liftname), **paper_plot_kwds)
Ghl = (plot_2d_diagram_additive_domain_sans_limits(hl, show_function=False) + plot_function_at_borders(hl, color='red', thickness=2) + bb)
## (Gh + Ghl).save(destdir + '{}{}-lifted-2d_diagram_sans_limits'.format(name, liftname) + ftype, figsize=igp.show_plots_figsize, **paper_plot_kwds)
save_graphics(Ghl, '{}{}-lift1-2d_diagram_sans_limits'.format(name, liftname), **paper_plot_kwds)

Gh += plot_2d_diagram_with_cones(h, show_function=False)
save_graphics(Gh, '{}{}-2d_diagram_plus_limits'.format(name, liftname), **paper_plot_kwds)
Ghl += plot_2d_diagram_with_cones(hl, show_function=False)
## (Gh + Ghl).save(destdir + '{}{}-lifted-2d_diagram_plus_limits'.format(name, liftname) + ftype, figsize=igp.show_plots_figsize, **paper_plot_kwds)
save_graphics(Ghl, '{}{}-lift1-2d_diagram_plus_limits'.format(name, liftname), **paper_plot_kwds)


####

# figures for full linewidth
igp.show_plots_figsize = 12
paper_plot_kwds['fontsize'] = 16
igp.plot_function_at_borders_kwds = { 'thickness': 2 }
#igp.plot_limit_cone_arrow_distance = 8.0/igp.show_plots_figsize * plot_limit_cone_arrow_distance
igp.plot_limit_cone_arrow_length = 8.0/igp.show_plots_figsize * plot_limit_cone_arrow_length

g = plot_covered_intervals(h, thickness=2, **ticks_keywords(h))
save_graphics(g, "{}-only-function".format(name), aspect_ratio=0.3, **paper_plot_kwds)

ftype = ".png"

h = kzh_minimal_has_only_crazy_perturbation_1()
x = h.end_points()
l=x[17]; u = x[18]; ll = x[19]; uu=x[20]; f=x[37]
assert l + uu == f and u + ll == f
w = u - l
color1 = 'lightblue'
color2 = 'mediumslateblue'
color3 = 'darkblue'

igp.plot_2d_complex_continuous_kwds = igp.plot_2d_complex_discontinuous_kwds = {'color': 'grey', 'alpha': 0.3, 'zorder': -2}   # 'linestyle': 'dotted',  'zorder': -10

g = plot_2d_complex(h)
g_at_borders = plot_function_at_borders(h, covered_components=generate_covered_components(h))
g += g_at_borders
tk = [0, h.a0, h.a1, h.a2, l, u, ll, uu, f-h.a2, f-h.a1, l+u, f, 1];
tkfx=[0,"", "$a_1$", "", "$l$","$u$","$f-u$","","","$f-a_1$","","$f$","$1$"]  # omit f-l (doesn't fit)
tkfy=[0,"$a_0$", "", "$a_2$", "$l$","$u$","$f-u$","$f-l$","$f-a_2$", "", "$f-a_0$","$f$",""]  # omit 1 (doesn't fit)
save_kwds = {'ticks': [tk, tk], 'tick_formatter': [tkfx, tkfy], 'show_legend': False,
             'figsize': igp.show_plots_figsize}
save_kwds.update(paper_plot_kwds)

all_special_faces = set(generate_faces_with_projections_intersecting(h, h.special_intervals, break_symmetry=False, halfopen_cover_only=True))
for F in all_special_faces:
    if F.is_2D():
        if number_of_projections_intersecting(F, h.special_intervals) == 1:
            g += F.plot(zorder=-5, fill_color=color1)
        else:
            g += F.plot(zorder=-5, fill_color=color2)
    else:
        if is_additive_face_sans_limits(h, F):
            g += F.plot(zorder=10, rgbcolor=additive_color, edge_thickness=1)
        else:
            g += F.plot(zorder=10, rgbcolor=color3, edge_thickness=1)
for z in tk:
    g += line([(0, z), (z, 0)], linestyle="dotted", thickness=0.5, zorder=9, color='black')
    g += line([(1, z), (z, 1)], linestyle="dotted", thickness=0.5, zorder=9, color='black')

## Vertices from the published proof in Equi VI:
bkpt = h.end_points()
vs6 = [(h.ucl, h.ucr, bkpt[31], 1, -1, 0),
       (bkpt[10], (h._f - h.ucl) - bkpt[10], (h._f - h.ucl), 0, -1, -1)]
g_vs6 = sum(plot_limit_cone_of_vertex(x, y, epstriple_to_cone((xeps, yeps, zeps)), color='red') for (x, y, z, xeps, yeps, zeps) in vs6)
#g += g_vs6

save_graphics(g, '2d_crazy_nf_with_func')



# figures for a poster-sized foldout
ftype = ".pdf"
save_kwds['figsize'] = igp.show_plots_figsize = 36
igp.plot_limit_cone_arrow_rel_distance = 0

h = kzh_minimal_has_only_crazy_perturbation_1()
h._vertices_used = []     # otherwise not recorded
#g = plot_2d_diagram(h, colorful=True, show_projections=False)
g = plot_2d_diagram_additive_domain_sans_limits(h, show_function=False, edge_thickness=4, vertex_size=45) + plot_2d_diagram_with_cones(h, show_function=False) + g_at_borders
try:
    logging.info("Running facet_test")
    facet_test(h, known_extreme=True)
except NotImplementedError as e:
    logging.info("Exception: {}  (this is normal)".format(e))
    pass
vs = [ v for v in h._vertices_used ]

##### Trying to reproduce the same vertices as in appendix of Equi VI -- but they are not the same
##### check the exact script (...strategic?...)
## h = kzh_minimal_has_only_crazy_perturbation_1()
## components6 = generate_covered_intervals(h)
## igp.generate_symbolic_two_sided_discontinuous_basis_functions = ('slopes', 'jumps')   # Restore classic behavior
## symbolic6 = generate_symbolic(h, components6)
## M6, vs6 = generate_additivity_equations(h, symbolic6, reduce_system=True, return_vertices=True)
## vs6 = [ v for v in vs6 if isinstance(v, tuple) ]   # omit special labels 'f', '1'
## g += sum(plot_limit_cone_of_vertex(x, y, epstriple_to_cone((xeps, yeps, zeps)), color='red') for (x, y, z, xeps, yeps, zeps) in vs6)

# plot black arrows slightly larger and below
g += sum(plot_limit_cone_of_vertex(x, y, epstriple_to_cone((xeps, yeps, zeps)),
                                   color='black', zorder=9,
                                   width=3, vertex_size=30)
         for (x, y, z, xeps, yeps, zeps) in vs)
#g += g_vs6

save_graphics(g, 'kzh_crazy_2d_with_func')

## igp.strategical_covered_components = False
## h = kzh_minimal_has_only_crazy_perturbation_1()
## plot_background = polygon2d([[0,0], [0,1], [1,1], [1,0]], fill=False, color='grey') + g_at_borders
## generate_directed_move_composition_completion(h, max_num_rounds=0, error_if_max_num_rounds_exceeded=False, plot_background=plot_background)
## g = h._completion.plot()
## save_graphics(g, 'kzh_crazy_moves_with_func')

## igp.strategical_covered_components = True
## h = kzh_minimal_has_only_crazy_perturbation_1()
## generate_directed_move_composition_completion(h, max_num_rounds=0, error_if_max_num_rounds_exceeded=False, plot_background=plot_background)
## g = h._completion.plot()
## save_graphics(g, 'kzh_crazy_moves_strategic_with_func')

## h = kzh_minimal_has_only_crazy_perturbation_1()
## generate_directed_move_composition_completion(h, plot_background=plot_background)
## g = h._completion.plot()
## save_graphics(g, 'kzh_crazy_completion_with_func')
