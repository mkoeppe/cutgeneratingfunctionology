################################
###  Software-paper
################################

import cutgeneratingfunctionology.igp as igp
from cutgeneratingfunctionology.igp import *

destdir = "survey_graphics/software_paper_graphics/"
ftype = ".pdf"

logging.disable(logging.INFO)

def plot_figure_1():
    def rainbow_brg(num):
        return ['blue','red','green'][:num]
    igp.rainbow = rainbow_brg
    h1 = gmic()
    g1 = plot_with_colored_slopes(h1)
    g1.save(destdir+"gmic_plot"+ftype, figsize=4, show_legend=False)
    igp.rainbow = rainbow
    h2 = equiv5_random_discont_1()
    g2 = plot_with_colored_slopes(h2)
    g2.save(destdir+"discont_pwl_plot"+ftype, figsize=4, show_legend=False)

def plot_figure_2():
    face = Face(([0.2, 0.3], [0.75, 0.85], [1, 1.2]))
    g = face.plot()
    g += line([(0.2, 0), (0.2, 1)], color='grey')
    g += line([(0.3, 0), (0.3, 1)], color='grey')
    g += line([(0, 0.75), (1, 0.75)], color='grey')
    g += line([(0, 0.85), (1, 0.85)], color='grey')
    g += line([(0, 1), (1, 0)], color='grey')
    g += line([(0.2, 1), (1, 0.2)], color='grey')
    g += line([(0, 1), (1, 1)], color='black')
    g += line([(1, 0), (1, 1)], color='black')

    kwds = { 'alpha': proj_plot_alpha, 'zorder': -10 }
    IJK_kwds = [ kwds for i in range(3) ]
    for i in range(3):
        IJK_kwds[i]['color'] = proj_plot_colors[i]
        plot_kwds_hook(IJK_kwds[i])
    g += plot_projections_of_one_face(face, IJK_kwds)
    g += line([(1, 0.15), (0.15, 1)], linestyle=':', color='grey')
    for v in face.vertices:
        g += point(v, color='blue', zorder=10)
    for v in [(0.2, 0.75),(0.2, 1),(0.3, 0.9),(0.35, 0.85),(0.45, 0.75)]:
        g += point(v, color='red', zorder=10)
    g.save(destdir+"construct_a_face"+ftype)

def plot_figure_3():
    h = not_minimal_2()
    g = plot_2d_diagram_with_cones(h)
    deco = line([(1/5, 3/5), (1/5, 1+3/10*h(1/5))], color='cyan', thickness=2, zorder=-1)+line([(1/5, 3/5), (-3/10*h(3/5), 3/5)], color='cyan', thickness=2, zorder=-1)+ line([(1/5, 3/5), (4/5, 0)], color='cyan', thickness=2, zorder=-1)+line([(4/5, 0), (4/5, 1+3/10)], color='cyan', thickness=2, zorder=-1)+text(r"$\Delta\pi(x,y)$", (1/5, 3/5), axis_coords=False, vertical_alignment='bottom', horizontal_alignment='left', color='black', fontsize=15)+text(r"$\pi(x)$", (1/5, 1+3/10*h(1/5)+0.02), axis_coords=False, vertical_alignment='bottom', horizontal_alignment='center', color='black', fontsize=15)+text(r"$\pi(y)$", (-3/10*h(3/5)-0.02, 3/5), axis_coords=False, vertical_alignment='center', horizontal_alignment='right', color='black', fontsize=15)+text(r"$\pi(x+y)$", (4/5, 1+3/10+0.02), axis_coords=False, vertical_alignment='bottom', horizontal_alignment='center', color='black', fontsize=15)
    g.save(destdir+"2d_diagram_with_cones_continuous"+ftype, show_legend=False)
    g += deco
    g.save(destdir+"2d_diagram_with_cones_continuous_deco"+ftype, show_legend=False)

    h = equiv5_random_discont_1()
    g = plot_2d_diagram_with_cones(h)
    g.save(destdir+"2d_diagram_with_cones_discontinuous"+ftype, show_legend=False)

def plot_figure_4():
    # cherry-pick [49d8f23]
    h = example7slopecoarse2()
    gcone=plot_2d_diagram_with_cones(h, conesize=100)
    gface=plot_2d_diagram(h)
    gcone.save(destdir+"2d_with_cones_example7slopecoarse2"+ftype, figsize=10, show_legend=False)
    gface.save(destdir+"2d_with_faces_example7slopecoarse2"+ftype, figsize=10, show_legend=False)

def plot_figure_5():
    h = hildebrand_discont_3_slope_1()
    gcone=plot_2d_diagram_with_cones(h)
    gface=plot_2d_diagram(h)
    gcone.save(destdir+"2d_with_cones_discontinuous"+ftype,show_legend=False)
    gface.save(destdir+"2d_with_faces_discontinuous"+ftype,show_legend=False)

def plot_figure_6():
    h = gj_2_slope(f=3/5, lambda_1=1/3)
    def rainbow_fig_6(num):
        return ['#cc00ff','orange','#cc00ff','cyan','red','blue','orange'][:num]
    igp.show_plots_figsize = 6
    igp.rainbow=rainbow_fig_6
    igp.remove_duplicate_or_empty_components = False
    fname=destdir+"gj2slope-%s"+ftype
    generate_covered_components_strategically(h, show_plots=fname) 
    igp.show_plot_figsize = 10
    igp.rainbow = rainbow
    igp.remove_duplicate_or_empty_components = True

def plot_figure_7():
    h = hildebrand_discont_3_slope_1()
    def rainbow_fig_7(num):
        return ['red','red','red','blue','mediumspringgreen','#cc00ff','turquoise'][:num]
    igp.show_plots_figsize = 6
    igp.rainbow=rainbow_fig_7
    igp.remove_duplicate_or_empty_components = False
    fname=destdir+"animation_2d_diagram_disc-%s"+ftype
    generate_covered_components_strategically(h, show_plots=fname) 
    igp.show_plot_figsize = 10
    igp.rainbow = rainbow
    igp.remove_duplicate_or_empty_components = True

def plot_figures_in_table_3():
    h = drlm_backward_3_slope(f=1/12, bkpt=4/12)
    extremality_test(h)
    perturb = h._perturbations[0]
    epsilon_interval = find_epsilon_interval(h, perturb)
    epsilon = min(abs(epsilon_interval[0]), epsilon_interval[1])
    g = plot_perturbation_diagram(h, perturb)
    g.save(destdir+"not-extreme-perturbation-1"+ftype, figsize=4)

def plot_figures_in_table_4():
    def rainbow_brg(num):
        return ['blue', 'red','green'][:num]
    igp.rainbow = rainbow_brg
    h = gmic(f=4/5)
    g = plot_with_colored_slopes(h)
    g.save(destdir+"gmic45"+ftype, figsize=3, show_legend=False)
    hr = restrict_to_finite_group(h)
    gr = plot_with_colored_slopes(hr)
    gr.save(destdir+"gmic45-restricted"+ftype, figsize=3, show_legend=False)
    ha = automorphism(hr, 2)
    ga = plot_with_colored_slopes(ha)
    ga.save(destdir+"gmic45-auto-rest"+ftype, figsize=3, show_legend=False)
    hi = interpolate_to_infinite_group(ha)
    gi = plot_with_colored_slopes(hi)
    gi.save(destdir+"gmic45-auto-interpol"+ftype, figsize=3, show_legend=False)
    igp.rainbow = rainbow

#### MAIN ####
plot_figure_1()
plot_figure_2()
plot_figure_3()
plot_figure_4()
plot_figure_5()
plot_figure_6()
plot_figure_7()
plot_figures_in_table_3()
plot_figures_in_table_4()

os.system("cd %s && (pdflatex -synctex=1 -src-specials -interaction=nonstopmode software_paper_graphics)" % (destdir,)) 
