
import cutgeneratingfunctionology.igp as igp
from cutgeneratingfunctionology.igp import *

destdir = "survey_graphics/crazy_perturbation_graphics/"
ftype = ".png"

logging.disable(logging.INFO)


# plot_covered_intervals_with_markers branch
h = zhou_two_sided_discontinuous_cannot_assume_any_continuity()
g = plot_covered_intervals(h)
g.show(figsize=4, show_legend=False)
g.save(destdir + "zhou_two_sided_discontinuous_cannot_assume_any_continuity-covered_intervals.png", figsize=4, show_legend=False)
extremality_test(h)
# INFO: 2016-04-18 03:27:22,579 Finding epsilon interval for perturbation... done.  Interval is [-4/9, 4/3]
pert = h._perturbations[0] * (4/9)

default_rainbow = rainbow

def black_rainbow(num):
        return ['black', 'black', 'black'][:num]

igp.rainbow = black_rainbow
g = plot(pert, color='magenta')
g += plot(h+pert, color='blue')
g += plot(h-pert, color='red')
g += plot_covered_intervals(h)
g.show(figsize=4, show_legend=False)
g.save(destdir + "zhou_two_sided_discontinuous_cannot_assume_any_continuity-perturbation-1.png", figsize=4, show_legend=False)
igp.rainbow = default_rainbow

## # bhkish branch
## #graphics for the paper entitled 'crazy example'
## # note bhk_discontinuous is replaced by kzh_discontinuous_bhk_irrational

## sage: import igp
## sage: from cutgeneratingfunctionology.igp import *
## sage: load("bhk.sage")
## sage: delta =(1/199, sqrt(2)/199)
## sage: h = bhk_irrational(delta =(1/199, sqrt(2)/199))
## #sage: g = plot_with_colored_slopes(h)
## sage: g = plot_covered_intervals(h)
## # bhk_irrational
## sage: g.show(figsize=20, aspect_ratio=1/4)
## sage: g.save("bhk_irrational.png", figsize=20, aspect_ratio=1/4, show_legend=False)
## # bhk_irrational_zoom_1
## #sage: g.show(xmin=3/20-1/60, xmax=13/80+1/60, ymin=1/4, ymax=1/3, show_legend=False)
## sage: g.show(xmin=29/200, xmax=34/200, ymin=1/4, ymax=28/100, show_legend=False)
## sage: g.save('bhk_irrational_zoom_1.png',xmin=29/200, xmax=34/200, ymin=1/4, ymax=28/100, show_legend=False)
## # bhk_irrational_2d
## sage: g2d = plot_2d_diagram(h, colorful=True)
## sage: g2d.show(figsize=40)
## sage: g2d.save("bhk_irrational_2d.png", figsize=40, show_legend=False)

## # bhk_irrational_2d_zoom_1
## sage: g2d.show(xmin=3/10, xmax=7/20, ymin=29/200, ymax=34/200, show_legend=False)

## sage: bkpt = h.end_points()
## sage: generators = [bkpt[4]-bkpt[2], bkpt[6]-bkpt[4]]
## sage: pwl = piecewise_function_from_breakpoints_and_slopes([0,1],[0])
## INFO: 2016-04-10 15:00:45,739 Rational case.
## sage: crazy_piece_1 = CrazyPiece((bkpt[8], bkpt[9]), generators, [(bkpt[8], 1), (bkpt[9], -1)])
## sage: crazy_piece_2 = CrazyPiece((bkpt[10], bkpt[11]), generators, [(bkpt[10], 1), (bkpt[11], -1)])
## sage: cp = CrazyPerturbation(pwl, [crazy_piece_1, crazy_piece_2])
## sage: find_epsilon_for_crazy_perturbation(h, cp, True)
## Launched png viewer for Graphics object consisting of 1304 graphics primitives
## 0
## sage: f = 4/5; d1 = 3/5; d2 = 1/10; a0 = 15/100;
## sage: y1 = 0; y2 = 0;
## sage: d3 = f - d1 - d2; c2 = 0;
## sage: bb = 1/16
## sage: rnf = f.parent().fraction_field()
## sage: c3 = -1/(1-f);
## sage: c1 = (1-d2*c2-d3*c3+2*(y1+y2))/d1
## sage: b = (y1 + bb) / (c2 - c3)
## sage: j = b*(c1-c3)
## sage: b_shift_pieces = [singleton_piece(rnf(0), 0), \
## ....:                       open_piece((rnf(0),j), (b, 0)), \
## ....:                       singleton_piece(b, j),
## ....:                       open_piece((b, 0), (f-b, 0)), \
## ....:                       singleton_piece(f-b, -j),
## ....:                       open_piece((f-b, 0), (f, -j)), \
## ....:                       singleton_piece(f,0), \
## ....:                       open_piece((f, -j), (f+b, 0)), \
## ....:                       singleton_piece(f+b, -j),
## ....:                       open_piece((f+b, 0), (1-b, 0)), \
## ....:                       singleton_piece(1-b, j),
## ....:                       open_piece((1-b, 0), (rnf(1), j)), \
## ....:                       singleton_piece(rnf(1), 0) ]
## sage: b_shift = FastPiecewise(b_shift_pieces)
## sage: plot(b_shift)
## Launched png viewer for Graphics object consisting of 49 graphics primitives
## # bhk_irrational_and_disc_b
## sage: gcovered = plot_covered_intervals(h)
## sage: (gcovered + plot(b_shift,color = 'magenta'))
## sage: (gcovered + plot(b_shift,color = 'magenta')).save('bhk_irrational_and_disc_b.png',show_legend=False)

## sage: hb = h + b_shift
## sage: gb = plot_with_colored_slopes(hb)
## sage: gbcovered = plot_covered_intervals(hb)
## # bhk_irrational_disc_b
## sage: gb.show(show_legend=False)
## sage: gbcovered.save('bhk_irrational_disc_b.png', show_legend=False)


## sage: hd = bhk_discontinuous_finite_dim_5(f,d1,d2,a0, delta,b=b)
## sage: t1=delta[0]; t2=delta[0]+delta[1];
## sage: a1 = a0+t1; a2 = a0+t2
## sage: l0 = line([(2/15, h(2/15)), (a0-b, h(a0-b)), (a0, 1/4+j), (a0+b, h(a0+t2+b)), (13/80+1/60, h(13/80+1/60+t2))], color='magenta', linestyle='--', zorder=-1) + line([(a0-b, h(a0-b)), (a0, 1/4), (a0+b, h(a0+t2+b))], color='magenta', linestyle='--', zorder=-1) + line([(a0, 1/4+j), (a0, 1/4)], color='magenta', linestyle=':', zorder=-1)
## sage: l1 = line([(2/15+t1, h(2/15)), (a1-b, h(a0-b)), (a1, 1/4+j), (a1+b, h(a0+t2+b)), (13/80+1/60+t1, h(13/80+1/60+t2))], color='red', linestyle='--', zorder=-1) + line([(a1-b, h(a0-b)), (a1, 1/4), (a1+b, h(a0+t2+b))], color='red', linestyle='--', zorder=-1) + line([(a1, 1/4+j), (a1, 1/4)], color='red', linestyle=':', zorder=-1)
## sage: l2 = line([(2/15+t2, h(2/15)), (a2-b, h(a0-b)), (a2, 1/4+j), (a2+b, h(a0+t2+b)), (13/80+1/60+t2, h(13/80+1/60+t2))], color='orange', linestyle='--', zorder=-1) + line([(a2-b, h(a0-b)), (a2, 1/4), (a2+b, h(a0+t2+b))], color='orange', linestyle='--', zorder=-1) + line([(a2, 1/4+j), (a2, 1/4)], color='orange', linestyle=':', zorder=-1)
## # bhk_disc_b_zoom.png
## sage: gzigzag =  (l0+l1+l2)+ plot(h, color='grey')+ plot(hd, color='black')
## sage: gzigzag.show(xmin=2/15, xmax=13/80+1/60, ymin=1/4, ymax=3/8, show_legend=False, gridlines=None)
## sage: gzigzag.save('bhk_disc_b_zoom.png', xmin=2/15, xmax=13/80+1/60, ymin=1/4, ymax=3/8, show_legend=False, gridlines=None)

## # 'bhk_disc_b_modified.png'
## sage: gd = plot_covered_intervals(hd)
## sage: gd.save('bhk_disc_b_modified.png', figsize=10, show_legend=False, gridlines=None)

## # bhk_disc_too_many_additivities
## sage: find_epsilon_for_crazy_perturbation(hd, cp, True)

## # change in bhk_discontinuous() #return h_pwl, y_shift, b_shift, h_temp, h_below, h_above, move_shift,  (a1-a0, a2-a0), (bkpts[7], bkpts[8])
## sage: h_pwl, y_shift, b_shift, h_temp, h_below, h_above, move_shift,(t1,t2), (ucl, ucr)=bhk_discontinuous()

## # bhk_y_shift_disc_b.png
## sage: hy = h_pwl+y_shift
## sage: hp = b_shift - h_below + h_above + move_shift
## sage: gyp = plot(hy)+plot(hp, color='magenta')
## sage: gyp.save('bhk_y_shift_disc_b.png')

## # bhk_disc_not_min.png
## sage: h_org = hy+hp
## sage: ghorg =  plot_with_colored_slopes(h_org)
## sage: ghorg.save('bhk_disc_not_min.png', show_legend=False, gridlines=None)

## # bhk_disc_min.png
## #sage: h_org, (t1,t2), (ucl, ucr) = bhk_discontinuous()
## INFO: 2016-04-11 01:09:41,099 Coerced into real number field: Real Number Field in `a` with defining polynomial y^2 - 2
## sage: finite_dimensional_extremality_test(h_org)
## False
## sage: h_pert = h_org._perturbations[0]
## sage: ghorgpert = plot(h_pert, color='magenta')
## sage: (ghorg+ghorgpert).save('bhk_disc_not_min_pert.png', show_legend=False, gridlines=None)
## sage: hmin = lift(h_org)
## sage: ghmin = plot_with_colored_slopes(hmin)
## sage: ghmin.save('bhk_disc_min.png', show_legend=False, gridlines=None)


## sage: minimality_test(hmin)
## INFO: 2016-04-11 01:23:42,962 pi(0) = 0
## INFO: 2016-04-11 01:23:42,963 pi is subadditive.
## INFO: 2016-04-11 01:23:42,965 pi is symmetric.
## INFO: 2016-04-11 01:23:42,965 Thus pi is minimal.
## True
## sage: generate_uncovered_intervals(hmin)
## INFO: 2016-04-11 01:24:02,695 Computing maximal additive faces...
## INFO: 2016-04-11 01:24:06,096 Computing maximal additive faces... done
## INFO: 2016-04-11 01:24:06,485 Completing 33 functional directed moves and 2 covered components...
## INFO: 2016-04-11 01:24:06,851 Completing 29 functional directed moves and 2 covered components...
## INFO: 2016-04-11 01:24:06,939 New dense move from strip lemma: [<Int(0.2737500000000000?, 0.3362500000000000?)>, <Int(0.4637500000000000?, 0.526250000000000?)>]
## INFO: 2016-04-11 01:24:06,960 Completing 0 functional directed moves and 3 covered components...
## INFO: 2016-04-11 01:24:06,961 Completion finished.  Found 0 directed moves and 3 covered components.
## []


## sage: y1 = (hmin(ucl)-hmin.limit(ucl,1)).list()[0]
## sage: y2 = (hmin.limit(ucr,-1)-hmin(ucr)).list()[0]
## sage: c2 = ((hmin.limit(ucr,-1) - hmin.limit(ucl,1))/(ucr-ucl)).list()[0]
## sage: c3 = -5
## sage: b = hmin.end_points()[1].list()[0]
## sage: bb = b*(c2-c3)-y1
## sage: hmin = bhk_discontinuous(f=4/5, d1=3/5, d2=5/40, a0=19/100, delta_ratio=sqrt(2)/3, bb=19/23998, c2=5/11999, y1=185/1846, y2=240/11999, field=None)



## sage: finite_dimensional_extremality_test(hmin,show_all_perturbations=True)
## INFO: Finite dimensional test: Solution space has dimension 5
## False
## sage: len(hmin._perturbations)
## 5
## sage: colours = rainbow(5)
## sage: bkpt = h.end_points()
## sage: generators = [bkpt[4]-bkpt[2], bkpt[6]-bkpt[4]]
## sage: pwl = piecewise_function_from_breakpoints_and_slopes([0,1],[0])
## INFO: 2016-04-10 15:00:45,739 Rational case.
## sage: crazy_piece_1 = CrazyPiece((bkpt[8], bkpt[9]), generators, [(bkpt[8], 1), (bkpt[9], -1)])
## sage: crazy_piece_2 = CrazyPiece((bkpt[10], bkpt[11]), generators, [(bkpt[10], 1), (bkpt[11], -1)])
## sage: cp = CrazyPerturbation(pwl, [crazy_piece_1, crazy_piece_2])
## sage: find_epsilon_for_crazy_perturbation(h, cp, True)
## Launched png viewer for Graphics object consisting of 1304 graphics primitives
## 0
## sage: f = 4/5; d1 = 3/5; d2 = 1/10; a0 = 15/100;
## sage: y1 = 0; y2 = 0;
## sage: d3 = f - d1 - d2; c2 = 0;
## sage: bb = 1/16
## sage: rnf = f.parent().fraction_field()
## sage: c3 = -1/(1-f);
## sage: c1 = (1-d2*c2-d3*c3+2*(y1+y2))/d1
## sage: b = (y1 + bb) / (c2 - c3)
## sage: j = b*(c1-c3)
## sage: b_shift_pieces = [singleton_piece(rnf(0), 0), \
## ....:                       open_piece((rnf(0),j), (b, 0)), \
## ....:                       singleton_piece(b, j),
## ....:                       open_piece((b, 0), (f-b, 0)), \
## ....:                       singleton_piece(f-b, -j),
## ....:                       open_piece((f-b, 0), (f, -j)), \
## ....:                       singleton_piece(f,0), \
## ....:                       open_piece((f, -j), (f+b, 0)), \
## ....:                       singleton_piece(f+b, -j),
## ....:                       open_piece((f+b, 0), (1-b, 0)), \
## ....:                       singleton_piece(1-b, j),
## ....:                       open_piece((1-b, 0), (rnf(1), j)), \
## ....:                       singleton_piece(rnf(1), 0) ]
## sage: b_shift = FastPiecewise(b_shift_pieces)
## sage: plot(b_shift)
## Launched png viewer for Graphics object consisting of 49 graphics primitives
## # bhk_irrational_and_disc_b
## sage: gcovered = plot_covered_intervals(h)
## sage: (gcovered + plot(b_shift,color = 'magenta'))
## sage: (gcovered + plot(b_shift,color = 'magenta')).save('bhk_irrational_and_disc_b.png',show_legend=False)

## sage: hb = h + b_shift
## sage: gb = plot_with_colored_slopes(hb)
## sage: gbcovered = plot_covered_intervals(hb)
## # bhk_irrational_disc_b
## sage: gb.show(show_legend=False)
## sage: gbcovered.save('bhk_irrational_disc_b.png', show_legend=False)


## sage: hd = bhk_discontinuous_finite_dim_5(f,d1,d2,a0, delta,b=b)
## sage: t1=delta[0]; t2=delta[0]+delta[1];
## sage: a1 = a0+t1; a2 = a0+t2
## sage: l0 = line([(2/15, h(2/15)), (a0-b, h(a0-b)), (a0, 1/4+j), (a0+b, h(a0+t2+b)), (13/80+1/60, h(13/80+1/60+t2))], color='magenta', linestyle='--', zorder=-1) + line([(a0-b, h(a0-b)), (a0, 1/4), (a0+b, h(a0+t2+b))], color='magenta', linestyle='--', zorder=-1) + line([(a0, 1/4+j), (a0, 1/4)], color='magenta', linestyle=':', zorder=-1)
## sage: l1 = line([(2/15+t1, h(2/15)), (a1-b, h(a0-b)), (a1, 1/4+j), (a1+b, h(a0+t2+b)), (13/80+1/60+t1, h(13/80+1/60+t2))], color='red', linestyle='--', zorder=-1) + line([(a1-b, h(a0-b)), (a1, 1/4), (a1+b, h(a0+t2+b))], color='red', linestyle='--', zorder=-1) + line([(a1, 1/4+j), (a1, 1/4)], color='red', linestyle=':', zorder=-1)
## sage: l2 = line([(2/15+t2, h(2/15)), (a2-b, h(a0-b)), (a2, 1/4+j), (a2+b, h(a0+t2+b)), (13/80+1/60+t2, h(13/80+1/60+t2))], color='orange', linestyle='--', zorder=-1) + line([(a2-b, h(a0-b)), (a2, 1/4), (a2+b, h(a0+t2+b))], color='orange', linestyle='--', zorder=-1) + line([(a2, 1/4+j), (a2, 1/4)], color='orange', linestyle=':', zorder=-1)
## # bhk_disc_b_zoom.png
## sage: gzigzag =  (l0+l1+l2)+ plot(h, color='grey')+ plot(hd, color='black')
## sage: gzigzag.show(xmin=2/15, xmax=13/80+1/60, ymin=1/4, ymax=3/8, show_legend=False, gridlines=None)
## sage: gzigzag.save('bhk_disc_b_zoom.png', xmin=2/15, xmax=13/80+1/60, ymin=1/4, ymax=3/8, show_legend=False, gridlines=None)

## # 'bhk_disc_b_modified.png'
## sage: gd = plot_covered_intervals(hd)
## sage: gd.save('bhk_disc_b_modified.png', figsize=10, show_legend=False, gridlines=None)

## # bhk_disc_too_many_additivities
## sage: find_epsilon_for_crazy_perturbation(hd, cp, True)

## # change in bhk_discontinuous() #return h_pwl, y_shift, b_shift, h_temp, h_below, h_above, move_shift,  (a1-a0, a2-a0), (bkpts[7], bkpts[8])
## sage: h_pwl, y_shift, b_shift, h_temp, h_below, h_above, move_shift,(t1,t2), (ucl, ucr)=bhk_discontinuous()

## # bhk_y_shift_disc_b.png
## sage: hy = h_pwl+y_shift
## sage: hp = b_shift - h_below + h_above + move_shift
## sage: gyp = plot(hy)+plot(hp, color='magenta')
## sage: gyp.save('bhk_y_shift_disc_b.png')

## # bhk_disc_not_min.png
## sage: h_org = hy+hp
## sage: ghorg =  plot_with_colored_slopes(h_org)
## sage: ghorg.save('bhk_disc_not_min.png', show_legend=False, gridlines=None)

## # bhk_disc_min.png
## #sage: h_org, (t1,t2), (ucl, ucr) = bhk_discontinuous()
## INFO: 2016-04-11 01:09:41,099 Coerced into real number field: Real Number Field in `a` with defining polynomial y^2 - 2
## sage: finite_dimensional_extremality_test(h_org)
## False
## sage: h_pert = h_org._perturbations[0]
## sage: ghorgpert = plot(h_pert, color='magenta')
## sage: (ghorg+ghorgpert).save('bhk_disc_not_min_pert.png', show_legend=False, gridlines=None)
## sage: hmin = lift(h_org)
## sage: ghmin = plot_with_colored_slopes(hmin)
## sage: ghmin.save('bhk_disc_min.png', show_legend=False, gridlines=None)


## sage: minimality_test(hmin)
## INFO: 2016-04-11 01:23:42,962 pi(0) = 0
## INFO: 2016-04-11 01:23:42,963 pi is subadditive.
## INFO: 2016-04-11 01:23:42,965 pi is symmetric.
## INFO: 2016-04-11 01:23:42,965 Thus pi is minimal.
## True
## sage: generate_uncovered_intervals(hmin)
## INFO: 2016-04-11 01:24:02,695 Computing maximal additive faces...
## INFO: 2016-04-11 01:24:06,096 Computing maximal additive faces... done
## INFO: 2016-04-11 01:24:06,485 Completing 33 functional directed moves and 2 covered components...
## INFO: 2016-04-11 01:24:06,851 Completing 29 functional directed moves and 2 covered components...
## INFO: 2016-04-11 01:24:06,939 New dense move from strip lemma: [<Int(0.2737500000000000?, 0.3362500000000000?)>, <Int(0.4637500000000000?, 0.526250000000000?)>]
## INFO: 2016-04-11 01:24:06,960 Completing 0 functional directed moves and 3 covered components...
## INFO: 2016-04-11 01:24:06,961 Completion finished.  Found 0 directed moves and 3 covered components.
## []


## sage: y1 = (hmin(ucl)-hmin.limit(ucl,1)).list()[0]
## sage: y2 = (hmin.limit(ucr,-1)-hmin(ucr)).list()[0]
## sage: c2 = ((hmin.limit(ucr,-1) - hmin.limit(ucl,1))/(ucr-ucl)).list()[0]
## sage: c3 = -5
## sage: b = hmin.end_points()[1].list()[0]
## sage: bb = b*(c2-c3)-y1
## sage: hmin, (t1,t2), (ucl, ucr) = bhk_discontinuous(f=4/5, d1=3/5, d2=5/40, a0=19/100, delta_ratio=sqrt(2)/3, bb=19/23998, c2=5/11999, y1=185/1846, y2=240/11999, field=None)



## sage: finite_dimensional_extremality_test(hmin,show_all_perturbations=True)
## INFO: Finite dimensional test: Solution space has dimension 5
## False
## sage: len(hmin._perturbations)
## 5

## # bhk_disc_min_pert.png
## sage: colours = rainbow(5)
## sage: ghminpert = plot(hmin, color='black')
## sage: for i in range(5):
## ....:     ghminpert += plot(hmin._perturbations[i], color=colours[i])
## ....:     
## sage: ghminpert.save('bhk_disc_min_pert.png', aspect_ratio=1/4)


## sage: perturbs = hmin._perturbations
## sage: len(perturbs)
## 5
## sage: P = perturbation_polyhedron(hmin, perturbs)


## sage: import igp
## sage: from cutgeneratingfunctionology.igp import *
## sage: load("bhk.sage")
## sage: hmin, (t1,t2), (ucl, ucr) = bhk_discontinuous(f=4/5, d1=3/5, d2=5/40, a0=19/100, delta_ratio=sqrt(2)/3, bb=19/23998, c2=5/11999, y1=185/1846, y2=240/11999, field=None)
## INFO: 2016-04-11 02:40:27,881 Coerced into real number field: Real Number Field in `a` with defining polynomial y^2 - 2
## sage: minimality_test(hmin)
## INFO: 2016-04-11 02:40:35,110 pi(0) = 0
## INFO: 2016-04-11 02:40:35,328 pi is subadditive.
## INFO: 2016-04-11 02:40:35,329 pi is symmetric.
## INFO: 2016-04-11 02:40:35,329 Thus pi is minimal.
## True
## sage: finite_dimensional_extremality_test(hmin,show_all_perturbations=True)
## INFO: 2016-04-11 02:40:40,182 Computing maximal additive faces...
## INFO: 2016-04-11 02:40:43,966 Computing maximal additive faces... done
## INFO: 2016-04-11 02:40:44,398 Completing 33 functional directed moves and 2 covered components...
## INFO: 2016-04-11 02:40:44,783 Completing 29 functional directed moves and 2 covered components...
## INFO: 2016-04-11 02:40:44,886 New dense move from strip lemma: [<Int(0.2737500000000000?, 0.3362500000000000?)>, <Int(0.4637500000000000?, 0.526250000000000?)>]
## INFO: 2016-04-11 02:40:44,909 Completing 0 functional directed moves and 3 covered components...
## INFO: 2016-04-11 02:40:44,911 Completion finished.  Found 0 directed moves and 3 covered components.
## INFO: 2016-04-11 02:40:49,148 Finite dimensional test: Solution space has dimension 5
## INFO: 2016-04-11 02:40:49,410 Finding epsilon interval for perturbation...
## INFO: 2016-04-11 02:40:49,556 Finding epsilon interval for perturbation... done.  Interval is [-0.03781372208195416?, 0.05286190515876323?]
## INFO: 2016-04-11 02:40:49,559 Thus the function is NOT extreme.
## INFO: 2016-04-11 02:40:49,821 Finding epsilon interval for perturbation...
## INFO: 2016-04-11 02:40:49,970 Finding epsilon interval for perturbation... done.  Interval is [-0.0702422082017725?, 0.03781372208195416?]
## INFO: 2016-04-11 02:40:50,225 Finding epsilon interval for perturbation...
## INFO: 2016-04-11 02:40:50,381 Finding epsilon interval for perturbation... done.  Interval is [-0.0787640170331994?, 0.006378390149643534?]
## INFO: 2016-04-11 02:40:50,635 Finding epsilon interval for perturbation...
## INFO: 2016-04-11 02:40:50,784 Finding epsilon interval for perturbation... done.  Interval is [-0.03781372208195416?, 0.0702422082017725?]
## INFO: 2016-04-11 02:40:51,042 Finding epsilon interval for perturbation...
## INFO: 2016-04-11 02:40:51,192 Finding epsilon interval for perturbation... done.  Interval is [-0.0424012958362057?, 0.03781372208195416?]
## False
## sage: hlift = lift_until_extreme(hmin)

## sage: extremality_test(hlift)
## True

# has_crazy_perturbation.png
hlift = kzh_minimal_has_only_crazy_perturbation_1()

gcrazy = plot_with_colored_slopes(hlift)
gcrazy.show(aspect_ratio=1/4, figsize=10, show_legend=False, gridlines=None)

h = hlift
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
    

(gcrazy+plot(cp*0.1)).save(destdir + 'has_crazy_perturbation.png', aspect_ratio=1/4, figsize=10, show_legend=False, gridlines=False, ticks=[[],[]], tick_formatter=[[],[]])


## ### logging.debug proof ####
## sage: igp.strategical_covered_components = True
## sage: logging.getLogger().setLevel(logging.DEBUG)
## sage: h = hildebrand_2_sided_discont_2_slope_1()
## sage: extremality_test(h, True)

## sage: h = hildebrand_2_sided_discont_2_slope_1() #disc ext
## sage: h = drlm_backward_3_slope(1/8, 2/8) #cont ext
## sage: h = drlm_backward_3_slope(1/8, 3/8) #cont not ext
## sage: h = example7slopecoarse2() # cont not ext
## sage: h = dr_projected_sequential_merge_3_slope() #cont ext
## sage: h = zhou_two_sided_discontinuous_cannot_assume_any_continuity() #disc not ext
## sage: h = rlm_dpl1_extreme_3a() # disc ext
## sage: h = kzh_minimal_has_only_crazy_perturbation_1() #disc crazy
## sage: extremality_test(h)
## #sage: extremality_test(h, True)


# to show plot in covered interval computation, call ci = generate_covered_components_strategically(h, True)

# sage: igp.strategical_covered_components = False
# sage: logging.getLogger().setLevel(logging.INFO)


####### uniform continuous lemma ######
### right-angle ####
u = 0; v = 0; eta = 1; delta = 0.25;
face = Face(([u,u+eta],[v,v+eta],[u+v,u+v+2*eta]))
t1 = 0.05; t2 = 0.2
x1 = u+t1; x2 = u+t2
y1 = v + delta;
y2 = y1 + (t2 - t1)
y3 = y2 + (t2 - t1)
yn1 = y3 + delta
yn = yn1 + (t2 - t1)
g = face.plot(alpha=0.1, zorder=-10)
g += disk((u,v), 0.1, (0, pi/2), color="mediumspringgreen", zorder=-5)
g += line([(x1, y2), (x2, y1)], color='black', linestyle=':') + line([(x1, y3), (x2, y2)], color='black', linestyle=':') +line([(x1, yn), (x2, yn1)], color='black', linestyle=':')
g += line([(0, y1), (x2, y1)], color='black', linestyle=':')+line([(0, y2), (x2, y2)], color='black', linestyle=':')+line([(0, y3), (x2, y3)], color='black', linestyle=':')+line([(0, yn1), (x2, yn1)], color='black', linestyle=':')+line([(0, yn), (x2, yn)], color='black', linestyle=':')+line([(x1, 0), (x1, yn)], color='black', linestyle=':')+line([(x2, 0), (x2, yn)], color='black', linestyle=':')
g += point([(u,v),(x1,y2),(x2,y1),(x1,y3),(x2,y2),(x1,yn),(x2,yn1)], color='black',size=30,zorder=10)
ticks = [[u, x1, x2, u+delta, u+eta], [v, y1, y2, y3, yn1, yn, v+eta]]
tick_formatter = [["$u$", "$x$", "$x'$", "", r"$u+\eta$"],["$v$","$y_1$","$y_2$","$y_3$","$y_{N-1}$","$y_N$",r"$v+\eta$"]]
g += text("$(u,v)$", (u-0.01,v-0.01), axis_coords=False, vertical_alignment='top', horizontal_alignment='right',color='black',fontsize=20)
g += text(r"$u+\delta$", (u+delta, -0.01), axis_coords=False, vertical_alignment='top', horizontal_alignment='left',color='black',fontsize=20)
g.fontsize(20)
g.show(ticks=ticks, tick_formatter=tick_formatter)
g.save(destdir + "proof_uniform_cont_case_1.png", ticks=ticks, tick_formatter=tick_formatter)

### obtuse-angle ####
u = 1.1; v = 0.1; eta = 1; delta = 0.25;
face = Face(([u-eta,u+eta],[v,v+eta],[u+v,u+v+2*eta]))
t1 = -0.1; t2 = 0.05
x1 = u+t1; x2 = u+t2
y1 = v + delta;
y2 = y1 + (t2 - t1)
y3 = y2 + (t2 - t1)
yn1 = y3 + delta
yn = yn1 + (t2 - t1)
g = face.plot(alpha=0.1, zorder=-10)
g += line([(u,v), (u-eta, v+eta)], color='black')+line([(u,v), (u+eta, v)], color='black')
g += disk((u,v), 0.1, (0, pi*3/4), color="mediumspringgreen", zorder=-5)
g += line([(x1, y2), (x2, y1)], color='black', linestyle=':') + line([(x1, y3), (x2, y2)], color='black', linestyle=':') +line([(x1, yn), (x2, yn1)], color='black', linestyle=':')
g += line([(0, y1), (x2, y1)], color='black', linestyle=':')+line([(0, y2), (x2, y2)], color='black', linestyle=':')+line([(0, y3), (x2, y3)], color='black', linestyle=':')+line([(0, yn1), (x2, yn1)], color='black', linestyle=':')+line([(0, yn), (x2, yn)], color='black', linestyle=':')+line([(x1, 0), (x1, yn)], color='black', linestyle=':')+line([(x2, 0), (x2, yn)], color='black', linestyle=':')
g += point([(u,v),(x1,y2),(x2,y1),(x1,y3),(x2,y2),(x1,yn),(x2,yn1)], color='black',size=30,zorder=10)
ticks = [[u, x1, x2, u-delta, u-eta, u+delta, u+eta], [v, y1, y2, y3, yn1, yn, v+eta]]
tick_formatter = [["$u$", "$x$", "", "", r"$u-\eta$", "", r"$u+\eta$"],["$v$","$y_1$","$y_2$","$y_3$","$y_{N-1}$","$y_N$",r"$v+\eta$"]]
g += text("$(u,v)$", (u-0.02,v), axis_coords=False, vertical_alignment='center', horizontal_alignment='right',color='black',fontsize=20)
g += text(r"$u-\delta$", (u-delta, -0.01), axis_coords=False, vertical_alignment='top', horizontal_alignment='right',color='black',fontsize=20)
g += text(r"$u+\delta$", (u+delta, -0.01), axis_coords=False, vertical_alignment='top', horizontal_alignment='left',color='black',fontsize=20)
g += text("$x'$", (x2, -0.01), axis_coords=False, vertical_alignment='top', horizontal_alignment='left',color='black',fontsize=20)
g.fontsize(20)
g.show(ticks=ticks, tick_formatter=tick_formatter)
g.save(destdir + "proof_uniform_cont_case_2.png", ticks=ticks, tick_formatter=tick_formatter)

### sharp-angle ####
u = 1; v = 0.1; eta = 1; delta = 0.25;
face = Face(([u-eta,u],[v,v+eta],[u+v,u+v+eta]))
t1 = -0.2; t2 = -0.05
x1 = u+t1; x2 = u+t2
y1 = v + delta;
y2 = y1 + (t2 - t1)
y3 = y2 + (t2 - t1)
yn1 = y3 + delta
yn = yn1 + (t2 - t1)
g = face.plot(alpha=0.1, zorder=-10)
g += line([(u,v), (u-eta, v+eta)], color='black')+line([(u,v), (u, v+eta)], color='black')
g += disk((u,v), 0.1, (pi/2, pi*3/4), color="mediumspringgreen", zorder=-5)
g += line([(x1, y2), (x2, y1)], color='black', linestyle=':') + line([(x1, y3), (x2, y2)], color='black', linestyle=':') +line([(x1, yn), (x2, yn1)], color='black', linestyle=':')
g += line([(0, y1), (x2, y1)], color='black', linestyle=':')+line([(0, y2), (x2, y2)], color='black', linestyle=':')+line([(0, y3), (x2, y3)], color='black', linestyle=':')+line([(0, yn1), (x2, yn1)], color='black', linestyle=':')+line([(0, yn), (x2, yn)], color='black', linestyle=':')+line([(x1, 0), (x1, yn)], color='black', linestyle=':')+line([(x2, 0), (x2, yn)], color='black', linestyle=':')
g += point([(u,v),(x1,y2),(x2,y1),(x1,y3),(x2,y2),(x1,yn),(x2,yn1)], color='black',size=30,zorder=10)
ticks = [[u, x1, x2, u-delta, u-eta], [v, y1, y2, y3, yn1, yn, v+eta]]
tick_formatter = [["$u$", "$x$", "$x'$", "", r"$u-\eta$"],["$v$","$y_1$","$y_2$","$y_3$","$y_{N-1}$","$y_N$",r"$v+\eta$"]]
#g += text("$v$", (-0.01,v), axis_coords=False, vertical_alignment='center', horizontal_alignment='right',color='black',fontsize=20)
g += text("$(u,v)$", (u,v), axis_coords=False, vertical_alignment='top', horizontal_alignment='left',color='black',fontsize=20)
g += text(r"$u-\delta$", (u-delta, -0.01), axis_coords=False, vertical_alignment='top', horizontal_alignment='right',color='black',fontsize=20)
g.fontsize(20)
g.show(ticks=ticks, tick_formatter=tick_formatter)
g.save(destdir + "proof_uniform_cont_case_3.png", ticks=ticks, tick_formatter=tick_formatter)

## ### sharp-angle consider p2 p3 ####
## u = 0; v = 0; eta = 1
## face = Face(([u,u+eta],[v,v+eta],[u+v,u+v+eta]))
## g = face.plot()
## g += point([(u+eta,v)], color='black',size=30,zorder=10)
## g += text("$(u,v)$", (u+eta,v-0.01), axis_coords=False, vertical_alignment='top', horizontal_alignment='center',color='black',fontsize=20)
## g.show(axes=False)
## g.save("proof_uniform_cont_case_4.png", axes=False)

### V for case 1 and 2 ####
u = 1.1; v = 0.1; delta = eta = 1;
face = Face(([u-delta,u+delta],[v,v+delta],[u+v,u+v+2*delta]))
t1 = 0.3; t2 = 0.75;
x2 = u+t1; x1 = u+t2;
y1 = v+t1; y2 = v+t2;
g = face.plot(alpha=0.1, zorder=-10)
g += line([(u,v), (u-eta, v+eta)], color='black')+line([(u,v), (u+eta, v)], color='black') + line([(u,v), (u, v+eta)], color='black')
g += line([(x1, y1), (x2, y2)], color='black', linestyle=':')
g += line([(0, y1), (x1, y1)], color='black', linestyle=':')+line([(0, y2), (x2, y2)], color='black', linestyle=':')+line([(x1, 0), (x1, y1)], color='black', linestyle=':')+line([(x2, 0), (x2, y2)], color='black', linestyle=':')
g += point([(u,v),(x1,y1),(x2,y2)], color='black',size=30,zorder=10)
g += text("$(u,v)$", (u-0.02,v), axis_coords=False, vertical_alignment='center', horizontal_alignment='right',color='black',fontsize=20)
ticks = [[u, x1, x2, u-delta, u+delta], [v, y1, y2, v+delta]]
tick_formatter = [["$u$", "$x$", "$x'$", r"$u-\delta$", r"$u+\delta$"],["$v$","$y$","$y'$",r"$v+\delta$"]]
g.fontsize(20)
g.show(ticks=ticks, tick_formatter=tick_formatter)
g.save(destdir + "proof_uniform_V_case12.png", ticks=ticks, tick_formatter=tick_formatter)

### V for case 3 ####
u = 1.1; v = 0.1; eta = delta = 1;
face = Face(([u-eta,u],[v,v+eta],[u+v-eta,u+v+eta]))
t1 = 0.3; t2 = 0.75
y1 = v + t1; y2 = v+t2
x1 = u - t1/2; x2 = x1+y1-y2
g = face.plot(alpha=0.1, zorder=-10)
g += line([(u,v), (u-eta, v+eta)], color='black')+line([(u,v), (u, v+eta)], color='black')+line([(u,v), (u-eta, v)], color='black')
g += line([(x1, y1), (x2, y2)], color='black', linestyle=':')
g += line([(0, y1), (u, y1)], color='black', linestyle=':')+line([(0, y2), (x2, y2)], color='black', linestyle=':')+line([(x1, 0), (x1, y1)], color='black', linestyle=':')+line([(x2, 0), (x2, y2)], color='black', linestyle=':')
g += point([(u,v),(x1,y1),(x2,y2)], color='black',size=30,zorder=10)
ticks = [[u, x1, x2, u-delta], [v, y1, y2, v+delta]]
tick_formatter = [["$u$", "$x$", "$x'$",r"$u-\delta$"],["$v$","$y$","$y'$",r"$v+\delta$"]]
g += text("$(u,v)$", (u,v), axis_coords=False, vertical_alignment='top', horizontal_alignment='left',color='black',fontsize=20)
g += text("=",(u-t1/4,v+t1+0.01), rotation="vertical", axis_coords=False, vertical_alignment='center', horizontal_alignment='center',color='blue', fontsize=15)+text("=",(u-t1*3/4,v+t1+0.01), rotation="vertical",axis_coords=False, vertical_alignment='center', horizontal_alignment='center',color='blue', fontsize=15)
g.fontsize(20)
g.show(ticks=ticks, tick_formatter=tick_formatter,figsize=5.9)
g.save(destdir + "proof_uniform_V_case3.png", ticks=ticks, tick_formatter=tick_formatter,figsize=5.9)


### W for case 1 and 2 ####
u = 0; v = 0.5; delta = eta = 1;
face = Face(([u-delta,u+delta],[v,v+delta],[u+v,u+v+2*delta]))
t1 = 0.3; t2 = 0.75;
z1 = u+v+t1; z2 = u+v+t2;
y = v+t1/2
x1 = z1-y; x2 = z2-y;
g = face.plot(alpha=0.1, zorder=-10)
g += line([(u,v), (u-eta, v+eta)], color='black')+line([(u,v), (u+eta, v)], color='black') + line([(u,v), (u, v+eta)], color='black')
g += line([(u, y), (x2, y)], color='black', linestyle=':')
g += line([(u, v), (u+v, 0)], color='black', linestyle=':')+line([(u, v+t1), (z1, 0)], color='black', linestyle=':')+line([(x2, y), (z2, 0)], color='black', linestyle=':')+line([(u+delta, v), (u+v+delta, 0)], color='black', linestyle=':')
g += point([(u,v),(x1,y),(x2,y)], color='black',size=30,zorder=10) #(u+delta,v),(u, v+delta)
g += text("$(u,v)$", (u-0.02,v), axis_coords=False, vertical_alignment='center', horizontal_alignment='right',color='black',fontsize=20)+text(r"$(u+\delta,v)$", (u+delta+0.02,v), axis_coords=False, vertical_alignment='center', horizontal_alignment='left',color='black',fontsize=20)+text(r"$(u,v+\delta)$", (u,v+delta+0.02), axis_coords=False, vertical_alignment='bottom', horizontal_alignment='center',color='black',fontsize=20)+text(r"$(x,y)$", (x1+0.02,y+0.02), axis_coords=False, vertical_alignment='bottom', horizontal_alignment='left',color='black',fontsize=20)+text("$(x',y)$", (x2+0.02,y+0.02), axis_coords=False, vertical_alignment='bottom', horizontal_alignment='left',color='black',fontsize=20)
g += arrow((0.3,0), (u+v+delta+0.4,0), color='black',width=1)
zticks = [u+v, z1, z2, u+v+delta+0.05]
zlocs = [u+v, z1, z2, u+v+delta]
ztick_formatter = ["$u+v$", "$z$", "$z'$", r"$u+v+\delta$"]
for i in range(4):
    g += text(ztick_formatter[i], (zticks[i], -0.01), axis_coords=False, vertical_alignment='top', horizontal_alignment='center',color='black',fontsize=20)
    g += line([(zlocs[i],0),(zlocs[i],0.03)],color='black',thickness=0.5)
g += text("=",(u+t1/4,v+t1/2+0.01), rotation="vertical", axis_coords=False, vertical_alignment='center', horizontal_alignment='center',color='blue', fontsize=15)+text("=",(u,v+t1/4), axis_coords=False, vertical_alignment='center', horizontal_alignment='center',color='blue', fontsize=15)+text("=",(u,v+t1*3/4), axis_coords=False, vertical_alignment='center', horizontal_alignment='center',color='blue', fontsize=15)
# unsaved change:
#g += text("=",(u,v+t1/4), axis_coords=False, vertical_alignment='center', horizontal_alignment='center',color='blue', fontsize=15)+text("=",(u,v+t1*3/4), axis_coords=False, vertical_alignment='center', horizontal_alignment='center',color='blue', fontsize=15) #+ text("=",(u+t1/4,v+t1/2+0.01), rotation="vertical", axis_coords=False, vertical_alignment='center', horizontal_alignment='center',color='blue', fontsize=15)+
g.fontsize(20)
g.show(axes=False)
g.save(destdir + "proof_uniform_W_case12.png",axes=False)

### W for case 3a ####
u = eta; v = 0.5; eta =1; delta = 1;
face = Face(([u-eta,u],[v,v+eta],[u+v,u+v+eta]))
t1 = 0.3; t2 = 0.75
z1 = u+v+t1; z2 = u+v+t2;
y = v+delta
x1 = z1-y; x2 = z2-y;
g = face.plot(alpha=0.1, zorder=-10)
g += line([(u,v), (u-eta, v+eta)], color='black')+ line([(u,v), (u, v+eta)], color='black')
g += line([(u, y), (0, y)], color='black', linestyle=':')
g += line([(u, v), (u+v, 0)], color='black', linestyle=':')+line([(x1, y), (z1, 0)], color='black', linestyle=':')+line([(x2, y), (z2, 0)], color='black', linestyle=':')+line([(u, y), (u+v+delta, 0)], color='black', linestyle=':')
g += point([(u,v),(x1,y),(x2,y)], color='black',size=30,zorder=10) #(u+delta,v),(u, v+delta)
g += text("$(u,v)$", (u-0.02,v), axis_coords=False, vertical_alignment='center', horizontal_alignment='right',color='black',fontsize=20)+text("$(x,y)$", (x1,y+0.02), axis_coords=False, vertical_alignment='bottom', horizontal_alignment='center',color='black',fontsize=20)+text("$(x',y)$", (x2,y+0.02), axis_coords=False, vertical_alignment='bottom', horizontal_alignment='center',color='black',fontsize=20)
g += text(r"$y=v+\delta$", (u+0.02,v+delta), axis_coords=False, vertical_alignment='center', horizontal_alignment='left',color='black',fontsize=20)
g += arrow((0.3,0), (u+v+delta+0.4,0), color='black',width=1)
zticks = [u+v, z1, z2, u+v+delta+0.05]
zlocs = [u+v, z1, z2, u+v+delta]
ztick_formatter = ["$u+v$", "$z$", "$z'$", r"$u+v+\delta$"]
for i in range(4):
    g += text(ztick_formatter[i], (zticks[i], -0.01), axis_coords=False, vertical_alignment='top', horizontal_alignment='center',color='black',fontsize=20)
    g += line([(zlocs[i],0),(zlocs[i],0.03)],color='black',thickness=0.5)
g.fontsize(20)
g.show(axes=False)
g.save(destdir + "proof_uniform_W_case3a.png",axes=False)


### W for case 3b ####
u = eta; v = 0.5; eta =1; delta = 1;
face = Face(([u-eta,u],[v,v+eta],[u+v-eta,u+v+eta]))
t1 = -0.2; t2 = 0.25
z1 = u+v+t1; z2 = u+v+t2;
y = v+delta/2
x1 = z1-y; x2 = z2-y;
g = face.plot(alpha=0.1, zorder=-10)
g += line([(u,v), (u-eta, v)], color='black')+ line([(u,v), (u, v+eta)], color='black')
g += line([(u-delta, y), (u, y)], color='black', linestyle=':')
g += line([(u, v), (u+v, 0)], color='black', linestyle=':')+line([(x1, y), (z1, 0)], color='black', linestyle=':')+line([(x2, y), (z2, 0)], color='black', linestyle=':')+line([(u, y), (u+y, 0)], color='black', linestyle=':')+line([(u-delta, y), (u-delta+y, 0)], color='black', linestyle=':')
g += point([(u,v),(x1,y),(x2,y)], color='black',size=30,zorder=10) #(u+delta,v),(u, v+delta)
g += text("$(u,v)$", (u+0.02,v), axis_coords=False, vertical_alignment='center', horizontal_alignment='left',color='black',fontsize=20)+text("$(x,y)$", (x1,y+0.02), axis_coords=False, vertical_alignment='bottom', horizontal_alignment='center',color='black',fontsize=20)+text("$(x',y)$", (x2,y+0.02), axis_coords=False, vertical_alignment='bottom', horizontal_alignment='center',color='black',fontsize=20)
g += text(r"$y=v+\frac{\delta}{2}$", (u+0.02,y), axis_coords=False, vertical_alignment='center', horizontal_alignment='left',color='black',fontsize=20)
g += arrow((0.3,0), (u+y+0.4,0), color='black',width=1)
zticks = [u+y-delta-0.05,u+v, z1, z2, u+y+0.05]
zlocs = [u+y-delta, u+v, z1, z2, u+y]
ztick_formatter = [r"$u+v-\frac{\delta}{2}$", "$u+v$", "$z$", "$z'$", r"$u+v+\frac{\delta}{2}$"]
for i in range(5):
    g += text(ztick_formatter[i], (zticks[i], -0.01), axis_coords=False, vertical_alignment='top', horizontal_alignment='center',color='black',fontsize=20)
    g += line([(zlocs[i],0),(zlocs[i],0.03)],color='black',thickness=0.5)
g.fontsize(20)
g.show(axes=False)
g.save(destdir + "proof_uniform_W_case3b.png",axes=False)


#### MAIN ####
os.system("cd %s && (pdflatex -synctex=1 -src-specials -interaction=nonstopmode crazy_perturbation_graphics)" % (destdir,)) 
