########## facet paper ###########

sage: h = hildebrand_discont_3_slope_1()
sage: hh = discontinuous_interpolation([0,1/2,5/8,7/8],[0,1,3/4,1/4],[0,1/2,3/4,1/4],[1/2,1,3/4,1/4])
#sage: h3 = restrict_to_finite_group(h, oversampling=3)
#sage: extremality_test(h3, True)
sage: g = plot_2d_diagram_additive_domain_sans_limits(h)
sage: g.show(show_legend=False)
sage: g.save('2d_simple_E_pi_extreme_not_facet.png', show_legend=False)
Launched png viewer for Graphics object consisting of 98 graphics primitives
sage: gg = plot_2d_diagram_additive_domain_sans_limits(hh)
sage: gg.show(show_legend=False)
sage: gg.save('2d_simple_E_pi_prime_extreme_not_facet.png', show_legend=False)

sage: fn = kzh_minimal_has_only_crazy_perturbation_1()
sage: gen = set(itertools.chain(generate_type_1_vertices(fn, operator.gt), \
                                generate_type_2_vertices(fn, operator.gt)))
sage: dp = [delta_pi_general(fn, x, y, (xeps, yeps, zeps)) for (x, y, z, xeps, yeps, zeps) in gen]
sage: igp.show_RNFElement_by_embedding=False
sage: assert min(dp) == 19/23998


sage: h = kzh_minimal_has_only_crazy_perturbation_1()
sage: x = h.end_points()
sage: l=x[17]; u = x[18]; ll = x[19]; uu=x[20]; f=x[37]
sage: g = plot_2d_complex(h)
sage: g += polygon([(0, l), (0, u), (1, u), (1, l)], color='yellow',fill=True, zorder=-5)
sage: g += polygon([(0, ll), (0, uu), (1, uu), (1, ll)], color='yellow',fill=True, zorder=-5)
sage: g += polygon([(ll,0), (uu,0), (uu,1), (ll,1)], color='yellow',fill=True, zorder=-5)
sage: g += polygon([(l,0), (u,0), (u,1), (l,1)], color='yellow',fill=True, zorder=-5)
sage: g += polygon([(l,0), (u,0), (0,u), (0,l)], color='yellow',fill=True, zorder=-5)
sage: g += polygon([(ll,0), (uu,0), (0,uu), (0,ll)], color='yellow',fill=True, zorder=-5)
sage: g += polygon([(l,1), (u,1), (1,u), (1,l)], color='yellow',fill=True, zorder=-5)
sage: g += polygon([(ll,1), (uu,1), (1,uu), (1,ll)], color='yellow',fill=True, zorder=-5)
sage: g += polygon([(l,1), (u,1), (u,1-w)], color='red',fill=True, zorder=-5)
sage: g += polygon([(ll,1), (uu,1), (uu,1-w)], color='red',fill=True, zorder=-5)
sage: g += polygon([(0,l), (0,u), (w,l)], color='red',fill=True, zorder=-5)
sage: g += polygon([(0,ll), (0,uu), (w,ll)], color='red',fill=True, zorder=-5)
sage: g += polygon([(1,l), (1,u), (1-w,u)], color='red',fill=True, zorder=-5)
sage: g += polygon([(1,ll), (1,uu), (1-w,uu)], color='red',fill=True, zorder=-5)
sage: g += polygon([(l,0), (u,0), (l,w)], color='red',fill=True, zorder=-5)
sage: g += polygon([(ll,0), (uu,0), (ll,w)], color='red',fill=True, zorder=-5)
sage: g += polygon([(l,l), (l,u), (u,u), (u,l)], color='red',fill=True, zorder=-5)
sage: g += polygon([(ll,l), (ll,u), (uu,u), (uu,l)], color='red',fill=True, zorder=-5)
sage: g += polygon([(l,ll), (l,uu), (u,uu), (u,ll)], color='red',fill=True, zorder=-5)
sage: g += polygon([(ll,ll), (ll,uu), (uu,uu), (uu,ll)], color='red',fill=True, zorder=-5)
sage: g += polygon([(l,ll-l), (l,uu-l), (u,uu-u), (u,ll-u)], color='red',fill=True, zorder=-5)
sage: g += polygon([(ll-l,l), (uu-l,l), (uu-u,u), (ll-u,u)], color='red',fill=True, zorder=-5)
sage: g += polygon([(ll,1+l-ll), (ll,1+u-ll), (uu,1+u-uu), (uu,1+l-uu)], color='red',fill=True, zorder=-5)
sage: g += polygon([(1+l-ll,ll), (1+u-ll,ll), (1+u-uu,uu), (1+l-uu,uu)], color='red',fill=True, zorder=-5)
sage: tk=[0,l,u,ll,uu,f,1];
sage: tkf=[0,"$l$","$u$","$f-l$","$f-u$","$f$",1]
sage: g.show(ticks=[tk,tk], tick_formatter=[tkf, tkf], show_legend=False)
sage: g.save('2d_crazy_nf.png',ticks=[tk,tk], tick_formatter=[tkf, tkf], show_legend=False)
