"""
Automatic verification of the paper "All Cyclic Group Facets Inject".

We check cases (a') and (b') of the subadditivity proof
in the paper :cite:`koeppe-zhou:cyclic-group-facets-inject`
through symbolic computation.

::

    sage: import cutgeneratingfunctionology.igp.procedures.injective_2_slope_fill_in_proof as proof
    sage: from cutgeneratingfunctionology.igp.procedures.injective_2_slope_fill_in_proof import *
    sage: logging.disable(logging.INFO) # Suppress output in automatic tests.

    sage: setup_case_aprime1_MMM_type_I()
    sage: setup_case_aprime2_MMM_type_II()
    sage: setup_case_aprime3_MWW()

    sage: setup_case_bprime2_MMM()
    sage: setup_case_bprime3_MWW_type_I()
    sage: setup_case_bprime3_MWW_type_II()
    sage: setup_case_bprime4_WMW()

"""

from cutgeneratingfunctionology.igp import *

### Library.

MW = [None, 'M', 'W']

def construct_phi(pi, signs):
    global IJK_plusminus
    IJK_plusminus = [None, [None, None, None], [None, None, None]]
    phi = [None, None, None]
    for i in range(3):
        end_points = pi[i].end_points()
        x = end_points[0]
        fx = pi[i](end_points[0])
        b = [x]
        v = [fx]
        for j in range(1, len(end_points)):
            for k in [signs[i], -signs[i]]:
                old_x = x
                x  += d_plusminus[k][i]
                fx += d_plusminus[k][i] * s_plusminus[k]
                b.append(x)
                v.append(fx)
                IJK_plusminus[k][i] = [min(old_x, x), max(old_x, x)]
            assert x == end_points[j]
            assert fx == pi[i](x)
        phi[i] = piecewise_function_from_breakpoints_and_values(b, v,
                                                                merge=False, field=ParametricRealField)
    return phi

def plot_case(phi, pts):
    g = F.plot()
    g += line([[x, 0], [x, 1]], color='grey', legend_label='Complex Delta P_phi')
    g += line([[0, y], [1, y]], color='grey')
    g += line([(0, z), (z, 0)], color='grey')
    g += points([ p[0:2] for p in pts], zorder=10)
    return g

def setup_case(sign1, sign2, sign3, show_plots=False):
    global phi, x, y, z, P12, P13, P23

    signs = (sign1, sign2, sign3)

    ## with K.off_the_record():
    ##     plot_fun_IJK(pi, phi).show(figsize=[8,2.5])

    phi = construct_phi(pi, signs)

    (x, y, z) = [ phi_i.end_points()[1] for phi_i in phi ]

    P12 = (x, y, x + y)
    P13 = (x, z - x, z)
    P23 = (z - y, y, z)

    if show_plots:
        plot_case(phi, [P12, P13, P23]).show(legend_title='Signs {}'.format([MW[s] for s in signs]))


def plot_fun_IJK(pi, phi):
    plots = [ pi[i].plot(color='black', **ticks_keywords(pi[i]))
              + phi[i].plot(color='red', **ticks_keywords(phi[i]))
              for i in range(3) ]

    return graphics_array(plots)

def delta_IJK(pi, xy):
    x = xy[0]
    y = xy[1]
    if len(xy) == 3:
        assert x + y == xy[2]
    return pi[0](x) + pi[1](y) - pi[2](x+y)

### Cases.

########################################
## (a') Right angle, bottom left.
########################################

def setup_pi_case_aprime():
    global KK, inv_mq, u, v, w, pi_u, pi_v, pi_w, s_1, s_2, s_3, s_p, s_m
    KK = ParametricRealField([QQ('1/6'), QQ('0'), QQ('0'), QQ('1/4'), QQ('1/8'),
                              QQ('1/3'), QQ('1/8'), QQ('-1/8'), QQ('1/2'), QQ('-1/4')],
                             names=['inv_mq', 'u', 'v', 'pi_u', 'pi_v',
                                    's_1', 's_2', 's_3', 's_p', 's_m'],
                             mutable_values=True, big_cells=True, allow_refinement=False)
    inv_mq, u, v, pi_u, pi_v, s_1, s_2, s_3, s_p, s_m = KK.gens()

    assert 0 < inv_mq < 1
    assert s_p > s_m

    global s
    s = [s_1, s_2, s_3]
    w = u + v
    pi_w = pi_u + pi_v

    global I, J, K, IJK, F, pi
    I = [u, u + inv_mq]
    J = [v, v + inv_mq];
    K = [w + inv_mq, w + 2 * inv_mq]; #K2 = [w] + K
    IJK = (I, J, K)
    F = Face([I, J, K])

    global pi
    pi = [piecewise_function_from_breakpoints_and_values(I,
                                                         [ pi_u + (x - u) * s[0] for x in I ],
                                                         field=ParametricRealField),
          piecewise_function_from_breakpoints_and_values(J,
                                                         [ pi_v + (y - v) * s[1] for y in J ],
                                                         field=ParametricRealField),
          piecewise_function_from_breakpoints_and_values(K,
                                                         [ pi_w + (z - w) * s[2] for z in K ],
                                                         field=ParametricRealField, merge=False)]

    assert delta_IJK(pi, (u+inv_mq, v)) >= 0   # s_1 >= s_3
    assert delta_IJK(pi, (u, v+inv_mq)) >= 0   # s_2 >= s_3

    global d1_plus, d2_plus, d3_plus, d_plus, d1_minus, d2_minus, d3_minus, d_minus, d_plusminus, s_plusminus
    d1_plus, d2_plus, d3_plus = d_plus  = [ inv_mq * (s_i - s_m) / (s_p - s_m) for s_i in s ]
    d1_minus, d2_minus, d3_minus = d_minus = [ inv_mq * (s_p - s_i) / (s_p - s_m) for s_i in s ]
    assert all(d_plus[i] > 0 and d_minus[i] > 0 for i in range(3))

    KK.freeze()

    assert delta_IJK(pi, (u+inv_mq, v)) == inv_mq * (s_1 - s_3) >= 0
    assert delta_IJK(pi, (u, v+inv_mq)) == inv_mq * (s_2 - s_3) >= 0

    assert all(d_plus[i] + d_minus[i] == inv_mq for i in range(3))

    d_plusminus = [ None, d_plus, d_minus ]   # index by 1, -1.
    s_plusminus = [ None, s_p, s_m ]

####

def setup_case_aprime_MMM(show_plots=False):
    setup_pi_case_aprime()
    setup_case(+1, +1, +1, show_plots=show_plots)
    assert P13 in F
    assert P23 in F

def setup_case_aprime1_MMM_type_I(show_plots=False):
    logging.info("a'1 MMM Type I")
    setup_case_aprime_MMM()
    KK.change_values(s_3=1/10)
    if show_plots:
        plot_case(phi, [P12, P13, P23]).show(legend_title="Case a' MMM Type I")
    with KK.unfrozen():
        assert P12[2] <= z
    with KK.temporary_assumptions():
        with KK.unfrozen():
            assert K[0] < P12[2]
        assert P12 in F
        assert delta_IJK(phi, P12) == d_minus[2] * (s_p - s_m) == inv_mq * (s_p - s_3) >= 0
    assert delta_IJK(phi, P13) == delta_IJK(phi, P23) == inv_mq * ((s_1 - s_3) + (s_2 - s_3)) >= 0

def setup_case_aprime2_MMM_type_II(show_plots=False):
    logging.info("a'2 MMM Type II")
    setup_case_aprime_MMM(show_plots=show_plots)
    with KK.unfrozen():
        assert P12[2] > z
    assert P12 in F
    assert delta_IJK(phi, P12) == inv_mq * ((s_1 - s_3) + (s_2 - s_3)) >= 0
    assert delta_IJK(phi, P13) == delta_IJK(phi, P23) == inv_mq * (s_p - s_3) >= 0

def setup_case_aprime3_MWW(show_plots=False):
    logging.info("a'3 MWW")
    setup_pi_case_aprime()
    setup_case(+1, -1, -1)
    KK.change_values(s_3=1/10)
    if show_plots:
        plot_case(phi, [P12, P13, P23]).show(legend_title="Case a' MWW")

    assert P12[2] <= z

    with KK.temporary_assumptions(case_id="a' MWW P12"):
        with KK.unfrozen():
            assert P12[2] > K[0]

        assert delta_IJK(phi, P12) == inv_mq * (s_1 - s_3) >= 0

    with KK.temporary_assumptions(case_id="a' MWW P13"):
        with KK.unfrozen():
            assert P13[1] < J[1]
        assert x - u > z - u - v - inv_mq
        assert is_pt_in_interval(IJK_plusminus[+1][1], P13[1])
        assert delta_IJK(phi, P13) == inv_mq * ((s_p - s_3) + (s_2 - s_3)) >= 0

    assert P23[0] >= I[1]

########################################
## (b') Acute angle, top left.
########################################

def setup_pi_case_bprime():
    global KK, inv_mq, u, v, w, pi_u, pi_v, pi_w, s_1, s_2, s_3, s_p, s_m
    KK = ParametricRealField([QQ('1/6'), QQ('0'), QQ('1'), QQ('1/4'), QQ('1/8'),
                              QQ('1/4'), QQ('1/16'), QQ('1/8'), QQ('1/2'), QQ('-1/4')],
                             names=['inv_mq', 'u', 'v', 'pi_u', 'pi_v',
                                    's_1', 's_2', 's_3', 's_p', 's_m'],
                             mutable_values=True, big_cells=True, allow_refinement=False)
    inv_mq, u, v, pi_u, pi_v, s_1, s_2, s_3, s_p, s_m = KK.gens()

    assert inv_mq > 0

    assert s_p > s_m

    global s
    s = [s_1, s_2, s_3]

    w = u + v
    pi_w = pi_u + pi_v

    global I, J, K, IJK, F, pi
    I = [u, u + inv_mq]
    J = [v - 2*inv_mq, v - inv_mq]; #J2 = [v - 2*inv_mq, v - inv_mq, v]
    K = [w - inv_mq, w]
    IJK = (I, J, K)
    F = Face([I, J, K])

    global pi
    pi = [piecewise_function_from_breakpoints_and_values(I,
                                                         [ pi_u + (x - u) * s[0] for x in I ],
                                                         field=ParametricRealField),
          piecewise_function_from_breakpoints_and_values(J,
                                                         [ pi_v + (y - v) * s[1] for y in J ],
                                                         field=ParametricRealField, merge=False),
          piecewise_function_from_breakpoints_and_values(K,
                                                         [ pi_w + (z - w) * s[2] for z in K ],
                                                         field=ParametricRealField)]

    assert delta_IJK(pi, (u, v-inv_mq)) >= 0          # s_2 <= s_1
    assert delta_IJK(pi, (u+inv_mq, v-inv_mq)) >= 0   # s_2 <= s_3

    global d1_plus, d2_plus, d3_plus, d_plus, d1_minus, d2_minus, d3_minus, d_minus, d_plusminus, s_plusminus
    d1_plus, d2_plus, d3_plus = d_plus  = [ inv_mq * (s_i - s_m) / (s_p - s_m) for s_i in s ]
    d1_minus, d2_minus, d3_minus = d_minus = [ inv_mq * (s_p - s_i) / (s_p - s_m) for s_i in s ]
    assert all(d_plus[i] > 0 and d_minus[i] > 0 for i in range(3))

    KK.freeze()

    assert delta_IJK(pi, (u, v-inv_mq)) == inv_mq * (s_3 - s_2) >= 0
    assert delta_IJK(pi, (u+inv_mq, v-inv_mq)) == inv_mq * (s_1 - s_2) >= 0

    assert all(d_plus[i] + d_minus[i] == inv_mq for i in range(3))

    d_plusminus = [ None, d_plus, d_minus ]   # index by 1, -1.
    s_plusminus = [ None, s_p, s_m ]

    # Trivial consequences of the slope inequalities
    assert d2_plus <= d1_plus and d2_plus <= d3_plus
    assert d2_minus >= d1_minus and d2_minus >= d3_minus

####

def setup_case_bprime2_MMM(show_plots=False):
    logging.info("b'2 MMM")
    setup_pi_case_bprime()
    setup_case(+1, +1, +1, show_plots=show_plots)

    assert P23[0] >= I[1]             # not in interior of F

    assert P13[1] >= y
    with KK.temporary_assumptions(case_id="b'2 MMM P13"):
        with KK.unfrozen():
            assert P13[1] < J[1]
        #phi[1](P13[1])
        assert phi[0](P13[0]) - phi[0](u) == d1_plus * s_p
        assert phi[2](P13[2]) - phi[2](w) == d3_minus * (- s_m)

        assert delta_IJK(phi, P13) >= 0

    #assert P12[0] + P12[1] < K[0]  # unknown
    with KK.temporary_assumptions(case_id="b'2 MMM P12"):
        with KK.unfrozen():
            assert P12[2] > K[0]
        assert phi[2].which_function(P12[2])._slope == s_p
        assert delta_IJK(phi, P12) >= 0

####

def setup_case_bprime3_MWW(show_plots=False):
    setup_pi_case_bprime()
    setup_case(+1, -1, -1, show_plots=show_plots)
    assert P12 in F
    assert P23 in F

def setup_case_bprime3_MWW_type_I(show_plots=False):
    logging.info("b'3 MWW Type I")
    setup_case_bprime3_MWW()
    with KK.unfrozen():
        assert P23[0] > x
    assert P13[1] > y
    assert P12[2] < z

    assert delta_IJK(phi, P12) >= 0

    assert delta_IJK(phi, P23) >= 0

    with KK.temporary_assumptions():
        with KK.unfrozen():
            assert P13[1] < J[1]
        assert F.interior_contains(P13)
        assert phi[1].which_function(P13[1])._slope == s_p
        assert delta_IJK(phi, P13) == inv_mq * (s_p - s_2) >= 0

def setup_case_bprime3_MWW_type_II(show_plots=False):
    logging.info("Case b'3 MWW Type II")
    setup_case_bprime3_MWW()
    KK.change_values(s_2=QQ('-1/5'))
    if show_plots:
        plot_case(phi, [P12, P13, P23]).show(legend_title="Case b'3 MWW Type II")
    with KK.unfrozen():
        assert P23[0] <= x
    assert z - y == P23[0] <= x
    assert P13[1] <= y
    assert P12[2] >= z

    #assert phi[2].which_function(P12[2])._slope == s_p      # does not work at breakpoints!
    assert is_pt_in_interval(IJK_plusminus[+1][2], P12[2])
    assert delta_IJK(phi, P12) >= 0

    assert P13 in F
    #assert phi[1].which_function(P13[1])._slope == s_m
    assert is_pt_in_interval(IJK_plusminus[-1][1], P13[1])
    assert delta_IJK(phi, P13) >= 0

    #assert phi[0].which_function(P23[0])._slope == s_p
    assert is_pt_in_interval(IJK_plusminus[+1][0], P23[0])
    assert delta_IJK(phi, P23) >= 0

####

def setup_case_bprime4_WMW(show_plots=False):
    logging.info("Case b'4 WMW")
    setup_pi_case_bprime()
    setup_case(-1, +1, -1, show_plots=show_plots)

    assert not F.interior_contains(P12)

    with KK.changed_values(s_3=QQ('1/3')):
        with KK.temporary_assumptions(case_id="b'4 WMW P23"):
            with KK.unfrozen():
                assert F.interior_contains(P23)
            assert phi[0].which_function(P23[0])._slope == s_p
            assert delta_IJK(phi, P23) >= 0
        with KK.temporary_assumptions(case_id="b'4 WMW P13"):
            with KK.unfrozen():
                assert F.interior_contains(P13)
            assert phi[1].which_function(P13[1])._slope == s_m
            assert delta_IJK(phi, P13) >= 0
