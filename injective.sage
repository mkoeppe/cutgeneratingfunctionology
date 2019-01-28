from cutgeneratingfunctionology.igp import *

### Library.

def construct_phi(pi, signs):
    phi = [None, None, None]
    for i in range(3):
        end_points = pi[i].end_points()
        x = end_points[0]
        fx = pi[i](end_points[0])
        b = [x]
        v = [fx]
        for j in range(1, len(end_points)):
            for k in [signs[i], -signs[i]]:
                x  += d_plusminus[k][i]
                fx += d_plusminus[k][i] * s_plusminus[k]
                b.append(x)
                v.append(fx)
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


def plot_fun_IJK(pi, phi):
    plots = [ pi[i].plot(color='black', **ticks_keywords(pi[i]))
              + phi[i].plot(color='red', **ticks_keywords(phi[i]))
              for i in range(3) ]

    return graphics_array(plots)

def delta_IJK(pi, xy):
    x = xy[0]
    y = xy[1]
    return pi[0](x) + pi[1](y) - pi[2](x+y)

### Cases.

KK.<inv_mq, u_prime, v_prime, pi_u_prime, pi_v_prime, s_1, s_2, s_3, s_p, s_m> = ParametricRealField([1/6, 0, 1, 1/4, 1/8, 1/4, 1/16, 1/8, 1/2, -1/4], mutable_values=True)

assert inv_mq > 0

s = [s_1, s_2, s_3]
w_prime = u_prime + v_prime
pi_w_prime = pi_u_prime + pi_v_prime

## Acute angle, top left.

assert s_2 <= s_1
assert s_2 <= s_3

I = [u_prime, u_prime + inv_mq]
J = [v_prime - 2*inv_mq, v_prime - inv_mq]; JJ = [v_prime - 2*inv_mq, v_prime - inv_mq, v_prime]
K = [w_prime - inv_mq, w_prime]
IJK = (I, J, K)
F = Face([I, J, K])

pi = [piecewise_function_from_breakpoints_and_values(I,
                                                     [ pi_u_prime + (x - u_prime) * s[0] for x in I ],
                                                     field=ParametricRealField),
      piecewise_function_from_breakpoints_and_values(JJ,
                                                     [ pi_v_prime + (y - v_prime) * s[1] for y in JJ ],
                                                     field=ParametricRealField, merge=False),
      piecewise_function_from_breakpoints_and_values(K,
                                                     [ pi_w_prime + (z - w_prime) * s[2] for z in K ],
                                                     field=ParametricRealField)]

d1_plus, d2_plus, d3_plus = d_plus  = [ inv_mq * (s_i - s_m) / (s_p - s_m) for s_i in s ]
d1_minus, d2_minus, d3_minus = d_minus = [ inv_mq * (s_p - s_i) / (s_p - s_m) for s_i in s ]
assert all(d_plus[i] > 0 and d_minus[i] > 0 and d_plus[i] + d_minus[i] == inv_mq for i in range(3))

KK.freeze()

d_plusminus = [ None, d_plus, d_minus ]   # index by 1, -1.
s_plusminus = [ None, s_p, s_m ]

# Trivial consequences of the slope inequalities
assert d2_plus <= d1_plus and d2_plus <= d3_plus
assert d2_minus >= d1_minus and d2_minus >= d3_minus

signs = (+1, +1, +1)

## with K.off_the_record():
##     plot_fun_IJK(pi, phi).show(figsize=[8,2.5])

phi = construct_phi(pi, signs)

(x, y, z) = [ phi_i.end_points()[1] for phi_i in phi ]

P12 = (x, y)
P13 = (x, z - x)
P23 = (z - y, y)

plot_case(phi, [P12, P13, P23]).show(legend_title='Signs {}'.format([print_sign(s) for s in signs]))

assert P23[0] > I[1]             # outside F

assert P13[1] >= y
with KK.temporary_assumptions():
    with KK.unfrozen():
        assert P13[1] <= J[1]
    #phi[1](P13[1])
    assert phi[0](P13[0]) - phi[0](u_prime) == d1_plus * s_p
    assert phi[2](P13[0] + P13[1]) - phi[2](u_prime + v_prime) == d3_minus * (- s_m)

    assert delta_IJK(phi, P13) >= 0

#assert P12[0] + P12[1] < K[0]  # unknown
with KK.temporary_assumptions():
    with KK.unfrozen():
        assert P12[0] + P12[1] >= K[0]
    assert delta_IJK(phi, P12) >= 0

####
signs = (+1, -1, -1)
phi = construct_phi(pi, signs)

(x, y, z) = [ phi_i.end_points()[1] for phi_i in phi ]

P12 = (x, y, x + y)
P13 = (x, z - x, z)
P23 = (z - y, y, z)

plot_case(phi, [P12, P13, P23]).show(legend_title='Signs {}'.format([print_sign(s) for s in signs]))

# There are 2 combinatorial types of the 2d diagram.
with KK.temporary_assumptions():
    with KK.unfrozen():
        assert P23[0] > x
    assert P13[1] > y
    assert P12[2] < z

    assert P12 in F
    assert delta_IJK(phi, P12) >= 0

    assert P23 in F
    assert delta_IJK(phi, P23) >= 0

    with KK.temporary_assumptions():
        with KK.unfrozen():
            assert P13 in F #assert P13[1] <= J[1]
        assert delta_IJK(phi, P13) >= 0

with KK.temporary_assumptions():
    with KK.changed_values(s_2=-1/5):
        plot_case(phi, [P12, P13, P23]).show(legend_title='Signs {}'.format([print_sign(s) for s in signs]))
        with KK.unfrozen():
            assert P23[0] < x
        assert P13[1] < y
        assert P12[2] > z

        assert P12 in F
        assert delta_IJK(phi, P12) >= 0

        assert P13 in F
        assert delta_IJK(phi, P13) >= 0

        assert P23 in F
        assert delta_IJK(phi, P23) >= 0

####
signs = (-1, +1, -1)
phi = construct_phi(pi, signs)

(x, y, z) = [ phi_i.end_points()[1] for phi_i in phi ]

P12 = (x, y, x + y)
P13 = (x, z - x, z)
P23 = (z - y, y, z)

assert P12 not in F

with KK.changed_values(s_3=1/3):
    with KK.temporary_assumptions():
        with KK.unfrozen():
            assert P23 in F
        assert delta_IJK(phi, P23) >= 0
    with KK.temporary_assumptions():
        with KK.unfrozen():
            assert P13 in F
        assert delta_IJK(phi, P13) >= 0
