# The following bug example of lift_extreme_function_for_finite_group_to_infinite_group()
# is now solved by generate_extreme_lifted_function_equiv()

h = FastPiecewise([[(QQ(0), 1/18), FastLinearFunction(QQ(18), QQ(0))], [(1/18, 1/9), FastLinearFunction(-126/11, 18/11)], [(1/9, 1/6), FastLinearFunction(18/55, 18/55)], [(1/6, 2/9), FastLinearFunction(342/55, -36/55)], [(2/9, 5/18), FastLinearFunction(-126/11, 36/11)], [(5/18, 1/3), FastLinearFunction(666/55, -36/11)], [(1/3, 7/18), FastLinearFunction(-306/55, 144/55)], [(7/18, 4/9), FastLinearFunction(18/55, 18/55)], [(4/9, 1/2), FastLinearFunction(342/55, -126/55)], [(1/2, 5/9), FastLinearFunction(-126/11, 72/11)], [(5/9, 11/18), FastLinearFunction(342/55, -36/11)], [(11/18, 2/3), FastLinearFunction(18/55, 18/55)], [(2/3, 13/18), FastLinearFunction(-306/55, 234/55)], [(13/18, 7/9), FastLinearFunction(666/55, -468/55)], [(7/9, 5/6), FastLinearFunction(-126/11, 108/11)], [(5/6, 8/9), FastLinearFunction(342/55, -54/11)], [(8/9, 17/18), FastLinearFunction(18/55, 18/55)], [(17/18, QQ(1)), FastLinearFunction(-126/11, 126/11)]])



# However, new bug example appears: h is not extreme, it has pert, but this pert is not a vertex of perturbation_lim_slopes_mip.  <--- Bug example solved by setting for example perturbation_components = [fn._stability_orbits[0], fn._stability_orbits[3], fn._stability_orbits[5]]. When perturbing all uncovered intervals (perturbation_components = fn._stability_orbits), setting n>1 and "if (sa == 0 and sb != 0) or (sa != 0 and sb == 0):" fails to give subadditive function, because not only slope values but also slope lengths are related, maybe another MIP for lengths?


h = FastPiecewise([[(QQ(0), 1/18), FastLinearFunction(QQ(18), QQ(0))], [(1/18, 1/9), FastLinearFunction(-QQ(14), 16/9)], [(1/9, 2/9), FastLinearFunction(QQ(2), QQ(0))], [(2/9, 5/18), FastLinearFunction(-QQ(2), 8/9)], [(5/18, 1/3), FastLinearFunction(QQ(6), -4/3)], [(1/3, 7/18), FastLinearFunction(-QQ(6), 8/3)], [(7/18, 4/9), FastLinearFunction(QQ(6), -QQ(2))], [(4/9, 1/2), FastLinearFunction(-QQ(6), 10/3)], [(1/2, 5/9), FastLinearFunction(QQ(6), -8/3)], [(5/9, 11/18), FastLinearFunction(-QQ(6), QQ(4))], [(11/18, 2/3), FastLinearFunction(QQ(6), -10/3)], [(2/3, 13/18), FastLinearFunction(-QQ(6), 14/3)], [(13/18, 7/9), FastLinearFunction(QQ(6), -QQ(4))], [(7/9, 5/6), FastLinearFunction(-QQ(2), 20/9)], [(5/6, 17/18), FastLinearFunction(QQ(2), -10/9)], [(17/18, QQ(1)), FastLinearFunction(-QQ(14), QQ(14))]])

pert = FastPiecewise([[(QQ(0), 1/9), FastLinearFunction(QQ(0), QQ(0))], [left_open_interval(1/9, 5/36), FastLinearFunction(-QQ(4), 4/9)], [left_open_interval(5/36, 7/36), FastLinearFunction(QQ(4), -2/3)], [left_open_interval(7/36, 2/9), FastLinearFunction(-QQ(4), 8/9)], [left_open_interval(2/9, 1/3), FastLinearFunction(QQ(0), QQ(0))], [left_open_interval(1/3, 10/27), FastLinearFunction(QQ(4), -4/3)], [left_open_interval(10/27, 7/18), FastLinearFunction(-QQ(8), 28/9)], [left_open_interval(7/18, 4/9), FastLinearFunction(QQ(0), QQ(0))], [left_open_interval(4/9, 7/15), FastLinearFunction(QQ(12), -16/3)], [left_open_interval(7/15, 1/2), FastLinearFunction(-QQ(8), QQ(4))], [left_open_interval(1/2, 5/9), FastLinearFunction(QQ(0), QQ(0))], [left_open_interval(5/9, 53/90), FastLinearFunction(-QQ(8), 40/9)], [left_open_interval(53/90, 11/18), FastLinearFunction(QQ(12), -22/3)], [left_open_interval(11/18, 2/3), FastLinearFunction(QQ(0), QQ(0))], [left_open_interval(2/3, 37/54), FastLinearFunction(-QQ(8), 16/3)], [left_open_interval(37/54, 13/18), FastLinearFunction(QQ(4), -26/9)], [left_open_interval(13/18, 5/6), FastLinearFunction(QQ(0), QQ(0))], [left_open_interval(5/6, 31/36), FastLinearFunction(-QQ(4), 10/3)], [left_open_interval(31/36, 11/12), FastLinearFunction(QQ(4), -32/9)], [left_open_interval(11/12, 17/18), FastLinearFunction(-QQ(4), 34/9)], [left_open_interval(17/18, QQ(1)), FastLinearFunction(QQ(0), QQ(0))]])

# Bug example 2017-8-14 [lifting 4dc8d6e]: generate_extreme_lifted_function_equiv() can not lift the following function h with "if sa * sb < 0." Solved by changing to  "sa * sb != 0".

h = FastPiecewise([[(QQ(0), 1/18), FastLinearFunction(QQ(18), QQ(0))], [(1/18, 1/6), FastLinearFunction(-414/67, 90/67)], [(1/6, 2/9), FastLinearFunction(-90/67, 36/67)], [(2/9, 1/3), FastLinearFunction(234/67, -36/67)], [(1/3, 7/18), FastLinearFunction(-414/67, 180/67)], [(7/18, 1/2), FastLinearFunction(234/67, -72/67)], [(1/2, 5/9), FastLinearFunction(-414/67, 252/67)], [(5/9, 2/3), FastLinearFunction(234/67, -108/67)], [(2/3, 13/18), FastLinearFunction(-414/67, 324/67)], [(13/18, 5/6), FastLinearFunction(234/67, -144/67)], [(5/6, 8/9), FastLinearFunction(-90/67, 126/67)], [(8/9, QQ(1)), FastLinearFunction(-414/67, 414/67)]])

# infinite lifting sequence by lift_until_extreme with igp.perturbation_style = 'use_pwl_template'
# lift to extreme in 2 rounds (3 equiv perturbations) by lift_until_extreme with igp.perturbation_style = 'slopes_proportional_to_limiting_slopes_for_positive_epsilon'

perturbation = FastPiecewise([[(QQ(0), 1/6), FastLinearFunction(QQ(0), QQ(0))], [left_open_interval(1/6, 8/45), FastLinearFunction(1296/67, -216/67)], [left_open_interval(8/45, 2/9), FastLinearFunction(-324/67, 72/67)], [left_open_interval(2/9, 5/18), FastLinearFunction(QQ(0), QQ(0))], [left_open_interval(5/18, 3/10), FastLinearFunction(-648/67, 180/67)], [left_open_interval(3/10, 29/90), FastLinearFunction(972/67, -306/67)], [left_open_interval(29/90, 1/3), FastLinearFunction(-648/67, 216/67)], [left_open_interval(1/3, 13/18), FastLinearFunction(QQ(0), QQ(0))], [left_open_interval(13/18, 11/15), FastLinearFunction(-648/67, 468/67)], [left_open_interval(11/15, 34/45), FastLinearFunction(972/67, -720/67)], [left_open_interval(34/45, 7/9), FastLinearFunction(-648/67, 504/67)], [left_open_interval(7/9, 5/6), FastLinearFunction(QQ(0), QQ(0))], [left_open_interval(5/6, 79/90), FastLinearFunction(-324/67, 270/67)], [left_open_interval(79/90, 8/9), FastLinearFunction(1296/67, -1152/67)], [left_open_interval(8/9, QQ(1)), FastLinearFunction(QQ(0), QQ(0))]])

#lift to extreme in 5 rounds (many equiv perturbations) by lift_until_extreme with igp.perturbation_style = 'slopes_proportional_to_limiting_slopes_for_negative_epsilon'

perturbation = FastPiecewise([[(QQ(0), 1/6), FastLinearFunction(QQ(0), QQ(0))], [left_open_interval(1/6, 7/36), FastLinearFunction(-324/67, 54/67)], [left_open_interval(7/36, 37/180), FastLinearFunction(1296/67, -261/67)], [left_open_interval(37/180, 2/9), FastLinearFunction(-324/67, 72/67)], [left_open_interval(2/9, 5/18), FastLinearFunction(QQ(0), QQ(0))], [left_open_interval(5/18, 13/45), FastLinearFunction(972/67, -270/67)], [left_open_interval(13/45, 14/45), FastLinearFunction(-648/67, 198/67)], [left_open_interval(14/45, 19/60), FastLinearFunction(972/67, -306/67)], [left_open_interval(19/60, 59/180), FastLinearFunction(-648/67, 207/67)], [left_open_interval(59/180, 1/3), FastLinearFunction(972/67, -324/67)], [left_open_interval(1/3, 13/18), FastLinearFunction(QQ(0), QQ(0))], [left_open_interval(13/18, 131/180), FastLinearFunction(972/67, -702/67)], [left_open_interval(131/180, 133/180), FastLinearFunction(-648/67, 477/67)], [left_open_interval(133/180, 67/90), FastLinearFunction(972/67, -720/67)], [left_open_interval(67/90, 23/30), FastLinearFunction(-648/67, 486/67)], [left_open_interval(23/30, 7/9), FastLinearFunction(972/67, -756/67)], [left_open_interval(7/9, 5/6), FastLinearFunction(QQ(0), QQ(0))], [left_open_interval(5/6, 17/20), FastLinearFunction(-324/67, 270/67)], [left_open_interval(17/20, 31/36), FastLinearFunction(1296/67, -1107/67)], [left_open_interval(31/36, 8/9), FastLinearFunction(-324/67, 288/67)], [left_open_interval(8/9, QQ(1)), FastLinearFunction(QQ(0), QQ(0))]])

# Observe that sa * sb > 0

#lift_extreme_function_for_finite_group_to_infinite_group() fails.

# [lifting fcc6896] fix the bug example by checking "sa * sb != 0" only.

# Bug example 2017-8-14 [lifting 21b1551]

# "The lifted function has uncovered intervals" <-- normal, since we set slopes on some uncovered intervals to 0.
h = FastPiecewise([[(QQ(0), 1/18), FastLinearFunction(QQ(18), QQ(0))], [(1/18, 1/9), FastLinearFunction(-414/31, 54/31)], [(1/9, 1/6), FastLinearFunction(234/31, -18/31)], [(1/6, 2/9), FastLinearFunction(-90/31, 36/31)], [(2/9, 1/3), FastLinearFunction(72/31, QQ(0))], [(1/3, 1/2), FastLinearFunction(-90/31, 54/31)], [(1/2, 5/9), FastLinearFunction(234/31, -108/31)], [(5/9, 13/18), FastLinearFunction(-90/31, 72/31)], [(13/18, 5/6), FastLinearFunction(72/31, -45/31)], [(5/6, 8/9), FastLinearFunction(-90/31, 90/31)], [(8/9, 17/18), FastLinearFunction(234/31, -198/31)], [(17/18, QQ(1)), FastLinearFunction(-414/31, 414/31)]])

perturbation = FastPiecewise([[(0, 1/9), FastLinearFunction(QQ(0), QQ(0))], [left_open_interval(1/9, 5/36), FastLinearFunction(324/31, -36/31)], [left_open_interval(5/36, 1/6), FastLinearFunction(-324/31, 54/31)], [left_open_interval(1/6, 1/2), FastLinearFunction(QQ(0), QQ(0))], [left_open_interval(1/2, 55/108), FastLinearFunction(324/31, -162/31)], [left_open_interval(55/108, 14/27), FastLinearFunction(-648/31, 333/31)], [left_open_interval(14/27, 29/54), FastLinearFunction(324/31, -171/31)], [left_open_interval(29/54, 59/108), FastLinearFunction(-648/31, 351/31)], [left_open_interval(59/108, 5/9), FastLinearFunction(324/31, -180/31)], [left_open_interval(5/9, 8/9), FastLinearFunction(QQ(0), QQ(0))], [left_open_interval(8/9, 11/12), FastLinearFunction(-324/31, 288/31)], [left_open_interval(11/12, 17/18), FastLinearFunction(324/31, -306/31)], [left_open_interval(17/18, QQ(1)), FastLinearFunction(QQ(0), QQ(0))]])


#################

q = 18
for f in range(1, q/2+1):
    print f
    for h in generate_extreme_functions_for_finite_group(q,f):
        if extremality_test(h):
            continue
        gen = generate_extreme_lifted_function_equiv(h)
        try:
            hl = gen.next()
        except:
            print "can not lift the function", sage_input(h)
            raise ValueError
