sage: logging.disable(logging.info)
sage: h = drlm_not_extreme_1()
sage: hh = lift_until_extreme(h) # same for other igp.perturbation_style
5/7
1
#sage: hh = lift_extreme_function_for_finite_group_to_infinite_group(h)
#sage: len(hh)
#2
sage: hh = lift_on_uncovered_components(h)
sage: len(hh)
4

sage: h = example7slopecoarse2()
sage: hh = lift_until_extreme(h) # same for other igp.perturbation_style
2/3
1
sage: igp.perturbation_style = 'slopes_proportional_to_limiting_slopes_for_positive_ epsilon'
sage: hh = lift_until_extreme(h)
2/3
1
sage: igp.perturbation_style = 'slopes_proportional_to_limiting_slopes_for_negative_epsilon'
sage: hh = lift_until_extreme(h)
2/3
1
#sage: hh = lift_extreme_function_for_finite_group_to_infinite_group(h)
#sage: len(hh)
#12
sage: hh = lift_on_uncovered_components(h)
sage: len(hh)
12

sage: bkpts = [0, 1/13, 3/13, 7/26, 4/13,5/13,21/52,23/52,6/13,8/13,33/52,35/52,9/13,10/13,21/26,11/13,1]
sage: values = [0,1,3/14,5/7,3/4,5/14,55/112,33/112,3/7,4/7,79/112,57/112,9/14,1/4,2/7,11/14,0]
sage: h = piecewise_function_from_breakpoints_and_values(bkpts, values)
sage: hh = lift_until_extreme(h)
12/13
25/26
51/52
103/104
207/208
415/416
831/832
1663/1664
3327/3328
6655/6656
13311/13312
26623/26624
53247/53248
106495/106496
212991/212992
425983/425984
851967/851968
sage: igp.perturbation_style = 'slopes_proportional_to_limiting_slopes_for_positive_
....: epsilon'
sage: hh = lift_until_extreme(h)
12/13
25/26
1sage: igp.perturbation_style = 'slopes_proportional_to_limiting_slopes_for_negative_
....: epsilon'
sage: hh = lift_until_extreme(h)
12/13
25/26
1

#sage: hh = lift_extreme_function_for_finite_group_to_infinite_group(h)
#sage: len(hh)
#2
#sage: extremality_test(hh[0])
sage: hh = lift_on_uncovered_components(h)   # long time. 30 mins?
sage: len(hh)                                # long time
128
sage: extremality_test(hh[0])
True

        
# The following bug example of lift_extreme_function_for_finite_group_to_infinite_group()

        sage: h = FastPiecewise([[(QQ(0), 1/18), FastLinearFunction(QQ(18), QQ(0))], [(1/18, 1/9), FastLinearFunction(-126/11, 18/11)], [(1/9, 1/6), FastLinearFunction(18/55, 18/55)], [(1/6, 2/9), FastLinearFunction(342/55, -36/55)], [(2/9, 5/18), FastLinearFunction(-126/11, 36/11)], [(5/18, 1/3), FastLinearFunction(666/55, -36/11)], [(1/3, 7/18), FastLinearFunction(-306/55, 144/55)], [(7/18, 4/9), FastLinearFunction(18/55, 18/55)], [(4/9, 1/2), FastLinearFunction(342/55, -126/55)], [(1/2, 5/9), FastLinearFunction(-126/11, 72/11)], [(5/9, 11/18), FastLinearFunction(342/55, -36/11)], [(11/18, 2/3), FastLinearFunction(18/55, 18/55)], [(2/3, 13/18), FastLinearFunction(-306/55, 234/55)], [(13/18, 7/9), FastLinearFunction(666/55, -468/55)], [(7/9, 5/6), FastLinearFunction(-126/11, 108/11)], [(5/6, 8/9), FastLinearFunction(342/55, -54/11)], [(8/9, 17/18), FastLinearFunction(18/55, 18/55)], [(17/18, QQ(1)), FastLinearFunction(-126/11, 126/11)]])
        sage: lift_extreme_function_for_finite_group_to_infinite_group(h)
        []

#This is also a bug example for lift_until_extreme(h).

# It produces an infinite sequence of lifting using igp.perturbation_style = default 'use_pwl_template'. Print covered interval length in each round.

sage: hh = lift_until_extreme(h)
1/3
17/18
35/36
71/72
143/144
287/288
575/576
1151/1152
2303/2304
4607/4608
9215/9216
18431/18432
36863/36864
73727/73728
147455/147456
294911/294912
589823/589824
1179647/1179648
2359295/2359296
4718591/4718592
9437183/9437184
18874367/18874368
37748735/37748736
75497471/75497472
150994943/150994944
301989887/301989888
603979775/603979776
1207959551/1207959552
2415919103/2415919104
4831838207/4831838208
9663676415/9663676416
19327352831/19327352832
38654705663/38654705664
77309411327/77309411328
154618822655/154618822656
309237645311/309237645312
618475290623/618475290624
1236950581247/1236950581248

# Lifting does not terminiate in 25 rounds with igp.perturbation_style = 'slopes_proportional_to_limiting_slopes_for_positive_epsilon'

sage: igp.perturbation_style ='slopes_proportional_to_limiting_slopes_for_positive_epsilon'
sage: hh = lift_until_extreme(h)
1/3
29/45
173/225
371/450
793/900
26473/28800
7603/8064
89021/94080
9371137/9878400
4377053/4609920
65673767/69148800
357607123/376476800
22530629281/23718038400
157720837763/166026268800
1104061216129/1162183881600
7728466045067/8135287171200
54099359188561/56947010198400
378695773492403/398629071388800
2650871129070649/2790403499721600
18556099921578827/19532824498051200
129892705253730241/136729771486358400
909248953680998243/957108400404508800
6364742725475349769/6699758802831561600
44553199225439939387/46898311619820931200
311872395015391858321/328288181338746518400
2183106766411629475283/2298017269371225628800
^C---------------------------------------------------------------------------

# Same thing for igp.perturbation_style = 'slopes_proportional_to_limiting_slopes_for_negative_epsilon'
sage: igp.perturbation_style ='slopes_proportional_to_limiting_slopes_for_negative_epsilon'
sage: hh = lift_until_extreme(h)
1/3
2/3
206/225
659/675
44/45
1321/1350
881/900
5287/5400
47/48
21151/21600
14101/14400
84607/86400
11281/11520
338431/345600
25069/25600
1353727/1382400
180497/184320
5414911/5529600
3609941/3686400
21659647/22118400
962651/983040
86638591/88473600
57759061/58982400
346554367/353894400
46207249/47185920
1386217471/1415577600
308048327/314572800
5544869887/5662310400
147863197/150994944
22179479551/22649241600
14786319701/15099494400
88717918207/90596966400
1314339529/1342177280

#This bug example is now solved by generate_extreme_lifted_function_equiv()

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
# Drawback: too many vertices. SLOW!

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
