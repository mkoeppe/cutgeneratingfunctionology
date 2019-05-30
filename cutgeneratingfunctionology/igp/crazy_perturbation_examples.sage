from six.moves import range
def kzh_discontinuous_bhk_irrational(f=4/5, d1=3/5, d2=5/40, a0=19/100, delta_ratio=sqrt(2)/3, bb=1/1000, c2=0, y1=1/10, y2=1/50, field=None):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: hmin = kzh_discontinuous_bhk_irrational(f=4/5, d1=3/5, d2=5/40, a0=19/100, delta_ratio=sqrt(2)/3, bb=19/23998, c2=5/11999, y1=185/1846, y2=240/11999, field=None)
        sage: minimality_test(hmin)
        True
    """
    [f, d1, d2, a0, delta_ratio, bb, c2, y1, y2] = nice_field_values([f, d1, d2, a0, delta_ratio, bb, c2, y1, y2])
    rnf = f.parent().fraction_field()
    d3 = f - d1 - d2
    c3 = -rnf(1)/(1-f)
    c1 = (1-d2*c2-d3*c3+2*(y1+y2))/d1
    d21 = d2 / 2
    d31 = ((c1 - c2) * d21 + 2*y2)  / (c1 - c3)
    d11 = a0 - d31
    d13 = a0 - d21
    d12 = (d1 - d13)/2 - d11
    d32 = d3/2 - d31
    d32a = delta_ratio * d32
    d32b = d32 - d32a
    d12a = (c2 -c3) / (c1 - c2) * d32a
    d12b = (c2 -c3) / (c1 - c2) * d32b
    d12c = d12 - d12a - d12b
    b = (y1 + bb) / (c2 - c3)

    slopes_pwl = [c1, c3]*3 + [c1, c2, c1, c2, c1] + [c3, c1]*3 + [c3]
    intervals_pwl_left = [d11, d31, d12a, d32a, d12b, d32b, d12c, d21]
    interval_lengths_pwl = intervals_pwl_left + [d13] + intervals_pwl_left[::-1] + [1-f]
    h_pwl = piecewise_function_from_interval_lengths_and_slopes(interval_lengths_pwl, slopes_pwl, field=rnf) #, merge=False)
    j = b*(c1-c3)
    bkpts = [sum(interval_lengths_pwl[0:i]) for i in range(19)]
    disc_1 = bkpts[7];
    disc_2 = bkpts[8]
    y_shift_pieces = [closed_piece((rnf(0), 0), (disc_1, 0)),\
                      open_piece((disc_1, -y1), (disc_2, -y1)), \
                      closed_piece((disc_2, -y1-y2), (f-disc_2, -y1-y2)), \
                      open_piece((f-disc_2, -y1-2*y2), (f-disc_1, -y1-2*y2)), \
                      closed_piece((f-disc_1, -2*y1-2*y2), (rnf(1), -2*y1-2*y2))] 
    y_shift = FastPiecewise(y_shift_pieces)
    b_shift_pieces = [singleton_piece(rnf(0), 0), \
                      open_piece((rnf(0),j), (b, 0)), \
                      singleton_piece(b, j),
                      open_piece((b, 0), (f-b, 0)), \
                      singleton_piece(f-b, -j),
                      open_piece((f-b, 0), (f, -j)), \
                      singleton_piece(f,0), \
                      open_piece((f, -j), (f+b, 0)), \
                      singleton_piece(f+b, -j),
                      open_piece((f+b, 0), (1-b, 0)), \
                      singleton_piece(1-b, j),
                      open_piece((1-b, 0), (rnf(1), j)), \
                      singleton_piece(rnf(1), 0) ]
    b_shift = FastPiecewise(b_shift_pieces)
    h_temp = h_pwl + y_shift + b_shift
    h_below = piecewise_function_from_breakpoints_and_slopes([rnf(0)]+bkpts[1:8]+bkpts[10:17]+[rnf(1)], [0, c3, c1, c3, c1, c3, c1, 0, c1, c3, c1, c3, c1, c3, 0])
    zigzag_lengths =  [bkpts[2]-bkpts[1]-b, bkpts[3]-bkpts[2], bkpts[4]-bkpts[3], bkpts[5]-bkpts[4], bkpts[6]-bkpts[5], b-bkpts[3]+bkpts[2]-bkpts[5]+bkpts[4], b-bkpts[4]+bkpts[3]-bkpts[6]+bkpts[5], bkpts[3]-bkpts[2], bkpts[4]-bkpts[3], bkpts[5]-bkpts[4], bkpts[6]-bkpts[5], bkpts[7]-bkpts[6]-b]
    h_above = piecewise_function_from_interval_lengths_and_slopes([bkpts[1]] + zigzag_lengths + [bkpts[10]-bkpts[7]] + zigzag_lengths[::-1]+[rnf(1)-bkpts[16]], [0] + [c3, c1] * 6 + [0] + [c1, c3] * 6 + [0])
    a1 = a0+d12a+d32a
    a2 = a0+d12a+d32a+d12b+d32b
    move_pieces = [right_open_piece((rnf(0), 0), (a0, 0)), \
                   singleton_piece(a0, h_below(a0)-h_above(a0)), \
                   open_piece((a0, 0), (a1, 0)), \
                   singleton_piece(a1, h_below(a1)-h_above(a1)), \
                   open_piece((a1, 0), (a2, 0)), \
                   singleton_piece(a2, h_below(a2)-h_above(a2)), \
                   open_piece((a2, 0), (f-a2, 0)), \
                   singleton_piece(f-a2, h_below(f-a2)-h_above(f-a2)), \
                   open_piece((f-a2, 0), (f-a1, 0)), \
                   singleton_piece(f-a1, h_below(f-a1)-h_above(f-a1)), \
                   open_piece((f-a1, 0), (f-a0, 0)), \
                   singleton_piece(f-a0, h_below(f-a0)-h_above(f-a0)), \
                   left_open_piece((f-a0, 0), (rnf(1), 0))]
    move_shift = FastPiecewise(move_pieces)
    return h_temp - h_below + h_above + move_shift

def kzh_minimal_has_only_crazy_perturbation_1():
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = kzh_minimal_has_only_crazy_perturbation_1()

    On one hand, normal perturbation does not exist::

        sage: extremality_test(h, crazy_perturbations=False)
        True

    On the other hand, there exists crazy perturbation, such as the one we construct below::

        sage: bkpts = h.end_points()
        sage: t1 = bkpts[10]-bkpts[6]
        sage: t2 = bkpts[13]-bkpts[6]
        sage: f = bkpts[37]
        sage: ucl = bkpts[17]
        sage: ucr = bkpts[18]
        sage: generators = [t1, t2]
        sage: pwl = piecewise_function_from_breakpoints_and_slopes([0,1],[0])
        sage: crazy_piece_1 = CrazyPiece((ucl, ucr), generators, [(ucl, 1), (ucr, -1)])
        sage: crazy_piece_2 = CrazyPiece((f-ucr, f-ucl), generators, [(f-ucr, 1), (f-ucl, -1)])
        sage: cp = PiecewiseCrazyFunction(pwl, [crazy_piece_1, crazy_piece_2])
        sage: cp == kzh_minimal_has_only_crazy_perturbation_1_perturbation()  # not tested - equality is not implemented
        True

    This crazy perturbation is valid, since it has positive epsilon::

        sage: find_epsilon_for_crazy_perturbation(h, cp)
        0.0003958663221935161?

    Therefore, the function ``kzh_minimal_has_only_crazy_perturbation_1()`` is not extreme.

    Note:

        This example is obtained by the following code::

            sage: hmin = kzh_discontinuous_bhk_irrational(f=4/5, d1=3/5, d2=5/40, a0=19/100, delta_ratio=sqrt(2)/3, bb=19/23998, c2=5/11999, y1=185/1846, y2=240/11999, field=None) # long time
            sage: hlift = lift_until_extreme(hmin, use_all_perturbations=False, use_largest_absolute_epsilon=False) # long time
            sage: h == hlift  # long time
            True

        Without knowing the input values of ``kzh_discontinuous_bhk_irrational()`` that gives a minimal valid 'hmin', one could use the default values to construct a non-minimal discontinuous function, and then lift the function to minimal, as follows::

            sage: h_org = kzh_discontinuous_bhk_irrational() # long time
            sage: hmin_from_org = lift(h_org, phase_1=True, use_all_perturbations=False, use_largest_absolute_epsilon=False) # long time
            sage: hmin ==  hmin_from_org # long time
            True

        'hmin' is minimal but not extreme. The solution space of the finite dimensional test has dimension 5::

            sage: finite_dimensional_extremality_test(hmin,show_all_perturbations=True) # long time
            False
            sage: perturbs = hmin._perturbations # long time
            sage: len(perturbs) # long time
            5

        A more general way of lifting 'hmin' is to call the following generator. (Use Sage Polyhedron if 'use_polyhedron' is set to ``False``. Use LP if with random objective function if 'use_polyhedron' is set to ``True``.) Unfortunately, this method is too slow. With ``use_polyhedron=False``, it takes 15-20 mins to find a lifted function::

            sage: gen = generate_lifted_functions(hmin, perturbs=perturbs, use_polyhedron=False) # not tested
            sage: h = next(gen) # not tested
    """
    # The following numbers do not correspond to h in docstring. To check.
    [sqrt2, o] = nice_field_values([sqrt(2), 1])
    n = 41
    rnf = o.parent().fraction_field()
    bkpt_rat = [0, 101/5000, 60153/369200, 849/5000, 849/5000, 849/5000, 19/100, 281986521/1490645000, 40294/201875, 36999/184600, 19/100, 1051/5000, 1051/5000, 14199/64600, 1051/5000, 342208579/1490645000, 193799/807500, 219/800, 269/800, 371/800, 421/800, 452201/807500, 850307421/1490645000, 2949/5000, 37481/64600, 2949/5000, 2949/5000, 61/100, 110681/184600, 121206/201875, 910529479/1490645000, 61/100, 3151/5000, 3151/5000, 3151/5000, 235207/369200, 3899/5000, 4/5, 4101/5000, 4899/5000, 1]
    bkpt_irr = [0, 0, 0, 0, 1925/298129, 77/7752, 0, 77/22152, 0, 0, 77/7752, 0, 1925/298129, 0, 77/7752, 77/22152, 0, 0, 0, 0, 0, 0, -77/22152, -77/7752, 0, -1925/298129, 0, -77/7752, 0, 0, -77/22152, 0, -77/7752, -1925/298129, 0, 0, 0, 0, 0, 0, 0]
    value_rat = [0, 2727/13000, 421071/959920, 4851099/11999000, 4851099/11999000, 4851099/11999000, 18196/59995, 10467633/22933000, 795836841/1937838500, 975607/2399800, 18196/59995, 4291761/11999000, 4291761/11999000, 50943/167960, 4291761/11999000, 122181831/298129000, 187742/524875, 933/2080, 683/2080, 1397/2080, 1147/2080, 337133/524875, 175947169/298129000, 7707239/11999000, 117017/167960, 7707239/11999000, 7707239/11999000, 41799/59995, 1424193/2399800, 1142001659/1937838500, 12465367/22933000, 41799/59995, 7147901/11999000, 7147901/11999000, 7147901/11999000, 538849/959920, 10273/13000, 1, 9667/13000, 3333/13000, 0]
    value_irr = [0, 0, 0, -1925/71994, 67375/3875677, 2695/100776, 0, -385/22152, 0, 0, 385/93016248, -1925/71994, 67375/3875677, 0, 2695/100776, -385/22152, 0, 0, 0, 0, 0, 0, 385/22152, -2695/100776, 0, -67375/3875677, 1925/71994, -385/93016248, 0, 0, 385/22152, 0, -2695/100776, -67375/3875677, 1925/71994, 0, 0, 0, 0, 0, 0]
    lim_left_rat = [101/650, 707/13000, 421071/959920, 4851099/11999000, 4851099/11999000, 4851099/11999000, 275183/599950, 10467633/22933000, 848837/2099500, 975607/2399800, 275183/599950, 4291761/11999000, 4291761/11999000, 240046061/775135400, 4291761/11999000, 122181831/298129000, 187742/524875, 933/2080, 668809/1919840, 1397/2080, 96237/147680, 337133/524875, 175947169/298129000, 7707239/11999000, 535089339/775135400, 7707239/11999000, 7707239/11999000, 324767/599950, 1424193/2399800, 1250663/2099500, 12465367/22933000, 324767/599950, 7147901/11999000, 7147901/11999000, 7147901/11999000, 538849/959920, 12293/13000, 549/650, 899/1000, 101/1000, 101/650]
    lim_left_irr = [0, 0, 0, 0, 67375/3875677, 385/93016248, -1925/71994, -385/22152, 0, 0, -385/7752, 0, 67375/3875677, 192500/3875677, 385/93016248, -385/22152, 0, 0, 0, 0, 0, 0, 385/22152, -385/93016248, -192500/3875677, -67375/3875677, 0, 385/7752, 0, 0, 385/22152, 1925/71994, -385/93016248, -67375/3875677, 0, 0, 0, 0, 0, 0, 0]
    lim_right_rat = [101/650, 707/13000, 421071/959920, 4851099/11999000, 4851099/11999000, 4851099/11999000, 275183/599950, 10467633/22933000, 848837/2099500, 975607/2399800, 275183/599950, 4291761/11999000, 4291761/11999000, 240046061/775135400, 4291761/11999000, 122181831/298129000, 187742/524875, 51443/147680, 683/2080, 1251031/1919840, 1147/2080, 337133/524875, 175947169/298129000, 7707239/11999000, 535089339/775135400, 7707239/11999000, 7707239/11999000, 324767/599950, 1424193/2399800, 1250663/2099500, 12465367/22933000, 324767/599950, 7147901/11999000, 7147901/11999000, 7147901/11999000, 538849/959920, 12293/13000, 549/650, 899/1000, 101/1000, 101/650]
    lim_right_irr = [0, 0, 0, 0, 67375/3875677, 385/93016248, -1925/71994, -385/22152, 0, 0, -385/7752, 0, 67375/3875677, 192500/3875677, 385/93016248, -385/22152, 0, 0, 0, 0, 0, 0, 385/22152, -385/93016248, -192500/3875677, -67375/3875677, 0, 385/7752, 0, 0, 385/22152, 1925/71994, -385/93016248, -67375/3875677, 0, 0, 0, 0, 0, 0, 0]
    pieces = [singleton_piece(rnf(0), rnf(0))]
    for i in range(1,n):
        pieces.append(open_piece((bkpt_rat[i-1]*o+bkpt_irr[i-1]*sqrt2, lim_right_rat[i-1]*o+lim_right_irr[i-1]*sqrt2),(bkpt_rat[i]*o+bkpt_irr[i]*sqrt2, lim_left_rat[i]*o+lim_left_irr[i]*sqrt2)))
        pieces.append(singleton_piece(bkpt_rat[i]*o+bkpt_irr[i]*sqrt2, value_rat[i]*o+value_irr[i]*sqrt2))
    return FastPiecewise(pieces)

def kzh_minimal_has_only_crazy_perturbation_1_perturbation():
    r"""
    A crazy perturbation for ``kzh_minimal_has_only_crazy_perturbation_1``.
    """
    h = kzh_minimal_has_only_crazy_perturbation_1()
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
    return PiecewiseCrazyFunction(pwl, [crazy_piece_1, crazy_piece_2])

def kzh_minimal_has_only_crazy_perturbation_1_check_subadditivity_slacks():
    def generate_all_faces(fn):
        zero_fn = FastPiecewise([(I, FastLinearFunction(0, 0)) for I, f in fn.list()], merge=False)
        return generate_maximal_additive_faces_general(zero_fn)
    def number_of_projections_intersecting(F, real_set):
        return len([I for I in F.minimal_triple
                    if not real_set.is_disjoint_from(realset_from_interval(I))])
    h = kzh_minimal_has_only_crazy_perturbation_1()
    bkpts = h.end_points()
    ucl = bkpts[17]
    ucr = bkpts[18]
    f = bkpts[37]
    special_intervals = RealSet.open(ucl, ucr) + RealSet.open(f - ucr, f - ucl)
    s = delta_pi_general(h, bkpts[39], 1 + ucl - bkpts[39], (-1, 0, 0))
    assert s == 19/23998
    logging.debug("s = {}".format(s))
    for F in generate_all_faces(h):
        n_F = number_of_projections_intersecting(F, special_intervals)
        if n_F:
            logging.debug("{} n_F = {}".format(F, n_F))
            deltas = sorted(set( delta_pi_of_face(h, vertex[0], vertex[1], F)
                                 for vertex in F.vertices ))
            if deltas == [0]:
                logging.debug("... additive face")
            else:
                logging.debug("... Delta pi values {}".format(deltas))
                assert deltas[0] >= n_F * s
                if deltas[0] == n_F * s:
                    logging.debug("... tight")
                    assert len(deltas) > 1  # strict for at least one vertex
