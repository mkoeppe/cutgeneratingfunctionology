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

def kzh_minimal_has_only_crazy_perturbation_1(parametric=False, field=None, **parametric_kwds):
    r"""
    This function is a two-sided discontinuous piecewise linear minimal valid function
    introduced in :cite:`koeppe-zhou:crazy-perturbation`
    which is not extreme, but which is not a convex combination of other piecewise linear minimal valid
    functions.  It has two special intervals `(l, u)` and `(f-u, f-l)`, on which every 
    nonzero perturbation is microperiodic (invariant under the action of a dense additive group).

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = kzh_minimal_has_only_crazy_perturbation_1()
        sage: minimality_test(h)
        True

    On the one hand, there do not exist any nonzero piecewise continuous perturbations::

        sage: extremality_test(h, crazy_perturbations=False)
        True

    To attempt to reproduce the automatic proof of this fact that is given in
    :cite:`koeppe-zhou:crazy-perturbation`, Appendix C, we can do the following::

        sage: logging.getLogger().setLevel(logging.DEBUG)  # not tested - for interactive use
        sage: import cutgeneratingfunctionology.igp as igp 
        sage: igp.strategical_covered_components = True
        sage: igp.show_values_of_fastpiecewise = False
        sage: igp.show_RNFElement_by_embedding = False
        sage: h = kzh_minimal_has_only_crazy_perturbation_1()         # long time
        sage: extremality_test(h, crazy_perturbations=False)          # long time
        True
        sage: igp.strategical_covered_components = False
        sage: igp.show_values_of_fastpiecewise = True
        sage: igp.show_RNFElement_by_embedding = True

    Unfortunately, the details of the code have changed and it now makes different choices for
    the full-rank subsystem, so the subsystem shown in the publication cannot be reproduced
    with the current version of the code.

    Let's look a bit closer.  We take out the component that is only densely covered,
    consisting of the two special intervals::

        sage: proper_covered_components = generate_covered_components(h)[:2]

    ``symbolic`` is the template for the perturbation function `\bar\pi`; it is
    a general element of the vector space generated by the discontinuous piecewise
    functions with same slopes on the intervals of each covered component, defined
    outside of the special intervals::

        sage: symbolic = generate_symbolic(h, proper_covered_components, basis_functions=('midpoints', 'slopes'))
        sage: plot_symbolic(symbolic, ymin=-1.1, ymax=1.1).show(figsize=(4,40))      # not tested

    Here we used the midpoints-and-slopes basis, which is "local."  It makes sure that ``symbolic``
    is undefined on the special intervals::

        sage: symbolic((h.ucl + h.ucr) / 2)
        Traceback (most recent call last):
        ...
        ValueError: ...

    We generate equations from additive and limit-additive vertices, avoiding
    (by ``undefined_ok=True``) to write down any equation that touches the special intervals.
    (The solution set is a superset of the projection to the quotient by the subspace of 
    arbitrary real functions supported on the special intervals.)

        sage: M, vs = generate_additivity_equations(h, symbolic, reduce_system=True, return_vertices=True, undefined_ok=True)
        sage: perturbations = M.right_kernel().basis()
        sage: len(perturbations)
        0

    The solution set is trivial.  Thus, any perturbation is zero outside of the special intervals.

    On the other hand, there exist crazy perturbations, such as the one we construct below::

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

    In :cite:`koeppe-zhou:discontinuous-facets`, it shown that
    `\pi = ``kzh_minimal_has_only_crazy_perturbation_1()`` ` is a weak facet.
    In the proof, we take an arbitrary minimal valid function `\pi'` with
    `E(\pi) \subseteq E(\pi')`.  Because we have no control over the limit-additivities
    of the function `\pi'`, all of our arguments have to use addivities-sans-limit of `\pi`. 
    First we show that `\pi'` has to be affine linear on every non-special intervals of `\pi`.

        sage: h = kzh_minimal_has_only_crazy_perturbation_1()
        sage: additive_faces_sans_limits = list(generate_additive_faces_sans_limits(h))
        sage: covered_components = generate_covered_components_strategically(h, additive_faces=additive_faces_sans_limits)
        sage: uncovered_intervals = uncovered_intervals_from_covered_components(covered_components)
        sage: uncovered_intervals == [open_interval(h.ucl, h.ucr), open_interval(h._f - h.ucr, h._f - h.ucl)]
        True

    The above is also done by ``facet_test``::

        sage: h = kzh_minimal_has_only_crazy_perturbation_1()
        sage: facet_test(h)                                     # long time - 45s
        Traceback (most recent call last):
        ...
        NotImplementedError: facet_test does not know how to continue

    ``facet_test`` also sets up a finite system for a general partial function outside 
    the special intervals, again only using additivities-sans-limit.

        sage: len(list(b for t, b in h._facet_symbolic.basis if t == 'function value at' and b in h.end_points()))      # long time
        19
        sage: len(list(b for t, b in h._facet_symbolic.basis if t == 'function value at' and b not in h.end_points()))  # long time
        18

    Again this system has a full rank, and no nontrivial solution exists.

        sage: len(h._facet_solution_basis)      # long time
        0

    NOTE:

    If ``parametric=True``, the breakpoints and slopes are given symbolic names.
    This is useful for printing.

        sage: h = kzh_minimal_has_only_crazy_perturbation_1(parametric=True)
        sage: h.which_pair(1/100)
        (<Int(0, (x1)~)>, <FastLinearFunction ((c3)~)*x + (0.1553846153846154?)>)
        sage: h.which_pair(250/800)
        (<Int(l~, u~)>, <FastLinearFunction ((c2)~)*x + (0.3482269355779649?)>)
        sage: h.which_pair(1/2)
        (<Int((f - u)~, (f - l)~)>, <FastLinearFunction ((c2)~)*x + (0.6514397033086091?)>)

    Note:

        This example is obtained by the following code::

            sage: import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.WARN)
            sage: igp.generate_symbolic_two_sided_discontinuous_basis_functions = ('slopes', 'jumps')   # Restore classic behavior
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
    [sqrt2, o] = nice_field_values([sqrt(2), 1], field=field)
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
    h = FastPiecewise(pieces)
    def set_special_attributes(h):
        bkpts = h.end_points()
        h.ucl = bkpts[17]
        h.ucr = bkpts[18]
        h._f = bkpts[37]
        h.a0 = bkpts[6]
        h.a1 = bkpts[10]
        h.a2 = bkpts[13]
        h.t1 = h.a1 - h.a0
        h.t2 = h.a2 - h.a0
        from cutgeneratingfunctionology.spam.real_set import RealSet   # param-safe version of RealSet
        h.special_intervals = RealSet.open(h.ucl, h.ucr) + RealSet.open(h._f - h.ucr, h._f - h.ucl)
        h.s = delta_pi_general(h, bkpts[39], 1 + h.ucl - bkpts[39], (-1, 0, 0))
        assert h.s == 19/23998
    set_special_attributes(h)
    if parametric:
        bkpt_names_dict = { bkpt: 'x{}'.format(i) for i, bkpt in enumerate(h.end_points())
                            if bkpt not in (0, 1) }
        bkpt_names_dict[h.a0] = 'a0'
        bkpt_names_dict[h.a1] = 'a1'
        bkpt_names_dict[h.a2] = 'a2'
        bkpt_names_dict[h.ucl] = 'l'
        bkpt_names_dict[h.ucr] = 'u'
        bkpt_names_dict[h._f] = 'f'
        slope_names_dict = {35/13: 'c1', 5/11999: 'c2', -5: 'c3'}
        param_dict, K = param_dict_for_piecewise(h,
                                                 bkpt_names_dict=bkpt_names_dict,
                                                 slope_names_dict=slope_names_dict,
                                                 **parametric_kwds)
        param_dict[h._f - h.ucl] = param_dict[h._f] - param_dict[h.ucl]
        param_dict[h._f - h.ucr] = param_dict[h._f] - param_dict[h.ucr]
        for ai in (h.a0, h.a1, h.a2):
            param_dict[h._f - ai] = param_dict[h._f] - param_dict[ai]
            param_dict[h._f - ai] = param_dict[h._f] - param_dict[ai]
            param_dict[h._f - ai] = param_dict[h._f] - param_dict[ai]
        h, K = param_piecewise(h, param_dict=param_dict, field=K)
        set_special_attributes(h)
    return h

def kzh_minimal_has_only_crazy_perturbation_1_perturbation():
    r"""
    A crazy perturbation for ``kzh_minimal_has_only_crazy_perturbation_1``.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = kzh_minimal_has_only_crazy_perturbation_1()
        sage: cp = kzh_minimal_has_only_crazy_perturbation_1_perturbation()
        sage: find_epsilon_for_crazy_perturbation(h, cp)                       # long time - 70s
        0.0003958663221935161?

    """
    h = kzh_minimal_has_only_crazy_perturbation_1()
    generators = [h.t1, h.t2]
    pwl = piecewise_function_from_breakpoints_and_slopes([0, 1], [0])
    crazy_piece_1 = CrazyPiece((h.ucl, h.ucr), generators, [(h.ucl, 1), (h.ucr, -1)])
    crazy_piece_2 = CrazyPiece((h._f - h.ucr, h._f - h.ucl), generators, [(h._f - h.ucr, 1), (h._f - h.ucl, -1)])
    return PiecewiseCrazyFunction(pwl, [crazy_piece_1, crazy_piece_2])

def param_dict_for_piecewise(h, bkpt_names_dict=None, slope_names_dict={}, extra_names=[], extra_values=[], name_sort_key=None, base_ring=None):
    bkpts = h.end_points()
    if base_ring is None:
        base_ring = bkpts[0].parent()
    if bkpt_names_dict is None:
        bkpt_names_dict = { bkpt: 'x{}'.format(i) for i, bkpt in enumerate(bkpts) if bkpt not in (0, 1) }
    def default_name_sort_key(name):
        if name == 'f':
            return ''   # comes first so that we see  f - a0 instead of -a0 + f
        else:
            return name
    if name_sort_key is None:
        name_sort_key = default_name_sort_key
    value_name_pairs = sorted(((value, name)
                               for value, name in itertools.chain(bkpt_names_dict.items(),
                                                                  slope_names_dict.items())),
                              key=lambda value_name: name_sort_key(value_name[1]))
    names  = [ name  for value, name in value_name_pairs ] + extra_names
    values = [ value for value, name in value_name_pairs ] + extra_values
    K = ParametricRealField(names=names, values=values, base_ring=base_ring)
    K._record = False   # FIXME: This should be part of public interface of ParametricRealField
    values_gens = list(zip(values, K.gens()))
    if extra_names:   # they don't take part in renaming!
        values_gens = values_gens[:-len(extra_names)]
    return dict(values_gens), K

def param_piecewise(h, param_dict=None, bkpt_names_dict=None, slope_names_dict={}, field=None, **param_dict_kwargs):
    r"""
    Replace the piecewise function h by one that renames all breakpoints, and optionally slopes,
    symbolically.  Useful for latexing or plotting functions that do not have convenient
    parametric descriptions.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)                   # disable output for automatic tests
        sage: h = gmic()
        sage: hp, K = param_piecewise(h, slope_names_dict={h.functions()[0]._slope: 'c1'})
        sage: hp.which_pair(1/10)
        ((0, (x1)~), <FastLinearFunction ((c1)~)*x + (0)>)
    """
    if param_dict is None:
        param_dict, field = param_dict_for_piecewise(h, bkpt_names_dict=bkpt_names_dict,
                                                     slope_names_dict=slope_names_dict,
                                                     **param_dict_kwargs)
    def param(x):
        return param_dict.get(x, field(x))
    def param_interval(interval):
        if len(interval) <= 2:
            # old fashioned interval
            return [ param(x) for x in interval ]
        else:
            # coho interval
            return closed_or_open_or_halfopen_interval(param(interval.a), param(interval.b),
                                                       interval.left_closed, interval.right_closed)
    def param_function(function):
        return FastLinearFunction(param(function._slope), field(function._intercept))
    param_pieces = [ (param_interval(interval), param_function(function)) for interval, function in h.list() ]
    return FastPiecewise(param_pieces), field

def number_of_projections_intersecting(F, real_set):
    """
    This is the number `n_F` from facets-paper.
    """
    return len([I for I in F.minimal_triple
                if not real_set.is_disjoint_from(realset_from_interval(interval_mod_1(I)))])

def generate_all_faces(fn):
    zero_fn = FastPiecewise([(I, FastLinearFunction(0, 0)) for I, f in fn.list()], merge=False)
    return generate_additive_faces_general(zero_fn)

def generate_intervals_and_two_sided_discontinuous_breakpoints(function, halfopen_cover_only=False):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)                   # disable output for automatic tests

        sage: h = kzh_minimal_has_only_crazy_perturbation_1(parametric=True)
        sage: list(generate_intervals_and_two_sided_discontinuous_breakpoints(h, halfopen_cover_only=False))
        [[0], [0, (x1)~], [(x1)~], [(x1)~, (x2)~], [(x2)~], [(x2)~, (x3)~], [(x3)~],
        ..., [(x39)~], [(x39)~, 1], [1]]
        sage: list(generate_intervals_and_two_sided_discontinuous_breakpoints(h, halfopen_cover_only=True))
        [[0], [0, (x1)~], [(x1)~], [(x1)~, (x2)~], [(x2)~, (x3)~], [(x3)~],
        ..., [(x39)~], [(x39)~, 1], [1]]

    """
    bkpt = function.end_points()
    n = len(bkpt)
    for i in range(n):
        l = function.limits(bkpt[i])
        left_continuous = l[-1] is None or l[-1] == l[0]
        right_continuous = l[1] is None or l[1] == l[0]
        if not halfopen_cover_only or (not left_continuous and not right_continuous):
            yield [bkpt[i]]
        if i+1 < n:
            yield [bkpt[i], bkpt[i+1]]

def generate_faces_with_projections_intersecting(function, real_set, break_symmetry=False, halfopen_cover_only=False):
    for triple in generate_triples_with_projections_intersecting(function, real_set, break_symmetry=break_symmetry, halfopen_cover_only=halfopen_cover_only):
        yield Face(triple)

def generate_triples_with_projections_intersecting(function, real_set, break_symmetry=False, halfopen_cover_only=False):
    r"""
    It breaks symmetry but produces duplicates.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)                   # disable output for automatic tests
        sage: h = kzh_minimal_has_only_crazy_perturbation_1(parametric=True)
        sage: faces = set(generate_faces_with_projections_intersecting(h, h.special_intervals, break_symmetry=False))    # long time - 90s
        sage: sum(F.plot(fill_color='lightblue', rgbcolor='blue') for F in faces).show(figsize=15, xmax=1) # not tested
    """
    # Adapted from generate_additive_faces_general

    def plus1(I):
        return [ x + 1 for x in I ]

    I_list = list(generate_intervals_and_two_sided_discontinuous_breakpoints(function, halfopen_cover_only=halfopen_cover_only))
    K_list = I_list + [ plus1(I) for I in I_list ]

    I_disjoint_list = [ I for I in I_list if real_set.is_disjoint_from(realset_from_interval(I)) ]
    I_nondisjoint_list = [ I for I in I_list if not real_set.is_disjoint_from(realset_from_interval(I)) ]
    K_nondisjoint_list = I_nondisjoint_list + [ plus1(I) for I in I_nondisjoint_list ]

    def check_face(I, J, K):
        IplusJ = interval_sum(I, J)
        if len(interval_intersection(IplusJ, K)) > 0:
            F = Face((I, J, K))
            n_F = number_of_projections_intersecting(F, real_set)
            if n_F:
                return F

    for I in I_nondisjoint_list:
        for J in I_list:
            for K in K_list:
                if check_face(I, J, K):
                    yield (I, J, K)
                    if not break_symmetry:
                        yield (J, I, K)

    for K in K_nondisjoint_list:
        for i, I in enumerate(I_disjoint_list):
            for J in I_disjoint_list[i:]:
                if check_face(I, J, K):
                    yield (I, J, K)
                    if not break_symmetry:
                        yield (J, I, K)

def kzh_minimal_has_only_crazy_perturbation_1_check_subadditivity_slacks(parametric=False):
    r"""
    Check a claim in the proof of the theorem that ``kzh_minimal_has_only_crazy_perturbation_1``
    is a weak facet (but not extreme, nor a facet).

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)                   # disable output for automatic tests
        sage: logging.getLogger().setLevel(logging.DEBUG)     # not tested - for interactive use
        sage: kzh_minimal_has_only_crazy_perturbation_1_check_subadditivity_slacks()

    """
    h = kzh_minimal_has_only_crazy_perturbation_1(parametric=parametric)
    logging.debug("s = {}".format(h.s))
    for F in generate_faces_with_projections_intersecting(h, h.special_intervals,
                                                          halfopen_cover_only=True):
        n_F = number_of_projections_intersecting(F, h.special_intervals)
        assert n_F > 0
        if n_F:
            logging.debug("{} n_F = {}".format(F, n_F))
            deltas = sorted(set( delta_pi_of_face(h, vertex[0], vertex[1], F)
                                 for vertex in F.vertices ))
            if deltas == [0]:
                logging.debug("... additive face")
            else:
                logging.debug("... Delta pi values {}".format(deltas))
                assert deltas[0] >= n_F * h.s
                if deltas[0] == n_F * h.s:
                    logging.debug("... tight")
                    assert len(deltas) > 1  # strict for at least one vertex
