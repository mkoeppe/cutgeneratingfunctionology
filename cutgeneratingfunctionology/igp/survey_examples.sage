from six.moves import range
## Various examples of functions that appear in the survey.

def not_minimal_1(): # was not_minimal.sage
    r"""
    A non-minimal function.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = not_minimal_1()
        sage: minimality_test(h, False)
        False
    """
    return piecewise_function_from_breakpoints_and_values([0, 1/5, 2/5, 4/5, 1], [0, 1/5, 3/5, 1, 0])

def not_minimal_2(): # was not_minimal_2.sage
    r"""
    A non-minimal function.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = not_minimal_2()
        sage: minimality_test(h, False)
        False
    """
    return piecewise_function_from_breakpoints_and_values([0, 1/5, 2/5, 3/5, 4/5, 1], [0, 1/5, 1/2, 4/5, 1, 0])

def not_extreme_1(): # was symmetric_rational_function1.sage
    r"""
    A non-extreme, minimal function.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = not_extreme_1()
        sage: minimality_test(h, False)
        True
        sage: extremality_test(h, False)
        False
    """
    slopes = [10/3,0,10/3,0,10/3,-10/3,0,-10/3,0,-10/3]
    interval_lengths = [1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10]
    return piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes)

def drlm_not_extreme_1():
    r"""Example from S. S. Dey, J.-P. P. Richard, Y. Li, and L. A. Miller,
    On the extreme inequalities of infinite group problems,
    Mathematical Programming 121 (2009), no. 1, 145-170,
    https://doi:10.1007/s10107-008-0229-6.
    Figure 1.

    It is the interpolation of the extreme function C13 for the finite
    group problem of order 7 from Araoz, Evans, Gomory and Johnson.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = drlm_not_extreme_1()
        sage: minimality_test(h, False)
        True
        sage: extremality_test(h, False)
        False
    """
    bkpt = [0, 2/7, 4/7, 6/7, 1]
    values = [0, 4/5, 1/5, 1, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def drlm_not_extreme_2():
    r"""
    Example from S. S. Dey, J.-P. P. Richard, Y. Li, and L. A. Miller,
    On the extreme inequalities of infinite group problems,
    Mathematical Programming 121 (2009), no. 1, 145-170,
    https://doi:10.1007/s10107-008-0229-6.
    Figure 3.

    Note: this is not any of ``drlm_2_slope_limit`` functions,
    since here ``s_positive = 3``, whereas in
    ``drlm_2_slope_limit(f=1/2, nb_pieces_left=2, nb_pieces_right=2)``,
    the ``s_positive`` has to be 4.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: h = drlm_not_extreme_2()
        sage: minimality_test(h, False)
        True
        sage: extremality_test(h, False)
        False
    """
    f1 = FastLinearFunction(QQ(3), 0)
    f2 = FastLinearFunction(QQ(0), 1/2)
    f3 = FastLinearFunction(QQ(3), -1/2)
    f4 = FastLinearFunction(QQ(3), -4/3)
    f6 = FastLinearFunction(QQ(3), -13/6)
    f7 = FastLinearFunction(QQ(0), 0)
    return FastPiecewise([[right_open_interval(0, 1/4), f1], \
                          [singleton_interval(1/4), f2], \
                          [left_open_interval(1/4, 1/2), f3], \
                          [open_interval(1/2, 3/4), f4], \
                          [singleton_interval(3/4), f2], \
                          [open_interval(3/4, 1),f6], \
                          [singleton_interval(1),f7]])

def phi_s_in_drlm_not_extreme_2(s=10):
    r"""Example from S. S. Dey, J.-P. P. Richard, Y. Li, and L. A. Miller,
    On the extreme inequalities of infinite group problems,
    Mathematical Programming 121 (2009), no. 1, 145-170,
    https://doi:10.1007/s10107-008-0229-6.
    Figure 2.

    `s` is an integer; `s > 2`.
    
    ``phi_s_in_drlm_not_extreme_2(s)`` is an extreme function.
    The pointwise limit as `s` tends to `\\infty` is not extreme. See ``drlm_not_extreme_2()``

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = phi_s_in_drlm_not_extreme_2()
        sage: extremality_test(h, False)
        True
    """
    return  gj_2_slope_repeat(f=1/2, s_positive=3, s_negative=-s, m=2, n=3)
#    f1(x) = 3*x
#    f2(x) = -s * (x - 1/4) + 1/2
#    f3(x) = 3*x - 1/2
#    f4(x) = -s*x + 1 + s/2
#    f5(x) = 3*x - 4/3
#    f6(x) = -s * (x - 3/4) + 1/2
#    f7(x) = 3*x - 13/6
#    f8(x) = -s * (x - 1)
#    a = (s/4 + 1/2)/(s + 3)
#    b = (s/4 + 1)/(s + 3)
#    c = 1/2
#    d = (3*s + 14)/6/(s + 3)
#    e = (9*s + 22)/12/(s + 3)
#    f = (9*s + 32)/12/(s + 3)
#    g = (12*s + 26)/12/(s + 3)
#    return FastPiecewise([[right_open_interval(0, a), f1], \
#                          [closed_interval(a, b), f2], \
#                          [left_open_interval(b, c), f3], \
#                          [left_open_interval(c, d), f4], \
#                          [open_interval(d, e), f5], \
#                          [closed_interval(e, f), f6], \
#                          [open_interval(f, g),f7], \
#                          [closed_interval(g, 1),f8]])

def bhk_irrational_extreme_limit_to_rational_nonextreme(n=Infinity):
    r"""
    A sequence of ``bhk_irrational`` functions, each extreme, indexed by n = 1, 2, ...
    whose limit (n = Infinity) is a ``bhk_irrational`` function with rational parameters, 
    and hence not extreme. 

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = bhk_irrational_extreme_limit_to_rational_nonextreme(1)
        sage: extremality_test(h, False)
        True
        sage: h = bhk_irrational_extreme_limit_to_rational_nonextreme(2)
        sage: extremality_test(h, False)
        True
        sage: h = bhk_irrational_extreme_limit_to_rational_nonextreme(17)
        sage: extremality_test(h, False)
        True
        sage: h = bhk_irrational_extreme_limit_to_rational_nonextreme()
        sage: extremality_test(h, False)
        False
    """
    del1 = 1/60
    if n != Infinity:
        del1 -= sqrt(2) / (90*n)
    del2 = 1/60
    return bhk_irrational(delta=(del1, del2))

def drlm_gj_2_slope_extreme_limit_to_nonextreme(s=Infinity):
    r"""
    A sequence of ``phi_s_in_drlm_not_extreme_2`` functions, each extreme,
    indexed by s, (where s is a real number, s = abs(negative_slope) and s > 2)
    whose limit (s = Infinity) is a ``drlm_not_extreme_2`` function which is not extreme.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = drlm_gj_2_slope_extreme_limit_to_nonextreme(3)
        sage: extremality_test(h, False)
        True
        sage: h = drlm_gj_2_slope_extreme_limit_to_nonextreme(17)
        sage: extremality_test(h, False)
        True
        sage: h = drlm_gj_2_slope_extreme_limit_to_nonextreme()
        sage: extremality_test(h, False)
        False
    """
    if s != Infinity:
        return phi_s_in_drlm_not_extreme_2(s=s)
    return drlm_not_extreme_2()

def drlm_2_slope_limit_1_1(f=1/2, nb_pieces_left=1, nb_pieces_right=1):
    r"""
    An iconic choice of parameters in ``drlm_2_slope_limit``.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: h = drlm_2_slope_limit_1_1()
        sage: extremality_test(h, False)
        True
    """
    return drlm_2_slope_limit(f=f, nb_pieces_left=nb_pieces_left, nb_pieces_right=nb_pieces_right)

def chen_3_slope_not_extreme(f=1/2, lam=8):
    r"""
    A continuous 3-slope function, constructed by K. Chen in his Ph.D. thesis [KChen_thesis].
    The function has non-degenerate intervals with a zero derivative.

    Chen claims that the function is extreme under certain conditions. But his proof on page 37
    made a mistake in setting up the equations:
    `-\\pi(C) + \\pi(CC) + \\pi(D) - \\pi(DD) = 0`. 
    This relation should not exist in `E(\\pi)`.
    Returned function can be extreme or NOT. Need to use ``extremality_test()`` to check.

    See the discussion in [KZh2015b, section 3].

    Parameters:
        * f (real) `\in (0,1)`;
        * lam (real): the first slope has length 1/lam.

    Requirement:
        * 0 < f <= 1/2;
        * lam > 3*max(1/f, 1 /(1-f))

    Examples:
        [KChen_thesis]  p.33, fig.7 NOT extreme::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: h = chen_3_slope_not_extreme(f=1/2, lam=8)
            sage: extremality_test(h, False)
            False

        Extreme example (with parameters NOT satisfying the requirement)::

            sage: h = chen_3_slope_not_extreme(f=2/3, lam=20)
            sage: extremality_test(h, False)
            True

    References:

    - [KChen_thesis]:  K. Chen, Topics in group methods for integer programming, Ph.D. thesis, Georgia Institute of Technology, June 2011.

    - [KZh2015b] M. Koeppe and Y. Zhou, An electronic compendium of extreme functions for the Gomory-Johnson infinite group problem, Operations Research Letters, 2015, http://dx.doi.org/10.1016/j.orl.2015.06.004
    """
    if not (bool(0 < f < 1) and (lam > 3*max(1/f, 1 /(1-f)))):
        raise ValueError("Bad parameters. Unable to construct the function.")
    alpha = f / 2 - 3 / (2*lam)
    beta = 1/2 - f/2 - 3 / (2*lam)
    bkpts = [0, 1/lam, 1/lam + alpha, 2/lam + alpha, 2/lam + 2*alpha, f, \
                4/lam + 2*alpha, 4/lam + 2*alpha + beta, 5/lam + 2*alpha + beta, 5/lam + 2*alpha + 2*beta, 1]
    values = [0, 2/3, 2/3, 1/3, 1/3, 1, 2/3, 2/3, 1/3, 1/3, 0]
    return  piecewise_function_from_breakpoints_and_values(bkpts, values)

def dr_projected_sequential_merge_3_slope(f=2/3, lambda_1=1/2, lambda_2=1/2, n=1):
    r"""
    Construct the one-dimensional projected sequential merge inequality: `h = g \lozenge_n^1 \xi`, where
        * `g =` ``multiplicative_homomorphism(gj_forward_3_slope(f=f, lambda_1=lambda_1, lambda_2=lambda_2),-1);``
        * `\xi =` ``gmic(f/n)``.

    See ``projected_sequential_merge()``.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = dr_projected_sequential_merge_3_slope()
        sage: extremality_test(h, False)
        True

    Reference:  p.311, fig.5,[39] SS Dey, JPP Richard, Relations between facets of low-and high-dimensional group problems,
                Mathematical programming 123 (2), 285-313.
    """
    g = multiplicative_homomorphism(gj_forward_3_slope(f=f, lambda_1=lambda_1, lambda_2=lambda_2),-1)
    h = projected_sequential_merge(g, n=n)
    return h

def gomory_fractional(f=4/5):
    r"""
    The Gomory fractional cut.  
    Not minimal.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = gomory_fractional(f=4/5)
        sage: minimality_test(h, f=4/5)
        False
    """
    h = FastPiecewise([right_open_piece((0, 0), (1, 1)), 
                       singleton_piece(1, 0)])
    # now normalize...
    gf = h * (1/f)
    # Remember what is f (mostly for plotting purposes)
    gf._f = f
    return gf

def ll_strong_fractional_bad_figure_3():
    r"""
    Corresponds to Figure 3 in Letchford-Lodi (2002); divided by its
    value at f=2/3 to normalize.

    This function is not subadditive, contrary to what is claimed.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = ll_strong_fractional_bad_figure_3()
        sage: minimality_test(h, False)
        False
    """
    f = 2/3
    h = FastPiecewise([closed_piece((0, 0), (2/3, 2/3)),
                       open_piece((2/3, 0), (1, 1/3)),
                       singleton_piece(1, 0)])
    return h * (1/f)

def ll_strong_fractional_bad_figure_3_corrected():
    r"""
    Corresponds to what Figure 3 in Letchford-Lodi (2002) should have
    looked like; divided by its value at f=2/3 to normalize.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = ll_strong_fractional_bad_figure_3_corrected()
        sage: extremality_test(h, False)
        True
        sage: h == ll_strong_fractional()
        True
    """
    f = 2/3
    h = FastPiecewise([closed_piece((0, 0), (2/3, 2/3)),
                       open_piece((2/3, 2/3 - 1/2), (1, 1 - 1/2)),
                       singleton_piece(1, 0)])
    return h * (1/f)

def california_ip():
    r"""
    The California Integer Programming cut, to be plotted with rgbcolor=(192, 54, 44).
    Not minimal.
    """
    x = var('x')
    scale = 1/10
    shift = 1
    g1(x) = scale * (x^2 + 4*x + 4) + shift
    g2(x) = scale * x^2 + shift
    g3(x) = scale * (x^2 - 4*x + 4) + shift
    ggb = FastPiecewise([([-2, -1], g1),
                         ([-1, 1], g2),
                        ([1, 2], g3)])
    pieces = [[(0, 5/14), FastLinearFunction(14/5, QQ(0))],
              [(6/7, QQ(1)), FastLinearFunction(-QQ(7), QQ(7))]]
    pieces += [ transform_piece_to_interval(piece, [5/14, 6/7], [-2, 2]) 
                for piece in ggb.list() ]
    logging.warning("This function is not piecewise linear; code for handling this function is not implemented.")
    return FastPiecewise(pieces)

def kzh_2q_example_1():
    r"""
    A continuous 4-slope non-extreme function, whose restriction to
    1/2q is extreme, thereby showing that an oversampling factor of 3
    is optimal.

    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    Example::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_2q_example_1()
        sage: extremality_test(h)
        False
        sage: h_1q = restrict_to_finite_group(h, oversampling=1)
        sage: extremality_test(h_1q)
        True
        sage: h_2q = restrict_to_finite_group(h, oversampling=2)
        sage: extremality_test(h_2q)
        True
        sage: h_3q = restrict_to_finite_group(h, oversampling=3)
        sage: extremality_test(h_3q)
        False

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].

    """
    bkpt = [0, 3/37, 4/37, 5/37, 9/37, 10/37, 11/37, 12/37, 13/37, 14/37, 15/37, 16/37,
            20/37, 21/37, 22/37, 25/37, 26/37, 28/37, 29/37, 33/37, 34/37, 36/37, 1]
    values = [0, 99/122, 21/122, 27/61, 19/61, 71/122, 15/61, 63/122, 59/122, 46/61, 51/122,
                42/61, 34/61, 101/122, 23/122, 1, 22/61, 18/61, 69/122, 53/122, 43/61, 39/61, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def zhou_two_sided_discontinuous_cannot_assume_any_continuity():
    r"""
    This function is two-sided discontinuous at the origin.  The
    extremality test then cannot assume any left or right continuity
    of the perturbation function at the breakpoints, even if such
    one-sided continuity is present in the function.

    This example has been constructed by Yuan Zhou (2014).

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = zhou_two_sided_discontinuous_cannot_assume_any_continuity()
        sage: extremality_test(h)
        False
    """
    return FastPiecewise([singleton_piece(0, 0),
                          open_piece((0, 2/3), (1/8, 3/4)),
                          closed_piece((1/8, 1/4), (3/8, 3/4)),
                          open_piece((3/8, 1/4), (4/8, 1/3)),
                          singleton_piece(4/8, 1),
                          open_piece((4/8, 1/3), (8/8, 2/3)),
                          singleton_piece(1, 0)])

### The sporadic extreme functions have been moved to 'extreme_functions_sporadic.sage'

## Example functions that appear in the paper Equivariant perturbation V.

def equiv5_random_discont_1():
    r"""
    A randomly generated discontinuous function that appears in Equiv V.
    """
    return discontinuous_interpolation([0, 1/5, 2/5, 3/5, 4/5],
                                       [0, 1, 2/5, 1/2, 3/5], [0, 1, 2/5, 3/5, 1], [0, 1, 0, 2/5, 3/5])

## Example functions that appear in the paper Equivariant perturbation VII.
def equiv7_example_1():
    r"""
    One-sided discontinuous minimal valid function that appears in Equiv VII example 7.8.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h1 = equiv7_example_1()
        sage: generate_uncovered_components(h1)
        [[<Int(0, 1/4)>, <Int(1/4, 1/2)>]]
        sage: finite_dimensional_extremality_test(h1)
        False
        sage: len(h1._perturbations)
        1
    """
    return FastPiecewise([singleton_piece(0, 0),open_piece((0, 1/2), (1/2, 1/2)),closed_piece((1/2, 1),(1,0))])

def equiv7_example_2_crazy_perturbation():
    r"""
    An effective perturbation function for the two-sided discontinuous function
    `\pi_2 =` ``minimal_no_covered_interval()`` that appears in Equiv VII example 7.9.

    This perturbation is bounded, but is not Lipschitz continuous on the interval
    `(0, 1/2)`. In fact, it is a highly discontinuous "locally microperiodic"
    perturbation, which does not have a limit at any point in `(0, 1/2)`.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = minimal_no_covered_interval()
        sage: cp = equiv7_example_2_crazy_perturbation()
        sage: find_epsilon_for_crazy_perturbation(h, cp)
        1/6
    """
    t1, t2 = nice_field_values([1, sqrt(2)])
    generators = [t1, t2]
    pwl = piecewise_function_from_breakpoints_and_slopes([0,1],[0])
    crazy_piece_1 = CrazyPiece((0, 1/4), generators, [(0, 1), (1/4, -1)])
    crazy_piece_2 = CrazyPiece((1/4, 1/2), generators, [(1/4, 1), (1/2, -1)])
    cp = PiecewiseCrazyFunction(pwl, [crazy_piece_1, crazy_piece_2])
    return cp

def equiv7_minimal_2_covered_2_uncovered():
    r"""
    A continuous minimal valid function that appears in Equiv VII.

    Five 0-slopes, all uncovered, translation and reflection on all five.
    (need one more equation)
    """
    # from functions_with_horizontal_slopes.sage example_7
    bkpts0 = [0,3,4,7,8,11,12,13,14,15,16,17,18,19,20,21,22,25,26,29,30,33,49]
    bkpts = [bkpt/49 for bkpt in bkpts0]
    values0 = [0,3,2,5,4,7,6,6,7,7,8,8,9,9,10,10,9,12,11,14,13,16,0]
    values = [value/16 for value in values0]
    h = piecewise_function_from_breakpoints_and_values(bkpts, values)
    return h

def equiv7_example_3():
    bkpt = [ x/10 for x in range(11) ]
    value = [0, 3/9, 3/9, 3/9, 4/9, 5/9, 6/9, 6/9, 6/9, 9/9, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, value)

def equiv7_example_xyz_1():
    bkpt = [0,1/12,3/12,4/12,7/12,8/12,10/12,11/12,12/12]
    value =[0,9/15,4/15,6/15,9/15,11/15,6/15,15/15,0]
    return piecewise_function_from_breakpoints_and_values(bkpt, value)

def equiv7_example_xyz_2():
    bkpt = [0,1/24,3/12,4/12,7/12,8/12,21/24,11/12,12/12]
    value =[0,29/40,4/15,6/15,9/15,11/15,11/40,15/15,0]
    return piecewise_function_from_breakpoints_and_values(bkpt, value)

def equiv7_example_xyz_3():
    bkpt = [0,1/36,3/12,4/12,7/12,8/12,32/36,11/12,12/12]
    value =[0,34/45,4/15,6/15,9/15,11/15,11/45,15/15,0]
    return piecewise_function_from_breakpoints_and_values(bkpt, value)

# Examples with posterior moves
def equiv7_example_post_1():
    bkpt = [0, 1/18, 1/9, 1/6, 5/18, 1/3, 4/9, 1/2, 11/18, 2/3, 13/18, 7/9, 5/6, 17/18, 1]
    value = [0, 3/4, 1/2, 3/4, 1/4, 1/2, 1/2, 3/4, 1/4, 1/2, 1/4, 1, 1/4, 3/4, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, value)

def equiv7_example_post_2():
    bkpt = [0, 1/9, 1/6, 2/9, 5/18, 1/3, 7/18, 4/9, 1/2, 5/9, 11/18, 13/18, 5/6, 8/9, 1]
    value = [0, 2/3, 2/3, 1/3, 1/3, 2/3, 1/3, 2/3, 2/3, 1/3, 1/3, 1, 1/3, 2/3, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, value)

def equiv7_example_post_3():
    bkpt = [0, 1/18, 1/6, 2/9, 5/18, 1/3, 7/18, 4/9, 1/2, 5/9, 2/3, 13/18, 5/6, 8/9, 1]
    value =[0, 2/3, 2/3, 1/3, 2/3, 2/3, 1/3, 1/3, 2/3, 1/3, 1/3, 1, 1/3, 2/3, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, value)

def equiv7_example_post_4():
    bkpt = [0, 1/9, 2/9, 1/3, 4/9, 5/9, 11/18, 2/3, 8/9, 17/18, 1]
    value =[0, 1/2, 1/4, 3/4, 1/2, 1, 1/2, 3/4, 1/4, 1/2, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, value)
