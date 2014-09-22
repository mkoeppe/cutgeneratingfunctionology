# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

## Various examples of functions that appear in the survey.

def not_minimal_1(): # was not_minimal.sage
    """
    A non-minimal function.

    EXAMPLES::
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = not_minimal_1()
        sage: minimality_test(h, False)
        False
    """
    return piecewise_function_from_breakpoints_and_values([0, 1/5, 2/5, 4/5, 1], [0, 1/5, 3/5, 1, 0])

def not_minimal_2(): # was not_minimal_2.sage
    """
    A non-minimal function.

    EXAMPLES::
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = not_minimal_2()
        sage: minimality_test(h, False)
        False
    """
    return piecewise_function_from_breakpoints_and_values([0, 1/5, 2/5, 3/5, 4/5, 1], [0, 1/5, 1/2, 4/5, 1, 0])

def not_extreme_1(): # was symmetric_rational_function1.sage
    """
    A non-extreme, minimal function.

    EXAMPLES::
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
    """Example from S. S. Dey, J.-P. P. Richard, Y. Li, and L. A. Miller,
    On the extreme inequalities of infinite group problems,
    Mathematical Programming 121 (2009), no. 1, 145–170,
    doi:10.1007/s10107-008-0229-6.
    Figure 1.

    It is the interpolation of the extreme function C13 for the finite
    group problem of order 7 from Araoz, Evans, Gomory and Johnson.

    EXAMPLES::
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = drlm_not_extreme_1()
        sage: minimality_test(h, False)
        True
        sage: extremality_test(h, False)
        False
    """
    return piecewise_function_from_robert_txt_file("dey-richard-not-extreme.txt")

def drlm_not_extreme_2():
    """Example from S. S. Dey, J.-P. P. Richard, Y. Li, and L. A. Miller,
    On the extreme inequalities of infinite group problems,
    Mathematical Programming 121 (2009), no. 1, 145–170,
    doi:10.1007/s10107-008-0229-6.
    Figure 3.
    Note: this is not any of drlm_2_slope_limit() funcitons,
          since here s_positive = 3, whereas in
          drlm_2_slope_limit(f=1/2, nb_pieces_left=2, nb_pieces_right=2),
          the s_positive has to be 4.

    EXAMPLES::
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
    """Example from S. S. Dey, J.-P. P. Richard, Y. Li, and L. A. Miller,
    On the extreme inequalities of infinite group problems,
    Mathematical Programming 121 (2009), no. 1, 145–170,
    doi:10.1007/s10107-008-0229-6.
    Figure 2.
    s is an integer;
    phi_s_in_drlm_not_extreme_2(s) is an extreme function.
    The pointwise limit as s tends to \infty is not extreme. see drlm_not_extreme_2()

    Note: s > 2

    EXAMPLES::
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
    """
    A sequence of `bhk_irrational` functions, each extreme, indexed by n = 1, 2, ...
    whose limit (n = Infinity) is a `bhk_irrational` function with rational parameters, 
    and hence not extreme. 

    EXAMPLES::
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
    """
    A sequence of `phi_s_in_drlm_not_extreme_2` functions, each extreme,
    indexed by s, (where s is a real number, s = abs(negative_slope) and s > 2)
    whose limit (s = Infinity) is a `drlm_not_extreme_2` function which is not extreme.

    EXAMPLES::
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
    """
    An iconic choice of parameters in drlm_2_slope_limit.

    EXAMPLES::
        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: h = drlm_2_slope_limit_1_1()
        sage: extremality_test(h, False)
        True
    """
    return drlm_2_slope_limit(f=f, nb_pieces_left=nb_pieces_left, nb_pieces_right=nb_pieces_right)

def hildebrand_5_slope_22_1():
    """
    One of Hildebrand's world-record 5-slope functions.

    EXAMPLES::
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = hildebrand_5_slope_22_1()
        sage: extremality_test(h, False)
        True
    """
    return piecewise_function_from_robert_txt_file("example5Slope22data.txt")

def hildebrand_5_slope_24_1():
    """
    One of Hildebrand's world-record 5-slope functions.

    EXAMPLES::
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = hildebrand_5_slope_24_1()
        sage: extremality_test(h, False)
        True
    """
    return piecewise_function_from_robert_txt_file("example5Slope24data.txt")

def hildebrand_5_slope_28_1():
    """
    One of Hildebrand's world-record 5-slope functions.

    EXAMPLES::
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = hildebrand_5_slope_28_1()
        sage: extremality_test(h, False)
        True
    """
    return piecewise_function_from_robert_txt_file("example5Slope28data.txt")

def chen_3_slope_not_extreme(f=1/2, lam=8):
    """
    A continuous 3-slope function that has non-degenerate intervals with a zero derivative.
    The paper claims that the function is extreme under conditions. BUT his proof on page 37
    made a mistake in setting up the equations:
    -\pi(C ) + \pi(CC) + \pi(D) - \pi(DD) = 0. This relation should not exist in E(\pi).
    Returned function can be extreme or NOT. Need to use extemality_test() to check.

    Parameters:
        f (real) \in (0,1);
        lam (real): the first slope has length 1/lam.

    Requirement:
        0 < f <= 1/2;
        lam > 3*max(1/f, 1 /(1-f))

    Examples:
        [KChen_thesis]  p.33, fig.7 NOT extreme::
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: h = chen_3_slope_not_extreme(f=1/2, lam=8)
            sage: extremality_test(h, False)
            False

        extreme example (with parameters not satisfying the requirement)::
            sage: h = chen_3_slope_not_extreme(f=2/3, lam=20)
            sage: extremality_test(h, False)
            True

    Reference:
        [KChen_thesis]:  K. Chen, Topics in group methods for integer programming,
                            Ph.D. thesis, Georgia Institute of Technology, June 2011.
    """
    if not (bool(0 < f < 1) and (lam > 3*max(1/f, 1 /(1-f)))):
        raise ValueError, "Bad parameters. Unable to construct the function."
    alpha = f / 2 - 3 / (2*lam)
    beta = 1/2 - f/2 - 3 / (2*lam)
    bkpts = [0, 1/lam, 1/lam + alpha, 2/lam + alpha, 2/lam + 2*alpha, f, \
                4/lam + 2*alpha, 4/lam + 2*alpha + beta, 5/lam + 2*alpha + beta, 5/lam + 2*alpha + 2*beta, 1]
    values = [0, 2/3, 2/3, 1/3, 1/3, 1, 2/3, 2/3, 1/3, 1/3, 0]
    return  piecewise_function_from_breakpoints_and_values(bkpts, values)

def dr_projected_sequential_merge_3_slope(f=2/3, lambda_1=1/4, lambda_2=1/4, n=1):
    """
    construct the one-dimensional projected sequential merge inequality: h = g @_n^1 xi, where
    g = multiplicative_homomorphism(gj_forward_3_slope(f=f, lambda_1=lambda_1, lambda_2=lambda_2),-1);
    xi = gmic(f/n).
    see projected_sequential_merge()

    EXAMPLES::
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

def singleton_piece(x, y):
    return (singleton_interval(x), FastLinearFunction(0, y))

def linear_function_through_points(p, q):
    slope = (q[1] - p[1]) / (q[0] - p[0])
    intercept = p[1] - slope * p[0]
    return FastLinearFunction(slope, intercept) 

def open_piece(p, q):
    return (open_interval(p[0], q[0]), linear_function_through_points(p, q))

def closed_piece(p, q):
    return (closed_interval(p[0], q[0]), linear_function_through_points(p, q))

def left_open_piece(p, q):
    return (left_open_interval(p[0], q[0]), linear_function_through_points(p, q))

def right_open_piece(p, q):
    return (right_open_interval(p[0], q[0]), linear_function_through_points(p, q))

def hildebrand_2_sided_discont_1_slope_1():
    """
    The first known example of function that is discontinuous on both
    sides of the origin but is also extreme.

    Constructed by Robert Hildebrand (previously unpublished). 

    EXAMPLES::
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = hildebrand_2_sided_discont_1_slope_1()
        sage: extremality_test(h, False)
        True
    """
    return FastPiecewise([singleton_piece(0, 0),
                          open_piece((0, 4/8), (1/8, 6/8)),
                          closed_piece((1/8, 2/8), (3/8, 6/8)),
                          open_piece((3/8, 2/8), (4/8, 4/8)),
                          singleton_piece(4/8, 1),
                          open_piece((4/8, 2/8), (5/8, 4/8)),
                          closed_piece((5/8, 2/8), (7/8, 6/8)),
                          open_piece((7/8, 4/8), (8/8, 6/8)),
                          singleton_piece(1, 0)])


def hildebrand_2_sided_discont_2_slope_1():
    """
    The second known example of function that is discontinuous on both
    sides of the origin but is also extreme.  This one has 2 slopes.

    Constructed by Robert Hildebrand (previously unpublished). 

    EXAMPLES::
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = hildebrand_2_sided_discont_2_slope_1()
        sage: extremality_test(h, False)
        True
    """
    return FastPiecewise([singleton_piece(0, 0),
                          open_piece((0, 4/8), (1/8, 2/8)),
                          singleton_piece(1/8, 6/8),
                          open_piece((1/8, 2/8), (3/8, 6/8)),
                          singleton_piece(3/8, 2/8),
                          open_piece((3/8, 6/8), (4/8, 4/8)),
                          singleton_piece(4/8, 1),
                          open_piece((4/8, 4/8), (5/8, 6/8)),
                          singleton_piece(5/8, 2/8),
                          open_piece((5/8, 6/8), (7/8, 2/8)),
                          singleton_piece(7/8, 6/8),
                          open_piece((7/8, 2/8), (1, 4/8)),
                          singleton_piece(1, 0)])

def hildebrand_discont_3_slope_1():
    """
    This is a very new discontinuous 3-slope function that is extreme.

    Constructed by Robert Hildebrand (previously unpublished).

    EXAMPLES::
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = hildebrand_discont_3_slope_1()
        sage: extremality_test(h, False)
        True
    """
    return FastPiecewise([right_open_piece((0, 0), (1/8, 6/8)),
                          closed_piece((1/8, 2/8), (3/8, 6/8)),
                          left_open_piece((3/8, 2/8), (4/8, 1)),
                          left_open_piece((4/8, 4/8), (5/8, 6/8)),
                          closed_piece((5/8, 6/8), (7/8, 2/8)), 
                          right_open_piece((7/8, 2/8), (1, 4/8)), 
                          singleton_piece(1, 0)])

def gomory_fractional(f=4/5):
    """
    The Gomory fractional cut.  
    Not minimal.

    EXAMPLES:
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = gomory_fractional(f=4/5)
        sage: minimality_test(h, f=4/5)
        False
    """
    h = FastPiecewise([right_open_piece((0, 0), (1, 1)), 
                       singleton_piece(1, 0)])
    # now normalize...
    return h * (1/f)

def ll_strong_fractional_bad_figure_3():
    """
    Corresponds to Figure 3 in Letchford-Lodi (2002); divided by its
    value at f=2/3 to normalize.

    This function is not subadditive, contrary to what is claimed.

    EXAMPLES::
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
    """
    Corresponds to what Figure 3 in Letchford-Lodi (2002) should have
    looked like; divided by its value at f=2/3 to normalize.

    EXAMPLES::
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
    """
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
    logging.warn("This function is not piecewise linear; code for handling this function is not implemented.")
    return FastPiecewise(pieces)
