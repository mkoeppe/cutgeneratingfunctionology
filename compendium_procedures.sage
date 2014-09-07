# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

def transform_coho_interval(interval, shift, divisor):
    """Shift `interval` by `shift` and then divide it by `divisor`.

    EXAMPLES:
    sage: transform_coho_interval([1, 3], 1, 2)
    [1, 2]
    sage: transform_coho_interval([1, 3], 1, -1)
    [-4, -2]
    sage: transform_coho_interval(right_open_interval(-1,1), 1, 2)==right_open_interval(0, 1)
    True
    sage: transform_coho_interval(right_open_interval(-1,1), 1, -1)==left_open_interval(-2, 0)
    True
    """
    # This really wants to be in a compose method of FastPiecewise.
    x, y = (interval[0 ] + shift) / divisor, (interval[1 ] + shift) / divisor
    if divisor > 0:
        if len(interval) == 2:         # old-fashioned closed interval
            return [x, y]
        else:
            return closed_or_open_or_halfopen_interval(x, y, interval.left_closed, interval.right_closed)
    else:
        if len(interval) == 2:         # old-fashioned closed interval
            return [y, x]
        else:
            return closed_or_open_or_halfopen_interval(y, x, interval.right_closed, interval.left_closed)

def multiplicative_homomorphism(function, multiplier):
    """
    Construct the function x -> function(multiplier * x). 
    `multiplier` must be a nonzero integer.
    """
    if not multiplier in ZZ or multiplier == 0:
        raise ValueError, "Bad parameter multiplier, needs to be a nonzero integer."
    elif multiplier > 0:
        i_range = range(multiplier)
    else: # multiplier < 0
        i_range = range(multiplier, 0)
    # This really wants to be in a compose method of FastLinear.
    new_pairs = [ (transform_coho_interval(interval, i, multiplier), \
                   FastLinearFunction(function._slope * multiplier, function._intercept - function._slope * i))
                  for (interval, function) in function.list() \
                  for i in i_range ]
    return FastPiecewise(new_pairs)

def projected_sequential_merge(g, n=1):
    """
    construct the one-dimensional projected sequential merge inequality: h = g(with f = nr) @_n^1 gmic(f = r).

    The notation "@_n^1" is given in [39] p.305, def.12 & eq.25 :
        h(x) = (xi(x) + n * g(fractional((n + 1) * x - r * xi(x)))) / (n + 1), where xi = gmic(r)

    Parameters:
        g (real) a valid inequality;
        n (integer).

    Function is known to be extreme under the conditions: c.f. [39] p.309, Lemma 1
        (1) g is a facet-defining inequality, with f=nr ;
        (2) E(g) is unique up to scaling, c.f.[39] p.289, def.6 & def.7 ;
        (3) [g]_{nr} is nondecreasing, c.f.[39] p.290, def.8 (where m=1).

    Note:
        g = gj_forward_3_slope() does not satisfy condition(3), but

        g = multiplicative_homomorphism(gj_forward_3_slope(),-1) satisfies (3).

    Examples:
        [39]  p.311, fig.5 ::

        sage: logging.disable(logging.INFO)
        sage: g = multiplicative_homomorphism(gj_forward_3_slope(f=2/3, lambda_1=1/4, lambda_2=1/4), -1)
        sage: h = projected_sequential_merge(g, n=1)

    Reference:
        [39] SS Dey, JPP Richard, Relations between facets of low-and high-dimensional group problems,
             Mathematical programming 123 (2), 285-313.
    """
    f = find_f(g)
    r = f / n
    xi = gmic(r)
    ith = bisect_left(g.end_points(), f)
    l = len(g.intervals())
    multiplier = (1 + n - f) / (1 - r)
    new_pairs = [ (transform_coho_interval(interval, 0, n), \
                   FastLinearFunction(function._slope * n, function._intercept))
                  for (interval, function) in g.list()[0:ith]]
    i = r * multiplier - f
    new_pairs = new_pairs + \
                [ (transform_coho_interval(interval, i, multiplier), \
                   FastLinearFunction(function._slope * multiplier, function._intercept - function._slope * i))
                  for (interval, function) in g.list()[ith:l] ]
    new_pairs = new_pairs + \
                [ (transform_coho_interval(interval, j + 1/(1-r), multiplier), \
                   FastLinearFunction(function._slope * multiplier, function._intercept - function._slope * (j + 1/(1-r))))
                  for (interval, function) in g.list() \
                  for j in range(n) ]
    new_g = FastPiecewise(new_pairs)
    h = (1 / (n+1)) * xi + (1 / (n+1)) * new_g
    return h

def finite_group_order_from_function_f_oversampling_order(fn, f=None, oversampling=None, order=None):
    """
    Determine a finite group order to use, based on the given parameters.

    EXAMPLES::
    sage: logging.disable(logging.INFO)
    sage: finite_group_order_from_function_f_oversampling_order(gmic(f=4/5), order=17)
    17
    sage: finite_group_order_from_function_f_oversampling_order(gmic(f=4/5))
    5
    sage: finite_group_order_from_function_f_oversampling_order(gj_forward_3_slope())
    45
    sage: finite_group_order_from_function_f_oversampling_order(gmic(f=4/5), oversampling=3)
    15
    """
    if not order is None:
        return order
    if not order is None and not oversampling is None:
        raise ValueError, "Only one of the `order` and `oversampling` parameters may be provided."
    if f is None:
        f = find_f(fn, no_error_if_not_minimal_anyway=True)
    bkpt_f = fn.end_points()
    if not f is None:
        bkpt_f = [f] + bkpt_f
    is_rational_bkpt_f, bkpt_f = is_all_QQ(bkpt_f)
    if is_rational_bkpt_f:
        bkpt_f_denominator = [denominator(x) for x in bkpt_f]
        q = lcm(bkpt_f_denominator)
        if oversampling is None:
            logging.info("Rational breakpoints and f; using group generated by them, (1/%s)Z" % q)
            return q
        else:
            grid_nb = oversampling * q
            logging.info("Rational breakpoints and f; using group generated by them, refined by the provided oversampling factor %s, (1/%s)Z" % (oversampling, grid_nb))
            if oversampling >= 3 and fn.is_continuous_defined():
                # TODO: Detect and handle the case that fn already was a finite group problem
                logging.info("This is a continuous function; because oversampling factor is >= 3, the extremality test on the restricted function will be equivalent to the extremality test on the original function.")
            return grid_nb
    else:
        raise ValueError, "This is a function with irrational breakpoints or f, so no natural finite group order is available; need to provide an `order` parameter."

def restrict_to_finite_group(function, f=None, oversampling=None, order=None):
    """Restrict the given `function` to the cyclic group of given `order`.
    
    If `order` is not given, it defaults to the group generated by the
    breakpoints of `function` and `f` if these data are rational.
    However, if `oversampling` is given, it must be a positive
    integer; then the group generated by the breakpoints of `function`
    and `f` will be refined by that factor.

    If `f` is not provided, uses the one found by `find_f`.

    Assume in the following that `f` and all breakpoints of `function`
    lie in the cyclic group and that `function` is continuous.

    Then the restriction is valid if and only if `function` is valid.
    The restriction is minimal if and only if `function` is minimal.
    The restriction is extreme if `function` is extreme. 
    FIXME: Add reference.

    If, in addition `oversampling` >= 3, then the following holds:
    The restriction is extreme if and only if `function` is extreme.
    This is Theorem 1.5 in [IR2].
    
    Reference:
        [IR2] A. Basu, R. Hildebrand, and M. Köppe, Equivariant perturbation in Gomory and Johnson’s infinite group problem.
                I. The one-dimensional case, Mathematics of Operations Research (2014), doi:10. 1287/moor.2014.0660

    EXAMPLES::
    sage: logging.disable(logging.WARN)
    sage: g = gj_2_slope()
    sage: extremality_test(g)
    True
    sage: gf = restrict_to_finite_group(g)
    sage: minimality_test(gf)
    True
    sage: finite_dimensional_extremality_test(gf)
    True
    sage: h = drlm_not_extreme_1()
    sage: extremality_test(h)
    False
    sage: h7 = restrict_to_finite_group(h)
    sage: h7.end_points()
    [0, 1/7, 2/7, 3/7, 4/7, 5/7, 6/7, 1]
    sage: minimality_test(h7)
    True
    sage: finite_dimensional_extremality_test(h7)
    True
    sage: h21 = restrict_to_finite_group(h, oversampling=3)
    sage: minimality_test(h21)
    True
    sage: finite_dimensional_extremality_test(h21)
    False
    sage: h28 = restrict_to_finite_group(h, order=28)
    sage: minimality_test(h28)
    True
    sage: finite_dimensional_extremality_test(h28)
    False
    sage: h5 = restrict_to_finite_group(h, order=5)
    sage: minimality_test(h5)
    False

    Reference:
        [IR2] A. Basu, R. Hildebrand, and M. Köppe, Equivariant perturbation in Gomory and Johnson’s infinite group problem.
                I. The one-dimensional case, Mathematics of Operations Research (2014), doi:10. 1287/moor.2014.0660
    """
    order = finite_group_order_from_function_f_oversampling_order(function, f, oversampling, order)
    pieces = [ (singleton_interval(x/order), FastLinearFunction(0, function(x/order))) for x in range(order+1) ]
    return FastPiecewise(pieces)

def interpolate_to_infinite_group(function):
    """Interpolate the given `function` to make it a function in the
infinite group problem.  

`function` may be a function of a finite (cyclic) group problem, 
represented as a `FastPiecewise` with singleton intervals within [0,1] as its parts.

(`function` is actually allowed, however, to be more general; it can be any `FastPiecewise`.)

    See `restrict_to_finite_group` for a discussion of the relation of
    the finite and infinite group problem.

    EXAMPLES::
    
    The same as restrict_to_finite_group(drlm_not_extreme_1()):

    sage: logging.disable(logging.WARN)
    sage: h7 = discrete_function_from_points_and_values([0/7, 1/7, 2/7, 3/7, 4/7, 5/7, 6/7, 7/7], [0/10, 4/10, 8/10, 5/10, 2/10, 6/10, 10/10, 0/10])
    sage: finite_dimensional_extremality_test(h7)
    True
    sage: h = interpolate_to_infinite_group(h7)
    sage: extremality_test(h)
    False
    sage: h21 = restrict_to_finite_group(h, oversampling=3)
    sage: finite_dimensional_extremality_test(h21)
    False
    sage: h28 = restrict_to_finite_group(h, oversampling=4)
    sage: finite_dimensional_extremality_test(h28)
    False
    sage: h14 = restrict_to_finite_group(h, oversampling=2)
    sage: finite_dimensional_extremality_test(h14)
    True

    """
    # TODO: Allow `function` to be just a list of (x, y) pairs.
    last_x, last_y = None, None
    pieces = []
    # We are not using
    # `piecewise_function_from_breakpoints_and_values` to allow to
    # interpolate some partially interpolated functions.
    for (interval, fn) in function.list():
        x = interval[0]
        y = fn(x)
        if last_x is not None and last_x < x:
            slope = (y - last_y) / (x - last_x)
            intercept = last_y - slope * last_x
            pieces.append(([last_x, x], FastLinearFunction(slope, intercept)))
        last_x = interval[1]
        last_y = fn(last_x)
    return FastPiecewise(pieces)
