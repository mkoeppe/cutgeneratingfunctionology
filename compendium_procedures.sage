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
    return FastPiecewise(new_pairs, merge=False)

def projected_sequential_merge(g=piecewise_function_from_breakpoints_and_values([0, 1/3, 5/12, 1/2, 5/6, 11/12, 1], \
                                                                                [0, 1, 1/2, 3/4, 1/4, 1/2, 0]), n=1):
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
    new_g = FastPiecewise(new_pairs, merge=False)
    h = (1 / (n+1)) * xi + (1 / (n+1)) * new_g
    return h
