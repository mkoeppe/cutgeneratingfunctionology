from six.moves import range
def cpl3_function(r0, z1, o1, o2):
    r"""
    Construct a CPL3= function.

    Parameters:
        * 0 < r0 (real) < 1, 
        * 0 < z1 (real) <= (1-r0)/4, 
        * 0 <= o1, o2 (real), o1 + o2 <= 1/2;
        * if z1 = (1-r0)/4, then o1 + o2 = 1/2.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO) # Suppress output in automatic tests.
        sage: p = cpl3_function(r0=1/7, z1=1/7, o1=1/4, o2=1/12)
        sage: p
        <FastPiecewise with 6 parts, 
         (0, 1/7)       <FastLinearFunction 0>   values: [0, 0]
         (1/7, 2/7)     <FastLinearFunction 7/4*x - 1/4>         values: [0, 1/4]
         (2/7, 3/7)     <FastLinearFunction 7/12*x + 1/12>       values: [1/4, 1/3]
         (3/7, 5/7)     <FastLinearFunction 7/6*x - 1/6>         values: [1/3, 2/3]
         (5/7, 6/7)     <FastLinearFunction 7/12*x + 1/4>        values: [2/3, 3/4]
         (6/7, 1)       <FastLinearFunction 7/4*x - 3/4>         values: [3/4, 1]>
        sage: q = plot(p)

    Reference:
        - [1] L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Infinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)) or (bool(z1 == (1-r0)/4) & bool(o1 + o2 != 1/2)):
        raise ValueError("Bad parameters. Unable to construct the function.")
    if not (bool(0 <= o1) & bool(0 <= o2) & bool(o1+o2 <= 1/2)):
        logging.info("Conditions for a CPL-3 function are NOT satisfied.")

    if z1 < (1-r0)/4:
        bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
        values = [0, 0, o1, o1+o2, 1-(o1+o2), 1-o1, 1]
    else:
        bkpt = [0, r0, r0+z1, 1-z1, 1]
        values = [0, 0, o1, 1-o1, 1]
    list_of_pairs = [[(bkpt[i], bkpt[i+1]), linear_function_through_points([bkpt[i], values[i]], [bkpt[i+1], values[i+1]])] for i in range(len(bkpt)-1)]
    return PiecewiseQuasiPeriodic(list_of_pairs)

def superadditive_lifting_function_from_group_function(fn, f=None):
    r"""
    Convert a standard representation 'phi' (a superadditive quasiperiodic function) from a group representation 'fn' (a subadditive periodic function).

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO) # Suppress output in automatic tests.
        sage: fn = mlr_cpl3_d_3_slope(r0=1/7, z1=1/7)
        sage: phi = superadditive_lifting_function_from_group_function(fn)
        sage: phi_expected = FastPiecewise([[(0, 1/7), FastLinearFunction(0,0)], [(1/7, 2/7), FastLinearFunction(7/4, -1/4)], [(2/7, 3/7), FastLinearFunction(7/12, 1/12)], [(3/7, 5/7), FastLinearFunction(7/6, -1/6)], [(5/7, 6/7), FastLinearFunction(7/12, 1/4)], [(6/7, 1), FastLinearFunction(7/4, -3/4)]])
        sage: phi == phi_expected
        True
        sage: q = plot(phi)
    """
    if f is None:
        f = find_f(fn)
    if not bool(0 < f < 1):
        raise ValueError("Bad parameter 'f'. Unable to construct the function.")
    
    bkpt = fn._end_points
    phi_at_bkpt = [bkpt[i]-fn(bkpt[i])*f for i in range(len(bkpt))]
    list_of_pairs = []
    for i in range(len(bkpt)-1):
        list_of_pairs.append([(bkpt[i], bkpt[i+1]), linear_function_through_points([bkpt[i], phi_at_bkpt[i]], [bkpt[i+1], phi_at_bkpt[i+1]])])
    phi = PiecewiseQuasiPeriodic(list_of_pairs)
    return phi

def group_function_from_superadditive_lifting_function(phi, f=None):
    r"""
    Convert a group representation 'fn' (a subadditive periodic function) from a standard representation 'phi' (a superadditive quasiperiodic function).

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO) # Suppress output in automatic tests.
        sage: phi = cpl3_function(r0=1/7, z1=1/7, o1=1/4, o2=1/12)
        sage: fn = group_function_from_superadditive_lifting_function(phi)
        sage: fn_expected = FastPiecewise([[(0, 1/7), FastLinearFunction(7, 0)], [(1/7, 2/7), FastLinearFunction(-21/4, 7/4)], [(2/7, 3/7), FastLinearFunction(35/12, -7/12)], [(3/7, 5/7), FastLinearFunction(-7/6, 7/6)], [(5/7, 6/7), FastLinearFunction(35/12, -7/4)], [(6/7, 1), FastLinearFunction(-21/4, 21/4)]])
        sage: fn == fn_expected
        True
        sage: q = plot(fn)
    """
    bkpt = phi._end_points
    if f is None:
        m = phi(bkpt[0])-bkpt[0]
        for i in range(1,len(bkpt)):
            if phi(bkpt[i])-bkpt[i] < m:

                m = phi(bkpt[i])-bkpt[i]
                f = bkpt[i]
    if not bool(0 < f < 1):
        raise ValueError("Bad parameter 'f'. Unable to construct the function.")

    fn_at_bkpt = [(bkpt[i]-phi(bkpt[i]))/f for i in range(len(bkpt))]
    list_of_pairs = []
    for i in range(len(bkpt)-1):
        list_of_pairs.append([(bkpt[i], bkpt[i+1]), linear_function_through_points([bkpt[i], fn_at_bkpt[i]], [bkpt[i+1], fn_at_bkpt[i+1]])])
    fn = FastPiecewise(list_of_pairs)
    return fn

def mlr_cpl3_a_2_slope(r0=3/13, z1=3/26, field=None, conditioncheck=True):
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = mlr_cpl3_a_2_slope()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt a.
        - Proven extreme p.188, thm.19.

    Parameters:
        * 0 < r0 (real) < 1, 
        * 0 < z1 (real) <= (1-r0)/4

    Function is known to be extreme under the conditions:
        * 0 < r0 < 1, 
        * 0 <= z1 < 1

    Note:
        Parameter z1 is not actually used.
        Same as ``gmic(f=r0)``.

    Examples:
        page 183, Fig 2, point a::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO) # Suppress output in automatic tests.
            sage: h = mlr_cpl3_a_2_slope(r0=3/13, z1=3/26)
            sage: extremality_test(h)
            True
            sage: phi = cpl3_function(r0=3/13, z1=3/26, o1=3/20, o2=3/20)
            sage: fn = group_function_from_superadditive_lifting_function(phi)
            sage: h == fn
            True
            sage: h == gmic(f=3/13)
            True

    Reference:
        - [1] L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Infinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not bool(0 < r0 < 1):
            raise ValueError("Bad parameters. Unable to construct the function.")
        if not bool(0 <= z1 < 1):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    bkpt = [0, r0, 1]
    slopes = [1/r0, 1/(r0-1)]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=field)

def mlr_cpl3_b_3_slope(r0=3/26, z1=1/13, field=None, conditioncheck=True):
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = mlr_cpl3_b_3_slope()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 3 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt b.
        - Proven extreme p.188, thm.19.

    Parameters:
        * 0 < r0 (real) < 1, 
        * 0 < z1 (real) <= (1-r0)/4.

    Function is known to be extreme under the conditions:
        * 3*r0 + 8*z1 <= 1

    Note:
        * When z1 = (1-r0)/4, the function is the same as ``gmic(f=r0)``.
        * ``mlr_cpl3_b_3_slope(r0,z1)`` is the same as ``drlm_backward_3_slope(f=r0,bkpt=r0+2*z1)``.

    Examples:
        page 183, Fig 2, point b::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO) # Suppress output in automatic tests.
            sage: h1 = mlr_cpl3_b_3_slope(r0=3/26, z1=1/13)
            sage: extremality_test(h1)
            True
            sage: phi = cpl3_function(r0=3/26, z1=1/13, o1=7/58, o2=7/58)
            sage: fn = group_function_from_superadditive_lifting_function(phi)
            sage: h1 == fn
            True
            sage: h2 = drlm_backward_3_slope(f=3/26,bkpt=7/26)
            sage: extremality_test(h2)
            True
            sage: h1 == h2
            True
            sage: h3 = mlr_cpl3_b_3_slope(r0=3/26, z1=23/104, conditioncheck=False)
            sage: h3 == gmic(f=3/26)
            True
   
    Reference:
        - [1] L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Infinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
            raise ValueError("Bad parameters. Unable to construct the function.")
        if not (bool(3*r0 + 8*z1 <= 1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    if z1 < (1-r0)/4:
        bkpt = [0, r0, r0+2*z1, 1-2*z1, 1]
        slopes = [1/r0, (2*z1-1)/(2*z1*(1+r0)), 1/(1+r0), (2*z1-1)/(2*z1*(1+r0))]
    else:
        bkpt = [0, r0, 1]
        slopes = [1/r0, (2*z1-1)/(2*z1*(1+r0))]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=field)
        
def mlr_cpl3_c_3_slope(r0=5/24, z1=1/12, field=None, conditioncheck=True):
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = mlr_cpl3_c_3_slope()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 3 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt c.
        - Proven extreme p.188, thm.19.

    Parameters:
        * 0 < r0 (real) < 1, 
        * 0 < z1 (real) <= (1-r0)/4

    Function is known to be extreme under the conditions:
        * 3*r0 + 4*z1 <= 1

    Note:
        * ``mlr_cpl3_c_3_slope(r0,z1)`` is the same as ``drlm_backward_3_slope(f=r0,bkpt=r0+z1)``.

    Examples:
        page 183, Fig 2, point c::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO) # Suppress output in automatic tests.
            sage: h1 = mlr_cpl3_c_3_slope(r0=5/24, z1=1/12)
            sage: extremality_test(h1)
            True
            sage: phi = cpl3_function(r0=5/24, z1=1/12, o1=7/29, o2=2/29)
            sage: fn = group_function_from_superadditive_lifting_function(phi)
            sage: h1 == fn
            True
            sage: h2 = drlm_backward_3_slope(f=5/24,bkpt=7/24)
            sage: extremality_test(h2)
            True
            sage: h1 == h2
            True

    Reference:
        - [1] L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Infinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
            raise ValueError("Bad parameters. Unable to construct the function.")
        if not bool(3*r0 + 4*z1 <= 1):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    bkpt = [0, r0, r0+z1, 1-z1, 1]
    slopes = [1/r0, (z1-1)/(z1*(1+r0)), 1/(1+r0), (z1-1)/(z1*(1+r0))]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=field)
   
def mlr_cpl3_d_3_slope(r0=1/6, z1=None, field=None, conditioncheck=True):
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = mlr_cpl3_d_3_slope()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 3 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt d.
        - Proven extreme p.188, thm.19.
       
    Parameters:
        * 0 < r0 (real) < 1, 
        * 0 < z1 (real) <= (1-r0)/4

    Function is known to be extreme under the conditions:
        * r0 = 2*z1, 
        * r0 + 8*z1 <= 1

    Note:
        * ``multiplicative_homomorphism(mlr_cpl3_d_3_slope(r0, z1=2*r0), -1)`` is the same as ``gj_forward_3_slope(f=1-r0, lambda_1=2*z1/(1-r0), lambda_2=z1/r0)``;
        * ``gj_forward_3_slope`` being extreme only requires:  r0 >= z1, r0 + 4*z1 <= 1.

    Examples:
        p.183, Fig 2, point d1::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO) # Suppress output in automatic tests.
            sage: h1 = mlr_cpl3_d_3_slope(r0=1/6, z1=1/12)
            sage: extremality_test(h1)
            True
            sage: phi = cpl3_function(r0=1/6, z1=1/12, o1=1/5, o2=0) 
            sage: fn1 = group_function_from_superadditive_lifting_function(phi)
            sage: h1 == fn1
            True
            sage: h2 = mlr_cpl3_d_3_slope(r0=1/6, z1=5/24, conditioncheck=False)
            sage: extremality_test(h2)
            False
            sage: phi = cpl3_function(r0=1/6, z1=5/24, o1=7/20, o2=3/20)
            sage: fn2 = group_function_from_superadditive_lifting_function(phi)
            sage: h2 == fn2
            True

    Reference:
        - [1] L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Infinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if z1 is None:
        z1 = r0/2
    if conditioncheck:
        if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
            raise ValueError("Bad parameters. Unable to construct the function.")
        if not (bool(r0 == 2*z1) & bool(r0+8*z1 <= 1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    if z1 < (1-r0)/4:
        bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
        slopes = [1/r0, (2*z1+1)/(2*z1*(r0-1)), (1-2*z1)/(2*z1*(1-r0)), 1/(r0-1), (1-2*z1)/(2*z1*(1-r0)), (2*z1+1)/(2*z1*(r0-1))]
    else:
        bkpt = [0, r0, r0+z1, 1-z1, 1]
        slopes = [1/r0, (2*z1+1)/(2*z1*(r0-1)), (1-2*z1)/(2*z1*(1-r0)), (2*z1+1)/(2*z1*(r0-1))]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=field)

def mlr_cpl3_f_2_or_3_slope(r0=1/6, z1=None, field=None, conditioncheck=True):
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = mlr_cpl3_f_2_or_3_slope()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 2 or 3; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt f.
        - Proven extreme p.188, thm.19.

    Parameters:
        * 0 < r0 (real) < 1, 
        * 0 < z1 (real) <= (1-r0)/4

    Function is known to be extreme under the conditions:
        * r0 <= z1, 
        * r0 + 5*z1 = 1

    Note:
        When z1 = (1-r0)/4, the function is the same as ``gmic(f=r0)``.

    Examples:
        page 184, Fig 3, point f1 and f2::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO) # Suppress output in automatic tests.
            sage: h1 = mlr_cpl3_f_2_or_3_slope(r0=1/6, z1=1/6)
            sage: extremality_test(h1, f=1/6)
            True
            sage: phi = cpl3_function(r0=1/6, z1=1/6, o1=1/3, o2=0) 
            sage: fn = group_function_from_superadditive_lifting_function(phi, f=1/6)
            sage: h1 == fn
            True
            sage: h2 = mlr_cpl3_f_2_or_3_slope(r0=3/23, z1=4/23)
            sage: extremality_test(h2)
            True
            sage: h3 = mlr_cpl3_f_2_or_3_slope(r0=3/23, z1=5/23, conditioncheck=False)
            sage: h3 == gmic(f=3/23)
            True

    Reference:
        - [1] L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Infinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if z1 is None:
        z1 = (1-r0)/5
    if conditioncheck:
        if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
            raise ValueError("Bad parameters. Unable to construct the function.")
        if not (bool(r0 <= z1) & bool(r0+5*z1 == 1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    if z1 < (1-r0)/4:
        bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
        slopes = [1/r0, (z1*r0+8*z1^2-r0-6*z1+1)/(z1*(r0^2+4*z1+8*z1*r0-1)), (r0+8*z1-2)/(r0^2+4*z1+8*z1*r0-1), (r0+8*z1-1)/(r0^2+4*z1+8*z1*r0-1), (r0+8*z1-2)/(r0^2+4*z1+8*z1*r0-1), (z1*r0+8*z1^2-r0-6*z1+1)/(z1*(r0^2+4*z1+8*z1*r0-1))]
    else:
        bkpt = [0, r0, 1]
        slopes = [1/r0, (z1*r0+8*z1^2-r0-6*z1+1)/(z1*(r0^2+4*z1+8*z1*r0-1))]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=field)

def mlr_cpl3_g_3_slope(r0=1/12, z1=None, field=None, conditioncheck=True):
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = mlr_cpl3_g_3_slope()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 3 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt g.
        - Proven extreme p.188, thm.19.

    Parameters:
        * 0 < r0 (real) < 1, 
        * 0 < z1 (real) <= (1-r0)/4

    Function is known to be extreme under the conditions:
        * r0 < z1, 
        * 2*r0 + 4*z1 = 1

    Examples:
        page 184, Fig 3, point g::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO) # Suppress output in automatic tests.
            sage: h1 = mlr_cpl3_g_3_slope(r0=1/12, z1=5/24)
            sage: extremality_test(h1)
            True
            sage: phi = cpl3_function(r0=1/12, z1=5/24, o1=5/18, o2=1/6)
            sage: fn1 = group_function_from_superadditive_lifting_function(phi)
            sage: h1 == fn1
            True
            sage: h2 = mlr_cpl3_g_3_slope(r0=1/12, z1=11/48, conditioncheck=False)
            sage: extremality_test(h2)
            False
            sage: phi = cpl3_function(r0=1/12, z1=11/48, o1=11/36, o2=7/36)
            sage: fn2 = group_function_from_superadditive_lifting_function(phi)
            sage: h2 == fn2
            True

    Reference:
        - [1] L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Infinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if z1 is None:
        z1 = (1-2*r0)/4
    if conditioncheck:
        if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
            raise ValueError("Bad parameters. Unable to construct the function.")
        if not (bool(r0 < z1) & bool(2*r0 + 4*z1 == 1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")
    if z1 < (1-r0)/4:
        bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
        slopes = [1/r0, 3/(3*r0-1), (1-3*z1)/(z1*(1-3*r0)), 3/(3*r0-1), (1-3*z1)/(z1*(1-3*r0)), 3/(3*r0-1)]
    else:
        bkpt = [0, r0, r0+z1, 1-z1, 1]
        slopes = [1/r0, 3/(3*r0-1), (1-3*z1)/(z1*(1-3*r0)), 3/(3*r0-1)]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=field)

def mlr_cpl3_h_2_slope(r0=1/4, z1=1/6, field=None, conditioncheck=True):
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = mlr_cpl3_h_2_slope()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt h.
        - Proven extreme p.188, thm.19.

    Parameters:
        * 0 < r0 (real) < 1, 
        * 0 < z1 (real) <= (1-r0)/4

    Function is known to be extreme under the conditions:
        * r0 + 4*z1 <= 1 < 2*r0 + 4*z1

    Note:
        When z1 = (1-r0)/4, the function is the same as ``gmic(f=r0)``.

    Examples:
        page 183, Fig 2, point h::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO) # Suppress output in automatic tests.
            sage: h1 = mlr_cpl3_h_2_slope(r0=1/4, z1=1/6)
            sage: extremality_test(h1)
            True
            sage: phi = cpl3_function(r0=1/4, z1=1/6, o1=1/4, o2=1/4) 
            sage: fn = group_function_from_superadditive_lifting_function(phi)
            sage: h1 == fn
            True
            sage: h2 = mlr_cpl3_h_2_slope(r0=1/4, z1=3/16)
            sage: h2 == gmic(f=1/4)
            True

    Reference:
        - [1] L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Infinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
            raise ValueError("Bad parameters. Unable to construct the function.")
        if not (bool(r0 + 4*z1 <= 1 < 2*r0 + 4*z1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")
    if z1 < (1-r0)/4:
        bkpt = [0, r0, r0+2*z1, 1-2*z1, 1]
        slopes = [1/r0, (4*z1-1)/(4*r0*z1), 1/r0, (4*z1-1)/(4*r0*z1)]
    else:
        bkpt = [0, r0, 1]
        slopes = [1/r0, (4*z1-1)/(4*r0*z1)]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=field)

def mlr_cpl3_k_2_slope(r0=7/27, z1=4/27, field=None, conditioncheck=True):
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = mlr_cpl3_k_2_slope()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt k.
        - Proven extreme p.188, thm.19.

    Parameters:
        * 0 < r0 (real) < 1, 
        * 0 < z1 (real) <= (1-r0)/4

    Function is known to be extreme under the conditions:
        * r0 <= 2*z1, 
        * r0 + 5*z1 = 1

    Examples:
        page 185, Fig 4, point k1::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO) # Suppress output in automatic tests.
            sage: h = mlr_cpl3_k_2_slope(r0=7/27, z1=4/27)
            sage: extremality_test(h)
            True
            sage: phi = cpl3_function(r0=7/27, z1=4/27, o1=1/3, o2=0) 
            sage: fn = group_function_from_superadditive_lifting_function(phi)
            sage: h == fn
            True

    Reference:
        - [1] L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Infinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275

    """
    if z1 is None:
        z1 = (1-r0)/5
    if conditioncheck:
        if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
            raise ValueError("Bad parameters. Unable to construct the function.")
        if not (bool(r0 <= 2*z1) & bool(r0+5*z1==1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")
    if z1 < (1-r0)/4:
        bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
        slopes = [1/r0, (3*z1-1)/(3*z1*r0), 1/r0, (2-3*r0-12*z1)/(3*r0*(1-r0-4*z1)), 1/r0, (3*z1-1)/(3*z1*r0)]
    else:
        bkpt = [0, r0, r0+z1, 1-z1, 1]
        slopes = [1/r0, (3*z1-1)/(3*z1*r0), 1/r0, (3*z1-1)/(3*z1*r0)]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=field)
        
def mlr_cpl3_l_2_slope(r0=8/25, z1=None, field=None, conditioncheck=True):
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = mlr_cpl3_l_2_slope()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt l.
        - Proven extreme p.188, thm.19.

    Parameters:
        * 0 < r0 (real) < 1, 
        * 0 < z1 (real) <= (1-r0)/4

    Function is known to be extreme under the conditions:
        * r0 = 2*z1, 
        * 6*z1 <= 1 < 7*z1
 
    Note:
        There is a typo in one of the given slopes [1] p.179, Table 3, Ext. pnt l, s3.
        The given slope is -4*z1/(2*z1-4*z1*r0) while the correct slope is (1-4*z1)/(2*z1-4*z1*r0).

    Examples:
        page 185, Fig 4, point l::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO) # Suppress output in automatic tests.
            sage: h1 = mlr_cpl3_l_2_slope(r0=8/25, z1=4/25)
            sage: extremality_test(h1)
            True 
            sage: phi = cpl3_function(8/25,4/25,4/9,0)
            sage: fn1 = group_function_from_superadditive_lifting_function(phi)
            sage: h1 == fn1
            True
            sage: h2 = mlr_cpl3_l_2_slope(r0=8/25, z1=17/100, conditioncheck=False)
            sage: extremality_test(h2)
            False
            sage: phi = cpl3_function(r0=8/25, z1=17/100, o1=17/36, o2=1/36)
            sage: fn2 = group_function_from_superadditive_lifting_function(phi)
            sage: h2 == fn2
            True

    Reference:
        - [1] L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Infinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if z1 is None:
        z1 = r0/2
    if conditioncheck:
        if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
            raise ValueError("Bad parameters. Unable to construct the function.")
        if not (bool(r0 == 2*z1) & bool(6*z1<=1<7*z1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")
    if z1 < (1-r0)/4:
        bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
        slopes = [1/r0, 2/(2*r0-1), (1-4*z1)/(2*z1-4*z1*r0), 2/(2*r0-1), (1-4*z1)/(2*z1-4*z1*r0), 2/(2*r0-1)]
    else:
        bkpt = [0, r0, r0+z1, 1-z1, 1]
        slopes = [1/r0, 2/(2*r0-1), (1-4*z1)/(2*z1-4*z1*r0), 2/(2*r0-1)]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=field)

def mlr_cpl3_n_3_slope(r0=9/25, z1=2/25, field=None, conditioncheck=True):
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = mlr_cpl3_n_3_slope()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 3 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt n.
        - Proven extreme p.188, thm.19.

    Parameters:
        * 0 < r0 (real) < 1, 
        * 0 < z1 (real) <= (1-r0)/4

    Function is known to be extreme under the conditions:
        * r0 > 2*z1, 
        * r0 + 8*z1 <= 1

    Note:
        * ``multiplicative_homomorphism( mlr_cpl3_n_3_slope(r0, z1), -1) == gj_forward_3_slope(f=1-r0, lambda_1=4*z1/(1-r0), lambda_2=2*z1/r0)``;
        * ``gj_forward_3_slope`` being extreme only requires:  r0 >= z1, r0 + 4*z1 <= 1.

    Examples:
        page 185, Fig 4, point n2::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO) # Suppress output in automatic tests.
            sage: h1 = mlr_cpl3_n_3_slope(r0=9/25, z1=2/25)
            sage: extremality_test(h1)
            True
            sage: phi = cpl3_function(r0=9/25, z1=2/25, o1=1/4, o2=0) 
            sage: fn1 = group_function_from_superadditive_lifting_function(phi)
            sage: h1 == fn1
            True
            sage: h2 = mlr_cpl3_n_3_slope(r0=9/25, z1=4/25, conditioncheck=False)
            sage: extremality_test(h2)
            True
            sage: phi = cpl3_function(r0=9/25, z1=4/25, o1=1/2, o2=0)
            sage: fn2 = group_function_from_superadditive_lifting_function(phi)
            sage: h2 == fn2
            True

    Reference:
        - [1] L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Infinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
            raise ValueError("Bad parameters. Unable to construct the function.")
        if not (bool(r0 > 2*z1) & bool(r0+8*z1<=1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    if z1 < (1-r0)/4:
        bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
        slopes = [1/r0, (1+r0)/(r0*(r0-1)), 1/r0, 1/(r0-1), 1/r0, (1+r0)/(r0*(r0-1))]
    else:
        bkpt = [0, r0, r0+z1, 1-z1, 1]
        slopes = [1/r0, (1+r0)/(r0*(r0-1)), 1/r0, (1+r0)/(r0*(r0-1))]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=field)

def mlr_cpl3_o_2_slope(r0=3/8, z1=None, field=None, conditioncheck=True):
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = mlr_cpl3_o_2_slope()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt o.
        - Proven extreme p.188, thm.19.

    Parameters:
        * 0 < r0 (real) < 1, 
        * 0 < z1 (real) <= (1-r0)/4

    Function is known to be extreme under the conditions:
        * r0 >= 2*z1, 
        * 2*r0 + 2*z1 = 1

    Examples:
        page 186, Fig 5, point o::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO) # Suppress output in automatic tests.
            sage: h1 = mlr_cpl3_o_2_slope(r0=3/8, z1=1/8)
            sage: extremality_test(h1, f=3/8)
            True
            sage: phi = cpl3_function(r0=3/8, z1=1/8, o1=1/2, o2=0) 
            sage: fn1 = group_function_from_superadditive_lifting_function(phi)
            sage: h1 == fn1
            True
            sage: h2 = mlr_cpl3_o_2_slope(r0=2/8, z1=3/16,conditioncheck=False)
            sage: extremality_test(h2)
            False
            sage: phi = cpl3_function(r0=2/8, z1=3/16, o1=3/8, o2=1/8)
            sage: fn2 = group_function_from_superadditive_lifting_function(phi)
            sage: h2 == fn2
            True

    Reference:
        - [1] L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Infinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if z1 is None:
        z1 = (1-2*r0)/2
    if conditioncheck:
        if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
            raise ValueError("Bad parameters. Unable to construct the function.")
        if not (bool(r0 >= 2*z1) & bool(2*r0+2*z1==1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")
    if z1 < (1-r0)/4:
        bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
        slopes = [1/r0, 2/(2*r0-1), (4*z1-4*z1*r0-1+2*r0)/(2*r0*z1*(1-2*r0)), 1/r0, (4*z1-4*z1*r0-1+2*r0)/(2*r0*z1*(1-2*r0)), 2/(2*r0-1)]
    else:
        bkpt = [0, r0, r0+z1, 1-z1, 1]
        slopes = [1/r0, 2/(2*r0-1), (4*z1-4*z1*r0-1+2*r0)/(2*r0*z1*(1-2*r0)), 2/(2*r0-1)]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=field)
        
def mlr_cpl3_p_2_slope(r0=5/12, z1=None, field=None, conditioncheck=True):
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = mlr_cpl3_p_2_slope()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt p.
        - Proven extreme p.188, thm.19.

    Parameters:
        * 0 < r0 (real) < 1, 
        * 0 < z1 (real) <= (1-r0)/4

    Function is known to be extreme under the conditions:
        * r0 > 2*z1, 
        * 2*r0 + 2*z1 = 1

    Note:
        There is a typo in one of the given slopes [1] p.179, Table 3, Ext. pnt p, s4.
        The given slope is (2*z1-10*z1*r0+r0)/(r0*(1-2*r0)*(4*z1-1+r0)) while the correct slope is (-r0+2*r0^2-2*z1+8*z1*r0)/(r0*(1-2*r0)*(1-r0-4*z1)).

    Examples:
        page 186, Fig 5, point p1 and p2::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO) # Suppress output in automatic tests.
            sage: h1 = mlr_cpl3_p_2_slope(r0=5/12, z1=1/12)
            sage: extremality_test(h1, f=5/12)
            True
            sage: phi = cpl3_function(r0=5/12, z1=1/12, o1=1/2, o2=0)
            sage: fn = group_function_from_superadditive_lifting_function(phi)
            sage: h1 == fn
            True
            sage: h2 = mlr_cpl3_p_2_slope(r0=7/21, z1=1/6, conditioncheck=False)
            sage: extremality_test(h2, f=1/3)
            True
            sage: phi = cpl3_function(r0=7/21, z1=1/6, o1=1/2, o2=0)
            sage: fn2 = group_function_from_superadditive_lifting_function(phi)
            sage: h2 == fn2
            True

    Reference:
        - [1] L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Infinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if z1 is None:
        z1 = (1-2*r0)/2
    if conditioncheck:
        if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
            raise ValueError("Bad parameters. Unable to construct the function.")
        if not (bool(r0 > 2*z1) & bool(2*r0+2*z1==1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")
    if z1 < (1-r0)/4:
        bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
        slopes = [1/r0, 2/(2*r0-1), 1/r0, (-r0+2*r0^2-2*z1+8*z1*r0)/(r0*(1-2*r0)*(1-r0-4*z1)), 1/r0, 2/(2*r0-1)]
    else:
        bkpt = [0, r0, r0+z1, 1-z1, 1]
        slopes = [1/r0, 2/(2*r0-1), 1/r0, 2/(2*r0-1)]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=field)
        
def mlr_cpl3_q_2_slope(r0=5/12, z1=3/24, field=None, conditioncheck=True):
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = mlr_cpl3_q_2_slope()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt q.
        - Proven extreme p.188, thm.19.


    Parameters:
        * 0 < r0 (real) < 1, 
        * 0 < z1 (real) <= (1-r0)/4

    Function is known to be extreme under the conditions:
        * r0 > 2*z1, 
        * r0+4*z1 <= 1 < r0+5*z1

    Examples:
        page 186, Fig 5, point q::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO) # Suppress output in automatic tests.
            sage: h1 = mlr_cpl3_q_2_slope(r0=5/12, z1=3/24)
            sage: extremality_test(h1)
            True
            sage: phi = cpl3_function(r0=5/12, z1=3/24, o1=3/8, o2=0) 
            sage: fn1 = group_function_from_superadditive_lifting_function(phi)
            sage: h1 == fn1
            True
            sage: h2 = mlr_cpl3_q_2_slope(r0=5/12, z1=7/48,conditioncheck=False)
            sage: extremality_test(h2)
            True
            sage: phi = cpl3_function(r0=5/12, z1=7/48, o1=1/2, o2=0)
            sage: fn2 = group_function_from_superadditive_lifting_function(phi)
            sage: h2 == fn2
            True

    Reference:
        - [1] L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Infinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
            raise ValueError("Bad parameters. Unable to construct the function.")
        if not (bool(r0 > 2*z1) & bool(r0+4*z1 <= 1 < r0+5*z1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")
    if z1 < (1-r0)/4:
        bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
        slopes = [1/r0,(r0+2*z1)/(r0*(-1+r0+2*z1)), 1/r0, (r0+2*z1)/(r0*(-1+r0+2*z1)), 1/r0, (r0+2*z1)/(r0*(-1+r0+2*z1))]
    else:
        bkpt = [0, r0, r0+z1, 1-z1, 1]
        slopes = [1/r0, (r0+2*z1)/(r0*(-1+r0+2*z1)), 1/r0, (r0+2*z1)/(r0*(-1+r0+2*z1))]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=field)

def mlr_cpl3_r_2_slope(r0=3/7, z1=1/7, field=None, conditioncheck=True):
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = mlr_cpl3_r_2_slope()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt r.
        - Proven extreme p.188, thm.19.

    Parameters:
        * 0 < r0 (real) < 1, 
        * 0 < z1 (real) <= (1-r0)/4
    Function is known to be extreme under the conditions:
        * r0 > 2*z1, 
        * r0+4*z1 <= 1 <= 2*r0+2*z1

    Examples:
        page 185, Fig , point r::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO) # Suppress output in automatic tests.
            sage: h = mlr_cpl3_r_2_slope(r0=3/7, z1=1/7)
            sage: extremality_test(h)
            True 
            sage: phi = cpl3_function(r0=3/7, z1=1/7, o1=1/2, o2=0) 
            sage: fn = group_function_from_superadditive_lifting_function(phi)
            sage: h == fn
            True

    Reference:
        - [1] L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Infinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
            raise ValueError("Bad parameters. Unable to construct the function.")
        if not (bool(r0 > 2*z1) & bool(r0+4*z1<=1<=2*r0+2*z1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    bkpt = [0, r0, r0+z1, 1-z1, 1]
    slopes = [1/r0, (2*z1-1)/(2*r0*z1), 1/r0, (2*z1-1)/(2*r0*z1)]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=field)

def cpl3_0(f=3/5, z1=1/25, z2=1/125, field=None, conditioncheck=True):
    # This is the gmic function.
    claimed_parameter_attribute = None
    if conditioncheck:
        if all(bool(l < 0) for l in [-z2, -z1, -f, f + 2*z1 + 2*z2 - 1]):
            logging.info("Conditions for extremality are satisfied.")
            claimed_parameter_attribute = 'extreme'
        else:
            logging.info("Conditions for extremality are NOT satisfied.")
            claimed_parameter_attribute = 'not_extreme'
    theta = ((-z1)/(f - 1), (-z2)/(f - 1))
    h = cpl_n_group_function(3,False)(f, (z1,z2), theta, field=field)
    h._claimed_parameter_attribute = claimed_parameter_attribute
    return h

def cpl3_1(f=291403/720000, z1=114649/720000, z2=8501/144000, field=None, conditioncheck=True):
    # 2-slope, new?
    # different from (but looks like) automorphism(dg_2_step_mir(1-f, (1-f-z1)/2))
    claimed_parameter_attribute = None
    if conditioncheck:
        if all(bool(l < 0) for l in [-z2, -f + 2*z2, -f - 3*z1 - 2*z2 + 1, f + 2*z1 + 2*z2 - 1]):
            logging.info("Conditions for extremality are satisfied.")
            claimed_parameter_attribute = 'extreme'
        else:
            logging.info("Conditions for extremality are NOT satisfied.")
            claimed_parameter_attribute = 'not_extreme'
    theta = ((-z1)/(f + 2*z2 - 1), 0)
    h = cpl_n_group_function(3,False)(f, (z1,z2), theta, field=field)
    h._claimed_parameter_attribute = claimed_parameter_attribute
    return h

def cpl3_2(f=8011/24000, z1=1969/48000, z2=401/1920, field=None, conditioncheck=True):
    # 3-slope, is kzh_3_slope_param_extreme_1(f,z1,1-2*(z1+z2)), same extremality condition.
    claimed_parameter_attribute = None
    if conditioncheck:
        if all(bool(l < 0) for l in [-z1, -2*f - 2*z1 - 2*z2 + 1, f + 2*z1 + 2*z2 - 1, 2*f + 3*z1 + z2 - 1]):
            logging.info("Conditions for extremality are satisfied.")
            claimed_parameter_attribute = 'extreme'
        else:
            logging.info("Conditions for extremality are NOT satisfied.")
            claimed_parameter_attribute = 'not_extreme'
    theta = ((f + 3/2*z1 + z2 - 1/2)/(2*f + 3*z1 + 3*z2 - 1), \
             z2/(4*f + 6*z1 + 6*z2 - 2))
    h = cpl_n_group_function(3,False)(f, (z1,z2), theta, field=field)
    h._claimed_parameter_attribute = claimed_parameter_attribute
    return h

def cpl3_3(f=9469/66000, z1=2869/66000, z2=24709/66000, field=None, conditioncheck=True):
    # 3-slope, new.
    claimed_parameter_attribute = None
    if conditioncheck:
        if all(bool(l < 0) for l in [-f, -f - 3*z1 - 2*z2 + 1, f + 2*z1 + 2*z2 - 1, 2*f + 2*z1 - z2]):
            logging.info("Conditions for extremality are satisfied.")
            claimed_parameter_attribute = 'extreme'
        else:
            logging.info("Conditions for extremality are NOT satisfied.")
            claimed_parameter_attribute = 'not_extreme'
    theta = ((-f*z1 - z1^2)/(f^2 + f*z1 + 2*f*z2 - f - z1), \
             (-z1*z2)/(f^2 + f*z1 + 2*f*z2 - f - z1))
    h = cpl_n_group_function(3,False)(f, (z1,z2), theta, field=field)
    h._claimed_parameter_attribute = claimed_parameter_attribute
    return h

def cpl3_4(f=8497/126000, z1=499/42000, z2=27863/126000, field=None, conditioncheck=True):
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = cpl3_4()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    The function is extreme under the condtions:
    f + 4*z1 + 4*z2 - 1 <= 0 && -z1 < 0 && -f < 0 && 2*f + 2*z1 - z2 <= 0

    TESTS::

        sage: from cutgeneratingfunctionology.igp import *
        sage: f=1/18; z1=1/25; z2=7/36
        sage: h = cpl3_4(f,z1,z2)
        sage: extremality_test(h)
        True
    """
    # 4-slope, new.
    # looks like automorphism(gj_forward_3_slope), but has 2 positive slopes.
    claimed_parameter_attribute = None
    if conditioncheck:
        if all(bool(l < 0) for l in [-z1, -f, f + 4*z1 + 4*z2 - 1, 2*f + 2*z1 - z2]):
            logging.info("Conditions for extremality are satisfied.")
            claimed_parameter_attribute = 'extreme'
        else:
            logging.info("Conditions for extremality are NOT satisfied.")
            claimed_parameter_attribute = 'not_extreme'
    theta = ((-f*z1 - z1^2 - f*z2 - z1*z2)/(f^2 + f*z1 + f*z2 - f - z1 - z2), \
             (-z1*z2 - z2^2)/(f^2 + f*z1 + f*z2 - f - z1 - z2))
    h = cpl_n_group_function(3,False)(f, (z1,z2), theta, field=field)
    h._claimed_parameter_attribute = claimed_parameter_attribute
    return h

def cpl3_5(f=3/5, z1=1/25, z2=1/125, field=None, conditioncheck=True):
    # 2-slope, is automorphism(gj_2_slope(1-f, (1-f-2*z1)/f))
    claimed_parameter_attribute = None
    if conditioncheck:
        if all(bool(l < 0) for l in [-z2, -z1, -2*f - 2*z1 + 1, f + 2*z1 + 2*z2 - 1]):
            logging.info("Conditions for extremality are satisfied.")
            claimed_parameter_attribute = 'extreme'
        else:
            logging.info("Conditions for extremality are NOT satisfied.")
            claimed_parameter_attribute = 'not_extreme'
    theta = (1/2, 0)
    h = cpl_n_group_function(3,False)(f, (z1,z2), theta, field=field)
    h._claimed_parameter_attribute = claimed_parameter_attribute
    return h

def cpl3_6(f=1/2000, z1=499/4000, z2=1/4, field=None, conditioncheck=True):
    # 3-slope, is drlm_backward_3_slope(f,f+z1)
    claimed_parameter_attribute = None
    if conditioncheck:
        if all(bool(l < 0) for l in [-z2, -z1, -f, f + 2*z1 + 2*z2 - 1, 3*f + 4*z1 - 1]):
            logging.info("Conditions for extremality are satisfied.")
            claimed_parameter_attribute = 'extreme'
        else:
            logging.info("Conditions for extremality are NOT satisfied.")
            claimed_parameter_attribute = 'not_extreme'
    theta = ((f + z1)/(f + 1), z2/(f + 1))
    h = cpl_n_group_function(3,False)(f, (z1,z2), theta, field=field)
    h._claimed_parameter_attribute = claimed_parameter_attribute
    return h

def cpl3_7(f=3/5, z1=1/25, z2=1/125, field=None, conditioncheck=True):
     # 2-slope, is automorphism(gj_2_slope(1-f, (1-f-2*(z1+z2))/f))
    claimed_parameter_attribute = None
    if conditioncheck:
        if all(bool(l < 0) for l in [-z2, -z1, -2*f - 2*z1 - 2*z2 + 1, f + 2*z1 + 2*z2 - 1]):
            logging.info("Conditions for extremality are satisfied.")
            claimed_parameter_attribute = 'extreme'
        else:
            logging.info("Conditions for extremality are NOT satisfied.")
            claimed_parameter_attribute = 'not_extreme'
    theta = (z1/(2*z1 + 2*z2), z2/(2*z1 + 2*z2))
    h = cpl_n_group_function(3,False)(f, (z1,z2), theta, field=field)
    h._claimed_parameter_attribute = claimed_parameter_attribute
    return h

def cpl3_8(f=3/5, z1=1/25, z2=1/125, field=None, conditioncheck=True):
     # 3-slope, is automorphism(gj_forward_3_slope(1-f, 2*(z1+z2)/(1-f), 2*z2/f))
    claimed_parameter_attribute = None
    if conditioncheck:
        if all(bool(l < 0) for l in [-z2, -z1, -f + 2*z2, f + 4*z1 + 4*z2 - 1]):
            logging.info("Conditions for extremality are satisfied.")
            claimed_parameter_attribute = 'extreme'
        else:
            logging.info("Conditions for extremality are NOT satisfied.")
            claimed_parameter_attribute = 'not_extreme'
    theta = ((-z1 - z2)/(f - 1), 0)
    h = cpl_n_group_function(3,False)(f, (z1,z2), theta, field=field)
    h._claimed_parameter_attribute = claimed_parameter_attribute
    return h

def cpl3_9(f=20489/330000, z1=61603/660000, z2=72581/660000, field=None, conditioncheck=True):
     # 3-slope, is drlm_backward_3_slope(f,f+z1+z2)
    claimed_parameter_attribute = None
    if conditioncheck:
        if all(bool(l < 0) for l in [-z2, -z1, -f, 3*f + 4*z1 + 4*z2 - 1]):
            logging.info("Conditions for extremality are satisfied.")
            claimed_parameter_attribute = 'extreme'
        else:
            logging.info("Conditions for extremality are NOT satisfied.")
            claimed_parameter_attribute = 'not_extreme'
    theta = ((f*z1 + z1^2 + z1*z2)/(f*z1 + f*z2 + z1 + z2), \
             (f*z2 + z1*z2 + z2^2)/(f*z1 + f*z2 + z1 + z2))
    h = cpl_n_group_function(3,False)(f, (z1,z2), theta, field=field)
    h._claimed_parameter_attribute = claimed_parameter_attribute
    return h

##### One linear equation on f, z1, z2 #####

def cpl3_10(f=14/33, z1=17/264, z2=None, field=None, conditioncheck=True):
    # 2-slope, new?
    if z2 is None:
        z2 = f / 2
    claimed_parameter_attribute = None
    if conditioncheck:
        if bool(-f + 2*z2 == 0) and all(bool(l < 0) for l in [2*f + 2*z1 - 1, -2*f - 3*z1 + 1, -f]):
            logging.info("Conditions for extremality are satisfied.")
            claimed_parameter_attribute = 'extreme'
        else:
            logging.info("Conditions for extremality are NOT satisfied.")
            claimed_parameter_attribute = 'not_extreme'
    theta = ((-z1)/(2*f - 1), 0)
    h = cpl_n_group_function(3,False)(f, (z1,z2), theta, field=field)
    h._claimed_parameter_attribute = claimed_parameter_attribute
    return h

def cpl3_11(f=73/160, z1=None, z2=27/160, field=None, conditioncheck=True):
    # 2-slope, is multipicative_homomorphism(gmic(f), 2)
    if z1 is None:
        z1 = -f + 1/2
    claimed_parameter_attribute = None
    if conditioncheck:
        if bool(-2*f - 2*z1 + 1 == 0) and all(bool(l < 0) for l in [-f + 2*z2, 2*f - 1, -z2]):
            logging.info("Conditions for extremality are satisfied.")
            claimed_parameter_attribute = 'extreme'
        else:
            logging.info("Conditions for extremality are NOT satisfied.")
            claimed_parameter_attribute = 'not_extreme'
    theta = (1/2, 0)
    h = cpl_n_group_function(3,False)(f, (z1,z2), theta, field=field)
    h._claimed_parameter_attribute = claimed_parameter_attribute
    return h

def cpl3_12(f=5/42, z1=61/336, z2=None, field=None, conditioncheck=True):
    # 3-slope, new?
    if z2 is None:
        z2 = -2*f - 3*z1 + 1
    claimed_parameter_attribute = None
    if conditioncheck:
        if bool(-2*f - 3*z1 - z2 + 1 == 0) and all(bool(l < 0) for l in [-f, -3*f - 4*z1 + 1, 3*f + 3*z1 - 1]):
            logging.info("Conditions for extremality are satisfied.")
            claimed_parameter_attribute = 'extreme'
        else:
            logging.info("Conditions for extremality are NOT satisfied.")
            claimed_parameter_attribute = 'not_extreme'
    theta = (2*z1/(9*z1 + 3*z2 - 1), (3*z1 + 3*z2 - 1)/(9*z1 + 3*z2 - 1))
    h = cpl_n_group_function(3,False)(f, (z1,z2), theta, field=field)
    h._claimed_parameter_attribute = claimed_parameter_attribute
    return h


def cpl3_13(f=1/14, z1=5/28, z2=None, field=None, conditioncheck=True):
    # 3-slope, is automorphism(gj_forward_3_slope(1-f, 2*(z1+z2)/(1-f), 2*z2/f))
    if z2 is None:
        z2 = f/2
    claimed_parameter_attribute = None
    if conditioncheck:
        if bool(-f + 2*z2 == 0) and all(bool(l < 0) for l in [-f, 3*f + 4*z1 - 1, -z1]):
            logging.info("Conditions for extremality are satisfied.")
            claimed_parameter_attribute = 'extreme'
        else:
            logging.info("Conditions for extremality are NOT satisfied.")
            claimed_parameter_attribute = 'not_extreme'
    theta = ((-f - 2*z1)/(2*f - 2), 0)
    h = cpl_n_group_function(3,False)(f, (z1,z2), theta, field=field)
    h._claimed_parameter_attribute = claimed_parameter_attribute
    return h

def cpl3_14(f=1/24, z1=9/40, z2=None, field=None, conditioncheck=True):
    # 3-slope, new?
    if z2 is None:
        z2 = (-f - 3*z1 + 1) / 2
    claimed_parameter_attribute = None
    if conditioncheck:
        if bool(-f - 3*z1 - 2*z2 + 1 == 0) and all(bool(l < 0) for l in [-f, -3*f - 5*z1 + 1, 3*f + 3*z1 - 1]):
            logging.info("Conditions for extremality are satisfied.")
            claimed_parameter_attribute = 'extreme'
        else:
            logging.info("Conditions for extremality are NOT satisfied.")
            claimed_parameter_attribute = 'not_extreme'
    theta = (z1/(9*z1 + 6*z2 - 2), (3*z1 + 3*z2 - 1)/(9*z1 + 6*z2 - 2))
    h = cpl_n_group_function(3,False)(f, (z1,z2), theta, field=field)
    h._claimed_parameter_attribute = claimed_parameter_attribute
    return h
