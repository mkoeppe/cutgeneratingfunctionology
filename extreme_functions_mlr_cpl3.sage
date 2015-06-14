# Make sure current directory is in path.
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

def cpl3_function(r0, z1, o1, o2):
    """
    Construct a CPL3 function.

    Parameters:
        0 < r0 (real) < 1, 0 < z1 (real) <= (1-r0)/4, 0 <= o1, o2 (real)

    EXAMPLE::

        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
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
    if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4) & bool(0 <= o1) & bool(0 <= o2)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if not bool(((1-r0)/2-2*z1 >= 0) and bool(1/2-o1-o2 >= 0)):
        logging.info("Conditions for a CPL-3 function are NOT satisfied.")

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    if r0+2*z1 < 1-2*z1:
        assert all(bkpt[i] < bkpt[i+1] for i in xrange(len(bkpt)-1)), "Constructed breakpoints are not in consecutive order"
        values = [0, 0, o1, o1+o2, 1-(o1+o2), 1-o1, 1]
        list_of_pairs = [[(bkpt[i], bkpt[i+1]), linear_function_through_points([bkpt[i], values[i]], [bkpt[i+1], values[i+1]])] for i in range(len(bkpt)-1)]
        return PiecewiseQuasiPeriodic(list_of_pairs)
    else:
        assert 0 < r0 < r0+z1 < r0+2*z1 == 1-2*z1 < 1-z1 < 1, "Constructed breakpoints are not in consecutive order"
        bkpt = [0, r0, r0+z1, r0+2*z1, 1-z1, 1]
        values = [0, 0, o1, o1+o2, 1-o1, 1]
        list_of_pairs = [[(bkpt[i], bkpt[i+1]), linear_function_through_points([bkpt[i], values[i]], [bkpt[i+1], values[i+1]])] for i in range(len(bkpt)-1)]
        return PiecewiseQuasiPeriodic(list_of_pairs)

def superadditive_lifting_function_from_group_function(fn, f=None):
    """
    Convert a standard representation 'phi' (a superadditive quasiperiodic function) from a group representation 'fn' (a subadditive periodic function).

    EXAMPLE::
        
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
        raise ValueError, "Bad parameter 'f'. Unable to construct the function."
    
    bkpt = fn._end_points
    phi_at_bkpt = [bkpt[i]-fn(bkpt[i])*f for i in range(len(bkpt))]

    list_of_pairs = []

    for i in range(len(bkpt)-1):
        list_of_pairs.append([(bkpt[i], bkpt[i+1]), linear_function_through_points([bkpt[i], phi_at_bkpt[i]], [bkpt[i+1], phi_at_bkpt[i+1]])])

    phi = PiecewiseQuasiPeriodic(list_of_pairs)

    return phi

def group_function_from_superadditive_lifting_function(phi, f=None):
    """
    Convert a group representation 'fn' (a subadditive periodic function) from a standard representation 'phi' (a superadditive quasiperiodic function).

    EXAMPLE::

        sage: phi = cpl3_function(r0=1/7, z1=1/7, o1=1/4, o2=1/12)
        sage: fn = group_function_from_superadditive_lifting_function(phi)
        sage: fn_expected = FastPiecewise([[(0, 1/7), FastLinearFunction(7, 0)], [(1/7, 2/7), FastLinearFunction(-21/4, 7/4)], [(2/7, 3/7), FastLinearFunction(35/12, -7/12)], [(3/7, 5/7), FastLinearFunction(-7/6, 7/6)], [(5/7, 6/7), FastLinearFunction(35/12, -7/4)], [(6/7, 1), FastLinearFunction(-21/4, 21/4)]])
        sage: fn == fn_expected
        True
        sage: q = plot(fn)
    """
    bkpt = phi._end_points
    if f is None:
        Min = phi(bkpt[0])-bkpt[0]
        for i in range(1,len(bkpt)):
            if phi(bkpt[i])-bkpt[i] < Min:

                Min = phi(bkpt[i])-bkpt[i]
                f = bkpt[i]

    if not bool(0 < f < 1):
        raise ValueError, "Bad parameter 'f'. Unable to construct the function."

    fn_at_bkpt = [(bkpt[i]-phi(bkpt[i]))/f for i in range(len(bkpt))]
    list_of_pairs = []
    for i in range(len(bkpt)-1):
        list_of_pairs.append([(bkpt[i], bkpt[i+1]), linear_function_through_points([bkpt[i], fn_at_bkpt[i]], [bkpt[i+1], fn_at_bkpt[i+1]])])

    fn = FastPiecewise(list_of_pairs)

    return fn

def mlr_cpl3_a_2_slope(r0=3/13, z1=3/26, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt a.
        - Proven extreme p.188, thm.19.

    Parameters:
        0 < r0 (real) < 1, 0 < z1 (real) <= (1-r0)/4
    Function is known to be extreme under the conditions:
        0 < r0 < 1, 0 <= z1 < 1

    Note:
        Parameter z1 is not actually used.
        Same as gmic(f=r0).

    Examples:
        page 183, Fig 2, point a::

            sage: logging.disable(logging.INFO) # Suppress output
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
    if not bool(0 < r0 < 1):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if conditioncheck:
        if not bool(0 <= z1 < 1):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    bkpt = [0, r0, 1]
    slopes = [1/r0, 1/(r0-1)]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_b_3_slope(r0=3/26, z1=1/13, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 3 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt b.
        - Proven extreme p.188, thm.19.

    Parameters:
        0 < r0 (real) < 1, 0 < z1 (real) <= (1-r0)/4
    Function is known to be extreme under the conditions:
        3*r0 + 8*z1 <= 1

    Note:
        When z1 = (1-r0)/4, the function is the same as gmic(f=r0).
        mlr_cpl3_b_3_slope(r0,z1) is the same as drlm_backward_3_slope(f=r0,bkpt=r0+2*z1).

    Examples:
        page 183, Fig 2, point b::

            sage: logging.disable(logging.INFO) # Suppress output
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
    if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if conditioncheck:
        if not (bool(3*r0 + 8*z1 <= 1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    bkpt = [0, r0, r0+2*z1, 1-2*z1, 1]
    if all(bkpt[i] < bkpt[i+1] for i in xrange(len(bkpt)-1)):
        slopes = [1/r0, (2*z1-1)/(2*z1*(1+r0)), 1/(1+r0), (2*z1-1)/(2*z1*(1+r0))]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)
    else:
        assert 0 < r0 < r0+2*z1 == 1-2*z1 < 1, "Constructed breakpoints are not in consecutive order"
        bkpt = [0, r0, 1]
        slopes = [1/r0, (2*z1-1)/(2*z1*(1+r0))]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)
        
def mlr_cpl3_c_3_slope(r0=5/24, z1=1/12, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 3 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt c.
        - Proven extreme p.188, thm.19.

    Parameters:
        0 < r0 (real) < 1, 0 < z1 (real) <= (1-r0)/4
    Function is known to be extreme under the conditions:
        3*r0 + 4*z1 <= 1

    Note:
        mlr_cpl3_c_3_slope(r0,z1) is the same as drlm_backward_3_slope(f=r0,bkpt=r0+z1).

    Examples:
        page 183, Fig 2, point c::

            sage: logging.disable(logging.INFO) # Suppress output
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
    if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if conditioncheck:
        if not bool(3*r0 + 4*z1 <= 1):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    bkpt = [0, r0, r0+z1, 1-z1, 1]
    slopes = [1/r0, (z1-1)/(z1*(1+r0)), 1/(1+r0), (z1-1)/(z1*(1+r0))]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)
   
def mlr_cpl3_d_3_slope(r0=1/6, z1=None, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 3 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt d.
        - Proven extreme p.188, thm.19.
       
    Parameters:
        0 < r0 (real) < 1, 0 < z1 (real) <= (1-r0)/4
    Function is known to be extreme under the conditions:
        r0 = 2*z1, r0 + 8*z1 <= 1

    Examples:
        p.183, Fig 2, point d1::

            sage: logging.disable(logging.INFO) # Suppress output
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
    if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if conditioncheck:
        if not (bool(r0 == 2*z1) & bool(r0+8*z1 <= 1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    if all(bkpt[i] < bkpt[i+1] for i in xrange(len(bkpt)-1)):
        slopes = [1/r0, (2*z1+1)/(2*z1*(r0-1)), (1-2*z1)/(2*z1*(1-r0)), 1/(r0-1), (1-2*z1)/(2*z1*(1-r0)), (2*z1+1)/(2*z1*(r0-1))]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)
    else:
        bkpt = [0, r0, r0+z1, 1-z1, 1]
        slopes = [1/r0, (2*z1+1)/(2*z1*(r0-1)), (1-2*z1)/(2*z1*(1-r0)), (2*z1+1)/(2*z1*(r0-1))]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_f_2_or_3_slope(r0=1/6, z1=None, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 2 or 3; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt f.
        - Proven extreme p.188, thm.19.

    Parameters:
        0 < r0 (real) < 1, 0 < z1 (real) <= (1-r0)/4
    Function is known to be extreme under the conditions:
        r0 <= z1, r0 + 5*z1 = 1

    Note:
        When z1 = (1-r0)/4, the function is the same as gmic(f=r0).

    Examples:
        page 184, Fig 3, point f1 and f2::

            sage: logging.disable(logging.INFO) # Suppress output
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
    if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if conditioncheck:
        if not (bool(r0 <= z1) & bool(r0+5*z1 == 1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    if all(bkpt[i] < bkpt[i+1] for i in xrange(len(bkpt)-1)):
        slopes = [1/r0, (z1*r0+8*z1^2-r0-6*z1+1)/(z1*(r0^2+4*z1+8*z1*r0-1)), (r0+8*z1-2)/(r0^2+4*z1+8*z1*r0-1), (r0+8*z1-1)/(r0^2+4*z1+8*z1*r0-1), (r0+8*z1-2)/(r0^2+4*z1+8*z1*r0-1), (z1*r0+8*z1^2-r0-6*z1+1)/(z1*(r0^2+4*z1+8*z1*r0-1))]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)
    else:
        bkpt = [0, r0, 1]
        slopes = [1/r0, (z1*r0+8*z1^2-r0-6*z1+1)/(z1*(r0^2+4*z1+8*z1*r0-1))]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_g_3_slope(r0=1/12, z1=None, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 3 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt g.
        - Proven extreme p.188, thm.19.

    Parameters:
        0 < r0 (real) < 1, 0 < z1 (real) <= (1-r0)/4
    Function is known to be extreme under the conditions:
        r0 < z1, 2*r0 + 4*z1 = 1

    Examples:
        page 184, Fig 3, point g::

            sage: logging.disable(logging.INFO) # Suppress output
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
    if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if conditioncheck:
        if not (bool(r0 < z1) & bool(2*r0 + 4*z1 == 1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    if all(bkpt[i] < bkpt[i+1] for i in xrange(len(bkpt)-1)):
        slopes = [1/r0, 3/(3*r0-1), (1-3*z1)/(z1*(1-3*r0)), 3/(3*r0-1), (1-3*z1)/(z1*(1-3*r0)), 3/(3*r0-1)]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)
    else:
        bkpt = [0, r0, r0+z1, 1-z1, 1]
        slopes = [1/r0, 3/(3*r0-1), (1-3*z1)/(z1*(1-3*r0)), 3/(3*r0-1)]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_h_2_slope(r0=1/4, z1=1/6, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt h.
        - Proven extreme p.188, thm.19.

    Parameters:
        0 < r0 (real) < 1, 0 < z1 (real) <= (1-r0)/4
    Function is known to be extreme under the conditions:
        r0 + 4*z1 <= 1 < 2*r0 + 4*z1

    Note:
        When z1 = (1-r0)/4, the function is the same as gmic(f=r0).

    Examples:
        page 183, Fig 2, point h::

            sage: logging.disable(logging.INFO) # Suppress output
            sage: h1 = mlr_cpl3_h_2_slope(r0=1/4, z1=1/6)
            sage: extremality_test(h1)
            True
            sage: phi = cpl3_function(r0=1/4, z1=1/6, o1=1/4, o2=1/4) 
            sage: fn = group_function_from_superadditive_lifting_function(phi)
            sage: h1 == fn
            True
            sage: h2 = mlr_cpl3_h_2_slope(r0=1/4, z1=3/16, conditioncheck=False)
            sage: h2 == gmic(f=1/4)
            True

    Reference:
        - [1] L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Infinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if conditioncheck:
        if not (bool(r0 + 4*z1 <= 1 < 2*r0 + 4*z1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    bkpt = [0, r0, r0+2*z1, 1-2*z1, 1]
    if all(bkpt[i] < bkpt[i+1] for i in xrange(len(bkpt)-1)):
        slopes = [1/r0, (4*z1-1)/(4*r0*z1), 1/r0, (4*z1-1)/(4*r0*z1)]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)
    else:
        bkpt = [0, r0, 1]
        slopes = [1/r0, (4*z1-1)/(4*r0*z1)]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_k_2_slope(r0=7/27, z1=4/27, field=None, conditioncheck=True): # phi_end_points for group function
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt k.
        - Proven extreme p.188, thm.19.

    Parameters:
        0 < r0 (real) < 1, 0 < z1 (real) <= (1-r0)/4
    Function is known to be extreme under the conditions:
        r0 <= 2*z1, r0 + 5*z1 = 1

    Examples:
        page 185, Fig 4, point k1::

            sage: logging.disable(logging.INFO) # Suppress output
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
    if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if conditioncheck:
        if not (bool(r0 <= 2*z1) & bool(r0+5*z1==1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    print bkpt
    if all(bkpt[i] < bkpt[i+1] for i in xrange(len(bkpt)-1)):
        slopes = [1/r0, (3*z1-1)/(3*z1*r0), 1/r0, (2-3*r0-12*z1)/(3*r0*(1-r0-4*z1)), 1/r0, (3*z1-1)/(3*z1*r0)]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)
    else:
        assert 0 < r0 < r0+z1 < r0+2*z1 == 1-2*z1 < 1-z1 < 1
        bkpt = [0, r0, r0+z1, 1-z1, 1]
        slopes = [1/r0, (3*z1-1)/(3*z1*r0), 1/r0, (3*z1-1)/(3*z1*r0)]
        ## print bkpt
        ## print slopes
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)
        
def mlr_cpl3_l_2_slope(r0=8/25, z1=None, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt l.
        - Proven extreme p.188, thm.19.

    Parameters:
        0 < r0 (real) < 1, 0 < z1 (real) <= (1-r0)/4
    Function is known to be extreme under the conditions:
        r0 = 2*z1, 6*z1 <= 1 < 7*z1
 
    Note:
        There is a typo in one of the given slopes [1] p.179, Table 3, Ext. pnt l, s3.
        The given slope is -4*z1/(2*z1-4*z1*r0) while the correct slope is (1-4*z1)/(2*z1-4*z1*r0).

    Examples:
        page 185, Fig 4, point l::

            sage: logging.disable(logging.INFO) # Suppress output
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
    if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if conditioncheck:
        if not (bool(r0 == 2*z1) & bool(6*z1<=1<7*z1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    if all(bkpt[i] < bkpt[i+1] for i in xrange(len(bkpt)-1)):
        slopes = [1/r0, 2/(2*r0-1), (1-4*z1)/(2*z1-4*z1*r0), 2/(2*r0-1), (1-4*z1)/(2*z1-4*z1*r0), 2/(2*r0-1)]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)
    else:
        bkpt = [0, r0, r0+z1, 1-z1, 1]
        slopes = [1/r0, 2/(2*r0-1), (1-4*z1)/(2*z1-4*z1*r0), 2/(2*r0-1)]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_n_3_slope(r0=16/22, z1=1/22, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 3 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt n.
        - Proven extreme p.188, thm.19.

    Parameters:
        0 < r0 (real) < 1, 0 < z1 (real) <= (1-r0)/4
    Function is known to be extreme under the conditions:
        r0 > 2*z1, r0 + 8*z1 < 1

    Examples:
        page 185, Fig 4, point n2::

            sage: logging.disable(logging.INFO) # Suppress output
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
    if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if conditioncheck:
        if not (bool(r0 > 2*z1) & bool(r0+8*z1<1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    if all(bkpt[i] < bkpt[i+1] for i in xrange(len(bkpt)-1)):
        slopes = [1/r0, (1+r0)/(r0*(r0-1)), 1/r0, 1/(r0-1), 1/r0, (1+r0)/(r0*(r0-1))]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)
    else:
        bkpt = [0, r0, r0+z1, 1-z1, 1]
        slopes = [1/r0, (1+r0)/(r0*(r0-1)), 1/r0, (1+r0)/(r0*(r0-1))]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_o_2_slope(r0=3/8, z1=None, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt o.
        - Proven extreme p.188, thm.19.

    Parameters:
        0 < r0 (real) < 1, 0 < z1 (real) <= (1-r0)/4
    Function is known to be extreme under the conditions:
        r0 >= 2*z1, 2*r0 + 2*z1 = 1

    Examples:
        page 186, Fig 5, point o::

            sage: logging.disable(logging.INFO) # Suppress output
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
    if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if conditioncheck:
        if not (bool(r0 >= 2*z1) & bool(2*r0+2*z1==1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    if all(bkpt[i] < bkpt[i+1] for i in xrange(len(bkpt)-1)):
        slopes = [1/r0, 2/(2*r0-1), (4*z1-4*z1*r0-1+2*r0)/(2*r0*z1*(1-2*r0)), 1/r0, (4*z1-4*z1*r0-1+2*r0)/(2*r0*z1*(1-2*r0)), 2/(2*r0-1)]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)
    else:
        bkpt = [0, r0, r0+z1, 1-z1, 1]
        slopes = [1/r0, 2/(2*r0-1), (4*z1-4*z1*r0-1+2*r0)/(2*r0*z1*(1-2*r0)), 2/(2*r0-1)]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)
        
def mlr_cpl3_p_2_slope(r0=5/12, z1=None, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt p.
        - Proven extreme p.188, thm.19.

    Parameters:
        0 < r0 (real) < 1, 0 < z1 (real) <= (1-r0)/4
    Function is known to be extreme under the conditions:
        r0 > 2*z1, 2*r0 + 2*z1 = 1

    Note:
        There is a typo in one of the given slopes [1] p.179, Table 3, Ext. pnt p, s4.
        The given slope is (2*z1-10*z1*r0+r0)/(r0*(1-2*r0)*(4*z1-1+r0)) while the correct slope is (-r0+2*r0^2-2*z1+8*z1*r0)/(r0*(1-2*r0)*(1-r0-4*z1)).

    Examples:
        page 186, Fig 5, point p1 and p2::

            sage: logging.disable(logging.INFO) # Suppress output
            sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
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
    if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if conditioncheck:
        if not (bool(r0 > 2*z1) & bool(2*r0+2*z1==1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    if all(bkpt[i] < bkpt[i+1] for i in xrange(len(bkpt)-1)):
        slopes = [1/r0, 2/(2*r0-1), 1/r0, (-r0+2*r0^2-2*z1+8*z1*r0)/(r0*(1-2*r0)*(1-r0-4*z1)), 1/r0, 2/(2*r0-1)]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)
    else:
        bkpt = [0, r0, r0+z1, 1-z1, 1]
        slopes = [1/r0, 2/(2*r0-1), 1/r0, 2/(2*r0-1)]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)
        
def mlr_cpl3_q_2_slope(r0=5/12, z1=3/24, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt q.
        - Proven extreme p.188, thm.19.


    Parameters:
        0 < r0 (real) < 1, 0 < z1 (real) <= (1-r0)/4
    Function is known to be extreme under the conditions:
        r0 > 2*z1, r0+4*z1 <= 1 < r0+5*z1

    Examples:
        page 186, Fig 5, point q::

            sage: logging.disable(logging.INFO) # Suppress output
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
    if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if conditioncheck:
        if not (bool(r0 > 2*z1) & bool(r0+4*z1 <= 1 < r0+5*z1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    if all(bkpt[i] < bkpt[i+1] for i in xrange(len(bkpt)-1)):
        slopes = [1/r0,(r0+2*z1)/(r0*(-1+r0+2*z1)), 1/r0, (r0+2*z1)/(r0*(-1+r0+2*z1)), 1/r0, (r0+2*z1)/(r0*(-1+r0+2*z1))]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)
    else:
        bkpt = [0, r0, r0+z1, 1-z1, 1]
        slopes = [1/r0, (r0+2*z1)/(r0*(-1+r0+2*z1)), 1/r0, (r0+2*z1)/(r0*(-1+r0+2*z1))]
        return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_r_2_slope(r0=3/7, z1=1/7, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Infinity; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered [1] p.179, Table 3, Ext. pnt r.
        - Proven extreme p.188, thm.19.

    Parameters:
        0 < r0 (real) < 1, 0 < z1 (real) <= (1-r0)/4
    Function is known to be extreme under the conditions:
        r0 > 2*z1, r0+4*z1 <= 1 <= 2*r0+2*z1

    Examples:
        page 185, Fig , point r::

            sage: logging.disable(logging.INFO) # Suppress output
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
    if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if conditioncheck:
        if not (bool(r0 > 2*z1) & bool(r0+4*z1<=1<=2*r0+2*z1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    bkpt = [0, r0, r0+z1, 1-z1, 1]
    slopes = [1/r0, (2*z1-1)/(2*r0*z1), 1/r0, (2*z1-1)/(2*r0*z1)]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)


