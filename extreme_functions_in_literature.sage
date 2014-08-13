import logging

logging.basicConfig(format='%(levelname)s: %(asctime)s %(message)s', level=logging.INFO)


def gmic(f=4/5):
    """
    Summary: 
        Name: GMIC
        Infinite (or Finite); Dim = 1; Slopes = 2; Continuous; Analysis of subadditive polytope method; 
        Discovered [55] p.7-8, Eq.8; 
        Proven (Infinite) [60] p.377, thm.3.3
                (Finite)  [57] p.514, Appendix 3.
    Parameters:
        f (real) \in (0,1)

    Examples: 
        [61] p.343, Fig. 1, Example 1 ::

            sage: h = gmic(4/5)	
		
    Reference: 
        [55]: R.E. Gomory, An algorithm for the mixed integer problem, Tech. Report RM-2597, RAND Corporation, 1960.

        [57]: R.E. Gomory, Some polyhedra related to combinatorial problems, Linear Algebra and its Application 2 (1969) 451-558.

        [60]: R.E. Gomory and E.L. Johnson, Some continuous functions related to corner polyhedra, part II, Mathematical Programming 3 (1972) 359–389.

        [61]: R.E. Gomory and E.L. Johnson, T-space and cutting planes, Mathematical Programming 96 (2003) 341–375.	
    """
    if not bool(0 < f < 1):
        raise ValueError, "Bad parameters. Unable to construct the function."
    gmi_bkpt = [0,f,1]
    gmi_values = [0,1,0]
    return piecewise_function_from_breakpoints_and_values(gmi_bkpt, gmi_values)


def gj_2_slope(f=3/5, lambda_1=1/6):
    """
    Summary: 
        Name: GJ's 2-Slope
        Infinite (or Finite); Dim = 1; Slopes = 2; Continuous; Analysis of subadditive polytope method; 
        Discovered [61] p.352, Fig.5, construction 1;  
        Proven (Infinite) [60] p.377, thm.3.3
                          [61] p.352, thm.4; p.354, thm.5
    Parameters:
        f (real) \in (0,1)
        lambda_1 (real) in (0,1]
		
    Function is known to be extreme under the conditions: 
        0 < lambda_1 <=1, lambda_1 < f/(1 - f)

    Examples:
        [61] p.354, Fig.6 ::

            sage: h = gj_2_slope(f=3/5, lambda_1=1/6)
            sage: h = gj_2_slope(f=3/5, lambda_1=1/2)
            sage: h = gj_2_slope(f=3/5, lambda_1=1)

    Reference: 
        [60]: R.E. Gomory and E.L. Johnson, Some continuous functions related to corner polyhedra, part II, Mathematical Programming 3 (1972) 359–389.

        [61]: R.E. Gomory and E.L. Johnson, T-space and cutting planes, Mathematical Programming 96 (2003) 341–375.	
    """
    if not (bool(0 < f < 1) & bool(0 < lambda_1 < f/(1 - f))):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if not (bool(lambda_1 <= 1)):
        logging.info("Conditions for extremality are NOT satisfied.")
    else:
        logging.info("Conditions for extremality are satisfied.")
    bkpts = [0, (f - lambda_1*(1 - f))/2, (f + lambda_1*(1 - f))/2, f, 1]
    values = [0, (1 + lambda_1)/2, (1 - lambda_1)/2, 1, 0]
    return piecewise_function_from_breakpoints_and_values(bkpts, values)


def gj_2_slope_repeat(f=3/5, s_positive=4, s_negative=-5, m=4, n=3):
    """
    Summary: 
        Name: GJ's 2-Slope-repeat
        Infinite (or Finite); Dim = 1; Slopes = 2; Continuous; Analysis of subadditive polytope method; 
        Discovered [61] p.354, Fig.7, construction 2; 
        Proven (Infinite) [60] p.377, thm.3.3
                          [61] p.354, thm.5; p.355, thm.6	
    Parameters:
        f (real) \in (0,1)
        s_positive, s_negative (real)
        m, n >= 2 (integer)

    Function is known to be extreme under the conditions: 
        0 < f < 1; 
        s_positive > 1/f;  s_negative < 1/(f - 1); 
        m >= (s_positive - s_positive*s_negative*f) / (s_positive - s_negative);  
        n >= (- s_negative + s_positive*s_negative*(f - 1)) / (s_positive - s_negative);

    Examples:
        [61] p.354, Fig.7 ::

            sage: h = gj_2_slope_repeat(f=3/5, s_positive=4, s_negative=-5, m=4, n=3)

    Reference: 
        [60]: R.E. Gomory and E.L. Johnson, Some continuous functions related to corner polyhedra, part II, Mathematical Programming 3 (1972) 359–389.

        [61]: R.E. Gomory and E.L. Johnson, T-space and cutting planes, Mathematical Programming 96 (2003) 341–375.	
    """
    if not (bool(0 < f < 1) & (m >= 2) & (n >= 2) & bool (s_positive > 1 / f) & bool(s_negative < 1/(f - 1))):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if not (bool(m >= (s_positive - s_positive*s_negative*f) / (s_positive - s_negative)) & bool(n >= (- s_negative + s_positive*s_negative*(f - 1)) / (s_positive - s_negative))):
        logging.info("Conditions for extremality are NOT satisfied.")
    else:
        logging.info("Conditions for extremality are satisfied.")
    len1_positive = (1 - s_negative*f) / (s_positive - s_negative) / m
    len1_negative = (f - m*len1_positive) / (m - 1)
    len2_negative = (1 - s_positive*(f - 1)) / (s_positive - s_negative) / n
    len2_positive = (1 - f - n*len2_negative) / (n - 1)
    interval_lengths = [len1_positive, len1_negative] * (m - 1) + [len1_positive, len2_negative] + [len2_positive, len2_negative]*(n - 1)
    slopes = [s_positive, s_negative]*(m + n - 1)
    return piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes)


def two_step_mir(f=4/5, alpha=3/10):
    """
    Summary: 
        Name: 2-Step MIR
        Infinite (or Finite); Dim = 1; Slopes = 2; Continuous; Simple sets method; 
        Discovered [33]  p.39 def.8, Fig.5
        Proven (Infinite) [60] p.377, thm.3.3

    Parameters:
        f (real) \in (0,1)
        alpha (real) \in (0,f)

	Function is known to be extreme under the conditions: 
        0 < alpha < f < 1; 
        f / alpha < ceil(f / alpha) <= 1 / alpha

    Examples:
        [33] p.40, Fig.5 ::

            sage: h = two_step_mir(f=4/5, alpha=3/10)

    Reference: 
        [33]: S. Dash and O. G¨unl¨uk, Valid inequalities based on simple mixed-integer sets., 
                Proceedings 10th Conference on Integer Programming and Combinatorial Optimization 
                (D. Bienstock and G. Nemhauser, eds.), Springer-Verlag, 2004, pp. 33–45.

        [60]: R.E. Gomory and E.L. Johnson, Some continuous functions related to corner polyhedra, part II, Mathematical Programming 3 (1972) 359–389.	
    """
    if not (bool(0 < alpha < f < 1) & bool(f / alpha < ceil(f / alpha) <= 1 / alpha)):
        logging.info("Conditions for extremality are NOT satisfied.")
    else:
        logging.info("Conditions for extremality are satisfied.")
    rho = f - alpha * floor(f / alpha)
    tau = ceil(f / alpha)
    s_positive = (1 - rho*tau) / (rho*tau*(1 - f))
    s_negative = - 1/(1 - f)
    interval_lengths = [rho, alpha - rho] * tau
    interval_lengths[-1] = 1 - f
    slopes = [s_positive, s_negative] * tau
    return piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes)


def interval_length_n_step_mir(n, m, a, b):
    if m == n:
        return [b[n - 1], a[n - 1] - b[n - 1]]
    else:
        l = interval_length_n_step_mir(n, m + 1, a, b)
        result = l * ceil(b[m - 1] / a[m]) 
        result[-1] = a[m - 1] - b[m - 1]
        return result


def n_step_mir(f=4/5, n=3, a=[1, 3/10, 8/100]):
    """
    Summary: 
        Name: n-Step MIR
        Infinite (or Finite); Dim = 1; Slopes = 2; Continuous; Simple sets method; 
        Discovered [74]  p.328, def.3, thm.2
        Proven (Infinite) [60] p.377, thm.3.3

    Parameters:
        f (real) \in (0,1)
        n (integer)
        a (list of reals, with length = n) \in (0,f)

	Function is known to be extreme under the conditions: 
        0 < a[1] < f < 1 == a[0];
        a[i] > 0, for i = 0, 1, ... , n-1. 
        b[i - 1] / a[i] < ceil(b[i - 1] / a[i]) <= a[i - 1] / a[i],  for i = 1, 2, ... , n-1.
        where,
        b[0] = f,
        b[i] = b[i - 1] - a[i] * floor(b[i - 1] / a[i]),  for i = 1, 2, ... , n-1.

    Examples:
        [74] p.333 - p.335, Fig.1 - Fig.6 ::
            sage: h = n_step_mir(f=4/5, n=1, a=[1])
            sage: h = n_step_mir(f=4/5, n=2, a=[1, 3/10])
            sage: h = n_step_mir(f=4/5, n=3, a=[1, 3/10, 8/100])
            sage: h = n_step_mir(f=4/5, n=4, a=[1, 3/10, 8/100, 3/100])
            sage: h = n_step_mir(f=4/5, n=5, a=[1, 45/100, 2/10, 558/10000, 11/1000])
            sage: h = n_step_mir(f=4/5, n=6, a=[1, 48/100, 19/100, 8/100, 32/1000, 12/1000])

    Reference: 
        [60]: R.E. Gomory and E.L. Johnson, Some continuous functions related to corner polyhedra, part II, Mathematical Programming 3 (1972) 359–389.	

        [74]: K. Kianfar and Y. Fathi, Generalized mixed integer rounding valid inequalities: 
                Facets for infinite group polyhedra, Mathematical Programming 120 (2009) 313–346.
    """
    if not bool(0 < f < 1 == a[0]):
        raise ValueError, "Bad parameters. Unable to construct the function."
    b = []
    b.append(f)
    t = True
    for i in range(1, n):
        b.append(b[i - 1] - a[i] * floor(b[i - 1] / a[i]))
        if not (bool(a[i] > 0) & bool(b[i - 1] / a[i] < ceil(b[i - 1] / a[i]) <= a[i - 1] / a[i])):
            t = False
            logging.info("Conditions for extremality are NOT satisfied.")
    if t:
        logging.info("Conditions for extremality are satisfied.") 
    interval_lengths =  interval_length_n_step_mir(n, 1, a, b)
    nb_interval = len(interval_lengths)
    interval_length_positive = sum(interval_lengths[i] for i in range(0, nb_interval, 2))
    interval_length_negative = sum(interval_lengths[i] for i in range(1, nb_interval, 2))
    s_negative = a[0] /(b[0] - a[0])
    s_positive = - s_negative * interval_length_negative / interval_length_positive 
    slopes = [s_positive, s_negative] * (nb_interval // 2)
    return piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes)


def forward_3_slope(f=4/5, lambda_1=2/9, lambda_2=1/3):
    """
    Summary: 
        Name: Forward 3-Slope
        Infinite (or Finite); Dim = 1; Slopes = 3; Continuous; Analysis of subadditive polytope method; 
        Discovered [61] p.359, Construction.3, Fig.8; 
        Proven [61] p.359, thm.8.

    Parameters:
        f (real) \in (0,1)
        lambda_1, lambda_2 (real) \in (0,1)

	Function is known to be extreme under the conditions: 
        0 <= lambda_1 <= 1/2   and   0 <= lambda_2 <= 1  (in literature)

    Note: 
        Since the domain and range are in [0,1]. I think the conditions should be:
        (0 <= lambda_1 <= 1/2)  &  (0 <= lambda_2 <= 1 - lambda_1) & (0 < lambda_1 * f + lambda_2 * (f - 1) < lambda_1 * f < f / 2)

    Examples: 
        [61] p.360, Fig.8 ::
            sage: h = forward_3_slope(f=4/5, lambda_1=2/9, lambda_2=1/6)
            sage: h = forward_3_slope(f=4/5, lambda_1=2/9, lambda_2=1/3)
            sage: h = forward_3_slope(f=4/5, lambda_1=2/9, lambda_2=1/2)
        Try irrational case ::
            sage: h = forward_3_slope(f=sqrt(17)/5, lambda_1=sqrt(5)/9, lambda_2=1/sqrt(10))

    Reference: 
        [61]: R.E. Gomory and E.L. Johnson, T-space and cutting planes, Mathematical Programming 96 (2003) 341–375.	
    """
    if not bool(0 < f < 1):
        raise ValueError, "Bad parameters. Unable to construct the function."
    a = lambda_1 * f
    a1 = a + lambda_2 * (f - 1)
    if not bool(0 < a1 < a < f / 2):
        raise ValueError, "Bad parameters. Unable to construct the function."
    # note the discrepancy with the published literature
    if not (bool(0 <= lambda_1 <= 1/2) & bool(0 <= lambda_2 <= 1 - lambda_1)):
        logging.info("Conditions for extremality are NOT satisfied.")
    else:
        logging.info("Conditions for extremality are satisfied.") 
    bkpts = [0, a1, a, f - a, f - a1, f, 1]
    values = [0, lambda_1 + lambda_2, lambda_1, 1 - lambda_1, 1 - lambda_1 - lambda_2, 1, 0]
    return piecewise_function_from_breakpoints_and_values(bkpts, values)


def backward_3_slope(f=1/12, bkpt=2/12):
    """
    Summary: 
        Name: Backward 3-Slope
        Infinite; Dim = 1; Slopes = 3; Continuous; Group relations method; 
        Discovered [40] p.154 eq.5; 
        Proven [40] p.153 thm.6.

    Parameters:
        f, bkpt (rational number? or real) \in (0,1)

    Function is known to be extreme under the conditions: 
        f, bkpt are rational numbers(*) such that   
        f < bkpt < (1+f)/4 < 1

    Note: 
        In [61] p.374, Appendix C, p.360. Fig.10, they claim that 
        the function (named pi3(u)) is facet (= extreme?), without assuming f, bkpt being rational numbers. 
        Also, one can add a perturbation (zigzag) on the third slope as shown in Fig.10.
        Proofs are not provided, however.

    Examples:
        Finite group --> Example 3.8 in [8] p.386, 
        Infinite group --> Interpolation using Equation 5 from [40] p.154 ::

        sage: h = backward_3_slope(f=1/12, bkpt=2/12)

        sage: h = backward_3_slope(f=1/12, bkpt=3/12)

    Reference: 
        [40]: S.S. Dey, J.-P.P. Richard, Y. Li, and L.A. Miller, Extreme inequalities for infinite group problems., 
                Mathematical Programming 121 (2010) 145–170.

        [61]: R.E. Gomory and E.L. Johnson, T-space and cutting planes, Mathematical Programming 96 (2003) 341–375.	
    """
    # FIXME: (*) Must f and bkpt be rational numbers?  
    #       The proof is based on interpolation of finite cyclic group extreme functions(cf. [8]). Need rational number. Can we relax them to irrational numbers?
    if not bool(0 < f < bkpt < 1 + f - bkpt < 1):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if not ((f in QQ) & (bkpt in QQ) & bool(0 < f < bkpt < ((1 + f)/4) < 1)):
        logging.info("Conditions for extremality are NOT satisfied.")
    else:
        logging.info("Conditions for extremality are satisfied.")
    bkpts = [0, f, bkpt, 1 + f - bkpt, 1]
    # values = [0, 1, bkpt/(1 + f), (1 + f - bkpt)/(1 + f),0]
    slopes = [1/f, (1 + f - bkpt)/(1 + f)/(f - bkpt), 1/(1 + f), (1 + f - bkpt)/(1 + f)/(f - bkpt)]
    return piecewise_function_from_breakpoints_and_slopes(bkpts, slopes)


def gj_2_slope_limit(f=0.6, nb_pieces_left=3, nb_pieces_right=4):
    """
    Summary: 
        Name: GJ's 2-Slope Limit
        Infinite; Dim = 1; Slopes = 1; Discontinuous; Group relations method; 
        Discovered [40] p.158 def.10; 
        Proven [40] p.159 thm.8.

    Parameters:
        f (real) \in (0,1)
        nb_pieces_left (positive integer) : number of linear pieces to the left of f; 
        nb_pieces_right (positive integer) : number of linear pieces to the right of f;

    Function is known to be extreme under the conditions: 
        nb_pieces_left * (1-f) <= nb_pieces_right * f;

    Examples:
        [40] p.159 Fig.4 ::
            sage: h = gj_2_slope_limit(f=0.6, nb_pieces_left=3, nb_pieces_right=4)

    Reference:
        [40]: S.S. Dey, J.-P.P. Richard, Y. Li, and L.A. Miller, Extreme inequalities for infinite group problems., 
                Mathematical Programming 121 (2010) 145–170.	
    """
    m = nb_pieces_left
    d = nb_pieces_right	
    if not ((m in ZZ) & (d in ZZ) & (m >= 1) & (d >= 1) & bool(0 < f < 1)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if not bool(m*(1 - f) <= d*f):
        logging.info("Conditions for extremality are NOT satisfied.")
    else:
        logging.info("Conditions for extremality are satisfied.")
    s = (m + d)/((d + 1)*f - (m - 1)*(1 - f))
    delta_2 = (s - s*f + 1)/(d + 1)
    if m == 1:
        delta_1 = 0
    else:
        delta_1 = (s*f - 1)/(m - 1)
    def psi(x):
        if (x >= 0) & (x <= 1):
            if (x <= f):
                if (m * x / f) in ZZ:
                    return x / f
                else:
                    return s*x - floor(m * x / f)*delta_1
            else:
                if (d*(1 - x)/(1 - f)) in ZZ:
                    return (1 - x)/(1 - f)
                else:
                    return s*x - s*f + 1 - (d - floor(d*(1 - x)/(1 - f)))*delta_2
        else:
            raise ValueError, "outside domain"
    # FIXME: Need to develop code for discontinuous functions.
    logging.warn("This is a discontinuous function; code for handling discontinuous functions is not implemented yet.")
    return psi		
