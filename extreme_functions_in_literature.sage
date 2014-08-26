import logging

logging.basicConfig(format='%(levelname)s: %(asctime)s %(message)s', level=logging.INFO)


def gmic(f=4/5):
    """
    Summary:
        - Name: GMIC;
        - Infinite (or Finite); Dim = 1; Slopes = 2; Continuous; Analysis of subadditive polytope method;
        - Discovered [55] p.7-8, Eq.8;
        - Proven (Infinite) [60] p.377, thm.3.3; (Finite) [57] p.514, Appendix 3.

    Parameters:
        f (real) \in (0,1).

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
        - Name: GJ's 2-Slope;
        - Infinite (or Finite); Dim = 1; Slopes = 2; Continuous; Analysis of subadditive polytope method;
        - Discovered [61] p.352, Fig.5, construction 1;
        - Proven (Infinite) [60] p.377, thm.3.3; [61] p.352, thm.4; p.354, thm.5.

    Parameters:
        f (real) \in (0,1);
        lambda_1 (real) in (0,1].

    Function is known to be extreme under the conditions:
        0 < lambda_1 <=1, lambda_1 < f/(1 - f).

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
        - Name: GJ's 2-Slope-repeat;
        - Infinite (or Finite); Dim = 1; Slopes = 2; Continuous; Analysis of subadditive polytope method;
        - Discovered [61] p.354, Fig.7, construction 2;
        - Proven (Infinite) [60] p.377, thm.3.3; [61] p.354, thm.5; p.355, thm.6.

    Parameters:
        f (real) \in (0,1);
        s_positive, s_negative (real);
        m, n >= 2 (integer).

    Function is known to be extreme under the conditions:
        0 < f < 1;
        s_positive > 1/f;  s_negative < 1/(f - 1);
        m >= (s_positive - s_positive*s_negative*f) / (s_positive - s_negative);
        n >= (- s_negative + s_positive*s_negative*(f - 1)) / (s_positive - s_negative).

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
        - Name: 2-Step MIR;
        - Infinite (or Finite); Dim = 1; Slopes = 2; Continuous; Simple sets method;
        - Discovered [33]  p.39 def.8, Fig.5;
        - Proven (Infinite) [60] p.377, thm.3.3.

    Parameters:
        f (real) \in (0,1);
        alpha (real) \in (0,f).

    Function is known to be extreme under the conditions:
        0 < alpha < f < 1;
        f / alpha < ceil(f / alpha) <= 1 / alpha.

    Examples:
        [33] p.40, Fig.5 ::

            sage: h = two_step_mir(f=4/5, alpha=3/10)

    Reference:
        [33]: S. Dash and O. G¨unl¨uk, Valid inequalities based on simple mixed-integer sets.,
                Proceedings 10th Conference on Integer Programming and Combinatorial Optimization
                (D. Bienstock and G. Nemhauser, eds.), Springer-Verlag, 2004, pp. 33–45.

        [60]: R.E. Gomory and E.L. Johnson, Some continuous functions related to corner polyhedra, part II, Mathematical Programming 3 (1972) 359–389.
    """
    if not (bool(0 < alpha < f < 1) & bool(f / alpha < ceil(f / alpha))):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if not bool(ceil(f / alpha) <= 1 / alpha):
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


def n_step_mir(f=4/5, a=[1, 3/10, 8/100]):
    """
    Summary:
        - Name: n-Step MIR;
        - Infinite (or Finite); Dim = 1; Slopes = 2; Continuous; Simple sets method;
        - Discovered [74]  p.328, def.3, thm.2;
        - Proven (Infinite) [60] p.377, thm.3.3.

    Parameters:
        f (real) \in (0,1);
        a (list of reals, with length = n) \in (0,f).

    Function is known to be extreme under the conditions:
        0 < a[1] < f < 1 == a[0];
        a[i] > 0, for i = 0, 1, ... , n-1;
        b[i - 1] / a[i] < ceil(b[i - 1] / a[i]) <= a[i - 1] / a[i],  for i = 1, 2, ... , n-1;
        where,
        b[0] = f;
        b[i] = b[i - 1] - a[i] * floor(b[i - 1] / a[i]),  for i = 1, 2, ... , n-1.

    Note:
        if a[i] > b[i-1] for some i, then the n_step_mir function degenerates, i.e.
        n_step_mir(f, [a[0], .. , a[n - 1]]) = n_step_mir(f, [a[0], .. a[i - 1], a[i + 1], ... , a[n - 1]])

    Examples:
        [74] p.333 - p.335, Fig.1 - Fig.6 ::
        
            sage: h = n_step_mir(f=4/5, a=[1])
            sage: h = n_step_mir(f=4/5, a=[1, 3/10])
            sage: h = n_step_mir(f=4/5, a=[1, 3/10, 8/100])
            sage: h = n_step_mir(f=4/5, a=[1, 3/10, 8/100, 3/100])
            sage: h = n_step_mir(f=4/5, a=[1, 45/100, 2/10, 558/10000, 11/1000])
            sage: h = n_step_mir(f=4/5, a=[1, 48/100, 19/100, 8/100, 32/1000, 12/1000])

    Reference:
        [60]: R.E. Gomory and E.L. Johnson, Some continuous functions related to corner polyhedra, part II, Mathematical Programming 3 (1972) 359–389.

        [74]: K. Kianfar and Y. Fathi, Generalized mixed integer rounding valid inequalities:
                Facets for infinite group polyhedra, Mathematical Programming 120 (2009) 313–346.
    """
    if (a == []) | (not bool(0 < f < 1 == a[0])):
        raise ValueError, "Bad parameters. Unable to construct the function."
    b = []
    b.append(f)
    n = len(a)
    t = True
    for i in range(1, n):
        b.append(b[i - 1] - a[i] * floor(b[i - 1] / a[i]))
        if not (bool(0 < a[i]) & bool(b[i - 1] / a[i] < ceil(b[i - 1] / a[i]))):
            raise ValueError, "Bad parameters. Unable to construct the function."
        if not bool(ceil(b[i - 1] / a[i]) <= a[i - 1] / a[i]):
            t = False
    if t:
        logging.info("Conditions for extremality are satisfied.")
    else:
        logging.info("Conditions for extremality are NOT satisfied.")
    interval_lengths =  interval_length_n_step_mir(n, 1, a, b)
    nb_interval = len(interval_lengths)
    interval_length_positive = sum(interval_lengths[i] for i in range(0, nb_interval, 2))
    interval_length_negative = sum(interval_lengths[i] for i in range(1, nb_interval, 2))
    s_negative = a[0] /(b[0] - a[0])
    s_positive = - s_negative * interval_length_negative / interval_length_positive 
    slopes = [s_positive, s_negative] * (nb_interval // 2)
    return piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes)


def gj_forward_3_slope(f=4/5, lambda_1=2/9, lambda_2=1/3):
    """
    Summary: 
        - Name: GJ's Forward 3-Slope;
        - Infinite (or Finite); Dim = 1; Slopes = 3; Continuous; Analysis of subadditive polytope method;
        - Discovered [61] p.359, Construction.3, Fig.8;
        - Proven [61] p.359, thm.8.

    Parameters:
        f (real) \in (0,1);
        lambda_1, lambda_2 (real) \in (0,1).

    Function is known to be extreme under the conditions:
        0 <= lambda_1 <= 1/2;
        0 <= lambda_2 <= 1  (in literature).

    Note: 
        Since the domain and range are in [0,1]. I think the conditions should be:
        (0 <= lambda_1 <= 1/2)  &  (0 <= lambda_2 <= 1 - lambda_1) & (0 < lambda_1 * f + lambda_2 * (f - 1) < lambda_1 * f < f / 2).

    Examples:
        [61] p.360, Fig.8 ::
            sage: h = gj_forward_3_slope(f=4/5, lambda_1=2/9, lambda_2=1/6)
            sage: h = gj_forward_3_slope(f=4/5, lambda_1=2/9, lambda_2=1/3)
            sage: h = gj_forward_3_slope(f=4/5, lambda_1=2/9, lambda_2=1/2)
        Try irrational case ::
            sage: h = gj_forward_3_slope(f=sqrt(17)/5, lambda_1=sqrt(5)/9, lambda_2=1/sqrt(10))

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


def drlm_backward_3_slope(f=1/12, bkpt=2/12):
    """
    Summary:
        - Name: drlm's Backward 3-Slope;
        - Infinite; Dim = 1; Slopes = 3; Continuous; Group relations method;
        - Discovered [40] p.154 eq.5;
        - Proven [40] p.153 thm.6.

    Parameters:
        f, bkpt (real) \in (0,1).

    Function is known to be extreme under the conditions:
        f < bkpt < (1+f)/4 < 1.

    Note:
        In [40], they require that f, bkpt are rational numbers.
        The proof is based on interpolation of finite cyclic group extreme functions(cf. [8]), so it needs rational numbers.
        But in fact, by analysing covered intervals and using the condition f < bkpt < (1+f)/4 < 1,
        one can prove that the function is extreme without assuming f, bkpt being rational numbers.
        
        In [61] p.374, Appendix C, p.360. Fig.10, they consider real number f, bkpt, and claim (without proof) that:
        1) the function (named pi3(u)) is facet (thus extreme);
        2) can add a perturbation (zigzag) on the third slope as shown in Fig.10;

    Examples:
        Finite group --> Example 3.8 in [8] p.386,
        Infinite group --> Interpolation using Equation 5 from [40] p.154 ::

        sage: h = drlm_backward_3_slope(f=1/12, bkpt=2/12)

        sage: h = drlm_backward_3_slope(f=1/12, bkpt=3/12)

    Reference:
        [8]: J. Ar´aoz, L. Evans, R.E. Gomory, and E.L. Johnson, Cyclic groups and knapsack facets,
                Mathematical Programming 96 (2003) 377–408.
                
        [40]: S.S. Dey, J.-P.P. Richard, Y. Li, and L.A. Miller, Extreme inequalities for infinite group problems.,
                Mathematical Programming 121 (2010) 145–170.

        [61]: R.E. Gomory and E.L. Johnson, T-space and cutting planes, Mathematical Programming 96 (2003) 341–375.
    """
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

def dg_2_step_mir_limit(f=3/5, d=3):
    """
    Summary:
        - Name: DG-2-Step MIR Limit;
        - Infinite; Dim = 1; Slopes = 1; Discontinuous; Simple sets method;
        - Discovered [33] p.41, def.12;
        - Proven [33] p.43, lemma 14.

    Parameters:
        f (real) \in (0,1);
        d (positive interger): number of slopes on [0,f).

    Function is known to be extreme under the conditions:
        0 < f < 1;
        d >= ceil(1 / (1 - f)) - 1.

    Note:
        This is the limit function as alpha in two_step_mir()
        tends (from left) to f/d, where d is integer;
        c.f. [33] p.42, lemma 13.

        It's a special case of drlm_2_slope_limit(),
        dg_2_step_mir_limit(f, d) =
        multiplicative_homomorphism(drlm_2_slope_limit(f=1-f, nb_pieces_left=1, nb_pieces_right=d), -1)

    Examples:
        [33] p.42, Fig.6 ::
            sage: h = dg_2_step_mir_limit(f=3/5, d=3)

    Reference:
        [33]: S. Dash and O. G¨unl¨uk, Valid inequalities based on simple mixed-integer sets.,
                Proceedings 10th Conference on Integer Programming and Combinatorial Optimization
                (D. Bienstock and G. Nemhauser, eds.), Springer-Verlag, 2004, pp. 33–45.
    """
    if not (bool(0 < f < 1) & (d >= 1)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if not bool(d >= ceil(1 / (1 - f)) - 1):
        logging.info("Conditions for extremality are NOT satisfied.")
    else:
        logging.info("Conditions for extremality are satisfied.")
    f = nice_field_values([f])[0]
    field = f.parent()
    pieces = []
    for k in range(d):
        left_x = f * k / d
        right_x = f * (k + 1) / d
        pieces = pieces + \
                 [[singleton_interval(left_x), FastLinearFunction(field(0), left_x / f)], \
                  [open_interval(left_x, right_x), FastLinearFunction(1 / (f - 1), (k + 1)/(d + 1)/(1 - f))]]
    pieces.append([closed_interval(f, field(1)), FastLinearFunction(1 / (f - 1), 1 /(1 - f))])
    h = FastPiecewise(pieces, merge=False)
    # FIXME: Need to develop code for discontinuous functions.
    logging.warn("This is a discontinuous function; code for handling discontinuous functions is not implemented yet.")
    return h

def drlm_2_slope_limit(f=3/5, nb_pieces_left=3, nb_pieces_right=4):
    """
    Summary:
        - Name: drlm's 2-Slope Limit;
        - Infinite; Dim = 1; Slopes = 1; Discontinuous; Group relations method;
        - Discovered [40] p.158 def.10;
        - Proven [40] p.159 thm.8.

    Parameters:
        f (real) \in (0,1);
        nb_pieces_left (positive integer) : number of linear pieces to the left of f;
        nb_pieces_right (positive integer) : number of linear pieces to the right of f.

    Function is known to be extreme under the conditions:
        nb_pieces_left * (1-f) <= nb_pieces_right * f.

    Examples:
        [40] p.159 Fig.4 ::
            sage: h = drlm_2_slope_limit(f=3/5, nb_pieces_left=3, nb_pieces_right=4)

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
    # in irrational case, try to coerce to common number field
    [f, s, delta_1, delta_2, m, d] =  nice_field_values([f, s, delta_1, delta_2, m, d])
    pieces = []
    for k in range(m):
        pieces = pieces + \
                 [[singleton_interval(f * k / m), FastLinearFunction(0, k / m)], \
                  [open_interval(f * k / m, f * (k + 1) / m), FastLinearFunction(s, -k * delta_1)]]
    pieces.append([singleton_interval(f), FastLinearFunction(0, 1)])
    for k in range(d, 0, - 1):
        pieces = pieces + \
                 [[open_interval(1 - (1 - f)* k / d, 1 - (1 - f)*(k - 1)/d), FastLinearFunction(s, -s*f + 1 - (d - k + 1)*delta_2)], \
                  [singleton_interval(1 - (1 - f)*(k - 1)/d), FastLinearFunction(0, (k - 1) / d)]]
    psi = FastPiecewise(pieces, merge=False)    
    # FIXME: Need to develop code for discontinuous functions.
    logging.warn("This is a discontinuous function; code for handling discontinuous functions is not implemented yet.")
    return psi

def drlm_3_slope_limit(f=1/5):
    """
    Summary:
        - Name: drlm-3-Slope Limit;
        - Infinite; Dim = 1; Slopes = 2; Discontinuous; Group relations method;
        - Discovered [40] p.161 def.11;
        - Proven [40] p.161 thm.9.

    Parameters:
        f (real) \in (0,1);

    Function is known to be extreme under the conditions:
        0 < f < 1/3.

    Note:
        This is the limit function as bkpt tends to f in drlm_backward_3_slope(f, bkpt).

    Examples:
        [40] p.162 Fig.5 ::
            sage: h = drlm_3_slope_limit(f=1/5)

    Reference:
        [40]: S.S. Dey, J.-P.P. Richard, Y. Li, and L.A. Miller, Extreme inequalities for infinite group problems.,
                Mathematical Programming 121 (2010) 145–170.
    """
    if not bool(0 < f < 1):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if not bool(0 < f < 1/3):
        logging.info("Conditions for extremality are NOT satisfied.")
    else:
        logging.info("Conditions for extremality are satisfied.")
    f = nice_field_values([f])[0]
    field = f.parent()
    pieces = [[closed_interval(0, f), FastLinearFunction(1/f, 0, field=field)], \
              [open_interval(f, 1), FastLinearFunction(1/(f + 1), 0, field=field)], \
              [singleton_interval(field(1)), FastLinearFunction(field(0), 0, field=field)]]
    kappa = FastPiecewise(pieces, merge=False)
    # FIXME: Need to develop code for discontinuous functions.
    logging.warn("This is a discontinuous function; code for handling discontinuous functions is not implemented yet.")
    return kappa

def bccz_counterexample(f=2/3, q=4, eta=1/1000, maxiter=10000):
    """
    return function psi, a counterexample to GJ's conjecture;
    psi is a continuous facet (hence extreme), but is not piecewise linear. c.f.[IR1]

    Parameters:
        f (real) \in (0,1);
        q (real), q > 2 : ratio of the geometric series;
        eta (real), 0 <= eta < 1 : to control the serie sum;
        maxiter (integer) : maximum number of iterations;

    Note:
        psi is the uniform limit of the sequence of functions psi_n,
        generated by psi_n_in_bccz_counterexample_construction(f, [e[0], e[1], ..., e[n - 1]]).
        e is a geometric series with ratio q, such that:
        0 < ... < e[n] <= e[n - 1] <= ... <= e[1] <= e[0] <= 1 - f and \sum_{i = 0}^{\infty} {2^i * e[i]} <= f.
        The first n terms of e are generated by generate_example_e_for_psi_n(f, n, q, eta)
        
    See also:
        def generate_example_e_for_psi_n(f, n q, eta)
        def psi_n_in_bccz_counterexample_construction(f, e)

    Examples:
        quick evaluations::
            sage: bccz_counterexample(f=2/3, q=4, eta=0, maxiter=10000)(r=1/5)
            sage: bccz_counterexample(f=2/3, q=4, eta=0, maxiter=10000)(r=1/4)

        too many iterations::
            sage: bccz_counterexample(f=2/3, q=4, eta=0, maxiter=10000)(r=9/40)

    Reference:
        [IR1]:  A. Basu, M. Conforti, G. Cornuéjols, and G. Zambelli, A counterexample to a conjecture of Gomory and Johnson,
                    Mathematical Programming Ser. A 133 (2012), 25–38.
    """
    if not (bool(0 < f < 1) & bool(q > 2) & bool(0 <= eta < 1)):
        raise ValueError, "Bad parameters."
    def evaluate_psi_at_r(r):
        if r == 0:
            return 0
        if bool (f <= r <= 1):
            return (1 - r) / (1 - f)
        if bool (0 < r < f):
            z = (1 - eta)*(q - 2) / q * min(f, 1 - f)
            # or take z = min((1 - eta)*(q - 2)*f / q , 1 - f)
            n = 0
            x_left = 0
            x_right = f
            y_left = 0
            y_right = 1
            while not (bool((x_left + x_right - z)/2 <= r <= (x_left + x_right + z)/2)) and (n < maxiter):
                if bool(r < (x_left + x_right - z)/2):
                    x_right = (x_left + x_right - z)/2 
                    y_right = (y_left + y_right + z/(1 - f))/2
                else:
                    x_left = (x_left + x_right + z)/2
                    y_left = (y_left + y_right - z/(1 - f))/2
                z = z / q
                n += 1
            if n == maxiter:
                logging.warn("Reaching max number of iterations, return approximate psi(%s)" %r)
            return (y_left + y_right)/2 - (r - (x_left + x_right)/2) / (1 - f)
        else:
            raise ValueError, "outside domain"
    logging.warn("This function is not piecewise linear; code for handling this function is not implemented.")
    return evaluate_psi_at_r


def generate_example_e_for_psi_n(f=2/3, n=7, q=4, eta=1/1000):
    """
    return the first n terms of a geometric series e that satisfies 
    0 < ... < e[n] <= e[n - 1] <= ... <= e[1] <= e[0] <= 1 - f and \sum_{i = 0}^{\infty} {2^i * e[i]} <= f.
    This can be used in psi_n_in_bccz_counterexample_construction(f, [e[0],...,e[n-1]]), so that the function constructed is extreme.

    Parameters:
        f (real) \in (0,1);
        n (integer);
        q (real), q > 2 : ratio of the geometric series;
        eta (real), 0 <= eta < 1 : to control the serie sum, \sum_{i = 0}^{\infty} {2^i * e[i]} <= (1 - eta)*f.

    Note:
        If (eta == 0) and (f >= 1/2), then \sum_{i = 0}^{\infty} {2^i * e[i]} = f.
        This case is not mentioned in [IR1], but using a similar (or same) proof to [IR1], one can show that:
        1) psi_n still converges uniformly to psi;
        2) The limit funciton psi is a continuous facet (hence extreme);
        3) psi is not piecewise linear. 
        Also notice that:
        4) psi is not in W^{1,1}.

    Reference: 
        [IR1]:  A. Basu, M. Conforti, G. Cornuéjols, and G. Zambelli, A counterexample to a conjecture of Gomory and Johnson,
                    Mathematical Programming Ser. A 133 (2012), 25–38.
    """
    if n == 0:
        return []
    if not (bool(0 < f < 1) & bool(q > 2) & bool(0 <= eta < 1)):
        raise ValueError, "Bad parameters."
    x = (1 - eta)*(q - 2) / q * min(f, 1 - f)
    # or take x = min((1 - eta)*(q - 2)*f / q , 1 - f) 
    e = [x / q^i for i in range(n)]
    return e


def psi_n_in_bccz_counterexample_construction(f=2/3, e=[1/12, 1/24]):
    """
    Summary: 
        - Name: psi_n in the construction of BCCZ's counterexample to GJ's conjecture;
        - Infinite; Dim = 1; Slopes = 2; Continuous;  Analysis of subadditive polytope method.
        - Discovered [IR1]  p.30, section.3, fig.1;
        - Proven [IR1] p.35, thm.4.7.

    Note:
        The (uniform) limit \psi = \lim_{n to \infty} \psi_n is well defined if \sum_{i = 0}^{\infty} {2^i * e[i]} < f.
        The (unifrom) limit \psi is a continuous facet, but is not piecewise linear. A counterexample of GJ's conjecture.
        Could use the function generate_example_e_for_psi_n(f, n, q) to generate a sequence e that satisfies the conditions for extremality.

    Parameters:
        f (real) \in (0,1);
        e (list of reals, with length = n) \in (0,f).

    Function is known to be extreme under the conditions:
        0 < f < 1;
        0 < e[n - 1] <= e[n - 2] <= ... <= e[1] <= e[0] <= 1 - f;
        \sum_{i = 0}^{n - 1} 2^i * e[i] < f.

    Examples:
        [IR1]  p.30, fig.1::
            sage: h = psi_n_in_bccz_counterexample_construction(f=2/3, e=[1/12, 1/24])

            sage: h = psi_n_in_bccz_counterexample_construction(f=4/5, e=[1/5, 1/20, 1/80, 1/320, 1/1280])

            sage: h = psi_n_in_bccz_counterexample_construction(f=4/5, e=generate_example_e_for_psi_n(f=4/5, n=7, q=3, eta=0))

            sage: sum([plot(psi_n_in_bccz_counterexample_construction(e=generate_example_e_for_psi_n(n=n)), color=color, legend_label="psi_%d"%n) for n, color in zip(range(7),rainbow(7))])

    Reference: 
        [IR1]:  A. Basu, M. Conforti, G. Cornuéjols, and G. Zambelli, A counterexample to a conjecture of Gomory and Johnson,
                    Mathematical Programming Ser. A 133 (2012), 25–38.
    """
    if not bool(0 < f < 1):
        raise ValueError, "Bad parameters. Unable to construct the function."
    n = len(e)
    if n == 0:
        logging.info("Conditions for extremality are satisfied.")
        return piecewise_function_from_breakpoints_and_values([0,f,1], [0,1,0])
    a = [1]
    b = [f]
    sum_e = 0
    if not bool(0 < e[0]):
        raise ValueError, "Bad parameters. Unable to construct the function."
    t = bool(e[0] <= 1 - f)
    for i in range(0, n):
        a.append((b[i] + e[i]) / 2)
        b.append((b[i] - e[i]) / 2)
        sum_e = sum_e + (2^i) * e[i]
        if not (bool(e[i] > 0) & bool(sum_e < f)):
            raise ValueError, "Bad parameters. Unable to construct the function."
        if not (i == 0) | bool(e[i] <= e[i-1]):
            t = False
    if t:
        logging.info("Conditions for extremality are satisfied.")
    else:
        logging.info("Conditions for extremality are NOT satisfied.")
    interval_lengths =  interval_length_n_step_mir(n + 1, 1, a, b)
    nb_interval = len(interval_lengths)
    interval_length_positive = sum(interval_lengths[i] for i in range(0, nb_interval, 2))
    interval_length_negative = sum(interval_lengths[i] for i in range(1, nb_interval, 2))
    s_negative = a[0] /(b[0] - a[0])
    s_positive = - s_negative * interval_length_negative / interval_length_positive 
    slopes = [s_positive, s_negative] * (nb_interval // 2)
    return piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes)

def bhk_irrational(f=4/5, d1=3/5, d2=1/10, a0=15/100, delta=(1/200, sqrt(2)/200)):
    """
    Summary:
        - Name: BHK's irrational function.
        - Infinite; Dim = 1; Slopes = 3; Continuous;  Covered intervals and equivariant perturbation.
        - Discovered [IR2]  p.33, section.5.2, fig.9-10.
        - Proven [IR2] p.34, thm.5.3.

    Parameters:
        f (real) \in (0,1);
        d1 (real): length of the positive slope;
        d2 (real): length of the zero slopes;
        a0 (real): length of the first zig-zag;
        delta (n-tuple of reals): length of the extra zig-zags.

    Function is known to be extreme under the conditions:
        0 < f < 1;
        d1, d2, a0, delta > 0;
        d1 + d2 < f;
        len(delta) == 2
        sum(delta) <  d2 / 4; Weaker condition: 2*delta[0] + delta[1] < d2 / 2;
        the two components of delta are linearly independent over \Q.

    Relation between the code parameters and the paper parameters:
        t1 = delta[0], t2 = delta[0] + delta[1], ...
        a1 = a0 + t1, a2 = a0 + t2, ...
        A = f/2 - a0/2 - d2/4,
        A0 = f/2 - a0/2 + d2/4.

    Examples:
        [IR2]  p.34, thm.5.3::
            sage: h = bhk_irrational(f=4/5, d1=3/5, d2=1/10, a0=15/100, delta=(1/200, sqrt(2)/200))
            sage: h = bhk_irrational(f=4/5, d1=3/5, d2=1/10, a0=15/100, delta=(1/200, 6* sqrt(2)/200, 1/500))

    Reference:
        [IR2] A. Basu, R. Hildebrand, and M. Köppe, Equivariant perturbation in Gomory and Johnson’s infinite group problem.
                I. The one-dimensional case, Mathematics of Operations Research (2014), doi:10. 1287/moor.2014.0660
    """
    if not (bool(0 < f < 1) and bool(d1 > 0) and bool(d2 > 0) and bool(a0 > 0) and (len(delta) >= 2) \
            and bool(min(delta) > 0) and bool(d1 + d2 < f) and (sum(delta) < f/2 - d2/4 - 3*a0/2) ):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if len(delta) == 2:
        if is_QQ_linearly_independent(delta) and  2*delta[0] + delta[1] < d2 / 2:
            logging.info("Conditions for extremality are satisfied.")
        else:
            logging.info("Conditions for extremality are NOT satisfied.")
    else:
        logging.info("len(delta) >= 3, Conditions for extremality are unknown.")
    d3 = f - d1 - d2
    c2 = 0
    c3 = -1/(1-f)
    c1 = (1-d2*c2-d3*c3)/d1
    d21 = d2 / 2
    d31 = c1 / (c1 - c3) * d21
    d11 = a0 - d31
    d13 = a0 - d21
    d12 = (d1 - d13)/2 - d11
    d32 = d3/2 - d31
    zigzag_lengths = []
    zigzag_slopes = []
    delta_positive = 0
    delta_negative = 0
    for delta_i in delta:
        delta_i_negative = c1 * delta_i / (c1 - c3)
        delta_i_positive = delta_i - delta_i_negative
        delta_positive += delta_i_positive
        delta_negative += delta_i_negative
        zigzag_lengths = zigzag_lengths + [delta_i_positive, delta_i_negative]
        zigzag_slopes = zigzag_slopes + [c1, c3]
    d12new = d12 - delta_positive
    d32new = d32 - delta_negative

    slopes_left = [c1,c3] + zigzag_slopes + [c1,c3,c2]
    slopes = slopes_left + [c1] + slopes_left[::-1] + [c3]
    intervals_left = [d11,d31] + zigzag_lengths + [d12new,d32new,d21]
    interval_lengths = intervals_left + [d13] + intervals_left[::-1] + [1-f]
    return piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes)

def bhk_gmi_irrational(f=4/5, d1=3/5, d2=1/10, a0=15/100, delta=(1/200, sqrt(2)/200), alpha=95/100):
    """
    A version of the irrational function with non-zero second slope,
    obtained by forming a convex combination of a modified version of the irrational function with the GMI cut.
    """
    if not (bool(0 < f < 1) and bool(d1 > 0) and bool(d2 > 0) and bool(a0 > 0) and (len(delta) >= 2) \
            and bool(min(delta) > 0) and bool(d1 + d2 < f) and (sum(delta) < f/2 - d2/4 - 3*a0/2) \
            and bool(0 < alpha < 1)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    # FIXME: Extremality condition ?
    d3 = f - d1 - d2
    c2 = 0
    c3 = (-1/(1-f) - (1 - alpha)/f) / alpha
    c1 = (1-d2*c2-d3*c3)/d1
    d21 = d2 / 2
    d31 = c1 / (c1 - c3) * d21
    d11 = a0 - d31
    d13 = a0 - d21
    d12 = (d1 - d13)/2 - d11
    d32 = d3/2 - d31
    zigzag_lengths = []
    zigzag_slopes = []
    delta_positive = 0
    delta_negative = 0
    for delta_i in delta:
        delta_i_negative = c1 * delta_i / (c1 - c3)
        delta_i_positive = delta_i - delta_i_negative
        delta_positive += delta_i_positive
        delta_negative += delta_i_negative
        zigzag_lengths = zigzag_lengths + [delta_i_positive, delta_i_negative]
        zigzag_slopes = zigzag_slopes + [c1, c3]
    d12new = d12 - delta_positive
    d32new = d32 - delta_negative
    slopes_left = [c1,c3] + zigzag_slopes + [c1,c3,c2]
    slopes = slopes_left + [c1] + slopes_left[::-1] + [-1/(1 - f)]
    intervals_left = [d11,d31] + zigzag_lengths + [d12new,d32new,d21]
    interval_lengths = intervals_left + [d13] + intervals_left[::-1] + [1-f]

    bhk = piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes)
    gmi = gmic(f)
    return alpha * bhk + (1 - alpha) * gmi

def chen_3_slope(f=1/2, lam=8):
    """
    A continuous 3-slope extreme function that has non-degenerate intervals with a zero derivative.

    Parameters:
        f (real) \in (0,1);
        lam (real): the first slope has length 1/lam.

    Requirement:
        0 < f <= 1/2;
        lam > 3*max(1/f, 1 /(1-f))

    Examples:
        [KChen_thesis]  p.33, fig.7 ::
            sage: h = chen_3_slope(f=1/2, lam=8)

    Reference:
        [KChen_thesis]:  K. Chen, Topics in group methods for integer programming,
                            Ph.D. thesis, Georgia Institute of Technology, June 2011.
    """
    # FIXME:
    # (1) Is this example extreme?? extremality_test(h) returns False.
    # (2) Is f <= 1/2 necessary?  used in "All intervals are covered" ?
    # (3) generate_maximal_additive_faces() has a bug. try: chen_3_slope(f=2/3, lam=10)
    if not (bool(0 < f < 1) and (lam > 3*max(1/f, 1 /(1-f)))):
        raise ValueError, "Bad parameters. Unable to construct the function."
    alpha = f / 2 - 3 / (2*lam)
    beta = 1/2 - f/2 - 3 / (2*lam)
    bkpts = [0, 1/lam, 1/lam + alpha, 2/lam + alpha, 2/lam + 2*alpha, f, \
                4/lam + 2*alpha, 4/lam + 2*alpha + beta, 5/lam + 2*alpha + beta, 5/lam + 2*alpha + 2*beta, 1]
    values = [0, 2/3, 2/3, 1/3, 1/3, 1, 2/3, 2/3, 1/3, 1/3, 0]
    return  piecewise_function_from_breakpoints_and_values(bkpts, values)

def chen_4_slope(f=7/10, s_pos=2, s_neg=-4, lam1=1/4, lam2=1/4):
    """
    This 4-slope function is shown [KChen_thesis] to be extreme.
    This function can also be shown to be a facet.

    Parameters:
        f (real) \in (0,1);
        s_pos, s_neg (real): positive slope and negative slope
        lam1, lam2 (real).

    Requirement:
        1/2 <= f <= 1;
        s_pos > 1/f;
        s_neg < 1/(f - 1);
        0 < lam1 < min(1/2, (s_pos - s_neg) / s_pos / (1 - s_neg * f));
        f - 1 / s_pos < lam2 < min(1/2, (s_pos - s_neg) / s_neg / (s_pos * (f - 1) - 1)).

    Examples:
        [KChen_thesis]  p.38, fig.8 ::

        sage: h = chen_4_slope(f=7/10, s_pos=2, s_neg=-4, lam1=1/4, lam2=1/4)

        sage: h = chen_4_slope(f=1/2, s_pos=4, s_neg=-4, lam1=1/3, lam2=1/3)

        The following function's parameters do not satisfy the requirement, however it is extreme ::

        sage: h = chen_4_slope(f=1/2, s_pos=5, s_neg=-5, lam1=1/5, lam2=1/5)

    Reference:
        [KChen_thesis]:  K. Chen, Topics in group methods for integer programming,
                            Ph.D. thesis, Georgia Institute of Technology, June 2011.
    """
    if not (bool(0 < f < 1) and bool(s_pos > 1/f) and bool(s_neg < 1/(f - 1)) \
                            and bool(0 < lam1 < 1) and bool(0 < lam2 < 1)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if bool(1/2 <= f) and bool(lam1 < (s_pos - s_neg) / s_pos / (1 - s_neg * f)) and \
        bool (f - 1 / s_pos < lam2 < (s_pos - s_neg) / s_neg / (s_pos * (f - 1) - 1)):
        logging.info("Conditions for extremality are satisfied.")
    else:
        logging.info("Conditions for extremality are NOT satisfied.")
    slopes = [s_pos, s_neg, 1/f, s_neg, s_pos, s_neg, s_pos, 1/(f-1), s_pos, s_neg]
    aa = lam1 * (1 - s_neg * f) / 2 / (s_pos - s_neg)
    a = lam1 * f / 2
    b = f - a
    bb = f - aa
    c = 1 + lam2 * (f - 1) / 2
    d = 1 + f - c
    cc = 1 + (s_pos * lam2 * (f - 1) - lam2) / 2 / (s_pos - s_neg)
    dd = 1 + f - cc
    return piecewise_function_from_breakpoints_and_slopes([0, aa, a, b, bb, f, dd, d, c, cc, 1], slopes)
