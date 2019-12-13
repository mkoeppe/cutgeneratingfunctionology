from six.moves import range

from .fast_piecewise import PiecewiseLinearFunction_1d
from .parametric_family import ParametricFamily

class ParametricFamily_gmic(ParametricFamily):

    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = gmic()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - Name: GMIC (Gomory mixed integer cut);
        - Infinite (or Finite); Dim = 1; Slopes = 2; Continuous; Analysis of subadditive polytope method;
        - Discovered [55] p.7-8, Eq.8;
        - Proven extreme (for infinite group) [60] p.377, thm.3.3; (finite group) [57] p.514, Appendix 3.
        - (Although only extremality has been established in literature, the same proof shows that) gmic is a facet.

    Parameters:
        f (real) `\in (0,1)`.

    Examples: 
        [61] p.343, Fig. 1, Example 1 ::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: h = gmic(4/5)
            sage: extremality_test(h, False)
            True

    Reference: 
        [55]: R.E. Gomory, An algorithm for the mixed integer problem, Tech. Report RM-2597, RAND Corporation, 1960.

        [57]: R.E. Gomory, Some polyhedra related to combinatorial problems, Linear Algebra and its Application 2 (1969) 451-558.

        [60]: R.E. Gomory and E.L. Johnson, Some continuous functions related to corner polyhedra, part II, Mathematical Programming 3 (1972) 359-389.

        [61]: R.E. Gomory and E.L. Johnson, T-space and cutting planes, Mathematical Programming 96 (2003) 341-375.
    """

    def _construct_function(self, f=4/5, field=None, conditioncheck=True):
        if not bool(0 < f < 1):
            raise ValueError("Bad parameters. Unable to construct the function.")
        claimed_parameter_attribute = None
        if conditioncheck:
            claimed_parameter_attribute = 'extreme'
        gmi_bkpt = [0,f,1]
        gmi_values = [0,1,0]
        h = piecewise_function_from_breakpoints_and_values(gmi_bkpt, gmi_values, field=field)
        h._claimed_parameter_attribute = claimed_parameter_attribute
        return h

gmic = ParametricFamily_gmic()

class ParametricFamily_gj_2_slope(ParametricFamily):

    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = gj_2_slope()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - Name: Gomory--Johnson's 2-Slope;
        - Infinite (or Finite); Dim = 1; Slopes = 2; Continuous; Analysis of subadditive polytope method;
        - Discovered [61] p.352, Fig.5, construction 1;
        - Proven extreme (infinite group) [60] p.377, thm.3.3; [61] p.352, thm.4; p.354, thm.5.
        - gj_2_slope is a facet.

    Parameters:
        * f (real) `\in (0,1)`;      
        * lambda_1 (real) in (0,1].

    Function is known to be extreme under the conditions:
        * 0 < lambda_1 <=1, 
        * lambda_1 < f/(1 - f).

    Examples: [61] p.354, Fig.6 ::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = gj_2_slope(f=3/5, lambda_1=1/6)
        sage: extremality_test(h, False)
        True
        sage: h = gj_2_slope(f=3/5, lambda_1=1/2)
        sage: extremality_test(h, False)
        True
        sage: h = gj_2_slope(f=3/5, lambda_1=1)
        sage: extremality_test(h, False, f=3/5)         # Provide f to suppress warning
        True

    Reference:
        [60]: R.E. Gomory and E.L. Johnson, Some continuous functions related to corner polyhedra, part II, Mathematical Programming 3 (1972) 359-389.

        [61]: R.E. Gomory and E.L. Johnson, T-space and cutting planes, Mathematical Programming 96 (2003) 341-375.
    """

    def _construct_function(self, f=3/5, lambda_1=1/6, field=None, conditioncheck=True):
        if not (bool(0 < f < 1) & bool(0 < lambda_1 < f/(1 - f))):
            raise ValueError("Bad parameters. Unable to construct the function.")
        claimed_parameter_attribute = None
        if conditioncheck:
            if not (bool(lambda_1 <= 1)):
                logging.info("Conditions for extremality are NOT satisfied.")
                claimed_parameter_attribute = 'constructible'
            else:
                logging.info("Conditions for extremality are satisfied.")
                claimed_parameter_attribute = 'extreme'
        bkpts = [0, (f - lambda_1*(1 - f))/2, (f + lambda_1*(1 - f))/2, f, 1]
        values = [0, (1 + lambda_1)/2, (1 - lambda_1)/2, 1, 0]
        h = piecewise_function_from_breakpoints_and_values(bkpts, values, field=field)
        h._claimed_parameter_attribute = claimed_parameter_attribute
        return h

gj_2_slope = ParametricFamily_gj_2_slope()

class ParametricFamily_gj_2_slope_repeat(ParametricFamily):

    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = gj_2_slope_repeat()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - Name: Gomory--Johnson's 2-Slope-repeat;
        - Infinite (or Finite); Dim = 1; Slopes = 2; Continuous; Analysis of subadditive polytope method;
        - Discovered [61] p.354, Fig.7, construction 2;
        - Proven extreme (for infinite group) [60] p.377, thm.3.3; [61] p.354, thm.5; p.355, thm.6.
        - gj_2_slope_repeat is a facet.

    Parameters:
        * f (real) `\in (0,1)`;       
        * s_positive, s_negative (real);
        * m, n >= 2 (integer).

    Function is known to be extreme under the conditions:
        * `0 < f < 1`;
        * s_positive > 1/f;  s_negative < 1/(f - 1);
        * m >= (s_positive - s_positive*s_negative*f) / (s_positive - s_negative);
        * n >= (- s_negative + s_positive*s_negative*(f - 1)) / (s_positive - s_negative).

    Examples: [61] p.354, Fig.7 ::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = gj_2_slope_repeat(f=3/5, s_positive=4, s_negative=-5, m=4, n=3)
        sage: extremality_test(h, False)
        True

    Reference:
        [60]: R.E. Gomory and E.L. Johnson, Some continuous functions related to corner polyhedra, part II, Mathematical Programming 3 (1972) 359-389.

        [61]: R.E. Gomory and E.L. Johnson, T-space and cutting planes, Mathematical Programming 96 (2003) 341-375.
    """

    def _construct_function(self, f=3/5, s_positive=4, s_negative=-5, m=4, n=3, field=None, conditioncheck=True):
        if not (bool(0 < f < 1) & (m >= 2) & (n >= 2) & bool (s_positive > 1 / f) & bool(s_negative < 1/(f - 1))):
            raise ValueError("Bad parameters. Unable to construct the function.")
        claimed_parameter_attribute = None
        if conditioncheck:
            if not (bool(m >= (s_positive - s_positive*s_negative*f) / (s_positive - s_negative)) & bool(n >= (- s_negative + s_positive*s_negative*(f - 1)) / (s_positive - s_negative))):
                logging.info("Conditions for extremality are NOT satisfied.")
                claimed_parameter_attribute = 'constructible'
            else:
                logging.info("Conditions for extremality are satisfied.")
                claimed_parameter_attribute = 'extreme'
        len1_positive = (1 - s_negative*f) / (s_positive - s_negative) / m
        len1_negative = (f - m*len1_positive) / (m - 1)
        len2_negative = (1 - s_positive*(f - 1)) / (s_positive - s_negative) / n
        len2_positive = (1 - f - n*len2_negative) / (n - 1)
        interval_lengths = [len1_positive, len1_negative] * (m - 1) + [len1_positive, len2_negative] + [len2_positive, len2_negative]*(n - 1)
        slopes = [s_positive, s_negative]*(m + n - 1)
        h = piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes, field=field)
        h._claimed_parameter_attribute = claimed_parameter_attribute
        return h

gj_2_slope_repeat = ParametricFamily_gj_2_slope_repeat()

class ParametricFamily_dg_2_step_mir(ParametricFamily):

    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = dg_2_step_mir()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - Name: Dash-Gunluk's 2-Step MIR;
        - Infinite (or Finite); Dim = 1; Slopes = 2; Continuous; Simple sets method;
        - Discovered [33]  p.39 def.8, Fig.5;
        - Proven extreme (for infinite group) [60] p.377, thm.3.3.
        - dg_2_step_mir is a facet.

    Parameters:
        * f (real) in (0,1);
        * alpha (real) in (0,f).

    Function is known to be extreme under the conditions:
        * 0 < alpha < f < 1;
        * f / alpha < ceil(f / alpha) <= 1 / alpha.

    Examples: [33] p.40, Fig.5 ::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = dg_2_step_mir(f=4/5, alpha=3/10)
        sage: extremality_test(h, False)
        True

    Reference:
        [33]: S. Dash and O. Gunluk, Valid inequalities based on simple mixed-integer sets.,
            Proceedings 10th Conference on Integer Programming and Combinatorial Optimization
            (D. Bienstock and G. Nemhauser, eds.), Springer-Verlag, 2004, pp. 33-45.

        [60]: R.E. Gomory and E.L. Johnson, Some continuous functions related to corner polyhedra, part II, Mathematical Programming 3 (1972) 359-389.
    """

    def _claimed_parameter_attribute(self, f, alpha, **kwargs):
        if not (bool(0 < alpha < f < 1) & bool(f / alpha < ceil(f / alpha))):
            return 'not_constructible'
        if not bool(ceil(f / alpha) <= 1 / alpha):
            return 'constructible'
        else:
            return 'extreme'

    def _construct_function(self, f=4/5, alpha=3/10, field=None, conditioncheck=True):
        if not (bool(0 < alpha < f < 1) & bool(f / alpha < ceil(f / alpha))):
            raise ValueError("Bad parameters. Unable to construct the function.")
        rho = f - alpha * floor(f / alpha)
        tau = ceil(f / alpha)
        s_positive = (1 - rho*tau) / (rho*tau*(1 - f))
        s_negative = - 1/(1 - f)
        interval_lengths = [rho, alpha - rho] * tau
        interval_lengths[-1] = 1 - f
        slopes = [s_positive, s_negative] * tau
        return piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes, field=field)

dg_2_step_mir = ParametricFamily_dg_2_step_mir()

def interval_length_n_step_mir(n, m, a, b):
    if m == n:
        return [b[n - 1], a[n - 1] - b[n - 1]]
    else:
        l = interval_length_n_step_mir(n, m + 1, a, b)
        result = l * ceil(b[m - 1] / a[m]) 
        result[-1] = a[m - 1] - b[m - 1]
        return result

class ParametricFamily_kf_n_step_mir (ParametricFamily):

    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kf_n_step_mir()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - Name: Kianfar-Fathi's n-Step MIR;
        - Infinite (or Finite); Dim = 1; Slopes = 2; Continuous; Simple sets method;
        - Discovered [74]  p.328, def.3, thm.2;
        - Proven extreme (for infinite group) [60] p.377, thm.3.3.
        - (Although only extremality has been established in literature, the same proof shows that) ``kf_n_step_mir`` is a facet.

    Parameters:
        * f (real) `\in (0,1)`;
        * a (list of reals, with length = n) `\in (0,f)`.

    Function is known to be extreme under the conditions:
        * 0 < a[1] < f < 1 == a[0];
        * a[i] > 0, for i = 0, 1, ... , n-1;
        * b[i - 1] / a[i] < ceil(b[i - 1] / a[i]) <= a[i - 1] / a[i],  for i = 1, 2, ... , n-1;
    where,
        * b[0] = f;
        * b[i] = b[i - 1] - a[i] * floor(b[i - 1] / a[i]),  for i = 1, 2, ... , n-1.

    Note:
        if a[i] > b[i-1] for some i, then the kf_n_step_mir function degenerates, i.e.
        kf_n_step_mir(f, [a[0], .. , a[n - 1]]) = kf_n_step_mir(f, [a[0], .. a[i - 1], a[i + 1], ... , a[n - 1]])

    Examples: [74] p.333 - p.335, Fig.1 - Fig.6 ::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = kf_n_step_mir(f=4/5, a=[1])
        sage: extremality_test(h, False)
        True
        sage: h = kf_n_step_mir(f=4/5, a=[1, 3/10])
        sage: extremality_test(h, False)
        True
        sage: h = kf_n_step_mir(f=4/5, a=[1, 3/10, 8/100])
        sage: extremality_test(h, False)
        True
        sage: h = kf_n_step_mir(f=4/5, a=[1, 3/10, 8/100, 3/100])
        sage: extremality_test(h, False)
        True
        sage: h = kf_n_step_mir(f=4/5, a=[1, 45/100, 2/10, 558/10000, 11/1000])
        sage: extremality_test(h, False)
        True
        sage: h = kf_n_step_mir(f=4/5, a=[1, 48/100, 19/100, 8/100, 32/1000, 12/1000])
        sage: extremality_test(h, False)
        True

    Reference:
        [60]: R.E. Gomory and E.L. Johnson, Some continuous functions related to corner polyhedra, part II, Mathematical Programming 3 (1972) 359-389.

        [74]: K. Kianfar and Y. Fathi, Generalized mixed integer rounding valid inequalities:
                Facets for infinite group polyhedra, Mathematical Programming 120 (2009) 313-346.
    """

    def _claimed_parameter_attribute(self, f, a, **kwargs):
        if (not a) | (not bool(0 < f < 1 == a[0])):
            return 'not_constructible'
        b = []
        b.append(f)
        n = len(a)
        for i in range(1, n):
            b.append(b[i - 1] - a[i] * floor(b[i - 1] / a[i]))
            if not (bool(0 < a[i]) & bool(b[i - 1] / a[i] < ceil(b[i - 1] / a[i]))):
                return 'not_constructible'
        for i in range(1, n):
            if not bool(ceil(b[i - 1] / a[i]) <= a[i - 1] / a[i]):
                return 'constructible'
        return 'extreme'

    def _construct_function(self, f=4/5, a=(1, 3/10, 8/100), field=None, conditioncheck=True):
        b = []
        b.append(f)
        n = len(a)
        for i in range(1, n):
            b.append(b[i - 1] - a[i] * floor(b[i - 1] / a[i]))
        interval_lengths =  interval_length_n_step_mir(n, 1, a, b)
        nb_interval = len(interval_lengths)
        interval_length_positive = sum(interval_lengths[i] for i in range(0, nb_interval, 2))
        interval_length_negative = sum(interval_lengths[i] for i in range(1, nb_interval, 2))
        s_negative = a[0] /(b[0] - a[0])
        s_positive = - s_negative * interval_length_negative / interval_length_positive
        slopes = [s_positive, s_negative] * (nb_interval // 2)
        return piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes, field=field)

kf_n_step_mir = ParametricFamily_kf_n_step_mir()


class ParametricFamily_gj_forward_3_slope(ParametricFamily):

    # FIXME: What is the relation between the parameters shown in the figure (taken from param graphics) and the construction parameters?
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = gj_forward_3_slope()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), ticks=[h.end_points(),[]], tick_formatter=[["$0$","$a'$","$a$","$b$","$b'$","$f$","$1$"], []], thickness=2)
        slope_formatter = ["$s^+$", "$s^-$", "$\\frac{1}{f}$", "$s^-$", "$s^+$", "$s^-$"]
        for i in range(len(h.end_points())-1):
            x = (h.end_points()[i] + h.end_points()[i+1])/2 -1/50
            y = h(x) +1/10
            g += text(slope_formatter[i], (x, y), axis_coords=False, vertical_alignment='bottom',horizontal_alignment='left',color='black')
        sphinx_plot(g)

    Summary:
        - Name: Gomory--Johnson' Forward 3-Slope;
        - Infinite (or Finite); Dim = 1; Slopes = 3; Continuous; Analysis of subadditive polytope method;
        - Discovered [61] p.359, Construction.3, Fig.8;
        - Proven extreme [61] p.359, thm.8.
        - gj_forward_3_slope is a facet.

    Parameters:
        * f (real) `\in (0,1)`;
        * lambda_1, lambda_2 (real) `\in (0,1)`.

    Function is known to be extreme under the conditions:
        * 0 <= lambda_1 <= 1/2;
        * 0 <= lambda_2 <= 1  (in literature).

    Note: 
        Since the domain and range are in [0,1], I think the conditions for a three-slope extreme function should be:

        (0 <= lambda_1 <= 1/2)  &  (0 <= lambda_2 <= 1) & (0 < lambda_1 * f + lambda_2 * (f - 1) < lambda_1 * f).

    Examples:
        [61] p.360, Fig.8 ::

            sage: from cutgeneratingfunctionology.igp import *
            sage: h = gj_forward_3_slope(f=4/5, lambda_1=4/9, lambda_2=1/3)
            sage: extremality_test(h, False)
            True
            sage: h = gj_forward_3_slope(f=4/5, lambda_1=4/9, lambda_2=2/3)
            sage: extremality_test(h, False)
            True
            sage: h = gj_forward_3_slope(f=4/5, lambda_1=4/9, lambda_2=1)
            sage: extremality_test(h, False)
            True

        Try irrational case ::

            sage: h = gj_forward_3_slope(f=sqrt(17)/5, lambda_1=2*sqrt(5)/9, lambda_2=2/sqrt(10))
            sage: extremality_test(h, False)
            True

    Reference:
        [61]: R.E. Gomory and E.L. Johnson, T-space and cutting planes, Mathematical Programming 96 (2003) 341-375.
    """

    def _construct_function(self, f=4/5, lambda_1=4/9, lambda_2=2/3, field=None, conditioncheck=True):
        if not bool(0 < f < 1):
            raise ValueError("Bad parameters. Unable to construct the function.")
        a = lambda_1 * f / 2
        a1 = a + lambda_2 * (f - 1) / 2
        if not bool(0 < a1 < a < f / 2):
            raise ValueError("Bad parameters. Unable to construct the function.")
        claimed_parameter_attribute = None
        if conditioncheck:
            # note the discrepancy with the published literature
            if not (bool(0 <= lambda_1 <= 1/2) & bool(0 <= lambda_2 <= 1)):
                logging.info("Conditions for extremality are NOT satisfied.")
                claimed_parameter_attribute = 'constructible'
            else:
                logging.info("Conditions for extremality are satisfied.")
                claimed_parameter_attribute = 'extreme'
        bkpts = [0, a1, a, f - a, f - a1, f, 1]
        values = [0, (lambda_1 + lambda_2)/2, lambda_1 / 2, 1 - lambda_1 / 2, 1 - (lambda_1 + lambda_2)/2, 1, 0]
        h = piecewise_function_from_breakpoints_and_values(bkpts, values, field=field)
        h._claimed_parameter_attribute = claimed_parameter_attribute
        return h

gj_forward_3_slope = ParametricFamily_gj_forward_3_slope()

class ParametricFamily_drlm_backward_3_slope(ParametricFamily):

    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = drlm_backward_3_slope()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - Name: Dey--Richard--Li--Miller's Backward 3-Slope;
        - Infinite; Dim = 1; Slopes = 3; Continuous; Group relations method;
        - Discovered [40] p.154 eq.5;
        - Proven [40] p.153 thm.6.
        - (Although only extremality has been established in literature, the same proof shows that) drlm_backward_3_slope is a facet.

    Parameters:
        f, bkpt (real) `\in (0,1)`.

    Function is known to be extreme under the conditions:
        f < bkpt < (1+f)/4 < 1.

    Note:
        In [40], they require that f, bkpt are rational numbers.
        The proof is based on interpolation of finite cyclic group extreme functions(cf. [8]), so it needs rational numbers.
        But in fact, by analysing covered intervals and using the condition f < bkpt <= (1+f)/4 < 1,
        one can prove that the function is extreme without assuming f, bkpt being rational numbers.
        
        In [61] p.374, Appendix C, p.360. Fig.10, they consider real number f, bkpt, and claim (without proof) that:

        (1) the function (named pi3(u)) is facet (thus extreme);
        (2) can add a perturbation (zigzag) on the third slope as shown in Fig.10;

        An extremality proof for the general (not necessarily rational) case appears in [KZh2015b, section 4].

    Examples:
        - Finite group --> Example 3.8 in [8] p.386,
        - Infinite group --> Interpolation using Equation 5 from [40] p.154 ::

            sage: from cutgeneratingfunctionology.igp import *
            sage: h = drlm_backward_3_slope(f=1/12, bkpt=2/12)
            sage: extremality_test(h, False)
            True
            sage: h = drlm_backward_3_slope(f=1/12, bkpt=3/12)
            sage: extremality_test(h, False)
            True

    References:

    - [8] J. Araoz, L. Evans, R.E. Gomory, and E.L. Johnson, Cyclic groups and knapsack facets,
      Mathematical Programming 96 (2003) 377-408.
                
    - [40] S.S. Dey, J.-P.P. Richard, Y. Li, and L.A. Miller, On the extreme inequalities of infinite group problems,
      Mathematical Programming 121 (2010) 145-170.

    - [61] R.E. Gomory and E.L. Johnson, T-space and cutting planes, Mathematical Programming 96 (2003) 341-375.

    - [KZh2015b] M. Koeppe and Y. Zhou, An electronic compendium of extreme functions for the
      Gomory-Johnson infinite group problem, Operations Research Letters, 2015,
      http://dx.doi.org/10.1016/j.orl.2015.06.004
    """

    def _construct_function(self, f=1/12, bkpt=2/12, field=None, conditioncheck=True):
        if not bool(0 < f < bkpt < 1 + f - bkpt < 1):
            raise ValueError("Bad parameters. Unable to construct the function.")
        claimed_parameter_attribute = None
        if conditioncheck:
            #if not ((f in QQ) & (bkpt in QQ) & bool(0 < f < bkpt < ((1 + f)/4) < 1)):
            if not bool(0 < f < bkpt <= ((1 + f)/4) < 1):
                logging.info("Conditions for extremality are NOT satisfied.")
                claimed_parameter_attribute = 'constructible'
            else:
                logging.info("Conditions for extremality are satisfied.")
                claimed_parameter_attribute = 'extreme'
        bkpts = [0, f, bkpt, 1 + f - bkpt, 1]
        # values = [0, 1, bkpt/(1 + f), (1 + f - bkpt)/(1 + f),0]
        slopes = [1/f, (1 + f - bkpt)/(1 + f)/(f - bkpt), 1/(1 + f), (1 + f - bkpt)/(1 + f)/(f - bkpt)]
        h = piecewise_function_from_breakpoints_and_slopes(bkpts, slopes, field=field)
        h._claimed_parameter_attribute = claimed_parameter_attribute
        return h

drlm_backward_3_slope = ParametricFamily_drlm_backward_3_slope()

class ParametricFamily_dg_2_step_mir_limit(ParametricFamily):

    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = dg_2_step_mir_limit()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - Name: Dash-Gunluk 2-Step MIR Limit;
        - Infinite; Dim = 1; Slopes = 1; Discontinuous; Simple sets method;
        - Discovered [33] p.41, def.12;
        - Proven extreme [33] p.43, lemma 14.
        - dg_2_step_mir_limit is a facet.

    Parameters:
        * f (real) `\in (0,1)`;
        * d (positive integer): number of slopes on [0,f).

    Function is known to be extreme under the conditions:
        * 0 < f < 1;
        * d >= ceil(1 / (1 - f)) - 1.

    Note:
        This is the limit function as alpha in dg_2_step_mir()
        tends (from left) to f/d, where d is integer;
        cf. [33] p.42, lemma 13.

        It's a special case of ``drlm_2_slope_limit()``:

        ``dg_2_step_mir_limit(f, d) = multiplicative_homomorphism(drlm_2_slope_limit(f=1-f, nb_pieces_left=1, nb_pieces_right=d), -1)``.

    Examples: [33] p.42, Fig.6 ::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: h = dg_2_step_mir_limit(f=3/5, d=3)
        sage: extremality_test(h, False)
        True

    Reference:
        [33]: S. Dash and O. Gunluk, Valid inequalities based on simple mixed-integer sets.,
                Proceedings 10th Conference on Integer Programming and Combinatorial Optimization
                (D. Bienstock and G. Nemhauser, eds.), Springer-Verlag, 2004, pp. 33-45.
    """

    def _claimed_parameter_attribute(self, f, d, **kwargs):
        if not (bool(0 < f < 1) & (d >= 1)):
            return 'not_constructible'
        if not bool(d >= ceil(1 / (1 - f)) - 1):
            return 'constructible'
        else:
            return 'extreme'

    def _construct_function(self, f=3/5, d=3, field=None, **kwds):
        if not (bool(0 < f < 1) & (d >= 1)):
            raise ValueError("Bad parameters. Unable to construct the function.")
        f = nice_field_values([f], field)[0]
        field = f.parent()
        pieces = []
        for k in range(d):
            left_x = f * k / d
            right_x = f * (k + 1) / d
            pieces = pieces + \
                     [[singleton_interval(left_x), FastLinearFunction(field(0), left_x / f)], \
                      [open_interval(left_x, right_x), FastLinearFunction(1 / (f - 1), (k + 1)/(d + 1)/(1 - f))]]
        pieces.append([closed_interval(f, field(1)), FastLinearFunction(1 / (f - 1), 1 /(1 - f))])
        return PiecewiseLinearFunction_1d(pieces)

dg_2_step_mir_limit = ParametricFamily_dg_2_step_mir_limit()

class ParametricFamily_drlm_2_slope_limit(ParametricFamily):

    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = drlm_2_slope_limit()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - Name: Dey--Richard--Li--Miller's 2-Slope Limit;
        - Infinite; Dim = 1; Slopes = 1; Discontinuous; Group relations method;
        - Discovered [40] p.158 def.10;
        - Proven extreme [40] p.159 thm.8.
        - (Although only extremality has been established in literature, the same proof shows that) drlm_2_slope_limit is a facet.

    Parameters:
        * f (real) `\in (0,1)`;
        * nb_pieces_left (positive integer) : number of linear pieces to the left of f;
        * nb_pieces_right (positive integer) : number of linear pieces to the right of f.

    Function is known to be extreme under the conditions:
        nb_pieces_left * (1-f) <= nb_pieces_right * f.

    Examples:
        [40] p.159 Fig.4 ::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
            sage: h = drlm_2_slope_limit(f=3/5, nb_pieces_left=3, nb_pieces_right=4)
            sage: extremality_test(h, False)
            True

    Reference:
        [40]: S.S. Dey, J.-P.P. Richard, Y. Li, and L.A. Miller, On the extreme inequalities of infinite group problems,
                Mathematical Programming 121 (2010) 145-170.
    """

    def _construct_function(self, f=3/5, nb_pieces_left=3, nb_pieces_right=4, field=None, conditioncheck=True):
        m = nb_pieces_left
        d = nb_pieces_right
        if not ((m in ZZ) & (d in ZZ) & (m >= 1) & (d >= 1) & bool(0 < f < 1)):
            raise ValueError("Bad parameters. Unable to construct the function.")
        claimed_parameter_attribute = None
        if conditioncheck:
            if not bool(m*(1 - f) <= d*f):
                logging.info("Conditions for extremality are NOT satisfied.")
                claimed_parameter_attribute = 'constructible'
            else:
                logging.info("Conditions for extremality are satisfied.")
                claimed_parameter_attribute = 'extreme'
        s = (m + d)/((d + 1)*f - (m - 1)*(1 - f))
        delta_2 = (s - s*f + 1)/(d + 1)
        if m == 1:
            delta_1 = 0
        else:
            delta_1 = (s*f - 1)/(m - 1)
        # in irrational case, try to coerce to common number field
        [f, s, delta_1, delta_2, m, d] =  nice_field_values([f, s, delta_1, delta_2, m, d], field)
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
        psi = FastPiecewise(pieces)
        psi._claimed_parameter_attribute = claimed_parameter_attribute
        return psi

drlm_2_slope_limit = ParametricFamily_drlm_2_slope_limit()

class ParametricFamily_drlm_3_slope_limit(ParametricFamily):

    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = drlm_3_slope_limit()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - Name: Dey--Richard--Li--Miller's 3-Slope Limit;
        - Infinite; Dim = 1; Slopes = 2; Discontinuous; Group relations method;
        - Discovered [40] p.161 def.11;
        - Proven extreme [40] p.161 thm.9.
        - (Although only extremality has been established in literature, the same proof shows that) drlm_3_slope_limit is a facet.

    Parameters:
        f (real) `\in (0,1)`;

    Function is known to be extreme under the conditions:
        0 < f < 1/3.

    Note:
        This is the limit function as bkpt tends to f in drlm_backward_3_slope(f, bkpt).

    Examples:
        [40] p.162 Fig.5 ::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
            sage: h = drlm_3_slope_limit(f=1/5)
            sage: extremality_test(h, False)
            True

    Reference:
        [40]: S.S. Dey, J.-P.P. Richard, Y. Li, and L.A. Miller, On the extreme inequalities of infinite group problems,
                Mathematical Programming 121 (2010) 145-170.
    """

    def _construct_function(self, f=1/5, field=None, conditioncheck=True):
        if not bool(0 < f < 1):
            raise ValueError("Bad parameters. Unable to construct the function.")
        claimed_parameter_attribute = None
        if conditioncheck:
            if not bool(0 < f < 1/3):
                logging.info("Conditions for extremality are NOT satisfied.")
                claimed_parameter_attribute = 'constructible'
            else:
                logging.info("Conditions for extremality are satisfied.")
                claimed_parameter_attribute = 'extreme'
        f = nice_field_values([f], field)[0]
        field = f.parent()
        pieces = [[closed_interval(0, f), FastLinearFunction(1/f, 0)], \
                  [open_interval(f, 1), FastLinearFunction(1/(f + 1), 0)], \
                  [singleton_interval(field(1)), FastLinearFunction(field(0), 0)]]
        kappa = FastPiecewise(pieces)
        kappa._claimed_parameter_attribute = claimed_parameter_attribute
        return kappa

drlm_3_slope_limit = ParametricFamily_drlm_3_slope_limit()

class ParametricFamily_bccz_counterexample(ParametricFamily):

    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = bccz_counterexample()
        h = psi_n_in_bccz_counterexample_construction(e=generate_example_e_for_psi_n(n=7))
        g = plot(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=1, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    return function psi, a counterexample to Gomory--Johnson's conjecture
    constructed by Basu--Conforti--Cornuejols--Zambelli [IR1].

    psi is a continuous facet (hence extreme), but is not piecewise linear. cf. [IR1]

    It can be considered as an absolutely continuous, measurable, non-piecewise linear
    "2-slope function".  A separate case with different parameters, which gives rise to
    a continuous "1-slope function", is discussed in [KZh2015b, section 5].

    Parameters:

    - f (real) `\in (0,1)`;
    - q (real), q > 2: ratio of the geometric series;
    - eta (real), 0 <= eta < 1: to control the series sum;
    - maxiter (integer): maximum number of iterations;

    Note:
        psi is the uniform limit of the sequence of functions psi_n,
        generated by psi_n_in_bccz_counterexample_construction(f, [e[0], e[1], ..., e[n - 1]]).
        e is a geometric series with ratio q, such that:
        0 < ... < e[n] <= e[n - 1] <= ... <= e[1] <= e[0] <= 1 - f and \sum_{i = 0}^{\infty} {2^i * e[i]} <= f.
        The first n terms of e are generated by generate_example_e_for_psi_n(f, n, q, eta)
        
    See also::

        def generate_example_e_for_psi_n(f, n q, eta)
        def psi_n_in_bccz_counterexample_construction(f, e)

    Examples:
        quick exact evaluations::

            sage: from cutgeneratingfunctionology.igp import *
            sage: bccz_counterexample(f=2/3, q=4, eta=0, maxiter=10000)(r=1/5)
            21/40
            sage: bccz_counterexample(f=2/3, q=4, eta=0, maxiter=10000)(r=1/4)
            3/4

        too many iterations::

            sage: bccz_counterexample(f=2/3, q=4, eta=0, maxiter=10000)(r=9/40) # doctest: +SKIP

    References:

    - [IR1]:  A. Basu, M. Conforti, G. Cornuejols, and G. Zambelli, A counterexample to a conjecture of Gomory and Johnson,
      Mathematical Programming Ser. A 133 (2012), 25-38.

    - [KZh2015b] M. Koeppe and Y. Zhou, An electronic compendium of extreme functions for the
      Gomory-Johnson infinite group problem, Operations Research Letters, 2015,
      http://dx.doi.org/10.1016/j.orl.2015.06.004
    """

    def _construct_function(self, f=2/3, q=4, eta=1/1000, maxiter=10000):
        if not (bool(0 < f < 1) & bool(q > 2) & bool(0 <= eta < 1)):
            raise ValueError("Bad parameters.")
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
                    logging.warning("Reaching max number of iterations, return approximate psi(%s)" %r)
                return (y_left + y_right)/2 - (r - (x_left + x_right)/2) / (1 - f)
            else:
                raise ValueError("outside domain")
        logging.warning("This function is not piecewise linear; code for handling this function is not implemented.")
        return evaluate_psi_at_r

bccz_counterexample = ParametricFamily_bccz_counterexample()

def generate_example_e_for_psi_n(f=2/3, n=7, q=4, eta=1/1000):
    r"""
    Return the first n terms of a geometric series e that satisfies 
    0 < ... < e[n] <= e[n - 1] <= ... <= e[1] <= e[0] <= 1 - f and \sum_{i = 0}^{\infty} {2^i * e[i]} <= f.

    This can be used in psi_n_in_bccz_counterexample_construction(f, [e[0],...,e[n-1]]), so that the function constructed is extreme.

    Parameters:
        - f (real) `\in (0,1)`;
        - n (integer);
        - q (real), q > 2: ratio of the geometric series;
        - eta (real), 0 <= eta < 1: to control the series sum, \sum_{i = 0}^{\infty} {2^i * e[i]} <= (1 - eta)*f.

    Note:
        If (eta == 0) and (f >= 1/2), then \sum_{i = 0}^{\infty} {2^i * e[i]} = f.

    This case is not mentioned in [IR1], but using a similar proof, one can show that:
        (1) psi_n still converges uniformly to psi;
        (2) The limit funciton psi is a continuous facet (hence extreme);
        (3) psi is not piecewise linear. 

    Also notice that: psi is not in W^{1,1}. See [KZh2015b, section 5].

    References:

    - [IR1] A. Basu, M. Conforti, G. Cornuejols, and G. Zambelli, A counterexample to a conjecture of Gomory and Johnson,
      Mathematical Programming Ser. A 133 (2012), 25-38.

    - [KZh2015b] M. Koeppe and Y. Zhou, An electronic compendium of extreme functions for the
      Gomory-Johnson infinite group problem, Operations Research Letters, 2015,
      http://dx.doi.org/10.1016/j.orl.2015.06.004
    """
    if n == 0:
        return []
    if not (bool(0 < f < 1) & bool(q > 2) & bool(0 <= eta < 1)):
        raise ValueError("Bad parameters.")
    x = (1 - eta)*(q - 2) / q * min(f, 1 - f)
    # or take x = min((1 - eta)*(q - 2)*f / q , 1 - f) 
    e = [x / q^i for i in range(n)]
    return e

class ParametricFamily_psi_n_in_bccz_counterexample_construction(ParametricFamily):

    r"""
    Summary: 
        - Name: psi_n in the construction of BCCZ's counterexample to GJ's conjecture;
        - Infinite; Dim = 1; Slopes = 2; Continuous;  Analysis of subadditive polytope method.
        - Discovered [IR1]  p.30, section.3, fig.1;
        - Proven extreme [IR1] p.35, thm.4.7.

    Note:
        The (uniform) limit `\psi = \lim_{n to \infty} \psi_n` is well defined if `\sum_{i = 0}^{\infty} {2^i * e[i]} < f`.

        The (uniform) limit \psi is a continuous facet, but is not piecewise linear. A counterexample of GJ's conjecture.

        Could use the function generate_example_e_for_psi_n(f, n, q) to generate a sequence e that satisfies the conditions for extremality.

        psi_n_in_bccz_counterexample_construction() is a special case of kf_n_step_mir(), with f=f, a = [1, (f + e[0])/2, (f - e[0] + 2*e[1])/4, ...]

    Parameters:
        * f (real) `\in (0,1)`;
        * e (list of reals, with length = n) `\in (0,f)`.

    Function is known to be extreme under the conditions:
        * 0 < f < 1;
        * 0 < e[n - 1] <= e[n - 2] <= ... <= e[1] <= e[0] <= 1 - f;
        * \sum_{i = 0}^{n - 1} 2^i * e[i] < f.

    Examples:
        [IR1]  p.30, fig.1::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: h = psi_n_in_bccz_counterexample_construction(f=2/3, e=[1/12, 1/24])
            sage: extremality_test(h, False)
            True
            sage: h = psi_n_in_bccz_counterexample_construction(f=4/5, e=[1/5, 1/20, 1/80, 1/320, 1/1280])
            sage: extremality_test(h, False, f=4/5)         # Suppress warning about non-unique f
            True
            sage: h = psi_n_in_bccz_counterexample_construction(f=4/5, e=generate_example_e_for_psi_n(f=4/5, n=7, q=3, eta=0))
            sage: extremality_test(h, False)
            True
            sage: sum([plot(psi_n_in_bccz_counterexample_construction(e=generate_example_e_for_psi_n(n=n)), color=color, legend_label="psi_%d"%n) for n, color in zip(range(7),rainbow(7))]) # not tested

    Reference: 
        [IR1]:  A. Basu, M. Conforti, G. Cornuejols, and G. Zambelli, A counterexample to a conjecture of Gomory and Johnson,
                    Mathematical Programming Ser. A 133 (2012), 25-38.
    """

    def _construct_function(self, f=2/3, e=(1/12, 1/24), field=None, conditioncheck=True):
        if not bool(0 < f < 1):
            raise ValueError("Bad parameters. Unable to construct the function.")
        claimed_parameter_attribute = None
        n = len(e)
        if n == 0:
            if conditioncheck:
                logging.info("Conditions for extremality are satisfied.")
                claimed_parameter_attribute = 'extreme'
            h = piecewise_function_from_breakpoints_and_values([0,f,1], [0,1,0])
            h._claimed_parameter_attribute = claimed_parameter_attribute
            return h
        a = [1]
        b = [f]
        sum_e = 0
        if not bool(0 < e[0]):
            raise ValueError("Bad parameters. Unable to construct the function.")
        t = bool(e[0] <= 1 - f)
        for i in range(0, n):
            a.append((b[i] + e[i]) / 2)
            b.append((b[i] - e[i]) / 2)
            sum_e = sum_e + (2^i) * e[i]
            if not (bool(e[i] > 0) & bool(sum_e < f)):
                raise ValueError("Bad parameters. Unable to construct the function.")
            if not (i == 0) | bool(e[i] <= e[i-1]):
                t = False
        if conditioncheck:
            if t:
                logging.info("Conditions for extremality are satisfied.")
                claimed_parameter_attribute = 'extreme'
            else:
                logging.info("Conditions for extremality are NOT satisfied.")
                claimed_parameter_attribute = 'constructible'
        interval_lengths =  interval_length_n_step_mir(n + 1, 1, a, b)
        nb_interval = len(interval_lengths)
        interval_length_positive = sum(interval_lengths[i] for i in range(0, nb_interval, 2))
        interval_length_negative = sum(interval_lengths[i] for i in range(1, nb_interval, 2))
        s_negative = a[0] /(b[0] - a[0])
        s_positive = - s_negative * interval_length_negative / interval_length_positive 
        slopes = [s_positive, s_negative] * (nb_interval // 2)
        h = piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes, field=field)
        h._claimed_parameter_attribute = claimed_parameter_attribute
        return h

psi_n_in_bccz_counterexample_construction = ParametricFamily_psi_n_in_bccz_counterexample_construction()

class ParametricFamily_bhk_irrational(ParametricFamily):

    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = bhk_irrational()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=1, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Summary:
        - Name: Basu-Hildebrand-Koeppe's irrational function.
        - Infinite; Dim = 1; Slopes = 3; Continuous;  Covered intervals and equivariant perturbation.
        - Discovered :cite:`basu-hildebrand-koeppe:equivariant` p.33, section.5.2, fig.9-10.
        - Proven extreme :cite:`basu-hildebrand-koeppe:equivariant` p.34, thm.5.3.
        - (Although only extremality has been established in literature, the same proof shows that), bhk_irrational is a facet.

    Parameters:
        * f (real) `\in (0,1)`;
        * d1 (real): length of the positive slope;
        * d2 (real): length of the zero slopes;
        * a0 (real): length of the first zig-zag;
        * delta (n-tuple of reals): length of the extra zig-zags.

    Function is known to be extreme under the conditions:
        * 0 < f < 1;
        * d1, d2, a0, delta > 0;
        * d1 + d2 < f;
        * len(delta) == 2
        * sum(delta) <  d2 / 4; Weaker condition: 2*delta[0] + delta[1] < d2 / 2;
        * the two components of delta are linearly independent over \Q.

    Relation between the code parameters and the paper parameters:
        * t1 = delta[0], t2 = delta[0] + delta[1], ...
        * a1 = a0 + t1, a2 = a0 + t2, ...
        * A = f/2 - a0/2 - d2/4,
        * A0 = f/2 - a0/2 + d2/4.

    Examples:
        :cite:`basu-hildebrand-koeppe:equivariant`  p.34, thm.5.3::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.NOTSET) # enable INFO messages disabled by other doctests
            sage: h = bhk_irrational(f=4/5, d1=3/5, d2=1/10, a0=15/100, delta=(1/200, sqrt(2)/200))
            INFO: ...
            sage: extremality_test(h, False)
            INFO: ...
            True

        :cite:`basu-hildebrand-koeppe:equivariant`  thm 5.4: Not extreme for rational data::

            sage: h = bhk_irrational(delta=[1/200, 3/200])
            INFO: ...
            sage: extremality_test(h, False)
            INFO: ...
            INFO: ... Total: 1 stability orbit...
            False

        A generalization with 3 zigzags instead of 2 as in :cite:`basu-hildebrand-koeppe:equivariant`::

            sage: h = bhk_irrational(f=4/5, d1=3/5, d2=1/10, a0=15/100, delta=(1/200, 6* sqrt(2)/200, 1/500))
            INFO: ...
            sage: extremality_test(h, False) # not tested - takes 20min on Macbook Pro 2.7Ghz Core i7
            INFO: ...
            INFO: ... Total: 3 stability orbits...
            False

        Verify that p constructed below is an effective perturbation of h::

            sage: sqrt2 = h.end_points()[0].parent().gen()
            sage: pb = ((0, 0), (343/1000, -3/100), (1331/4000, -9/400), (1249/4000, -3/400), (151/500, 0), (303/1000, 0), (1171/4000, 3/400), (1089/4000, 9/400), (131/500, 3/100), (87/250, -3/100), (1351/4000, -9/400), (1269/4000, -3/400), (307/1000, 0), (77/250, 0), (1191/4000, 3/400), (1109/4000, 9/400), (267/1000, 3/100), (353/1000, -3/100), (1371/4000, -9/400), (1289/4000, -3/400), (39/125, 0), (313/1000, 0), (1211/4000, 3/400), (1129/4000, 9/400), (34/125, 3/100), (179/500, -3/100), (1391/4000, -9/400), (1309/4000, -3/400), (317/1000, 0), (159/500, 0), (1231/4000, 3/400), (1149/4000, 9/400), (277/1000, 3/100), (363/1000, -3/100), (1411/4000, -9/400), (1329/4000, -3/400), (161/500, 0), (323/1000, 0), (1251/4000, 3/400), (1169/4000, 9/400), (141/500, 3/100), (46/125, -3/100), (1431/4000, -9/400), (1349/4000, -3/400), (327/1000, 0), (41/125, 0), (1271/4000, 3/400), (1189/4000, 9/400), (287/1000, 3/100), (373/1000, -3/100), (1451/4000, -9/400), (1369/4000, -3/400), (83/250, 0), (333/1000, 0), (1291/4000, 3/400), (1209/4000, 9/400), (73/250, 3/100), (189/500, -3/100), (1471/4000, -9/400), (1389/4000, -3/400), (337/1000, 0), (169/500, 0), (1311/4000, 3/400), (1229/4000, 9/400), (297/1000, 3/100), (383/1000, -3/100), (1491/4000, -9/400), (1409/4000, -3/400), (171/500, 0), (343/1000, 0), (1331/4000, 3/400), (1249/4000, 9/400), (151/500, 3/100), (97/250, -3/100), (1511/4000, -9/400), (1429/4000, -3/400), (347/1000, 0), (87/250, 0), (1351/4000, 3/400), (1269/4000, 9/400), (307/1000, 3/100), (493/1000, -3/100), (1931/4000, -9/400), (1849/4000, -3/400), (113/250, 0), (453/1000, 0), (1771/4000, 3/400), (1689/4000, 9/400), (103/250, 3/100), (249/500, -3/100), (1951/4000, -9/400), (1869/4000, -3/400), (457/1000, 0), (229/500, 0), (1791/4000, 3/400), (1709/4000, 9/400), (417/1000, 3/100), (503/1000, -3/100), (1971/4000, -9/400), (1889/4000, -3/400), (231/500, 0), (463/1000, 0), (1811/4000, 3/400), (1729/4000, 9/400), (211/500, 3/100), (127/250, -3/100), (1991/4000, -9/400), (1909/4000, -3/400), (467/1000, 0), (117/250, 0), (1831/4000, 3/400), (1749/4000, 9/400), (427/1000, 3/100), (513/1000, -3/100), (2011/4000, -9/400), (1929/4000, -3/400), (59/125, 0), (473/1000, 0), (1851/4000, 3/400), (1769/4000, 9/400), (54/125, 3/100), (259/500, -3/100), (2031/4000, -9/400), (1949/4000, -3/400), (477/1000, 0), (239/500, 0), (1871/4000, 3/400), (1789/4000, 9/400), (437/1000, 3/100), (523/1000, -3/100), (2051/4000, -9/400), (1969/4000, -3/400), (241/500, 0), (483/1000, 0), (1891/4000, 3/400), (1809/4000, 9/400), (221/500, 3/100), (66/125, -3/100), (2071/4000, -9/400), (1989/4000, -3/400), (487/1000, 0), (61/125, 0), (1911/4000, 3/400), (1829/4000, 9/400), (447/1000, 3/100), (533/1000, -3/100), (2091/4000, -9/400), (2009/4000, -3/400), (123/250, 0), (493/1000, 0), (1931/4000, 3/400), (1849/4000, 9/400), (113/250, 3/100), (269/500, -3/100), (2111/4000, -9/400), (2029/4000, -3/400), (497/1000, 0), (249/500, 0), (1951/4000, 3/400), (1869/4000, 9/400), (457/1000, 3/100), (1, 0))
            sage: pv = ((0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0))
            sage: bkpts = [b[0]+b[1]*sqrt2 for b in pb]
            sage: values = [v[0]+v[1]*sqrt2 for v in pv]
            sage: p = piecewise_function_from_breakpoints_and_values(bkpts, values)
            sage: h1 = h + 8/10000*p
            sage: minimality_test(h1)
            INFO: ...
            True
            sage: h2 = h - 8/10000*p
            sage: minimality_test(h2)
            INFO: ...
            True
    """

    def _construct_function(self, f=4/5, d1=3/5, d2=1/10, a0=15/100, delta=(1/200, sqrt(2)/200), field=None):
        if not (bool(0 < f < 1) and bool(d1 > 0) and bool(d2 > 0) and bool(a0 > 0) 
                and all(bool(deltai > 0) for deltai in delta) and bool(d1 + d2 < f) and (sum(delta) < f/2 - d2/4 - 3*a0/2) ):
            raise ValueError("Bad parameters. Unable to construct the function.")
        if len(delta) < 2:
            logging.info("Conditions for extremality are NOT satisfied.")
        elif len(delta) == 2:
            if is_QQ_linearly_independent(*delta) and  2*delta[0] + delta[1] < d2 / 2:
                logging.info("Conditions for extremality are satisfied if it is a minimal function.")
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
        return piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes, field=field)

bhk_irrational = ParametricFamily_bhk_irrational()

class ParametricFamily_bhk_slant_irrational(ParametricFamily):

    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = bhk_slant_irrational()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=1, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    A version of the irrational function with non-zero second slope

    Parameters:
        * f (real) `\in (0,1)`;
        * d1 (real): length of the positive slopes;
        * d2 (real): length of the slant (c2) slopes;
        * a0 (real): length of the first zig-zag;
        * delta (n-tuple of reals): length of the extra zig-zags.
        * c2 (real): slant slope, c2 = 0 in bhk_irrational() 

    Function is known to be extreme under the conditions:
        * 0 < f < 1;
        * d1, d2, a0, delta > 0;
        * d1 + d2 < f;
        * len(delta) == 2
        * sum(delta) <  d2 / 4; Weaker condition: 2*delta[0] + delta[1] < d2 / 2;
        * the two components of delta are linearly independent over \Q.
        * (?? sufficient ??) -1 / (1 - f) <= c2 <= (1 - d1 - d2) / (d1 + d2) / (1 - f);
        
    Also needs:
        * sum(delta) < (f/2 - d2/4 - 3*a0/2) * c1 / (c1-c2);
        * 5*a0 > 2*(c1-c2)/(c1-c3) * d2 + d2 / 2 + d1.
        * d?? > 0
        * ...

    Relation between the code parameters and the paper parameters:
        * t1 = delta[0], t2 = delta[0] + delta[1], ...
        * a1 = a0 + t1, a2 = a0 + t2, ...
        * A = f/2 - a0/2 - d2/4,
        * A0 = f/2 - a0/2 + d2/4.

    Examples::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = bhk_slant_irrational(f=4/5, d1=3/5, d2=1/10, a0=15/100, delta=(1/200, sqrt(2)/200), c2=1/16)
        sage: h2 = bhk_gmi_irrational(f=4/5, d1=3/5, d2=1/10, a0=15/100, delta=(1/200, sqrt(2)/200), alpha=95/100)
        sage: h == h2
        True

    Bug example (function is not minimal)::

        sage: h3 = bhk_irrational(f=4/5, d1=3/5, d2=1/10, a0=40/180, delta=(1/400, sqrt(2)/400))
        sage: minimality_test(h3)
        False
    """

    def _construct_function(self, f=4/5, d1=3/5, d2=1/10, a0=15/100, delta=(1/200, sqrt(2)/200), c2=0, field=None):
        if not (bool(0 < f < 1) and bool(d1 > 0) and bool(d2 > 0) and bool(a0 > 0)
                and all(bool(deltai > 0) for deltai in delta) and bool(d1 + d2 < f)):
            raise ValueError("Bad parameters. Unable to construct the function.")

        d3 = f - d1 - d2
        c3 = -1/(1-f)
        c1 = (1-d2*c2-d3*c3)/d1
        d21 = d2 / 2
        d31 = (c1 - c2) / (c1 - c3) * d21
        d11 = a0 - d31
        d13 = a0 - d21
        d12 = (d1 - d13)/2 - d11
        d32 = d3/2 - d31

        if bool(d32 < 0):
            raise ValueError("Bad parameters. Unable to construct the function. ")

        a0_max = min(f - d2/2, d1 + d2/2 + 2*d31) / 3  # since A > 0 and d12 > 0

        if not bool(d31 < a0 < a0_max):
            raise ValueError("Bad parameters. %s < a0 < %s is not satisfied. Unable to construct the function." % (d31, a0_max))

        sumdelta_max = (f/2 - d2/4 - 3*a0/2) * c1 / (c1 - c2) # since d32 > 0
        if not bool(sum(delta) < sumdelta_max):
            raise ValueError("Bad parameters. sum(delta) < %s is not satisfied. Unable to construct the function." % sumdelta_max)

        a0_min = (2*(c1-c2)/(c1-c3) * d2 + d2 / 2 + d1) / 5 # since d11 >= d12
        c2_min = -1 / (1 - f) # since c2 >= c3
        c2_max = (1 - d1 - d2) / (d1 + d2) / (1 - f) # since c2 <= c1
        if bool(a0 < a0_min):
            logging.info("Conditions for extremality are NOT satisfied. Minimality requires a0 >= %s" % a0_min)
        elif not bool(c2_min <= c2 <= c2_max):
            logging.info("Conditions for extremality are NOT satisfied. Need %s < c2 < %s" % (c2_min, c2_max))
        elif len(delta) < 2:
            logging.info("Conditions for extremality are NOT satisfied.")
        elif len(delta) == 2:
            if is_QQ_linearly_independent(*delta) and  2*delta[0] + delta[1] < d2 / 2:
                logging.info("Conditions for extremality are satisfied if it is a minimal function.")
            else:
                logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("len(delta) >= 3, Conditions for extremality are unknown.")

        zigzag_lengths = []
        zigzag_slopes = []
        delta_positive = 0
        delta_negative = 0
        for delta_i in delta:
            delta_i_negative = (c1 - c2) * delta_i / (c1 - c3)
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
        return piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes, field=field)

bhk_slant_irrational = ParametricFamily_bhk_slant_irrational()

class ParametricFamily_bhk_gmi_irrational(ParametricFamily):

    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = bhk_gmi_irrational()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=1, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    A version of the irrational function with non-zero second slope,
    obtained by forming a convex combination of a modified version of the irrational function with the GMI cut.
    Constructed by Chun Yu Hong, 2013.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = bhk_gmi_irrational()
        sage: extremality_test(h, False)
        True
    """

    def _construct_function(self, f=4/5, d1=3/5, d2=1/10, a0=15/100, delta=(1/200, sqrt(2)/200), alpha=95/100, field=None):
        if not (bool(0 < f < 1) and bool(d1 > 0) and bool(d2 > 0) and bool(a0 > 0) and (len(delta) >= 2) \
                and bool(min(delta) > 0) and bool(d1 + d2 < f) and (sum(delta) < f/2 - d2/4 - 3*a0/2) \
                and bool(0 < alpha < 1)):
            raise ValueError("Bad parameters. Unable to construct the function.")
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

        bhk = piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes, field=field)
        gmi = gmic(f, field=field)
        return alpha * bhk + (1 - alpha) * gmi

bhk_gmi_irrational = ParametricFamily_bhk_gmi_irrational()

class ParametricFamily_chen_4_slope(ParametricFamily):

    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = chen_4_slope()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    This 4-slope function is shown [KChen_thesis] to be a facet.

    Parameters:
        * f (real) `\in (0,1)`;
        * s_pos, s_neg (real): positive slope and negative slope
        * lam1, lam2 (real).

    Function is claimed to be extreme under the following conditions,
    according to literature [KChen_thesis]:
        * 1/2 <= f < 1;
        * s_pos >= 1/f;
        * s_neg <= 1/(f - 1);
        * 0 <= lam1 < min(1/2, (s_pos - s_neg) / s_pos / (1 - s_neg * f));
        * f - 1 / s_pos < lam2 < min(1/2, (s_pos - s_neg) / s_neg / (s_pos * (f - 1) - 1)).

    Note:
        The lower bound 0 < lam1 claimed in [KChen_thesis] is not sufficient for extremality.
        The following example satisfies the extremality conditions according to [KChen_thesis],
        however it violates the subadditivity, and thus is not an extreme function::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: h = chen_4_slope(f=7/10, s_pos=2, s_neg=-4, lam1=1/100, lam2=49/100, condition_according_to_literature=True)
            sage: h._claimed_parameter_attribute
            'extreme'
            sage: extremality_test(h, False)
            False

        On the other hand, the hypotheses stated by Chen are also not necessary for extremality.
        For example, the following function does not satisfy the hypotheses, however it is extreme::

            sage: h = chen_4_slope(f=7/10, s_pos=2, s_neg=-4, lam1=1/10, lam2=1/10, condition_according_to_literature=True)
            sage: h._claimed_parameter_attribute
            'constructible'
            sage: extremality_test(h, False)
            True

    We propose to revised the conditions for extremality, as follows.
        s_pos >= 1/f;
        s_neg <= 1/(f - 1);
        lam1 <= 1/2;
        lam2 <= 1/2;
        ((f*s_pos-1) * (1-f*s_neg)) * lam1 <= (s_pos-s_neg) * lam2;
        ((f-1)*s_neg-1) * (1+(1-f)*s_pos) * lam2 <= (s_pos-s_neg) * lam1.

    Examples:
        [KChen_thesis]  p.38, fig.8::

            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: h = chen_4_slope(f=7/10, s_pos=2, s_neg=-4, lam1=1/4, lam2=1/4)
            sage: extremality_test(h, False)
            True

            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: h = chen_4_slope(f=1/2, s_pos=4, s_neg=-4, lam1=1/3, lam2=1/3)
            sage: extremality_test(h, False)
            True

        The following parameters do not satisfy the requirement, however the function is extreme::

            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: h = chen_4_slope(f=1/2, s_pos=5, s_neg=-5, lam1=1/5, lam2=1/5)
            sage: extremality_test(h, False)
            True

    Reference:
        [KChen_thesis]:  K. Chen, Topics in group methods for integer programming,
                            Ph.D. thesis, Georgia Institute of Technology, June 2011.
    """

    def _construct_function(self, f=7/10, s_pos=2, s_neg=-4, lam1=1/4, lam2=1/4, field=None, conditioncheck=True, condition_according_to_literature=False, merge=True):     
        if not (bool(0 < f < 1) and bool(s_pos >= 1/f) and bool(s_neg <= 1/(f - 1)) \
                                and bool(0 <= lam1 <= 1) and bool(0 <= lam2 <= 1)):
            raise ValueError("Bad parameters. Unable to construct the function.")
        claimed_parameter_attribute = None
        if conditioncheck:
            if condition_according_to_literature:
                if bool(1/2 <= f) and bool(lam1 < 1/2) and bool(lam2 < 1/2) and \
                   bool(0 <= lam1  < (s_pos - s_neg) / s_pos / (1 - s_neg * f)) and \
                   bool (f - 1 / s_pos < lam2 < (s_pos - s_neg) / s_neg / (s_pos * (f - 1) - 1)):
                    claimed_parameter_attribute = 'extreme'
                else:
                    claimed_parameter_attribute = 'constructible'
            else:
                if bool(lam1 <= 1/2) and bool(lam2 <= 1/2) \
               and bool(-f^2*s_pos*s_neg*lam2 + 2*f*s_pos*s_neg*lam2 + f*s_pos*lam2 + f*s_neg*lam2 - s_pos*s_neg*lam2 - s_pos*lam1 + s_neg*lam1 - s_pos*lam2 - s_neg*lam2 - lam2 <= 0) \
               and bool(-f^2*s_pos*s_neg*lam1 + f*s_pos*lam1 + f*s_neg*lam1 - s_pos*lam2 + s_neg*lam2 - lam1 <= 0):
                    claimed_parameter_attribute = 'extreme'
                else:
                    claimed_parameter_attribute = 'constructible'
            if claimed_parameter_attribute == 'extreme':
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
        h = piecewise_function_from_breakpoints_and_slopes([0, aa, a, b, bb, f, dd, d, c, cc, 1], slopes, field=field, merge=merge)
        h._claimed_parameter_attribute = claimed_parameter_attribute
        return h

chen_4_slope = ParametricFamily_chen_4_slope()

class ParametricFamily_chen_4_slope_reworded(ParametricFamily):

    def _construct_function(self, f=7/10, aa=19/240, a=7/80, c=77/80, cc=29/30,field=None, conditioncheck=True, condition_according_to_literature=False, merge=False):
        if not (bool(0 < aa < a < f/2) and bool((1+f)/2 < c < cc < 1)):
            raise ValueError("Bad parameters. Unable to construct the function.")
        claimed_parameter_attribute = None
        v = (a*aa*cc - a*aa - (a*aa*cc - aa^2 - (a*aa - aa^2)*c)*f)/(((a - aa)*c - a*cc + aa)*f^2 - ((a - aa)*c - a*cc + aa)*f)
        w = (a*cc^2 + a*c - (a*c + a)*cc - (a*cc^2 + (a - aa)*c - ((a - aa)*c + a + aa)*cc + aa)*f)/(((a - aa)*c - a*cc + aa)*f^2 - ((a - aa)*c - a*cc + aa)*f)
        lam1 = 2*a/f
        lam2 = 2*(1-c)/(1-f)
        s_pos = v/aa
        s_neg = -w/(1-cc)
        if conditioncheck:
            if condition_according_to_literature:
                if bool(1/2 <= f) and bool(lam1 < 1/2) and bool(lam2 < 1/2) and \
                   bool(0 <= lam1  < (s_pos - s_neg) / s_pos / (1 - s_neg * f)) and \
                   bool (f - 1 / s_pos < lam2 < (s_pos - s_neg) / s_neg / (s_pos * (f - 1) - 1)):
                    claimed_parameter_attribute = 'extreme'
                else:
                    claimed_parameter_attribute = 'constructible'
            else:
                if bool(lam1 <= 1/2) and bool(lam2 <= 1/2) \
               and bool(-f^2*s_pos*s_neg*lam2 + 2*f*s_pos*s_neg*lam2 + f*s_pos*lam2 + f*s_neg*lam2 - s_pos*s_neg*lam2 - s_pos*lam1 + s_neg*lam1 - s_pos*lam2 - s_neg*lam2 - lam2 <= 0) \
               and bool(-f^2*s_pos*s_neg*lam1 + f*s_pos*lam1 + f*s_neg*lam1 - s_pos*lam2 + s_neg*lam2 - lam1 <= 0):
                    claimed_parameter_attribute = 'extreme'
                else:
                    claimed_parameter_attribute = 'constructible'
            if claimed_parameter_attribute == 'extreme':
                logging.info("Conditions for extremality are satisfied.")
            else:
                logging.info("Conditions for extremality are NOT satisfied.")
        slopes = [s_pos, s_neg, 1/f, s_neg, s_pos, s_neg, s_pos, 1/(f-1), s_pos, s_neg]
        values = [0, v, a/f, 1-a/f, 1-v, 1, 1-w, 1-(1-c)/(1-f), (1-c)/(1-f), w, 0]
        b = f - a
        bb = f - aa
        d = 1 + f - c
        dd = 1 + f - cc
        h = piecewise_function_from_breakpoints_and_values([0, aa, a, b, bb, f, dd, d, c, cc, 1], values, merge=merge, field=field)
        h._claimed_parameter_attribute = claimed_parameter_attribute
        return h

chen_4_slope_reworded = ParametricFamily_chen_4_slope_reworded()

class ParametricFamily_rlm_dpl1_extreme_3a(ParametricFamily):
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = rlm_dpl1_extreme_3a()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    From Richard--Li--Miller [RLM2009].

    For 0 < f < 1/3, by Thm.28, the DPL1 function \phi (whose corresponding h is shown on p.273, Fig.3-lowerleft) is "extreme" (not the usual definition).

    See def.19 for the definition of DPLn, which is a special family of discontinuous piecewise linear functions.
    See Prop.18 and Fig 1 for relation between the DPLn representation `\phi` and the group representation `\pi`, where `\pi(u)` is called `f(u)`, and `f` is called `r_0` throughout this paper.

    All we know from the paper is that `\pi` on p.273, Fig.3-lowerleft is subadditive. However, the extremality is unknown (see discussion after thm.28 p.272).

    Indeed, the function rlm_dpl1_fig3_lowerleft(f) is a facet (and thus extreme) for any 0 < f < 1/3. This can be verified using the covered components and the additivity equations. (Specifically, 2 * \pi(f+) = \pi(2f+) and 2* \pi((1+f) / 2 +) = \pi(f+))

    This is worked out in [KZh2015b, section 2].

    Example: p.273, Fig.3-lowerleft ::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: h = rlm_dpl1_extreme_3a(f=1/4)
        sage: extremality_test(h, False)
        True

    All other 3 functions (corresponding to \phi from the DPL1 family) shown in Fig.3 are proven to be extreme.
    They are covered by ``drlm_3_slope_limit`` and ``drlm_2_slope_limit`` classes:
        * upper-left:  drlm_3_slope_limit(1/3)
        * upper-right: drlm_2_slope_limit(f=3/5, nb_pieces_left=1, nb_pieces_right=1)
        * lower-right: drlm_2_slope_limit(f=3/5, nb_pieces_left=1, nb_pieces_right=2)

    Reference:

    .. [RLM2009] J.-P. P. Richard, Y. Li, and L. A. Miller, Valid inequalities for MIPs and group polyhedra
       from approximate liftings, Mathematical Programming 118 (2009), no. 2, 253-277, doi:10.1007/s10107-007-0190-9

    .. [KZh2015b] M. Koeppe and Y. Zhou, An electronic compendium of extreme functions for the
       Gomory-Johnson infinite group problem, Operations Research Letters, 2015,
       http://dx.doi.org/10.1016/j.orl.2015.06.004
    """

    def _construct_function(self, f=1/4, field=None, conditioncheck=True):
        if not bool(0 < f < 1):
            raise ValueError("Bad parameters. Unable to construct the function.")
        claimed_parameter_attribute = None
        if conditioncheck:
            if bool(f < 1/3):
                pass # is the fig3_lowerleft case
            else:
                pass # is not the fig3_lowerleft case
            claimed_parameter_attribute = 'extreme'
        f = nice_field_values([f], field)[0]
        field = f.parent()
        pieces = [[closed_interval(field(0), f), FastLinearFunction(1/f, 0)], \
                  [open_interval(f, (1 + f)/2), FastLinearFunction(2/(1 + 2*f), 0)], \
                  [singleton_interval((1 + f)/2), FastLinearFunction(field(0), 1/2)], \
                  [open_interval((1 + f)/2, 1), FastLinearFunction(2/(1 + 2*f), -1/(1 + 2*f))], \
                  [singleton_interval(field(1)), FastLinearFunction(field(0), 0)]]
        h = FastPiecewise(pieces)
        h._claimed_parameter_attribute = claimed_parameter_attribute
        return h

rlm_dpl1_extreme_3a = ParametricFamily_rlm_dpl1_extreme_3a()

class ParametricFamily_ll_strong_fractional(ParametricFamily):

    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = ll_strong_fractional()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Letchford--Lodi's strong fractional cut.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = ll_strong_fractional(f=2/3)
        sage: extremality_test(h, False)
        True
        sage: h = ll_strong_fractional(f=2/7)
        sage: minimality_test(h, False)
        False

    Reference:
        [78] Letchford-Lodi (2002) Thm. 2, Fig. 3 (but note this figure shows the wrong function; 
             see ``ll_strong_fractional_bad_figure_3`` and ``ll_strong_fractional_bad_figure_3_corrected``)

        [33] S. Dash and O. Gunluk (2004) Thm. 16

    Remarks:
        Discontinuous, 1-slope;

        For f >= 1/2, this function is facet (extreme), and is identical to
        ``drlm_2_slope_limit(f=f, nb_pieces_left=1, nb_pieces_right=1)``.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: f=2/3
        sage: l = ll_strong_fractional(f)
        sage: d = drlm_2_slope_limit(f=f, nb_pieces_left=1, nb_pieces_right=ceil(1/f)-1)
        sage: dg = automorphism(dg_2_step_mir_limit(f=1-f, d=ceil(1/f)-1))
        sage: show(plot(l, color='red', legend_label='ll_strong_fractional')) # not tested
        sage: show(plot(d, color='blue', legend_label='drlm_2_slope_limit')) # not tested
        sage: show(plot(dg, color='green', legend_label='automorphism(dg_2_step_mir_limit)')) # not tested
        sage: l == d == dg
        True

    Remarks:
        The function is NOT minimal for 0 < f < 1/2.  It equals
        ``drlm_2_slope_limit(f=f, nb_pieces_left=1, nb_pieces_right=ceil(1/f)-1)``,
        except for limits at breakpoints.

    EXAMPLES::

        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: f=1/3
        sage: l = ll_strong_fractional(f)
        sage: d = drlm_2_slope_limit(f=f, nb_pieces_left=1, nb_pieces_right=ceil(1/f)-1)
        sage: dg = automorphism(dg_2_step_mir_limit(f=1-f, d=ceil(1/f)-1))
        sage: show(plot(l, color='red', legend_label='ll_strong_fractional')) # not tested
        sage: show(plot(d, color='blue', legend_label='drlm_2_slope_limit')) # not tested
        sage: show(plot(dg, color='green', legend_label='automorphism(dg_2_step_mir_limit)')) # not tested
        sage: d == dg
        True
        sage: l == d
        False
    """

    def _claimed_parameter_attribute(self, f, **kwargs):
        if not bool(0 < f < 1):
            return 'not_constructible'
        if not bool(1/2 <= f < 1):
            return 'constructible'
        else:
            return 'extreme'

    def _construct_function(self, f=2/3, field=None, conditioncheck=True):
        if not bool(0 < f < 1):
            raise ValueError("Bad parameters. Unable to construct the function.")
        f = nice_field_values([f], field)[0]
        field = f.parent()
        k = ceil(1/f) -1
        pieces = [[closed_interval(0,f), FastLinearFunction(1/f, 0)]]
        for p in range(k-1):
            pieces.append([left_open_interval(f + (1 - f)* p / k, f + (1 - f)*(p + 1)/k), FastLinearFunction(1/f, -(p + 1)/f/(k + 1))])
        p = k - 1
        pieces.append([open_interval(f + (1 - f)* p / k, f + (1 - f)*(p + 1)/k), FastLinearFunction(1/f, -(p + 1)/f/(k + 1))])
        pieces.append([singleton_interval(1), FastLinearFunction(0, 0)])
        return PiecewiseLinearFunction_1d(pieces)

ll_strong_fractional = ParametricFamily_ll_strong_fractional()

class ParametricFamily_bcdsp_arbitrary_slope(ParametricFamily):

    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = bcdsp_arbitrary_slope()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    A family of extreme functions with an arbitrary number `k` of slopes. (k >= 2)

    Function is known to be extreme under the condition:
        0 < f <= 1/2.

    Tests show that the function is also extreme when f <= 4/5.

    Examples::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = bcdsp_arbitrary_slope(f=1/2, k=2)
        sage: h == gmic(f=1/2)
        True
        sage: h = bcdsp_arbitrary_slope(f=1/2, k=3)
        sage: h == gj_forward_3_slope(f=1/2, lambda_1=1/2, lambda_2=1/4)
        True
        sage: h = bcdsp_arbitrary_slope(f=1/2, k=4)
        sage: number_of_slopes(h)
        4
        sage: extremality_test(h)
        True
        sage: h = bcdsp_arbitrary_slope(f=4/5, k=10)
        sage: number_of_slopes(h)
        10
        sage: extremality_test(h)
        True

    Reference:
         [arbitrary_num_slopes] A. Basu, M. Conforti, M. Di Summa, and J. Paat, Extreme Functions with an Arbitrary Number of Slopes, 2015, http://www.ams.jhu.edu/~abasu9/papers/infinite-slopes.pdf, to appear in Proceedings of IPCO 2016.
    """

    def _construct_function(self, f=1/2, k=4, field=None, conditioncheck=True):
        if not bool(0 < f < 1) or k not in ZZ or k < 2:
            raise ValueError("Bad parameters. Unable to construct the function.")
        claimed_parameter_attribute = None
        if conditioncheck:
            if not bool(0 < f <= 1/2):
                logging.info("Conditions for extremality are NOT satisfied.")
                claimed_parameter_attribute = 'constructible'
            else:
                logging.info("Conditions for extremality are satisfied.")
                claimed_parameter_attribute = 'extreme'
        f = nice_field_values([f], field)[0]
        field = f.parent()
        bkpts = [field(0)]
        slopes = []
        for i in range(k-2, 0, -1):
            bkpts += [f/(8^i), 2*f/(8^i)]
            slopes += [(2^i - f)/f/(1-f), 1/(f-1)]
        bkpts = bkpts + [f - x for x in bkpts[::-1]] + [field(1)]
        slopes = slopes + [1/f] + slopes[::-1] + [1/(f-1)]
        h = piecewise_function_from_breakpoints_and_slopes(bkpts, slopes, field=field)
        h._claimed_parameter_attribute = claimed_parameter_attribute
        return h

bcdsp_arbitrary_slope = ParametricFamily_bcdsp_arbitrary_slope()

extreme_function_with_world_record_number_of_slopes = bcdsp_arbitrary_slope

class bcds_discontinuous_everywhere:
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = bcds_discontinuous_everywhere()
        g = h.plot(show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, ticks=[[QQ('1/2'), 1],[1]], tick_formatter=[[r'$f=\frac{1}{2}$', "$1$"], ["$1$"]])
        sphinx_plot(g)

    An extreme function whose graph is dense in R \times [0,1].

    Reference:
        [bcds_discontinous_everywhere] Amitabh Basu, Michele Conforti, Marco Di Summa, An extreme function which is nonnegative and discontinuous everywhere, arXiv:1802.01499 [math.OC]
    """
    def __init__(self):
        self._f = 1/2

    def __call__(self, x):
        """
        Examples::

            sage: from cutgeneratingfunctionology.igp import bcds_discontinuous_everywhere, delta_pi
            sage: h = bcds_discontinuous_everywhere()
            sage: delta_pi(h, 1/5+sqrt(3), 3/7+sqrt(1/3)) >= 0
            True
            sage: delta_pi(h, -13/9+sqrt(17), 3/7-3*sqrt(17))>=0
            True
            sage: delta_pi(h, 19/8, 3/7-3*sqrt(101))>=0
            True
            sage: delta_pi(h, 19/8+sqrt(101/9), 3/7-1/3*sqrt(101))>=0
            True
        """
        try:
            xx=AA(x)
        except TypeError:
            raise NotImplementedError("Not implemented for non-algebraic numbers.")
        p=xx.minpoly()
        if p.degree()>2:
            raise NotImplementedError("Not implemented for algebraic numbers with degree greater than 2.") 
        elif p.degree()==2:
            s=-p.coefficients(sparse=False)[1]/2
        else:
            s=xx
        ss=fractional(s)
        if ss<=1/2:
            return 2*ss
        else:
            return 2-2*ss
               
    def plot(self, xmin=0, xmax=1, ymin=0, ymax=1, color=None, rgbcolor=None,
             aspect_ratio='automatic', thickness=None, **kwds):
        ymin = max(0, ymin)
        ymax = min(1, ymax)
        color = color or rgbcolor or 'blue'
        if xmax >= xmin and ymax >= ymin:
            return polygon(((xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)),
                           fill=True, alpha=0.2, thickness=0, color=color,
                           aspect_ratio=aspect_ratio,
                           **kwds)
        else:
            return Graphics()

class kzh_extreme_and_weak_facet_but_not_facet:

    """
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_extreme_and_weak_facet_but_not_facet()
        g = h.plot(show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h.pi))
        sphinx_plot(g)

    An extreme function and weak facet that is not a facet.

    This factory class generates a fresh function from an infinite family.
    The function depends on the sequence of its evaluations.

    TESTS::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = kzh_extreme_and_weak_facet_but_not_facet()
        sage: hm = kzh_minimal_has_only_crazy_perturbation_1()+1/1000000*kzh_minimal_has_only_crazy_perturbation_1_perturbation()
        sage: minimality_test_randomized(h, kzh_minimal_has_only_crazy_perturbation_1(), hm,
        ....:                            limits=False, extra_cosets=[sqrt(3), -sqrt(3)],
        ....:                            max_iterations=10)
        True

    Reference:
        Matthias Koeppe, Yuan Zhou. On the notions of facets, weak facets, and extreme functions of the Gomory-Johnson infinite group problem.
    """
    def __init__(self, pi=None):
        self.Cplus = set()
        if pi is None:
            pi = kzh_minimal_has_only_crazy_perturbation_1()
        self.pi = pi

    def __call__(self, x):
        """
        Examples::

            sage: from cutgeneratingfunctionology.igp import kzh_extreme_and_weak_facet_but_not_facet, delta_pi
            sage: h = kzh_extreme_and_weak_facet_but_not_facet()
            sage: bool(delta_pi(h, 1/5+sqrt(3), 3/7+sqrt(1/3)) >= 0)
            True
            sage: bool(delta_pi(h, -13/9, 3/7)>=0)
            True
            sage: bool(delta_pi(h, 19/8, 3/7-3*sqrt(101))>=0)
            True
            sage: bool(delta_pi(h, 19/8+sqrt(101/9), 3/7-1/3*sqrt(101))>=0)
            True
        """

        x = fractional(x)
        f = self.pi._f
        l = self.pi.ucl
        u = self.pi.ucr
        t1 = self.pi.t1
        t2 = self.pi.t2
        s = self.pi.s
        pi=self.pi
        if x<l or u<x<f-u or x>f-l:
            return pi(x)
        if l<x<u and self.is_in_C(x):
            return pi(x)
        if f-u<x<f-l and self.is_in_C(f-x):
            return pi(x)
        if l<x<u and self.is_in_Cplus(x):
            return pi(x)+s
        if f-u<x<f-l and self.is_in_Cplus(f-x):
            return pi(x)+s
        return pi(x)-s

    def plot(self, *args, **kwds):
        return plot_with_colored_slopes(self.pi, **kwds)

    def is_in_T(self, x):
        try:
            t1,t2,xx=nice_field_values([self.pi.t1, self.pi.t2, x])
            coe=xx.list()
            # make sure x is properly embedded.
            if not len(coe)==2:
                return False
            xxx=nice_field_values([coe[0]+coe[1]*sqrt(2)])[0]
            if not xxx==xx:
                return False
            c1=coe[0]/self.pi.t2
            c2=coe[1]/(self.pi.t1/sqrt(2))
            if c1 in ZZ and c2 in ZZ:
                return True
            else: 
                return False
        except (TypeError,ValueError):
            return False
 
    def is_in_C(self, x):  
        f=self.pi._f
        l=self.pi.ucl
        u=self.pi.ucr
        t1=self.pi.t1
        t2=self.pi.t2
        if self.is_in_T(x-(l+u)/2) or self.is_in_T(x-(l+u-t1)/2) or self.is_in_T(x-(l+u-t2)/2):
            return True
        else:
            return False
                      
    def is_in_Cplus(self,x): 
        f=self.pi._f
        l=self.pi.ucl
        u=self.pi.ucr
        t1=self.pi.t1
        t2=self.pi.t2
        phi_x=l+u-x
        if self.is_in_T(x-phi_x):
            self.Cplus.add(x)
            return True
        for v in self.Cplus:
            if self.is_in_T(x-v):
                return True
            if self.is_in_T(phi_x-v):
                return False
        self.Cplus.add(x)
        return True

class ParametricFamily_kzh_3_slope_param_extreme_1(ParametricFamily):

    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_3_slope_param_extreme_1()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    New extreme function discovered by computer based search followed by parametric search.
    It has 3 slopes in the general case.

    Parameters: real numbers `f`, `a`, `b`, where `a` is the length of the first interval right to `f` and `b` is the length of interval centered at `(1+f)/2`.

    Function is known to be extreme under the conditions:
        0 <= a and  0 <= b <= f and 3*f+4*a-b-1 <= 0

    Examples::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = kzh_3_slope_param_extreme_1(f=6/19, a=1/19, b=5/19)
        sage: extremality_test(h)
        True
        sage: h = kzh_3_slope_param_extreme_1(f=1/3, a=1/18, b=1/4)
        sage: extremality_test(h)
        True
    """

    def _construct_function(self, f=6/19, a=1/19, b=5/19, field=None, conditioncheck=True):
        if not bool(0 < f < f+a < (1+f-b)/2 < (1+f+b)/2 < 1-a < 1):
            raise ValueError("Bad parameters. Unable to construct the function.")
        claimed_parameter_attribute = None
        if conditioncheck:
            if not bool(0 <= a and  0 <= b <= f and 3*f+4*a-b-1 <= 0):
                logging.info("Conditions for extremality are NOT satisfied.")
                claimed_parameter_attribute = 'constructible'
            else:
                logging.info("Conditions for extremality are satisfied.")
                claimed_parameter_attribute = 'extreme'
        v = (f*f+f*a-3*f*b-3*a*b+b)/(f*f+f-3*f*b)
        bkpts = [0, f, f+a, (1+f-b)/2, (1+f+b)/2, 1-a, 1]
        values = [0, 1, v, (f-b)/2/f, (f+b)/2/f, 1-v, 0]
        h = piecewise_function_from_breakpoints_and_values(bkpts, values, field=field)
        h._claimed_parameter_attribute = claimed_parameter_attribute
        return h

kzh_3_slope_param_extreme_1 = ParametricFamily_kzh_3_slope_param_extreme_1()

class ParametricFamily_kzh_3_slope_param_extreme_2(ParametricFamily):

    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_3_slope_param_extreme_2()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    New extreme function discovered by computer based search followed by parametric search.
    The function looks like ``gj_forward_3_slope`` + ``drlm_backward_3_slope``.

    Parameters: real numbers `f`, `a`, `b`, where `a` is the length of interval centered at `f/2`, and `b` is the length of interval centered at `(1+f)/2`.

    Function is known to be extreme under the (sufficient) conditions:
        0 < a < f < 1; 2*b - a <= f <= a + b and f <= (1+a-b)/2

    Examples::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = kzh_3_slope_param_extreme_2(f=5/9, a=3/9, b=2/9)
        sage: extremality_test(h)
        True
        sage: h = kzh_3_slope_param_extreme_2(f=4/9, a=2/9, b=3/9)
        sage: extremality_test(h) # Claimed conditions for extremality are NOT satisfied.
        True
        sage: h = kzh_3_slope_param_extreme_2(5/9, 2/5, 7/20)
        sage: extremality_test(h)
        False
    """

    def _construct_function(self, f=5/9, a=3/9, b=2/9, field=None, conditioncheck=True):
        if not bool(0 < a < f < 1 and 0 < b < 1-f):
            raise ValueError("Bad parameters. Unable to construct the function.")
        claimed_parameter_attribute = None
        if conditioncheck:
            if not bool(2*b - a <= f <= a + b and f <= (1+a-b)/2):
                logging.info("Conditions for extremality are NOT satisfied.")
                claimed_parameter_attribute = 'constructible'
            else:
                logging.info("Conditions for extremality are satisfied.")
                claimed_parameter_attribute = 'extreme'
        v = (f*(f-a+b-2)-a*b+2*a)/(f+b-1)/f/4;
        bkpts = [0, (f-a)/4, (f-a)/2, (f+a)/2, f-(f-a)/4, f, \
                 (1+f-b)/2, (1+f+b)/2, 1]
        values = [0, v, (f-a)/f/2, (f+a)/f/2, 1-v, 1, (f-b)/f/2, (f+b)/f/2, 0]
        h = piecewise_function_from_breakpoints_and_values(bkpts, values, field=field)
        h._claimed_parameter_attribute = claimed_parameter_attribute
        return h

kzh_3_slope_param_extreme_2 = ParametricFamily_kzh_3_slope_param_extreme_2()

class ParametricFamily_kzh_4_slope_param_extreme_1(ParametricFamily):
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_4_slope_param_extreme_1()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    New extreme function discovered by computer based search followed by parametric search.
    The function looks like ``gj_forward_3_slope`` + ``kzh_3_slope_param_extreme_1``.

    Parameters: real numbers `f`, `a`, `b`, where `a` is the length of interval centered at `f/2`, and `b` is the first bkpt, which is also the length of interval centered at `(1+f)/2`.

    Function is known to be extreme under the (sufficient) conditions:
        a > max{3*f-2, f/2}; b > 0; a + 3*b > 2*f - 1; 2*a + 3*b < 3*f - 1.

    Examples::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = kzh_4_slope_param_extreme_1(f=13/18, a=7/18, b=1/18)
        sage: extremality_test(h)
        True
        sage: h = kzh_4_slope_param_extreme_1(f=13/18, a=14/37, b=1/19)
        sage: extremality_test(h)
        True
    """

    def _construct_function(self, f=13/18, a=7/18, b=1/18, field=None, conditioncheck=True):
        w = (f+a)/f/4;
        v = (3*b+1)/f/3;
        c = (1/4*f - 1/4*a - 1/2*b)
        if not bool(0 < b < (f-a)/2 < f < (1+f-b)/2-c < (1+f-b)/2 < 1):
            raise ValueError("Bad parameters. Unable to construct the function.")
        claimed_parameter_attribute = None
        if conditioncheck:
            if not bool(b > 0 and -3*f + 2*a + 3*b + 1 < 0 and  f - 2*a < 0 \
                        and  2*f - a - 3*b - 1 < 0 and  3*f - a - 2 < 0):
                logging.info("Conditions for extremality are NOT satisfied.")
                claimed_parameter_attribute = 'constructible'
            else:
                logging.info("Conditions for extremality are satisfied.")
                claimed_parameter_attribute = 'extreme'
        bkpts = [0, b, (f-a)/2, (f+a)/2, f-b, f,\
                 (1+f-b)/2-c, (1+f-b)/2, (1+f+b)/2, (1+f+b)/2+c, 1]
        values = [0, v, (f-a)/f/2, (f+a)/f/2, 1-v, 1, \
                  w, (1-v)/2, (1+v)/2, 1-w, 0]
        h = piecewise_function_from_breakpoints_and_values(bkpts, values, field=field)
        h._claimed_parameter_attribute = claimed_parameter_attribute
        return h

kzh_4_slope_param_extreme_1 = ParametricFamily_kzh_4_slope_param_extreme_1()
