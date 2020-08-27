"""
Piecewise linear functions of one real variable
"""
from __future__ import division, print_function, absolute_import

from bisect import bisect_left

from sage.structure.element import Element, ModuleElement
from sage.structure.richcmp import richcmp, op_NE, op_EQ
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from .intervals import *
from . import FastLinearFunction
from cutgeneratingfunctionology.spam.parametric_real_field_element import is_parametric_element

class PiecewiseLinearFunction_1d (ModuleElement):
    r"""
    Returns a piecewise function from a list of (interval, function)
    pairs.

    Uses binary search to allow for faster function evaluations
    than the standard class ``PiecewisePolynomial``.

    The parent, an instance of ``PiecewiseLinearFunctionsSpace``, has
    the concept of a ``domain``, a ``codomain``, and a ``base_ring``.

    For this element class, ``domain`` must be an exact real field.
    All breakpoints lie in ``domain``, but the function can be called
    with any element of a real field that has a common pushout.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = FastPiecewise([[(3/10, 15/40), FastLinearFunction(1, 0)], [(13/40, 14/40), FastLinearFunction(1, 0)]], merge=True)
        sage: len(h.intervals())
        1
        sage: h.intervals()[0][0], h.intervals()[0][1]
        (3/10, 3/8)
        sage: h = FastPiecewise([[(3/10, 15/40), FastLinearFunction(1, 0)],
        ....:                    [(13/40, 14/40), FastLinearFunction(1, 0)],
        ....:                    [(17,18), FastLinearFunction(77,78)]], merge=True)
        sage: len(h.intervals())
        2
        sage: h.intervals()[0][0], h.intervals()[0][1]
        (3/10, 3/8)

    By default, function values are cached but this can be turned off to reduce memory use::

        sage: h = FastPiecewise([[(3/10, 15/40), FastLinearFunction(1, 0)]])
        sage: hasattr(h, '_call_cache')
        True
        sage: FastPiecewise.cache = False
        sage: h = FastPiecewise([[(3/10, 15/40), FastLinearFunction(1, 0)]])
        sage: hasattr(h, '_call_cache')
        False
        sage: FastPiecewise.cache = True
        sage: h = FastPiecewise([[(3/10, 15/40), FastLinearFunction(1, 0)]])
        sage: hasattr(h, '_call_cache')
        True
        sage: h = FastPiecewise([[(3/10, 15/40), FastLinearFunction(1, 0)]], cache=False)
        sage: hasattr(h, '_call_cache')
        False

    """

    cache = True

    def __init__(self, list_of_pairs, var=None, periodic_extension=True, merge=True, cache=None, parent=None):
        r"""
        ``list_of_pairs`` is a list of pairs (I, fcn), where
        fcn is a Sage function (such as a polynomial over RR, or functions
        using the lambda notation), and I is an interval such as I = (1,3).

        If the optional ``var`` is specified, this variable will be used for printing purposes.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: f1(x) = 1
            sage: f2(x) = 1 - x
            sage: f = FastPiecewise([[(0,1),f1],[(1,2),f2]])
            sage: f.list()
            [[(0, 1), x |--> 1], [(1, 2), x |--> -x + 1]]
            sage: f.length()
            2
        """
        # Sort intervals according to their left endpoints; In case of equality, place single point before interval.
        list_of_pairs = sorted(list_of_pairs, key = lambda i_f: coho_interval_left_endpoint_with_epsilon(i_f[0]))
        if merge:
            merged_list_of_pairs = []
            intervals_to_scan = []
            singleton = None
            common_f = None
            for (i, f) in list_of_pairs:
                if len(i) == 1:
                    i = singleton_interval(i[0])            # upgrade to coho interval
                if common_f == f:
                    intervals_to_scan.append(i)
                    singleton = None
                elif common_f is not None and singleton is not None and common_f(singleton) == f(singleton):
                    intervals_to_scan.append(i)
                    singleton = None
                    common_f = f
                elif i[0] == i[1] and common_f is not None and common_f(i[0]) == f(i[0]):
                    intervals_to_scan.append(i)
                else:
                    merged_intervals = union_of_coho_intervals_minus_union_of_coho_intervals([[interval] for interval in intervals_to_scan], [],
                                                                                             old_fashioned_closed_intervals=True)
                    for merged_interval in merged_intervals:
                        merged_list_of_pairs.append((merged_interval, common_f))
                    intervals_to_scan = [i]
                    if i[0] == i[1]:
                        singleton = i[0]
                    else:
                        singleton = None
                    common_f = f
            merged_intervals = union_of_coho_intervals_minus_union_of_coho_intervals([[interval] for interval in intervals_to_scan], [],
                                                                                     old_fashioned_closed_intervals=True)
            for merged_interval in merged_intervals:
                merged_list_of_pairs.append((merged_interval, common_f))
            list_of_pairs = merged_list_of_pairs

        self._length = len(list_of_pairs)
        if parent is None:
            from sage.structure.element import get_coercion_model
            from sage.rings.rational_field import QQ
            cm = get_coercion_model()
            enclosing_domain = cm.common_parent(QQ(1), *(a for interval, func in list_of_pairs
                                                         for a in (interval[0], interval[1])))
            from sage.categories.pushout import pushout
            def codomain_from_func(func):
                if isinstance(func, FastLinearFunction):
                    return cm.common_parent(func._slope, func._intercept)
                return QQ
            codomains = [codomain_from_func(func) for interval, func in list_of_pairs]
            if codomains:
                codomain = pushout(enclosing_domain, cm.common_parent(*codomains))
            else:
                codomain = enclosing_domain
            parent = PiecewiseLinearFunctionsSpace(enclosing_domain, codomain)
        else:
            enclosing_domain = parent.domain()
            codomain = parent.codomain()

        ModuleElement.__init__(self, parent)

        intervals = self._intervals = [convert_interval(x[0], enclosing_domain)
                                       for x in list_of_pairs]
        functions = self._functions = [x[1] for x in list_of_pairs]
        self._list = [[self._intervals[i], self._functions[i]] for i in range(self._length)]

        # end_points are distinct.
        end_points = []
        # ith_at_end_points records in which interval the end_point first appears as a left_end or right_end.
        ith_at_end_points = []
        # record the value at each end_point, value=None if end_point is not in the domain.
        values_at_end_points = []
        # record function values at [x, x+, x-] for each endpoint x.
        limits_at_end_points = []
        left_limit = None
        for i in range(len(intervals)):
            left_value = None
            if len(intervals[i]) <= 2 or intervals[i].left_closed:
                left_value = functions[i](intervals[i][0])
            if intervals[i][0] != intervals[i][1]:
                right_limit = functions[i](intervals[i][0])
            else:
                right_limit = None
            if (end_points == []) or (end_points[-1] != intervals[i][0]):
                end_points.append(intervals[i][0])
                ith_at_end_points.append(i)
                values_at_end_points.append(left_value)
                if limits_at_end_points != []:
                    limits_at_end_points[-1][1]= None
                limits_at_end_points.append([left_value, right_limit, None])
            else:
                if left_value is not None:
                    values_at_end_points[-1] = left_value
                    limits_at_end_points[-1][0] = left_value
                limits_at_end_points[-1][1] = right_limit
            right_value = None
            if len(intervals[i]) <= 2 or intervals[i].right_closed:
                right_value = functions[i](intervals[i][1])
            if intervals[i][0] != intervals[i][1]:
                left_limit = functions[i](intervals[i][1])
                end_points.append(intervals[i][1])
                ith_at_end_points.append(i)
                values_at_end_points.append(right_value)
                limits_at_end_points.append([right_value, None, left_limit])
            elif right_value is not None:
                values_at_end_points[-1] = right_value
        if periodic_extension and limits_at_end_points != []:
            #if values_at_end_points[0] != values_at_end_points[-1]:
            #    logging.warning("Function is actually not periodically extendable.")
            #    periodic_extension = False
            #else:
                limits_at_end_points[0][-1] = limits_at_end_points[-1][-1]
                limits_at_end_points[-1][1] = limits_at_end_points[0][1]
        self._end_points = end_points
        self._ith_at_end_points = ith_at_end_points
        self._values_at_end_points = values_at_end_points
        self._limits_at_end_points = limits_at_end_points
        self._periodic_extension = periodic_extension
        if not any(is_parametric_element(x) for x in end_points):
            # We do not wish to keep a call cache when the breakpoints are
            # parametric, because that may create `mysterious' proof cells.
            if cache is None:
                cache = self.__class__.cache
            if cache:
                self._call_cache = dict()
        is_continuous = True
        if len(end_points) == 1 and end_points[0] is None:
            is_continuous = False
        elif len(end_points)>= 2:
            [m0, r0, l0] = limits_at_end_points[0]
            [m1, r1, l1] = limits_at_end_points[-1]
            if m0 is None or r0 is None or  m0 != r0 or l1 is None or m1 is None or l1 != m1:
                is_continuous = False
            else:
                for i in range(1, len(end_points)-1):
                    [m, r, l] = limits_at_end_points[i]
                    if l is None or m is None or r is None or not(l == m == r):
                        is_continuous = False
                        break
        self._is_continuous = is_continuous
        self._is_two_sided_discontinuous = not ( is_continuous or \
                                                 limits_at_end_points[0][0] == limits_at_end_points[0][1] or
                                                 limits_at_end_points[-1][-1] == limits_at_end_points[-1][0] )

    # The following makes this class hashable and thus enables caching
    # of the above functions; but we must promise not to modify the
    # contents of the instance.
    def __hash__(self):
        return id(self)

    def _richcmp_(self, other, op):
        r"""
        Comparison.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.WARN) # Suppress output in automatic tests.
            sage: f = FastPiecewise([[open_interval(1,3), FastLinearFunction(0,3)]])
            sage: g =  FastPiecewise([[open_interval(1,2), FastLinearFunction(0,3)], [right_open_interval(2,3), FastLinearFunction(0,3)]], merge=False)
            sage: f == g
            True
            sage: g == f
            True
            sage: bkpt = [0, 1/8, 3/8, 1/2, 5/8, 7/8, 1]
            sage: limits = [(0, 0, 1/2), (1/4, 1/4, 3/4), (3/4, 1/4, 3/4), (1, 1/2, 1), (3/4, 3/4, 3/4), (1/4, 1/4, 1/4), (0, 0, 1/2)]
            sage: h = piecewise_function_from_breakpoints_and_limits(bkpt, limits)
            sage: h == hildebrand_discont_3_slope_1()
            True
            sage: f = FastPiecewise([[open_interval(1,2), FastLinearFunction(1,3)], [right_open_interval(3,4), FastLinearFunction(2,3)]])
            sage: g = FastPiecewise([[open_interval(1,2), FastLinearFunction(1,3)], [right_open_interval(3,4), FastLinearFunction(2,3)]])
            sage: f == g
            True
            sage: g == f
            True

        Test for discrete_function_from_points_and_values::

            sage: h1 = discrete_function_from_points_and_values([i/11 for i in range(12)],[0, 5/6, 1/5, 2/3, 2/5, 1/2, 3/5, 1/3, 4/5, 1/6, 1, 0])
            sage: h2 = discrete_function_from_points_and_values([i/11 for i in range(12)],[0, 5/6, 4/9, 2/3, 5/18, 1/2, 13/18, 1/3, 5/9, 1/6, 1, 0])
            sage: h1 == h2
            False
        """
        if not isinstance(other, FastPiecewise):
            return NotImplemented
        def is_equal():
            if self._periodic_extension != other._periodic_extension:
                return False
            domain = union_of_coho_intervals_minus_union_of_coho_intervals([self.intervals()], [], old_fashioned_closed_intervals=True)
            if union_of_coho_intervals_minus_union_of_coho_intervals([other.intervals()],[], old_fashioned_closed_intervals=True) != domain:
                return False
            difference = self - other
            return (difference.intervals() == domain) and all(fn == FastLinearFunction(0,0) for fn in difference.functions())
        if op == op_EQ:
            return is_equal()
        elif op == op_NE:
            return not is_equal()
        else:
            # FIXME: Implement partial order of pointwise comparison
            return NotImplemented

    def is_continuous(self):
        r"""
        return if the function is continuous
        """
        return self._is_continuous

    def is_two_sided_discontinuous(self):
        r"""
        return if the function is discontinuous at 0+ and at 1-.
        """
        return self._is_two_sided_discontinuous

    def is_discrete(self):
        r"""
        Return if the function is discrete, i.e., all pieces are singletons
        """
        return all(interval_length(interval) == 0 for interval in self.intervals())

    def end_points(self):
        r"""
        Returns a list of all interval endpoints for this function.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: f1(x) = 1
            sage: f2(x) = 2
            sage: f3(x) = 1-x
            sage: f4(x) = x^2-5
            sage: f = FastPiecewise([[open_interval(0,1),f1],[singleton_interval(1),f2],[open_interval(1,2),f3],[(2,3),f4]])
            sage: f.end_points()
            [0, 1, 2, 3]
            sage: f = FastPiecewise([[open_interval(0,1),f1],[open_interval(2,3),f3]])
            sage: f.end_points()
            [0, 1, 2, 3]
        """
        return self._end_points

    def values_at_end_points(self):
        r"""
        Returns a list of function values at all endpoints for this function.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = 4
            sage: f5(x) = sin(2*x)
            sage: f6(x) = x-3
            sage: f7(x) = 7
            sage: f = FastPiecewise([[right_open_interval(0,1),f1],
            ....:                    [right_open_interval(1,2),f2],
            ....:                    [open_interval(2,3),f3],
            ....:                    [singleton_interval(3),f4],
            ....:                    [left_open_interval(3,6),f5],
            ....:                    [open_interval(6,7),f6],
            ....:                    [(9,10),f7]])
            sage: f.values_at_end_points()
            [1, 0, None, 4, sin(12), None, 7, 7]
        """
        return self._values_at_end_points

    def limits_at_end_points(self):
        r"""
        Returns a list of 3-tuples [function value, right limit, left limit] at all endpoints for this function.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = 4
            sage: f5(x) = sin(2*x)
            sage: f6(x) = x-3
            sage: f7(x) = 7
            sage: f = FastPiecewise([[right_open_interval(0,1),f1],
            ....:                    [right_open_interval(1,2),f2],
            ....:                    [open_interval(2,3),f3],
            ....:                    [singleton_interval(3),f4],
            ....:                    [left_open_interval(3,6),f5],
            ....:                    [open_interval(6,7),f6],
            ....:                    [(9,10),f7]], periodic_extension= False)
            sage: f.limits_at_end_points()
            [[1, 1, None], [0, 0, 1], [None, e^2, -1], [4, sin(6), e^3], [sin(12), 3, sin(12)], [None, None, 4], [7, 7, None], [7, None, 7]]
        """
        return self._limits_at_end_points

    def functions(self):
        """
        Returns the list of functions (the "pieces").

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp.fast_piecewise import *
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = FastPiecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            sage: f.functions()
            [x |--> 1, x |--> -x + 1, x |--> e^x, x |--> sin(2*x)]
        """
        return self._functions

    def intervals(self):
        """
        A piecewise non-polynomial example.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp.fast_piecewise import *
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = FastPiecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            sage: f.intervals()
            [(0, 1), (1, 2), (2, 3), (3, 10)]
        """
        return self._intervals

    def list(self):
        """
        Returns the pieces of this function as a list of functions.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp.fast_piecewise import *
            sage: f1(x) = 1
            sage: f2(x) = 1 - x
            sage: f = FastPiecewise([[(0,1),f1],[(1,2),f2]])
            sage: f.list()
            [[(0, 1), x |--> 1], [(1, 2), x |--> -x + 1]]
        """
        return self._list

    def length(self):
        """
        Returns the number of pieces of this function.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp.fast_piecewise import *
            sage: f1(x) = 1
            sage: f2(x) = 1 - x
            sage: f = FastPiecewise([[(0,1),f1],[(1,2),f2]])
            sage: f.length()
            2
        """
        return self._length



    def which_function(self, x0):
        r"""
        Returns the function piece used to evaluate self at `x_0`.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = FastPiecewise([[(0,1),f1],
            ....:                    [(1,2),f2],
            ....:                    [(2,3),f3],
            ....:                    [(3,10),f4]])
            sage: f.which_function(0.5) is f1
            True
            sage: f.which_function(1) in [f1, f2]
            True
            sage: f.which_function(5/2) is f3
            True
            sage: f.which_function(3) in [f3, f4]
            True
            sage: f.which_function(-1)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point -1, outside of domain.
            sage: f.which_function(11)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 11, outside of domain.
            sage: f = FastPiecewise([[right_open_interval(0,1),f1],
            ....:                    [right_open_interval(1,2),f2],
            ....:                    [right_open_interval(2,3),f3],
            ....:                    [closed_interval(3,10),f4]])
            sage: f.which_function(0.5) is f1
            True
            sage: f.which_function(1) is f2
            True
            sage: f.which_function(5/2) is f3
            True
            sage: f.which_function(3) is f4
            True
            sage: f = FastPiecewise([[open_interval(0,1),f1],
            ....:                    [right_open_interval(2,3),f3]])
            sage: f.which_function(0)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 0, outside of domain.
            sage: f.which_function(0.5) is f1
            True
            sage: f.which_function(1)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 1, outside of domain.
            sage: f.which_function(3/2)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 3/2, outside of domain.
            sage: f.which_function(2) is f3
            True
            sage: f.which_function(5/2) is f3
            True
            sage: f.which_function(3)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 3, outside of domain.

        Test that it works correctly for functions set up with merge=False::

            sage: f = FastPiecewise([[right_open_interval(0,1), FastLinearFunction(0, 0)],
            ....:                    [right_open_interval(1,2), FastLinearFunction(1, -1)]],
            ....:                   merge=False)
            sage: f.which_function(1)
            <FastLinearFunction x - 1>

        """
        return self.which_pair(x0)[1]

    def which_pair(self, x0):
        endpts = self.end_points()
        ith = self._ith_at_end_points
        i = bisect_left(endpts, x0)
        if i >= len(endpts):
            raise ValueError("Value not defined at point %s, outside of domain." % x0)
        if x0 == endpts[i]:
            if self._values_at_end_points[i] is not None:
                if is_pt_in_interval(self.intervals()[ith[i]], x0):
                    return self.intervals()[ith[i]], self.functions()[ith[i]]
                else:
                    return self.intervals()[ith[i]+1], self.functions()[ith[i]+1]
            else:
                raise ValueError("Value not defined at point %s, outside of domain." % x0)
        if i == 0:
            raise ValueError("Value not defined at point %s, outside of domain." % x0)
        if is_pt_in_interval(self._intervals[ith[i]],x0):
            return self.intervals()[ith[i]], self.functions()[ith[i]]
        raise ValueError("Value not defined at point %s, outside of domain." % x0)

    def __call__(self,x0):
        r"""
        Evaluates self at `x_0`.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: f1(x) = 1
            sage: f2(x) = 1-x

            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = FastPiecewise([[(0,1),f1],
            ....:                    [(1,2),f2],
            ....:                    [(2,3),f3],
            ....:                    [(3,10),f4]])
            sage: f(0.5)
            1
            sage: f(1)
            0
            sage: f(5/2)
            e^(5/2)
            sage: f(3)
            sin(6)
            sage: f(-1)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point -1, outside of domain.
            sage: f(11)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 11, outside of domain.
            sage: f = FastPiecewise([[right_open_interval(0,1),f1],
            ....:                    [right_open_interval(1,2),f2],
            ....:                    [right_open_interval(2,3),f3],
            ....:                    [closed_interval(3,10),f4]])
            sage: f(0.5)
            1
            sage: f(1)
            0
            sage: f(5/2)
            e^(5/2)
            sage: f(3)
            sin(6)
            sage: f = FastPiecewise([[open_interval(0,1),f1],
            ....:                    [right_open_interval(2,3),f3]])
            sage: f(0)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 0, outside of domain.
            sage: f(0.5)
            1
            sage: f(1)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 1, outside of domain.
            sage: f(3/2)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 3/2, outside of domain.
            sage: f(2)
            e^2
            sage: f(5/2)
            e^(5/2)
            sage: f(3)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 3, outside of domain.

        Big cell parametrics::

            sage: K.<f> = ParametricRealField([4/5], big_cells=True)
            sage: h = gmic(f)
            sage: h(4/5)
            (4/5/f)~
            sage: with K.frozen():
            ....:     f == 4/5
            Traceback (most recent call last):
            ...
            ParametricRealFieldFrozenError...

        """
        from . import ParametricRealField
        # fast path
        if hasattr(self, '_call_cache'):
            result = self._call_cache.get(x0)
            if result is not None:
                return result
        elif isinstance(self.parent().domain(), ParametricRealField):
            # Big Cells!
            with self.parent().domain().off_the_record():
                interval, function = self.which_pair(x0)
            if not is_pt_in_interval(interval, x0):
                raise AssertionError
            return function(x0)

        # Remember that intervals are sorted according to their left endpoints; singleton has priority.
        endpts = self.end_points()
        ith = self._ith_at_end_points
        i = bisect_left(endpts, x0)
        if i >= len(endpts):
            raise ValueError("Value not defined at point %s, outside of domain." % x0)
        if x0 == endpts[i]:
            if self._values_at_end_points[i] is not None:
                result = self._values_at_end_points[i]
                if hasattr(self, '_call_cache'):
                    self._call_cache[x0] = result
                return result
            else:
                raise ValueError("Value not defined at point %s, outside of domain." % x0)
        if i == 0:
            raise ValueError("Value not defined at point %s, outside of domain." % x0)
        if is_pt_in_interval(self._intervals[ith[i]],x0):
            result = self.functions()[ith[i]](x0)
            if hasattr(self, '_call_cache'):
                self._call_cache[x0] = result
            return result
        raise ValueError("Value not defined at point %s, outside of domain." % x0)

    def limits(self, x0):
        r"""
        returns [function value at `x_0`, function value at `x_0^+`, function value at `x_0^-`].

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp.fast_piecewise import *
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = 4
            sage: f5(x) = sin(2*x)
            sage: f6(x) = x-3
            sage: f7(x) = 7
            sage: f = FastPiecewise([[right_open_interval(0,1),f1],
            ....:                    [right_open_interval(1,2),f2],
            ....:                    [open_interval(2,3),f3],
            ....:                    [singleton_interval(3),f4],
            ....:                    [left_open_interval(3,6),f5],
            ....:                    [open_interval(6,7),f6],
            ....:                    [(9,10),f7]], periodic_extension=False)
            sage: f.limits(1/2)
            [1, 1, 1]
            sage: f.limits(1)
            [0, 0, 1]
            sage: f.limits(2)
            [None, e^2, -1]
            sage: f.limits(3)
            [4, sin(6), e^3]
            sage: f.limits(6)
            [sin(12), 3, sin(12)]
            sage: f.limits(7)
            [None, None, 4]
            sage: f.limits(8)
            [None, None, None]
            sage: f.limits(9)
            [7, 7, None]
        """
        endpts = self.end_points()
        ith = self._ith_at_end_points
        i = bisect_left(endpts, x0)
        if i >= len(endpts):
            return [None, None, None]
        if x0 == endpts[i]:
            return self.limits_at_end_points()[i]
        if i == 0:
            return [None, None, None]
        if is_pt_in_interval(self._intervals[ith[i]],x0):
            result = self.functions()[ith[i]](x0)
            return [result, result, result]
        return [None, None, None]

    def limit(self, x0, epsilon):
        r"""
        returns limit (from right if `\epsilon > 0`, from left if `\epsilon < 0`) value at `x_0`;
        if `\epsilon = 0`, returns value at `x_0`.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp.fast_piecewise import *
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = 4
            sage: f5(x) = sin(2*x)
            sage: f6(x) = x-3
            sage: f7(x) = 7
            sage: f = FastPiecewise([[right_open_interval(0,1),f1],
            ....:                    [right_open_interval(1,2),f2],
            ....:                    [open_interval(2,3),f3],
            ....:                    [singleton_interval(3),f4],
            ....:                    [left_open_interval(3,6),f5],
            ....:                    [open_interval(6,7),f6],
            ....:                    [(9,10),f7]], periodic_extension=False)
            sage: f.limit(1,0)
            0
            sage: f.limit(1,1)
            0
            sage: f.limit(2,-1)
            -1
            sage: f.limit(2,0)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 2, outside of domain.
            sage: f.limit(7,1)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 7+, outside of domain.
            sage: f.limit(8,-1)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 8-, outside of domain.
        """
        result =self.limits(x0)[epsilon]
        if result is None:
            from . import print_sign
            raise ValueError("Value not defined at point %s%s, outside of domain." % (x0, print_sign(epsilon)))
        return result

    def limit_within_relint(self, x0, interval):
        a, b = interval_to_endpoints(interval)
        if x0 == a < b:
            return self.limit(a, +1)
        elif a < b == x0:
            return self.limit(b, -1)
        elif a <= x0 <= b:
            return self(x0)
        else:
            raise ValueError("outside of face")

    def slope(self, x0, epsilon=1):
        endpts = self.end_points()
        limits = self.limits_at_end_points()
        i = bisect_left(endpts, x0)
        if i >= len(endpts):
            raise ValueError("Slope at %s is not defined, outside of domain." % x0)
        if x0 == endpts[i]:
            if epsilon > 0:
                if limits[i][epsilon] == limits[i][0]:
                    return (limits[i+1][-epsilon]-limits[i][0])/(endpts[i+1]-endpts[i])
                elif limits[i][epsilon] > limits[i][0]:
                    return +Infinity
                else:
                    return -Infinity
            elif epsilon < 0:
                if limits[i][epsilon] == limits[i][0]:
                    return (limits[i-1][-epsilon]-limits[i][0])/(endpts[i-1]-endpts[i])
                elif limits[i][epsilon] < limits[i][0]:
                    return +Infinity
                else:
                    return -Infinity
            else:
                raise ValueError("Slope at %s is not defined." % x0)
        if i == 0:
            raise ValueError("Slope at %s is not defined, outside of domain." % x0)
        return  (limits[i][-1]-limits[i-1][1])/(endpts[i]-endpts[i-1])

    def which_function_on_interval(self, interval):
        x = (interval[0] + interval[1]) / 2
        # FIXME: This should check that the given interval is contained in the defining interval!
        # This could be implemented by refactoring which_function using new function which_function_index.
        return self.which_function(x)

    def _add_(self,other):
        r"""
        Add self and another piecewise function.

        In contrast to ``PiecewisePolynomial.__add__``, this does not do zero extension of domains.
        Rather, the result is only defined on the intersection of the domains.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: f = FastPiecewise([[singleton_interval(1), FastLinearFunction(0,17)]])
            sage: g = FastPiecewise([[[0,2], FastLinearFunction(0,2)]])
            sage: (f+g).list()
            [[<Int{1}>, <FastLinearFunction 19>]]
            sage: h = FastPiecewise([[open_interval(1,3), FastLinearFunction(0,3)]])
            sage: (g+h).list()
            [[<Int(1, 2]>, <FastLinearFunction 5>]]
            sage: j = FastPiecewise([[open_interval(0,1), FastLinearFunction(0,1)], [[1, 3], FastLinearFunction(0, 5)]])
            sage: (g+j).list()
            [[<Int(0, 1)>, <FastLinearFunction 3>], [(1, 2), <FastLinearFunction 7>]]
        """
        intervals = intersection_of_coho_intervals([self.intervals(), other.intervals()])
        return self.parent()([ (interval, self.which_function_on_interval(interval) + other.which_function_on_interval(interval))
                               for interval in intervals ], merge=True)

    def _neg_(self):
        return self.parent()([[interval, -f] for interval,f in self.list()], merge=True)

    def _rmul_(self, scalar):
        r"""
        Multiply self by a scalar (element of base).
        """
        return self.parent()([[interval, scalar * f] for interval,f in self.list()])

    _lmul_ = _rmul_

    def _acted_upon_(self, actor, self_on_left):
        r"""
        Multiply self by a matrix from left or right.
        """
        if self_on_left:
            return FastPiecewise([[interval, f * actor] for interval,f in self.list()])
        else:
            return FastPiecewise([[interval, actor * f] for interval,f in self.list()])

    def __truediv__(self, other):
        r"""
        Divide by a scalar.
        """
        return self * (QQ(1) / other)

    __div__ = __truediv__

    def _sub_(self, other):
        return self + (-other)

    ## Following just fixes a bug in the plot method in piecewise_old.py
    ## (see doctests below).  Also adds plotting of single points
    ## and discontinuity markers.
    def plot(self, *args, **kwds):
        r"""
        Returns the plot of self.

        Keyword arguments are passed onto the plot command for each piece
        of the function. E.g., the ``plot_points`` keyword affects each
        segment of the plot.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = FastPiecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            sage: P = f.plot(rgbcolor=(0.7,0.1,0), plot_points=40)
            sage: P
            Graphics object...

        Remember: to view this, type ``show(P)`` or ``P.save("path/myplot.png")``
        and then open it in a graphics viewer such as GIMP.

        TESTS:

        We should not add each piece to the legend individually, since
        this creates duplicates (:trac:`12651`). This tests that only
        one of the graphics objects in the plot has a non-``None``
        ``legend_label``::

            sage: f1(x) = sin(x)
            sage: f2(x) = cos(x)
            sage: f = FastPiecewise([[(-1,0), f1],[(0,1), f2]])
            sage: p = f.plot(legend_label='$f(x)$')
            sage: lines = [
            ....:   line
            ....:   for line in p._objects
            ....:   if line.options()['legend_label'] is not None ]
            sage: len(lines)
            1

        The implementation of the plot method in Sage 5.11 piecewise.py
        is incompatible with the use of the xmin and xmax arguments.  Test that
        this has been fixed::

            sage: q = f.plot(xmin=0, xmax=3)
            sage: q = plot(f, xmin=0, xmax=3)
            sage: q = plot(f, 0, 3)
            sage: q = plot(f, 0, 3, color='red')

        The implementation should crop according to the given xmin, xmax::

            sage: q = plot(f, 1/2, 3)
            sage: q = plot(f, 1, 2)
            sage: q = plot(f, 2, 3)

        Also the following plot syntax should be accepted::

            sage: q = plot(f, [2, 3])

        """
        from sage.plot.all import plot, Graphics, point
        from . import delete_one_time_plot_kwds
        g = Graphics()
        if 'rgbcolor' in kwds:
            color=kwds['rgbcolor']
        elif 'color' in kwds:
            color=kwds['color']
        else:
            color = 'blue'
        if not 'plot_points' in kwds:
            plot_pts = 200
        else:
            plot_pts = kwds['plot_points']
        ### Code duplication with xmin/xmax code in plot.py.
        n = len(args)
        xmin = None
        xmax = None
        if n == 0:
            # if there are no extra args, try to get xmin,xmax from
            # keyword arguments
            xmin = kwds.pop('xmin', None)
            xmax = kwds.pop('xmax', None)
        elif n == 1:
            # if there is one extra arg, then it had better be a tuple
            xmin, xmax = args[0]
            args = []
            ## The case where the tuple is longer than 2 elements is for the
            ## case of symbolic expressions; it does not apply here.
            ## FIXME: We should probably signal an error.
        elif n == 2:
            # if there are two extra args, they should be xmin and xmax
            xmin = args[0]
            xmax = args[1]
            args = []
        ## The case with three extra args is for the case of symbolic
        ## expressions; it does not apply here.  FIXME: We should
        ## probably signal an error.
        point_kwds = dict()
        if 'alpha' in kwds:
            point_kwds['alpha'] = kwds['alpha']
        if 'legend_label' in kwds and self.is_discrete():
            point_kwds['legend_label'] = kwds['legend_label']
        # Whether to plot discontinuity markers
        discontinuity_markers = kwds.pop('discontinuity_markers', True)
        # record last right endpoint, then compare with next left endpoint to decide whether it needs to be plotted.
        last_end_point = []
        last_closed = True
        for (i, f) in self.list():
            a = i[0]
            b = i[1]
            left_closed = True
            right_closed = True
            if len(i) > 2: # coho interval
                left_closed = i.left_closed
                right_closed = i.right_closed
            # using the above data.
            if (xmin is not None) and (a < xmin):
                a = xmin
                left_closed = True
            if (xmax is not None) and (b > xmax):
                b = xmax
                right_closed = True
            if discontinuity_markers:
                # Handle open/half-open intervals here
                if a < b or (a == b and left_closed and right_closed):
                    if not (last_closed or last_end_point == [a, f(a)] and left_closed):
                        # plot last open right endpoint
                        g += point(last_end_point, color=color, pointsize=23, **point_kwds)
                        delete_one_time_plot_kwds(point_kwds)
                        g += point(last_end_point, rgbcolor='white', pointsize=10, **point_kwds)
                    if last_closed and last_end_point != [] and last_end_point != [a, f(a)] and not left_closed:
                        # plot last closed right endpoint
                        g += point(last_end_point, color=color, pointsize=23, **point_kwds)
                        delete_one_time_plot_kwds(point_kwds)
                    if not (left_closed or last_end_point == [a, f(a)] and last_closed):
                        # plot current open left endpoint
                        g += point([a, f(a)], color=color, pointsize=23, **point_kwds)
                        delete_one_time_plot_kwds(point_kwds)
                        g += point([a, f(a)], rgbcolor='white', pointsize=10, **point_kwds)
                    if left_closed and last_end_point != [] and last_end_point != [a, f(a)] and not last_closed:
                        # plot current closed left endpoint
                        g += point([a, f(a)], color=color, pointsize=23, **point_kwds)
                        delete_one_time_plot_kwds(point_kwds)
                    last_closed = right_closed
                    last_end_point = [b, f(b)]
            if a < b and (float(b) - float(a))/(plot_pts-1) != float(0):
                # We do not plot anything if (float(b) - float(a))/(plot_pts-1) == float(0) because
                # otherwise the plot method in src/plot/misc.py complains that
                # "start point and endpoint must be different"
                g += plot(f, *args, xmin=a, xmax=b, zorder=-1, **kwds)
                # If it's the first piece, pass all arguments. Otherwise,
                # filter out 'legend_label' so that we don't add each
                # piece to the legend separately (trac #12651).
                delete_one_time_plot_kwds(kwds)
                #delete_one_time_plot_kwds(point_kwds)
            elif a == b and left_closed and right_closed:
                g += point([a, f(a)], color=color, pointsize=23, **point_kwds)
                delete_one_time_plot_kwds(point_kwds)
        # plot open rightmost endpoint. minimal functions don't need this.
        if discontinuity_markers and not last_closed:
            g += point(last_end_point, color=color,pointsize=23, **point_kwds)
            delete_one_time_plot_kwds(point_kwds)
            g += point(last_end_point, rgbcolor='white', pointsize=10, **point_kwds)
        # For empty functions, if ticks were provided, use them (for uniformity).
        if not g:
            g._set_extra_kwds(kwds)
        return g

    def is_continuous_defined(self, xmin=0, xmax=1):
        r"""
        return ``True`` if self is defined on [xmin, xmax] and is continuous on [xmin, xmax].
        """
        bkpt = self._end_points
        if xmin < bkpt[0] or xmax > bkpt[-1]:
            return False
        if xmin == xmax:
            return (self(xmin) is not None)
        limits = self._limits_at_end_points
        i = 0
        while bkpt[i] < xmin:
            i += 1
        if bkpt[i] == xmin:
            if limits[i][0] is None or limits[i][1] is None or limits[i][0] != limits[i][1]:
                return False
            i += 1
        while bkpt[i] < xmax:
            if limits[i][-1] is None or limits[i][0] is None or limits[i][1] is None or \
                                        not (limits[i][-1] == limits[i][0] == limits[i][1]):
                return False
            i += 1
        if bkpt[i] == xmax:
            if limits[i][0] is None or limits[i][-1] is None or limits[i][0] != limits[i][-1]:
                return False
        return True

    def _repr_(self):
        from . import show_values_of_fastpiecewise
        rep = "<FastPiecewise with %s parts, " % len(self._functions)
        for interval, function in zip(self._intervals, self._functions):
            rep += "\n " + repr(interval) + "\t" + repr(function)
            if show_values_of_fastpiecewise:
                rep += "\t values: " + repr([function(interval[0]), function(interval[1])])
        rep += ">"
        return rep

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.
        """
        # FIXME: Add keyword arguments
        # FIXME: "sage_input(..., verify=True)" does not work yet
        # because of module trouble?
        return sib.name('FastPiecewise')(sib(self.list()))

    def sha1(self):
        r"""
        Return a SHA-1 hash of the function.

        The hash is intended to stay stable even when the code is updated.

        Merged and unmerged versions have the same hash.

        TESTS::

            sage: from cutgeneratingfunctionology.igp import *
            sage: h1 = piecewise_function_from_breakpoints_and_slopes([0, 1/4, 1/2, 1], [1, 1, -1], merge=False)
            sage: h1.sha1()
            'c562cf38581076609876b1c4fab604756690db7b'
            sage: h2 = piecewise_function_from_breakpoints_and_slopes([0, 1/4, 1/2, 1], [1, 1, -1], merge=True)
            sage: h2.sha1()
            'c562cf38581076609876b1c4fab604756690db7b'

        """
        from hashlib import sha1
        self_merged = self * ZZ(1)      # in case we were constructed with merge=False!
        data = list(zip(self_merged.end_points(), self_merged.limits_at_end_points()))
        from . import is_all_QQ
        from sage.misc.flatten import flatten
        try:
            from sage.misc import six
        except ImportError:    # Sage >= 9.2
            import six
        is_rational, _ = is_all_QQ(flatten(data))
        if not is_rational:
            logging.warning("For functions with non-rational data, cannot guarantee a stable SHA-1 hash.")
        stable_str = six.b(str(data))
        return sha1(stable_str).hexdigest()

    def _latex_(self, table=False, labels={}, name=r'\pi'):
        r"""
        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: f1(x) = 1
            sage: f2(x) = 1 - x
            sage: f = FastPiecewise([[(0,1),f1],[(1,2),f2]])
            sage: latex(f)
            \begin{cases}
            x \ {\mapsto}\ 1 &\text{on $(0, 1)$}\cr
            x \ {\mapsto}\ -x + 1 &\text{on $(1, 2)$}\cr
            \end{cases}

        ::

            sage: f(x) = sin(x*pi/2)
            sage: g(x) = 1-(x-1)^2
            sage: h(x) = -x
            sage: P = FastPiecewise([[(0,1), f], [(1,3),g], [(3,5), h]])
            sage: latex(P)
            \begin{cases}
            x \ {\mapsto}\ \sin\left(\frac{1}{2} \, \pi x\right) &\text{on $(0, 1)$}\cr
            x \ {\mapsto}\ -{\left(x - 1\right)}^{2} + 1 &\text{on $(1, 3)$}\cr
            x \ {\mapsto}\ -x &\text{on $(3, 5)$}\cr
            \end{cases}
        """
        if not table:
            # FIXME: This code from piecewise_old does not handle half-open intervals
            from sage.misc.latex import latex
            tex = ['\\begin{cases}\n']
            for (left, right), f in self.list():
                tex.append('%s &\\text{on $(%s, %s)$}\\cr\n' % (latex(f), left, right))
            tex.append(r'\end{cases}')
            return ''.join(tex)

        from sage.misc.latex import latex
        def labeled_latex(x):
            return labels.get(x, latex(x))
        latex.add_package_to_preamble_if_available("booktabs")
        s = []
        num_columns = 6
        s += [r'\begin{array}{*%sc}' % num_columns]
        s += [r'  \toprule']
        s += ['  ' + ' & '.join(['i',
                                 'x_i',
                                 name + '(x_i^-)',
                                 name + '(x_i)',
                                 name + '(x_i^+)',
                                 r'\text{slope}']) + r'\\']
        s += [r'  \midrule']
        end_points = self.end_points()
        for index, (bkpt, limits) in enumerate(zip(end_points, self.limits_at_end_points())):
            latex_limits = [ labeled_latex(x) for x in limits ]
            for eps in [-1, +1]:
                if limits[eps] == limits[0]:
                    latex_limits[eps] = ''
            slope = ''
            if index < len(end_points) - 1:
                slope = self.which_function((end_points[index] + end_points[index+1])/ 2)._slope
                slope = labeled_latex(slope)
            s += ['  ' + ' & '.join([labeled_latex(index),
                                     labeled_latex(bkpt),
                                     latex_limits[-1],
                                     latex_limits[0],
                                     latex_limits[1],
                                     r'\smash{\raisebox{-1.5ex}{$%s$}}' % slope]) + r'\\']
        s += [r'  \bottomrule']
        s += [r'\end{array}']
        return '\n'.join(s)

    ## @staticmethod
    ## def from_interval_lengths_and_slopes(cls, interval_lengths, slopes, field=None, merge=True):


FastPiecewise = PiecewiseLinearFunction_1d

from sage.structure.parent import Parent
from sage.categories.homset import Homset
from sage.structure.unique_representation import UniqueRepresentation

class PiecewiseFunctionsSpace(Homset):

    """
    Common base class of ``PiecewiseLinearFunctionsSpace`` and ``PiecewiseCrazyFunctionsSpace``.
    """

    def __init__(self, domain, codomain):
        from sage.categories.sets_cat import SetsWithPartialMaps
        from sage.categories.modules import Modules
        domain_cat = SetsWithPartialMaps()
        codomain_cat = codomain.category()
        if hasattr(codomain_cat, "WithBasis"):   # Adapted from Homset.
            # The above is a lame but fast check that category is a
            # subcategory of Modules(...).
            base_ring = codomain.base_ring()
        else:
            base_ring = codomain
        cat = domain_cat.Homsets() & Modules(base_ring)
        # Skip the Homset constructor because it uses the wrong category.
        # Following is copied from Homset.__init__
        self._domain = domain
        self._codomain = codomain
        self._Homset__category = domain_cat
        Parent.__init__(self, base=base_ring, category=cat)

    __contains__ = Parent.__contains__   # override Homset.__contains__, which is too strict

class PiecewiseLinearFunctionsSpace(UniqueRepresentation, PiecewiseFunctionsSpace):

    """
    Vector space of piecewise linear partial functions.

    EXAMPLES::

        sage: import logging; logging.disable(logging.WARN) # Suppress output in automatic tests.
        sage: from cutgeneratingfunctionology.igp import *
        sage: PWL = PiecewiseLinearFunctionsSpace(QQ, QQ)
        sage: PWL([((0, 1), FastLinearFunction(1, 1))])
        <FastPiecewise ...(0, 1)...x + 1...>

        sage: sqrt2, = nice_field_values([sqrt(2)])
        sage: PWL_NF = PiecewiseLinearFunctionsSpace(sqrt2.parent(), sqrt2.parent())
        sage: PWL_NF.base_ring()
        Real Number Field in `a` ...

    TESTS::

        sage: TestSuite(PWL).run(skip=['_test_zero'])

    Test that 'in' works::

        sage: gmic() in PiecewiseLinearFunctionsSpace(QQ, AA)
        True
        sage: hsqrt2 = gmic(sqrt(2)/2)
        sage: hsqrt2 in PiecewiseLinearFunctionsSpace(QQ, QQ)
        False
        sage: C = hsqrt2.parent().codomain(); C
        Real Number Field in `a` as the root of the defining polynomial y^2 - 2 near 1.414213562373095?
        sage: hsqrt2 in PiecewiseLinearFunctionsSpace(QQ, C)
        False
        sage: hsqrt2 in PiecewiseLinearFunctionsSpace(C, C)
        True

    """
    ## FIXME: We skip _test_zero because we need to improve equality testing.

    def _element_constructor_(self, x, *args, **kwds):
        """
        Construct a piecewise linear function.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: z = PiecewiseLinearFunctionsSpace(QQ, QQ)(0); z.list()
            [[(0, 1), <FastLinearFunction 0>]]
            sage: PiecewiseLinearFunctionsSpace(AA, AA)(z).parent()
            Vector space of piecewise linear partial functions from Algebraic Real Field to Algebraic Real Field over Algebraic Real Field
        """
        element_class = self.element_class
        if isinstance(x, FastPiecewise):
            return element_class(x.list(),
                                 periodic_extension=x._periodic_extension,
                                 merge=False, parent=self)
        if x in self.codomain():
            return element_class([((QQ(0), QQ(1)), FastLinearFunction(self.codomain().zero(),
                                                                      self.codomain()(x)))],
                                 periodic_extension=True, parent=self)
        return element_class(x, parent=self, *args, **kwds)

    def _repr_(self):
        return "Vector space of piecewise linear partial functions from {} to {} over {}".format(
            self.domain(), self.codomain(), self.base_ring())

    Element = FastPiecewise

    def _coerce_map_from_(self, S):
        """
        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: PiecewiseLinearFunctionsSpace(AA, AA).has_coerce_map_from(PiecewiseLinearFunctionsSpace(QQ, QQ))
            True

        TESTS::

            sage: alpha, = nice_field_values([sqrt(2)])
            sage: (alpha * gmic() + gmic()).parent()
            Vector space of piecewise linear partial functions from Rational Field to Real Number Field in `a` as the root of the defining polynomial y^2 - 2 near 1.414213562373095? over Real Number Field in `a` as the root of the defining polynomial y^2 - 2 near 1.414213562373095?

        """
        return isinstance(S, PiecewiseLinearFunctionsSpace) and self.codomain().has_coerce_map_from(S.codomain())

    def _an_element_(self):
        return self([((QQ(0), QQ(1)), FastLinearFunction(QQ(0), QQ(0))),
                     ((QQ(1), QQ(2)), FastLinearFunction(QQ(1), QQ(-1)))])

    # FIXME: Implement _get_action_ to get action of matrix on FastPiecewise right
