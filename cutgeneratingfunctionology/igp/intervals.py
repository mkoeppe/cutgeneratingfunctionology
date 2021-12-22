## 
## A lightweight representation of closed bounded intervals, possibly empty or degenerate.
##

def interval_sum(int1, int2):
    r"""
    Return the sum of two intervals.
    """
    if len(int1) == 0 or len(int2) == 0:
        return []
    if len(int1) == 1 and len(int2) == 1:
        return [int1[0]+int2[0]]
    elif len(int1) == 2 and len(int2) == 1:
        return [int1[0]+int2[0],int1[1]+int2[0]]
    elif len(int1) == 1 and len(int2) == 2:
        return [int1[0]+int2[0],int1[0]+int2[1]]
    else:    
        return [int1[0]+int2[0],int1[1]+int2[1]]

def interval_intersection(int1, int2):
    r"""
    Return the intersection of two intervals.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: interval_intersection([1], [2])
        []
        sage: interval_intersection([1,3], [2,4])
        [2, 3]
        sage: interval_intersection([1,3], [2])
        [2]
        sage: interval_intersection([2], [2,3])
        [2]
        sage: interval_intersection([1,3], [3, 4])
        [3]
        sage: interval_intersection([1,3], [4, 5])
        []
        sage: interval_intersection([], [4, 5])
        []
    """
    if len(int1) == 0 or len(int2) == 0:
        return []
    if len(int1) == 1 and len(int2) == 1:
        if int1[0] == int2[0]:
            return [int1[0]]
        else:
            return []
    elif len(int1) == 2 and len(int2) == 1:
        if int1[0] <= int2[0] <= int1[1]:
            return [int2[0]]
        else:
            return []
    elif len(int1) == 1 and len(int2) == 2:
        if int2[0] <= int1[0] <= int2[1]:
            return [int1[0]]
        else:
            return []    
    else:        
        max0 = max(int1[0],int2[0])
        min1 = min(int1[1],int2[1])
        if max0 > min1:
            return []
        elif max0 == min1:
            return [max0]
        else:
            return [max0,min1]

def interval_empty(interval):
    r"""
    Determine whether an interval is empty.            
    """
    if len(interval) == 0:
        return True
    else:
        return False

def interval_equal(int1, int2):
    r"""
    Determine whether two intervals are equal.
    (This ignores whether the intervals are represented as tuples or lists.)
    """
    return tuple(int1) == tuple(int2)

def element_of_int(x,int):
    r"""
    Determine whether value `x` is inside the interval int.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: element_of_int(1, [])
        False
        sage: element_of_int(1, [1])
        True
        sage: element_of_int(1, [2])
        False
        sage: element_of_int(1, [0,2])
        True
        sage: element_of_int(1, [1,2])
        True
        sage: element_of_int(2, [3,4])
        False
    """
    if len(int) == 0:
        return False
    elif len(int) == 1:
        if x == int[0]:
            return True
        else:
            return False
    elif int[0] <= x <= int[1]:
        return True
    else:
        return False

def interval_to_endpoints(int):
    r"""
    Convert a possibly degenerate interval to a pair (a,b)
    of its endpoints, suitable for specifying pieces of a ``FastPiecewise``.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: interval_to_endpoints([1])
        (1, 1)
        sage: interval_to_endpoints([1,3])
        (1, 3)
    """
    if len(int) == 0:
        raise ValueError("An empty interval does not have a pair representation")
    elif len(int) == 1:
        return (int[0], int[0])
    elif len(int) == 2:
        return (int[0], int[1])
    else:
        raise ValueError("Not an interval: %s" % (int,))

# Assume the lists of intervals are sorted.                
def find_interior_intersection(list1, list2):
    r"""
    Tests whether list1 and list2 contain a pair of intervals
    whose interiors intersect.

    Assumes both lists are sorted.
    
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: find_interior_intersection([[1, 2], [3, 4]], [[2, 3], [4, 5]])
        False
        sage: find_interior_intersection([[1, 2], [3, 5]], [[2, 4]])
        True
    """
    i=0
    j=0
    while i < len(list1) and j < len(list2):
        if len(interval_intersection(list1[i], list2[j])) == 2:
            return True
        else:
            if list1[i][0] < list2[j][0]:
                i = i + 1
            else:
                j = j + 1
    return False

def interval_minus_union_of_intervals(interval, remove_list):
    r"""Compute a list of intervals that represent the
    set difference of interval and the union of the 
    intervals in remove_list.

    Assumes remove_list is sorted (and pairwise essentially
    disjoint), and returns a sorted list.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: interval_minus_union_of_intervals([0, 10], [[-1, 0], [2, 3], [9,11]]) 
        [[0, 2], [3, 9]]
        sage: interval_minus_union_of_intervals([0, 10], [[-1, 0], [2, 3]]) 
        [[0, 2], [3, 10]]
        sage: interval_minus_union_of_intervals([0, 10], [[-1, 0], [2, 3], [9,11], [13, 17]])
        [[0, 2], [3, 9]]
    """
    scan = scan_union_of_coho_intervals_minus_union_of_coho_intervals([[interval]], [remove_list])
    return list(proper_interval_list_from_scan(scan))

##
## A representation of bounded intervals, closed, open, or half-open
##
        
import collections
_closed_or_open_or_halfopen_interval = collections.namedtuple('Interval', ['a', 'b', 'left_closed', 'right_closed'])

class closed_or_open_or_halfopen_interval (_closed_or_open_or_halfopen_interval):
    def __repr__(self):
        if self.a == self.b and self.left_closed and self.right_closed:
            r = "{" + repr(self.a) + "}"
        else:
            r = ("[" if self.left_closed else "(") \
                + repr(self.a) + ", " + repr(self.b) \
                + ("]" if self.right_closed else ")")
        return "<Int" + r + ">"

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.
        """
        if self.a == self.b and self.left_closed and self.right_closed:
            return sib.name('singleton_interval')(sib(self.a))
        else:
            if self.left_closed and self.right_closed:
                name = 'closed_interval'
            elif self.left_closed and not self.right_closed:
                name = 'right_open_interval'
            elif not self.left_closed and self.right_closed:
                name = 'left_open_interval'
            else:
                name = 'open_interval'
            return sib.name(name)(sib(self.a), sib(self.b))

def closed_interval(a, b):
    return closed_or_open_or_halfopen_interval(a, b, True, True)

def open_interval(a, b):
    return closed_or_open_or_halfopen_interval(a, b, False, False)

def singleton_interval(a):
    return closed_or_open_or_halfopen_interval(a, a, True, True)

def left_open_interval(a, b):
    return closed_or_open_or_halfopen_interval(a, b, False, True)

def right_open_interval(a, b):
    return closed_or_open_or_halfopen_interval(a, b, True, False)

def coho_interval_from_interval(int):
    if len(int) == 0:
        raise ValueError("An empty interval does not have a coho_interval representation")
    elif len(int) == 1:
        return singleton_interval(int[0])
    elif len(int) == 2:
        return closed_interval(int[0], int[1])
    else:
        raise ValueError("Not an interval: %s" % (int,))

def realset_from_interval(int):
    from cutgeneratingfunctionology.spam.real_set import RealSet, InternalRealInterval
    if len(int) == 0:
        return RealSet()
    elif len(int) == 1:
        return RealSet.point(int[0])
    elif len(int) == 2:
        return RealSet.closed(int[0], int[1])
    else:
        return RealSet(InternalRealInterval(int.a, int.left_closed, int.b, int.right_closed))

def interval_length(interval):
    r"""
    Determine the length of the given interval.

    interval can be old-fashioned or coho.
    """
    if len(interval) <= 1:
        return 0
    elif interval[1] >= interval[0]:
        return interval[1] - interval[0]
    else:
        return 0

def is_pt_in_interval(i, x0):
    r"""
    retrun whether the point x0 is contained in the (ordinary or coho) interval i.
    """
    if len(i) == 0:
        return False
    if len(i) == 1:
        return i[0] == x0
    if len(i) == 2:
        return bool(i[0] <= x0 <= i[1])
    else:  
        if i.left_closed and i.right_closed:
            return bool(i.a <= x0 <= i.b)
        if i.left_closed and not i.right_closed:
            return bool(i.a <= x0 < i.b)
        if not i.left_closed and i.right_closed:
            return bool(i.a < x0 <= i.b)
        if not i.left_closed and not i.right_closed:
            return bool(i.a < x0 < i.b)

def coho_interval_left_endpoint_with_epsilon(i, closure=False):
    r"""Return (x, epsilon)
    where x is the left endpoint
    and epsilon is 0 if the interval is left closed and 1 otherwise.
    """
    if len(i) == 0:
        raise ValueError("An empty interval does not have a left endpoint.")
    elif len(i) <= 2:
        # old-fashioned closed interval or singleton
        return i[0], 0 # Scanning from the left, turn on at left endpoint.
    else:
        # coho interval
        return i.a, 0 if (i.left_closed or closure) else 1

def coho_interval_right_endpoint_with_epsilon(i, closure=False):
    r"""Return (x, epsilon)
    where x is the right endpoint
    and epsilon is 1 if the interval is right closed and 0 otherwise.
    """
    if len(i) == 0:
        raise ValueError("An empty interval does not have a right endpoint.")
    elif len(i) == 1:
        # old-fashioned singleton
        return i[0], 1 # Scanning from the left, turn off at that point plus epsilon
    elif len(i) == 2:
        # old-fashioned proper closed interval
        return i[1], 1 # Scanning from the left, turn off at right endpoint plus epsilon
    else:
        # coho interval
        return i.b, 1 if (i.right_closed or closure) else 0

def coho_interval_contained_in_coho_interval(I, J):
    I = (coho_interval_left_endpoint_with_epsilon(I), coho_interval_right_endpoint_with_epsilon(I))
    J = (coho_interval_left_endpoint_with_epsilon(J), coho_interval_right_endpoint_with_epsilon(J))
    return J[0] <= I[0] and I[1] <= J[1]

def scan_coho_interval_left_endpoints(interval_list, tag=None, closure=False):
    r"""Generate events of the form ``(x, epsilon), delta=-1, tag``.

    This assumes that interval_list is sorted from left to right,
    and that the intervals are pairwise disjoint.
    """
    for i in interval_list:
        yield coho_interval_left_endpoint_with_epsilon(i, closure), -1, tag

def scan_coho_interval_right_endpoints(interval_list, tag=None, closure=False):
    r"""Generate events of the form ``(x, epsilon), delta=+1, tag``.

    This assumes that interval_list is sorted from left to right,
    and that the intervals are pairwise disjoint.
    """
    for i in interval_list:
        yield coho_interval_right_endpoint_with_epsilon(i, closure), +1, tag

def scan_coho_interval_list(interval_list, tag=None, closure=False):
    r"""Generate events of the form ``(x, epsilon), delta, tag``.

    This assumes that interval_list is sorted, and 
    that the intervals are pairwise disjoint. (disjoint needed?)

    delta is -1 for the beginning of an interval ('on').
    delta is +1 for the end of an interval ('off'). 

    This is so that the events sort lexicographically in a way that if
    we have intervals whose closures intersect in one point, such as
    [a, b) and [b, c], we see first the 'on' event and then the 'off'
    event.  In this way consumers of the scan can easily implement merging 
    of such intervals. 

    if closure is ``True``, considers intervals as closed.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: list(scan_coho_interval_list([closed_or_open_or_halfopen_interval(1, 2, True, False), closed_or_open_or_halfopen_interval(2, 3, True, True)]))
        [((1, 0), -1, None), ((2, 0), -1, None), ((2, 0), 1, None), ((3, 1), 1, None)]
    """
    return merge(scan_coho_interval_left_endpoints(interval_list, tag, closure),
                 scan_coho_interval_right_endpoints(interval_list, tag, closure))

## def scan_set_difference(a, b):
##     r"""`a` and `b` should be event generators."""

from heapq import *

def scan_union_of_coho_intervals_minus_union_of_coho_intervals(interval_lists, remove_lists, remove_closure=False):
    # Following uses the lexicographic comparison of the tuples.
    scan = merge(merge(*[scan_coho_interval_list(interval_list, True) for interval_list in interval_lists]),
                 merge(*[scan_coho_interval_list(remove_list, False, remove_closure) for remove_list in remove_lists]))
    interval_indicator = 0
    remove_indicator = 0
    on = False
    for ((x, epsilon), delta, tag) in scan:
        was_on = on
        if tag:                                       # interval event
            interval_indicator -= delta
            assert(interval_indicator) >= 0, "interval_indicator should stay non-negative; required sort order must be violated"
        else:                                           # remove event
            remove_indicator -= delta
            assert(remove_indicator) >= 0, "remove_indicator should stay non-negative; required sort order must be violated"
        now_on = interval_indicator > 0 and remove_indicator == 0
        if not was_on and now_on: # switched on
            yield (x, epsilon), -1, None
        elif was_on and not now_on: # switched off
            yield (x, epsilon), +1, None
        on = now_on
    # No unbounded intervals:
    assert interval_indicator == 0, "interval_indicator should be zero at the end of the scan"
    assert remove_indicator == 0, "remove_indicator should be zero at the end of the scan"

def intersection_of_coho_intervals(interval_lists):
    r"""Compute the intersection of the union of intervals. 
    Actually returns a generator.
    
    Each interval_list must be sorted, but intervals may overlap.  In
    this case, the output is broken into non-overlapping intervals at
    the points where the overlap multiplicity changes.
    
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: list(intersection_of_coho_intervals([[[1,2]], [[2,3]]]))
        [<Int{2}>]
        sage: list(intersection_of_coho_intervals([[[1,2], [2,3]], [[0,4]]]))
        [<Int[1, 2)>, <Int{2}>, <Int(2, 3]>]
        sage: list(intersection_of_coho_intervals([[[1,3], [2,4]], [[0,5]]]))
        [<Int[1, 2)>, <Int[2, 3]>, <Int(3, 4]>]
        sage: list(intersection_of_coho_intervals([[[1,2], left_open_interval(2,3)], [[0,4]]]))
        [<Int[1, 2]>, <Int(2, 3]>]
        sage: list(intersection_of_coho_intervals([[[1,3]], [[2,4]]]))
        [<Int[2, 3]>]
    """
    scan = merge(*[scan_coho_interval_list(interval_list, tag=index) for index, interval_list in enumerate(interval_lists)])
    interval_indicators = [ 0 for interval_list in interval_lists ]
    (on_x, on_epsilon) = (None, None)
    for ((x, epsilon), delta, index) in scan:
        was_on = all(on > 0 for on in interval_indicators)
        interval_indicators[index] -= delta
        assert interval_indicators[index] >= 0
        now_on = all(on > 0 for on in interval_indicators)
        if was_on: 
            assert on_x is not None
            assert on_epsilon >= 0
            assert epsilon >= 0
            if (on_x, on_epsilon) < (x, epsilon):
                yield closed_or_open_or_halfopen_interval(on_x, x,
                                                          on_epsilon == 0, epsilon > 0)
        if now_on:
            (on_x, on_epsilon) = (x, epsilon)
        else:
            (on_x, on_epsilon) = (None, None)
    assert all(on == 0 for on in interval_indicators) # no unbounded intervals

def coho_intervals_intersecting(a, b):
    r"""
    Determine if the two intervals intersect in at least 1 point.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: coho_intervals_intersecting(singleton_interval(1), singleton_interval(1))
        True
        sage: coho_intervals_intersecting(singleton_interval(1), singleton_interval(2))
        False
        sage: coho_intervals_intersecting(singleton_interval(1), open_interval(1,2))
        False
        sage: coho_intervals_intersecting(singleton_interval(1), right_open_interval(1,2))
        True
    """
    intervals = list(intersection_of_coho_intervals([[a], [b]]))
    assert len(intervals) <= 1
    return len(intervals) == 1

def coho_intervals_intersecting_full_dimensionally(a, b):
    r"""
    Determine if the two intervals intersect in a proper interval.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: coho_intervals_intersecting_full_dimensionally(singleton_interval(1), singleton_interval(1))
        False
        sage: coho_intervals_intersecting_full_dimensionally(singleton_interval(1), singleton_interval(2))
        False
        sage: coho_intervals_intersecting_full_dimensionally(singleton_interval(1), open_interval(1,2))
        False
        sage: coho_intervals_intersecting_full_dimensionally(singleton_interval(1), right_open_interval(1,2))
        False
        sage: coho_intervals_intersecting_full_dimensionally(open_interval(0,2), right_open_interval(1,3))
        True
    """
    intervals = list(intersection_of_coho_intervals([[a], [b]]))
    assert len(intervals) <= 1
    return len(intervals) == 1 and interval_length(intervals[0]) > 0

def coho_interval_list_from_scan(scan, old_fashioned_closed_intervals=False):
    r"""Actually returns a generator."""
    indicator = 0
    (on_x, on_epsilon) = (None, None)
    for ((x, epsilon), delta, tag) in scan:
        was_on = indicator > 0
        indicator -= delta
        assert indicator >= 0
        now_on = indicator > 0
        if not was_on and now_on:                        # switched on
            (on_x, on_epsilon) = (x, epsilon)
        elif was_on and not now_on:                     # switched off
            assert on_x is not None
            assert on_epsilon >= 0
            assert epsilon >= 0
            if (on_x, on_epsilon) < (x, epsilon):
                left_closed = on_epsilon == 0
                right_closed = epsilon > 0
                if old_fashioned_closed_intervals and left_closed and right_closed and on_x < x:
                    yield (on_x, x)
                else:
                    yield closed_or_open_or_halfopen_interval(on_x, x, left_closed, right_closed)
            (on_x, on_epsilon) = (None, None)
    assert indicator == 0

def union_of_coho_intervals_minus_union_of_coho_intervals(interval_lists, remove_lists, old_fashioned_closed_intervals=False, remove_closure=False):
    r"""Compute a list of closed/open/half-open intervals that represent
    the set difference of interval and the union of the intervals in
    remove_lists.

    Assume each of the lists in interval_lists and remove_lists are sorted (and
    each pairwise disjoint).  Returns a sorted list.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: union_of_coho_intervals_minus_union_of_coho_intervals([[[0,10]]], [[[2,2], [3,4]]])
        [<Int[0, 2)>, <Int(2, 3)>, <Int(4, 10]>]
        sage: union_of_coho_intervals_minus_union_of_coho_intervals([[[0, 10]]], [[[1, 7]], [[2, 5]]])
        [<Int[0, 1)>, <Int(7, 10]>]
        sage: union_of_coho_intervals_minus_union_of_coho_intervals([[[0,10], closed_or_open_or_halfopen_interval(10, 20, False, True)]], [])
        [<Int[0, 20]>]
        sage: union_of_coho_intervals_minus_union_of_coho_intervals([[[0,10], closed_or_open_or_halfopen_interval(10, 20, False, True)]], [], old_fashioned_closed_intervals=True)
        [(0, 20)]
    """
    gen = coho_interval_list_from_scan(scan_union_of_coho_intervals_minus_union_of_coho_intervals(interval_lists, remove_lists, remove_closure), old_fashioned_closed_intervals)
    return list(gen)

def proper_interval_list_from_scan(scan):
    r"""Return a generator of the proper intervals [a, b], a<b, in the scan.

    This ignores whether intervals are open/closed/half-open.
    """
    indicator = 0
    (on_x, on_epsilon) = (None, None)
    for ((x, epsilon), delta, tag) in scan:
        was_on = indicator > 0
        indicator -= delta
        assert indicator >= 0
        now_on = indicator > 0
        if not was_on and now_on:                        # switched on
            (on_x, on_epsilon) = (x, epsilon)
        elif was_on and not now_on:                     # switched off
            assert on_x is not None
            assert on_epsilon >= 0
            assert epsilon >= 0
            if on_x < x:
                yield [on_x, x]
            (on_x, on_epsilon) = (None, None)
    assert indicator == 0

def convert_interval(interval, ring):
    a = ring(interval[0])
    b = ring(interval[1])
    if len(interval) == 2:  # old fashioned closed interval
        return (a, b)
    return closed_or_open_or_halfopen_interval(a, b, interval.left_closed, interval.right_closed)
