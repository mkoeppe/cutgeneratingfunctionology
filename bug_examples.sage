# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

def not_minimal_3(): # this was a bug
    """
    sage: logging.disable(logging.INFO); 
    sage: h = not_minimal_3()
    sage: minimality_test(h, False)
    False
    """
    return piecewise_function_from_breakpoints_and_values([0,1/5,4/5,1],[0,1/2,1,0])

def not_minimal_wrong_range():
    """
    sage: logging.disable(logging.INFO); 
    sage: h = not_minimal_wrong_range()
    sage: minimality_test(h, False)
    False
    """
    return piecewise_function_from_breakpoints_and_values([0,1/2,1], [0,2,0])

def fake_f():
    """
    sage: logging.disable(logging.INFO); 
    sage: h = fake_f()
    sage: minimality_test(h, f=4/5)
    False
    sage: minimality_test(h, f=1/5)
    False
    """
    return piecewise_function_from_breakpoints_and_values([0,1/5,3/5,4/5,1],[0,1,0,1,0])

def limits_out_of_range():                                  # plotting bug
    """
    sage: logging.disable(logging.INFO); 
    sage: h = limits_out_of_range()
    sage: minimality_test(h, False)
    False
    """
    return FastPiecewise([[singleton_interval(0), FastLinearFunction(0,0)], [open_interval(0, 1/2), FastLinearFunction(6, -1)], [closed_interval(1/2,1), FastLinearFunction(-2, 2)]], merge=False)

def chen_tricky_uncovered_intervals():
    """
    sage: logging.disable(logging.INFO); 
    sage: h = chen_tricky_uncovered_intervals()
    sage: extremality_test(h, False)
    False
    """
    return chen_3_slope_not_extreme(f=1/sqrt(3), lam=10)    

def minimal_no_covered_interval():
    """
    sage: logging.disable(logging.WARN) 
    sage: h = minimal_no_covered_interval()
    sage: extremality_test(h, False)
    False
    """
    return FastPiecewise([[singleton_interval(0), FastLinearFunction(0, 0)], \
                          [open_interval(0, 1/2), FastLinearFunction(0, 1/2)], \
                          [singleton_interval(1/2), FastLinearFunction(0, 1)], \
                          [open_interval(1/2, 1), FastLinearFunction(0, 1/2)], \
                          [singleton_interval(1), FastLinearFunction(0, 0)]], merge=True)

def minimal_has_uncovered_interval():
    """
    sage: logging.disable(logging.WARN)
    sage: h = minimal_has_uncovered_interval()
    sage: extremality_test(h, False)
    False
    """
    return FastPiecewise([[singleton_interval(0), FastLinearFunction(0, 0)], \
                          [open_interval(0, 1/8), FastLinearFunction(0, 3/4)],\
                          [singleton_interval(1/8), FastLinearFunction(0, 1/2)], \
                          [open_interval(1/8, 1/4), FastLinearFunction(0, 1/4)], \
                          [singleton_interval(1/4), FastLinearFunction(0, 1)], \
                          [open_interval(1/4, 1), FastLinearFunction(0, 1/2)], \
                          [singleton_interval(1), FastLinearFunction(0,0)]], merge=True)

def lift_of_minimal_has_uncovered_interval():
    """
    This function was obtained by:
    sage_input(lift(minimal_has_uncovered_interval()))

    The function has 3 slopes and discontinuities.

    sage: logging.disable(logging.WARN)
    sage: h = lift_of_minimal_has_uncovered_interval()
    sage: len(generate_covered_intervals(h) + generate_uncovered_intervals(h)) >= 2
    True
    """
    return FastPiecewise([[singleton_interval(QQ(0)), FastLinearFunction(QQ(0), QQ(0))], [open_interval(0, 1/8), FastLinearFunction(QQ(0), 3/4)], [singleton_interval(1/8), FastLinearFunction(QQ(0), 1/2)], [open_interval(1/8, 1/4), FastLinearFunction(QQ(0), 1/4)], [singleton_interval(1/4), FastLinearFunction(QQ(0), QQ(1))], [left_open_interval(1/4, 1/2), FastLinearFunction(QQ(0), 1/2)], [left_open_interval(1/2, 9/16), FastLinearFunction(QQ(4), -3/2)], [left_open_interval(9/16, 11/16), FastLinearFunction(-QQ(4), QQ(3))], [left_open_interval(11/16, 3/4), FastLinearFunction(QQ(4), -5/2)], [open_interval(3/4, 1), FastLinearFunction(QQ(0), 1/2)], [singleton_interval(1), FastLinearFunction(QQ(0), QQ(0))]])

def lift_of_minimal_no_covered_interval():
    """
    sage: logging.disable(logging.WARN)
    sage: h = lift_of_minimal_no_covered_interval()
    sage: extremality_test(h)
    False
    """
    return FastPiecewise([[singleton_interval(QQ(0)), FastLinearFunction(-QQ(2), QQ(0))], [left_open_interval(0, 1/16), FastLinearFunction(-QQ(2), 1/2)], [left_open_interval(1/16, 3/32), FastLinearFunction(QQ(10), -1/4)], [left_open_interval(3/32, 1/8), FastLinearFunction(-QQ(14), QQ(2))], [left_open_interval(1/8, 3/8), FastLinearFunction(QQ(2), QQ(0))], [left_open_interval(3/8, 13/32), FastLinearFunction(-QQ(14), QQ(6))], [left_open_interval(13/32, 7/16), FastLinearFunction(QQ(10), -15/4)], [open_interval(7/16, 1/2), FastLinearFunction(-QQ(2), 3/2)], [singleton_interval(1/2), FastLinearFunction(-QQ(2), QQ(2))], [left_open_interval(1/2, 5/8), FastLinearFunction(QQ(2), -1/2)], [left_open_interval(5/8, 7/8), FastLinearFunction(-QQ(2), QQ(2))], [open_interval(7/8, 1), FastLinearFunction(QQ(2), -3/2)], [singleton_interval(1), FastLinearFunction(QQ(2), -QQ(2))]])

def example7slopecoarse2():
    """
    sage: logging.disable(logging.INFO); 
    sage: h = example7slopecoarse2()
    sage: extremality_test(h, False)
    False
    """
    return piecewise_function_from_robert_txt_file("example7slopeCoarse2.txt")
