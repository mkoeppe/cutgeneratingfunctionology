def not_minimal_3(): # this was a bug
    return piecewise_function_from_breakpoints_and_values([0,1/5,4/5,1],[0,1/2,1,0])

def not_minimal_wrong_range():
    return piecewise_function_from_breakpoints_and_values([0,1/2,1], [0,2,0])

def fake_f():
    return piecewise_function_from_breakpoints_and_values([0,1/5,3/5,4/5,1],[0,1,0,1,0])

def limits_out_of_range():                                  # plotting bug
    return FastPiecewise([[singleton_interval(0), FastLinearFunction(0,0)], [open_interval(0, 1/2), FastLinearFunction(6, -1)], [closed_interval(1/2,1), FastLinearFunction(-2, 2)]], merge=False)

def chen_tricky_uncovered_intervals():
    return chen_3_slope_not_extreme(f=1/sqrt(3), lam=10)    

def minimal_uncovered_interval_example()
    return FastPiecewise([[singleton_interval(0), FastLinearFunction(0,0)], [open_interval(0, 1/2), FastLinearFunction(0, 1/2)], [singleton_interval(1/2), FastLinearFunction(0,1)], [open_interval(1/2, 1), FastLinearFunction(0, 1/2)],[singleton_interval(1), FastLinearFunction(0,0)]], merge=True)
