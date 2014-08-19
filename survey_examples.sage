## Various examples of functions that appear in the survey.

def not_minimal_1(): # was not_minimal.sage
    return piecewise_function_from_breakpoints_and_values([0, 1/5, 2/5, 4/5, 1], [0, 1/5, 3/5, 1, 0])

def not_minimal_2(): # was not_minimal_2.sage
    return piecewise_function_from_breakpoints_and_values([0, 1/5, 2/5, 3/5, 4/5, 1], [0, 1/5, 1/2, 4/5, 1, 0])

def not_extreme_1(): # was symmetric_rational_function1.sage
    slopes = [10/3,0,10/3,0,10/3,-10/3,0,-10/3,0,-10/3]
    interval_lengths = [1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10]
    return piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes)

