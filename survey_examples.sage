## Various examples of functions that appear in the survey.

def not_minimal_1(): # was not_minimal.sage
    return piecewise_function_from_breakpoints_and_values([0, 1/5, 2/5, 4/5, 1], [0, 1/5, 3/5, 1, 0])

def not_minimal_2(): # was not_minimal_2.sage
    return piecewise_function_from_breakpoints_and_values([0, 1/5, 2/5, 3/5, 4/5, 1], [0, 1/5, 1/2, 4/5, 1, 0])

def not_extreme_1(): # was symmetric_rational_function1.sage
    slopes = [10/3,0,10/3,0,10/3,-10/3,0,-10/3,0,-10/3]
    interval_lengths = [1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10]
    return piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes)

def drlm_not_extreme_1():
    """Example from S. S. Dey, J.-P. P. Richard, Y. Li, and L. A. Miller,
    On the extreme inequalities of infinite group problems,
    Mathematical Programming 121 (2009), no. 1, 145–170,
    doi:10.1007/s10107-008-0229-6.
    Figure 1.

    It is the interpolation of the extreme function C13 for the finite
    group problem of order 7 from Araoz, Evans, Gomory and Johnson.
    """
    return piecewise_function_from_robert_txt_file("dey-richard-not-extreme.txt")

def bhk_irrational_extreme_limit_to_rational_nonextreme(n=Infinity):
    """
    A sequence of `bhk_irrational` functions, each extreme, indexed by n = 1, 2, ...
    whose limit (n = Infinity) is a `bhk_irrational` function with rational parameters, 
    and hence not extreme. 
    """
    del1 = 1/40 
    if n != Infinity:
        del1 -= sqrt(2) / (70*n)
    return the_irrational_function_t1_t2(del1=del1, del2=del2)

