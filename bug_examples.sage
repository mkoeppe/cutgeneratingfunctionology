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
    sage: simple_finite_dimensional_extremality_test(h, oversampling=4)
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
    bkpt = [0, 1/24, 1/12, 1/8, 1/6, 5/24, 7/24, 1/3, 3/8, 5/12, 11/24, 1/2, \
            13/24, 7/12, 5/8, 2/3, 5/6, 7/8, 11/12, 23/24, 1]
    values = [0, 3/4, 1/4, 3/4, 1/2, 3/4, 1/4, 1/2, 1/4, 3/4, 1/4, 1, \
              1/4, 1/2, 1/4, 1/2, 1/2, 3/4, 1/2, 3/4, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def example7slopecoarse2_lifted():
    """
    obtained via: 
    h = example7slopecoarse2(); lift_until_extreme(h, finite_dimensional_test_first=True); 
    hl = last_lifted(h).

    sage: logging.disable(logging.INFO)
    sage: h = example7slopecoarse2_lifted()
    sage: extremality_test(h, False)
    False
    """
    return FastPiecewise([[(QQ(0), 1/24), FastLinearFunction(QQ(18), QQ(0))], [left_open_interval(1/24, 1/16), FastLinearFunction(-QQ(6), QQ(1))], [left_open_interval(1/16, 1/12), FastLinearFunction(-QQ(18), 7/4)], [left_open_interval(1/12, 5/48), FastLinearFunction(QQ(18), -5/4)], [left_open_interval(5/48, 1/8), FastLinearFunction(QQ(6), QQ(0))], [left_open_interval(1/8, 1/6), FastLinearFunction(-QQ(6), 3/2)], [left_open_interval(1/6, 5/24), FastLinearFunction(QQ(6), -1/2)], [left_open_interval(5/24, 7/24), FastLinearFunction(-QQ(6), QQ(2))], [left_open_interval(7/24, 1/3), FastLinearFunction(QQ(6), -3/2)], [left_open_interval(1/3, 3/8), FastLinearFunction(-QQ(6), 5/2)], [left_open_interval(3/8, 19/48), FastLinearFunction(QQ(6), -QQ(2))], [left_open_interval(19/48, 5/12), FastLinearFunction(QQ(18), -27/4)], [left_open_interval(5/12, 7/16), FastLinearFunction(-QQ(18), 33/4)], [left_open_interval(7/16, 11/24), FastLinearFunction(-QQ(6), QQ(3))], [left_open_interval(11/24, 1/2), FastLinearFunction(QQ(18), -QQ(8))], [left_open_interval(1/2, 13/24), FastLinearFunction(-QQ(18), QQ(10))], [left_open_interval(13/24, 7/12), FastLinearFunction(QQ(6), -QQ(3))], [left_open_interval(7/12, 5/8), FastLinearFunction(-QQ(6), QQ(4))], [left_open_interval(5/8, 2/3), FastLinearFunction(QQ(6), -7/2)], [left_open_interval(2/3, 5/6), FastLinearFunction(QQ(0), 1/2)], [left_open_interval(5/6, 7/8), FastLinearFunction(QQ(6), -9/2)], [left_open_interval(7/8, 11/12), FastLinearFunction(-QQ(6), QQ(6))], [left_open_interval(11/12, 23/24), FastLinearFunction(QQ(6), -QQ(5))], [left_open_interval(23/24, QQ(1)), FastLinearFunction(-QQ(18), QQ(18))]])

def gmic_disjoint(f=4/5):
    """
    sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
    sage: h = gmic_disjoint(4/5)
    sage: extremality_test(h, False)
    True
    """        
    pieces = [[right_open_interval(0, f), FastLinearFunction(1/f, 0)],
              [[f, 1], FastLinearFunction(-1/(1-f), 1/(1-f))]]
    return FastPiecewise(pieces, merge=False)

def gmic_disjoint_with_singletons(f=4/5):
    """
    sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
    sage: h = gmic_disjoint_with_singletons(4/5)
    sage: extremality_test(h, False)
    True
    """        
    pieces = [singleton_piece(0, 0), 
              [open_interval(0, f), FastLinearFunction(1/f, 0)],
              [right_open_interval(f, 1), FastLinearFunction(-1/(1-f), 1/(1-f))],
              singleton_piece(1, 0)]
    return FastPiecewise(pieces, merge=False)

def bhk_raises_plotting_error():
    """
    There were some problems when plotting the 2d-diagram of functions having
    irrational data coerced into a (non-quadratic) RealNumberField.

    This was a plotting regression introduced by plotting function and shadows
    at the borders of 2d diagrams.

    TESTS::
    
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = bhk_raises_plotting_error()
        sage: g = plot_2d_diagram(h)
    """
    ## sage: extremality_test(h,True)
    ## ...
    ## verbose 0 (2716: plot.py, generate_plot_points) WARNING: When plotting,
    ## failed to evaluate function at 200 points.
    ## verbose 0 (2716: plot.py, generate_plot_points) Last error message:
    ## 'unsupported operand parent(s) for '-': 'Number Field in a with defining
    ## polynomial y^4 - 4*y^2 + 1' and '<type 'float'>''
    ## ...
    ## TypeError: unsupported operand parent(s) for '+': 'Number Field in a
    ## with defining polynomial y^4 - 4*y^2 + 1' and 'Real Field with 53 bits
    ## of precision'

    ## The problems concern some strange behaviors of RealNumberField_absolute
    ## sage: [x]=nice_field_values([2^(1/3)])
    ## INFO: ... Coerced into real number field: Number
    ## Field in a with defining polynomial y^3 - 2
    ## sage: x+float(1)
    ## TypeError: unsupported operand parent(s) for '+': 'Number Field in a
    ## with defining polynomial y^3 - 2' and '<type 'float'>'
    ## sage: x+RealField(53)(1)
    ## TypeError: unsupported operand parent(s) for '+': 'Number Field in a
    ## with defining polynomial y^3 - 2' and 'Real Field with 53 bits of precision'
    ## sage: RR(x)
    ## RuntimeError: maximum recursion depth exceeded while calling a Python object
    
    ## The following operations work for RealNumberField_absolute
    ## sage: x+RealIntervalField(53)(1)
    ## 2.259921049894873?
    ## sage: x+1
    ## RNF2.259921049894873?
    ## sage: x+1/10
    ## RNF1.359921049894873?
    
    ## When dealing with RealNumberField_quadratic, everything works well.
    ## sage: [y]=nice_field_values([2^(1/2)])
    ## INFO: ... Coerced into real number field: Number
    ## Field in a with defining polynomial y^2 - 2
    ## sage: y+float(1)
    ## 2.414213562373095

    ## Related issues:
    ## sage: K.<s> = NumberField(x^2 - 2)
    ## These operations give error: sage: s+float(1); s+RealField(53)(1); RR(s)

    ## If we define
    ## sage: K.<s> = QuadraticField(2)
    ## or
    ## sage: K.<s> = NumberField(x^2 - 2, embedding=1.4)
    ## then everything goes well.

    ## However,
    ## sage: K.<s> = NumberField(x^2 - 2, embedding=RIF(1.4))
    ## ValueError: 1.4000000000000000? is not a root of the defining polynomial
    ## of Number Field in s with defining polynomial x^2 - 2
    return bhk_irrational(f=4/5, d1=3/5, d2=1/10, a0=15/100, delta=(1/200, sqrt(2)/200, sqrt(3)/200), field=None)

def minimal_has_uncovered_breakpoints():
    """
    sage: logging.disable(logging.INFO)
    sage: h = minimal_has_uncovered_breakpoints()
    sage: finite_dimensional_extremality_test(h,show_all_perturbations = True)
    False
    sage: len(h._perturbations)
    3

    ..NOTE:
    If symbolic function were defined as a piecewise linear function with endpoints of \{covered_intervals\} being its breakpoints,
    then finite_dimensional_extremality_test(h,show_all_perturbations = True)
    would find a solution space that only has dimension 1,
    The finite dimensional perturbation \bar\pi_{\T} would be
    [<FastPiecewise with 4 parts,
    (0, 1/9) <FastLinearFunction x> values: [0, 1/9]
    (1/9, 7/9) <FastLinearFunction -1/3*x + 4/27> values: [1/9, -1/9]
    (7/9, 8/9) <FastLinearFunction x - 8/9> values: [-1/9, 0]
    (8/9, 1) <FastLinearFunction 0> values: [0, 0]>]
    Any equivariant perturbation \bar\pi_{\text{zero}(\T)} vanishes at breakpoint of h.
    Thus, any perturbation \bar\pi generated by extremality_test(h,show_all_perturbations = True) would satisfy \bar\pi(1/3) = (1/3)*\bar\pi(1/9).
    However this does not cover the entire perturbation space, since perturbtions \bar\pi such as
    <FastPiecewise with 5 parts,
    (0, 2/9) <FastLinearFunction 0> values: [0, 0]
    (2/9, 1/3) <FastLinearFunction x - 2/9> values: [0, 1/9]
    (1/3, 5/9) <FastLinearFunction -x + 4/9> values: [1/9, -1/9]
    (5/9, 2/3) <FastLinearFunction x - 2/3> values: [-1/9, 0]
    (2/3, 1) <FastLinearFunction 0> values: [0, 0]>
    where \bar\pi(1/3) \neq (1/3)*\bar\pi(1/9) will be missing.
    """
    bkpts = [0, 1/9, 2/9, 3/9, 5/9, 6/9, 7/9, 8/9, 1]
    values = [0, 5/8, 5/8, 4/8, 4/8, 3/8, 3/8, 1, 0]
    h = piecewise_function_from_breakpoints_and_values(bkpts, values)
    return h

def bhk_discontinuous_irrational(f=4/5, d1=3/5, d2=1/10, a0=15/100, delta=(1/199, sqrt(2)/199), shift_scale=1, field=None):
    """
    There exists crazy fat perturbation \tilde\pi in this example:
    \tilde\pi(x)
        = epsion if x \in (3/10, 7/20) such that x - 3/10 \in <t_a, t_b>_Z;
        = -epsion if x \in (3/10, 7/20) such that x - 7/20 \in <t_a, t_b>_Z;
        = epsion if x \in (9/20, 1/2) such that x - 9/20 \in <t_a, t_b>_Z;
        = -epsion if x \in (9/20, 1/2) such that x - 9/20 \in <t_a, t_b>_Z;
        = 0 otherwise,
    where t_a = delta[0], t_b = delta[0]+delta[1] and epsilon is a very small positive number, for example epsilon=0.001.
    (requires that d2/2 is not in <t_a, t_b>_Z;)
    """
    h_pwl = bhk_irrational(f=f, d1=d1, d2=d2, a0=a0, delta=delta, field=field)
    bkpt = h_pwl.end_points()
    v = h_pwl.values_at_end_points()
    f = bkpt[-2]
    c1= v[1]/bkpt[1]
    c3= 1/(f-1)
    n = len(delta)
    b = min([bkpt[i+2]-bkpt[i+1] for i in range(2*n+2)]) * shift_scale / 2
    jump = b*(c1-c3)
    disc_pts = [bkpt[2*i] for i in range(n+2)]+\
               [bkpt[-2*i] for i in range(n+2,0,-1)]+[bkpt[-1]]
    discp = []
    for i in range(len(disc_pts)-1):
        discp.append(singleton_piece(disc_pts[i],0))
        if i < n+2:
            j = jump
        else:
            j = -jump
        discp.append(open_piece((disc_pts[i],j), (disc_pts[i]+b,0)))
        discp.append(singleton_piece(disc_pts[i]+b,j))
        discp.append(open_piece((disc_pts[i]+b,0), (disc_pts[i+1]-b,0)))
        if (i+1 < n+2) or (i+1 == len(disc_pts)-1):
            j = jump
        else:
            j = -jump
        discp.append(singleton_piece(disc_pts[i+1]-b,j))
        discp.append(open_piece((disc_pts[i+1]-b,0), (disc_pts[i+1],j)))
    discp.append(singleton_piece(disc_pts[-1], 0))
    h_shift = FastPiecewise(discp)
    return h_pwl+h_shift


def discontinuous_bhk_irrational_dense_move_not_affine():
    """
    This example shows that when \pi is two sided discontinuous at the origin,
    the perturbation function \tilde\pi may be not affine linear on the dense intervals from the Strip Lemma.

    Lemma 9.7 implies \pi is affine linear on the intervals (3/10, 7/20) and (9/20, 1/2) where the moves are dense. We also know that \tilde\pi is affine linear over each coset of <t_a, t_b>_Z on these two intervals, where t_a, t_b are irrational translation moves satifying the hypotheses of the Strip Lemma.
    But as \pi is two-sided discontinuous at 0, Corollary 8.3 (or even Lemma 8.2) doesn't apply, thus \tilde\pi could be discontinuous on the dense intervals.
    On the 2d diagram, each face is either additive or strictly subadditive (every vertice of the face $F$ is strictly subaddtive). Then the case \delta\pi_F(x,y) > 0 and S != empty never happens in the proof of Lemma 8.7. Hence any bounded \bar\bi is an effective perturbation.

    There exists crazy fat perturbation \tilde\pi in this example:
    \tilde\pi(x)
        = epsion if x \in (3/10, 7/20) such that x - 3/10 \in <t_a, t_b>_Z;
        = -epsion if x \in (3/10, 7/20) such that x - 7/20 \in <t_a, t_b>_Z;
        = epsion if x \in (9/20, 1/2) such that x - 9/20 \in <t_a, t_b>_Z;
        = -epsion if x \in (9/20, 1/2) such that x - 9/20 \in <t_a, t_b>_Z;
        = 0 otherwise,
    where epsion is a very small positive number. 
    (requires that d2/2 is not in <t_a, t_b>_Z;)

    EXAMPLE:
    sage: logging.disable(logging.INFO)
    sage: dbhk =  discontinuous_bhk_irrational_dense_move_not_affine()
    sage: plot_with_colored_slopes(dbhk).show(figsize=20) #not tested
    sage: plot_2d_diagram(dbhk).show(figsize=20) #not tested
    """
    bhk = bhk_irrational(delta=(1/199, sqrt(2)/199))
    b=bhk.end_points()
    v=bhk.values_at_end_points()
    vm = v[0:8]+ [v[8]+1/800, v[9]+1/800, v[10]-1/800, v[11]-1/800] + v[12::]
    vl = [v[0]+1/500, v[1], v[2]+1/500, v[3]-1/500, v[4]+1/500, v[5]-1/500, v[6]+1/500, v[7], v[8], v[9]+1/800, v[10], v[11]-1/800, v[12], \
          v[13]-1/500,v[14]+1/500,v[15]-1/500,v[16]+1/500,v[17]-1/500,v[18],v[19]-1/500]
    vr = [v[1], v[2]+1/500, v[3]-1/500, v[4]+1/500, v[5]-1/500, v[6]+1/500, v[7], v[8]+1/800, v[9], v[10]-1/800, v[11], v[12], \
          v[13]-1/500,v[14]+1/500,v[15]-1/500,v[16]+1/500,v[17]-1/500,v[18],v[19]-1/500,v[20]+1/500]
    discp = []
    for i in range(20):
        discp.append(singleton_piece(b[i],vm[i]))
        discp.append(open_piece((b[i],vl[i]), (b[i+1],vr[i])))
    discp.append(singleton_piece(b[20],vm[20]))
    dbhk = FastPiecewise(discp)
    return dbhk
