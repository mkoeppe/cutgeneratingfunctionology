def extreme_5slope_no_0d_1d_1():
    """
    5-slope extreme function without any 0-d or 1-d maximal additive faces
    except for the symmetry reflection or x=0 or y=0.

    LP model: 5slope_q37_f25_m12_fulldim.lp
    q = 37; f = 25/37; given number of slopes = 5; maxstep = 1;
    slope gap >= q/12; subadditive slack >= q/12;
    obj = max_slope_slack. MIPFocus = 1

    Gurobi run on my laptop,
    Interrupt at the 2nd feasible solution (662s, obj = 22.81667)
    Solution file: solution_5slope_no_0d_1d_1.sol
    vertex function is a 5-slope extreme function

    This function is obtained by:
    sage: faces, fn = painted_faces_and_funciton_from_solution('.../solution_5slope_no_0d_1d_1.sol', 37, showplots=False)
    sage: h_list = investigate_faces_solution(37, 25/37, faces)
    sage: extreme_5slope_no_0d_1d_1 = h_list[0]

    EXAMPLES::

        sage: h = extreme_5slope_no_0d_1d_1()
        sage: extremality_test(h)
        True
    """
    bkpt = [0, 2/37, 3/37, 5/37, 7/37, 10/37, 11/37, 12/37, 13/37, 14/37, 15/37, \
            18/37, 20/37, 22/37, 23/37, 25/37, 27/37, 28/37, 29/37, 33/37, 34/37, 35/37, 1]
    values = [0, 59/90, 7/9, 64/135, 4/9, 73/90, 43/54, 29/45, 16/45, 11/54, 17/90, \
              5/9, 71/135, 2/9, 31/90, 1, 19/45, 73/270, 23/90, 67/90, 197/270, 26/45, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def extreme_5slope_no_0d_1d_2():
    """
    5-slope extreme function without any 0-d or 1-d maximal additive faces
    except for the symmetry reflection or x=0 or y=0.

    LP model: 5slope_q37_f25_m12_fulldim.lp
    q = 37; f = 25/37; given number of slopes = 5; maxstep = 1;
    slope gap >= q/12; subadditive slack >= q/12;
    obj = max_slope_slack. MIPFocus = 1

    Gurobi run on my laptop,
    Interrupt at the 4nd feasible solution (1285s, obj = 28.08205)
    Solution file: solution_5slope_no_0d_1d_2.sol
    vertex function is a 5-slope extreme function

    This function is obtained by:
    sage: faces, fn = painted_faces_and_funciton_from_solution('.../solution_5slope_no_0d_1d_2.sol', 37, showplots=False)
    sage: h_list = investigate_faces_solution(37, 25/37, faces)
    sage: extreme_5slope_no_0d_1d_2 = h_list[0]

    EXAMPLES::

        sage: h = extreme_5slope_no_0d_1d_2()
        sage: extremality_test(h)
        True
    """
    bkpt = [0, 1/37, 2/37, 5/37, 7/37, 10/37, 12/37, 13/37, 15/37, 18/37, 20/37, \
            23/37, 24/37, 25/37, 26/37, 27/37, 28/37, 29/37, 33/37, 34/37, 35/37, 36/37, 1]
    values = [0, 8/15, 4/13, 10/13, 47/78, 152/195, 239/390, 151/390, 43/195, 31/78, \
              3/13, 9/13, 7/15, 1, 151/195, 539/780, 121/260, 149/390, 241/390, 139/260, 241/780, 44/195, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def extreme_5slope_no_0d_1d_3():
    """
    5-slope extreme function without any 0-d or 1-d maximal additive faces
    except for the symmetry reflection or x=0 or y=0.

    LP model: 5slope_q37_f25_m12_fulldim.lp
    Solution file: solution_5slope_no_0d_1d_3and4.sol
    the first vertex function is a 5-slope extreme function
    """
    bkpt = [0, 1/37, 3/37, 5/37, 7/37, 10/37, 12/37, 13/37, 15/37, 18/37, 20/37, \
            22/37, 24/37, 25/37, 26/37, 27/37, 28/37, 29/37, 33/37, 34/37, 35/37, 36/37, 1]
    values = [0, 8/15, 11/30, 49/60, 13/20, 77/100, 181/300, 119/300, 23/100, 7/20, \
              11/60, 19/30, 7/15, 1, 119/150, 71/100, 151/300, 21/50, 29/50, 149/300, 29/100, 31/150, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def extreme_5slope_no_0d_1d_4():
    """
    5-slope extreme function without any 0-d or 1-d maximal additive faces
    except for the symmetry reflection or x=0 or y=0.

    LP model: 5slope_q37_f25_m12_fulldim.lp
    Solution file: solution_5slope_no_0d_1d_3and4.sol
    the second vertex function is a 5-slope extreme function
    """
    bkpt = [0, 1/37, 3/37, 5/37, 7/37, 10/37, 12/37, 13/37, 15/37, 18/37, 20/37, \
            22/37, 24/37, 25/37, 26/37, 27/37, 28/37, 29/37, 33/37, 34/37, 35/37, 36/37, 1]
    values = [0, 41/63, 61/126, 191/252, 149/252, 197/252, 155/252, 97/252, 55/252, \
              103/252, 61/252, 65/126, 22/63, 1, 97/126, 173/252, 115/252, 47/126, 79/126, 137/252, 79/252, 29/126, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def extreme_5slope_no_0d_1d_5():
    """
    5-slope extreme function without any 0-d or 1-d maximal additive faces
    except for the symmetry reflection or x=0 or y=0.

    LP model: 5slope_q37_f25_m12_fulldim.lp
    Solution file: solution_5slope_no_0d_1d_5.sol
    vertex function is a 5-slope extreme function
    """
    bkpt = [0, 1/37, 3/37, 5/37, 7/37, 10/37, 12/37, 13/37, 15/37, 18/37, 20/37, \
            22/37, 24/37, 25/37, 28/37, 29/37, 33/37, 34/37, 1]
    values = [0, 145/221, 6/13, 10/13, 127/221, 347/442, 261/442, 181/442, 95/442, \
              94/221, 3/13, 7/13, 76/221, 1, 101/221, 159/442, 283/442, 120/221, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def extreme_5slope_no_transrefl():
    """
    5-slope extreme function without any translation/reflection
    except for the symmetry reflection.

    LP model: 5slope_q22_f10_step1_m4.lp
    q = 22; f = 10/22; given number of slopes = 5; maxstep = 1;
    slope gap >= q/4; subadditive slack >= q/4;
    obj = max_slope_slack.

    Gurobi run on my laptop, (51.83s, obj = 30.25)
    Optimal solution: solution_5slope_no_transrefl.sol
    vertex function is a 5-slope extreme function

    This function is obtained by:
    sage: faces, fn = painted_faces_and_funciton_from_solution('.../solution_5slope_no_transrefl.sol', 22, showplots=False)
    sage: h_list = investigate_faces_solution(22, 10/22, faces)
    sage: new_5slope_1 = h_list[0]

    EXAMPLES::

        sage: h = extreme_5slope_no_transrefl()
        sage: extremality_test(h)
        True
        sage: uncovered_intervals_from_covered_intervals(generate_directly_covered_intervals(h))
        []
    """
    bkpt = [0, 1/22, 1/11, 3/22, 2/11, 3/11, 7/22, 4/11, 9/22, 5/11,\
             1/2, 6/11, 13/22, 7/11, 9/11, 19/22, 10/11, 21/22, 1]
    values = [0, 23/32, 3/4, 7/16, 13/16, 3/16, 9/16, 1/4, 9/32, 1, \
                11/32, 3/8, 3/4, 7/16, 9/16, 1/4, 5/8, 21/32, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def extreme_5slope_no_transrefl_2():
    """
    5-slope extreme function without any translation/reflection
    except for the symmetry reflection.

    LP model: 6slope_q39_f16_step1_m12.lp
    Solution file: solution_5slope_transrefl_2.sol
    The third vertex function is a 5-slope extreme function
    """
    bkpt = [0, 1/39, 1/13, 4/39, 2/13, 7/39, 3/13, 10/39, 4/13, 1/3, 5/13, 16/39, \
            17/39, 7/13, 22/39, 8/13, 25/39, 9/13, 28/39, 10/13, 31/39, 11/13, 34/39, 38/39, 1]
    values = [0, 7/8, 11/24, 19/24, 3/8, 17/24, 7/24, 5/8, 5/24, 13/24, 1/8, 1, 1/4, \
              1/2, 5/6, 5/12, 3/4, 1/3, 2/3, 1/4, 7/12, 1/6, 1/2, 3/4, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

destdir = "/media/sf_dropbox/2q_mip/"

def fn_variable(q, x):
    return 'fn_%s' % int(x*q)

def face_variable(q, face):
    """
    EXAMPLES::
    
        sage: q=7;
        sage: face_variable(q,  Face(([1/7, 2/7], [1/7], [2/7, 3/7])))
        'h_1_1'
        sage: face_variable(q,  Face(([1/7], [1/7, 2/7], [2/7, 3/7])))
        'v_1_1'
        sage: face_variable(q,  Face(([1/7, 2/7], [1/7, 2/7], [3/7])))
        'd_1_1'
        sage: face_variable(q,  Face(([1/7, 2/7], [1/7, 2/7], [2/7, 3/7])))
        'l_1_1'
        sage: face_variable(q,  Face(([1/7, 2/7], [1/7, 2/7], [3/7, 4/7])))
        'u_1_1'
        sage: face_variable(q,  Face(([2/7], [2/7], [4/7])))
        'p_2_2'
    """
    if face.is_0D():
        return 'p_%s_%s' % (int(face.vertices[0][0]*q), int(face.vertices[0][1]*q))
    elif face.is_2D():
        if face.minimal_triple[2][0] == face.minimal_triple[0][0] + face.minimal_triple[1][0]:
            l_or_u = 'l'
        else:
            l_or_u = 'u'
        return '%s_%s_%s' %(l_or_u, int(face.minimal_triple[0][0]*q), int(face.minimal_triple[1][0]*q))
    elif face.is_horizontal():
        return 'h_%s_%s' %(int(face.minimal_triple[0][0]*q), int(face.minimal_triple[1][0]*q))
    elif face.is_vertical():
        return 'v_%s_%s' %(int(face.minimal_triple[0][0]*q), int(face.minimal_triple[1][0]*q))
    elif face.is_diagonal():
        return 'd_%s_%s' %(int(face.minimal_triple[0][0]*q), int(face.minimal_triple[1][0]*q))
         
def variable_face(q, s):
    """
    EXAMPLES::

        sage: q=7;
        sage: variable_face(q, 'h_1_1')
        <Face ([1/7, 2/7], [1/7], [2/7, 3/7])>
        sage: variable_face(q, 'v_1_1')
        <Face ([1/7], [1/7, 2/7], [2/7, 3/7])>
        sage: variable_face(q, 'd_1_1')
        <Face ([1/7, 2/7], [1/7, 2/7], [3/7])>
        sage: variable_face(q, 'l_1_1')
        <Face ([1/7, 2/7], [1/7, 2/7], [2/7, 3/7])>
        sage: variable_face(q, 'u_1_1')
        <Face ([1/7, 2/7], [1/7, 2/7], [3/7, 4/7])>
        sage: variable_face(q, 'p_2_2')
        <Face ([2/7], [2/7], [4/7])>
    """
    assert s[0] in set(['h','v','d','l','u','p']), "varialbe %s is not a face" % s
    index1 = s.find('_')
    index2 = s.rfind('_')
    x = int(s[index1+1: index2])/q
    y = int(s[index2+1::])/q
    if s[0] == 'h':
        return Face(([x, x+1/q], [y], [x+y, x+y+1/q]))
    elif s[0] == 'v':
        return Face(([x], [y, y+1/q], [x+y, x+y+1/q]))
    elif s[0] == 'd':
        return Face(([x, x+1/q], [y, y+1/q], [x+y+1/q]))
    elif s[0] == 'l':
        return Face(([x, x+1/q], [y, y+1/q], [x+y, x+y+1/q]))
    elif s[0] == 'u':
        return Face(([x, x+1/q], [y, y+1/q], [x+y+1/q, x+y+2/q]))
    elif s[0] == 'p':
        return Face(([x], [y], [x+y]))

def vertex_variable(q, v):
    """
    EXAMPLES::

        sage: vertex_variable(7, (1/7, 2/7))
        'p_1_2'
    """
    return 'p_%s_%s' % (int(v[0]*q), int(v[1]*q))

def variable_vertex(q, s):
    """
    EXAMPLES::

        sage: variable_vertex(7, 'p_1_2')
        (1/7, 2/7)
    """
    assert s[0] == 'p', "varialbe %s is not a vertex" % s
    index1 = s.find('_')
    index2 = s.rfind('_')
    x = int(s[index1+1: index2])/q
    y = int(s[index2+1::])/q
    return (x, y)

def print_logical_constraints(filename, q, face):
    """
    EXAMPLES::

        sage: import sys;
        sage: print_logical_constraints(sys.stdout, 7, Face(([1/7, 2/7], [1/7, 2/7], [2/7, 3/7])))
        l_1_1 - p_1_1 >= 0
        l_1_1 - p_1_2 >= 0
        l_1_1 - p_2_1 >= 0
        l_1_1 - p_1_1 - p_1_2 - p_2_1 <= 0     
    """
    for v in face.vertices:
        print >> filename, '%s - %s >= 0' %(face_variable(q, face), vertex_variable(q, v))
    v = face.vertices[0]
    print >> filename, '%s' % face_variable(q, face),
    for v in face.vertices:
        print >> filename, '- %s' % vertex_variable(q, v),
    print >> filename, '<= 0'
   
def print_xy_swapped_constraints(filename, q, face):
    """
    EXAMPLES::

        sage: import sys;
        sage: face =Face(([1/7], [2/7], [3/7]))
        sage: print_xy_swapped_constraints(sys.stdout, 7, face)
        p_1_2 - p_2_1 = 0
    """
    print >> filename, '%s - %s = 0' %(face_variable(q, face), face_variable(q, x_y_swapped_face(face)))


def print_fn_bounds(filename, q):
    """
    EXAMPLES::

        sage: print_fn_bounds(sys.stdout, 3)
        0 <= fn_0 <= 1
        0 <= fn_1 <= 1
        0 <= fn_2 <= 1
        0 <= fn_3 <= 1
    """
    bkpt = [x/q for x in range(q+1)]
    for x in bkpt:
        print >> filename, '0 <= %s <= 1' % fn_variable(q, x)

def print_fn_minimality_test(filename, q, f, m=0):
    """
    EXAMPLES::

        sage: print_fn_minimality_test(sys.stdout, 3, 2/3)
        fn_0 = 0
        fn_0 + fn_2 = 1
        fn_1 + fn_1 = 1
        fn_2 + fn_3 = 1
        fn_1 + fn_1 - fn_2 >= 0
        fn_1 + fn_1 - fn_2 - 2 p_1_1 <= 0
        fn_1 + fn_2 - fn_3 >= 0
        fn_1 + fn_2 - fn_3 - 2 p_1_2 <= 0
        fn_2 + fn_2 - fn_1 >= 0
        fn_2 + fn_2 - fn_1 - 2 p_2_2 <= 0
    """
    # fn(0) = 0
    print >> filename, '%s = 0' % fn_variable(q, 0)
    # fn(f) = 1
    #print >> filename, '%s = 1' % fn_variable(q, f)
    bkpt = [x/q for x in range(q+1)]
    # symmetric conditions
    x = 0
    while x <= f/2:
        print >> filename, '%s + %s = 1' % (fn_variable(q, x), fn_variable(q, f - x))
        x += 1/q
    x = f
    while x <= (1+f)/2:
        print >> filename, '%s + %s = 1' % (fn_variable(q, x), fn_variable(q, 1 + f - x))
        x += 1/q 
    # strict-subadditivity and additivity conditions
    for i in range(1, q):
        for j in range(i, q):
            x = bkpt[i]
            y = bkpt[j]
            z = fractional(x + y)
            if m == 0:
                print >> filename, '%s + %s - %s >= 0' %(fn_variable(q, x), fn_variable(q, y), fn_variable(q, z))
            else:
                print >> filename, '%s + %s - %s - %s %s >= 0' %(fn_variable(q, x), fn_variable(q, y), fn_variable(q, z), \
                                                                 RR(1 / m), vertex_variable(q, (x, y)))
            print >> filename, '%s + %s - %s - 2 %s <= 0' %(fn_variable(q, x), fn_variable(q, y), \
                                                            fn_variable(q, z), vertex_variable(q, (x, y)))

def print_trivial_additive_points(filename, q, f):
    """
    EXAMPLES::

        sage: print_trivial_additive_points(sys.stdout, 3, 2/3)
        p_0_0 = 0
        p_0_1 = 0
        p_0_2 = 0
        p_0_3 = 0
        p_1_3 = 0
        p_2_3 = 0
        p_3_3 = 0
        p_0_2 = 0
        p_1_1 = 0
        p_2_0 = 0
        p_2_3 = 0
        p_3_2 = 0
    """
    bkpt = [x/q for x in range(q+1)]
    # border x = 0 and border y = 0 are green
    for x in bkpt:
        print >> filename, '%s = 0' % vertex_variable(q, (0, x))
    for x in bkpt[1::]:
        print >> filename, '%s = 0' % vertex_variable(q, (x, 1))
    # diagonals corresponding to f
    for x in bkpt:
        if x < f:
            print >> filename, '%s = 0' % vertex_variable(q, (x, f - x))
        elif x == f:
            print >> filename, '%s = 0' % vertex_variable(q, (x, f - x))
            print >> filename, '%s = 0' % vertex_variable(q, (x, f - x + 1))
        elif x > f:
            print >> filename, '%s = 0' % vertex_variable(q, (x, f - x + 1))

    #b = f - a
    #print >> filename, '%s = 0' % vertex_variable(q, (b - a + 1/q, a - 1/q))
    #print >> filename, '%s = 0' % vertex_variable(q, (b - a + 1/q, a))
    
def covered_interval_variable(q, x):
    """
    EXAMPLES::

        sage: covered_interval_variable(7, 3/7)
        'c_3'
    """
    return 'c_%s' % int(x * q)

def variable_covered_interval(q, s):
    """
    EXAMPLES::

        sage: variable_covered_interval(7, 'c_3')
        [3/7, 4/7]
    """
    assert s[0] == 'c', "varialbe %s is not a covered interval" % s
    index1 = s.rfind('_')
    x = int(s[index1+1::])/q
    return [x, x + 1/q]

def covered_i_variable(q, z, i):
    """
    EXAMPLES::

        sage: covered_i_variable(7, 3/7, 2)
        'c_3_2'
    """
    return 'c_%s_%s' % (int(z * q), i)

def translation_i_variable(q, x, z, i):
    """
    EXAMPLES::

        sage: translation_i_variable(7, 1/7, 3/7, 2)
        't_1_3_2'
    """
    return 't_%s_%s_%s' % (int(x * q), int(z * q), i)

def reflection_i_variable(q, x, z, i):
    """
    EXAMPLES::

        sage: reflection_i_variable(7, 1/7, 3/7, 2)
        'r_1_3_2'
    """   
    return 'r_%s_%s_%s' % (int(x * q), int(z * q), i)

def print_directly_covered_constraints(filename, q, z):
    """
    EXAMPLES::

        sage: print_directly_covered_constraints(sys.stdout, 3, 1/3)
        c_1_0 - l_1_0 <= 0
        c_1_0 - u_1_0 <= 0
        c_1_0 - l_1_1 <= 0
        c_1_0 - u_1_1 <= 0
        c_1_0 - l_1_2 <= 0
        c_1_0 - u_1_2 <= 0
        c_1_0 - l_0_1 <= 0
        c_1_0 - l_1_0 <= 0
        c_1_0 - l_2_2 <= 0
        c_1_0 - u_0_0 <= 0
        c_1_0 - u_1_2 <= 0
        c_1_0 - u_2_1 <= 0
        c_1_0 - l_1_0 - u_1_0 - l_1_1 - u_1_1 - l_1_2 - u_1_2 - l_0_1 - l_1_0 - l_2_2 - u_0_0 - u_1_2 - u_2_1 >= -11
    """
    bkpt = [x/q for x in range(q)]
    c_z_0 = covered_i_variable(q, z, 0)
    variable_list = []
    for y in bkpt:
        # I projection: l_zz,yy and u_zz,yy
        l = face_variable(q, Face(([z, z+1/q], [y, y+1/q], [z + y, z + y + 1/q])))
        print >> filename, '%s - %s <= 0' % (c_z_0, l)
        u = face_variable(q, Face(([z, z+1/q], [y, y+1/q], [z + y + 1/q, z + y + 2/q])))
        print >> filename, '%s - %s <= 0' % (c_z_0, u)
        variable_list += [l, u]
    for x in bkpt:
        # K projection: l_xx,zz-xx
        y = z - x
        if y < 0:
            y += 1
        l = face_variable(q, Face(([x, x+1/q], [y, y+1/q], [x + y, x + y + 1/q])))
        print >> filename, '%s - %s <= 0' % (c_z_0, l)
        variable_list += [l]
    for x in bkpt:
        # K projection: u_xx,zz-xx-1
        y = z - x - 1/q
        if y < 0:
            y += 1
        u = face_variable(q, Face(([x, x+1/q], [y, y+1/q], [x + y + 1/q, x + y + 2/q])))
        print >> filename, '%s - %s <= 0' % (c_z_0, u)
        variable_list += [u]
    assert len(variable_list) == 4*q
    print >> filename, '%s' % c_z_0,
    for lu in variable_list:
        print >> filename, '- %s' % lu,
    print >> filename, '>= %s' % (1 - 4*q)

    ## want: I,J projections don't intersect with each other
    ## problem: this model has no feasible solution.
    #for lu in variable_list[0: 2*q]:
    #    print >> filename, '+ %s' % lu,
    #print >> filename, '>= %s' % (2*q -1)

def print_translation_i_constraints(filename, q, x, z, i):
    """
    EXAMPLES::

        sage: print_translation_i_constraints(sys.stdout, 7, 1/7, 3/7, 2)
        c_1_1 - t_1_3_2 <= 0
        h_1_2 - t_1_3_2 <= 0
        t_1_3_2 - c_1_1 - h_1_2 <= 0
        sage: print_translation_i_constraints(sys.stdout, 7, 3/7, 1/7, 2)
        c_3_1 - t_3_1_2 <= 0
        h_1_2 - t_3_1_2 <= 0
        t_3_1_2 - c_3_1 - h_1_2 <= 0
    """
    c_x_last = covered_i_variable(q, x, i-1)
    t_x_z_i =  translation_i_variable(q, x, z, i)
    print >> filename, '%s - %s <= 0' % (c_x_last, t_x_z_i)
    if x <= z:
        move = Face(([x, x + 1/q], [z - x], [z, z + 1/q])) # h_xx,zz-xx
    else:
        move = Face(([z, z + 1/q], [x - z], [x, x + 1/q])) # h_zz,xx-zz
    print >> filename, '%s - %s <= 0' % (face_variable(q, move), t_x_z_i)
    print >> filename, '%s - %s - %s <= 0' % (t_x_z_i, c_x_last, face_variable(q, move))
    
def print_reflection_i_constraints(filename, q, x, z, i):
    """
    EXAMPLES::

        sage: print_reflection_i_constraints(sys.stdout, 7, 1/7, 3/7, 2)
        c_1_1 - r_1_3_2 <= 0
        d_1_3 - r_1_3_2 <= 0
        r_1_3_2 - c_1_1 - d_1_3 <= 0
    """
    c_x_last = covered_i_variable(q, x, i-1)
    r_x_z_i =  reflection_i_variable(q, x, z, i)
    print >> filename, '%s - %s <= 0' % (c_x_last, r_x_z_i)
    move = Face(([x, x + 1/q], [z, z + 1/q], [x + z + 1/q])) # d_xx,zz
    print >> filename, '%s - %s <= 0' % (face_variable(q, move), r_x_z_i)
    print >> filename, '%s - %s - %s <= 0' % (r_x_z_i, c_x_last, face_variable(q, move))
    
def print_undirectly_covered_i_constraints(filename, q, z, i):
    """
    EXAMPLES::

        sage: print_undirectly_covered_i_constraints(sys.stdout, 3, 1/3, 2)
        c_1_2 - c_1_1 <= 0
        c_1_2 - t_0_1_2 <= 0
        c_1_2 - r_0_1_2 <= 0
        c_1_2 - t_2_1_2 <= 0
        c_1_2 - r_2_1_2 <= 0
        c_1_2 - c_1_1 - t_0_1_2 - r_0_1_2 - t_2_1_2 - r_2_1_2 >= -4
    """
    bkpt = [x/q for x in range(q)]
    c_z_now = covered_i_variable(q, z, i)
    c_z_last = covered_i_variable(q, z, i-1)
    print >> filename, '%s - %s <= 0' % (c_z_now, c_z_last)
    variable_list = [c_z_last]
    for x in bkpt:
        if x != z:
            t_x_z_i = translation_i_variable(q, x, z, i)
            r_x_z_i = reflection_i_variable(q, x, z, i)
            print >> filename, '%s - %s <= 0' % (c_z_now, t_x_z_i)
            print >> filename, '%s - %s <= 0' % (c_z_now, r_x_z_i)
            variable_list += [t_x_z_i, r_x_z_i]
    assert len(variable_list) == 2 * q - 1
    print >> filename, '%s' % c_z_now,
    for v in variable_list:
         print >> filename, '- %s' % v,
    print >> filename, '>= %s' % (2 - 2 * q)     

def print_obj_max_slope_slack(filename, nums):
    """
    slope_slack = s_0 - s_(nums-1)

    EXAMPLES::

        sage: print_obj_max_slope_slack(sys.stdout, 5)
        s_0 - s_4
    """
    print >> filename, 's_0 - s_%s' % (nums - 1),

def print_obj_max_subadd_slack(filename, q, weight=1): #is a constant!
    """
    subadd_slack = q * (\sum_x fn(x)) 
                 = q * (q - 1)

    EXAMPLES::

        sage: print_obj_max_subadd_slack(sys.stdout, 3)
        1 fn_0 + 1 fn_1 + 1 fn_2 
    """
    bkpt = [x/q for x in range(q)]
    print >> filename, '%s %s' % (weight, fn_variable(q, bkpt[0])),
    for x in bkpt[1::]:
        print >> filename, '+ %s %s' % (weight, fn_variable(q, x)),

def print_obj_min_undirectly_covered_times(filename, q, step=None, weight=1):
    """
    EXAMPLES::

        sage: print_obj_min_undirectly_covered_times(sys.stdout, 2)
        + 1 t_0_1_1 + 1 r_0_1_1 + 1 t_1_0_1 + 1 r_1_0_1
    """
    if step is None:
        step = q - 1
    bkpt = [x/q for x in range(q)]
    if step > 0:
        for x in bkpt:
            for z in bkpt:
                if x != z:
                    print >> filename, '+ %s %s' % (weight, translation_i_variable(q, x, z, step)),
                    print >> filename, '+ %s %s' % (weight, reflection_i_variable(q, x, z, step)),

def print_obj_min_covered_times_max_subadd_slack(filename, q, maxstep=None):
    """
    EXAMPLES::

        sage: print_obj_min_covered_times_max_subadd_slack(sys.stdout, 2)
        1 fn_0 + 1 fn_1 + 1 l_0_0 + 1 u_0_0 + 1 l_0_1 + 1 u_0_1 + 1 l_1_0 + 1 u_1_0 + 1 l_1_1 + 1 u_1_1 + 
        1 t_0_1_1 + 1 r_0_1_1 + 1 t_1_0_1 + 1 r_1_0_1
    """
    if maxstep is None:
        maxstep = q
    print_obj_max_subadd_slack(filename, q, weight = 1) # should weight q instead of 1.
    print_obj_min_directly_covered_times(filename, q, weight = 1)
    print_obj_min_undirectly_covered_times(filename, q, step = maxstep - 1, weight = 1)

def print_obj_5slope22(filename, q, weight=1):
    h = piecewise_function_from_robert_txt_file("/media/sf_dropbox/data/example5Slope22data.txt")
    bkpt = [x/q for x in range(q + 1)]
    m = 0
    for x in bkpt:
        for y in bkpt:
            if x <= y and h(x) + h(y) != h(fractional(x + y)):
                print >> filename, '+ %s %s' % ( weight, vertex_variable(q, (x, y)) ),
                #m += 1
    #print m

def print_obj_min_add_points(filename, q, weight=1):
    bkpt = [x/q for x in range(q + 1)]
    for x in bkpt:
        for y in bkpt:
            if x <= y:
                print >> filename, '+ %s %s' % ( weight, vertex_variable(q, (x, y)) ),

def print_obj_min_directly_covered_times(filename, q, weight=1):
    """
    EXAMPLES::

        sage: print_obj_min_directly_covered_times(sys.stdout, 2)
        + 1 l_0_0 + 1 u_0_0 + 1 l_0_1 + 1 u_0_1 + 1 l_1_0 + 1 u_1_0 + 1 l_1_1 + 1 u_1_1 
    """
    bkpt = [x/q for x in range(q)]
    for x in bkpt:
        for y in bkpt:
            print >> filename, '+ %s %s + %s %s' % ( weight, face_variable(q, Face(([x, x+1/q], [y, y+1/q], [x + y, x + y + 1/q]))), \
                                                     weight, face_variable(q, Face(([x, x+1/q], [y, y+1/q], [x + y + 1/q, x + y + 2/q]))) ),
def all_faces(q):
    faces_2d = []
    faces_diag = []
    faces_hor = []
    faces_ver = []
    faces_0d = []
    for xx in range(q):
        for yy in range(q):
            faces_2d.append( Face(([xx/q, (xx+1)/q], [yy/q, (yy+1)/q], [(xx+yy)/q, (xx+yy+1)/q])) )
            faces_2d.append( Face(([xx/q, (xx+1)/q], [yy/q, (yy+1)/q], [(xx+yy+1)/q, (xx+yy+2)/q])) )
            faces_diag.append( Face(([xx/q, (xx+1)/q], [yy/q, (yy+1)/q], [(xx+yy+1)/q])) )

    for xx in range(q):
        for yy in range(q+1):
            faces_hor.append( Face(([xx/q, (xx+1)/q], [yy/q], [(xx+yy)/q, (xx+yy+1)/q])) )
            faces_ver.append( Face(([yy/q], [xx/q, (xx+1)/q], [(xx+yy)/q, (xx+yy+1)/q])) )

    for xx in range(q+1):
        for yy in range(q+1):
            faces_0d.append( Face(([xx/q], [yy/q], [(xx+yy)/q])) )
    return faces_2d, faces_diag, faces_hor, faces_ver, faces_0d

def write_lpfile(q, f, nums, maxstep=None, m=0):
    """
    EXAMPLES:

        sage: write_lpfile(22, 10/22, 5, 1, 4)
    """
    if maxstep is None:
        maxstep = q
    filename = open(destdir + "%sslope_q%s_f%s_m%s_fulldim.lp" % (nums, q, int(f*q), m), "w")

    faces_2d, faces_diag, faces_hor, faces_ver, faces_0d = all_faces(q)

    print >> filename, '\ MIP model with q = %s, f = %s, num of slopes = %s, small_m = %s, without non-trivial 0d and 1d maximal faces' % (q, f, nums, m)

    print >> filename, 'Maximize'
    #print >> filename, 0
    print_obj_max_slope_slack(filename, nums)
    #print_obj_max_subadd_slack(filename, q) # is a constant!
    #print_obj_min_directly_covered_times(filename, q)
    #print_obj_min_undirectly_covered_times(filename, q)
    #print_obj_min_covered_times_max_subadd_slack(filename, q, maxstep=maxstep)
    #print_obj_5slope22(filename, q, weight=1)
    #print_obj_min_add_points(filename, q, weight=1)
    print >> filename

    print >> filename, 'Subject to'
    for face in faces_2d + faces_diag + faces_hor + faces_ver:
        #if face.minimal_triple[0][0] <= face.minimal_triple[1][0]:
            print_logical_constraints(filename, q, face)
    
    for face in faces_0d:
        if face.minimal_triple[0][0] < face.minimal_triple[1][0]:
            print_xy_swapped_constraints(filename, q, face)

    # Additive domain is the union of fall-dimensional convex sets,
    # except for the trivial additive points: symmetry reflection and x=0 and y=0.
    print_no_maximal_faces_diag(filename, q, f, faces_diag)
    print_no_maximal_faces_hor(filename, q, f, faces_hor)
    print_no_maximal_faces_ver(filename, q, f, faces_ver)
    print_no_maximal_faces_0d(filename, q, f, faces_0d)

    print_fn_minimality_test(filename, q, f, m)

    print_trivial_additive_points(filename, q, f)

    for zz in range(q):
        for xx in range(q):
            x = xx / q
            z = zz / q
            if x != z:
                for step in range(1, maxstep):
                    print_translation_i_constraints(filename, q, x, z, step)
                    print_reflection_i_constraints(filename, q, x, z, step)

    for zz in range(q):
        z = zz / q
        print_directly_covered_constraints(filename, q, z)
        for step in range(1, maxstep):
            print_undirectly_covered_i_constraints(filename, q, z, step)

    #for zz in range(q):
    #    z = zz / q
    #    print >> filename, '%s = 0' % covered_i_variable(q, z, maxstep - 1)

    # set maxstep = 1, search for 5 slope function without any
    # translation/reflection except for the symmetry reflection.

    z = 0
    while z < f/2:
        print >> filename, '%s + %s <= 1' % (covered_i_variable(q, z, maxstep - 1), \
                                             covered_i_variable(q, f - z - 1/q, maxstep - 1))
        z += 1/q
    z = f
    while z < (1+f)/2:
        print >> filename, '%s + %s <= 1' % (covered_i_variable(q, z, maxstep - 1), \
                                             covered_i_variable(q, 1 + f - z - 1/q, maxstep - 1))
        z += 1/q

    print_slope_constraints(filename, q, nums, m)
          
    print >> filename, 'Bounds'
    print_fn_bounds(filename, q)
    print_slope_bounds(filename, q, nums)

    print >> filename, 'Binary'
    for face in faces_2d + faces_diag + faces_hor + faces_ver + faces_0d :
        print >> filename, face_variable(q, face),

    for z in range(q):
        for step in range(maxstep):
            print >> filename, 'c_%s_%s' % (z, step),
    for z in range(q):
        for x in range(q):
            if x != z:
                for step in range(1, maxstep):
                    print >> filename, 't_%s_%s_%s' % (x, z, step),
                    print >> filename, 'r_%s_%s_%s' % (x, z, step),

    for k in range(nums):
        for j in range(q):
            print >> filename, '%s' % interval_slope_variable(j, k),
        
    print >> filename
    print >> filename, 'End'
    filename.close()

def painted_faces_and_funciton_from_solution(filename, q, showplots=True):
    """
    Read the solution file, draw 2d complex and plot fn.
    EXAMPLES::

        sage: faces, fn = painted_faces_and_funciton_from_solution(
        ...                 '/media/sf_dropbox/2q_mip/2q_13_9_4.sol', 13)
    """
    faces = []
    bkpt = [x/q for x in range(q+1)]
    values = [0 for x in range(q+1)]
    with open(filename) as sol_file:
        for line in sol_file:
            i = line.find(' ')
            s = line[0:i]
            if s[0] == 'f':
                j = s.find('_')
                k = int(s[(j + 1)::])
                v = eval(line[(i + 1)::])
                values[k] = v
            elif s[0] in set(['h','v','d','l','u','p']):
                face = variable_face(q, s)
                v = eval(line[(i + 1)::])
                if v == 0:
                    faces.append(face)
    fn = piecewise_function_from_breakpoints_and_values(bkpt, values)
    if showplots:
        plot_painted_faces(q, faces).show(show_legend=False)
        plot_with_colored_slopes(fn).show()
    return faces, fn

def investigate_faces_solution(q, f, faces):
    """
    Check the vertex-function corresponding to given painted faces.
    EXAMPLES::

        sage: faces, fn = painted_faces_and_funciton_from_solution( \
        ...               '/media/sf_dropbox/2q_mip/5slope_22_10_m0_min_add_point.sol',\
        ...               22, showplots=False)
        sage: investigate_faces_solution(22, 10/22, faces)
    """
    components = generate_covered_intervals_from_faces(faces)
    additive_vertices = generate_additive_vertices_from_faces(q, faces)
    fn_sym = generate_symbolic_continuous(None, components, field=QQ)
    ff = int(f * q)
    h_list = []
    for h in generate_vertex_function(q, ff, fn_sym, additive_vertices):
        #print h
        #extremality_test(h,True)
        h_list.append(h)
    return h_list

def slope_variable(k):
    """
    EXAMPLES::

        sage: slope_variable(3)
        's_3'
    """
    return 's_%s' % k

def interval_slope_variable(j, k):
    """
    EXAMPLES::

        sage: interval_slope_variable(7, 3)
        'i_7_s_3'
    """
    return 'i_%s_s_%s' % (j, k)

def print_slope_constraints(filename, q, nums, m=0):
    """
    EXAMPLES::

        sage: print_slope_constraints(sys.stdout, 3, 3)
        s_0 - s_1 > 0.0333333333333333
        s_1 - s_2 > 0.0333333333333333
        s_0 - 3 fn_1 = 0
        i_0_s_0 = 1
        s_2 + 3 fn_2 = 0
        i_2_s_2 = 1
        s_0 + 3 fn_1 - 3 fn_2 + 6 i_1_s_0 <= 6
        s_0 + 3 fn_1 - 3 fn_2 - 6 i_1_s_0 >= -6
        s_1 + 3 fn_1 - 3 fn_2 + 6 i_1_s_1 <= 6
        s_1 + 3 fn_1 - 3 fn_2 - 6 i_1_s_1 >= -6
        s_2 + 3 fn_1 - 3 fn_2 + 6 i_1_s_2 <= 6
        s_2 + 3 fn_1 - 3 fn_2 - 6 i_1_s_2 >= -6
        + i_0_s_0 + i_0_s_1 + i_0_s_2 = 1
        + i_1_s_0 + i_1_s_1 + i_1_s_2 = 1
        + i_2_s_0 + i_2_s_1 + i_2_s_2 = 1
        + i_0_s_0 + i_1_s_0 + i_2_s_0 >= 1
        + i_0_s_1 + i_1_s_1 + i_2_s_1 >= 1
        + i_0_s_2 + i_1_s_2 + i_2_s_2 >= 1
    """
    # s_0 > s_1 > ... > s_nums-1
    for k in range(0, nums - 1):
        if m == 0:
            print >> filename, '%s - %s >= 0' % (slope_variable(k), slope_variable(k+1))
        else:
            print >> filename, '%s - %s >= %s' % (slope_variable(k), slope_variable(k+1), RR(q/m))

    # first interval has the largest positive slope s_0
    print >> filename, 's_0 - %s fn_1 = 0' % q
    print >> filename, 'i_0_s_0 = 1'
    # last interval has slope s_nums-1
    print >> filename, 's_%s + %s fn_%s = 0' % (nums - 1, q, q - 1)
    print >> filename, 'i_%s_s_%s = 1' % (q - 1, nums - 1)
    # Condition: s_k + q(fn_j - fn_(j+1)) = 0 iff i_j_s_k = 1
    # ==> 1) s_k + q * fn_j - q * fn_(j+1) <= 2*q * (1 - i_j_s_k)
    # ==> 2) s_k + q * fn_j - q * fn_(j+1) >= - 2*q * (1 - i_j_s_k)
    # ==> 3) sum i_j_s_k over k = 1
    for j in range(1, q-1):
        for k in range(nums):       
            print >> filename, 's_%s + %s fn_%s - %s fn_%s + %s %s <= %s' % (k, q, j, q, j + 1, 2*q, interval_slope_variable(j, k), 2*q)
            print >> filename, 's_%s + %s fn_%s - %s fn_%s - %s %s >= %s' % (k, q, j, q, j + 1, 2*q, interval_slope_variable(j, k), -2*q)
    for j in range(q):
        for k in range(nums):
            print >> filename, '+ %s' % interval_slope_variable(j, k),
        print >> filename, '= 1'
    # Condition: sum i_j_s_k over j >= 1
    for k in range(nums):
        for j in range(q):
            print >> filename, '+ %s' % interval_slope_variable(j, k),
        print >> filename, '>= 1'

def print_slope_bounds(filename, q, nums):
    """
    EXAMPLES::

        sage: print_slope_bounds(sys.stdout, 3, 3)
        -3 <= s_0 <= 3
        -3 <= s_1 <= 3
        -3 <= s_2 <= 3
    """
    for k in range(nums):
        print >> filename, '%s <= %s <= %s' % (-q, slope_variable(k), q)

def print_no_maximal_faces_diag(filename, q, f, faces_diag):
    """
    EXAMPLES::

        sage: q = 3; f = 2/3;
        sage: faces_2d, faces_diag, faces_hor, faces_ver, faces_0d = all_faces(q)
        sage: print_no_maximal_faces_diag(sys.stdout, q, f, faces_diag)
        l_0_0 + u_0_0 - d_0_0 <= 1
        l_0_2 + u_0_2 - d_0_2 <= 1
        l_1_1 + u_1_1 - d_1_1 <= 1
        l_1_2 + u_1_2 - d_1_2 <= 1
        l_2_0 + u_2_0 - d_2_0 <= 1
        l_2_1 + u_2_1 - d_2_1 <= 1
    """
    for face in faces_diag:
        (i, j, k) = face.minimal_triple
        if k[0] != f and k[0] != 1 + f:
            face_l = Face(( i, j, [k[0] - 1/q, k[0]] ))
            face_u = Face(( i, j, [k[0], k[0] + 1/q] ))
            var_d = face_variable(q, face)
            var_l = face_variable(q, face_l)
            var_u = face_variable(q, face_u)
            # var_d == 0 iff var_l == 0 or var_u == 0,
            # i.e. var_d + 1 >= var_l + var_u
            print >> filename, '%s + %s - %s <= 1' % (var_l, var_u, var_d)

def print_no_maximal_faces_hor(filename, q, f, faces_hor):
    """
    EXAMPLES::

        sage: q = 3; f = 2/3;
        sage: faces_2d, faces_diag, faces_hor, faces_ver, faces_0d = all_faces(q)
        sage: print_no_maximal_faces_hor(sys.stdout, q, f, faces_hor)
        l_0_1 + u_0_0 - h_0_1 <= 1
        l_0_2 + u_0_1 - h_0_2 <= 1
        l_1_1 + u_1_0 - h_1_1 <= 1
        l_1_2 + u_1_1 - h_1_2 <= 1
        l_2_1 + u_2_0 - h_2_1 <= 1
        l_2_2 + u_2_1 - h_2_2 <= 1
    """
    for face in faces_hor:
        (i, j, k) = face.minimal_triple
        var_h = face_variable(q, face)
        if j[0] != 0 and j[0] != 1:
            face_l = Face(( i, [j[0], j[0] + 1/q], k ))
            face_u = Face(( i, [j[0] - 1/q, j[0]], k ))
            var_l = face_variable(q, face_l)
            var_u = face_variable(q, face_u)
            # var_h == 0 only if var_l == 0 or var_u == 0,
            # i.e. var_h + 1 >= var_l + var_u
            print >> filename, '%s + %s - %s <= 1' % (var_l, var_u, var_h)

def print_no_maximal_faces_ver(filename, q, f, faces_ver):
    """
    EXAMPLES::

        sage: q = 3; f = 2/3;
        sage: faces_2d, faces_diag, faces_hor, faces_ver, faces_0d = all_faces(q)
        sage: print_no_maximal_faces_ver(sys.stdout, q, f, faces_ver)
        l_1_0 + u_0_0 - v_1_0 <= 1
        l_2_0 + u_1_0 - v_2_0 <= 1
        l_1_1 + u_0_1 - v_1_1 <= 1
        l_2_1 + u_1_1 - v_2_1 <= 1
        l_1_2 + u_0_2 - v_1_2 <= 1
        l_2_2 + u_1_2 - v_2_2 <= 1
    """
    for face in faces_ver:
        (i, j, k) = face.minimal_triple
        var_v = face_variable(q, face)
        if i[0] != 0 and i[0] != 1:
            face_l = Face(( [i[0], i[0] + 1/q], j, k ))
            face_u = Face(( [i[0] - 1/q, i[0]], j, k ))
            var_l = face_variable(q, face_l)
            var_u = face_variable(q, face_u)
            # var_v == 0 only if var_l == 0 or var_u == 0,
            # i.e. var_v + 1 >= var_l + var_u
            print >> filename, '%s + %s - %s <= 1' % (var_l, var_u, var_v)

def print_no_maximal_faces_0d(filename, q, f, faces_0d):
    """
    EXAMPLES::

        sage: q = 3; f = 2/3;
        sage: faces_2d, faces_diag, faces_hor, faces_ver, faces_0d = all_faces(q)
        sage: print_no_maximal_faces_0d(sys.stdout, q, f, faces_0d)
        h_0_1 + h_1_1 + v_1_0 + v_1_1 + d_0_1 + d_1_0 - p_1_1 <= 5
        h_0_2 + h_1_2 + v_1_1 + v_1_2 + d_0_2 + d_1_1 - p_1_2 <= 5
        h_1_1 + h_2_1 + v_2_0 + v_2_1 + d_1_1 + d_2_0 - p_2_1 <= 5
        h_1_2 + h_2_2 + v_2_1 + v_2_2 + d_1_2 + d_2_1 - p_2_2 <= 5
    """
    for face in faces_0d:
        (x, y) = face.vertices[0]
        if x != 0 and x != 1 and y != 0 and y != 1:
            xx = int(x * q)
            yy = int(y * q)
            p  = 'p_%s_%s' % (xx, yy)
            h1 = 'h_%s_%s' % (xx - 1, yy)
            h2 = 'h_%s_%s' % (xx, yy)
            v1 = 'v_%s_%s' % (xx, yy - 1)
            v2 = 'v_%s_%s' % (xx, yy)
            d1 = 'd_%s_%s' % (xx - 1, yy)
            d2 = 'd_%s_%s' % (xx, yy - 1)
            # p == 0 only if at least one of h1, h2, v1, v2, d1, d2 is 0
            # i.e. p + 5 >= h1 + h2 + v1 + v2 + d1 + d2
            print >> filename, '%s + %s + %s + %s + %s + %s - %s <= 5' % (h1, h2, v1, v2, d1, d2, p)
 
