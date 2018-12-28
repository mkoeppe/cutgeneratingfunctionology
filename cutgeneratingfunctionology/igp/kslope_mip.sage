from six.moves import range

def fn_variable(q, x):
    return 'fn_%s' % int(x*q)

def face_variable(q, face):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
    assert s[0] in set(['h','v','d','l','u','p']), "variable %s is not a face" % s
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
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: vertex_variable(7, (1/7, 2/7))
        'p_1_2'
    """
    return 'p_%s_%s' % (int(v[0]*q), int(v[1]*q))

def variable_vertex(q, s):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: variable_vertex(7, 'p_1_2')
        (1/7, 2/7)
    """
    assert s[0] == 'p', "variable %s is not a vertex" % s
    index1 = s.find('_')
    index2 = s.rfind('_')
    x = int(s[index1+1: index2])/q
    y = int(s[index2+1::])/q
    return (x, y)

def print_logical_constraints(filename, q, face):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: import sys;
        sage: print_logical_constraints(sys.stdout, 7, Face(([1/7, 2/7], [1/7, 2/7], [2/7, 3/7])))
        l_1_1 - p_1_1 >= 0
        l_1_1 - p_1_2 >= 0
        l_1_1 - p_2_1 >= 0
        l_1_1 - p_1_1 - p_1_2 - p_2_1 <= 0     
    """
    for v in face.vertices:
        print('%s - %s >= 0' %(face_variable(q, face), vertex_variable(q, v)), file=filename)
    v = face.vertices[0]
    print('%s' % face_variable(q, face), end=' ', file=filename)
    for v in face.vertices:
        print('- %s' % vertex_variable(q, v), end=' ', file=filename)
    print('<= 0', file=filename)
   
def print_xy_swapped_constraints(filename, q, face):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: import sys;
        sage: face =Face(([1/7], [2/7], [3/7]))
        sage: print_xy_swapped_constraints(sys.stdout, 7, face)
        p_1_2 - p_2_1 = 0
    """
    print('%s - %s = 0' %(face_variable(q, face), face_variable(q, x_y_swapped_face(face))), file=filename)


def print_fn_bounds(filename, q):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: print_fn_bounds(sys.stdout, 3)
        0 <= fn_0 <= 1
        0 <= fn_1 <= 1
        0 <= fn_2 <= 1
        0 <= fn_3 <= 1
    """
    bkpt = [x/q for x in range(q+1)]
    for x in bkpt:
        print('0 <= %s <= 1' % fn_variable(q, x), file=filename)

def print_fn_minimality_test(filename, q, f, m=0):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
    print('%s = 0' % fn_variable(q, 0), file=filename)
    # fn(f) = 1
    #print >> filename, '%s = 1' % fn_variable(q, f)
    bkpt = [x/q for x in range(q+1)]
    # symmetric conditions
    x = 0
    while x <= f/2:
        print('%s + %s = 1' % (fn_variable(q, x), fn_variable(q, f - x)), file=filename)
        x += 1/q
    x = f
    while x <= (1+f)/2:
        print('%s + %s = 1' % (fn_variable(q, x), fn_variable(q, 1 + f - x)), file=filename)
        x += 1/q 
    # strict-subadditivity and additivity conditions
    for i in range(1, q):
        for j in range(i, q):
            x = bkpt[i]
            y = bkpt[j]
            z = fractional(x + y)
            if m == 0:
                print('%s + %s - %s >= 0' %(fn_variable(q, x), fn_variable(q, y), fn_variable(q, z)), file=filename)
            else:
                print('%s + %s - %s - %s %s >= 0' %(fn_variable(q, x), fn_variable(q, y), fn_variable(q, z), \
                                                                 RR(1 / m), vertex_variable(q, (x, y))), file=filename)
            print('%s + %s - %s - 2 %s <= 0' %(fn_variable(q, x), fn_variable(q, y), \
                                                            fn_variable(q, z), vertex_variable(q, (x, y))), file=filename)

def print_trivial_additive_points(filename, q, f):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
        print('%s = 0' % vertex_variable(q, (0, x)), file=filename)
    for x in bkpt[1::]:
        print('%s = 0' % vertex_variable(q, (x, 1)), file=filename)
    # diagonals corresponding to f
    for x in bkpt:
        if x < f:
            print('%s = 0' % vertex_variable(q, (x, f - x)), file=filename)
        elif x == f:
            print('%s = 0' % vertex_variable(q, (x, f - x)), file=filename)
            print('%s = 0' % vertex_variable(q, (x, f - x + 1)), file=filename)
        elif x > f:
            print('%s = 0' % vertex_variable(q, (x, f - x + 1)), file=filename)

    #b = f - a
    #print >> filename, '%s = 0' % vertex_variable(q, (b - a + 1/q, a - 1/q))
    #print >> filename, '%s = 0' % vertex_variable(q, (b - a + 1/q, a))

def move_variable(q, x, z):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: move_variable(7, 1/7, 3/7)
        'm_1_3'
    """
    return 'm_%s_%s' % (int(x * q), int(z * q))

def covered_i_variable(q, z, i):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: covered_i_variable(7, 3/7, 2)
        'c_3_2'
    """
    return 'c_%s_%s' % (int(z * q), i)

def translation_i_variable(q, x, z, i):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: translation_i_variable(7, 1/7, 3/7, 2)
        't_1_3_2'
    """
    return 't_%s_%s_%s' % (int(x * q), int(z * q), i)

def reflection_i_variable(q, x, z, i):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: reflection_i_variable(7, 1/7, 3/7, 2)
        'r_1_3_2'
    """   
    return 'r_%s_%s_%s' % (int(x * q), int(z * q), i)

def print_directly_covered_constraints(filename, q, z):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
        print('%s - %s <= 0' % (c_z_0, l), file=filename)
        u = face_variable(q, Face(([z, z+1/q], [y, y+1/q], [z + y + 1/q, z + y + 2/q])))
        print('%s - %s <= 0' % (c_z_0, u), file=filename)
        variable_list += [l, u]
    for x in bkpt:
        # K projection: l_xx,zz-xx
        y = z - x
        if y < 0:
            y += 1
        l = face_variable(q, Face(([x, x+1/q], [y, y+1/q], [x + y, x + y + 1/q])))
        print('%s - %s <= 0' % (c_z_0, l), file=filename)
        variable_list += [l]
    for x in bkpt:
        # K projection: u_xx,zz-xx-1
        y = z - x - 1/q
        if y < 0:
            y += 1
        u = face_variable(q, Face(([x, x+1/q], [y, y+1/q], [x + y + 1/q, x + y + 2/q])))
        print('%s - %s <= 0' % (c_z_0, u), file=filename)
        variable_list += [u]
    assert len(variable_list) == 4*q
    print('%s' % c_z_0, end=' ', file=filename)
    for lu in variable_list:
        print('- %s' % lu, end=' ', file=filename)
    print('>= %s' % (1 - 4*q), file=filename)

    ## want: I,J projections don't intersect with each other
    ## problem: this model has no feasible solution.
    #for lu in variable_list[0: 2*q]:
    #    print >> filename, '+ %s' % lu,
    #print >> filename, '>= %s' % (2*q -1)

def print_move_constraints(filename, q, x, z):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: print_move_constraints(sys.stdout, 7, 1/7, 3/7)
        m_1_3 - h_1_2 <= 0
        m_1_3 - v_5_3 <= 0
        h_1_2 + v_5_3 - m_1_3 <= 1
        sage: print_move_constraints(sys.stdout, 7, 3/7, 1/7)
        m_3_1 - h_3_5 <= 0
        m_3_1 - v_2_1 <= 0
        h_3_5 + v_2_1 - m_3_1 <= 1
    """
    m_x_z = move_variable(q, x, z)
    if z >= x:
        y = z - x
    else:
        y = z - x + 1
    # Maps a horizontal edge to a forward translation.
    forward_move = Face(([x, x + 1/q], [y], [x + y, x + y + 1/q])) # h_xx,(zz-xx)%q
    # Maps a vertical edge to a backward translation.
    backward_move = Face(([1 - y], [z, z + 1/q], [1 - y + z, 1 - y + z + 1/q])) # v_(xx-zz)%q,zz
    # m_x_z = min(forward_move, backward_move)
    # m_x_z == 0 only if forward_move == 0 or backward_move == 0,
    # i.e. m_x_z + 1 >= forward_move + backward_move
    print('%s - %s <= 0' % (m_x_z, face_variable(q, forward_move)), file=filename)
    print('%s - %s <= 0' % (m_x_z, face_variable(q, backward_move)), file=filename)
    print('%s + %s - %s <= 1' % (face_variable(q, forward_move), \
                                              face_variable(q, backward_move), m_x_z), file=filename)

def print_translation_i_constraints(filename, q, x, z, i):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: print_translation_i_constraints(sys.stdout, 7, 1/7, 3/7, 2)
        c_1_1 - t_1_3_2 <= 0
        m_1_3 - t_1_3_2 <= 0
        t_1_3_2 - c_1_1 - m_1_3 <= 0
    """
    c_x_last = covered_i_variable(q, x, i-1)
    m_x_z = move_variable(q, x, z)
    t_x_z_i =  translation_i_variable(q, x, z, i)
    print('%s - %s <= 0' % (c_x_last, t_x_z_i), file=filename)
    print('%s - %s <= 0' % (m_x_z, t_x_z_i), file=filename)
    print('%s - %s - %s <= 0' % (t_x_z_i, c_x_last, m_x_z), file=filename)
    
def print_reflection_i_constraints(filename, q, x, z, i):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: print_reflection_i_constraints(sys.stdout, 7, 1/7, 3/7, 2)
        c_1_1 - r_1_3_2 <= 0
        d_1_3 - r_1_3_2 <= 0
        r_1_3_2 - c_1_1 - d_1_3 <= 0
    """
    c_x_last = covered_i_variable(q, x, i-1)
    r_x_z_i =  reflection_i_variable(q, x, z, i)
    print('%s - %s <= 0' % (c_x_last, r_x_z_i), file=filename)
    move = Face(([x, x + 1/q], [z, z + 1/q], [x + z + 1/q])) # d_xx,zz
    print('%s - %s <= 0' % (face_variable(q, move), r_x_z_i), file=filename)
    print('%s - %s - %s <= 0' % (r_x_z_i, c_x_last, face_variable(q, move)), file=filename)
    
def print_undirectly_covered_i_constraints(filename, q, z, i):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
    print('%s - %s <= 0' % (c_z_now, c_z_last), file=filename)
    variable_list = [c_z_last]
    for x in bkpt:
        if x != z:
            t_x_z_i = translation_i_variable(q, x, z, i)
            r_x_z_i = reflection_i_variable(q, x, z, i)
            print('%s - %s <= 0' % (c_z_now, t_x_z_i), file=filename)
            print('%s - %s <= 0' % (c_z_now, r_x_z_i), file=filename)
            variable_list += [t_x_z_i, r_x_z_i]
    assert len(variable_list) == 2 * q - 1
    print('%s' % c_z_now, end=' ', file=filename)
    for v in variable_list:
         print('- %s' % v, end=' ', file=filename)
    print('>= %s' % (2 - 2 * q), file=filename)     

def print_obj_max_slope_slack(filename, kslopes):
    r"""
    slope_slack = s_0 - s_(kslopes-1)

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: print_obj_max_slope_slack(sys.stdout, 5)
        s_0 - s_4
    """
    print('s_0 - s_%s' % (kslopes - 1), end=' ', file=filename)

def print_obj_max_slope_slack_with_weights(filename, kslopes, weights):
    r"""
    slope_slack_with_weights = \sum_{0 <= i <= kslopes-2} {weights[i] * (s_i - s_(i+1))}

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: print_obj_max_slope_slack_with_weights(sys.stdout, 5, [1, 1, 1, 1])
        + 1 s_0 - 1 s_1 + 1 s_1 - 1 s_2 + 1 s_2 - 1 s_3 + 1 s_3 - 1 s_4
        sage: print_obj_max_slope_slack_with_weights(sys.stdout, 5, [1, 2, 2, 1])
        + 1 s_0 - 1 s_1 + 2 s_1 - 2 s_2 + 2 s_2 - 2 s_3 + 1 s_3 - 1 s_4
    """
    if not weights:
        print('s_0 - s_%s' % (kslopes - 1), end=' ', file=filename)
    else:
        for i in range(kslopes - 1):
            print('+ %s s_%s - %s s_%s' % (weights[i], i, weights[i], i+1), end=' ', file=filename)

def print_obj_max_subadd_slack(filename, q, weight=1): #is a constant!
    r"""
    subadd_slack = q * (\sum_x {fn(x)}) = q * (q - 1)

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: print_obj_max_subadd_slack(sys.stdout, 3)
        1 fn_0 + 1 fn_1 + 1 fn_2 
    """
    bkpt = [x/q for x in range(q)]
    print('%s %s' % (weight, fn_variable(q, bkpt[0])), end=' ', file=filename)
    for x in bkpt[1::]:
        print('+ %s %s' % (weight, fn_variable(q, x)), end=' ', file=filename)

def print_obj_min_undirectly_covered_times(filename, q, step=None, weight=1):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
                    print('+ %s %s' % (weight, translation_i_variable(q, x, z, step)), end=' ', file=filename)
                    print('+ %s %s' % (weight, reflection_i_variable(q, x, z, step)), end=' ', file=filename)

def print_obj_min_covered_times_max_subadd_slack(filename, q, maxstep=None):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
                print('+ %s %s' % ( weight, vertex_variable(q, (x, y)) ), end=' ', file=filename)
                #m += 1
    #print m

def print_obj_min_add_points(filename, q, weight=1):
    bkpt = [x/q for x in range(q + 1)]
    for x in bkpt:
        for y in bkpt:
            if x <= y:
                print('+ %s %s' % ( weight, vertex_variable(q, (x, y)) ), end=' ', file=filename)

def print_obj_min_directly_covered_times(filename, q, weight=1):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: print_obj_min_directly_covered_times(sys.stdout, 2)
        + 1 l_0_0 + 1 u_0_0 + 1 l_0_1 + 1 u_0_1 + 1 l_1_0 + 1 u_1_0 + 1 l_1_1 + 1 u_1_1 
    """
    bkpt = [x/q for x in range(q)]
    for x in bkpt:
        for y in bkpt:
            print('+ %s %s + %s %s' % ( weight, face_variable(q, Face(([x, x+1/q], [y, y+1/q], [x + y, x + y + 1/q]))), \
                                                     weight, face_variable(q, Face(([x, x+1/q], [y, y+1/q], [x + y + 1/q, x + y + 2/q]))) ), end=' ', file=filename)
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

def write_lpfile(q, f, kslopes, maxstep=None, m=0, type_cover=None, weights=[]):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: write_lpfile(37, 25/37, 5, m=12, type_cover='fulldim') # not tested
    """
    if maxstep is None:
        maxstep = q
    if type_cover == 'fulldim' or type_cover == 'fulldim_covers':
        maxstep = 1
    else:
       type_cover = 'step%s' % maxstep

    if weights:
        strw = '_' + str(weights)
    else:
        strw = ''

    destdir = output_dir+"2q_mip/"
    mkdir_p(destdir)
    filename = open(destdir + "%sslope_q%s_f%s_%s_m%s%s.lp" % (kslopes, q, int(f*q), type_cover, m, strw), "w")

    faces_2d, faces_diag, faces_hor, faces_ver, faces_0d = all_faces(q)

    print(r'\ MIP model with q = %s, f = %s, num of slopes = %s, small_m = %s' % (q, f, kslopes, m), file=filename)
    if type_cover == 'fulldim':
        print(r'\ without non-trivial 0d and 1d maximal faces.', file=filename)
    elif type_cover == 'fulldim_covers':
        print(r'\ without any translation/reflection except for the symmetry reflection.', file=filename)
    else:
        print(r'\ maximum number of steps in covering = %s.' % maxstep, file=filename)

    if weights:
        print(r'\obj = max_slope_slack with weights %s.' % weights, file=filename)

    print('Maximize', file=filename)
    #print >> filename, 0
    print_obj_max_slope_slack_with_weights(filename, kslopes, weights)
    #print_obj_max_slope_slack(filename, kslopes)
    #print_obj_max_subadd_slack(filename, q) # is a constant!
    #print_obj_min_directly_covered_times(filename, q)
    #print_obj_min_undirectly_covered_times(filename, q)
    #print_obj_min_covered_times_max_subadd_slack(filename, q, maxstep=maxstep)
    #print_obj_5slope22(filename, q, weight=1)
    #print_obj_min_add_points(filename, q, weight=1)
    print(file=filename)

    print('Subject to', file=filename)
    for face in faces_2d + faces_diag + faces_hor + faces_ver:
        #if face.minimal_triple[0][0] <= face.minimal_triple[1][0]:
            print_logical_constraints(filename, q, face)
    
    for face in faces_0d:
        if face.minimal_triple[0][0] < face.minimal_triple[1][0]:
            print_xy_swapped_constraints(filename, q, face)

    if type_cover == 'fulldim':
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
                if maxstep > 1:
                    print_move_constraints(filename, q, x, z)
                for step in range(1, maxstep):
                    print_translation_i_constraints(filename, q, x, z, step)
                    print_reflection_i_constraints(filename, q, x, z, step)

    for zz in range(q):
        z = zz / q
        print_directly_covered_constraints(filename, q, z)
        for step in range(1, maxstep):
            print_undirectly_covered_i_constraints(filename, q, z, step)

    # all intervals are covered (consider symmetry).
    z = 0
    while z < f/2:
        print('%s + %s <= 1' % (covered_i_variable(q, z, maxstep - 1), \
                                             covered_i_variable(q, f - z - 1/q, maxstep - 1)), file=filename)
        z += 1/q
    z = f
    while z < (1+f)/2:
        print('%s + %s <= 1' % (covered_i_variable(q, z, maxstep - 1), \
                                             covered_i_variable(q, 1 + f - z - 1/q, maxstep - 1)), file=filename)
        z += 1/q

    print_slope_constraints(filename, q, kslopes, m)
          
    print('Bounds', file=filename)
    print_fn_bounds(filename, q)
    print_slope_bounds(filename, q, kslopes)

    print('Binary', file=filename)
    for face in faces_2d + faces_diag + faces_hor + faces_ver + faces_0d :
        print(face_variable(q, face), end=' ', file=filename)

    for z in range(q):
        for step in range(maxstep):
            print('c_%s_%s' % (z, step), end=' ', file=filename)
    for z in range(q):
        for x in range(q):
            if x != z:
                if maxstep > 1:
                    print('m_%s_%s' % (x, z), end=' ', file=filename)
                for step in range(1, maxstep):
                    print('t_%s_%s_%s' % (x, z, step), end=' ', file=filename)
                    print('r_%s_%s_%s' % (x, z, step), end=' ', file=filename)

    for k in range(kslopes):
        for j in range(q):
            print('%s' % interval_slope_variable(j, k), end=' ', file=filename)
        
    print(file=filename)
    print('End', file=filename)
    filename.close()

def painted_faces_and_funciton_from_solution(filename, q):
    r"""
    Read the solution file, return colored faces and the function fn (inexact floating point).

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: faces, fn = painted_faces_and_funciton_from_solution('solution_2q_example_m4.sol', 37)    # not tested
    """
    faces = []
    bkpt = [x/q for x in range(q+1)]
    values = [0 for x in range(q+1)]

    destdir = output_dir+"2q_mip/"
    mkdir_p(destdir)
    with open(destdir + filename) as sol_file:
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
    return faces, fn

def refind_function_from_lpsolution(filename, q, f):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h_list = refind_function_from_lpsolution('solution_5slope_fulldim_1.sol', 37, 25/37) # not tested
        sage: h_list[0]==kzh_5_slope_fulldim_1() # not tested
        True                                     # not tested
    """
    
    faces, fn = painted_faces_and_funciton_from_solution(filename, q)
    components = generate_covered_intervals_from_faces(faces)
    additive_vertices = generate_additive_vertices_from_faces(q, faces)  
    fn_sym = generate_symbolic_continuous(None, components, field=QQ)
    ff = int(f * q)
    h_list = []
    for h in generate_vertex_function(q, ff, fn_sym, additive_vertices):
        h_list.append(h)
    return h_list

def slope_variable(k):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: slope_variable(3)
        's_3'
    """
    return 's_%s' % k

def interval_slope_variable(j, k):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: interval_slope_variable(7, 3)
        'i_7_s_3'
    """
    return 'i_%s_s_%s' % (j, k)

def print_slope_constraints(filename, q, kslopes, m=0):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: print_slope_constraints(sys.stdout, 3, 3, m=0)
        s_0 - s_1 >= 0
        s_1 - s_2 >= 0
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
    # s_0 > s_1 > ... > s_kslopes-1
    for k in range(0, kslopes - 1):
        if m == 0:
            print('%s - %s >= 0' % (slope_variable(k), slope_variable(k+1)), file=filename)
        else:
            print('%s - %s >= %s' % (slope_variable(k), slope_variable(k+1), RR(q/m)), file=filename)

    # first interval has the largest positive slope s_0
    print('s_0 - %s fn_1 = 0' % q, file=filename)
    print('i_0_s_0 = 1', file=filename)
    # last interval has slope s_kslopes-1
    print('s_%s + %s fn_%s = 0' % (kslopes - 1, q, q - 1), file=filename)
    print('i_%s_s_%s = 1' % (q - 1, kslopes - 1), file=filename)
    # Condition: s_k + q(fn_j - fn_(j+1)) = 0 iff i_j_s_k = 1
    # ==> 1) s_k + q * fn_j - q * fn_(j+1) <= 2*q * (1 - i_j_s_k)
    # ==> 2) s_k + q * fn_j - q * fn_(j+1) >= - 2*q * (1 - i_j_s_k)
    # ==> 3) sum i_j_s_k over k = 1
    for j in range(1, q-1):
        for k in range(kslopes):
            print('s_%s + %s fn_%s - %s fn_%s + %s %s <= %s' % (k, q, j, q, j + 1, 2*q, interval_slope_variable(j, k), 2*q), file=filename)
            print('s_%s + %s fn_%s - %s fn_%s - %s %s >= %s' % (k, q, j, q, j + 1, 2*q, interval_slope_variable(j, k), -2*q), file=filename)
    for j in range(q):
        for k in range(kslopes):
            print('+ %s' % interval_slope_variable(j, k), end=' ', file=filename)
        print('= 1', file=filename)
    # Condition: sum i_j_s_k over j >= 1
    for k in range(kslopes):
        for j in range(q):
            print('+ %s' % interval_slope_variable(j, k), end=' ', file=filename)
        print('>= 1', file=filename)

def print_slope_bounds(filename, q, kslopes):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: print_slope_bounds(sys.stdout, 3, 3)
        -3 <= s_0 <= 3
        -3 <= s_1 <= 3
        -3 <= s_2 <= 3
    """
    for k in range(kslopes):
        print('%s <= %s <= %s' % (-q, slope_variable(k), q), file=filename)

def print_no_maximal_faces_diag(filename, q, f, faces_diag):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
            print('%s + %s - %s <= 1' % (var_l, var_u, var_d), file=filename)

def print_no_maximal_faces_hor(filename, q, f, faces_hor):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
            print('%s + %s - %s <= 1' % (var_l, var_u, var_h), file=filename)

def print_no_maximal_faces_ver(filename, q, f, faces_ver):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
            print('%s + %s - %s <= 1' % (var_l, var_u, var_v), file=filename)

def print_no_maximal_faces_0d(filename, q, f, faces_0d):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
            print('%s + %s + %s + %s + %s + %s - %s <= 5' % (h1, h2, v1, v2, d1, d2, p), file=filename)


########### Copy from 2q_search.sage #############
def plot_painted_faces(q, faces):
    r"""
    (Ignore please.)
    Return a plot of the 2d complex of q-grid with green faces.
    """
    points = [x/q for x in range(q+1)]
    values = [0 for x in range(q+1)]
    function = discrete_function_from_points_and_values(points, values)

    p = Graphics()
    p.set_legend_options(loc='upper right')
    p += plot_2d_complex(function)

    kwds = { 'alpha': proj_plot_alpha, 'zorder': -10, 'color': 'grey'}
    IJK_kwds = [ kwds for i in range(3) ]
    for face in faces:
        p += face.plot()
        p += plot_projections_of_one_face(face, IJK_kwds)
    return p

def generate_additive_vertices_from_faces(q, faces):
    r"""
    Return the set of points (x<=y) on q*q 2d-grid that are covered by faces.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: fn = gmic(1/3); q=3;
        sage: faces = generate_maximal_additive_faces(fn)
        sage: generate_additive_vertices_from_faces(q, faces)
        {(0, 0), (0, 1/3), (0, 2/3), (0, 1), (1/3, 1), (2/3, 2/3), (2/3, 1), (1, 1)}
    """
    additive_vertices = set()
    for face in faces:
        if face.is_2D():
            i, j, k = face.minimal_triple
            for x in range(i[0]*q, i[1]*q + 1):
                for y in range(j[0]*q, j[1]*q + 1):
                    if x <= y and k[0] <= (x + y) / q <= k[1]:
                        additive_vertices.add((x/q, y/q))
        elif face.is_0D():
            additive_vertices.add(face.vertices[0])
        else:
            i, j, k = face.minimal_triple
            if face.is_horizontal():
                for x in range(i[0]*q, i[1]*q + 1):
                    if x/q <= j[0]:
                        additive_vertices.add((x/q, j[0]))
            elif face.is_vertical():
                for y in range(j[0]*q, j[1]*q + 1):
                    if i[0] <= y/q:
                        additive_vertices.add((i[0], y/q))
            else:
                for x in range(i[0]*q, i[1]*q + 1):
                    if 2*x/q <= k[0]:
                        additive_vertices.add((x/q, k[0] - x/q))   
    return additive_vertices
    
def generate_ieqs_and_eqns(q, ff, fn_sym, additive_vertices):
    r"""
    Return the equalities (by additivity) and inequalities (by subadditivity)
    that the slope variables must satisfy.

    Inputs:
        - q, ff are integers.
        - fn_sym is the symbolic function generated by considering covered_intervals, as follows.
          intervals in the same component of a function must have the same slope value.
          Take slope of the i-th component as i-th unit vector. 
          Then fn_sym maps x (in [0,1]) to a vector of dim = number of components.
        - additive_vertices are the green points on the 2d-grid, where additivity is attained by fn_sym.

    Output:
        - ieqdic, eqndic are dictionaries that maps ieq/eqn to the 2d-complex-vertices from which it comes, where key = ieq or eqn, value = set of 2d-complex-vertices.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: q=3; ff=1;
        sage: components=[[[0, 1/3]], [[1/3, 1]]]
        sage: fn_sym = generate_symbolic_continuous(None, components, field=QQ)
        sage: additive_vertices = {(0, 0), (0, 1/3), (0, 2/3), (0, 1), (1/3, 1), (2/3, 2/3), (2/3, 1), (1, 1)}
        sage: ieqdic, eqndic = generate_ieqs_and_eqns(q, ff, fn_sym, additive_vertices)
        sage: ieqdic
        {(0, 0, 0): set(),
        (0, 1/3, -1/3): {(1/3, 1/3), (1/3, 2/3)},
        (0, 1/3, 0): set(),
        (0, 1/3, 1/3): set(),
        (1, -1/3, -1/3): set(),
        (1, -1/3, 0): set(),
        (1, 0, 0): set()}
        sage: eqndic
        {(-1, 1/3, 0): set(),
        (0, 0, 0): {(0, 0), (0, 1/3), (0, 2/3), (0, 1)},
        (0, 1/3, 2/3): {(1/3, 1), (2/3, 2/3), (2/3, 1), (1, 1)}}
    """
    ieqdic = {}
    eqndic = {}
    for x in range(q):
        v = fn_sym(x/q)
        ieq = tuple([0]) + tuple(v)  # fn(x/q) >= 0
        if not ieq in ieqdic:
            ieqdic[ieq]=set([])
        ieq = tuple([1]) + tuple([-w for w in v]) #fn(x/q) <=1
        if not ieq in ieqdic:
            ieqdic[ieq]=set([])
    # fn(0) = 0
    eqn = tuple([0]) + tuple(fn_sym(0))
    if not eqn in eqndic:
        eqndic[eqn] = set([]) # or = [(0,0)]?
    # fn(1) = 0
    eqn = tuple([0]) + tuple(fn_sym(1))
    if not eqn in eqndic:
        eqndic[eqn] = set([])
    # fn(f) = 1
    eqn = tuple([-1]) + tuple(fn_sym(ff/q))
    if not eqn in eqndic:
        eqndic[eqn] = set([])
    # Note: If only do this for bkpts, some implied additive points on the grid
    # (whose x or y coordinate lies in between two bkpts) will be missing!
    # FIXME: Use maximal_additive_faces, don't need inside additve_vertices
    for x in range(q+1):
        for y in range(x, q+1):
            v = tuple([0]) + tuple(delta_pi(fn_sym, x/q, y/q))
            if (x/q, y/q) in additive_vertices:
                if v in eqndic:
                    eqndic[v].add((x/q, y/q))
                else:
                    eqndic[v] = set([(x/q, y/q)])
            else:
                if v in ieqdic:
                    ieqdic[v].add((x/q, y/q))
                else:
                    ieqdic[v] = set([(x/q, y/q)])
    return ieqdic, eqndic

def generate_vertex_function(q, ff, fn_sym, additive_vertices, kslopes=3):
    r"""
    Generate real valued functions which correspond to vertices 
    of the polytope defined by ``[ieqs, eqns] = generate_ieqs_and_eqns(..)``

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: q=3; ff=1;
        sage: components=[[[0, 1/3]], [[1/3, 1]]]
        sage: fn_sym = generate_symbolic_continuous(None, components, field=QQ)
        sage: additive_vertices = {(0, 0), (0, 1/3), (0, 2/3), (0, 1), (1/3, 1), (2/3, 2/3), (2/3, 1), (1, 1)}
        sage: h = next(generate_vertex_function(q, ff, fn_sym, additive_vertices, kslopes=2))
        sage: h == gmic(1/3)
        True
    """
    ieqdic, eqndic = generate_ieqs_and_eqns(q, ff, fn_sym, additive_vertices)
    p = Polyhedron(ieqs = list(ieqdic.keys()), eqns = list(eqndic.keys()))
    if not p.is_empty():
        #if p.n_vertices() > 1:
        #    print p.Vrepresentation()
        ## print "boundedness is %s" % p.is_compact()
        for x in p.vertices():
            k = len(set(x))
            if k >= kslopes:
                #print "%s gives a %s-slope function h =" % (x, k)
                v = vector(QQ,x)
                yield v * fn_sym
            ####### debug...
            #else:
            #    print "%s gives a %s-slope function. Ignore" % (x, k)
            #   # this will probably print many lines
            ####### ...
    #else:
    #    # this shouldn't happen!
    #    print "p.is_empty() is True"
