destdir = "/media/sf_dropbox/2q_mip/"

y_grid = 4

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
        print >> filename, '0 <= %s <= %s' % (fn_variable(q, x), y_grid)

def print_fn_minimality_test(filename, q, f):
    """
    EXAMPLES::

        sage: print_fn_minimality_test(sys.stdout, 3, 2/3)
        fn_0 = 0
        fn_0 + fn_2 = 1
        fn_1 + fn_1 = 1
        fn_2 + fn_3 = 1
        fn_1 + fn_1 - fn_2 - 1/10000 p_1_1 >= 0
        fn_1 + fn_1 - fn_2 - 2 p_1_1 <= 0
        fn_1 + fn_2 - fn_3 - 1/10000 p_1_2 >= 0
        fn_1 + fn_2 - fn_3 - 2 p_1_2 <= 0
        fn_2 + fn_2 - fn_1 - 1/10000 p_2_2 >= 0
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
        print >> filename, '%s + %s = %s' % (fn_variable(q, x), fn_variable(q, f - x), y_grid)
        x += 1/q
    x = f
    while x <= (1+f)/2:
        print >> filename, '%s + %s = %s' % (fn_variable(q, x), fn_variable(q, 1 + f - x), y_grid)
        x += 1/q 
    # strict-subadditivity and additivity conditions
    #small_m = 0
    small_m = 1
    for i in range(1, q):
        for j in range(i, q):
            x = bkpt[i]
            y = bkpt[j]
            z = fractional(x + y)
            #print >> filename, '%s + %s - %s >= 0' %(fn_variable(q, x), fn_variable(q, y), fn_variable(q, z))
            print >> filename, '%s + %s - %s - %s %s >= 0' %(fn_variable(q, x), fn_variable(q, y), fn_variable(q, z), \
                                                             small_m, vertex_variable(q, (x, y)))
            print >> filename, '%s + %s - %s - %s %s <= 0' %(fn_variable(q, x), fn_variable(q, y), \
                                                            fn_variable(q, z), 2*y_grid, vertex_variable(q, (x, y)))

def print_trivial_additive_points(filename, q, f):
    """
    EXAMPLES::

        sage: print_trivial_additive_points(sys.stdout, 3, 2/3, 1/3)
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
        p_1_0 = 0
        p_1_1 = 0
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

def translation_variable(q, x, z):
    """
    EXAMPLES::

        sage: translation_variable(7, 1/7, 3/7)
        't_1_3'
    """
    return 't_%s_%s' % (int(x * q), int(z * q))

def reflection_variable(q, x, z):
    """
    EXAMPLES::

        sage: reflection_variable(7, 1/7, 3/7)
        'r_1_3'
    """
    return 'r_%s_%s' % (int(x * q), int(z * q))

def print_translation_constraints(filename, q, x, z):
    """
    EXAMPLES::

        sage: print_translation_constraints(sys.stdout, 7, 1/7, 3/7)
        c_1 - t_1_3 <= 0
        h_1_2 - t_1_3 <= 0
        t_1_3 - c_1 - h_1_2 <= 0
        sage: print_translation_constraints(sys.stdout, 7, 3/7, 1/7)
        c_3 - t_3_1 <= 0
        h_1_2 - t_3_1 <= 0
        t_3_1 - c_3 - h_1_2 <= 0
    """
    print >> filename, '%s - %s <= 0' % (covered_interval_variable(q, x), translation_variable(q, x, z))
    if x <= z:
        move = Face(([x, x + 1/q], [z - x], [z, z + 1/q])) # h_xx,zz-xx
    else:
        move = Face(([z, z + 1/q], [x - z], [x, x + 1/q])) # h_zz,xx-zz
    print >> filename, '%s - %s <= 0' % (face_variable(q, move), translation_variable(q, x, z))
    print >> filename, '%s - %s - %s <= 0' % (translation_variable(q, x, z), \
                                              covered_interval_variable(q, x), \
                                              face_variable(q, move))

def print_reflection_constraints(filename, q, x, z):
    """
    EXAMPLES::

        sage: print_reflection_constraints(sys.stdout, 7, 1/7, 3/7)
        c_1 - r_1_3 <= 0
        d_1_3 - r_1_3 <= 0
        r_1_3 - c_1 - d_1_3 <= 0
    """
    print >> filename, '%s - %s <= 0' % (covered_interval_variable(q, x), reflection_variable(q, x, z))
    move = Face(([x, x + 1/q], [z, z + 1/q], [x + z + 1/q])) # d_xx,zz
    print >> filename, '%s - %s <= 0' % (face_variable(q, move), reflection_variable(q, x, z))
    print >> filename, '%s - %s - %s <= 0' % (reflection_variable(q, x, z), \
                                              covered_interval_variable(q, x), \
                                              face_variable(q, move))

################# TRY ######################
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

def print_obj_max_subadd_slack(filename, q, weight=1):
    """
    subadd_slack = q * (\sum_x fn(x))

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
                m += 1
    print m
#########################################

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
    
def write_lpfile(q, f, maxstep=None):
    """
    EXAMPLES:

        sage: write_lpfile(8, 7/8, 2/8)
    """
    if maxstep is None:
        maxstep = q
    #filename = open(destdir + "2q_%s_%s_%s_%s.lp" % (q, int(f*q), int(a*q), maxstep), "w")
    filename = open(destdir + "5slope_%s_%s.lp" % (q, int(f*q)), "w")
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

    print >> filename, '\ MIP model'# '\ MIP model for 2q_search with q = %s, f = %s, a = %s' % (q, f, a)

    print >> filename, 'Maximize'
    #print_obj_max_subadd_slack(filename, q)
    #print_obj_min_directly_covered_times(filename, q)
    #print_obj_min_undirectly_covered_times(filename, q)
    #print_obj_min_covered_times_max_subadd_slack(filename, q, maxstep=maxstep)
    print_obj_5slope22(filename, q, weight=1)
    print >> filename

    print >> filename, 'Subject to'
    for face in faces_2d + faces_diag + faces_hor + faces_ver:
        #if face.minimal_triple[0][0] <= face.minimal_triple[1][0]:
            print_logical_constraints(filename, q, face)
    
    for face in faces_0d:
        if face.minimal_triple[0][0] < face.minimal_triple[1][0]:
            print_xy_swapped_constraints(filename, q, face)

    print_fn_minimality_test(filename, q, f)

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

    for zz in range(q):
        z = zz / q
        #if z != a - 1/q and z != f - a:
        #    print >> filename, '%s = 0' % covered_i_variable(q, z, maxstep - 1)
        #else:
        #    print >> filename, '%s = 1' % covered_i_variable(q, z, maxstep - 1)
        print >> filename, '%s = 0' % covered_i_variable(q, z, maxstep - 1)
          
    print >> filename, 'Bounds'
    print_fn_bounds(filename, q)

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
    print >> filename

    print >> filename, 'Generals'
    for xx in range(q):
        print >> filename, '%s' % fn_variable(q, xx/q),
    print >> filename, '%s' % fn_variable(q, 1)

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

def investigate_faces_solution(q, f, a, faces):
    """
    Check if the vertex-function corresponding to given painted faces is a good 2q-example.
    EXAMPLES::

        sage: faces, fn = painted_faces_and_funciton_from_solution(
        ...                 '/media/sf_dropbox/2q_mip/2q_13_9_4.sol', 13, showplots=False)
        sage: investigate_faces_solution(13, 9/13, 4/13, faces)
    """
    covered_intervals = generate_covered_intervals_from_faces(faces)
    additive_vertices = generate_additive_vertices_from_faces(q, faces)
    final_uncovered = [[a - 1/q, a], [f - a, f - a + 1 / q]]
    components = covered_intervals  + [final_uncovered]
    fn_sym = generate_symbolic_continuous(None, components, field=QQ)
    ff = int(f * q)
    for h in generate_vertex_function(q, ff, fn_sym, additive_vertices):
        print h
        if not extremality_test(h): # h is not extreme
            h_2q = restrict_to_finite_group(h, f=f, oversampling=2, order=None)
            if extremality_test(h_2q): # but h restricted to 1/2q is extreme
                print "h is a valid example!"
            else:
                print "h restricted to 1/2q is not extreme. Not a good example"
        else:
            print "h is extreme. Not a good example"

#http://www.gurobi.com/documentation/5.6/reference-manual/lp_format
#\ LP format example

# Gurobi command
# m = read('2q_8_7_2.lp')
# m.optimize()
# m.write('2q_8_7_2.sol')

# sage: faces, fn = painted_faces_and_funciton_from_solution('/media/sf_dropbox/2q_mip/2q_8_7_2.sol', 8)

##  model.printAttr('x')          # all non-zero solution values
##  model.printAttr('lb', 'x*')   # bounds for vars whose names begin with 'x'
##  model.printAttr(['lb', 'ub']) # lower and upper bounds
# v = m.getVars()
## gurobi> print v[0].varName, v[0].x
# filename = open("2q_8_7_2.output", "w")
# for i in range(len(v)):
#     if v[i].x == 0:
#         print >> filename, '%s = %s' %(v[i].varName, v[i].x)
# filename.close()
#
#q = 35;
#for ff in range(q-1, 2, -1):
#    for aa in range(int((ff - 1)/2), 1, -1):
#        write_lpfile(q, ff/q, aa/q)
# fname = '2q_%s_%s_%s_%s.lp' % (q, ff, aa, q)
# m = read(fname)
# m.optimize()
