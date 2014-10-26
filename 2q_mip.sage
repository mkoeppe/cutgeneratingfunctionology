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

def print_fn_minimality_test(filename, q, f):
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
    # subadditivity and additivity conditions 
    for i in range(1, q):
        for j in range(i, q):
            x = bkpt[i]
            y = bkpt[j]
            z = fractional(x + y)
            print >> filename, '%s + %s - %s >= 0' %(fn_variable(q, x), fn_variable(q, y), fn_variable(q, z))
            print >> filename, '%s + %s - %s - 2 %s <= 0' %(fn_variable(q, x), fn_variable(q, y), \
                                                            fn_variable(q, z), vertex_variable(q, (x, y)))

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

def print_interval_covered_constraints(filename, q, z):
    """
    EXAMPLES::

        sage: print_interval_covered_constraints(sys.stdout, 3, 1/3)
        c_1 = 0
        + l_1_0 + u_1_0 + l_1_1 + u_1_1 + l_1_2 + u_1_2 + l_0_1 + l_1_0 + l_2_2 + u_0_0 + u_1_2 + u_2_1 + t_0_1 + r_0_1 + t_2_1 + r_2_1 <= 15
    """
    print >> filename, '%s = 0' % covered_interval_variable(q, z)
    m = 0
    bkpt = [x/q for x in range(q)]
    for y in bkpt: 
        # I projection: l_zz,yy and u_zz,yy
        print >> filename, '+ %s + %s' % ( face_variable(q, Face(([z, z+1/q], [y, y+1/q], [z + y, z + y + 1/q]))), \
                                           face_variable(q, Face(([z, z+1/q], [y, y+1/q], [z + y + 1/q, z + y + 2/q]))) ),
        m += 2
    for x in bkpt:
        # K projection: l_xx,zz-xx
        y = z - x
        if y < 0:
            y += 1
        print >> filename, '+ %s' % face_variable(q, Face(([x, x+1/q], [y, y+1/q], [x + y, x + y + 1/q]))),
        m += 1
    for x in bkpt:
        # K projection: u_xx,zz-xx-1
        y = z - x - 1/q
        if y < 0:
            y += 1
        print >> filename, '+ %s' % face_variable(q, Face(([x, x+1/q], [y, y+1/q], [x + y + 1/q, x + y + 2/q]))),
        m += 1
    for x in bkpt:
        if x != z:
            # t_xx,zz and r_xx,zz
            print >> filename, '+ %s + %s' % (translation_variable(q, x, z), reflection_variable(q, x, z)),
            m += 2
    m = m - 1
    assert m == 6*q - 3
    print >> filename, '<= %s' %m 

def print_interval_not_covered_constraints(filename, q, z):
    """
    EXAMPLES::

        sage: print_interval_not_covered_constraints(sys.stdout, 3, 1/3)
        c_1 = 1
        l_1_0 >= 1
        u_1_0 >= 1
        l_1_1 >= 1
        u_1_1 >= 1
        l_1_2 >= 1
        u_1_2 >= 1
        l_0_1 >= 1
        l_1_0 >= 1
        l_2_2 >= 1
        u_0_0 >= 1
        u_1_2 >= 1
        u_2_1 >= 1
        t_0_1 >= 1
        r_0_1 >= 1
        t_2_1 >= 1
        r_2_1 >= 1
    """
    print >> filename, '%s = 1' % covered_interval_variable(q, z)
    bkpt = [x/q for x in range(q)]
    for y in bkpt: 
        # I projection: l_zz,yy and u_zz,yy
        print >> filename, '%s >= 1' % face_variable(q, Face(([z, z+1/q], [y, y+1/q], [z + y, z + y + 1/q])))
        print >> filename, '%s >= 1' % face_variable(q, Face(([z, z+1/q], [y, y+1/q], [z + y + 1/q, z + y + 2/q])))
    for x in bkpt:
        # K projection: l_xx,zz-xx
        y = z - x
        if y < 0:
            y += 1
        print >> filename, '%s >= 1' % face_variable(q, Face(([x, x+1/q], [y, y+1/q], [x + y, x + y + 1/q])))
    for x in bkpt:
        # K projection: u_xx,zz-xx-1
        y = z - x - 1/q
        if y < 0:
            y += 1
        print >> filename, '%s >= 1' % face_variable(q, Face(([x, x+1/q], [y, y+1/q], [x + y + 1/q, x + y + 2/q])))
    for x in bkpt:
        if x != z:
            # t_xx,zz and r_xx,zz
            print >> filename, '%s >= 1' % translation_variable(q, x, z)
            print >> filename, '%s >= 1' % reflection_variable(q, x, z)

def print_objective(filename, q):
    """
    EXAMPLES::

        sage: print_objective(sys.stdout, 2)
        + l_0_0 + u_0_0 + l_0_1 + u_0_1 + l_1_0 + u_1_0 + l_1_1 + u_1_1 
    """
    bkpt = [x/q for x in range(q)]
    for x in bkpt:
        for y in bkpt:
            print >> filename, '+ %s + %s' % ( face_variable(q, Face(([x, x+1/q], [y, y+1/q], [x + y, x + y + 1/q]))), \
                                               face_variable(q, Face(([x, x+1/q], [y, y+1/q], [x + y + 1/q, x + y + 2/q]))) ),
    print >> filename
    
def write_lpfile(q, f, a):
    """
    EXAMPLES:

        sage: write_lpfile(8, 7/8, 2/8)
    """
    filename = open(destdir + "2q_%s_%s_%s.lp" % (q, int(f*q), int(a*q)), "w")
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

    print >> filename, '\ MIP model for 2q_search with q = %s, f = %s, a = %s' % (q, f, a)
    print >> filename, 'Maximize'
    print_objective(filename, q)
    print >> filename, 'Subject to'
    for face in faces_2d + faces_diag + faces_hor + faces_ver:
        #if face.minimal_triple[0][0] <= face.minimal_triple[1][0]:
            print_logical_constraints(filename, q, face)
    
    for face in faces_0d:
        if face.minimal_triple[0][0] < face.minimal_triple[1][0]:
            print_xy_swapped_constraints(filename, q, face)

    print_fn_minimality_test(filename, q, f)

    for zz in range(q):
        for xx in range(q):
            x = xx / q
            z = zz / q
            if x != z:
                print_translation_constraints(filename, q, x, z)
                print_reflection_constraints(filename, q, x, z)

    for zz in range(q):
        z = zz / q
        if z != a - 1/q and z != f - a:
            print_interval_covered_constraints(filename, q, z)
        else:
            print_interval_not_covered_constraints(filename, q, z)
          
    print >> filename, 'Bounds'
    print_fn_bounds(filename, q)
    
    #print >> filename, 'Generals'
    #for xx in range(q):
    #    print >> filename, '%s' % fn_variable(q, xx/q),
    #print >> filename, '%s' % fn_variable(q, 1)
    
    print >> filename, 'Binary'
    for face in faces_2d + faces_diag + faces_hor + faces_ver + faces_0d :
        print >> filename, face_variable(q, face),
    for z in range(q):
        print >> filename, 'c_%s' % z,
    for z in range(q):
        for x in range(q):
            if x != z:
                print >> filename, 't_%s_%s' % (x, z),
                print >> filename, 'r_%s_%s' % (x, z),
    print >> filename
    print >> filename, 'End'
    filename.close()

#http://www.gurobi.com/documentation/5.6/reference-manual/lp_format
#\ LP format example
#Maximize
#  x + y + z
#Subject To
#  c0: x + y = 1
#  c1: x + 5 y + 2 z <= 10
#  qc0: x + y + [ x ^ 2 - 2 x * y + 3 y ^ 2 ] <= 5
#Bounds
#  0 <= x <= 5
#  z >= 2
#Generals
#  x y z
#Binary
#  x y z
#End
