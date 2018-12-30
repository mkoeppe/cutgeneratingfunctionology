##########################################
# MIP approach to search for 2q example
# use code in kslope_mip.sage
##########################################

from six.moves import range

def print_trivial_additive_points_2q(filename, q, f, a):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: print_trivial_additive_points_2q(sys.stdout, 3, 2/3, 1/3)
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

    b = f - a
    print('%s = 0' % vertex_variable(q, (b - a + 1/q, a - 1/q)), file=filename)
    print('%s = 0' % vertex_variable(q, (b - a + 1/q, a)), file=filename)

def write_lpfile_2q(q, f, a, kslopes, maxstep=None, m=0):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: write_lpfile_2q(37, 25/37, 11/37, 4, maxstep=2, m=4) # not tested
    """
    if maxstep is None:
        maxstep = q

    destdir = output_dir+"2q_mip/"
    mkdir_p(destdir)
    filename = open(destdir + "mip_q%s_f%s_a%s_%sslope_%smaxstep_m%s.lp" % (q, int(f*q), int(a*q), kslopes, maxstep, m), "w")
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

    print(r'\ MIP model with q = %s, f = %s, a = %s, num of slopes = %s, maxstep of tran/refl = %s, small_m = %s' % (q, f, a, kslopes, maxstep, m), file=filename)

    print('Maximize', file=filename)
    #print >> filename, 0
    print_obj_max_slope_slack(filename, kslopes)
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

    # if no 1/m trick for subadditivity, set m=0 in print_fn_minimality_test.
    print_fn_minimality_test(filename, q, f, m)

    print_trivial_additive_points_2q(filename, q, f, a)

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

    for zz in range(q):
        z = zz / q
        if z != a - 1/q and z != f - a:
            print('%s = 0' % covered_i_variable(q, z, maxstep - 1), file=filename)
        else:
            print('%s = 1' % covered_i_variable(q, z, maxstep - 1), file=filename)

    print_slope_constraints_2q(filename, q, f, a, kslopes, m)
          
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

def print_slope_constraints_2q(filename, q, f, a, kslopes, m=0):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: print_slope_constraints_2q(sys.stdout, 4, 3/4, 2/4, 3, m=0)
        s_0 - s_1 >= 0
        s_1 - s_2 >= 0
        s_0 - 4 fn_1 = 0
        i_0_s_0 = 1
        s_2 + 4 fn_3 = 0
        i_3_s_2 = 1
        s_0 + 4 fn_1 - 4 fn_2 + 8 i_1_s_0 <= 8
        s_0 + 4 fn_1 - 4 fn_2 - 8 i_1_s_0 >= -8
        s_1 + 4 fn_1 - 4 fn_2 + 8 i_1_s_1 <= 8
        s_1 + 4 fn_1 - 4 fn_2 - 8 i_1_s_1 >= -8
        s_2 + 4 fn_1 - 4 fn_2 + 8 i_1_s_2 <= 8
        s_2 + 4 fn_1 - 4 fn_2 - 8 i_1_s_2 >= -8
        s_0 + 4 fn_2 - 4 fn_3 + 8 i_2_s_0 <= 8
        s_0 + 4 fn_2 - 4 fn_3 - 8 i_2_s_0 >= -8
        s_1 + 4 fn_2 - 4 fn_3 + 8 i_2_s_1 <= 8
        s_1 + 4 fn_2 - 4 fn_3 - 8 i_2_s_1 >= -8
        s_2 + 4 fn_2 - 4 fn_3 + 8 i_2_s_2 <= 8
        s_2 + 4 fn_2 - 4 fn_3 - 8 i_2_s_2 >= -8
        + i_0_s_0 + i_0_s_1 + i_0_s_2 = 1
        + i_1_s_0 + i_1_s_1 + i_1_s_2 = 1
        + i_2_s_0 + i_2_s_1 + i_2_s_2 = 1
        + i_3_s_0 + i_3_s_1 + i_3_s_2 = 1
        + i_0_s_0 + i_1_s_0 + i_2_s_0 + i_3_s_0 >= 1
        + i_0_s_1 + i_1_s_1 + i_2_s_1 + i_3_s_1 >= 1
        + i_0_s_2 + i_1_s_2 + i_2_s_2 + i_3_s_2 >= 1
        i_1_s_0 + i_0_s_0 <= 1
        i_1_s_0 - i_1_s_0 = 0
        i_1_s_0 + i_2_s_0 <= 1
        i_1_s_0 + i_3_s_0 <= 1
        i_1_s_1 + i_0_s_1 <= 1
        i_1_s_1 - i_1_s_1 = 0
        i_1_s_1 + i_2_s_1 <= 1
        i_1_s_1 + i_3_s_1 <= 1
        i_1_s_2 + i_0_s_2 <= 1
        i_1_s_2 - i_1_s_2 = 0
        i_1_s_2 + i_2_s_2 <= 1
        i_1_s_2 + i_3_s_2 <= 1
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
    # Two special intervals have the same slope value,
    # which is different from the slope values of other intervals
    ja = int(a * q) - 1
    jb = int((f - a) * q)
    for k in range(kslopes):
        for j in range(q):
            if j == jb:
                print('%s - %s = 0' % ( interval_slope_variable(ja, k), \
                                                     interval_slope_variable(jb, k) ), file=filename)
            elif j != ja:
                print('%s + %s <= 1' % ( interval_slope_variable(ja, k), \
                                                      interval_slope_variable(j, k) ), file=filename)

def refind_function_from_lpsolution_2q(filename, q, f, a):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = refind_function_from_lpsolution_2q('solution_2q_example_m4.sol', 37, 25/37, 11/37) # not tested
    """
    faces, fn = painted_faces_and_funciton_from_solution(filename, q)
    covered_intervals = generate_covered_intervals_from_faces(faces)
    additive_vertices = generate_additive_vertices_from_faces(q, faces)
    final_uncovered = [[a - 1/q, a], [f - a, f - a + 1 / q]]
    components = covered_intervals  + [final_uncovered]
    fn_sym = generate_symbolic_continuous(None, components, field=QQ)
    ff = int(f * q)
    for h in generate_vertex_function(q, ff, fn_sym, additive_vertices):
        if not extremality_test(h): # h is not extreme
            logging.info("Testing the extremality of pi restricted to 1/2q.")
            h_2q = restrict_to_finite_group(h, f=f, oversampling=2, order=None)
            if extremality_test(h_2q): # but h restricted to 1/2q is extreme
                logging.info("Find a valid 2q_example!")
                return h
# Comments:
# def kzh_2q_example_1() is obtained by
#    LP model: mip_q37_f25_a11_4slope_2maxstep_m4.lp
#    q = 37; f = 25/37; a = 11/37;
#    given number of slopes = 4; maxstep = 2;
#    slope gap >= q/4; subadditivity slack = 1/4;
#    obj = max_slope_slack.
#
#    Gurobi run on my laptop (276s, obj = 33.66393)
#    Optimal solution: solution_2q_example_m4.sol
#
#    ( Or
#    (1) LP model: mip_q37_f25_a11_4slope_2maxstep_m12.lp
#    slope gap >= q/12; subadditivity slack = 1/12;
#    Gurobi run on my laptop (1080s, obj = 33.66393)
#    Optimal solution: solution_2q_example_m12.sol
#
#    (2) LP model: mip_q37_f25_a11_4slope_2maxstep_m12slopegaponly.lp
#    slope gap >= q/12; no 1/m trick on subadditivity slack;
#    Gurobi run on point.math.ucdavis.edu (7412s, obj = 33.66393)
#    Optimal solution: solution_2q_example_m12slopegaponly.sol )
#
#    The vertex function is a 4-slope non-extreme function,
#    whose restriction to 1/2q is extreme.
#
#    This two_q_example is obtained by:
#    sage: h = refind_function_from_lpsolution_2q('solution_2q_example_m4.sol', 37, 25/37, 11/37) # not tested
