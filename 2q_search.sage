def covered_intervals_by_adding_face(face, edges, already_covered):
    """
    Return the covered_intervals by adding a new face to some existing faces.
    Inputs:
        face: the new face to add
        edges: a list of minimal_triple of the existing 1d faces
        already_covered: the covered_intervals induced by the existing faces
    """
    directly_covered_intervals = generate_directly_covered_by_faces([face], already_covered)
    if face.is_1D():
        new_edge = [face.minimal_triple]
    else:
        new_edge = []
    covered_intervals = generate_indirectly_covered_by_edges(edges + new_edge, directly_covered_intervals)
    return covered_intervals
    
#def random_face(q):
#    ii_left = ZZ.random_element(q)
#    ii_right = ZZ.random_element(ii_left, q)
#    jj_left = ZZ.random_element(q)
#    jj_right = ZZ.random_element(jj_left, q)
#    kk_left = ZZ.random_element(ii_left + jj_left, ii_right + jj_right + 1)
#    kk_right = ZZ.random_element(kk_left, ii_right + jj_right + 1)
#    return Face(([ii_left / q, ii_right / q], [jj_left / q, jj_right / q], [kk_left / q, kk_right / q]))
    
def pick_random_face(q, ff, aa, bb, to_cover):
    """
    Return a random face whose I, J, K projections are subintervals of 
    intervals in the to_cover list, unless too many (> 500) unsuccesful tries. 
    """
    times = 0
    while True:
        times += 1
        x = pick_random_integer(q, to_cover)
        y = pick_random_integer(q, to_cover)
        if x > y:
            x, y = y, x
        z = ZZ.random_element(2)
        if (x, y + 1) == (ff, ff) or (x + 1, y) == (ff, ff):
            continue
        if not intersect_with_segments([to_cover], [(x + y + z) / q]) or \
                                      (x + z, y + z) == (ff, ff):
            z = 1 - z
            if (x + y + z == aa - 1) or (x + y + z == bb) or (x + z, y + z) == (ff, ff):
                continue
            if times < 500 and not intersect_with_segments([to_cover], [(x + y + z) / q]):
            # at first, be strict and cover each interval EXACTLY once, until too many tries.
                continue
        i = [x/q, (x+1)/q]
        j = [y/q, (y+1)/q]
        if z == 0:
            k = [(x+y)/q, (x+y+1)/q]
        else:
            k = [(x+y+1)/q, (x+y+2)/q]
        return Face((i,j,k))

def pick_random_integer(q, to_cover):
    """
    Return a random integer z,  such that
    [z/q, (z+1)/q] is a subinterval of one interval in to_cover list.
    """
    n = ZZ.random_element(len(to_cover))
    x = int(to_cover[n][0] * q)
    y = int(to_cover[n][1] * q)
    z = ZZ.random_element(x,y)
    return z

def generate_painted_complex(q, ff, aa):
    """
    Randomly paint triangles green (while trying to avoid over covering intervals) in a 2d-complex,
    until all intervals are covered except for [(aa-1)/q, aa/q] and its reflection.
    Return green_faces and covered_intervals.

    Inputs:
        q, ff, aa are integers.
        3 <= ff < q and aa < ff/2
        q = grid_nb = lcm of denomiators.
    """
    f = ff / q
    a = aa / q
    grid = 1 / q
    bb = ff - aa
    b = bb / q
    # suppose that aa < bb
    final_uncovered = [[a - grid, a], [b, b + grid]]
    # reflection
    faces_1d = [ Face(([0,f], [0,f], [f])), Face(([f,1], [f,1], [1+f]))]
    # translation 
    faces_1d += [ Face(([b - a + grid], [a - grid, a], [b, b + grid])), \
                  Face(([a - grid, a], [b - a + grid], [b, b + grid])) ]
    # f(0) = 0; f(1) = 0;
    faces_1d += [ Face(([0, 1], [0], [0, 1])), \
                  Face(([0], [0, 1], [0, 1])) ] #, \
    #              Face(([0, 1], [1], [1, 2])), \
    #              Face(([1], [0, 1], [1, 2]))  ]
    edges = [ face.minimal_triple for face in faces_1d ]
    faces = [] #2d-faces
    # first paint lower left and upper right corner?
    already_covered = []
    to_cover = interval_minus_union_of_intervals([0,1], final_uncovered)
    # logging.info("Need to cover %s" % to_cover)
    while to_cover:
        face = pick_random_face(q, ff, aa, bb, to_cover)
        #logging.info("Try %s" % face)
        new_covered = covered_intervals_by_adding_face(face, edges, already_covered)
        if intersect_with_segments(new_covered, [a - grid]):
        # FIXME: causes infinite loop!
            continue
        additive_vertices = generate_additive_vertices_from_faces(q, faces + faces_1d + [face, x_y_swapped_face(face)])
        new_additive_vertices = generate_induced_additive_vertices(q, ff, additive_vertices, new_covered)
        if not new_additive_vertices:
        # FIXME: causes infinite loop!
            continue
        faces = generate_faces_from_vertices(q, new_additive_vertices)
        print "pick %s" % face
        plot_painted_faces(q, faces+faces_1d).show(show_legend=False)
        new_covered = generate_directly_covered_by_faces(faces, [])
        already_covered = generate_indirectly_covered_by_edges(edges, new_covered)
        to_cover = uncovered_intervals_from_covered_intervals(already_covered + [final_uncovered])
        print "to_cover = %s" % to_cover
    green_faces = faces + faces_1d
    covered_intervals = [merge_within_comp(component, one_point_overlap_suffices=True) for component in already_covered] 
    return green_faces, covered_intervals

def intersect_with_segments(covered_interval, segments_left_endpoints):
    """
    Return whether the covered_interval (inner) intersects with any given length-one segments.
    """
    for intervals in covered_interval:
        for i in intervals:
            for j in segments_left_endpoints:
                if i[0] <= j < i[1]:
                    return True
    return False

def generate_induced_additive_vertices(q, ff, additive_vertices, covered_intervals):
    uncovered_intervals = uncovered_intervals_from_covered_intervals(covered_intervals)
    uncovered_components = []
    for i in uncovered_intervals:
        x = i[0]
        while x < i[1]:
            uncovered_components.append([[x, x + 1/q]])
            x += 1/q
    components = covered_intervals  + uncovered_components
    fn_sym = generate_symbolic_continuous(None, components, field=QQ)
    #h_list = list(generate_vertex_function(q, ff, fn_sym, additive_vertices))
    found = False
    common_additive_vertices = set([(x/q, y/q) for x in range(q+1) for y in range(q+1) if x <= y])
    for h in generate_vertex_function(q, ff, fn_sym, additive_vertices, in_paint_phase=True):
        h_1q = restrict_to_finite_group(h, f=ff/q, order=q)
        h_additive_vertices = {(x,y) for (x, y, z, xeps, yeps, zeps) in generate_additive_vertices(h_1q)}
        common_additive_vertices.intersection_update(h_additive_vertices)
        found = True
    if not found:
        return False
    else:
        return common_additive_vertices

def plot_painted_faces(q, faces):
    """
    Return a plot of the 2d complex of q-grid with green faces.
    """
    points = [x/q for x in range(q+1)]
    values = [0 for x in range(q+1)]
    function = discrete_function_from_points_and_values(points, values)

    p = Graphics()
    #p.set_legend_options(handlelength = 0)
    p.set_legend_options(loc='upper right')
    p += plot_2d_complex(function)

    kwds = { 'alpha': proj_plot_alpha, 'zorder': -10, 'color': 'grey'}
    IJK_kwds = [ kwds for i in range(3) ]
    for face in faces:
        p += face.plot()
        p += plot_projections_of_one_face(face, IJK_kwds)
    return p

def generate_faces_from_vertices(q, vertices):
    xy_swapped_vertices = set([(y, x) for (x, y) in vertices])
    vertices.update(xy_swapped_vertices)
    faces = []
    for x in range(q):
        for y in range(q):
            if (x/q, (y+1)/q) in vertices and ((x+1)/q, y/q) in vertices:
                if (x/q, y/q) in vertices:
                    faces.append(Face(([x/q, (x+1)/q], [y/q, (y+1)/q], [(x+y)/q, (x+y+1)/q])))
                if ((x+1)/q, (y+1)/q) in vertices:
                    faces.append(Face(([x/q, (x+1)/q], [y/q, (y+1)/q], [(x+y+1)/q, (x+y+2)/q])))
    return faces

def generate_additive_vertices_from_faces(q, faces):
    """
    Return the set of points on 2d-grid that are covered by faces.
    Assume that all given 2d-faces are unit triangles.
    """
    additive_vertices = set()
    for face in faces:
        if face.is_2D():
            # suppose every 2d-face is a unit triangle
            for (x, y) in face.vertices:
                if x <= y:
                    additive_vertices.add((x, y))
        elif face.is_0D():
            additive_vertices.add(face.vertices)
        else:
            i, j, k = face.minimal_triple
            if face.is_horizontal():
                for x in range(i[0]*q, i[1]*q + 1):
                    if x/q <= j[0]:
                        additive_vertices.add((x/q, j[0]))
            elif face.is_vertical():
                for y in range(j[0]*q, j[1]*q + 1):
                    if i[0] <= y/q:
                        additive_vertices.add((j[0], y/q))
            else:
                for x in range(i[0]*q, i[1]*q + 1):
                    if 2*x/q <= k[0]:
                        additive_vertices.add((x/q, k[0] - x/q))   
    return additive_vertices
    
def generate_ieqs_and_eqns(q, ff, fn_sym, additive_vertices):
    """
    Return the equalities (by additivity) and inequalities (by subadditivity)
    that the slope variables must satisfy.

    Inputs:
        q, ff are integers.

        fn_sym is the symbolic function generated by considering covered_intervals:
        intervals in the same component of a function must have the same slope value.
        Take slope of the i-th component as i-th unit vector. 
        Then fn_sym maps x (\in [0,1]) to a vector of dim = number of components.

        additive_vertices are the green points on the 2d-grid, 
        where additivity is attained by fn_sym.
    """
    ieqs = []
    eqns = []
    for x in range(q):
        v = fn_sym(x/q)
        ieq_0 = [0] + list(v)  # fn(x/q) >= 0
        ieq_1 = [1] + [-w for w in v] #fn(x/q) <=1
        ieqs += [ieq_0, ieq_1]
    eqns += [[0] + list(fn_sym(0)), [0] + list(fn_sym(1)), [-1] + list(fn_sym(ff/q))] #fn(0) = 0; fn(1) = 0; fn(f) = 1;
    bkpt = fn_sym.end_points()
    for x in bkpt:
        for y in bkpt:
            if x <= y:
                v = [0] + list(delta_pi(fn_sym,x,y))
                if (x, y) in additive_vertices:
                    eqns.append(v)
                else:
                    ieqs.append(v)
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]
    for x in bkpt:
        for z in bkpt2:
            if x < z < 1+x:
                y = z - x
                v = [0] + list(delta_pi(fn_sym,x,y))
                if (x, y) in additive_vertices or (y, x) in additive_vertices:
                    eqns.append(v)
                else:
                    ieqs.append(v)
    return ieqs, eqns

def generate_vertex_function(q, ff, fn_sym, additive_vertices, in_paint_phase=False):
    """
    Generate real valued functions which correspond to vertices 
    of the polytope defined by [ieqs, eqns] = generate_ieqs_and_eqns(..)
    """
    ieqs, eqns = generate_ieqs_and_eqns(q, ff, fn_sym, additive_vertices)
    p = Polyhedron(ieqs = ieqs, eqns = eqns)
    if not p.is_empty():
        #if p.n_vertices() > 1:
        #    print p.Vrepresentation()
        ## print "boundedness is %s" % p.is_compact()
        for x in p.vertices():
            if in_paint_phase:
                yield vector(QQ,x) * fn_sym
            else:
                k = len(set(x))
                if k > 2: # can even write k > 3
                    print "%s gives a %s-slope function h =" % (x, k) 
                    v = vector(QQ,x)
                    yield v * fn_sym
                #else:
                #    print "%s gives a %s-slope function. Ignore" % (x, k)

def random_2q_example(q, ff=None, aa=None):
    """
    q, ff, aa are integers.
    3 <= ff < q and aa < ff/2
    want to find a function h with 2 uncovered segments:
    [(aa-1)/q, aa/q], [bb/q, (bb+1)/q] which are related by translation and reflection,
    such that:
        h is not extreme (i.e h restricted to 1/3q is not extreme), but
        h restricted to 1/2q is extreme

    EXAMPLES::

        sage: logging.disable(logging.INFO)
        sage: random_2q_example(13, 10, 4)
        sage: logging.disable()
    """
    attempts = 0
    random_ff = (ff is None)
    random_aa = (ff is None) or (aa is None)
    while attempts < 1:
    #while attempts >=0: 
        if random_ff:
            ff = ZZ.random_element(int(q/2),q)
        if random_aa:
            aa = ZZ.random_element(int(ff/4), int((ff+1)/2))
        # print "attempt #%s with q = %s, ff = %s, aa = %s" % (attempts, q, ff, aa)
        green_faces, covered_intervals = generate_painted_complex(q, ff, aa)
        if len(covered_intervals) <= 2: # <=1? too few slopes
            continue
        additive_vertices = generate_additive_vertices_from_faces(q, green_faces)
        f = ff / q
        final_uncovered = [[(aa - 1) / q, aa / q], [(ff - aa) / q, (ff - aa + 1) / q]]
        components = covered_intervals  + [final_uncovered]
        fn_sym = generate_symbolic_continuous(None, components, field=QQ)
        for h in generate_vertex_function(q, ff, fn_sym, additive_vertices):
            print h
            print "in attempt #%s with q = %s, ff = %s, aa = %s" % (attempts, q, ff, aa)
            # should check if final_nucovered is indeed uncovered. 
            # often it is not the case. then h is not a good example.
            if not extremality_test(h): # h is not extreme
                h_2q = restrict_to_finite_group(h, f=f, oversampling=2, order=None)
                if extremality_test(h_2q): # but h restricted to 1/2q is extreme
                    return h
                else:
                    print "h restricted to 1/2q is not extreme. Not a good example"
            else:
                print "h is extreme. Not a good example"
        attempts += 1
    print "Example function not found. Please try again."
