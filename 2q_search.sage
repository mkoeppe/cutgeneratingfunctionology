## initial parameter for paint_complex:
## reflection
# faces = [ Face(([0,f], [0,f], [f])), Face(([f,1], [f,1], [1+f]))]
## translation
# faces += [ Face(([b - a + grid], [a - grid, a], [b, b + grid])), \
#            Face(([a - grid, a], [b - a + grid], [b, b + grid])) ]
## f(0) = 0; f(1) = 0;
# faces += [ Face(([0, 1], [0], [0, 1])), \
#            Face(([0], [0, 1], [0, 1])) ]
# candidate_face_set = generate_candidate_face_set(q, aa, [])

def paint_complex(q, ff, aa, faces, candidate_face_set):
    """
    Randomly paint triangles green in a 2d-complex,
    until all intervals are covered except for [(aa-1)/q, aa/q] and its reflection.
    Return green_faces and covered_intervals.

    Inputs:
        q, ff, aa are integers.
        3 <= ff < q and aa < ff/2
        q = grid_nb = lcm of denomiators.
    """
    while candidate_face_set:
        face_picked = candidate_face_set.pop()
        covered_intervals = generate_covered_intervals_from_faces(faces + [face_picked])
        additive_vertices = generate_additive_vertices_from_faces(q, faces + [face_picked])
        additive_vertices = generate_implied_additive_vertices(q, ff, additive_vertices, covered_intervals)
        if additive_vertices:
            # otherwise, infeasible. continue to next face in candidate_face_set
            new_faces = generate_maximal_faces_from_vertices(q, additive_vertices)
            new_covered = generate_covered_intervals_from_faces(new_faces)
            if not intersect_with_segments(new_covered, [(aa - 1)/q]):
                # otherwise, infeasible. continue to next face in candidate_face_set
                new_candidate_face_set = generate_candidate_face_set(q, aa, new_covered)
                if not new_candidate_face_set:
                    # all covered, finish
                    return new_faces, new_covered
                result = paint_complex(q, ff, aa, new_faces, new_candidate_face_set)
                if not result is False:
                    # found a way of painting, return
                    return result
    # all possibilities are tired, none is feasible
    return False

def generate_candidate_face_set(q, aa, covered):
    to_cover = generate_to_cover_set(q, covered)
    to_cover.discard((aa - 1) / q)
    to_cover.discard(1 - aa / q)
    # NOTE: candidate_faces only takes faces with x <= y
    candidate_faces = [ Face(([x, x + 1/q], [y, y + 1/q], [x + y, x + y + 1/q])) \
                        for x in to_cover for y in to_cover \
                        if x <= y and (x + y) != (aa - 1)/q and (x + y) != (1 - aa / q)]
    candidate_faces += [ Face(([x, x + 1/q], [y, y + 1/q], [x + y + 1/q, x + y + 2/q])) \
                         for x in to_cover for y in to_cover \
                         if x <= y and (x + y + 1/q) != (aa - 1)/q and (x + y + 1/q) != (1 - aa / q)]
    return set(candidate_faces)

def generate_to_cover_set(q, covered):
    to_cover = set([x/q for x in range(q)])
    for component in covered:
        for i in component:
            for x in range(i[0]*q, i[1]*q):
                to_cover.discard(x/q)
    return to_cover

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

def generate_covered_intervals_from_faces(faces):
    #TODO
    pass
def generate_implied_additive_vertices(q, ff, additive_vertices, covered_intervals):
    #TODO
    pass

def generate_maximal_faces_from_vertices(q, additive_vertices):
    #TODO
    pass

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

def generate_additive_vertices_from_faces(q, faces):
    """
    Return the set of points on 2d-grid that are covered by faces.
    """
    additive_vertices = set()
    for face in faces:
        if face.is_2D():
            i, j, k = face.minimal_triple
            for x in range(i[0]*q, i[1]*q + 1):
                for y in range(j[0]*q, j[1]*q + 1):
                    if x <= y and k[0] <= (x + y) / q <= k[1]:
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
    # TODO: use a dictionary that maps ieq to the 2d-complex-vertices from which it comes;
    #       key = ieq, value = set of 2d-complex-vertices
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
    # TODO: Rewrite; use a dictionary that maps ieq to the 2d-complex-vertices from which it comes;
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
    # TODO: Rewrite
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
