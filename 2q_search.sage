def initial_faces_for_paint_complex(q, ff, aa):
    """
    Return inital 1d-faces corresponding to:
        translation between (aa-1)/q and (ff-aa)/q,
        reflection around f/2 and (1-f)/2, and
        f(0) = 0;

    EXAMPLES::

        sage: q=7; ff=10; aa=4
        sage: faces = initial_faces_for_paint_complex(q, ff, aa)
        sage: plot_painted_faces(q,faces)
    """
    f = ff / q
    a = aa / q
    b = f - a
    grid = 1 / q
    # reflection
    faces = [ Face(([0,f], [0,f], [f])), Face(([f,1], [f,1], [1+f]))]
    # translation
    faces += [ Face(([b - a + grid], [a - grid, a], [b, b + grid])), \
               Face(([a - grid, a], [b - a + grid], [b, b + grid])) ]
    # f(0) = 0, f(1) = 0
    faces += [ Face(([0, 1], [0], [0, 1])), \
               Face(([0], [0, 1], [0, 1])), \
               Face(([0, 1], [1], [1, 2])), \
               Face(([1], [0, 1], [1, 2])),  ]
    return faces

import random

def paint_complex(q, ff, aa, faces, candidate_face_set):
    """
    Randomly paint triangles green in a 2d-complex,
    until all intervals are covered except for [(aa-1)/q, aa/q] and its reflection.
    Return green_faces and covered_intervals if a possible way of painting is found,
    otherwise, return False

    Inputs:
        q, ff, aa are integers.
        3 <= ff < q and aa < ff/2
        q = grid_nb = lcm of denomiators.

    EXAMPLES::

        sage: q=7; ff=10; aa=4
        sage: faces = initial_faces_for_paint_complex(q, ff, aa)
        sage: candidate_face_set = generate_candidate_face_set(q, ff, aa, [])
        sage: paint_complex(q, ff, aa, faces, candidate_face_set)
        False
    """
    # FIXME: Perhaps it's better to use generator, "yield" all possible paintings
    while candidate_face_set:
        face_picked = random.sample(candidate_face_set, 1)[0]
        candidate_face_set.remove(face_picked)
        covered_intervals = generate_covered_intervals_from_faces(faces + [face_picked])
        additive_vertices = generate_additive_vertices_from_faces(q, faces + [face_picked])
        implied_additive_vertices = generate_implied_additive_vertices(q, ff, additive_vertices, covered_intervals)
        if not implied_additive_vertices is False:
            # implied_additive_vertices is False means infeasible. continue to next face in candidate_face_set
            ################
            if implied_additive_vertices:
                print implied_additive_vertices
            ################
            additive_vertices.update(implied_additive_vertices)
            new_faces = generate_faces_from_vertices(q, additive_vertices)
            new_covered = generate_covered_intervals_from_faces(new_faces)
            if not intersect_with_segments(new_covered, [(aa - 1)/q]):
                # otherwise, infeasible. continue to next face in candidate_face_set
                new_candidate_face_set = generate_candidate_face_set(q, ff, aa, new_covered)
                if not new_candidate_face_set:
                    # all covered, finish
                    #########
                    print "Painting exists for q = %s, ff = %s, aa = %s" % (q, ff, aa)
                    plot_painted_faces(q, new_faces).show(show_legend=False)
                    #########
                    return new_faces, new_covered
                result = paint_complex(q, ff, aa, new_faces, new_candidate_face_set)
                if not result is False:
                    # found a way of painting, return
                    return result
    # all possibilities are tired, none is feasible
    return False

def generate_candidate_face_set(q, ff, aa, covered):
    """
    Return a set of candidate_faces to paint in next step,
    whose I, J are currently uncovered,
    and are different from [(aa-1)/q, aa/q] and its reflection.

    EXAMPLES::

        sage: q=7; ff=10; aa=4
        sage: candidate_face_set = generate_candidate_face_set(q, ff, aa, [])
        sage: plot_painted_faces(q, candidate_face_set)
    """
    to_cover = generate_to_cover_set(q, covered)
    to_cover.discard((aa - 1) / q)
    to_cover.discard((ff - aa) / q)
    # NOTE: candidate_faces only takes faces with x <= y
    candidate_faces = [ Face(([x, x + 1/q], [y, y + 1/q], [x + y, x + y + 1/q])) \
                        for x in to_cover for y in to_cover \
                        if x <= y and (x + y) != (aa - 1)/q and (x + y) != ((ff - aa) / q)]
    candidate_faces += [ Face(([x, x + 1/q], [y, y + 1/q], [x + y + 1/q, x + y + 2/q])) \
                         for x in to_cover for y in to_cover \
                         if x <= y and (x + y + 1/q) != (aa - 1)/q and (x + y + 1/q) != ((ff - aa) / q)]
    return set(candidate_faces)

def generate_to_cover_set(q, covered):
    """
    Return a set {k/q | 0 <=k < q, [k/q, (k+1)/q] is uncovered} 
    """
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

def generate_implied_additive_vertices(q, ff, additive_vertices, covered_intervals):
    """
    Return new additive_vertices implied by the known additive_vertices
    """
    # FIXME: This function returns an empty set.
    # Hrepresentation is always minimal. Unable to track back to 
    # corrospending implied additive points on 2d-complex when an ieq is found to be eq.
    # We may need to use cddlib more directly
    to_cover_set = generate_to_cover_set(q, covered_intervals)
    uncovered_components = [[[x, x + 1/q]] for x in to_cover_set]
    components = covered_intervals  + uncovered_components
    fn_sym = generate_symbolic_continuous(None, components, field=QQ)
    ieqdic, eqndic = generate_ieqs_and_eqns(q, ff, fn_sym, additive_vertices)
    p = Polyhedron(ieqs = ieqdic.keys(), eqns = eqndic.keys())
    if p.is_empty():
        # infeasible
        return False

    n = p.n_vertices()
    implied_additive_vertices = set([])
    #for face_vertex in p.facial_incidences():
    #FIXME: p.facial_incidences() causes DeprecationWarning
    #    if len(face_vertex[1]) == n:
    #        i = face_vertex[0]
    #        fi = p.Hrepresentation(i)
    for fi in p.Hrepresentation():
        n_incident = len([v.index() for v in fi.incident()])
        if n_incident == n: # fi is an implied equality
            if tuple(fi) in ieqdic:
                implied_additive_vertices.update(ieqdic[fi])
    return implied_additive_vertices

def generate_faces_from_vertices(q, additive_vertices):
    # FIXME: Factor from generate_maximal_additive_faces_continuous(function)
    xy_swapped_vertices = set([(y, x) for (x, y) in additive_vertices])
    vertices = additive_vertices.union(xy_swapped_vertices)
    faces = []
    vertical_edge = set([])
    horizontal_edge = set([])
    pts = set([])
    for x in range(q):
        for y in range(x,q):
            pt_00 = (x/q, y/q)
            pt_10 = ((x+1)/q, y/q)
            pt_01 = (x/q, (y+1)/q)
            pt_11 = ((x+1)/q, (y+1)/q)
            i = [x/q, (x+1)/q]
            j = [y/q, (y+1)/q]
            k1 = [(x+y)/q, (x+y+1)/q]
            k2 = [(x+y+1)/q, (x+y+2)/q]
            if (pt_01 in vertices) and (pt_10 in vertices):
                diagonal_to_add = True
                if pt_00 in vertices:
                    faces.append(Face((i,j,k1)))
                    diagonal_to_add = False
                    vertical_edge.add(pt_00)
                    horizontal_edge.add(pt_00)
                    pts.add(pt_00)
                if pt_11 in vertices:
                    faces.append(Face((i,j,k2)))
                    diagonal_to_add = False
                    vertical_edge.add(pt_10)
                    horizontal_edge.add(pt_01)
                    pts.add(pt_11)
                if diagonal_to_add:
                    faces.append(Face((i,j,[(x+y+1)/q])))
                pts.add(pt_01)
                pts.add(pt_10)
            if (pt_00 in vertices) and (pt_01 in vertices) and not (pt_00 in vertical_edge):
                faces.append(Face(([x/q], j, k1)))
                vertical_edge.add(pt_00)
                pts.add(pt_00)
                pts.add(pt_01)
            if (pt_00 in vertices) and (pt_10 in vertices) and not (pt_00 in horizontal_edge):
                faces.append(Face((i, [y/q], k1)))
                horizontal_edge.add(pt_00)
                pts.add(pt_00)
                pts.add(pt_10)
            if (pt_00 in vertices) and not (pt_00 in pts):
                faces.append(Face(([x/q], [y/q], [(x+y)/q])))
                pts.add(pt_00)
        if (pt_01 in vertices) and (pt_11 in vertices) and not (pt_01 in horizontal_edge):
                faces.append(Face((i, [1], k2)))
                horizontal_edge.add(pt_01)
                pts.add(pt_01)
                pts.add(pt_11)
    faces += [ x_y_swapped_face(face) for face in faces \
                                       if face.minimal_triple[0] != face.minimal_triple[1] ]
    return faces

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
                        additive_vertices.add((x/q, y/q))
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

    Output:
        ieqdic, eqndic are dictionaries that maps ieq/eqn to 
        the 2d-complex-vertices from which it comes;
        key = ieq or eqn, value = set of 2d-complex-vertices.
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

    bkpt = fn_sym.end_points()
    for x in bkpt:
        for y in bkpt:
            if x <= y:
                v = tuple([0]) + tuple(delta_pi(fn_sym,x,y))
                if (x, y) in additive_vertices:
                    if v in eqndic:
                        eqndic[v].add((x,y))
                    else:
                        eqndic[v] = set([(x,y)])
                else:
                    if v in ieqdic:
                        ieqdic[v].add((x,y))
                    else:
                        ieqdic[v] = set([(x,y)])
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]
    for x in bkpt:
        for z in bkpt2:
            if x < z < 1+x:
                y = z - x
                if x > y:
                    x, y = y, x
                v = tuple([0]) + tuple(delta_pi(fn_sym,x,y))
                if (x, y) in additive_vertices:
                    if v in eqndic:
                        eqndic[v].add((x,y))
                    else:
                        eqndic[v] = set([(x,y)])
                else:
                    if v in ieqdic:
                        ieqdic[v].add((x,y))
                    else:
                        ieqdic[v] = set([(x,y)])
    return ieqdic, eqndic

def generate_vertex_function(q, ff, fn_sym, additive_vertices):
    """
    Generate real valued functions which correspond to vertices 
    of the polytope defined by [ieqs, eqns] = generate_ieqs_and_eqns(..)
    """
    ieqdic, eqndic = generate_ieqs_and_eqns(q, ff, fn_sym, additive_vertices)
    p = Polyhedron(ieqs = ieqdic.keys(), eqns = eqndic.keys())
    if not p.is_empty():
        #if p.n_vertices() > 1:
        #    print p.Vrepresentation()
        ## print "boundedness is %s" % p.is_compact()
        for x in p.vertices():
            k = len(set(x))
            if k > 2: # can even write k > 3
                print "%s gives a %s-slope function h =" % (x, k)
                v = vector(QQ,x)
                yield v * fn_sym
            else:
                print "%s gives a %s-slope function. Ignore" % (x, k)
    else:
        print "p.is_empty() is True"

def random_2q_example(q, ff=None, aa=None, max_attempts=None):
    """
    q, ff, aa are integers.
    3 <= ff < q and aa < ff/2
    want to find a function h with 2 uncovered segments:
    [(aa-1)/q, aa/q], [bb/q, (bb+1)/q] which are related by translation and reflection,
    such that:
        h is not extreme (i.e h restricted to 1/3q is not extreme), but
        h restricted to 1/2q is extreme

    EXAMPLES::

        sage: random_2q_example(13, 10, 4)
    """
    attempts = 0
    if ff is None:
        candidate_ff = range(3,q)
    else:
        candidate_ff = [ff]
    candidate_ff_aa = set([])
    if aa is None:
        for ff_0 in candidate_ff:
            for aa_0 in range(1, int((ff_0 + 1)/2)):
                candidate_ff_aa.add((ff_0, aa_0))
    else:
        for ff_0 in candidate_ff:
            if 1 <= aa < int((ff_0 + 1)/2):
                candidate_ff_aa.add((ff_0, aa))

    while candidate_ff_aa and ((max_attempts is None) or (attempts < max_attempts)): 
        (ff, aa) = random.sample(candidate_ff_aa, 1)[0]
        candidate_ff_aa.remove((ff, aa))
        # print "attempt #%s with q = %s, ff = %s, aa = %s" % (attempts, q, ff, aa)
        faces = initial_faces_for_paint_complex(q, ff, aa)
        candidate_face_set = generate_candidate_face_set(q, ff, aa, [])
        can_paint = paint_complex(q, ff, aa, faces, candidate_face_set)
        if can_paint is False:
            print "Painting doesn't exist for q = %s, ff = %s, aa = %s" % (q, ff, aa)
        elif len(can_paint[1]) >= 2:
            # otherwise, too few slopes in covered_intervals, continue to the next attempt.
            green_faces, covered_intervals = can_paint
            additive_vertices = generate_additive_vertices_from_faces(q, green_faces)
            final_uncovered = [[(aa - 1) / q, aa / q], [(ff - aa) / q, (ff - aa + 1) / q]]
            components = covered_intervals  + [final_uncovered]
            fn_sym = generate_symbolic_continuous(None, components, field=QQ)
            for h in generate_vertex_function(q, ff, fn_sym, additive_vertices):
                print h
                print "in attempt #%s with q = %s, ff = %s, aa = %s" % (attempts, q, ff, aa)
                if not extremality_test(h): # h is not extreme
                    h_2q = restrict_to_finite_group(h, f=f, oversampling=2, order=None)
                    if extremality_test(h_2q): # but h restricted to 1/2q is extreme
                        print "h is a valid example!"
                        return h
                    else:
                        print "h restricted to 1/2q is not extreme. Not a good example"
                else:
                    print "h is extreme. Not a good example"
        attempts += 1          
    print "Example function not found. Please try again."
