def initial_faces_for_paint_complex(q, ff, aa):
    """
    Return inital 1d-faces corresponding to:
        translation between (aa-1)/q and (ff-aa)/q,
        reflection around f/2 and (1-f)/2, and
        f(0) = 0;

    EXAMPLES::

        sage: q=8; ff=7; aa=2
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
    additive_vertices = generate_additive_vertices_from_faces(q, faces )
    new_faces = generate_faces_from_vertices(q, additive_vertices)
    return new_faces

import random

######
#global picked_face_set
#picked_face_set = set([])
######

def paint_complex(q, ff, aa, faces, candidate_face_set, last_non_candidate):
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

        sage: q=8; ff=7; aa=2
        sage: faces = initial_faces_for_paint_complex(q, ff, aa)
        sage: non_candidate = set([])
        sage: candidate_face_set = generate_candidate_face_set(q, ff, aa, [], None, non_candidate)
        sage: additive_vertices, green_faces, covered_intervals = 
        ...        paint_complex(q, ff, aa, faces, candidate_face_set, non_candidate).next()
        sage: plot_painted_faces(q, green_faces)
    """
    non_candidate = copy(last_non_candidate)
    while candidate_face_set:
        face_picked = random.sample(candidate_face_set, 1)[0]
        ###### debug...
        #print face_picked
        #picked_face_set.add(face_picked)
        ###### ...
        additive_vertices = generate_additive_vertices_from_faces(q, faces + [face_picked])
        new_faces = generate_faces_from_vertices(q, additive_vertices)
        covered_intervals = generate_covered_intervals_from_faces(new_faces)
        implied_additive_vertices = generate_implied_additive_vertices(q, ff, additive_vertices, covered_intervals)
        if not implied_additive_vertices is False:
            # implied_additive_vertices is False means infeasible. continue to next face in candidate_face_set
            if implied_additive_vertices:
                # implied_additive_vertices is not empty, update the following things
                additive_vertices.update(implied_additive_vertices)
                new_faces = generate_faces_from_vertices(q, additive_vertices)
                covered_intervals = generate_covered_intervals_from_faces(new_faces)
            if not intersect_with_segments(covered_intervals, [(aa - 1)/q]) and \
                   new_faces_are_legal(new_faces, non_candidate):
                new_candidate_face_set = generate_candidate_face_set(q, ff, aa, covered_intervals, face_picked, non_candidate)
                if not new_candidate_face_set:
                    # stop recursion
                    if generate_to_cover_set(q, covered_intervals) == set([(aa - 1) / q, (ff - aa) / q]):
                        # all covered, finish
                        ######### debug...
                        #print "    Painting exists for q = %s, ff = %s, aa = %s" % (q, ff, aa)
                        #plot_painted_faces(q, new_faces).show(show_legend=False)
                        #print picked_face_set
                        ######### ...
                        yield additive_vertices, new_faces, covered_intervals
                else:
                    for additive_vertices, new_faces, covered_intervals in \
                                paint_complex(q, ff, aa, new_faces, new_candidate_face_set, non_candidate):
                        yield additive_vertices, new_faces, covered_intervals

        ##### debug...
        #    else:
        #        # infeasible. continue to next face in candidate_face_set
        #        print "disselect this face %s" % face_picked
        #else:
        #    # infeasible. continue to next face in candidate_face_set
        #    print "disselect this face %s" % face_picked
        #picked_face_set.remove(face_picked)
        #### ...
        candidate_face_set.remove(face_picked)
        non_candidate.add(face_picked)

def new_faces_are_legal(new_faces, non_candidate):
    """
    check if none of implied additive faces is in non_candidate
    """
    for face in new_faces:
        if face in non_candidate:
            return False
    return True

def generate_candidate_face_set(q, ff, aa, covered, last_face=None, non_candidate=set([])):
    """
    Return a set of candidate_faces (lexicographically > last_face, not in non_candidate)
    to paint in next step, whose I, J are currently uncovered,
    and are different from [(aa-1)/q, aa/q] and its reflection.
    Note that candidate_faces only takes faces with I <= J.

    EXAMPLES::

        sage: q=8; ff=7; aa=2
        sage: candidate_face_set = generate_candidate_face_set(q, ff, aa, [])
        sage: plot_painted_faces(q, candidate_face_set)
    """
    to_cover = generate_to_cover_set(q, covered)
    to_cover.discard((aa - 1) / q)
    to_cover.discard((ff - aa) / q)
    # NOTE: candidate_faces only takes faces with x <= y
    candidate_faces = set([])
    for x in to_cover:
        for y in to_cover:
            if x <= y:
                if (x + y) != (aa - 1)/q and (x + y) != ((ff - aa) / q):
                    face = Face(([x, x + 1/q], [y, y + 1/q], [x + y, x + y + 1/q]))
                    if (last_face is None or face > last_face):
                        candidate_faces.add(face)
                if (x + y + 1/q) != (aa - 1)/q and (x + y + 1/q) != ((ff - aa) / q):
                    face = Face(([x, x + 1/q], [y, y + 1/q], [x + y + 1/q, x + y + 2/q]))
                    if (last_face is None or face > last_face):
                        candidate_faces.add(face)
    return candidate_faces - non_candidate

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
    to_cover_set = generate_to_cover_set(q, covered_intervals)
    uncovered_components = [[[x, x + 1/q]] for x in to_cover_set]
    components = covered_intervals  + uncovered_components
    fn_sym = generate_symbolic_continuous(None, components, field=QQ)
    ieqdic, eqndic = generate_ieqs_and_eqns(q, ff, fn_sym, additive_vertices)
    p = Polyhedron(ieqs = ieqdic.keys(), eqns = eqndic.keys())
    if p.is_empty():
        # infeasible
        return False

    implied_additive_vertices = set([])

    for ieq in ieqdic.keys():
        if is_implied_eq(ieq, p):
            implied_additive_vertices.update(ieqdic[ieq])
    return implied_additive_vertices

def is_implied_eq(ieq, p):
    """
    Return if ieq <type 'tuple'> must hold with equality
    on every vertices of polyhedron p
    """
    n = p.ambient_dim()
    for x in p.vertices():
        if sum([ x[i] * ieq[i + 1] for i in range(n) ]) != - ieq[0]:
            return False
    return True

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
            ####### debug...
            else:
                print "%s gives a %s-slope function. Ignore" % (x, k)
                # this will probably print many lines
            ####### ...
    else:
        # this shouldn't happen!
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

        sage: random_2q_example(8, 7, 2)
        sage: random_2q_example(8, None, None)
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
        print "attempt #%s with q = %s, ff = %s, aa = %s" % (attempts, q, ff, aa)
        faces = initial_faces_for_paint_complex(q, ff, aa)
        covered_intervals = generate_covered_intervals_from_faces(faces)
        candidate_face_set = generate_candidate_face_set(q, ff, aa, covered_intervals)
        for additive_vertices, green_faces, covered_intervals in paint_complex(q, ff, aa, faces, candidate_face_set, set([])):
            if len(covered_intervals) >= 2:
            # otherwise, too few slopes in covered_intervals, continue to the next attempt.
                final_uncovered = [[(aa - 1) / q, aa / q], [(ff - aa) / q, (ff - aa + 1) / q]]
                components = covered_intervals  + [final_uncovered]
                fn_sym = generate_symbolic_continuous(None, components, field=QQ)
                for h in generate_vertex_function(q, ff, fn_sym, additive_vertices):
                    print h
                    #print "in attempt #%s with q = %s, ff = %s, aa = %s" % (attempts, q, ff, aa)
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

def generate_covered_intervals_from_faces(faces):
    covered_intervals = generate_directly_covered_intervals_from_faces(faces)
    # debugging plot:
    # show(plot_covered_intervals(function, covered_intervals), \
    #      legend_fancybox=True, \
    #      legend_title="Directly covered, merged", \
    #      legend_loc=2) # legend in upper left

    edges = [ face.minimal_triple for face in faces if face.is_1D()]

    any_change = True
    ## FIXME: Here we saturate the covered interval components
    ## with the edge relations.  There should be a smarter way
    ## to avoid this while loop.  Probably by keeping track 
    ## of a set of non-covered components (connected by edges).
    ## --Matthias
    while any_change:
        any_change = False
        for edge in edges:
            intervals = []
            # 0 stands for I; 1 stands for J; 2 stands for K
            IJK = []
            for i in range(len(edge)):
                if len(edge[i]) == 2:
                    intervals.append(edge[i])
                    IJK.append(i)
            if edge_merge(covered_intervals,intervals,IJK):
                any_change = True

    covered_intervals = remove_empty_comp(covered_intervals)
    return covered_intervals

def generate_directly_covered_intervals_from_faces(faces):
    covered_intervals = []
    for face in faces:
        if face.is_2D():
            component = []
            for int1 in face.minimal_triple:
                component.append(interval_mod_1(int1))
            component.sort()
            component = merge_within_comp(component)
            covered_intervals.append(component)
            
    remove_duplicate(covered_intervals)
    
    #show(plot_covered_intervals(function, covered_intervals), xmax=1.5)

    for i in range(len(covered_intervals)):
        for j in range(i+1, len(covered_intervals)):
            if find_interior_intersection(covered_intervals[i], covered_intervals[j]):
                covered_intervals[j] = merge_two_comp(covered_intervals[i],covered_intervals[j])
                covered_intervals[i] = []
                    
    covered_intervals = remove_empty_comp(covered_intervals)
    return covered_intervals