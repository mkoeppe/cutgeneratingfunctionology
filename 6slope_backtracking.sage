# backtracking search for 6-slope extreme functions, 
# using incremental computation for additive_vertices, additive_faces and covered_intervals.
# do not consider translation/reflection other than the symmetry reflection regarding f. (so, no edge-merge)
# copy-paste from 2q_search.sage and edited.

import random

#global num_of_slopes
num_of_slopes = 6 # set to 6 as we are looking for 6-slope functions

#global to_check_implied
to_check_implied = 0

#global check_implied_period
check_implied_period = 2

def trivial_additive_vertices_for_paint_complex(q, f):
    """
    Return additive_vertices (x<=y) corresonding to
    fn(0) = 0, fn(1) = 0 and reflection around f/2 and (1-f)/2.

    Examples::

        sage: trivial_additive_vertices_for_paint_complex(5, 3/5)
        set([(0, 4/5), (0, 1), (1/5, 2/5), (0, 0), (2/5, 1), (3/5, 1), (4/5, 4/5), (0, 2/5), (1, 1), (0, 3/5), (0, 1/5), (1/5, 1), (4/5, 1)])
    """
    additive_vertices = set([])
    bkpt = [x/q for x in range(q+1)]
    # border x = 0 and border y = 0 are green
    for x in bkpt:
        additive_vertices.add((0, x))
    for x in bkpt[1::]:
        additive_vertices.add((x, 1))
    # diagonals corresponding to f
    for x in bkpt:
        if x <= f/2:
            additive_vertices.add((x, f - x))
        elif (x >= f) and (x <= f -x +1):
            additive_vertices.add((x, f - x + 1))
    return additive_vertices

def initial_vertices_faces_and_covered_intervals(q, f):
    """
    Return additive_vertices, additive_faces_set and covered_intervals,
    corresponding to fn(0) = 0, fn(1) = 0 and reflection around f/2 and (1-f)/2.

    EXAMPLES::

        sage: additive_vertices, faces_set, covered_intervals = initial_vertices_faces_and_covered_intervals(5, 3/5);
        sage: faces_set
        set([<Face ([3/5, 4/5], [4/5, 1], [8/5, 9/5])>, <Face ([4/5, 1], [4/5, 1], [9/5, 2])>, <Face ([0, 1/5], [0, 1/5], [0, 1/5])>, <Face ([4/5, 1], [4/5, 1], [8/5, 9/5])>, <Face ([0, 1/5], [2/5, 3/5], [2/5, 3/5])>])
        sage: covered_intervals
        [[[0, 1/5], [2/5, 3/5]], [[3/5, 4/5], [4/5, 1]]]
    """
    additive_vertices = trivial_additive_vertices_for_paint_complex(q, f)
    faces_set = set([])
    covered_intervals = []
    for xx in range(q):
        for yy in range(xx, q):
            for zz in range(2):
                face = Face(([xx/q, (xx+1)/q], [yy/q, (yy+1)/q], [(xx+yy+zz)/q, (xx+yy+zz+1)/q]))
                if additive_vertices.issuperset([(x,y) for (x,y) in face.vertices if x <= y]):
                    faces_set.add(face)
                    covered_intervals = directly_covered_by_adding_face(covered_intervals, face)
    return additive_vertices, faces_set, covered_intervals

def faces_around_vertex(q, v):
    """
    Given a grid vertex v (v[0] <= v[1]), return small triangle faces (I <= J) that are around v.

    EXAMPLES::

        sage: faces_around_vertex(5, (1/5, 2/5))
        [<Face ([0, 1/5], [2/5, 3/5], [2/5, 3/5])>, <Face ([0, 1/5], [2/5, 3/5], [3/5, 4/5])>, <Face ([1/5, 2/5], [2/5, 3/5], [3/5, 4/5])>, <Face ([0, 1/5], [1/5, 2/5], [2/5, 3/5])>, <Face ([1/5, 2/5], [1/5, 2/5], [2/5, 3/5])>, <Face ([1/5, 2/5], [1/5, 2/5], [3/5, 4/5])>]
        sage: faces_around_vertex(5, (1/5, 1/5))
        [<Face ([0, 1/5], [1/5, 2/5], [1/5, 2/5])>, <Face ([0, 1/5], [1/5, 2/5], [2/5, 3/5])>, <Face ([1/5, 2/5], [1/5, 2/5], [2/5, 3/5])>, <Face ([0, 1/5], [0, 1/5], [1/5, 2/5])>]
    """
    (x, y) = v
    if x > y:
        return []
    # Note: v[0]<=v[1]; only return faces with I <= J
    xl = x - 1/q
    xr = x + 1/q
    yl = y - 1/q
    yr = y + 1/q
    faces =[ \
           Face(([xl, x],[y, yr],[xl+y, xl+y+1/q])), \
           Face(([xl, x],[y, yr],[xl+y+1/q, xl+y+2/q])), \
           Face(([x, xr],[y, yr],[x+y, x+y+1/q])), \
           Face(([xl, x],[yl, y],[xl+yl+1/q, xl+yl+2/q])), \
           ]
    if x <= yl:
        faces += [ \
                 Face(([x, xr],[yl, y],[x+yl, xl+y+1/q])), \
                 Face(([x, xr],[yl, y],[x+yl+1/q, xl+y+2/q])), \
                 ]
    return faces

def directly_covered_by_adding_face(last_covered_intervals, face):
    """
    Compute incrementally new covered_intervals by adding a new face.
    Consider only directly covered and symmetry reflection regarding f.

    EXAMPLES::

        sage: last_covered_intervals = [[[0, 1/5], [2/5, 3/5]], [[3/5, 4/5], [4/5, 1]]]
        sage: face = Face(([1/5, 2/5], [1/5, 2/5], [3/5, 4,5]))
        sage: directly_covered_by_adding_face(last_covered_intervals, face)
        [[[0, 1/5], [2/5, 3/5]], [[1/5, 2/5], [3/5, 4/5], [4/5, 1]]]
    """
    covered_intervals = copy(last_covered_intervals)
    component = []
    for int1 in face.minimal_triple:
        [x, y] = interval_mod_1(int1)
        component.append([x, y])
        # consider the symmetry reflection.
        if y <= f:
            component.append([f-y, f-x])
        else:
            component.append([1+f-y, 1+f-x])
    component.sort()
    component = merge_within_comp(component)
    covered_intervals.append(component)
    #remove_duplicate(covered_intervals)
    for i in range(len(covered_intervals)-1):
        if find_interior_intersection(covered_intervals[i], covered_intervals[-1]):
                covered_intervals[-1] = merge_two_comp(covered_intervals[i],covered_intervals[-1])
                covered_intervals[i] = []
    covered_intervals = remove_empty_comp(covered_intervals)
    return covered_intervals

######
#global picked_face_set
#picked_face_set = set([])
#npick = -1
######

def paint_complex(q, f, last_additive_vertices, last_faces_set, last_covered_intervals, \
                  candidate_faces, last_non_candidate):
    """
    Randomly paint triangles green in a 2d-complex, until all intervals are covered.
    Return additive_vertices, green_faces and covered_intervals if a possible way of painting is found,
    otherwise, return False

    EXAMPLES::

        sage: q = 5; f = 3/5; num_of_slopes = 2;
        sage: additive_vertices, faces_set, covered_intervals = \
        ...     initial_vertices_faces_and_covered_intervals(q, f)
        sage: candidate_faces = generate_candidate_faces(q, f, covered_intervals, None)
        sage: additive_vertices, green_faces, covered_intervals = \
        ...     paint_complex(q, f, additive_vertices, faces_set, covered_intervals, candidate_faces, set([])).next()
        sage: plot_painted_faces(q, green_faces)
    """
    non_candidate = copy(last_non_candidate)
    num_candidate_faces = len(candidate_faces)
    n = 0
    ### debug
    #global npick
    #npick += 1
    ###
    while n < num_candidate_faces:
        face_picked = candidate_faces[n]
        n += 1
        ###### debug...
        #print face_picked
        #picked_face_set.add(face_picked)
        ###### ...
        faces_set = copy(last_faces_set)
        faces_set.add(face_picked)
        additive_vertices = copy(last_additive_vertices)
        covered_intervals = directly_covered_by_adding_face(last_covered_intervals, face_picked)
        legal_picked = True
        for v in face_picked.vertices:
            if (v[0] <= v[1]) and (not v in additive_vertices) and legal_picked:
                additive_vertices.add(v)
                for face in faces_around_vertex(q, v):
                    if not face in faces_set and \
                           additive_vertices.issuperset([(x,y) for (x,y) in face.vertices if x <= y]):
                        if face in non_candidate:
                            legal_picked = False
                            break
                        else:
                            faces_set.add(face)
                            covered_intervals = directly_covered_by_adding_face(covered_intervals, face)
        if not legal_picked:
            non_candidate.add(face_picked)
            continue

        global to_check_implied
        to_check_implied += 1
        if to_check_implied == check_implied_period:
            to_check_implied = 0
            ### debug
            #nface = len(faces_set)
            #ntocover =len(generate_to_cover(q, covered_intervals))
            ###
            implied_additive_vertices = generate_implied_additive_vertices(q, f, additive_vertices, covered_intervals)
            if implied_additive_vertices is False:
                # implied_additive_vertices is False means infeasible or too few slopes.
                # continue to next face in candidate_faces
                legal_picked = False
            else:
                for v in implied_additive_vertices:
                    if (v[0] <= v[1]) and (not v in additive_vertices) and legal_picked:
                        additive_vertices.add(v)
                        for face in faces_around_vertex(q, v):
                            if not face in faces_set and \
                                   additive_vertices.issuperset([(x,y) for (x,y) in face.vertices if x <= y]):
                                if face in non_candidate:
                                    legal_picked = False
                                    break
                                else:
                                    faces_set.add(face)
                                    covered_intervals = directly_covered_by_adding_face(covered_intervals, face)
            ###debug
            #for i in range(npick):
            #    print '',
            #if implied_additive_vertices is False:
            #    print 'stop',
            #else:
            #    print len(implied_additive_vertices),
            #print len(faces_set)-nface,
            #print ntocover - len(generate_to_cover(q, covered_intervals))
            ###
        # If infeasible or encounter non_candidate or
        # only few slopes are possible given current covered_intervals, stop recursion.
        if legal_picked and num_slopes_at_best(q, covered_intervals) >= num_of_slopes:
            new_candidate_faces= generate_candidate_faces(q, f, covered_intervals, face_picked)
            if not new_candidate_faces:
                # stop recursion
                if not generate_to_cover(q, covered_intervals):
                    # all covered, finish
                    ######### debug...
                    #print "    Painting exists for q = %s, f = %s" % (q, f)
                    #plot_painted_faces(q, faces_set).show(show_legend=False)
                    #print picked_face_set
                    ######### ...
                    yield additive_vertices, faces_set, covered_intervals
            else:
                for additive_vertices, faces_set, covered_intervals in \
                        paint_complex(q, f, additive_vertices, faces_set, covered_intervals, new_candidate_faces, non_candidate):
                    yield additive_vertices, faces_set, covered_intervals
        non_candidate.add(face_picked)
    ###debug
    #npick -= 1
    ###

def num_slopes_at_best(q, covered_intervals):
    """
    Return an upper bound on the final number of slopes given current covered_intervals.

    EXAMPLES::

        sage: covered_intervals = [[[0, 1/5], [2/5, 3/5]], [[3/5, 4/5], [4/5, 1]]]
        sage: num_slopes_at_best(5, covered_intervals)
        3
    """
    uncovered_num = q - sum([len(component) for component in covered_intervals])
    return uncovered_num + len(covered_intervals)

def generate_candidate_faces(q, f, covered, last_face=None):
    """
    Return a set of candidate_faces (lexicographically > last_face, not in non_candidate)
    to paint in next step, whose I, J are currently uncovered.
    Note that candidate_faces only takes faces with I <= J.

    EXAMPLES::

        sage: q = 5; f = 3/5; num_of_slopes = 2;
        sage: additive_vertices, faces_set, covered_intervals = \
        ...     initial_vertices_faces_and_covered_intervals(q, f)
        sage: candidate_faces = generate_candidate_faces(q, f, covered_intervals, None)
        sage: candidate_faces
        set([<Face ([1/5, 2/5], [1/5, 2/5], [2/5, 3/5])>, <Face ([1/5, 2/5], [1/5, 2/5], [3/5, 4/5])>])
    """
    to_cover = generate_to_cover(q, covered)
    # NOTE: candidate_faces only takes faces with x <= y
    candidate_faces = []
    for x in to_cover:
        for y in to_cover:
            if x <= y:
                #if (x + y) in to_cover:
                face = Face(([x, x + 1/q], [y, y + 1/q], [x + y, x + y + 1/q]))
                if (last_face is None or face > last_face):
                    candidate_faces.append(face)
                #if (x + y + 1/q) in to_cover:
                face = Face(([x, x + 1/q], [y, y + 1/q], [x + y + 1/q, x + y + 2/q]))
                if (last_face is None or face > last_face):
                    candidate_faces.append(face)
    return candidate_faces

def generate_to_cover(q, covered):
    """
    Return a set {k/q | 0 <=k < q, [k/q, (k+1)/q] is uncovered}

    EXAMPLES::

        sage: covered_intervals = [[[0, 1/5], [2/5, 3/5]], [[3/5, 4/5], [4/5, 1]]]
        sage: generate_to_cover(5, covered_intervals)
        [1/5]
    """
    to_cover = set([x/q for x in range(q)])
    for component in covered:
        for i in component:
            for x in range(i[0]*q, i[1]*q):
                to_cover.discard(x/q)
    return sorted(list(to_cover))

def generate_implied_additive_vertices(q, f, additive_vertices, covered_intervals):
    """
    Return new additive_vertices implied by the known additive_vertices

    EXAMPLES::

        sage: q = 5; f = 3/5; num_of_slopes = 2;
        sage: additive_vertices = set([(0, 4/5), (0, 1), (1/5, 2/5), (0, 0), (2/5, 1), (3/5, 1), \
        ...    (4/5, 1), (0, 2/5), (1, 1), (0, 3/5), (0, 1/5), (2/5, 2/5), (1/5, 1), (4/5, 4/5)])
        sage: covered_intervals = [[[0, 1/5], [2/5, 3/5]], [[1/5, 2/5], [3/5, 4/5], [4/5, 1]]]
        sage: generate_implied_additive_vertices(q, f, additive_vertices, covered_intervals)
        set([(2/5, 4/5)])
        sage: num_of_slopes = 3;
        sage: generate_implied_additive_vertices(q, f, additive_vertices, covered_intervals)
        False
    """
    to_cover = generate_to_cover(q, covered_intervals)
    uncovered_components = [[[x, x + 1/q]] for x in to_cover]
    components = covered_intervals  + uncovered_components
    # stop backtracking if there are too fews slopes.
    if len(components) < num_of_slopes:
        return False
    fn_sym = generate_symbolic_continuous(None, components, field=QQ)
    ieqdic, eqndic = generate_ieqs_and_eqns(q, f, fn_sym, additive_vertices)
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

def generate_ieqs_and_eqns(q, f, fn_sym, additive_vertices):
    """
    Return the equalities (by additivity) and inequalities (by subadditivity)
    that the slope variables must satisfy.

    Inputs:
        q, f*q are integers

        fn_sym is the symbolic function generated by considering covered_intervals:
        intervals in the same component of a function must have the same slope value.
        Take slope of the i-th component as the i-th unit vector. 
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
    eqn = tuple([-1]) + tuple(fn_sym(f))
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

def generate_vertex_function(q, f, fn_sym, additive_vertices):
    """
    Generate real valued functions which correspond to vertices 
    of the polytope defined by [ieqs, eqns] = generate_ieqs_and_eqns(..)
    """
    ieqdic, eqndic = generate_ieqs_and_eqns(q, f, fn_sym, additive_vertices)
    p = Polyhedron(ieqs = ieqdic.keys(), eqns = eqndic.keys())
    if not p.is_empty():
        #if p.n_vertices() > 1:
        #    print p.Vrepresentation()
        ## print "boundedness is %s" % p.is_compact()
        for x in p.vertices():
            k = len(set(x))
            if k >= num_of_slopes:
                print "%s gives a %s-slope function h =" % (x, k)
                v = vector(QQ,x)
                yield v * fn_sym
            ####### debug...
            #else:
            #    print "%s gives a %s-slope function. Ignore" % (x, k)
            #    # this will probably print many lines
            ####### ...
    #else:
        # this shouldn't happen! 
        # This can happen if we only check generate_implied_additive_vertices once in a while
        #print "p.is_empty() is True"

def search_6slope_example(q, f):
    last_additive_vertices, last_faces_set, last_covered_intervals = initial_vertices_faces_and_covered_intervals(q, f)
    candidate_faces = generate_candidate_faces(q, f, last_covered_intervals, None)
    found = False
    for additive_vertices, green_faces, covered_intervals in \
            paint_complex(q, f, last_additive_vertices, last_faces_set, last_covered_intervals, candidate_faces, set([])):
        if len(covered_intervals) >= num_of_slopes:
            # otherwise, too few slopes in covered_intervals, continue to the next attempt.
            fn_sym = generate_symbolic_continuous(None, covered_intervals, field=QQ)
            for h in generate_vertex_function(q, f, fn_sym, additive_vertices):
                print h
                found = True
                yield h #return h
    if not found:
        print "Example function not found. Please try again."
