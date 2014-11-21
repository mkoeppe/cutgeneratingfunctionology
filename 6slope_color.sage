# call search_6slope_example(q, f) in 6slope_color.sage
# find a fulldim_covers_6slope_extreme functions
def fulldim_covers_6slope_q25_1():
    # q = 25; f = 8
    bkpt = [0, 1/25, 2/25, 6/25, 7/25, 8/25, 9/25, 2/5, 11/25, 12/25, \
            13/25, 3/5, 16/25, 17/25, 18/25, 4/5, 21/25, 22/25, 23/25, 24/25, 1]
    values = [0, 17/36, 1/4, 3/4, 19/36, 1, 37/144, 5/24, 121/288, 107/288, 7/12, \
                35/72, 19/72, 53/72, 37/72, 5/12, 181/288, 167/288, 19/24, 107/144, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def fulldim_covers_6slope_q26_1():
    # q = 26, f = 13
    bkpt = [0, 1/26, 3/26, 5/26, 3/13, 7/26, 4/13, 5/13, 6/13, 1/2, \
            7/13, 8/13, 9/13, 19/26, 10/13, 21/26, 23/26, 25/26, 1]
    values = [0, 3/7, 5/7, 1/7, 2/7, 5/7, 6/7, 2/7, 4/7, 1, \
              4/7, 2/7, 6/7, 5/7, 2/7, 1/7, 5/7, 3/7, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

# Similar 2d-diagram
def fulldim_covers_6slope_q38_1():
    # q = 38, f = 19
    bkpt = [0, 1/38, 5/38, 7/38, 9/38, 5/19, 6/19, 7/19, 9/19, 1/2, \
            10/19, 12/19, 13/19, 14/19, 29/38, 31/38, 33/38, 37/38, 1]
    values = [0, 7/17, 11/17, 3/17, 5/17, 12/17, 14/17, 6/17, 10/17, 1, \
            10/17, 6/17, 14/17, 12/17, 5/17, 3/17, 11/17, 7/17, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

# backtracking search for 6-slope extreme functions, 
# using incremental computation for polytope, vertices_color, faces_color and covered_intervals.
# do not consider translation/reflection other than the symmetry reflection regarding f.

# polyhedral computatoin library:
# http://www.sagemath.org/doc/reference/libs/sage/libs/ppl.html#sage.libs.ppl.Polyhedron.minimize

# work with array and integer variables. 
# do NOT use class Face in functions.sage
# do NOT call functions, except for piecewise_function_from_breakpoints_and_values() and plot_..(), from functions.sage

# q and f are integers.
# vertices_color (q+1)*(q+1) 0-1 array, 0: green, 1: unknown, currently white, 2: non_candidate, must be white)
# faces_color q*q*2 0-1-2 array, 0: green, 1: unknown, currently white, 2: non-candidate, must be white.
# vertices have integer coordinates, k is integer in function value fn(k).
# face is represented by 3 integers (x, y, w), 0 <= x, y < q, w=0(lower triangle) or 1 (upper triangle)
# covered_intervals is a list of sets. For example: [set([0, 2]), set([3, 4])] means:
# [0,1] and [0,2] are covered with slope s1, [3,4] and [4,5] are covered with slope s2.
# only record vertex (x,y) and face = (x, y, w) with x <= y.
# polytope defines the feasible region of (\pi(0), \pi(1/q),..., \pi(1)).

from sage.libs.ppl import C_Polyhedron, Constraint, Constraint_System, Generator, Generator_System, Variable, point, Poly_Con_Relation
import numpy

#global num_of_slopes
num_of_slopes = 6 # set to 6 as we are looking for 6-slope functions

def initial_vertices_color(q, f):
    """
    paint green ( = 0) for vertices (x<=y) corresonding to
    \pi(0) = 0, \pi(1) = 0 and reflection around f/2 and (1+f)/2.

    EXAMPLES::

        sage: initial_vertices_color(5, 3)
        array([[0, 0, 0, 0, 0, 0],
               [1, 1, 0, 1, 1, 0],
               [1, 1, 1, 1, 1, 0],
               [1, 1, 1, 1, 1, 0],
               [1, 1, 1, 1, 0, 0],
               [1, 1, 1, 1, 1, 0]])
    """
    vertices_color = numpy.ones((q+1,q+1),int)
    # color: 1-white; 0-green
    # border x = 0 and border y = 0 are green
    for x in range(q+1):
        vertices_color[0, x] = 0
        vertices_color[x, q] = 0
    # diagonals corresponding to f
    for x in range(q+1):
        if x <= f/2:
            vertices_color[x, f - x] = 0
        elif (x >= f) and (x <= f - x + q):
            vertices_color[x, f - x + q] = 0
    # implied additive vertices
    cs = initial_cs(q, f, vertices_color)
    polytope = C_Polyhedron(cs)
    for x in range(1, q):
        for y in range(x, q):
            if (vertices_color[x, y] == 1):
                    z = (x + y) % q
                    if polytope.relation_with( Variable(x) + Variable(y) == Variable(z) ).implies(Poly_Con_Relation.is_included()):
                        # find implied additive vertices
                        vertices_color[x, y] = 0
    return vertices_color

def initial_faces_color_and_covered_intervals(q, f, vertices_color):
    """
    Return initial faces_color and covered_intervals,
    corresponding to \pi(0) = 0, \pi(1) = 0 and reflection around f/2 and (1+f)/2.

    EXAMPLES::

        sage: q=5; f=3;
        sage: vertices_color = initial_vertices_color(q, f);
        sage: faces_color, covered_intervals = initial_faces_color_and_covered_intervals(q, f, vertices_color)
        sage: faces_color[:,:,0]
        array([[0, 1, 0, 1, 1],
               [1, 1, 1, 1, 1],
               [1, 1, 1, 1, 1],
               [1, 1, 1, 1, 1],
               [1, 1, 1, 1, 0]])
        sage: faces_color[:,:,1]
        array([[1, 1, 1, 1, 1],
               [1, 1, 1, 1, 1],
               [1, 1, 1, 1, 1],
               [1, 1, 1, 1, 0],
               [1, 1, 1, 1, 0]])
        sage: covered_intervals
        [set([0, 2]), set([3, 4])]
    """
    faces_color = numpy.ones((q, q, 2),int)
    covered_intervals = []
    for x in range(q):
        for y in range(x, q):
            for w in range(2):
                if x < y:
                    vertices = [(x+1, y), (x, y+1), (x+w, y+w)]
                else:
                    vertices = [(x, y+1), (x+w, y+w)]
                if sum(vertices_color[v] for v in vertices) == 0 :
                    face = (x, y, w)
                    faces_color[face] = 0
                    covered_intervals = directly_covered_by_adding_face(covered_intervals, face, q, f)
    return faces_color, covered_intervals

def initial_cs(q, f, vertices_color):
    """
    Return the initial Constraint_System that defines 
    the feasible region of (\pi(0), \pi(1/q),...,\pi(1))

    EXAMPLES::

        sage: q=5; f=3;
        sage: vertices_color = initial_vertices_color(q, f);
        sage: cs = initial_cs(q, f, vertices_color)
        sage: cs
        Constraint_System {x0==0, x5==0, x3-1==0, x1>=0, -x1+1>=0,
                           x2>=0, -x2+1>=0, x3>=0, -x3+1>=0, x4>=0, -x4+1>=0,
                           2*x1-x2>=0, x1+x2-x3==0, x1+x3-x4>=0, -x0+x1+x4>=0,
                           2*x2-x4>=0, -x0+x2+x3>=0, -x1+x2+x4>=0, -x1+2*x3>=0,
                           -x2+x3+x4>=0, x3-2*x4==0}
        sage: C_Polyhedron(cs).minimized_generators()
        Generator_System {point(0/4, 3/4, 1/4, 4/4, 2/4, 0/4), 
                          point(0/6, 2/6, 4/6, 6/6, 3/6, 0/6)}
    """
    cs = Constraint_System()
    fn = [ Variable(i) for i in range(q+1) ]

    cs.insert(fn[0] == 0)
    cs.insert(fn[q] == 0)
    cs.insert(fn[f] == 1)
    for i in range(1,q):
        cs.insert(fn[i] >= 0)
        cs.insert(fn[i] <= 1)
    # symmetry is taken care by initial green vertices

    for x in range(1, q):
        for y in range(x, q):
            z = (x + y) % q
            if vertices_color[x, y] == 0:
                cs.insert(fn[x] + fn[y] == fn[z])
            else:
                cs.insert(fn[x] + fn[y] >= fn[z])
    return cs

def edges_around_vertex(q, v):
    """
    Given a grid vertex v (assume that v[0] <= v[1], v is not on the border),
    return a list of elements corresponding to the edges connected to v.
    Each element has the form (set([a, b]), v'). 
    v' is a grid vertex next to v. 
    If the edge (v, v') is green, then segments a and b are connected.

    EXAMPLES::

        sage: edges_around_vertex(5, (1, 3))
        [(set([0, 3]), (0, 3)), (set([3, 4]), (1, 4)), (set([1, 4]), (2, 3)), \
         (set([2, 3]), (1, 2)), (set([0, 3]), (0, 4)), (set([1, 2]), (2, 2))]
        sage: edges_around_vertex(5, (1, 2))
        [(set([0, 2]), (0, 2)), (set([2, 3]), (1, 3)), (set([1, 3]), (2, 2)), \
         (set([1, 2]), (1, 1)), (set([0, 2]), (0, 3))]
        sage: edges_around_vertex(5, (1, 1))
        [(set([0, 1]), (0, 1)), (set([1, 2]), (1, 2)), (set([0, 1]), (0, 2))]
    """
    (xx, yy) = v
    if xx > yy:
        return []
    # Note: v[0]<=v[1]; only return faces with I <= J
    xl = xx - 1
    xr = xx + 1
    yl = yy - 1
    yr = yy + 1
    if xx < yl:
        return [ \
           (set([ xl, (xl + yy) % q ]), (xl, yy)), \
           (set([ yy, (xx + yy) % q ]), (xx, yr)), \
           (set([ xx, (xx + yy) % q ]), (xr, yy)), \
           (set([ yl, (xx + yl) % q]), (xx, yl)), \
           (set([ xl, yy ]), (xl, yr)), \
           (set([ xx, yl ]), (xr, yl)) ]
    elif xx == yl:
        return [ \
           (set([ xl, (xl + yy) % q ]), (xl, yy)), \
           (set([ yy, (xx + yy) % q ]), (xx, yr)), \
           (set([ xx, (xx + yy) % q ]), (xr, yy)), \
           (set([ yl, (xx + yl) % q]), (xx, yl)), \
           (set([ xl, yy ]), (xl, yr)) ]
    else:
        return [ \
           (set([ xl, (xl + yy) % q ]), (xl, yy)), \
           (set([ yy, (xx + yy) % q ]), (xx, yr)), \
           (set([ xl, yy ]), (xl, yr)) ]

def faces_around_vertex(q, v):
    """
    Given a grid vertex v (assume that v[0] <= v[1], v is not on the border),
    return small triangle faces (only those with x <= y)
    and their vertices (only those with x <= y) that are around v.

    EXAMPLES::

        sage: faces_around_vertex(5, (1, 2))
        [((0, 2, 0), [(0, 2), (1, 2), (0, 3)]), ((0, 2, 1), [(1, 2), (0, 3), (1, 3)]), \
         ((1, 2, 0), [(1, 2), (2, 2), (1, 3)]), ((0, 1, 1), [(1, 1), (0, 2), (1, 2)]), \
         ((1, 1, 0), [(1, 1), (1, 2)]), ((1, 1, 1), [(1, 2), (2, 2)])]
        sage: faces_around_vertex(5, (1, 1))
        [((0, 1, 0), [(0, 1), (1, 1), (0, 2)]), ((0, 1, 1), [(1, 1), (0, 2), (1, 2)]), \
        ((1, 1, 0), [(1, 1), (1, 2)]), ((0, 0, 1), [(0, 1), (1, 1)])]
    """
    (xx, yy) = v
    if xx > yy:
        return []
    # Note: v[0]<=v[1]; only return faces with I <= J
    xl = xx - 1
    xr = xx + 1
    yl = yy - 1
    yr = yy + 1
    if xx < yl:
        return [ \
           ((xl, yy, 0), [(xl, yy), (xx, yy), (xl, yr)]), \
           ((xl, yy, 1), [(xx, yy), (xl, yr), (xx, yr)]), \
           ((xx, yy, 0), [(xx, yy), (xr, yy), (xx, yr)]), \
           ((xl, yl, 1), [(xx, yl), (xl, yy), (xx, yy)]), \
           ((xx, yl, 0), [(xx, yl), (xr, yl), (xx, yy)]), \
           ((xx, yl, 1), [(xx, yy),(xr, yl),(xr, yy)]) ]
    elif xx == yl:
        return [ \
           ((xl, yy, 0), [(xl, yy), (xx, yy), (xl, yr)]), \
           ((xl, yy, 1), [(xx, yy), (xl, yr), (xx, yr)]), \
           ((xx, yy, 0), [(xx, yy), (xr, yy), (xx, yr)]), \
           ((xl, yl, 1), [(xx, yl), (xl, yy), (xx, yy)]), \
           ((xx, yl, 0), [(xx, yl), (xx, yy)]), \
           ((xx, yl, 1), [(xx, yy),(xr, yy)]) ]
    else:
        return [ \
           ((xl, yy, 0), [(xl, yy), (xx, yy), (xl, yr)]), \
           ((xl, yy, 1), [(xx, yy), (xl, yr), (xx, yr)]), \
           ((xx, yy, 0), [(xx, yy), (xx, yr)]), \
           ((xl, yl, 1), [(xl, yy), (xx, yy)])]

def directly_covered_by_adding_face(last_covered_intervals, face, q, f):
    """
    Compute incrementally new covered_intervals by adding a new face.
    Consider only directly covered and symmetry reflection regarding f.

    EXAMPLES::

        sage: last_covered_intervals = [set([0,2]), set([3,4])]
        sage: face = (1, 1, 1)
        sage: directly_covered_by_adding_face(last_covered_intervals, face, 5, 3)
        [set([0, 2]), set([1, 3, 4])]
    """
    covered_intervals = []
    (x, y, w) = face
    z = (x + y + w) % q 
    newly_covered = set([ x, y, z, (f - x - 1) % q, (f - y - 1) % q, (f - z - 1) % q ])
    for component in last_covered_intervals:
        if component & newly_covered:
            newly_covered.update(component)
        else:
            covered_intervals.append(component)
    covered_intervals.append(newly_covered)
    return covered_intervals

def update_covered_uncovered_by_adding_face(last_covered_intervals, last_uncovered_intervals, face, q):
    """
    Compute incrementally new covered_intervals and new uncovered_intervals 
    by adding a new green triangle.

    EXAMPLES::

        sage: last_covered_intervals = [set([0,2]), set([3,4])]
        sage: last_uncovered_intervals = [set([1])]
        sage: face = (1, 1, 1)
        sage: update_covered_uncovered_by_adding_face(last_covered_intervals, last_uncovered_intervals, face, 5)
        ([set([0, 2]), set([1, 3, 4])], [])
    """
    covered_intervals = []
    uncovered_intervals = []
    (x, y, w) = face
    z = (x + y + w) %q
    newly_covered = set([x, y, z]) # note that reflection regarding f is considered in initial_covered_uncovered()
    for component in last_uncovered_intervals:
        if component & newly_covered:
            newly_covered.update(component)
        else:
            uncovered_intervals.append(component)
    for component in last_covered_intervals:
        if component & newly_covered:
            newly_covered.update(component)
        else:
            covered_intervals.append(component)
    covered_intervals.append(newly_covered)
    return covered_intervals, uncovered_intervals

def update_covered_uncovered_by_adding_edge(last_covered_intervals, last_uncovered_intervals, to_merge_set, q):
    """
    Compute incrementally new covered_intervals and new uncovered_intervals 
    by adding a new green edge, which implies that the elements in to_merge_set must be connected.

    EXAMPLES::

        sage: last_covered_intervals = [set([0,2]), set([3,4])]
        sage: last_uncovered_intervals = [set([1])]
        sage: to_merge_set = set([0, 1])
        sage: update_covered_uncovered_by_adding_edge(last_covered_intervals, last_uncovered_intervals, to_merge_set, 5)
        ([set([3, 4]), set([0, 1, 2])], [])
    """
    covered_intervals = []
    uncovered_intervals = []
    new_set_is_covered = False
    for component in last_uncovered_intervals:
        if component & to_merge_set:
            to_merge_set.update(component)
        else:
            uncovered_intervals.append(component)
    for component in last_covered_intervals:
        if component & to_merge_set:
            to_merge_set.update(component)
            new_set_is_covered = True
        else:
            covered_intervals.append(component)
    if new_set_is_covered:
        covered_intervals.append(to_merge_set)
    else:
        uncovered_intervals.append(to_merge_set)
    return covered_intervals, uncovered_intervals

def generate_to_cover(q, covered_intervals):
    """
    Return a sorted list {k | 0 <=k < q, [k, (k+1)] is uncovered}

    EXAMPLES::

        sage: covered_intervals = [set([0,2]), set([3,4])]
        sage: generate_to_cover(5, covered_intervals)
        [1]
    """
    to_cover = set(range(q))
    for component in covered_intervals:
        to_cover -= component
    return sorted(list(to_cover))

def generate_candidate_faces(q, f, covered_intervals, last_face=None):
    """
    Return a list of candidate_faces (lexicographically > last_face)
    to paint in next step, whose I, J are currently uncovered.
    Note that candidate_faces only takes faces with I <= J.

    EXAMPLES::

        sage: q = 5; f = 3; num_of_slopes = 2;
        sage: vertices_color = initial_vertices_color(q, f);
        sage: faces_color, covered_intervals = initial_faces_color_and_covered_intervals(q, f, vertices_color)
        sage: generate_candidate_faces(q, f, covered_intervals, None)
        [(1, 1, 0), (1, 1, 1)]
    """
    to_cover = generate_to_cover(q, covered_intervals)
    # NOTE: candidate_faces only takes faces with x <= y
    candidate_faces = []
    for x in to_cover:
        for y in to_cover:
            if x <= y:
                face = (x, y, 0)
                if (last_face is None or face > last_face):
                    candidate_faces.append(face)
                face = (x, y, 1)
                if (last_face is None or face > last_face):
                    candidate_faces.append(face)
    return candidate_faces

def num_slopes_at_best(q, covered_intervals, uncovered_intervals=None):
    """
    Return an upper bound on the final number of slopes,
    given current covered_intervals (and optionally, uncovered_intervals 
    that provides connected components of non-covered intervals).

    EXAMPLES::

        sage: covered_intervals = [set([0,2])]
        sage: uncovered_intervals = [set([1]), set([3,4])]
        sage: num_slopes_at_best(5, covered_intervals)
        4
        sage: num_slopes_at_best(5, covered_intervals, uncovered_intervals)
        3
    """
    if uncovered_intervals is None:
        uncovered_num = q - sum([len(component) for component in covered_intervals])
    else:
        uncovered_num = len(uncovered_intervals)
    return uncovered_num + len(covered_intervals)

def paint_complex_heuristic(q, f, last_vertices_color, last_faces_color, last_covered_intervals, candidate_faces, last_cs):
    """
    Paint triangles green in a 2d-complex, until all intervals are covered.
    Return the polytope which defines the feasible region of (\pi(0), \pi(1/q),...,\pi(1)) for that paint_complex_heuristic

    Heuristic: 
        I, J projections of the next face to paint are chosen from currently uncovered intervals.
        Stop painting immediately once everything is covered.
        Edges are not considered.

    EXAMPLES::

        sage: q = 5; f = 3; num_of_slopes = 2;
        sage: vertices_color = initial_vertices_color(q, f);
        sage: faces_color, covered_intervals = initial_faces_color_and_covered_intervals(q, f, vertices_color)
        sage: cs = initial_cs(q, f, vertices_color)
        sage: candidate_faces = generate_candidate_faces(q, f, covered_intervals, None)
        sage: for result_polytope in paint_complex_heuristic(q, f, vertices_color, faces_color, covered_intervals, candidate_faces, cs):
        ...       print result_polytope.minimized_generators()
        Generator_System {point(0/6, 2/6, 4/6, 6/6, 3/6, 0/6)}
        Generator_System {point(0/4, 3/4, 1/4, 4/4, 2/4, 0/4)}
    """
    num_candidate_faces = len(candidate_faces)
    n = 0
    while n < num_candidate_faces:
        faces_color = copy(last_faces_color)
        for i in range(n): 
            # set color =2 for non_candidate faces
            faces_color[candidate_faces[i]] = 2
        face_picked = candidate_faces[n]
        n += 1
        faces_color[face_picked] = 0
        vertices_color = copy(last_vertices_color)
        covered_intervals = directly_covered_by_adding_face(last_covered_intervals, face_picked, q, f)
        cs = copy(last_cs)

        legal_picked = True
        (x, y, w) = face_picked
        if x < y:
            vertices_picked = [(x+1, y), (x, y+1), (x+w, y+w)]
        else:
            vertices_picked = [(x, y+1), (x+w, y+w)]
        for (x, y) in vertices_picked:
            if (vertices_color[x, y] == 1) and legal_picked:
                # new green vertice
                vertices_color[x, y] = 0
                # update cs
                z = (x + y) % q
                cs.insert( Variable(x) + Variable(y) == Variable(z) )
                for (face, vertices) in faces_around_vertex(q, (x, y)):
                    if faces_color[face] != 0 and sum(vertices_color[v] for v in vertices) == 0:
                        # find new green face.
                        if faces_color[face] == 2: # face is in non_candidate
                            legal_picked = False
                            break
                        else:
                            faces_color[face] = 0
                            covered_intervals = directly_covered_by_adding_face(covered_intervals, face, q, f)
        polytope = C_Polyhedron(cs)
        if not legal_picked or num_slopes_at_best(q, covered_intervals) < num_of_slopes or polytope.is_empty():
            # If encounter non_candidate or too few slopes, stop recursion or infeasible.
            #faces_color[face_picked] = 2
            continue
        # look for implied additive vertices
        for x in range(1, q):
            for y in range(x, q):
                if legal_picked and (vertices_color[x, y] == 1):
                    z = (x + y) % q
                    if polytope.relation_with( Variable(x) + Variable(y) == Variable(z) ).implies(Poly_Con_Relation.is_included()):
                        # find implied additive vertices
                        vertices_color[x, y] = 0
                        for (face, vertices) in faces_around_vertex(q, (x, y)):
                            if faces_color[face] != 0 and sum(vertices_color[v] for v in vertices) == 0:
                                # find new green face.
                                if faces_color[face] == 2: # face is in non_candidate
                                    legal_picked = False
                                    break
                                else:
                                    faces_color[face] = 0
                                    covered_intervals = directly_covered_by_adding_face(covered_intervals, face, q, f)
        # If infeasible or encounter non_candidate or
        # only few slopes are possible given current covered_intervals, stop recursion.
        if legal_picked and num_slopes_at_best(q, covered_intervals) >= num_of_slopes:
            new_candidate_faces= generate_candidate_faces(q, f, covered_intervals, face_picked)
            if not new_candidate_faces:
                # stop recursion
                if not generate_to_cover(q, covered_intervals):
                    # all covered, finish
                    yield polytope
            else:
                for result_polytope in paint_complex_heuristic(q, f, vertices_color, faces_color, \
                        covered_intervals, new_candidate_faces, polytope.constraints()):
                    yield result_polytope
        #faces_color[face_picked] = 2

def paint_complex_fulldim_covers(q, f, last_vertices_color, last_faces_color, last_covered_intervals, (x, y, w), last_cs):
    """
    Paint triangles green in a 2d-complex, until all possibilities are tried.
    If all intervals are covered, 
    return the polytope which defines the feasible region of (\pi(0), \pi(1/q),...,\pi(1)) for that paint_complex_fulldim_covers

    fulldim_covers: 
        Allow painting any face, even if covered already.
        Edges are not considered.

    EXAMPLES::

        sage: q = 5; f = 3; num_of_slopes = 2;
        sage: vertices_color = initial_vertices_color(q, f);
        sage: faces_color, covered_intervals = initial_faces_color_and_covered_intervals(q, f, vertices_color)
        sage: cs = initial_cs(q, f, vertices_color)
        sage: for result_polytope in paint_complex_fulldim_covers(q, f, vertices_color, faces_color, covered_intervals, (0, 0, 0), cs):
        ...       print result_polytope.minimized_generators()
        Generator_System {point(0/6, 2/6, 4/6, 6/6, 3/6, 0/6)}
        Generator_System {point(0/4, 3/4, 1/4, 4/4, 2/4, 0/4)}
    """
    if (x, y, w) == (q - 1, q - 1, 1):
        # finish painting, check if all intervals are covered
        if not generate_to_cover(q, last_covered_intervals):
            # all covered, return valid painting
            polytope = C_Polyhedron(last_cs)
            if not polytope.is_empty():
                yield polytope
    else:
        # move to the next triangle
        w += 1
        if w == 2:
            w = 0
            y += 1
        if y == q:
            x += 1
            y = x
        # face_picked = (x, y, w)
        faces_color = copy(last_faces_color)
        if faces_color[(x, y, w)] == 1: # color is unkown
            vertices_color = copy(last_vertices_color)
            cs = copy(last_cs)
            # First, try out green triangle (x, y, w)
            faces_color [(x, y, w)] = 0
            covered_intervals = directly_covered_by_adding_face(last_covered_intervals, (x, y, w), q, f)
            legal_picked = True
            if x < y:
                vertices_picked = [(x+1, y), (x, y+1), (x+w, y+w)]
            else:
                vertices_picked = [(x, y+1), (x+w, y+w)]
            for (i, j) in vertices_picked:
                if (vertices_color[i, j] == 1) and legal_picked:
                    # new green vertice
                    vertices_color[i, j] = 0
                    # update cs
                    k = (i + j) % q
                    cs.insert( Variable(i) + Variable(j) == Variable(k) )
                    for (face, vertices) in faces_around_vertex(q, (i, j)):
                        if faces_color[face] != 0 and sum(vertices_color[v] for v in vertices) == 0:
                            # find new green face.
                            if faces_color[face] == 2: # face is in non_candidate
                                legal_picked = False
                                break
                            else:
                                faces_color[face] = 0
                                covered_intervals = directly_covered_by_adding_face(covered_intervals, face, q, f)
            if legal_picked and num_slopes_at_best(q, covered_intervals) >= num_of_slopes:
                # If encounter non_candidate or too few slopes, stop recursion
                polytope = C_Polyhedron(cs)
                if not polytope.is_empty():
                    # If infeasible, stop recursion
                    # look for implied additive vertices
                    for i in range(1, q):
                        for j in range(i, q):
                            if legal_picked and (vertices_color[i, j] == 1):
                                k = (i + j) % q
                                if polytope.relation_with( Variable(i) + Variable(j) == Variable(k) \
                                                            ).implies(Poly_Con_Relation.is_included()):
                                    # find implied additive vertices
                                    vertices_color[i, j] = 0
                                    for (face, vertices) in faces_around_vertex(q, (i, j)):
                                        if faces_color[face] != 0 and sum(vertices_color[v] for v in vertices) == 0:
                                            # find new green face.
                                            if faces_color[face] == 2: # face is in non_candidate
                                                legal_picked = False
                                                break
                                            else:
                                                faces_color[face] = 0
                                                covered_intervals = directly_covered_by_adding_face(covered_intervals, face, q, f)
                    # If infeasible or encounter non_candidate or
                    # only few slopes are possible after considering implied things, stop recursion.
                    if legal_picked and num_slopes_at_best(q, covered_intervals) >= num_of_slopes:
                        for result_polytope in paint_complex_fulldim_covers(q, f, \
                                vertices_color, faces_color, covered_intervals, (x, y, w), cs):
                            yield result_polytope
            # Now, try out white triangle (x, y, w)
            faces_color = copy(last_faces_color)
            faces_color [(x, y, w)] = 2
        for result_polytope in paint_complex_fulldim_covers(q, f, \
                last_vertices_color, faces_color, last_covered_intervals, (x, y, w), last_cs):
            yield result_polytope

def initial_covered_uncovered(q, f, vertices_color):
    """
    Return initial covered_intervals and uncovered_intervals,
    corresponding to \pi(0) = 0, \pi(1) = 0 and reflection around f/2 and (1+f)/2.

    EXAMPLES::

        sage: q=5; f=3;
        sage: vertices_color = initial_vertices_color(q, f);
        sage: initial_covered_uncovered(q, f, vertices_color)
        ([set([0, 2]), set([3, 4])], [set([1])])
    """
    # lower-left and upper-right green triangle
    covered_intervals = [set([0, f-1]), set([f, q-1])]
    # consider reflection regarding f.
    uncovered_intervals = [set([i, f - i - 1]) for i in range(1, (f + 1) // 2)]
    uncovered_intervals += [set([i, q + f - i - 1]) for i in range(f + 1, (f + q + 1) // 2)]
    for x in range(1, q):
        for y in range(x, q):
            if vertices_color[(x, y)] == 0:
                covered_intervals, uncovered_intervals = update_around_green_vertex(q, (x, y), \
                                        vertices_color, covered_intervals, uncovered_intervals)
    return covered_intervals, uncovered_intervals

def update_around_green_vertex(q, (x, y), vertices_color, covered_intervals, uncovered_intervals):
    """
    Painting vertex (x, y) from white to green induces some new green triangles and edges around this vertex.
    Return new covered_intervals and uncovered_intervals corresponding to these changes.

    EXAMPLES::

        sage: q = 5; f = 3;
        sage: vertices_color = initial_vertices_color(q, f);
        sage: (covered_intervals, uncovered_intervals) = ([set([0, 2]), set([3, 4])], [set([1])]);
        sage: x = 1; y = 1;  vertices_color[x, y] = 0;
        sage: update_around_green_vertex(q, (x, y), vertices_color, covered_intervals, uncovered_intervals)
        ([set([3, 4]), set([0, 1, 2])], [])
        sage: vertices_color = initial_vertices_color(q, f);
        sage: x = 2; y = 2; vertices_color[x,y]=0
        sage: update_around_green_vertex(q, (x, y), vertices_color, covered_intervals, uncovered_intervals)
        ([set([0, 2]), set([1, 3, 4])], [])
    """
    deja_v = set([])
    for (face, vertices) in faces_around_vertex(q, (x, y)):
        if sum(vertices_color[v] for v in vertices) == 0:
            # find new green face.
            deja_v.update(vertices)
            covered_intervals, uncovered_intervals = \
                update_covered_uncovered_by_adding_face(covered_intervals, uncovered_intervals, face, q)
    for (to_merge_set, v) in edges_around_vertex(q, (x, y)):
        # if an green edge is included in a green face, don't need to update_covered_uncovered again.
        if vertices_color[v] == 0 and not v in deja_v:
            # find new green edge.
            covered_intervals, uncovered_intervals = \
                update_covered_uncovered_by_adding_edge(covered_intervals, uncovered_intervals, to_merge_set, q)
    return covered_intervals, uncovered_intervals

def paint_complex_complete(q, f, last_vertices_color, last_covered_intervals, last_uncovered_intervals, (x, y), last_cs):
    """
    Paint vertices green in a 2d-complex, until all possibilities are tried.
    If all intervals are covered (consider both green triangles and edges),
    return the polytope which defines the feasible region of (\pi(0), \pi(1/q),...,\pi(1)) for that paint_complex_complete

    EXAMPLES::
        sage: q = 5; f = 3; num_of_slopes = 2;
        sage: vertices_color = initial_vertices_color(q, f);
        sage: cs = initial_cs(q, f, vertices_color)
        sage: covered_intervals, uncovered_intervals = initial_covered_uncovered(q, f, vertices_color)
        sage: for result_polytope in paint_complex_complete(q, f, vertices_color, covered_intervals, uncovered_intervals, (0, 0), cs):
        ...       print result_polytope.minimized_generators()
        Generator_System {point(0/6, 2/6, 4/6, 6/6, 3/6, 0/6)}
        Generator_System {point(0/4, 3/4, 1/4, 4/4, 2/4, 0/4)}
    """
    if (x, y) == (q - 1, q - 1):
        # finish painting, check if all intervals are covered
        if not last_uncovered_intervals:
            # all covered, return valid painting
            polytope = C_Polyhedron(last_cs)
            if not polytope.is_empty():
                yield polytope
    else:
        # move to the next vertices
        y += 1
        if y == q:
            x += 1
            y = x
        # vertex_picked = (x, y)
        vertices_color = copy(last_vertices_color)
        if vertices_color[(x, y)] == 1: # color is unkown
            cs = copy(last_cs)
            covered_intervals = copy(last_covered_intervals)
            uncovered_intervals = copy(last_uncovered_intervals)
            # First, try out green vertex (x, y)
            vertices_color[(x, y)] = 0
            z = (x + y) % q
            cs.insert( Variable(x) + Variable(y) == Variable(z) )
            covered_intervals, uncovered_intervals = update_around_green_vertex(q, (x, y), \
                                    vertices_color, covered_intervals, uncovered_intervals)
            # If too few slopes, stop recursion
            if num_slopes_at_best(q, covered_intervals, uncovered_intervals) >= num_of_slopes:
                polytope = C_Polyhedron(cs)
                # If infeasible, stop recursion
                if not polytope.is_empty():
                    legal_picked = True
                    # look for implied additive vertices
                    for i in range(1, q):
                        for j in range(i, q):
                            if legal_picked and (vertices_color[i, j] != 0):
                                k = (i + j) % q
                                if polytope.relation_with( Variable(i) + Variable(j) == Variable(k) \
                                                            ).implies(Poly_Con_Relation.is_included()):
                                    # find implied additive vertices
                                    if vertices_color[i, j] == 2:
                                        # encounter non_candidate vertex, stop recursion
                                        legal_picked = False
                                        break
                                    else:
                                        # paint this implied vertex green
                                        vertices_color[i, j] = 0
                                        covered_intervals, uncovered_intervals = update_around_green_vertex(q, (i, j), \
                                                                    vertices_color, covered_intervals, uncovered_intervals)
                                        if num_slopes_at_best(q, covered_intervals, uncovered_intervals) < num_of_slopes:
                                            legal_picked = False
                                            break
                    if legal_picked:
                        for result_polytope in paint_complex_complete(q, f, vertices_color, \
                                                covered_intervals, uncovered_intervals, (x, y), cs):
                            yield result_polytope
            # Now, try out white vertex (x, y)
            vertices_color = copy(last_vertices_color)
            vertices_color[(x, y)] = 2
        for result_polytope in paint_complex_complete(q, f, vertices_color, \
                                last_covered_intervals, last_uncovered_intervals, (x, y), last_cs):
            yield result_polytope

def plot_painted_faces(q, faces):
    """
    Return a plot of the 2d complex of q-grid with green faces.

    EXAMPLES::

        sage: q=5; faces = [(1, 1, 0), (1, 1, 1)]
        sage: plot_painted_faces(q, faces)
    """
    points = [x/q for x in range(q+1)]
    values = [0 for x in range(q+1)]
    function = discrete_function_from_points_and_values(points, values)

    p = Graphics()
    p.set_legend_options(loc='upper right')
    p += plot_2d_complex(function)

    kwds = { 'alpha': proj_plot_alpha, 'zorder': -10, 'color': 'grey'}
    IJK_kwds = [ kwds for i in range(3) ]
    for (x, y, w) in faces:
        face = Face(([x/q, (x+1)/q], [y/q, (y+1)/q], [(x+y+w)/q, (x+y+1+w)/q]))
        p += face.plot()
        p += plot_projections_of_one_face(face, IJK_kwds)
    return p

def generate_vertex_function(q, polytope, v_set=set([])):
    """
    Generate real valued functions corresponding to vertices of the polytope.
    Return those functions that have required number of slope values,
    but have not been found before (not in v_set).

    EXAMPLES::

        sage: q=5; f=3; num_of_slopes = 2;
        sage: vertices_color = initial_vertices_color(q, f);
        sage: polytope = C_Polyhedron(initial_cs(q, f, vertices_color))
        sage: h = generate_vertex_function(q, polytope).next()
    """
    bkpt = [i/q for i in range(q+1)]
    for v in polytope.minimized_generators():
        v_n = v.coefficients()
        v_d = v.divisor()
        num = len(set([v_n[i+1] - v_n[i] for i in range(q)]))
        if num >= num_of_slopes:
            values = [v_n[i] / v_d for i in range(q+1)]
            if not tuple(values) in v_set:
                v_set.add(tuple(values))
                h = piecewise_function_from_breakpoints_and_values(bkpt, values)
                yield h

def search_6slope_example(q, f, mode='heuristic', print_function=True):
    """
    Search for extreme functions that have required number of slope values.

    If `mode` is (defaut) 'heuristic', use paint_complex_heuristic() to paint;
    If `mode` is 'fulldim_covers', use paint_complex_fulldim_covers() to paint;
    If `mode` is 'complete', use paint_complex_complete() to paint;

    EXAMPLES::

        sage: q=5; f=3; num_of_slopes = 2;
        sage: h = search_6slope_example(q, f, mode='heuristic', print_function=True).next()
    """
    #initialization
    vertices_color = initial_vertices_color(q, f)
    cs = initial_cs(q, f, vertices_color)
    if mode == 'heuristic':
        faces_color, covered_intervals = initial_faces_color_and_covered_intervals(q, f, vertices_color)
        candidate_faces = generate_candidate_faces(q, f, covered_intervals, None)
        gen = paint_complex_heuristic(q, f, vertices_color, faces_color, covered_intervals, candidate_faces, cs)
    elif mode == 'fulldim_covers':
        faces_color, covered_intervals = initial_faces_color_and_covered_intervals(q, f, vertices_color)
        gen = paint_complex_fulldim_covers(q, f, vertices_color, faces_color, covered_intervals, (0, 0, 0), cs)
    elif mode == 'complete':
        covered_intervals, uncovered_intervals = initial_covered_uncovered(q, f, vertices_color)
        gen = paint_complex_complete(q, f, vertices_color, covered_intervals, uncovered_intervals, (1, 0), cs)
    else:
        raise ValueError, "mode must be one of {'heuristic', 'fulldim_covers', 'complete'."
    v_set = set([])
    for result_polytope in gen:
        for h in generate_vertex_function(q, result_polytope, v_set):
            if print_function:
                print h
            yield h #return h
    if not v_set:
        print "Example function not found. Please try again."

import time
def search_6slope_given_q(q, f_list=None, mode='heuristic', print_function=True):
    """
    EXAMPLES::

        sage: q = 12; num_of_slopes = 3;
        sage: logging.disable(logging.INFO)
        sage: h_list = search_6slope_given_q(q, None, mode='heuristic', print_function=False)
    """
    h_list = []
    n_sol = 0
    if f_list is None:
        f_list = range(1, (q // 2) + 1)
    start_cpu_t = time.clock()
    for f in f_list:
        t = time.localtime()
        cpu_t = time.clock()
        print "f = %s, start_wall_time = %s:%s" % (f, t[3], t[4])
        for h in search_6slope_example(q, f, mode, print_function):
              h_list.append(h)
              n_sol += 1
              print "Found solution No.%s" % n_sol
        print "q = %s, f = %s takes cpu_time = %s" % (q, f, time.clock() - cpu_t)
    t = time.localtime()
    print "finish_wall_time = %s:%s" % (t[3], t[4])
    print "Total cpu_time = %s" %(time.clock() - start_cpu_t)
    return h_list

