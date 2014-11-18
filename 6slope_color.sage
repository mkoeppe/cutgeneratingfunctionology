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
# vertices_color (q+1)*(q+1) 0-1 array, 0: green, 1: white.
# faces_color q*q*2 0-1-2 array, 0: green, 1: currently white, 2: non-candidate, must be white.
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
    \pi(0) = 0, \pi(1) = 0 and reflection around f/2 and (1-f)/2.

    EXAMPLES::

        sage: initial_vertices_color(5, 3/5)
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
    return vertices_color

def initial_faces_color_and_covered_intervals(q, f, vertices_color):
    """
    Return initial faces_color and covered_intervals,
    corresponding to \pi(0) = 0, \pi(1) = 0 and reflection around f/2 and (1-f)/2.

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

def faces_around_vertex(q, v):
    """
    Given a grid vertex v (assume that v[0] <= v[1]), 
    return small triangle faces (only those with x <= y)
    and their vertices (only those with x <= y) that are around v.

    EXAMPLES::

        sage: faces_around_vertex(5, (1, 2))
        [((0, 2, 0), [(0, 2), (1, 2), (0, 3)]), ((0, 2, 1), [(1, 2), (0, 3), (1, 3)]), \
         ((1, 2, 0), [(1, 2), (2, 2), (1, 3)]), ((0, 1, 1), [(1, 1), (0, 2), (1, 2)]), \
         ((1, 1, 0), [(1, 1), (1, 2)]), ((1, 1, 1), [(1, 2), (2, 2)])]
        sage: faces_around_vertex(5, (1, 1))
        [((0, 1, 0), [(0, 1), (1, 1), (0, 2)]), ((0, 1, 1), [(1, 1), (0, 2), (1, 2)]), \
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

def num_slopes_at_best(q, covered_intervals):
    """
    Return an upper bound on the final number of slopes given current covered_intervals.

    EXAMPLES::

        sage: covered_intervals = [set([0,2]), set([3,4])]
        sage: num_slopes_at_best(5, covered_intervals)
        3
    """
    uncovered_num = q - sum([len(component) for component in covered_intervals])
    return uncovered_num + len(covered_intervals)

def paint_complex(q, f, last_vertices_color, last_faces_color, last_covered_intervals, candidate_faces, last_cs):
    """
    Randomly paint triangles green in a 2d-complex, until all intervals are covered.
    Return the polytope which defines the feasible region of (\pi(0), \pi(1/q),...,\pi(1)) for that paint_complex

    EXAMPLES::

        sage: q = 5; f = 3; num_of_slopes = 2;
        sage: vertices_color = initial_vertices_color(q, f);
        sage: faces_color, covered_intervals = initial_faces_color_and_covered_intervals(q, f, vertices_color)
        sage: cs = initial_cs(q, f, vertices_color)
        sage: candidate_faces = generate_candidate_faces(q, f, covered_intervals, None)
        sage: for result_polytope in paint_complex(q, f, vertices_color, faces_color, covered_intervals, candidate_faces, cs):
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
                for result_polytope in paint_complex(q, f, vertices_color, faces_color, \
                        covered_intervals, new_candidate_faces, polytope.constraints()):
                    yield result_polytope
        #faces_color[face_picked] = 2

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

def generate_vertex_function(q, polytope):
    """
    Generate real valued functions corresponding to vertices of the polytope.
    Return those functions that have required number of slope values.

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
            h = piecewise_function_from_breakpoints_and_values(bkpt, values)
            yield h

def search_6slope_example(q, f):
    """
    Search for extreme functions that have required number of slope values.

    EXAMPLES::

        sage: q=5; f=3; num_of_slopes = 2;
        sage: h = search_6slope_example(q, f).next()
    """
    #initialization
    vertices_color = initial_vertices_color(q, f)
    faces_color, covered_intervals = initial_faces_color_and_covered_intervals(q, f, vertices_color)
    cs = initial_cs(q, f, vertices_color)
    candidate_faces = generate_candidate_faces(q, f, covered_intervals, None)

    found = False
    for result_polytope in \
            paint_complex(q, f, vertices_color, faces_color, covered_intervals, candidate_faces, cs):
        for h in generate_vertex_function(q, result_polytope):
            print h
            found = True
            yield h #return h
    if not found:
        print "Example function not found. Please try again."
