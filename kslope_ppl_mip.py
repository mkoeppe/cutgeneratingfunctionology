# Make sure current directory is in path.
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

# backtracking search for 6-slope extreme functions, 
# using incremental computation for polytope, vertices_color, faces_color and covered/uncovered_intervals.

# polyhedral computation library:
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

# Main: def search_kslope(k_slopes, q, f_list, mode, prep) and def measure_stats()

from sage.libs.ppl import C_Polyhedron, Constraint, Constraint_System, Generator, Generator_System, Variable, \
                          Poly_Con_Relation, MIP_Problem, Linear_Expression
## Can't import 'point' -- clashes with plot2d point
import numpy
from sage.numerical.mip import MixedIntegerLinearProgram, MIPVariable, MIPSolverException
import sage.numerical.backends.glpk_backend as backend

poly_is_included = Poly_Con_Relation.is_included()

q_threshold = 20
dim_threshold = 11
exp_dim_prep = 9;
exp_dim_lrs = 13;
# q_threshold was 20, dim_threshold was 8, changed 3/19/2015
# use glkp mip instead of ppl polyhedron to find implied things if q > q_threshold
# do preprocessing if exp_dim >= exp_dim_prep;
# vertex-enumeration using ppl if exp_dim < exp_dim_lrs
# vertex-enumeration using lrs if exp_dim >= exp_dim_lrs

# from http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
import os, errno
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

@cached_function
def additive_constraint(q, x, y):
    z = (x + y) % q
    return Variable(x) + Variable(y) == Variable(z)

@cached_function
def delta_expr(q, x, y):
    z = (x + y) % q
    return Variable(x) + Variable(y) - Variable(z)

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
        if x <= QQ(f)/2:
            vertices_color[x, f - x] = 0
        elif (x >= f) and (x <= f - x + q):
            vertices_color[x, f - x + q] = 0
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
                if all(vertices_color[v] == 0 for v in vertices):
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
            if vertices_color[x, y] == 0:
                cs.insert(delta_expr(q, x, y) == 0) # fn[x] + fn[y] == fn[z]
            else:
                cs.insert(delta_expr(q, x, y) >= 0) # fn[x] + fn[y] >= fn[z]
    return cs

def initial_cs_reordered(q, f, vertices_color):
    """
    symmetry + subadd + nonnegativity constraints
    """
    cs = Constraint_System()
    fn = [ Variable(i) for i in range(q+1) ]

    cs.insert(fn[0] == 0)
    cs.insert(fn[q] == 0)
    cs.insert(fn[f] == 1)
    for x in range(1, q):
        for y in range(x, q):
            if vertices_color[x, y] == 0:
                cs.insert(delta_expr(q, x, y) == 0) # fn[x] + fn[y] == fn[z]
    for x in range(1, q):
        for y in range(x, q):
            if vertices_color[x, y] != 0:
                cs.insert(delta_expr(q, x, y) >= 0) # fn[x] + fn[y] >= fn[z]
    for i in range(1,q):
        cs.insert(fn[i] >= 0)
        cs.insert(fn[i] <= 1)
    return cs

def initial_cs_matrix(q, f):
    n1 = int(f / 2)
    n2 = int((f + q) / 2) - f
    cs_matrix = matrix(QQ, 2 + n1 + n2, q, sparse = True)
    # fn[0] = 0, fn[f] = 1
    cs_matrix[0, 0] = 1
    cs_matrix[1, f] = 1
    # symmetry
    for i in range(1, n1 + 1):
        cs_matrix[1 + i, i] += 1
        cs_matrix[1 + i, f - i] += 1
        cs_matrix[1 + i, f] -= 1
    for i in range(1, n2 + 1):
        cs_matrix[1 + n1 + i, f + i] += 1
        cs_matrix[1 + n1 + i, q - i] += 1
        cs_matrix[1 + n1 + i, f] -= 1
    return cs_matrix

def initial_mip(q, f, vertices_color):
    global m
    global delta
    global var_id
    var_id = numpy.zeros((q, q), int)
    m = MixedIntegerLinearProgram(maximization=True, solver = "GLPK") # "GLPK" slow! "ppl" very slow! try solver = "Gurobi"
    fn = m.new_variable(real=True, nonnegative=True)
    for i in range(q + 1):
        # 0 <= fn[i] <= 1
        m.set_min(fn[i], 0)
        m.set_max(fn[i], 1)    
    m.set_max(fn[0], 0) # fn[0] == 0
    m.set_max(fn[q], 0) # fn[q] == 0
    m.set_min(fn[f], 1) # fn[f] == 1
    num_var = q + 1
    delta = m.new_variable(real=True, nonnegative=True) # delta[x, y] >= 0
    for x in range(1, q):
        for y in range(x, q):
            z = (x + y) % q
            # delta[x, y] = fn[x] + fn[y] - fn[z]
            m.add_constraint(fn[x] + fn[y] - fn[z] - delta[x, y], min=0, max=0)
            m.set_min(delta[x, y], 0)           # fn[x] + fn[y] >= fn[z]
            if vertices_color[x, y] == 0:
                m.set_max(delta[x, y], 0)       # fn[x] + fn[y] == fn[z]
            else:
                m.set_max(delta[x, y], None)
            var_id[x, y] = num_var
            num_var += 1
    m.set_objective(None)
    m.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
    m.solver_parameter("obj_upper_limit", 0.01)

@cached_function
def edges_around_vertex(q, v):
    """
    Given a grid vertex v (assume that v[0] <= v[1], v is not on the border),
    return a list of elements corresponding to the edges connected to v.
    Each element has the form (set([a, b]), v'). 
    v' is a grid vertex next to v. 
    If the edge (v, v') is green, then segments a and b are connected.

    EXAMPLES::

        sage: edges_around_vertex(5, (1, 3))
        [(set([0, 3]), (0, 3)), (set([3, 4]), (1, 4)), (set([1, 4]), (2, 3)), (set([2, 3]), (1, 2)), (set([0, 3]), (0, 4)), (set([1, 2]), (2, 2))]
        sage: edges_around_vertex(5, (1, 2))
        [(set([0, 2]), (0, 2)), (set([2, 3]), (1, 3)), (set([1, 3]), (2, 2)), (set([1, 2]), (1, 1)), (set([0, 2]), (0, 3))]
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

@cached_function
def faces_around_vertex(q, v):
    """
    Given a grid vertex v (assume that v[0] <= v[1], v is not on the border),
    return small triangle faces (only those with x <= y)
    and their vertices (only those with x <= y) that are around v.

    EXAMPLES::

        sage: faces_around_vertex(5, (1, 2))
        [((0, 2, 0), [(0, 2), (1, 2), (0, 3)]), ((0, 2, 1), [(1, 2), (0, 3), (1, 3)]), ((1, 2, 0), [(1, 2), (2, 2), (1, 3)]), ((0, 1, 1), [(1, 1), (0, 2), (1, 2)]), ((1, 1, 0), [(1, 1), (1, 2)]), ((1, 1, 1), [(1, 2), (2, 2)])]
        sage: faces_around_vertex(5, (1, 1))
        [((0, 1, 0), [(0, 1), (1, 1), (0, 2)]), ((0, 1, 1), [(1, 1), (0, 2), (1, 2)]), ((1, 1, 0), [(1, 1), (1, 2)]), ((0, 0, 1), [(0, 1), (1, 1)])]
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
    Compute incrementally new covered_intervals and new uncovered_intervals, resulting from
    adding a new green edge that connects the elements in to_merge_set.

    EXAMPLES::

        sage: last_covered_intervals = [set([0,2]), set([3,4])]
        sage: last_uncovered_intervals = [set([1])]
        sage: to_merge_set = set([0, 1])
        sage: update_covered_uncovered_by_adding_edge(last_covered_intervals, last_uncovered_intervals, to_merge_set, 5)
        ([set([3, 4]), set([0, 1, 2])], [])
    """
    covered_intervals = []
    uncovered_intervals = []
    new_set = copy(to_merge_set)
    new_set_is_covered = False
    for component in last_uncovered_intervals:
        if component & new_set:
            new_set.update(component)
        else:
            uncovered_intervals.append(component)
    for component in last_covered_intervals:
        if component & new_set:
            new_set.update(component)
            new_set_is_covered = True
        else:
            covered_intervals.append(component)
    if new_set_is_covered:
        covered_intervals.append(new_set)
    else:
        uncovered_intervals.append(new_set)
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

def generate_uncovered_set(q, uncovered_intervals):
    uncovered_set = set([])
    for component in uncovered_intervals:
        uncovered_set.update(component)
    return uncovered_set

def generate_candidate_faces(q, f, covered_intervals, last_face=None, faces_color=None, sym=False):
    """
    Return a list of candidate_faces (lexicographically > last_face)
    to paint in next step, whose I, J are currently uncovered.
    Note that candidate_faces only takes faces with I <= J.

    EXAMPLES::

        sage: q = 5; f = 3;
        sage: vertices_color = initial_vertices_color(q, f);
        sage: faces_color, covered_intervals = initial_faces_color_and_covered_intervals(q, f, vertices_color)
        sage: generate_candidate_faces(q, f, covered_intervals)
        [(1, 1, 0), (1, 1, 1)]
    """
    to_cover = generate_to_cover(q, covered_intervals)
    # NOTE: candidate_faces only takes faces with x <= y
    candidate_faces = []
    for x in to_cover:
        for y in to_cover:
            if (x <= y) and (not sym or x + y < q - 1) :
                face = (x, y, 0)
                if (last_face is None or face > last_face) and (faces_color is None or faces_color[face] == 1):
                    candidate_faces.append(face)
                face = (x, y, 1)
                if (last_face is None or face > last_face) and (faces_color is None or faces_color[face] == 1):
                    candidate_faces.append(face)
    return candidate_faces

def num_slopes_at_best(q, f, covered_intervals, uncovered_intervals=None):
    """
    Return an upper bound on the final number of slopes,
    given current covered_intervals (and optionally, uncovered_intervals 
    that provides connected components of non-covered intervals).

    EXAMPLES::

        sage: covered_intervals = [set([0,2])]
        sage: uncovered_intervals = [set([1]), set([3,4])]
        sage: num_slopes_at_best(5, 3, covered_intervals)
        3
        sage: num_slopes_at_best(5, 3, covered_intervals, uncovered_intervals)
        3
    """
    if uncovered_intervals is None:
        #uncovered_num = q - sum([len(component) for component in covered_intervals])
        # consider connected components in uncovered_intervals
        to_cover = set(range(0, (f + 1) / 2) + range(f, (f + q + 1) / 2))
        for component in covered_intervals:
            to_cover -= component
        uncovered_num = len(to_cover)
    else:

        uncovered_num = len(uncovered_intervals)
    return uncovered_num + len(covered_intervals)

def update_around_green_face(q, f, vertices_color, faces_color, covered_intervals, (x, y, w)):
    """
    Subfunction of paint_complex_heuristic() and paint_complex_fulldim_covers().
    Painting triangle (x, y, w) from white to green induces some new green triangles around it.
    Update vertices_color, faces_color, correspondingly.
    If there is non_candidate among implied green faces, return (False, None, changed_vertices, changed_faces).
    Otherwise, return (True, updated covered_intervals, changed_vertices, changed_faces).
    """
    changed_vertices = []
    changed_faces = []
    if x < y:
        vertices_picked = [(x+1, y), (x, y+1), (x+w, y+w)]
    else:
        vertices_picked = [(x, y+1), (x+w, y+w)]
    for (i, j) in vertices_picked:
        if (vertices_color[i, j] == 1):
            # new green vertice
            vertices_color[i, j] = 0
            changed_vertices.append((i, j))
            for (face, vertices) in faces_around_vertex(q, (i, j)):
                if faces_color[face] != 0 and all(vertices_color[v] == 0 for v in vertices):
                    # find new green face.
                    if faces_color[face] == 2: # face is in non_candidate
                        return False, None, changed_vertices, changed_faces
                    faces_color[face] = 0
                    changed_faces.append(face)
                    covered_intervals = directly_covered_by_adding_face(covered_intervals, face, q, f)
    return True, covered_intervals, changed_vertices, changed_faces

def update_implied_faces_pol(q, f, vertices_color, changed_vertices, faces_color, changed_faces, covered_intervals, polytope):
    """
    Subfunction of paint_complex_heuristic() and paint_complex_fulldim_covers().
    Look for implied additive vertices and faces, given polytope.
    Update vertices_color, changed_vertices, faces_color, changed_faces and polytope
    If there is non_candidate among implied green faces, return False.
    Otherwise, return True and updated covered_intervals.
    """
    for i in range(1, q):
        for j in range(i, q):
            if (vertices_color[i, j] == 1) and \
                    polytope.relation_with( additive_constraint(q, i, j) ).implies(poly_is_included):
                # find implied additive vertices
                vertices_color[i, j] = 0
                changed_vertices.append((i, j))
                for (face, vertices) in faces_around_vertex(q, (i, j)):
                    if faces_color[face] != 0 and all(vertices_color[v] == 0 for v in vertices):
                        # find new green face.
                        if faces_color[face] == 2: # face is in non_candidate
                            return False, None
                        faces_color[face] = 0
                        changed_faces.append(face)
                        covered_intervals = directly_covered_by_adding_face(covered_intervals, face, q, f)
    return True, covered_intervals

def update_implied_faces_mip(q, f, vertices_color, changed_vertices, faces_color, changed_faces, covered_intervals, sym=False):
    """
    Subfunction of paint_complex_combined_mip().
    Look for implied additive vertices and faces, given MILP m.
    Update vertices_color, changed_vertices, faces_color, changed_faces
    If there is non_candidate among implied green faces, return False.
    Otherwise, return True and updated covered_intervals.
    """
    for i in range(1, q):
        for j in range(i, q):
            if (vertices_color[i, j] == 1) and (not sym or i + j <= q):
                #m.set_objective(delta[i, j])
                m_glpk = m.get_backend()
                m_glpk.objective_coefficient(var_id[i, j], 1)
                obj_val = m.solve(objective_only=True)
                m_glpk.objective_coefficient(var_id[i, j], 0)
                if  obj_val == 0:
                    # find implied additive vertices
                    vertices_color[i, j] = 0
                    changed_vertices.append((i, j))
                    for (face, vertices) in faces_around_vertex(q, (i, j)):
                        if faces_color[face] != 0 and all(vertices_color[v] == 0 for v in vertices):
                            # find new green face.
                            if faces_color[face] == 2: # face is in non_candidate
                                return False, None
                            faces_color[face] = 0
                            changed_faces.append(face)
                            covered_intervals = directly_covered_by_adding_face(covered_intervals, face, q, f)
    return True, covered_intervals


def paint_complex_heuristic(k_slopes, q, f, vertices_color, faces_color, last_covered_intervals, candidate_faces, cs):
    """
    Paint triangles green in a 2d-complex, until all intervals are covered.
    Return the polytope which defines the feasible region of (\pi(0), \pi(1/q),...,\pi(1)) for that paint_complex_heuristic

    Heuristic: 
        I, J projections of the next face to paint are chosen from currently uncovered intervals.
        Stop painting immediately once everything is covered.
        Edges are not considered.

    EXAMPLES::

        sage: q = 5; f = 3; k_slopes = 2;
        sage: vertices_color = initial_vertices_color(q, f);
        sage: faces_color, covered_intervals = initial_faces_color_and_covered_intervals(q, f, vertices_color)
        sage: cs = initial_cs(q, f, vertices_color)
        sage: candidate_faces = generate_candidate_faces(q, f, covered_intervals)
        sage: for result_polytope in paint_complex_heuristic(k_slopes, q, f, vertices_color, faces_color, covered_intervals, candidate_faces, cs):
        ...       result_polytope.minimized_generators()
        Generator_System {point(0/6, 2/6, 4/6, 6/6, 3/6, 0/6)}
        Generator_System {point(0/4, 3/4, 1/4, 4/4, 2/4, 0/4)}
    """
    for (x, y, w) in candidate_faces:
        faces_color[(x, y, w)] = 0
        covered_intervals = directly_covered_by_adding_face(last_covered_intervals, (x, y, w), q, f)
        legal_picked, covered_intervals, changed_vertices, changed_faces = update_around_green_face( \
                    q, f, vertices_color, faces_color, covered_intervals, (x, y, w))
        # If encounter non_candidate or too few slopes, stop recursion
        if legal_picked and num_slopes_at_best(q, f, covered_intervals) >= k_slopes:
            polytope = C_Polyhedron(cs)
            # update polytope
            for (i, j) in changed_vertices:
                polytope.add_constraint( additive_constraint(q, i, j) )
            if not polytope.is_empty():
                # If infeasible, stop recursion
                # look for implied additive vertices and faces
                legal_picked, covered_intervals = update_implied_faces_pol(q, f, \
                        vertices_color, changed_vertices, faces_color, changed_faces, covered_intervals, polytope)
                # If encounter non_candidate or too few slopes, stop recursion.
                if legal_picked and num_slopes_at_best(q, f, covered_intervals) >= k_slopes:
                    new_candidate_faces= generate_candidate_faces(q, f, covered_intervals, last_face=(x, y, w))
                    if not new_candidate_faces:
                        # stop recursion
                        if not generate_to_cover(q, covered_intervals):
                            # all covered, finish
                            yield polytope
                    else:
                        for result_polytope in paint_complex_heuristic(k_slopes, q, f, vertices_color, faces_color, \
                                covered_intervals, new_candidate_faces, polytope.constraints()):
                            # Note: use minimized_constraints() in 'heuristic' mode takes longer. WHY??
                            yield result_polytope
        # Now, try out white triangle (x, y, w)
        for face in changed_faces:
            faces_color[face] = 1
        for v in changed_vertices:
            vertices_color[v] = 1
        faces_color[(x, y, w)] = 2
    for face in candidate_faces:
        faces_color[face] = 1

def paint_complex_fulldim_covers(k_slopes, q, f, vertices_color, faces_color, last_covered_intervals, (x, y, w), cs):
    """
    Paint triangles green in a 2d-complex, until all possibilities are tried.
    If all intervals are covered, 
    return the polytope which defines the feasible region of (\pi(0), \pi(1/q),...,\pi(1)) for that paint_complex_fulldim_covers

    fulldim_covers: 
        Allow painting any face, even if covered already.
        Edges are not considered.

    EXAMPLES::

        sage: q = 5; f = 3; k_slopes = 2;
        sage: vertices_color = initial_vertices_color(q, f);
        sage: faces_color, covered_intervals = initial_faces_color_and_covered_intervals(q, f, vertices_color)
        sage: cs = initial_cs(q, f, vertices_color)
        sage: for result_polytope in paint_complex_fulldim_covers(k_slopes, q, f, vertices_color, faces_color, covered_intervals, (0, 0, 0), cs):
        ...       result_polytope.minimized_generators()
        Generator_System {point(0/6, 2/6, 4/6, 6/6, 3/6, 0/6)}
        Generator_System {point(0/4, 3/4, 1/4, 4/4, 2/4, 0/4)}
    """
    picked_faces = []
    while (x, y, w) < (q - 1, q - 1, 1):
        # move to the next triangle
        w += 1
        if w == 2:
            w = 0
            y += 1
        if y == q:
            x += 1
            y = x
        if faces_color[(x, y, w)] == 1: # color is unkown
            picked_faces.append( (x, y, w) )
            # First, try out green triangle (x, y, w)
            faces_color [(x, y, w)] = 0
            covered_intervals = directly_covered_by_adding_face(last_covered_intervals, (x, y, w), q, f)
            legal_picked, covered_intervals, changed_vertices, changed_faces = update_around_green_face( \
                    q, f, vertices_color, faces_color, covered_intervals, (x, y, w))
            # If encounter non_candidate or too few slopes, stop recursion
            if legal_picked and num_slopes_at_best(q, f, covered_intervals) >= k_slopes:
                polytope = C_Polyhedron(cs)
                # update polytope
                for (i, j) in changed_vertices:
                    polytope.add_constraint( additive_constraint(q, i, j) )
                if not polytope.is_empty():
                    # If infeasible, stop recursion
                    # look for implied additive vertices and faces
                    legal_picked, covered_intervals = update_implied_faces_pol(q, f, \
                            vertices_color, changed_vertices, faces_color, changed_faces, covered_intervals, polytope)
                    # If encounter non_candidate or too few slopes, stop recursion
                    if legal_picked and num_slopes_at_best(q, f, covered_intervals) >= k_slopes:
                        for result_polytope in paint_complex_fulldim_covers(k_slopes, q, f, \
                                vertices_color, faces_color, covered_intervals, (x, y, w), polytope.minimized_constraints()):
                            yield result_polytope
            # Now, try out white triangle (x, y, w)
            for face in changed_faces:
                faces_color[face] = 1
            for v in changed_vertices:
                vertices_color[v] = 1
            faces_color [(x, y, w)] = 2
    if (x, y, w) == (q - 1, q - 1, 1):
        # finish painting, check if all intervals are covered
        if not generate_to_cover(q, last_covered_intervals):
            # all covered, return valid painting
            polytope = C_Polyhedron(cs)
            if not polytope.is_empty():
                yield polytope
    #recover picked_faces
    for face in picked_faces:
        faces_color[face] = 1

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
    uncovered_intervals = [set([i, f - i - 1]) for i in range(1, (f + 1) / 2)]
    uncovered_intervals += [set([i, q + f - i - 1]) for i in range(f + 1, (f + q + 1) / 2)]
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
    was_v = set([])
    for (face, vertices) in faces_around_vertex(q, (x, y)):
        if all(vertices_color[v] == 0 for v in vertices):
            # find new green face.
            was_v.update(vertices)
            covered_intervals, uncovered_intervals = \
                update_covered_uncovered_by_adding_face(covered_intervals, uncovered_intervals, face, q)
    for (to_merge_set, v) in edges_around_vertex(q, (x, y)):
        # if an green edge is included in a green face, don't need to update_covered_uncovered again.
        if vertices_color[v] == 0 and not v in was_v:
            # find new green edge.
            covered_intervals, uncovered_intervals = \
                update_covered_uncovered_by_adding_edge(covered_intervals, uncovered_intervals, to_merge_set, q)
    return covered_intervals, uncovered_intervals

def update_implied_vertices(k_slopes, q, f, vertices_color, changed_vertices, covered_intervals, uncovered_intervals, polytope):
    """
    Subfunction of paint_complex_complete().
    Look for implied additive vertices, given polytope.
    Update vertices_color, changed_vertices and polytope
    If there is non_candidate among implied green vertices, return False.
    Otherwise, return True, updated covered_intervals and uncovered_intervals.
    """
    for i in range(1, q):
        for j in range(i, q):
            if (vertices_color[i, j] != 0) and \
                    polytope.relation_with( additive_constraint(q, i, j) ).implies(poly_is_included):
                # find implied additive vertices
                if vertices_color[i, j] == 2:
                    # encounter non_candidate vertex, stop recursion
                    return False, None, None
                # paint this implied vertex green
                vertices_color[i, j] = 0
                changed_vertices.append((i, j))
                covered_intervals, uncovered_intervals = update_around_green_vertex(q, (i, j), \
                                        vertices_color, covered_intervals, uncovered_intervals)
                if num_slopes_at_best(q, f, covered_intervals, uncovered_intervals) < k_slopes:
                    return False, None, None
                #else: # This is not necessary. add_constraint only taks longer
                #    polytope.add_constraint( additive_constraint(q, i, j) )
    return True, covered_intervals, uncovered_intervals

def paint_complex_complete(k_slopes, q, f, vertices_color, last_covered_intervals, last_uncovered_intervals, (x, y), cs):
    """
    Paint vertices green in a 2d-complex, until all possibilities are tried.
    If all intervals are covered (consider both green triangles and edges),
    return the polytope which defines the feasible region of (\pi(0), \pi(1/q),...,\pi(1)) for that paint_complex_complete

    EXAMPLES::
        sage: q = 5; f = 3; k_slopes = 2;
        sage: vertices_color = initial_vertices_color(q, f);
        sage: cs = initial_cs(q, f, vertices_color)
        sage: covered_intervals, uncovered_intervals = initial_covered_uncovered(q, f, vertices_color)
        sage: for result_polytope in paint_complex_complete(k_slopes, q, f, vertices_color, covered_intervals, uncovered_intervals, (0, 0), cs):
        ...       result_polytope.minimized_generators()
        Generator_System {point(0/6, 2/6, 4/6, 6/6, 3/6, 0/6)}
        Generator_System {point(0/4, 3/4, 1/4, 4/4, 2/4, 0/4)}
    """
    picked_vertices = []
    while (x, y) < (q - 1, q - 1):
        # move to the next vertices
        y += 1
        if y == q:
            x += 1
            y = x
        if vertices_color[(x, y)] == 1: # color is unkown
            picked_vertices.append((x, y))
            changed_vertices = []
            # First, try out green vertex (x, y)
            vertices_color[(x, y)] = 0
            covered_intervals, uncovered_intervals = update_around_green_vertex(q, (x, y), \
                            vertices_color, last_covered_intervals, last_uncovered_intervals)
            # If too few slopes, stop recursion
            if num_slopes_at_best(q, f, covered_intervals, uncovered_intervals) >= k_slopes:
                polytope = C_Polyhedron(cs)
                polytope.add_constraint( additive_constraint(q, x, y) )
                # If infeasible, stop recursion
                if not polytope.is_empty():
                    # look for implied additive vertices
                    legal_picked, covered_intervals, uncovered_intervals = update_implied_vertices(k_slopes, q, f, \
                            vertices_color, changed_vertices, covered_intervals, uncovered_intervals, polytope)
                    if legal_picked:
                        for result_polytope in paint_complex_complete(k_slopes, q, f, vertices_color, \
                                covered_intervals, uncovered_intervals, (x, y), polytope.minimized_constraints()):
                            yield result_polytope
            # Now, try out white vertex (x, y)
            for (i, j) in changed_vertices:
                vertices_color[i, j] = 1
            vertices_color[(x, y)] = 2
    if (x, y) == (q - 1, q - 1):
        # finish painting, check if all intervals are covered
        if not last_uncovered_intervals:
            # all covered, return valid painting
            polytope = C_Polyhedron(cs)
            if not polytope.is_empty():
                yield polytope
    # recover picked_vertices
    for v in picked_vertices:
        vertices_color[v] = 1

def vertex_enumeration(polytope, prep=True, exp_dim=-1, vetime=False):
    if vetime:
        st = os.times();
    if not prep or (0 <= exp_dim < exp_dim_prep):
        # no preprocessing
        extreme_points = polytope.minimized_generators()
    elif exp_dim >= exp_dim_lrs:
        # preprocessing and vertex enumertation using redund + lrs
        cs = polytope.constraints()
        cs_prep_lrs_str = remove_redundancy_from_cs(cs, return_lrs=True)
        extreme_points = lrs_lrsinput_pploutput(cs_prep_lrs_str)
    else:
        # preprocessing and vertex enumertation using redund + ppl
        cs = polytope.constraints()
        cs_prep = remove_redundancy_from_cs(cs)
        polytope = C_Polyhedron(cs_prep)
        extreme_points = polytope.minimized_generators()
    if vetime:
        et = os.times(); 
        logging.info("user=%s, sys=%s, child user=%s, child sys=%s" %(et[0]-st[0], et[1]-st[1], et[2]-st[2], et[3]-st[3]))
        t = sum([et[i]-st[i] for i in range(4)]);
        logging.info("Vertex enumeration time = %s" % t)
    return extreme_points

def generate_vertex_values(k_slopes , q, polytope,  v_set=set([]), prep=True, exp_dim=-1, vetime=False):
    """
    Return vertices of the polytope, whose corresponding 
    piecewise_linear function h has at least  k_slopes.
    """
    if exp_dim == -1:
        exp_dim = int(q / 2) - 1
    extreme_points = vertex_enumeration(polytope, prep=prep, exp_dim=exp_dim, vetime=vetime)
    for v in extreme_points:
        v_n = v.coefficients()
        num = len(set([v_n[i+1] - v_n[i] for i in range(q)]))
        if num >= k_slopes:
            if not tuple(v_n) in v_set:
                v_set.add(tuple(v_n))
                yield v_n

def h_from_vertex_values(v_n):
    n = len(v_n)
    bkpt = [QQ(x) / (n-1) for x in range(n)]
    v_d = max(v_n)
    values = [QQ(y) / v_d for y in v_n]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def search_kslope_example(k_slopes, q, f, mode='combined', prep=True):
    """
    Search for extreme functions that have required number of slope values.

    If `mode` is 'heuristic', use paint_complex_heuristic() to paint;
    If `mode` is (defaut)'combined', use paint_complex_combined() to paint a few triangles, then enumerate vertex-functions;
    If `mode` is 'fulldim_covers', use paint_complex_fulldim_covers() to paint;
    If `mode` is 'complete', use paint_complex_complete() to paint;
    If `mode` is 'naive', enumerate vertex-fuctions and check whether all invervals are covered;

    EXAMPLES::

        sage: logging.disable(logging.INFO)
        sage: q=5; f=3; k_slopes = 2;
        sage: h = search_kslope_example(k_slopes, q, f, mode='heuristic').next()
        sage: h == piecewise_function_from_breakpoints_and_values([0, 3/5, 1], [0, 1, 0])
        True
    """
    #initialization
    vertices_color = initial_vertices_color(q, f)
    if mode == 'heuristic':
        cs = initial_cs(q, f, vertices_color)
        faces_color, covered_intervals = initial_faces_color_and_covered_intervals(q, f, vertices_color)
        candidate_faces = generate_candidate_faces(q, f, covered_intervals, last_face=None)
        gen = paint_complex_heuristic(k_slopes, q, f, vertices_color, faces_color, covered_intervals, candidate_faces, cs)
    elif mode == 'combined':      
        faces_color, covered_intervals = initial_faces_color_and_covered_intervals(q, f, vertices_color)
        candidate_faces = generate_candidate_faces(q, f, covered_intervals, last_face=None)
        cs_matrix = initial_cs_matrix(q, f)
        if q > q_threshold:
            initial_mip(q, f, vertices_color) # set initial global variable m and delta
            gen = paint_complex_combined_mip(k_slopes, q, f, vertices_color, faces_color, covered_intervals, candidate_faces, cs_matrix)
        else:
            cs = initial_cs(q, f, vertices_color)
            gen = paint_complex_combined_pol(k_slopes, q, f, vertices_color, faces_color, covered_intervals, candidate_faces, cs, cs_matrix)
    elif mode == 'sym':
        # note that q = 2*f
        vertices_color, faces_color, covered_intervals, candidate_faces, cs_matrix = initialization_sym(q, f)
        if not candidate_faces:
            # imposing too much initial green, candidate_faces is empty, but still have uncovered.
            raise ValueError, "imposing too much initial green, candidate_faces is empty, but still have uncovered."
            gen = gen_initial_polytope_sym(q, f, vertices_color, covered_intervals)
        else:
            gen = paint_complex_combined_mip(k_slopes, q, f, vertices_color, faces_color, covered_intervals, candidate_faces, cs_matrix, sym=True)
    elif mode == 'no_implied': # useless
        faces_color, covered_intervals = initial_faces_color_and_covered_intervals(q, f, vertices_color)
        candidate_faces = generate_candidate_faces(q, f, covered_intervals, last_face=None)
        cs_matrix = initial_cs_matrix(q, f)
        gen = paint_complex_no_implied(k_slopes, q, f, vertices_color, faces_color, covered_intervals, candidate_faces, cs_matrix)
    elif mode == 'fulldim_covers': #useless
        cs = initial_cs(q, f, vertices_color)
        faces_color, covered_intervals = initial_faces_color_and_covered_intervals(q, f, vertices_color)
        gen = paint_complex_fulldim_covers(k_slopes, q, f, vertices_color, faces_color, covered_intervals, (0, 0, 0), cs)
    elif mode == 'complete': #useless
        cs = initial_cs(q, f, vertices_color)
        covered_intervals, uncovered_intervals = initial_covered_uncovered(q, f, vertices_color)
        gen = paint_complex_complete(k_slopes, q, f, vertices_color, covered_intervals, uncovered_intervals, (1, 0), cs)
    elif mode == 'naive':
        cs = initial_cs(q, f, vertices_color)
        polytope = C_Polyhedron(cs)
    elif mode == 'sym_naive': #useless
        # note that q = 2*f
        #f2 = int(f / 2)
        #vertices_color[1, f2] = vertices_color[q - f2, q - 1] = 0 #???
        #vertices_color[f2, f2] = vertices_color[q - f2, q - f2] = 0 #???
        cs = initial_cs_sym(q, f, vertices_color)
        polytope = C_Polyhedron(cs)
    else:
        raise ValueError, "mode must be one of {'heuristic', 'combined', 'fulldim_covers', 'complete', 'naive', 'no_implied', 'sym', 'sym_naive'}."
    v_set = set([])
    if mode == 'combined' or mode == 'no_implied' or mode == 'sym':
        #global filename
        #filename = open(dir_math + "profiler/dim_threshold/%sslope_q%s_f%s_%sdim.txt" % (k_slopes, q, f, dim_threshold), "w")
        #print >> filename, "k_slope = %s, q = %s, f = %s, dim_threshold = %s" % (k_slopes, q, f, dim_threshold)
        for result_polytope, result_covered_intervals, exp_dim in gen:
            for values in generate_vertex_values(k_slopes, q, result_polytope, v_set, prep=prep, exp_dim=exp_dim):
                if all_intervals_covered(q, f, values, result_covered_intervals):
                    yield values
        #filename.close()
    elif mode == 'naive' or mode == 'sym_naive':
        for values in generate_vertex_values(k_slopes, q, polytope, v_set, prep=prep, vetime=True):
            if all_intervals_covered(q, f, values, []):
                yield values
    else:
        for result_polytope in gen:
            for values in generate_vertex_values(k_slopes, q, result_polytope, v_set, prep=prep):
                yield values
    if not v_set:
        logging.info("Example function not found. Please try again.")

import time
def search_kslope(k_slopes, q, f_list=None, mode='naive', prep=True, print_function=False):
    """
    EXAMPLES::

        sage: logging.disable(logging.INFO)
        sage: q = 12; k_slopes = 3;
        sage: h_list = search_kslope(k_slopes, q, None, mode='heuristic', print_function=False)
        sage: len(h_list)
        13
    """
    global h_list
    h_list = []
    n_sol = 0
    if f_list is None:
        f_list = range(1, (q / 2) + 1)
    time_ini = os.times()
    for f in f_list:
        time_start = os.times()
        logging.info( "Search for extreme funtions with q = %s, f = %s, k_slopes >= %s, mode = %s" % (q, f, k_slopes, mode) )
        for h in search_kslope_example(k_slopes, q, f, mode=mode, prep=prep):
            h_list.append(h)
            n_sol += 1
            if print_function:
                logging.info("h = %s" % (h,))
            time_end = os.times()
            t = sum([time_end[i]-time_start[i] for i in range(4)])
            logging.info("Found solution No.%s, in %s seconds" % (n_sol, t))
        time_end = os.times()
        t = sum([time_end[i]-time_start[i] for i in range(4)])
        logging.info("q = %s, f = %s takes cpu_time = %s" % (q, f, t))
    time_end = os.times()
    t = sum([time_end[i]-time_ini[i] for i in range(4)])
    logging.info("Total cpu_time = %s" % t)
    return h_list

def measure_stats(q, f_list, name=None, reordered_cs=False, prep=True):
    """
    Given (q, f), write {num_vertices, running times, 
    percentage and number of extreme functions, 
    percentage and number of fulldim_cover extreme functions,
    and number of k-slope extreme functions} to a file
    """
    if name is None:
        fout = sys.stdout
    else:
        fout = open("%s.txt" % name, "a")
    for f in f_list:
        cpu_t = time.clock()
        vertices_color = initial_vertices_color(q, f)
        if reordered_cs:
            cs = initial_cs_reordered(q, f, vertices_color)
        else:
            cs = initial_cs(q, f, vertices_color)
        if prep:
            cs = remove_redundancy_from_cs(cs)
        extreme_points = C_Polyhedron(cs).minimized_generators()
        tot_vertices = len(extreme_points)
        print >> fout, "q = %s, f = %s, num_vertices = %s, poly_gen_time = %s," % (q, f, tot_vertices, time.clock() - cpu_t)
        num_extreme = [0]*11
        num_full = [0]*11
        for v in extreme_points:
            v_n = v.coefficients()
            slopes = len(set([v_n[i+1] - v_n[i] for i in range(q)]))
            if slopes == 2:
                num_extreme[2] += 1
                num_full[2] += 1
            else:
                whether_covered = all_intervals_covered(q, f, v_n, [])
                if whether_covered:
                    num_extreme[slopes] += 1
                    if whether_covered == 'direct':
                        num_full[slopes] += 1
        tot_extreme = sum(num_extreme)
        tot_full = sum(num_full)
        print >> fout, "        perc_extreme = %s, perc_full = %s, cpu_time = %s," %(tot_extreme*100/tot_vertices, \
                                                                tot_full*100/tot_vertices, time.clock() - cpu_t)
        print >> fout, "        num_extreme = %s,  k-slope  %s" % (tot_extreme, num_extreme[2:])
        print >> fout, "        num_full = %s,  k-slope  %s" % (tot_full, num_full[2:])
    return

def dim_cs_matrix(q, changed_vertices, cs_matrix):
    """
    construct the new cs_matrix, and compute its rank
    """
    new_cs_matrix = cs_matrix
    for (i, j) in changed_vertices:
        new_cs_row = zero_vector(QQ, q)
        new_cs_row[i] += 1
        new_cs_row[j] += 1
        new_cs_row[(i + j) % q] -= 1
        new_cs_matrix = new_cs_matrix.stack(new_cs_row)
    d = q - new_cs_matrix.rank()
    return d, new_cs_matrix

def paint_complex_combined_pol(k_slopes, q, f, vertices_color, faces_color, last_covered_intervals, candidate_faces, cs, cs_matrix):
    """
    Combine 'heuristic' backracting search with vertex enumeration.
    Stop backtracting when q - rank(cs_matrix) <= dim_threshold
    """
    for (x, y, w) in candidate_faces:
        covered_intervals = directly_covered_by_adding_face(last_covered_intervals, (x, y, w), q, f)
        legal_picked, covered_intervals, changed_vertices, changed_faces = update_around_green_face( \
                    q, f, vertices_color, faces_color, covered_intervals, (x, y, w))
        # If encounter non_candidate or too few slopes, stop recursion
        if legal_picked and num_slopes_at_best(q, f, covered_intervals) >= k_slopes:
            polytope = C_Polyhedron(cs)
            # update polytope
            for (i, j) in changed_vertices:
                polytope.add_constraint( additive_constraint(q, i, j) )
            if not polytope.is_empty():
                # If infeasible, stop recursion
                the_last_changed = len(changed_vertices)
                # look for implied additive vertices and faces
                legal_picked, covered_intervals = update_implied_faces_pol(q, f, \
                        vertices_color, changed_vertices, faces_color, changed_faces, covered_intervals, polytope)
                # If encounter non_candidate or too few slopes, stop recursion.
                if legal_picked and num_slopes_at_best(q, f, covered_intervals) >= k_slopes:
                    new_candidate_faces= generate_candidate_faces(q, f, covered_intervals, last_face=(x, y, w))
                    # use changed_vertices[0:the_last_changed]. don't put implied equalities in cs_matrix.
                    exp_dim, new_cs_matrix = dim_cs_matrix(q, changed_vertices[0:the_last_changed], cs_matrix)
                    #exp_dim, new_cs_matrix = dim_cs_matrix(q, changed_vertices, cs_matrix)
                    if not new_candidate_faces or exp_dim <= dim_threshold:
                        # Suppose that k_slopes > 2. If k_slopes = 2, waist time on checking covered for 2-slope functions.
                        # stop recursion
                        yield polytope,  covered_intervals, exp_dim
                    else:
                        for result_polytope, result_covered_intervals, result_exp_dim in paint_complex_combined_pol(k_slopes, q, f, \
                                vertices_color, faces_color, covered_intervals, new_candidate_faces, polytope.constraints(), new_cs_matrix):
                            # Note: use minimized_constraints() in 'heuristic' mode takes longer. WHY??
                            yield result_polytope, result_covered_intervals, result_exp_dim
        # Now, try out white triangle (x, y, w)
        for face in changed_faces:
            faces_color[face] = 1
        for v in changed_vertices:
            vertices_color[v] = 1
        faces_color[(x, y, w)] = 2
    for face in candidate_faces:
        faces_color[face] = 1

def paint_complex_combined_mip(k_slopes, q, f, vertices_color, faces_color, last_covered_intervals, candidate_faces, cs_matrix, sym=False):
    """
    Combine 'heuristic' backracting search (using MILP library) with vertex enumeration.
    For q <= q_threshold, Polyhedron is faster. 
    If q - rank(cs_matrix) <= dim_threshold, stop backtracking search.
    Enumerate and check vertex functions then.
    if sym is True,
    restrict to the special 2d-diagram such as kzh_6_slope_fulldim_covers_2() and kzh_6_slope_fulldim_covers_3().
    q = 2f, f mod 2 = 1. symmetry between left and right parts of f. 2d-diagram is white in the middle.
    """
    the_last_constraint = m.number_of_constraints()
    for (x, y, w) in candidate_faces:
        covered_intervals = directly_covered_by_adding_face(last_covered_intervals, (x, y, w), q, f)
        legal_picked, covered_intervals, changed_vertices, changed_faces = update_around_green_face( \
                    q, f, vertices_color, faces_color, covered_intervals, (x, y, w))
        if sym and legal_picked:
            for (i, j) in changed_vertices:
                vertices_color[q - j, q - i] = vertices_color[i, j]
            for (x1, y1, w1) in changed_faces:
                (x2, y2, w2) = (q - y1 - 1, q - x1 - 1, 1 - w1)
                faces_color[(x2, y2, w2)] = faces_color[(x1, y1, w1)]
                covered_intervals = directly_covered_by_adding_face(covered_intervals, (x2, y2, w2), q, f)
        # If encounter non_candidate or too few slopes, stop recursion
        if legal_picked and num_slopes_at_best(q, f, covered_intervals) >= k_slopes:
            # update MILP m
            for (i, j) in changed_vertices:
                m.set_max(delta[i, j], 0)
            #m.set_objective(None)
            m.solver_parameter("primal_v_dual", "GLP_DUAL")
            try:
                m.solve(objective_only=True)
                is_feasible = True
            except MIPSolverException:
                # If infeasible, stop recursion
                is_feasible = False
            if is_feasible:
                # look for implied additive vertices and faces
                # add implied equalities to m
                the_last_changed_v = len(changed_vertices)
                the_last_changed_f = len(changed_faces)
                m.solver_parameter("primal_v_dual", "GLP_PRIMAL")
                legal_picked, covered_intervals = update_implied_faces_mip(q, f, \
                        vertices_color, changed_vertices, faces_color, changed_faces, covered_intervals, sym=sym)
                if sym and legal_picked:
                    for (i, j) in changed_vertices[the_last_changed_v::]:
                        vertices_color[q - j, q - i] = vertices_color[i, j]
                    for (x1, y1, w1) in changed_faces[the_last_changed_f::]:
                        (x2, y2, w2) = (q - y1 - 1, q - x1 - 1, 1 - w1)
                        faces_color[(x2, y2, w2)] = faces_color[(x1, y1, w1)]
                        covered_intervals = directly_covered_by_adding_face(covered_intervals, (x2, y2, w2), q, f)
                # If encounter non_candidate or too few slopes, stop recursion.
                if legal_picked and num_slopes_at_best(q, f, covered_intervals) >= k_slopes:
                    new_candidate_faces= generate_candidate_faces(q, f, covered_intervals, (x, y, w), faces_color, sym=sym)
                    for (i, j) in changed_vertices[the_last_changed_v::]:
                        m.set_max(delta[i, j], 0)
                    #exp_dim, new_cs_matrix = dim_cs_matrix(q, changed_vertices, cs_matrix)
                    exp_dim, new_cs_matrix = dim_cs_matrix(q, changed_vertices[0:the_last_changed_v], cs_matrix)
                    # ??? implied equalities are redundant in cs_matrix. don't put them in.
                    if not new_candidate_faces or exp_dim <= dim_threshold:
                        if sym:
                            cs = initial_cs_sym(q, f, vertices_color)
                        else:
                            cs = initial_cs(q, f, vertices_color)
                        polytope = C_Polyhedron(cs)
                        #print >> filename, "exp_dim = %s, got " % exp_dim,
                        #print >> filename, polytope
                        yield polytope, covered_intervals, exp_dim
                    else:
                        for result_polytope, result_covered_intervals, result_exp_dim in paint_complex_combined_mip(k_slopes, q, f, \
                                vertices_color, faces_color, covered_intervals, new_candidate_faces, new_cs_matrix, sym=sym):
                            yield result_polytope, result_covered_intervals, result_exp_dim
        # Now, try out white triangle (x, y, w)
        for (x1, y1, w1) in changed_faces:
            faces_color[(x1, y1, w1)] = 1
            if sym:
                faces_color[(q - y1 - 1, q - x1 - 1, 1 - w1)] = 1
        for (i, j) in changed_vertices:
            vertices_color[i, j] = 1
            if sym:
                vertices_color[q - j, q - i] = 1
            m.set_max(delta[(i, j)], None)
        faces_color[(x, y, w)] = 2
        set_subadd_lower_bound((x, y, w))
        if sym:
            faces_color[(q - y - 1, q - x - 1, 1 - w)] = 2
    for (x1, y1, w1) in candidate_faces:
        faces_color[(x1, y1, w1)] = 1
        if sym:
            faces_color[(q - y1 - 1, q - x1 - 1, 1 - w1)] = 1
    m.remove_constraints(range(the_last_constraint,  m.number_of_constraints()))

eps = QQ(1)/4

def set_subadd_lower_bound((x, y, w)):
    if x < y:
        m.add_constraint(delta[(x+1, y)] + delta[(x, y+1)] + delta[(x+w, y+w)], min = eps)
    else:
        m.add_constraint(2 * delta[(x, y+1)] + delta[(x+w, y+w)], min = eps)

def paint_complex_no_implied(k_slopes, q, f, vertices_color, faces_color, last_covered_intervals, candidate_faces, cs_matrix):
    """
    Do not consider implied green things.
    """
    for (x, y, w) in candidate_faces:
        covered_intervals = directly_covered_by_adding_face(last_covered_intervals, (x, y, w), q, f)
        legal_picked, covered_intervals, changed_vertices, changed_faces = update_around_green_face( \
                    q, f, vertices_color, faces_color, covered_intervals, (x, y, w))
        if legal_picked and num_slopes_at_best(q, f, covered_intervals) >= k_slopes:
            new_candidate_faces= generate_candidate_faces(q, f, covered_intervals, last_face=(x, y, w))
            exp_dim, new_cs_matrix = dim_cs_matrix(q, changed_vertices, cs_matrix)
            if not new_candidate_faces or exp_dim <= dim_threshold:
                cs = initial_cs(q, f, vertices_color)
                polytope = C_Polyhedron(cs)
                yield polytope, covered_intervals
            else:
                for result_polytope, result_covered_intervals in paint_complex_no_implied(k_slopes, q, f, \
                       vertices_color, faces_color, covered_intervals, new_candidate_faces, new_cs_matrix):
                   yield result_polytope, result_covered_intervals
        for face in changed_faces:
            faces_color[face] = 1
        for v in changed_vertices:
            vertices_color[v] = 1
        faces_color[(x, y, w)] = 2
    for face in candidate_faces:
       faces_color[face] = 1

def all_intervals_covered(q, f, values, last_covered_intervals):
    """
    Input:
        values: [pi(0), pi(1/q), .. pi(1)]
        last_covered_intervals: covered_intervals after painting complex

    Return whether all intervals are covered. incremental computation
    """
    to_cover = [i for i in generate_to_cover(q, last_covered_intervals) \
                if (1 <= i < (f + 1) / 2) or ((f + 1) <=  i < (f + q + 1) / 2)]
    if not to_cover:
        return True
    covered_intervals = copy(last_covered_intervals)
    uncovered_intervals = []
    for i in to_cover:
        uncovered_intervals.append(set([i, (f - i - 1) % q]))
    add_v = numpy.ones((q+1,q+1), bool)
    for x in range(q+1):
        for y in range(x, q+1):
            add_v[x, y] = add_v[y, x] = (values[x] + values[y] == values[(x + y) % q])
    was_connected = set([tuple(sym_pair) for sym_pair in uncovered_intervals]) #connected by symmetry
    was_face = set([])
    uncovered_set = set(to_cover)
    # directly covered
    for x in to_cover:
        if x in uncovered_set:
            covered_intervals, uncovered_intervals = update_directly_cover( \
                        q, x, add_v, was_face, was_connected, covered_intervals, uncovered_intervals)
            if not uncovered_intervals:
                return 'direct'
            uncovered_set = generate_uncovered_set(q, uncovered_intervals)
            if x in uncovered_set: # consider symmetry
                covered_intervals, uncovered_intervals = update_directly_cover( \
                        q, (f - x - 1) % q, add_v, was_face, was_connected, covered_intervals, uncovered_intervals)
                if not uncovered_intervals:
                    return 'direct'
                uncovered_set = generate_uncovered_set(q, uncovered_intervals)
                if x in uncovered_set:
                    #check white strip
                    if white_strip(q, f, x, add_v):
                        # x cannot be covered. not extreme function
                        return False

    to_cover = [i for i in uncovered_set if (1 <= i < (f + 1) / 2) or ((f + 1) <=  i < (f + q + 1) / 2)]
    # undirectly covered
    for x in to_cover:
        if x in uncovered_set:
            covered_intervals, uncovered_intervals = update_undirectly_cover( \
                        q, x, add_v, was_face, was_connected, covered_intervals, uncovered_intervals)
            if not uncovered_intervals:
                return 'undirect'
            uncovered_set = generate_uncovered_set(q, uncovered_intervals)
            if x in uncovered_set: # consider symmetry
                covered_intervals, uncovered_intervals = update_undirectly_cover( \
                        q, (f - x - 1) % q, add_v, was_face, was_connected, covered_intervals, uncovered_intervals)
                if not uncovered_intervals:
                    return 'undirect'
                uncovered_set = generate_uncovered_set(q, uncovered_intervals)
    return False

def update_directly_cover(q, x, add_v, was_face, was_connected, covered_intervals, uncovered_intervals):
    """
    Look for green triangles in 2d-complex that covers x.
    Bail out if x is covered.
    Update was_face, was_connected. 
    Return new covered_intervals and uncovered_intervals.
    """
    for y in range(q):
        for w in range(2):
            # I proj to x; K proj to x
            for face in [ (x, y, w), ((x - y) % q, (y - w) % q, w) ]:
                if all(add_v[v] for v in vertices_of_face(face)) and not face in was_face:
                    was_face.add(face)
                    was_connected.update(connected_pair_of_face(face, q))
                    covered_intervals, uncovered_intervals = update_covered_uncovered_by_adding_face( \
                                                     covered_intervals, uncovered_intervals, face, q)
                    return covered_intervals, uncovered_intervals
    return covered_intervals, uncovered_intervals

def update_undirectly_cover(q, x, add_v, was_face, was_connected, covered_intervals, uncovered_intervals):
    """
    Look for green edges in 2d-complex that connects x to others.
    Bail out if all covered.
    Update was_face, was_connected. 
    Return new covered_intervals and uncovered_intervals.
    """
    for y in range(1, q):
        # forward and backward translation to x.
        for u in [x, (x - y) % q]:
            connected = translation_pair(u, y, q)
            if add_v[u, y] and add_v[u + 1, y] and not connected in was_connected:
                was_connected.add(connected)
                covered_intervals, uncovered_intervals = update_covered_uncovered_by_adding_edge( \
                                   covered_intervals, uncovered_intervals, set(connected), q)
                if not uncovered_intervals:
                    return covered_intervals, uncovered_intervals
    for y in range(x) + range(x + 1, q):
        # reflection
        connected = sort_pair(x, y)
        if add_v[x, y + 1] and add_v[x + 1, y] and not connected in was_connected:
            was_connected.add(connected)
            covered_intervals, uncovered_intervals = update_covered_uncovered_by_adding_edge( \
                               covered_intervals, uncovered_intervals, set(connected), q)
            if not uncovered_intervals:
                return covered_intervals, uncovered_intervals
    return covered_intervals, uncovered_intervals

def white_strip(q, f, x, add_v):
    """
    Check if column and diagonal strip of x are all white.
    If so, x cannot be covered. Vertex function is not extreme.
    """
    for xx in [x, (f - x - 1) % q]:
        for y in range(1, q):
            # forward and backward translation to xx.
            for u in [xx, (xx - y) % q]:
                if add_v[u, y] and add_v[u + 1, y]:
                    return False
        for y in range(xx) + range(xx + 1, q):
            # reflection except for symmetry
            if add_v[xx, y + 1] and add_v[xx + 1, y] and ((xx + y + 1) % q != f):
                return False
    return True

#@cached_function
def vertices_of_face(face):
    """
    Return 3 vertices of the given triangle face.
    """
    (x, y, w) = face
    return [(x + w, y + w), (x + 1, y), (x, y + 1)]

#@cached_function
def connected_pair_of_face(face, q):
    """
    Return 2 translation and 1 reflection corresponding to 3 edges of the given triangle face.
    """
    (x, y, w) = face
    return [translation_pair(x, y + w, q), translation_pair(y, x + w, q), sort_pair(x, y)]

#@cached_function
def translation_pair(x, y, q):
    """
    Return the translation (2 components (first < second), they are connected)
    corresponding to the horizontal edge ( (x, y), (x + 1, y) )
    """
    return sort_pair(x, (x + y) % q)

def sort_pair(x, y):
    if x <= y:
        return (x, y)
    else:
        return (y, x)
    
def initialization_sym(q, f):
    f2 = int(f / 2)
    f3 = int(f / 3)
    vertices_color = initial_vertices_color(q, f)
    # impose some initial green vertices
    changed_vertices = [(1, f2), (f2, f2)]
    changed_vertices += [(2, f2 - 1), (2, f2), (f2 - 1, f2)]
    if f % 3 == 0: # Not very promising...
        changed_vertices += [(f3 - 1, f3 - 1), (f3 - 1, f3), (f3 - 1, f3 + 1), (f3 - 1, f3 + 2), \
                             (f3 - 2, f3 + 1), (f3, f3 + 1), (f3 + 1, f3 + 1)]
    elif f % 3 == 1:
        changed_vertices += [(f3 - 2, f3 + 1), (f3 - 2, f3 + 2), (f3 - 1, f3 + 1), (f3, f3), \
                             (f3, f3 + 1), (f3 + 1, f3 + 1), (f3 + 1, f3 + 2)]
    elif f % 3 == 2:
        changed_vertices += [(f3 - 1, f3 + 1), (f3 - 1, f3 + 2), (f3, f3 + 1), (f3 + 1, f3 + 1), (f3 + 1, f3 + 2)]
    # impose their symmetric vertices too
    for (i, j) in changed_vertices:
        vertices_color[i, j] = vertices_color[q - j, q - i] = 0
    # construct the initial constraint matrix
    cs_matrix = matrix(QQ, 1 + f2 + f, q, sparse = True)
    # fn[0] = 0, fn[f] = 1
    cs_matrix[0, 0] = 1
    cs_matrix[1, f] = 1
    # symmetry
    for i in range(1, f2 + 1):
        cs_matrix[1 + i, i] = 1
        cs_matrix[1 + i, f - i] = 1
    # fn[i] = fn[q - i]
    for i in range(1, f):
        cs_matrix[1 + f2 + i, i] = 1
        cs_matrix[1 + f2 + i, q - i] = -1
    # update cs_matrix according to the changed_vertices that we imposed 
    exp_dim, new_cs_matrix = dim_cs_matrix(q, changed_vertices, cs_matrix)
    # construct the initial LP problem
    global m
    global delta
    global var_id
    var_id = numpy.zeros((q, q), int)
    m = MixedIntegerLinearProgram(maximization=True, solver = "GLPK") # "GLPK" slow! "ppl" very slow! try solver = "Gurobi"
    fn = m.new_variable(real=True, nonnegative=True)
    for i in range(f + 1):
        # 0 <= fn[i] <= 1
        m.set_min(fn[i], 0)
        m.set_max(fn[i], 1)
    m.set_max(fn[0], 0) # fn[0] == 0
    #m.set_max(fn[q], 0) # fn[q] == 0
    m.set_min(fn[f], 1) # fn[f] == 1
    for i in range(f):
        m.add_constraint(fn[i] == fn[q - i])
    num_var = q + 1
    # only consider delta[x, y] where x <= y and  x + y <= q
    delta = m.new_variable(real=True, nonnegative=True) # delta[x, y] >= 0
    for x in range(1, f + 1):
        for y in range(x, q - x + 1):
            z = (x + y) % q # need this?
            # delta[x, y] = fn[x] + fn[y] - fn[z]
            m.add_constraint(fn[x] + fn[y] - fn[z] - delta[x, y], min=0, max=0)
            m.set_min(delta[x, y], 0)           # fn[x] + fn[y] >= fn[z]
            if vertices_color[x, y] == 0:
                m.set_max(delta[x, y], 0)       # fn[x] + fn[y] == fn[z]
            else:
                m.set_max(delta[x, y], None)
            var_id[x, y] = num_var
            num_var += 1
    m.set_objective(None)
    m.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
    m.solver_parameter("obj_upper_limit", 0.01)
    
    faces_color, covered_intervals = initial_faces_color_and_covered_intervals(q, f, vertices_color)
    # impose some initial non-candidate faces
    for x in range(f):
        for y in range(x, q - x):
            if x + y == q - 1:
                w_range = [0]
            else:
                w_range = [0, 1]
            for w in w_range:
                if x + y + w >= f:
                    faces_color[(x, y, w)] = faces_color[(q - y - 1, q - x - 1, 1 - w)] = 2
                    set_subadd_lower_bound((x, y, w))
                elif x >= 1 and (y >=  f2 + 1 or x + y + w <= f2):
                    faces_color[(x, y, w)] = faces_color[(q - y - 1, q - x - 1, 1 - w)] = 2
                    set_subadd_lower_bound((x, y, w))
    candidate_faces = generate_candidate_faces(q, f, covered_intervals, last_face=None, faces_color=faces_color, sym = True)
    return vertices_color, faces_color, covered_intervals, candidate_faces, new_cs_matrix

def initial_cs_sym(q, f, vertices_color):
    cs = Constraint_System()
    fn = [ Variable(i) for i in range(q+1) ]

    cs.insert(fn[0] == 0)
    cs.insert(fn[q] == 0)
    cs.insert(fn[f] == 1)
    for i in range(1, f):
        cs.insert(fn[i] >= 0)
        cs.insert(fn[i] <= 1)
        cs.insert(fn[i] - fn[q - i] == 0)
    for x in range(1, f + 1):
        for y in range(x, q - x + 1):
            if vertices_color[x, y] == 0:
                cs.insert(fn[x] + fn[y] == fn[x + y])
            else:
                cs.insert(fn[x] + fn[y] >= fn[x + y])
    return cs

def gen_initial_polytope_sym(q, f, vertices_color, covered_intervals):
    cs = initial_cs_sym(q, f, vertices_color)
    polytope = C_Polyhedron(cs)
    yield polytope, covered_intervals

output_dir = "./"
## override this in a file named "config.sage" if necessary

#### The following are no longer used.
## dir_math ="/homes/home02/y/yzh/Dropbox/group-relaxation-sage-code/"
## dir_yuan ="/media/sf_dropboxcode/"

def save_plot(q, hh, destdir = output_dir+"sym_mode_2d_diagrams/"):
    mkdir_p(destdir)
    logging.disable(logging.info)
    for i in range(len(hh)):
        v_n = hh[i]
        h = h_from_vertex_values(v_n)
        num = len(set([v_n[j+1] - v_n[j] for j in range(q)]))
        name = "%sq%s_%s" %(num, q, i)
        if not extremality_test(h):
            name += "_notextreme"
        name += ".png"
        g = plot_2d_diagram(h, colorful=True)
        g.save(destdir + name, figsize = 20, show_legend=False)
    logging.disable(logging.NOTSET)

def pattern_q_and_f(l, pattern):
    """
    pattern == 0 corresponds to the pattern described in the paper.
    In the paper, r is what is called l here.
    """
    if pattern <= 3:
        f = 18 * l + 11
    elif pattern == 4 or pattern == 5:
        f = 18 * l + 5
    elif pattern == 6:
        f = 18 * l + 13
    elif pattern == 7:
        f = 18 * l + 15
    q = 2 * f
    return q, f

def pattern_additive_vertices(l, pattern):
    """
    Return a list of pairs (i, j) of integers i <= j
    such that (i, j)/q is a prescribed additivity.

    These additivities form the additive triangles.

    pattern == 0 corresponds to the pattern described in the paper.
    Remaining patterns are undocumented.
    In the paper, r is what is called l here.
    """
    q, f = pattern_q_and_f(l, pattern)
    f2 = int(f / 2)
    f3 = int(f / 3)
    changed_vertices = []
    # impose some initial green vertices
    for k in range(f3 + 2, f2, 3):
        changed_vertices += [(k, k), (k, k + 1), (k, k + 2), (k - 1, k - 1), (k - 1, k), (k - 1, k + 2), (k + 1, k + 2)]
        if pattern <= 1 or pattern == 4 or pattern == 5 and k < f2 - 2 or pattern == 6 or pattern == 7:
            changed_vertices +=  [(k, k + 3)]
        if pattern == 3:
            changed_vertices += [(k - 1, k + 1), (k + 1, k + 1)]
    changed_vertices += [(1, f2), (f2, f2)]
    if pattern <= 3 or pattern == 6 or pattern == 7:
        changed_vertices += [(f2 - 1, f2 - 1), (f2 - 1, f2), (2, f2 - 1), (2, f2), (3, f2 - 1)]
    if pattern == 3:
        changed_vertices += [(f2 - 1, f2 + 1), (1, f2 - 1), (1, f2 + 1)]
    if pattern == 0:
        changed_vertices += [(f3, f3), (f3, f3 + 2)]
    if pattern == 6:
        changed_vertices += [(f3, f3), (f3, f3 + 1)]
    if pattern == 7:
        changed_vertices += [(f3 - 1, f3 - 1), (f3 - 1, f3), (f3 - 1, f3 + 1), (f3 - 1, f3 + 2), (f3, f3 + 1)]
    for k in range(1, l + 1):
        if pattern <= 3 or pattern == 6 or pattern == 7:
            i = 6 * k - 1
            j = f2 - 3 * k + 1
        elif pattern == 4 or pattern == 5:
            i = 6 * k - 3
            j = f2 - 3 * k + 2 
        changed_vertices += [(i - 1, j), (i - 1, j + 1), (i, j - 1), (i, j + 1), \
                         (i + 1, j - 2), (i + 1, j - 1), (i + 1, j), (i + 1, j + 1), \
                         (i + 2, j - 1), (i + 3, j - 2), (i + 3, j - 1), (i + 4, j - 2)]
        if pattern <= 1 or pattern == 4 or pattern == 5 and k > 1 or pattern == 6 or pattern == 7:
            changed_vertices += [(i - 1, j - 1), (i - 1, j + 2)]
        if pattern == 3:
            changed_vertices += [(i, j), (i + 2, j), (i + 2, j - 2)]
    return changed_vertices

def pattern_more_additive_vertices(l, pattern):
    """
    Extra additive points, to be added to pattern_additive_vertices
to reduce the dimension.
    """
    q, f = pattern_q_and_f(l, pattern)
    f2 = int(f / 2)
    f3 = int(f / 3)
    more_additivity = [(f3 + 2, q - 4), (f3 + 4, q - 10), (f2 - 2 , f2 - 2), (f3 + 4, q - 12), \
                       (f2 - 8, f2 - 8), (f2 - 17, f2 - 17), (f2 - 23, f2 - 23), \
                       (f3 + 8, q - 22), (f3 + 7, q - 19), (f2 - 14, f2 - 14)]
    ## components merge: l>=3, s1=s2; l>=5, s3=s4; l>=10, s6=s7; l>=13, s8=s9; l>=16, s12=s13.
    l_seuil = [0, 1, 3, 4, 5, 10, 13, 15, 15, 16]
    changed_vertices = [more_additivity[i] for i in range(len(l_seuil)) if l >= l_seuil[i]]
    return changed_vertices

def pattern_vertices_color(l, pattern=0, more_ini_additive=False):
    q, f = pattern_q_and_f(l, pattern)
    vertices_color = initial_vertices_color(q, f)
    changed_vertices = pattern_additive_vertices(l, pattern)
    if pattern == 0 and more_ini_additive:
        changed_vertices += pattern_more_additive_vertices(l, pattern)
    # impose their symmetric vertices too
    for (i, j) in changed_vertices:
        vertices_color[i, j] = vertices_color[q - j, q - i] = 0
    return vertices_color

def pattern_s(l, pattern):
    if pattern <= 1 or pattern == 6 or pattern == 7:
        s = [Variable(0), Variable(1)]
        for k in range(1, l + 1):
            s += [Variable(k), Variable(k + 1)] * 3
        if pattern == 0:
            s += [Variable(l + 1), Variable(l + 1), Variable(l + 1)]
        elif pattern == 1:
            s += [Variable(l + 1), Variable(l + 2), Variable(l + 1)]
        elif pattern == 6:
            s += [Variable(l + 1), Variable(l + 2), Variable(l + 2), Variable(l + 1)]
        elif pattern == 7:
            s += [Variable(l + 1), Variable(l + 2), Variable(l + 1), Variable(l + 2), Variable(l + 1)]
        for k in range(l, 0, -1):
            s += [Variable(k), Variable(k + 1), Variable(k)]
    elif pattern == 2:
        s = [Variable(0), Variable(1), Variable(1)]
        for k in range(1, l + 1):
            s += [Variable(2 * k), Variable(2 * k + 1), Variable(2 * k), Variable(2 * k + 1), Variable(2 * k), Variable(2 * k)]
        s += [Variable(2 * l + 2)]
        for k in range(l, 0, -1):
            s += [Variable(2 * k), Variable(2 * k + 1), Variable(2 * k)]
        s += [Variable(1)]
    elif pattern == 3:
        s = [Variable(0)] * 3
        for k in range(1, l + 1):
            s += [Variable(k)] * 6
        #s += [Variable(l + 1)]
        s += [-Variable(0)]
        for k in range(l, 0, -1):
            s += [Variable(k)] * 3
        s += [Variable(0)]
    elif pattern == 4:
        s = []
        for k in range(l):
            s += [Variable(k), Variable(k + 1)] * 3
        for k in range(l, 0, -1):
            s += [Variable(k), Variable(k + 1), Variable(k)]
        s += [Variable(0), Variable(1)]
    elif pattern == 5: # inverse
        s = [Variable(l + 2)]
        for k in range(l, 0, -1):
            s += [Variable(k), Variable(k + 1), Variable(k), Variable(k + 1), Variable(k), Variable(k)]
        s += [Variable(0)]
        for k in range(1, l+1):
            s += [Variable(k), Variable(k + 1), Variable(k)]
    return s + [s[0]] + s[-1::-1]

def pattern_fn(l, pattern):
    s = pattern_s(l, pattern)
    f = len(s)
    q = 2 * f
    fn = [Linear_Expression(0)] * (q + 1)
    for k in range(f):
        fn[k + 1] = fn[k] + s[k]
    for k in range(f):
        fn[q - k] = fn[k]
    return fn

def pattern_polytope(vertices_color, fn):
    q = len(fn) - 1
    f = int(q / 2)
    cs = Constraint_System()
    cs.insert(fn[0] == 0)
    cs.insert(fn[f] == 1)
    for i in range(1, f):
        cs.insert(fn[i] >= 0)
        cs.insert(fn[i] <= 1)
    for x in range(1, f + 1):
        for y in range(x, q - x + 1):
            if vertices_color[x, y] == 0:
                cs.insert(fn[x] + fn[y] == fn[x + y])
            else:
                cs.insert(fn[x] + fn[y] >= fn[x + y])
    polytope = C_Polyhedron(cs)
    return polytope

def pattern_extreme(l, k_slopes, pattern=0, show_plots=False,
                    test_extremality=False, polytope = None,
                    more_ini_additive=True, count_components=False, use_sha1=True):
    r"""
    Computes the functions corresponding to the extreme points of the
    polytope corresponding to subadditivities and prescribed additivities,
    according to `pattern` and size parameter `l` (called r in the paper).

    If `test_extremality` is True (default is False), check extremality.

    Prints lines:
    NUM-SLOPES  SV (input to pattern0_sym_fn)  DENOMINATOR

    Creates files in the directory 'sym_mode_2d_diagrams/patterns_PATTERN/'
    (this directory is created if it does not exist).

    Returns the max number of slopes found, regardless of extremality test.

    EXAMPLES::

        sage: pattern_extreme(3, 4, 0, test_extremality=True, count_components=True)
        #####  ...
        q = 130; f = 65; num_components = 8; num_slopes = 8; divisor = 205; extreme = True
        l = 3; sv = (49, 9, 9, 3, -19)
        h = pattern0_sym_fn(l, sv)
        #####  ...
        q = 130; f = 65; num_components = 6; num_slopes = 6; divisor = 33; extreme = True
        l = 3; sv = (9, 1, 1, 1, -3)
        h = pattern0_sym_fn(l, sv)
        8
    """
    q, f = pattern_q_and_f(l, pattern)
    fn = pattern_fn(l, pattern)
    if polytope is None:
        vertices_color = pattern_vertices_color(l, pattern, more_ini_additive=more_ini_additive)
        polytope = pattern_polytope(vertices_color, fn)
    v_set = set([])
    vv = []
    nn = []
    destdir = output_dir+"sym_mode_2d_diagrams/"+"patterns_%s/" % pattern
    mkdir_p(destdir)
    logging.disable(logging.info)
    for v in polytope.minimized_generators():
        #v.coefficients() is numerator of component's slope value
        v_n = [sum(p*q for p, q in zip(fn[i].coefficients(), v.coefficients())) for i in range(q+1)]
        num = len(set([v_n[i+1] - v_n[i] for i in range(q)]))
        if num >= k_slopes:
            if not tuple(v_n) in v_set:
                v_set.add(tuple(v_n))
                vv.append(v_n)
                nn.append(num)
                #print num, v.coefficients(), v.divisor()
                # Note: Imposing the invariance under x -> 1-x may
                # result in a non-extreme all covered vertex-function.
                h_is_extreme = 'untested'
                h = None
                if test_extremality or count_components or show_plots or use_sha1:
                    # all covered, continuous => extreme iff system of equations has unique solution
                    h = h_from_vertex_values(v_n)
                if test_extremality:
                    h_is_extreme = simple_finite_dimensional_extremality_test(h, show_plots=False, f=None, oversampling=None, order=q, show_all_perturbations=False)
                num_components = 'not computed'
                if count_components:
                    num_components = len(generate_covered_intervals(h))
                if use_sha1:
                    id = h.sha1()
                else:
                    id = len(vv)
                sage_name = "%sq%s_%s.sage" %(num, q, id)
                info = "q = %s; f = %s; num_components = %r; num_slopes = %s; divisor = %s; extreme = %r\n" % (q, f, num_components, num, v.divisor(), h_is_extreme)
                info += "l = %s; " % l
                info += "sv = %s\n" % (v.coefficients(),)
                if pattern != 0:
                    info += "v_n = %s\n" % (v_n,)
                    info += "h = h_from_vertex_values(v_n)\n"
                else:
                    info += "h = pattern%s_sym_fn(l, sv)\n" % pattern  # this function only exists for pattern=0
                print "##### ", destdir + sage_name
                print info,
                with open(destdir + sage_name, "w") as sage_file:
                    print >> sage_file, info
                if show_plots: # and h_is_extreme:
                    name = "%sq%s_%s.png" %(num, q, id)
                    g = plot_2d_diagram(h, colorful=True)
                    figsize = 10 * l
                    g.save(destdir + name, figsize = figsize, show_legend=False)
                    print "# Plot saved in", destdir + name
    logging.disable(logging.NOTSET)
    #return vv, nn
    if len(nn) > 0:
        return max(nn)
    else:
        return 0

def plot_pattern(l, more_ini_additive = True, show_plots=True):
    f = 18 * l + 11
    q = 2 * f
    s=[0, 1]
    for k in range(1, l + 1):
        s += [k, k + 1] * 3
    s += [l + 1] * 3
    for k in range(l, 0, -1):
        s += [k, k + 1, k]
    sk = s + [s[0]] + s[-1::-1]
    nc = max(sk)+1
    sc = [nc - i - 1 for i in sk] + [2*nc - i - 1 for i in sk]
    colors = rainbow(2*nc)
    vertices_color = pattern_vertices_color(l, pattern=0, more_ini_additive=more_ini_additive)
    for i in range(q):
        for j in range (i):
            vertices_color[i,j]=vertices_color[j,i]
    for i in range(q):
        vertices_color[q][i]=vertices_color[i][q]=0
    vertices_color[q][q]=0
    g = Graphics()

    vertical_edge = set([])
    horizontal_edge = set([])
    pts = set([])
    for i in range(q):
        for j in range(q):
            #if vertices_color[i,j]==0:
            #    g += points([(i,j)])
            if vertices_color[i+1,j]==0 and vertices_color[i,j+1]==0:
                if vertices_color[i,j]==0:
                    g += polygon([(i,j),(i+1,j),(i,j+1)], color=colors[sc[i]], fill=True)
                    pts.add((i,j));
                    vertical_edge.add((i,j))
                    horizontal_edge.add((i,j))
                if vertices_color[i+1, j+1]==0:
                    g += polygon([(i+1,j+1),(i+1,j),(i,j+1)], color=colors[sc[i]], fill=True)
                    pts.add((i+1, j+1))
                    vertical_edge.add((i+1,j))
                    horizontal_edge.add((i,j+1))
                if vertices_color[i,j] == 1 and vertices_color[i+1, j+1] == 1:
                    g += line([(i+1,j),(i,j+1)], color='springgreen',linestyle='-')
                pts.add((i+1,j));
                pts.add((i,j+1));
            if vertices_color[i,j]==0 and vertices_color[i,j+1]==0 and not ((i,j) in vertical_edge):
                g += line([(i,j),(i,j+1)], color='springgreen',linestyle='-')
                vertical_edge.add((i,j))
                pts.add((i,j))
                pts.add((i,j+1))
            if vertices_color[i,j]==0 and vertices_color[i+1,j]==0 and not ((i,j) in horizontal_edge):
                g += line([(i,j),(i+1,j)], color='springgreen',linestyle='-')
                horizontal_edge.add((i,j))
                pts.add((i,j))
                pts.add((i+1,j))
            if vertices_color[i,j]==0 and not ((i,j) in pts):
                g += points([(i,j)], color='springgreen')
                pts.add((i,j))
    if show_plots:
        g.show(gridlines=[range(q+1),range(q+1)],gridlinesstyle=dict(color="grey", linestyle=":"),figsize=10*l)
    else:
        return g

def write_panda_format_cs(cs, fname=None, newcode=True):
    if fname:
        dir = output_dir+"profiler/panda/"
        mkdir_p(dir)
        filename = open(dir+fname, "w")
    else:
        filename = sys.stdout
    if newcode:
        panda_string = convert_pplcs_to_panda(cs)
        print >> filename, panda_string
    else:
        print >> filename, 'Inequalities:'
        for c in cs:
            for x in c.coefficients():
                print >> filename, -x,
            print >> filename, -c.inhomogeneous_term()
            if c.is_equality():
                for x in c.coefficients():
                    print >> filename, x,
                print >> filename, c.inhomogeneous_term()
    if fname:
        filename.close()
    return

def convert_pplcs_to_panda(cs):
    s = "Names\n"
    q = cs.space_dimension() - 1;
    for i in range(1, q+2):
        s += "x%s " % i
    s += "\n"
    s += 'Equations\n'
    for c in cs:
        if c.is_equality():
            coefs = c.coefficients()
            for i in range(q+1):
                x = coefs[i]
                if x > 0:
                    s += '+%sx%s ' % (x, i+1)
                elif x < 0:
                    s += '%sx%s ' % (x, i+1)
            s += '= '
            s += '%s\n' % -c.inhomogeneous_term()
    s += 'Inequalities\n'
    for c in cs:
        if not c.is_equality():
            coefs = c.coefficients()
            for i in range(q+1):
                x = coefs[i]
                if x > 0:
                    s += '+%sx%s ' % (x, i+1)
                elif x < 0:
                    s += '%sx%s ' % (x, i+1)
            s += '>= '
            s += '%s\n' % -c.inhomogeneous_term()
    return s

def write_porta_ieq(q, f, destdir=None):
    vertices_color = initial_vertices_color(q, f);
    cs = initial_cs(q, f, vertices_color)
    if destdir is None:
        fname = None
    else:
        dir = output_dir+"profiler/porta/"
        mkdir_p(dir)
        fname = dir + "porta_q%sf%s.ieq" % (q, f)
    write_porta_format_cs(cs, q=q, f=f, fname=fname)
    return

def write_porta_format_cs(cs, q=None, f=None, fname=None):
    # q, f are used to find a valid point -- gmic function
    # this valid point is required by 'traf' as the 0 is not in the feasible region.
    if fname:
        filename = open(fname, "w")
    else:
        filename = sys.stdout
    porta_string = convert_pplcs_to_porta(cs, q=q, f=f)
    print >> filename, porta_string
    if fname:
        filename.close()
    return

def convert_pplcs_to_porta(cs, q=None, f=None):
    # q, f are used to find a valid point -- gmic function
    # this valid point is required by 'traf' as the 0 is not in the feasible region.
    s = ""
    s += "DIM = %s\n" % cs.space_dimension()
    s += "\n"
    if not q is None:
        s += 'VALID\n'
        s += '0 '
        for i in range(1, f):
            s += '%s/%s ' % (i, f)
        s += '1 '
        for i in range(q - f - 1, 0, -1):
            s += '%s/%s ' % (i, q - f)
        s += '0\n'
        s += '\n'
    s += 'INEQUALITIES_SECTION\n'
    for c in cs:
        coefs = c.coefficients()
        for i in range(q+1):
            x = coefs[i]
            if x > 0:
                s += '+%sx%s ' % (x, i+1)
            elif x < 0:
                s += '%sx%s ' % (x, i+1)
        if c.is_equality():
            s += '== '
        else:
            s += '>= '
        s += '%s\n' % -c.inhomogeneous_term()
    s += '\n'
    s += 'END\n'
    return s

def write_lrs_ine(q, f, destdir=None):
    vertices_color = initial_vertices_color(q, f);
    cs = initial_cs(q, f, vertices_color)
    if destdir is None:
        fname = None
    else:
        dir = output_dir+"profiler/lrs/"
        mkdir_p(dir)
        fname = dir + "lrs_q%sf%s.ine" % (q, f)
    write_lrs_format_cs(cs, fname=fname)
    return

def write_lrs_format_cs(cs, fname=None):
    if fname:
        filename = open(fname, "w")
    else:
        filename = sys.stdout
    lrs_string = convert_pplcs_to_lrs(cs, fname=fname)
    print >> filename, lrs_string
    if fname:
        filename.close()
    return

def convert_pplcs_to_lrs(cs, fname=None):
    if fname:
        s = fname + '\n'
    else:
        s = ""
    s += "H-representation" + '\n'
    m = len(cs)
    n = cs.space_dimension() + 1
    k = 0
    linearities = []
    for i in range(m):
        c = cs[i]
        if c.is_equality():
            k += 1
            linearities.append(i)
    if k > 0:
        s += "linearity %s " % k
    for i in linearities:
        s += repr(i + 1)  + ' '
    s += '\n'
    s += "begin\n"
    s += "%s %s rational\n" %(m, n)
    for c in cs:
        s += repr(c.inhomogeneous_term()) + ' '
        for x in c.coefficients():
            s += repr(x) + ' '
        s += '\n'
    s += 'end\n'
    return s

def read_lrs_to_cs_or_gs(fname):
    myfile = open(fname, "r")
    lrs_string = myfile.read()
    myfile.close()
    return convert_lrs_to_ppl(lrs_string)

def convert_lrs_to_ppl(lrs_string):
    """
    Convert lrs format H-representation to ppl cs;
    or lrs format V-representation (of a polytope) to ppl gs.

    COPY from src/geometry/polyhedron/backend_cdd.py and edit
    Polyhedron_cdd._init_from_cdd_output(self, cdd_output_string)
    """
    cddout=lrs_string.splitlines()

    # nested function
    def expect_in_cddout(expected_string):
        l = cddout.pop(0).strip()
        if l != expected_string:
            raise ValueError, ('Error while parsing cdd output: expected "'
                               +expected_string+'" but got "'+l+'".\n' )
    # nested function
    def cdd_linearities():
        l = cddout[0].split()
        if l[0] != "linearity":
            return []
        cddout.pop(0)
        assert len(l) == int(l[1])+2, "Not enough linearities given"
        return [int(i)-1 for i in l[2:]]  # make indices pythonic

    # nested function
    def cdd_convert(string, field=QQ):
        """
        Converts the cdd output string to a QQ numerical value.
        """
        return [field(x) for x in string.split()]

    # nested function
    def find_in_cddout(expected_string):
        """
        Find the expected string in a list of strings, and
        truncates ``cddout`` to start at that point. Returns
        ``False`` if search fails.
        """
        for pos in range(0,len(cddout)):
            l = cddout[pos].strip();
            if l==expected_string:
                # must not assign to cddout in nested function
                for i in range(0,pos+1):
                    cddout.pop(0)
                return True
        return False

    def lrs_row_to_linear_expression(l):
        rational_list = cdd_convert(l)
        num_list = [x.numerator() for x in rational_list]
        den_list = [x.denominator() for x in rational_list]
        common_den = lcm(den_list)
        ihom = int(common_den / den_list[0]) * num_list[0]
        coef = [int(common_den / den_list[i]) * num_list[i] for i in range(1, len(rational_list))]
        return Linear_Expression(coef, ihom)

    def lrs_row_to_point(l):
        rational_list = cdd_convert(l)
        num_list = [x.numerator() for x in rational_list]
        den_list = [x.denominator() for x in rational_list]
        common_den = lcm(den_list)
        coef = [int(common_den / den_list[i]) * num_list[i] for i in range(len(rational_list))]
        return sage.libs.ppl.point(Linear_Expression(coef, 0), common_den)

    if find_in_cddout('V-representation'):
        # Suppose it's the V-representation of a polytope.
        # Return the ppl generator_system that contains the vertices.
        # Sometimes the number of vertices is unknown from the input file, so read till 'end'
        # ex: in the output of lrs vertex enumertaion file.ext.
        #raise NotImplementedError, "V-representation Not implemented."
        gs = Generator_System()
        equations = cdd_linearities()
        expect_in_cddout('begin')
        l = cddout.pop(0).split()
        n = int(l[1]) - 1
        l = cddout.pop(0).strip()
        while l != 'end':
            l_type = l[0]
            l = l[1:]
            #if (i in equations) or l_type == '0':
            #    raise NotImplementedError, "V-representation of line or ray is NOT implemented."
            vertex = lrs_row_to_point(l)
            gs.insert(vertex)
            l = cddout.pop(0).strip()
        return gs

    if find_in_cddout('H-representation'):
        cs = Constraint_System()
        equations = cdd_linearities()
        expect_in_cddout('begin')
        l = cddout.pop(0).split()
        m = int(l[0])
        n = int(l[1]) - 1
        fn = [ Variable(i) for i in range(n) ]
        for i in range(m):
            l = cddout.pop(0).strip()
            lin_exp = lrs_row_to_linear_expression(l)
            if i in equations:
                cs.insert(lin_exp == 0)
            else:
                cs.insert(lin_exp >= 0)
        expect_in_cddout('end')
        return cs

from sage.misc.temporary_file import tmp_filename
from subprocess import Popen, PIPE

def lrs_redund(in_str, verbose=False):
        """
        To remove redundant inequalities from an H-representation or
        input points that are not vertices from a V-representation
        use the command 'redund' from lrslib.
        Input: lrs format in_str; Output: lrs format out_str;

        Copy and edit from def _volume_lrs(self, verbose=False),
        http://www.sagenb.org/src/geometry/polyhedron/base.py
        """
        #if is_package_installed('lrslib') != True:
        #    print 'You must install the optional lrs package ' \
        #          'for this function to work'
        #    raise NotImplementedError

        in_filename = tmp_filename()
        in_file = file(in_filename,'w')
        in_file.write(in_str)
        in_file.close()
        if verbose: print in_str

        redund_procs = Popen(['redund',in_filename],stdin = PIPE, stdout=PIPE, stderr=PIPE)
        out_str, err = redund_procs.communicate()
        if verbose:
            print out_str

        return out_str

def remove_redundancy_from_cs(cs, verbose=False, return_lrs=False):
    """
    Remove redundant inequalities from cs using the command 'redund' from lrslib.
    Return the new cs without redundancy
    in ppl format if return_lrs==False (default), or in lrs format if return_lrs==True 
    """
    in_str = convert_pplcs_to_lrs(cs, fname=None)
    out_str = lrs_redund(in_str, verbose=verbose)
    if return_lrs:
        return out_str
    else:
        return convert_lrs_to_ppl(out_str)

def lrs_lrs(in_str, verbose=False):
    """
    use the command 'lrs' from lrslib.
    Input: lrs format in_str; Output: lrs format out_str;
    """
    #if is_package_installed('lrslib') != True:
    #    print 'You must install the optional lrs package ' \
    #          'for this function to work'
    #    raise NotImplementedError
    in_filename = tmp_filename()
    in_file = file(in_filename,'w')
    in_file.write(in_str)
    in_file.close()
    if verbose: print in_str

    redund_procs = Popen(['lrs',in_filename],stdin = PIPE, stdout=PIPE, stderr=PIPE)
    out_str, err = redund_procs.communicate()
    if verbose:
        print out_str
    return out_str

def lrs_lrsinput_pploutput(in_str):
    """
    use the command 'lrs' from lrslib.
    Input: lrs format in_str; Output: ppl format extreme_points;
    
    EXAMPLES::

        sage: cube_in_str = "cube\n*cube of side 2 centred at origin\nH-representation\nbegin\n6  4 rational" + \
        ...                 "\n1 1 0 0\n1 0 1 0\n1 0 0 1\n1 -1 0 0\n1 0 -1 0\n1 0 0 -1\nend"
        sage: lrs_lrsinput_pploutput(cube_in_str)
        Generator_System {point(1/1, 1/1, 1/1), point(-1/1, 1/1, 1/1), point(1/1, -1/1, 1/1), point(-1/1, -1/1, 1/1), point(1/1, 1/1, -1/1), point(-1/1, 1/1, -1/1), point(1/1, -1/1, -1/1), point(-1/1, -1/1, -1/1)}
        sage: lrs_q5f3_str = "lrs_q5f3\nH-representation\nlinearity 5 1 2 3 13 21\nbegin\n21 7 rational" + \
        ...                  "\n0 1 0 0 0 0 0\n0 0 0 0 0 0 1\n-1 0 0 0 1 0 0\n0 0 1 0 0 0 0\n1 0 -1 0 0 0 0\n0 0 0 1 0 0 0\n1 0 0 -1 0 0 0\n0 0 0 0 1 0 0" + \
        ...                  "\n1 0 0 0 -1 0 0\n0 0 0 0 0 1 0\n1 0 0 0 0 -1 0\n0 0 2 -1 0 0 0\n0 0 1 1 -1 0 0\n0 0 1 0 1 -1 0\n0 -1 1 0 0 1 0" + \
        ...                  "\n0 0 0 2 0 -1 0\n0 -1 0 1 1 0 0\n0 0 -1 1 0 1 0\n0 0 -1 0 2 0 0\n0 0 0 -1 1 1 0\n0 0 0 0 1 -2 0\nend"
        sage: lrs_lrsinput_pploutput(lrs_q5f3_str)
        Generator_System {point(0/6, 2/6, 4/6, 6/6, 3/6, 0/6), point(0/4, 3/4, 1/4, 4/4, 2/4, 0/4)}    
    """
    v_lrs_str = lrs_lrs(in_str)
    extreme_points = convert_lrs_to_ppl(v_lrs_str)
    return extreme_points
   
def lcdd_rational(in_str, verbose=False):
    """
    use the command 'lcdd_gmp' from cddlib.
    Input: cdd format in_str; Output: cdd format out_str;
    """
    in_filename = tmp_filename()
    in_file = file(in_filename,'w')
    in_file.write(in_str)
    in_file.close()
    if verbose: 
        print in_str
    redund_procs = Popen(['lcdd_gmp',in_filename],stdin = PIPE, stdout=PIPE, stderr=PIPE)
    out_str, err = redund_procs.communicate()
    if verbose:
        print out_str
    return out_str

def measure_stats_detail(q, f, prep=True):
    vertices_color = initial_vertices_color(q, f);
    cs = initial_cs(q, f, vertices_color)
    polytope = C_Polyhedron(cs)
    exp_dim = int(q / 2) - 1
    time_start = os.times()
    extreme_points = vertex_enumeration(polytope, prep=prep, exp_dim=exp_dim)
    time_end = os.times()
    t = sum([time_end[i]-time_start[i] for i in range(4)])
    num = len(extreme_points)
    slope = []
    extreme = []
    vdenominator = []
    sdenominator = []
    for v in extreme_points:
        v_n = v.coefficients()
        s = set([v_n[i+1] - v_n[i] for i in range(q)])
        ns = len(s)
        slope.append(ZZ(ns))
        if ns == 2:
            extreme.append(ZZ(1))
        else:
            whether_covered = all_intervals_covered(q, f, v_n, [])
            if whether_covered:
                extreme.append(ZZ(1))
            else:
                extreme.append(ZZ(0))
        divisor = v.divisor()
        vdenominator.append(ZZ(divisor))
        gcd_s = gcd(list(s)+[divisor])
        sdenominator.append(ZZ(divisor / gcd_s))
    return t, ZZ(num), slope, extreme, vdenominator, sdenominator     
    
def write_stats_detail(q, prep=True, destdir=None):
    if destdir is None:
        filename = sys.stdout
    else:
        dir = output_dir+"profiler/sage_input/"
        mkdir_p(dir)
        filename = open(dir + "stats_q%s.sage" % q, "a")
    
    t_enumeration = []
    n_vertex = []
    n_slope = []
    is_extreme = []
    denominator_v = []
    denominator_s = []
    max_denominator_v = []
    max_denominator_s = []
    
    for f in range(1, int(q/2)+1):
        t, num, slope, extreme, vdenominator, sdenominator = measure_stats_detail(q, f, prep=prep)
        t_enumeration.append(t)
        n_vertex.append(num)
        n_slope.append(slope)
        is_extreme.append(extreme)
        denominator_v.append(vdenominator)
        denominator_s.append(sdenominator)
        max_denominator_v.append(max(vdenominator))
        max_denominator_s.append(max(sdenominator)) 

    print >> filename, "t_enumeration =",
    print >> filename, sage_input(t_enumeration)
    print >> filename
    print >> filename, "n_vertex =",
    print >> filename, sage_input(n_vertex)
    print >> filename
    print >> filename, "n_slope =",
    print >> filename, sage_input(n_slope)
    print >> filename
    print >> filename, "is_extreme =",
    print >> filename, sage_input(is_extreme)
    print >> filename
    print >> filename, "denominator_v =",
    print >> filename, sage_input(denominator_v)
    print >> filename
    print >> filename, "denominator_s =",
    print >> filename, sage_input(denominator_s)
    print >> filename
    print >> filename, "max_denominator_v =",
    print >> filename, sage_input(max_denominator_v)
    print >> filename
    print >> filename, "max_denominator_s =",
    print >> filename, sage_input(max_denominator_s)
    
    if destdir:
        filename.close()
    return

def time_to_find_first_extreme(k, q, f, mode='combined'):
    st = os.times(); 
    v = search_kslope_example(k, q, f, mode=mode, prep=True).next(); 
    et = os.times(); 
    t = sum([et[i]-st[i] for i in range(4)]);
    return t, v
    
def times_in_naive_search(k, q, f):
    n_sol = 0
    h_list = []
    logging.info( "Naive search for extreme funtions with q = %s, f = %s, k_slopes >= %s" % (q, f, k))
    time_start = os.times()
    for h in search_kslope_example(k, q, f, mode='naive', prep=True):
        n_sol += 1
        if n_sol == 1:
            time_end = os.times()
            t = sum([time_end[i]-time_start[i] for i in range(4)])
            logging.info("Time to first k-slope = %s" % t)
        h_list.append(h)
    time_end = os.times()
    t = sum([time_end[i]-time_start[i] for i in range(4)])
    logging.info("Total time of naive search = %s" % t)
    logging.info("Numer of k-slope functions = %s" % n_sol)
    return h_list
