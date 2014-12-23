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

# Main: def search_kslope(k_slopes, q, f_list, mode) and def measure_stats()

from sage.libs.ppl import C_Polyhedron, Constraint, Constraint_System, Generator, Generator_System, Variable, point, \
                          Poly_Con_Relation, MIP_Problem, Linear_Expression
import numpy
from sage.numerical.mip import MixedIntegerLinearProgram, MIPVariable, MIPSolverException
import sage.numerical.backends.glpk_backend as backend

poly_is_included = Poly_Con_Relation.is_included()

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
    Look for implied additive vertices and faces, given MIP_Problem m.
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

def generate_vertex_values(k_slopes , q, polytope,  v_set=set([])):
    """
    Return vertices of the polytope, whose corresponding 
    piecewise_linear function h has at least  k_slopes.
    """
    for v in polytope.minimized_generators():
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

def search_kslope_example(k_slopes, q, f, mode='heuristic'):
    """
    Search for extreme functions that have required number of slope values.

    If `mode` is (defaut) 'heuristic', use paint_complex_heuristic() to paint;
    If `mode` is 'combined', use paint_complex_combined() to paint a few triangles, then enumerate vertex-functions;
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
            gen = gen_initial_polytope_sym(q, f, vertices_color, covered_intervals)
        else:
            gen = paint_complex_combined_mip(k_slopes, q, f, vertices_color, faces_color, covered_intervals, candidate_faces, cs_matrix, sym=True)
    elif mode == 'no_implied':
        faces_color, covered_intervals = initial_faces_color_and_covered_intervals(q, f, vertices_color)
        candidate_faces = generate_candidate_faces(q, f, covered_intervals, last_face=None)
        cs_matrix = initial_cs_matrix(q, f)
        gen = paint_complex_no_implied(k_slopes, q, f, vertices_color, faces_color, covered_intervals, candidate_faces, cs_matrix)
    elif mode == 'fulldim_covers':
        cs = initial_cs(q, f, vertices_color)
        faces_color, covered_intervals = initial_faces_color_and_covered_intervals(q, f, vertices_color)
        gen = paint_complex_fulldim_covers(k_slopes, q, f, vertices_color, faces_color, covered_intervals, (0, 0, 0), cs)
    elif mode == 'complete':
        cs = initial_cs(q, f, vertices_color)
        covered_intervals, uncovered_intervals = initial_covered_uncovered(q, f, vertices_color)
        gen = paint_complex_complete(k_slopes, q, f, vertices_color, covered_intervals, uncovered_intervals, (1, 0), cs)
    elif mode == 'naive':
        cs = initial_cs(q, f, vertices_color)
        polytope = C_Polyhedron(cs)
    elif mode == 'sym_naive':
        # note that q = 2*f
        f2 = int(f / 2)
        vertices_color[1, f2] = vertices_color[q - f2, q - 1] = 0
        vertices_color[f2, f2] = vertices_color[q - f2, q - f2] = 0
        cs = initial_cs_sym(q, f, vertices_color)
        polytope = C_Polyhedron(cs)
    else:
        raise ValueError, "mode must be one of {'heuristic', 'combined', 'fulldim_covers', 'complete', 'naive', 'no_implied', 'sym', 'sym_naive'}."
    v_set = set([])
    if mode == 'combined' or mode == 'no_implied' or mode == 'sym':
        #global filename
        #filename = open(dir_math + "profiler/dim_threshold/%sslope_q%s_f%s_%sdim.txt" % (k_slopes, q, f, dim_threshold), "w")
        #print >> filename, "k_slope = %s, q = %s, f = %s, dim_threshold = %s" % (k_slopes, q, f, dim_threshold)
        for result_polytope, result_covered_intervals in gen:
            for values in generate_vertex_values(k_slopes, q, result_polytope, v_set):
                if all_intervals_covered(q, f, values, result_covered_intervals):
                    yield values
        #filename.close()
    elif mode == 'naive' or mode == 'sym_naive':
        for values in generate_vertex_values(k_slopes, q, polytope, v_set):
            if all_intervals_covered(q, f, values, []):
                yield values
    else:
        for result_polytope in gen:
            for values in generate_vertex_values(k_slopes, q, result_polytope, v_set):
                yield values
    if not v_set:
        logging.info("Example function not found. Please try again.")

import time
def search_kslope(k_slopes, q, f_list=None, mode='heuristic', print_function=False):
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
    start_cpu_t = time.clock()
    for f in f_list:
        cpu_t = time.clock()
        logging.info( "Search for extreme funtions with q = %s, f = %s, k_slopes >= %s, mode = %s" % (q, f, k_slopes, mode) )
        for h in search_kslope_example(k_slopes, q, f, mode):
            h_list.append(h)
            n_sol += 1
            if print_function:
                logging.info("h = %s" % (h,))
            logging.info("Found solution No.%s, in %s seconds" % (n_sol, time.clock() - cpu_t))
        logging.info("q = %s, f = %s takes cpu_time = %s" % (q, f, time.clock() - cpu_t))
    logging.info("Total cpu_time = %s" %(time.clock() - start_cpu_t))
    return h_list

def measure_stats(q, f_list, name=None):
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
        cs = initial_cs(q, f, vertices_color)
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

q_threshold = 20
dim_threshold = 8

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
                        yield polytope,  covered_intervals
                    else:
                        for result_polytope, result_covered_intervals in paint_complex_combined_pol(k_slopes, q, f, \
                                vertices_color, faces_color, covered_intervals, new_candidate_faces, polytope.constraints(), new_cs_matrix):
                            # Note: use minimized_constraints() in 'heuristic' mode takes longer. WHY??
                            yield result_polytope, result_covered_intervals
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
            # update constraint_system, set up MIP_Problem
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
                        yield polytope, covered_intervals
                    else:
                        for result_polytope, result_covered_intervals in paint_complex_combined_mip(k_slopes, q, f, \
                                vertices_color, faces_color, covered_intervals, new_candidate_faces, new_cs_matrix, sym=sym):
                            yield result_polytope, result_covered_intervals
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

dir_math ="/homes/home02/y/yzh/Dropbox/group-relaxation-sage-code/"
dir_yuan ="/media/sf_dropbox/sym_mode_2d_diagrams/"
def save_plot(q, hh, destdir = dir_math+"sym_mode_2d_diagrams/"):
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

def pattern_1(l, k_slopes, show_plots=False):
    f = int(18 * l + 11); 
    q = int(2 * f)
    f2 = int(9 * l + 5) # = int(f / 2)
    f3 = int(6 * l + 3) # = int(f / 3)
    vertices_color = initial_vertices_color(q, f)
    # impose some initial green vertices
    changed_vertices = [(f3 + 1, f3 + 1), (f2, f2)]
    for k in range(f3 + 2, f2, 3):
        changed_vertices += [(k, k), (k, k + 1), (k, k + 2), (k, k + 3), (k - 1, k + 2), (k + 1, k + 2), (k + 2, k + 2)]
    changed_vertices += [(1, f2), (2, f2 - 1), (2, f2), (3, f2 - 1)]
    for k in range(1, l + 1):
        i = 6 * k - 1
        j = f2 - 3 * k + 1
        changed_vertices += [(i - 1, j - 1), (i - 1, j), (i - 1, j + 1), (i - 1, j + 2), (i, j - 1), (i, j + 1), \
                             (i + 1, j - 2), (i + 1, j - 1), (i + 1, j), (i + 1, j + 1), \
                             (i + 2, j - 1), (i + 3, j - 2), (i + 3, j - 1), (i + 4, j - 2)]
    # impose their symmetric vertices too
    for (i, j) in changed_vertices:
        vertices_color[i, j] = vertices_color[q - j, q - i] = 0    
    #faces_color, covered_intervals = initial_faces_color_and_covered_intervals(q, f, vertices_color)
    cs = initial_cs_sym(q, f, vertices_color)
    polytope = C_Polyhedron(cs)
    v_set = set([])
    vv = []
    if show_plots == 'math':
        destdir = dir_math+"sym_mode_2d_diagrams/"+"patterns/"
    elif show_plots:
        destdir = dir_yuan+"patterns/"
    #logging.disable(logging.info)
    for v in polytope.minimized_generators():
        v_n = v.coefficients()
        num = len(set([v_n[i+1] - v_n[i] for i in range(q)]))
        if num >= k_slopes:
            if not tuple(v_n) in v_set:
                v_set.add(tuple(v_n))
                #if all_intervals_covered(q, f, v_n, covered_intervals):
                vv.append(v_n)
                if show_plots:
                        h = h_from_vertex_values(v_n)
                        name = "%sq%s_%s" %(num, q, len(vv))
                        #if not extremality_test(h):
                        #    name += "_notextreme"
                        #g2 =plot_with_colored_slopes(h);
                        #g2.save(destdir + name + "_plot.png" )
                        name += ".png"
                        g = plot_2d_diagram(h, colorful=True)
                        figsize = 12 * l + 8
                        g.save(destdir + name, figsize = figsize, show_legend=False)    
    #logging.disable(logging.NOTSET)
    return vv
