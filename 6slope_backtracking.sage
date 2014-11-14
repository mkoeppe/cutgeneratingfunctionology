# backtracking search for 6-slope extreme functions, 
# using incremental computation for vertices_color, additive_faces and covered_intervals.
# do not consider translation/reflection other than the symmetry reflection regarding f. (so, no edge-merge).
# incremental polyhedral computation for implied green vertices,
# using http://www.sagemath.org/doc/reference/libs/sage/libs/ppl.html#sage.libs.ppl.Polyhedron.minimize

# Note: vertices_color (q+1)*(q+1) 0-1 array, 0: green, 1: white.
# vertices have integer coordinates, k is integer in function value fn(k).
# faces, covered_intervals have rational coordinates.
# only record vertex (x,y) with x <=y and F(I,J,K) with I <= J.
# polytope defines the feasible region of (fn(0), fn(1/q),...,fn(1)).

from sage.libs.ppl import C_Polyhedron, Constraint, Constraint_System, Generator, Generator_System, Variable, point, Poly_Con_Relation
import numpy

#global num_of_slopes
num_of_slopes = 6 # set to 6 as we are looking for 6-slope functions

def initial_vertices_color(q, f):
    """
    paint green ( = 0) for vertices (x<=y) corresonding to
    fn(0) = 0, fn(1) = 0 and reflection around f/2 and (1-f)/2.

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
    ff = int(f*q)
    # border x = 0 and border y = 0 are green
    for x in range(q+1):
        vertices_color[0, x] = 0
        vertices_color[x, q] = 0
    # diagonals corresponding to f
    for x in range(q+1):
        if x <= ff/2:
            vertices_color[x, ff - x] = 0
        elif (x >= ff) and (x <= ff - x + q):
            vertices_color[x, ff - x + q] = 0
    return vertices_color

def initial_faces_and_covered_intervals(q, f, vertices_color):
    """
    Return additive_faces_set and covered_intervals,
    corresponding to fn(0) = 0, fn(1) = 0 and reflection around f/2 and (1-f)/2.

    EXAMPLES::

        sage: q=5; f=3/5;
        sage: vertices_color = initial_vertices_color(q, f);
        sage: faces_set, covered_intervals = initial_faces_and_covered_intervals(q, f, vertices_color)
        sage: faces_set
        set([<Face ([3/5, 4/5], [4/5, 1], [8/5, 9/5])>, <Face ([4/5, 1], [4/5, 1], [9/5, 2])>, 
             <Face ([0, 1/5], [0, 1/5], [0, 1/5])>, <Face ([4/5, 1], [4/5, 1], [8/5, 9/5])>, 
             <Face ([0, 1/5], [2/5, 3/5], [2/5, 3/5])>])
        sage: covered_intervals
        [[[0, 1/5], [2/5, 3/5]], [[3/5, 4/5], [4/5, 1]]]
    """
    faces_set = set([])
    covered_intervals = []
    for x in range(q):
        for y in range(x, q):
            for z in range(2):
                face = Face(([x/q, (x+1)/q], [y/q, (y+1)/q], [(x+y+z)/q, (x+y+z+1)/q]))
                if x < y:
                    vertices = [(x+1, y), (x, y+1), (x+z, y+z)]
                else:
                    vertices = [(x, y+1), (x+z, y+z)]
                if sum(vertices_color[v] for v in vertices) == 0 :
                    faces_set.add(face)
                    covered_intervals = directly_covered_by_adding_face(covered_intervals, face, f)
    return faces_set, covered_intervals

def initial_cs(q, f, vertices_color):
    """
    Return the initial Constraint_System that defines 
    the feasible region of (fn(0), fn(1/q),...,fn(1))

    EXAMPLES::

        sage: q=5; f=3/5;
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
    cs.insert(fn[int(f*q)] == 1)
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
    Given a grid vertex v (integers, v[0] <= v[1]), return small triangle faces (QQ, I <= J)
    and their vertices (integers, xx <= yy) that are around v.

    EXAMPLES::

        sage: faces_around_vertex(5, (1, 2))
        [(<Face ([0, 1/5], [2/5, 3/5], [2/5, 3/5])>, [(0, 2), (1, 2), (0, 3)]), 
         (<Face ([0, 1/5], [2/5, 3/5], [3/5, 4/5])>, [(1, 2), (0, 3), (1, 3)]), 
         (<Face ([1/5, 2/5], [2/5, 3/5], [3/5, 4/5])>, [(1, 2), (2, 2), (1, 3)]), 
         (<Face ([0, 1/5], [1/5, 2/5], [2/5, 3/5])>, [(1, 1), (0, 2), (1, 2)]), 
         (<Face ([1/5, 2/5], [1/5, 2/5], [2/5, 3/5])>, [(1, 1), (1, 2)]), 
         (<Face ([1/5, 2/5], [1/5, 2/5], [3/5, 4/5])>, [(1, 2), (2, 2)])]
        sage: faces_around_vertex(5, (1, 1))
        [(<Face ([0, 1/5], [1/5, 2/5], [1/5, 2/5])>, [(0, 1), (1, 1), (0, 2)]), 
         (<Face ([0, 1/5], [1/5, 2/5], [2/5, 3/5])>, [(1, 1), (0, 2), (1, 2)]), 
         (<Face ([1/5, 2/5], [1/5, 2/5], [2/5, 3/5])>, [(1, 1), (1, 2)]), 
         (<Face ([0, 1/5], [0, 1/5], [1/5, 2/5])>, [(0, 1), (1, 1)])]
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
           (Face(([xl/q, xx/q],[yy/q, yr/q],[(xl+yy)/q, (xl+yy+1)/q])), [(xl, yy), (xx, yy), (xl, yr)]), \
           (Face(([xl/q, xx/q],[yy/q, yr/q],[(xl+yy+1)/q, (xl+yy+2)/q])), [(xx, yy), (xl, yr), (xx, yr)]), \
           (Face(([xx/q, xr/q],[yy/q, yr/q],[(xx+yy)/q, (xx+yy+1)/q])), [(xx, yy), (xr, yy), (xx, yr)]), \
           (Face(([xl/q, xx/q],[yl/q, yy/q],[(xl+yl+1)/q, (xl+yl+2)/q])), [(xx, yl), (xl, yy), (xx, yy)]),\
           (Face(([xx/q, xr/q],[yl/q, yy/q],[(xx+yl)/q, (xl+yy+1)/q])), [(xx, yl), (xr, yl), (xx, yy)]), \
           (Face(([xx/q, xr/q],[yl/q, yy/q],[(xx+yl+1)/q, (xl+yy+2)/q])), [(xx, yy),(xr, yl),(xr, yy)]) ]
    elif xx == yl:
        return [ \
           (Face(([xl/q, xx/q],[yy/q, yr/q],[(xl+yy)/q, (xl+yy+1)/q])), [(xl, yy), (xx, yy), (xl, yr)]), \
           (Face(([xl/q, xx/q],[yy/q, yr/q],[(xl+yy+1)/q, (xl+yy+2)/q])), [(xx, yy), (xl, yr), (xx, yr)]), \
           (Face(([xx/q, xr/q],[yy/q, yr/q],[(xx+yy)/q, (xx+yy+1)/q])), [(xx, yy), (xr, yy), (xx, yr)]), \
           (Face(([xl/q, xx/q],[yl/q, yy/q],[(xl+yl+1)/q, (xl+yl+2)/q])), [(xx, yl), (xl, yy), (xx, yy)]),\
           (Face(([xx/q, xr/q],[yl/q, yy/q],[(xx+yl)/q, (xl+yy+1)/q])), [(xx, yl), (xx, yy)]), \
           (Face(([xx/q, xr/q],[yl/q, yy/q],[(xx+yl+1)/q, (xl+yy+2)/q])), [(xx, yy),(xr, yy)]) ]
    else:
        return [ \
           (Face(([xl/q, xx/q],[yy/q, yr/q],[(xl+yy)/q, (xl+yy+1)/q])), [(xl, yy), (xx, yy), (xl, yr)]), \
           (Face(([xl/q, xx/q],[yy/q, yr/q],[(xl+yy+1)/q, (xl+yy+2)/q])), [(xx, yy), (xl, yr), (xx, yr)]), \
           (Face(([xx/q, xr/q],[yy/q, yr/q],[(xx+yy)/q, (xx+yy+1)/q])), [(xx, yy), (xx, yr)]), \
           (Face(([xl/q, xx/q],[yl/q, yy/q],[(xl+yl+1)/q, (xl+yl+2)/q])), [(xl, yy), (xx, yy)])]

def directly_covered_by_adding_face(last_covered_intervals, face, f):
    """
    Compute incrementally new covered_intervals by adding a new face.
    Consider only directly covered and symmetry reflection regarding f.

    EXAMPLES::

        sage: last_covered_intervals = [[[0, 1/5], [2/5, 3/5]], [[3/5, 4/5], [4/5, 1]]]
        sage: face = Face(([1/5, 2/5], [1/5, 2/5], [3/5, 4,5]))
        sage: directly_covered_by_adding_face(last_covered_intervals, face, 3/5)
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

def paint_complex(q, f, last_vertices_color, last_faces_set, last_covered_intervals, \
                  candidate_faces, last_non_candidate, last_cs):
    """
    Randomly paint triangles green in a 2d-complex, until all intervals are covered.
    Return the polytope which defines the feasible region of (fn(0), fn(1/q),...,fn(1)) for that paint_complex

    EXAMPLES::
        sage: q = 5; f = 3/5; num_of_slopes = 2;
        sage: vertices_color = initial_vertices_color(q, f);
        sage: faces_set, covered_intervals = initial_faces_and_covered_intervals(q, f, vertices_color)
        sage: cs = initial_cs(q, f, vertices_color)
        sage: candidate_faces = generate_candidate_faces(q, f, covered_intervals, None)
        sage: for result_polytope in paint_complex(q, f, vertices_color, faces_set, covered_intervals, candidate_faces, set([]), cs):
        ...       print result_polytope.minimized_generators()
        Generator_System {point(0/6, 2/6, 4/6, 6/6, 3/6, 0/6)}
        Generator_System {point(0/4, 3/4, 1/4, 4/4, 2/4, 0/4)}
    """
    non_candidate = copy(last_non_candidate)
    num_candidate_faces = len(candidate_faces)
    n = 0
    while n < num_candidate_faces:
        face_picked = candidate_faces[n]
        n += 1
        faces_set = copy(last_faces_set)
        faces_set.add(face_picked)
        vertices_color = copy(last_vertices_color)
        covered_intervals = directly_covered_by_adding_face(last_covered_intervals, face_picked, f)
        cs = copy(last_cs)

        legal_picked = True
        for v_green in face_picked.vertices:
            x = int(v_green[0] * q)
            y = int(v_green[1] * q)
            if (x <= y) and (vertices_color[x, y] == 1) and legal_picked:
                # new green vertice
                vertices_color[x, y] = 0
                # update cs
                z = (x + y) % q
                cs.insert( Variable(x) + Variable(y) == Variable(z) )
                for (face, vertices) in faces_around_vertex(q, (x, y)):
                    if not face in faces_set and sum(vertices_color[v] for v in vertices) == 0:
                        # find new green face.
                        if face in non_candidate:
                            legal_picked = False
                            break
                        else:
                            faces_set.add(face)
                            covered_intervals = directly_covered_by_adding_face(covered_intervals, face, f)
        polytope = C_Polyhedron(cs)
        if polytope.is_empty() or not legal_picked or num_slopes_at_best(q, covered_intervals) < num_of_slopes:
            # If infeasible or encounter non_candidate or too few slopes, stop recursion.
            non_candidate.add(face_picked)
            continue
        # look for implied additive vertices
        for x in range(1, q):
            for y in range(x, q):
                if legal_picked and (vertices_color[x, y] == 1):
                    z = (x + y) % q
                    # TODO: compare "maximize" and "relation_with"
                    # try "maximize"
                    #if polytope.maximize(Variable(x) + Variable(y) - Variable(z))['sup_n'] == 0:
                    # try "relation_with"
                    if polytope.relation_with( Variable(x) + Variable(y) == Variable(z) ).implies(Poly_Con_Relation.is_included()):
                        # find implied additive vertices
                        vertices_color[x, y] = 0
                        polytope.add_constraint( Variable(x) + Variable(y) == Variable(z) ) #TODO: try without adding implied equality
                        for (face, vertices) in faces_around_vertex(q, (x, y)):
                            if not face in faces_set and sum(vertices_color[v] for v in vertices) == 0:
                                # find new green face.
                                if face in non_candidate:
                                    legal_picked = False
                                    break
                                else:
                                    faces_set.add(face)
                                    covered_intervals = directly_covered_by_adding_face(covered_intervals, face, f)
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
                for result_polytope in paint_complex(q, f, vertices_color, faces_set, \
                        covered_intervals, new_candidate_faces, non_candidate, polytope.minimized_constraints()): #TODO: try with polytope.constraints()
                    yield result_polytope
        non_candidate.add(face_picked)

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

def generate_candidate_faces(q, f, covered_intervals, last_face=None):
    """
    Return a set of candidate_faces (lexicographically > last_face, not in non_candidate)
    to paint in next step, whose I, J are currently uncovered.
    Note that candidate_faces only takes faces with I <= J.

    EXAMPLES::

        sage: q = 5; f = 3/5; num_of_slopes = 2;
        sage: vertices_color = initial_vertices_color(q, f);
        sage: faces_set, covered_intervals = initial_faces_and_covered_intervals(q, f, vertices_color)
        sage: generate_candidate_faces(q, f, covered_intervals, None)
        set([<Face ([1/5, 2/5], [1/5, 2/5], [2/5, 3/5])>, <Face ([1/5, 2/5], [1/5, 2/5], [3/5, 4/5])>])
    """
    to_cover = generate_to_cover(q, covered_intervals)
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

def generate_to_cover(q, covered_intervals):
    """
    Return a set {k/q | 0 <=k < q, [k/q, (k+1)/q] is uncovered}

    EXAMPLES::

        sage: covered_intervals = [[[0, 1/5], [2/5, 3/5]], [[3/5, 4/5], [4/5, 1]]]
        sage: generate_to_cover(5, covered_intervals)
        [1/5]
    """
    to_cover = set([x/q for x in range(q)])
    for component in covered_intervals:
        for i in component:
            for x in range(i[0]*q, i[1]*q):
                to_cover.discard(x/q)
    return sorted(list(to_cover))

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

def generate_vertex_function(q, polytope):
    """
    Generate real valued functions corresponding to vertices of the polytope.
    Return those functions that have required number of slope values.

    EXAMPLES::
        sage: q=5; f=3/5; num_of_slopes = 2;
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
        sage: q=5; f=3/5; num_of_slopes = 2;
        sage: h = search_6slope_example(q, f).next()
    """
    #initialization
    vertices_color = initial_vertices_color(q, f)
    faces_set, covered_intervals = initial_faces_and_covered_intervals(q, f, vertices_color)
    cs = initial_cs(q, f, vertices_color)
    candidate_faces = generate_candidate_faces(q, f, covered_intervals, None)

    found = False
    for result_polytope in \
            paint_complex(q, f, vertices_color, faces_set, covered_intervals, candidate_faces, set([]), cs):
        for h in generate_vertex_function(q, result_polytope):
            print h
            found = True
            yield h #return h
    if not found:
        print "Example function not found. Please try again."
