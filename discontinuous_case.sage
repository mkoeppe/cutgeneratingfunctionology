# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

########## Code for Discontinuous Case ###########

nonzero_eps = { (-1,-1,-1), (-1, 1,-1), (-1, 1, 1), (-1, 1, 0), (-1, 0,-1), ( 1,-1,-1), \
                ( 1,-1, 1), ( 1,-1, 0), ( 1, 1, 1), ( 1, 0, 1), ( 0,-1,-1), ( 0, 1, 1) }
continuous_xy_eps = { (-1,-1,-1), (1, 1, 1) }
type2_reduced_eps = { (0,-1,-1), (0, 1, 1), (1,-1,-1), (1, 1, 1), (-1,-1,-1), (-1, 1, 1), \
                      (1,-1, 0), (-1, 1, 0) }
dic_eps_to_cone = { (-1,-1,-1): [(-1, 0), (0, -1)], \
                    (-1, 1,-1): [(-1, 1), (-1, 0)], \
                    (-1, 1, 1): [(0, 1), (-1, 1)], \
                    (-1, 1, 0): [(-1, 1)], \
                    (-1, 0,-1): [(-1, 0)], \
                    ( 1,-1,-1): [(0, -1), (1, -1)], \
                    ( 1,-1, 1): [(1, -1), (1, 0)], \
                    ( 1,-1, 0): [(1, -1)], \
                    ( 1, 1, 1): [(1, 0), (0, 1)], \
                    ( 1, 0, 1): [(1, 0)], \
                    ( 0,-1,-1): [(0, -1)], \
                    ( 0, 1, 1): [(0, 1)], \
                    ( 0, 0, 0): [] \
                  }

def generate_type_1_vertices_general(fn, comparison, reduced=True, bkpt=None):
    """A generator of vertieces.

    ``..._general`` refers to the fact that it outputs 6-tuples (x,y,z,xeps,yeps,zeps).
    
    When reduced is:

    - ``True``: only outputs fewer triples satisfying ``comparison`` relation, for the purpose of ``minimality_test`` or setting up system of equations.
    - ``False``: outputs all triples satisfying ``comparison`` relation, for the purpose of plotting nonsubadditive or ``additive_limit_vertices``.
    """
    if bkpt is None:
        bkpt = fn.end_points()
    if not fn.is_continuous():
        #limits = fn.limits_at_end_points()
        limits = [fn.limits(x) for x in bkpt] # in case bkpt != fn.end_points()
    else:
        limits = [[fn(x)] for x in bkpt]
    for i in range(len(bkpt)):
        for j in range(i,len(bkpt)):
            x = bkpt[i]
            y = bkpt[j]
            z = fractional(x + y)
            if comparison(limits[i][0] + limits[j][0], fn(z)):
                yield (x, y, x+y, 0, 0, 0)
            if not fn.is_continuous():
                limits_x = limits[i]
                limits_y = limits[j]
                limits_z = fn.limits(z)
                if reduced and limits_x[0] == limits_x[1] == limits_x[-1] and limits_y[0] == limits_y[1] == limits_y[-1]:
                    eps_to_check = continuous_xy_eps # continuous at x and y
                else:
                    eps_to_check = nonzero_eps
                for (xeps, yeps, zeps) in eps_to_check:
                    if comparison(limits_x[xeps] + limits_y[yeps] - limits_z[zeps], 0):
                       yield (x, y, x+y, xeps, yeps, zeps)

def generate_type_2_vertices_general(fn, comparison, reduced=True, bkpt=None):
    """
    A generator of vertices.

    When reduced is:

    - ``True``: only outputs fewer triples satisfying ``comparison`` relation, for the purpose of ``minimality_test`` or setting up system of equations. Note: if `fn` is continuous at `y`, then use `fn(y^-) = fn(y) = fn(y^+)`.
    - ``False``: outputs all triples satisfying ``comparison`` relation, for the purpose of plotting nonsubadditive or ``additive_limit_vertices``.
    """
    if bkpt is None:
        bkpt = fn.end_points()
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]
    if not fn.is_continuous():
        #limits = fn.limits_at_end_points()
        limits = [fn.limits(x) for x in bkpt] # in case bkpt != fn.end_points()
    else:
        limits = [[fn(x)] for x in bkpt]
    for i in range(len(bkpt)):
        for k2 in range(i + 1, i + len(bkpt) - 1):
            # only need to check for 0 < y < 1. and note that bkpt2[i + len(bkpt) - 1] == bkpt[i] + 1.
            x = bkpt[i]
            z = bkpt2[k2]
            y = z - x
            if k2 < len(bkpt):
                k = k2
            else:
                k = k2 - len(bkpt) + 1
            if comparison(limits[i][0] + fn(y), limits[k][0]):
                yield (x, y, z, 0, 0, 0)
            if not fn.is_continuous():
                limits_x = limits[i]
                limits_z = limits[k]
                limits_y = fn.limits(y)
                # no trouble at 0- and 1+ since 0 < y < 1.
                if not (limits_y[0] == limits_y[1] == limits_y[-1]):
                    # then y is a in bkpt. this is done in type1check_general.
                    continue
                if reduced:
                    eps_to_check = type2_reduced_eps
                else:
                    eps_to_check = nonzero_eps
                for (xeps, yeps, zeps) in eps_to_check:
                    if comparison(limits_x[xeps] + limits_y[yeps] - limits_z[zeps], 0):
                       yield (x, y, x+y, xeps, yeps, zeps)

def generate_nonsymmetric_vertices_general(fn, f):
    """
    Generate vertices (x, y, xeps, yeps) that violate ``symmetric_test``.
    """
    bkpt = fn.end_points()
    limits = fn.limits_at_end_points()
    for i in range(len(bkpt)):
        x = bkpt[i]
        if x == f:
            continue
        if x < f:
            y = f - x
        else:
            y = 1 + f - x
        if limits[i][0] + fn(y) != 1:
            yield (x, y, 0, 0)
        if not fn.is_continuous():
            limits_x = limits[i]
            limits_y = fn.limits(y)
            if limits_x[-1] + limits_y[1] != 1:
                yield (x, y, -1, 1)
            if limits_x[1] + limits_y[-1] != 1:
                yield (x, y, 1, -1)

def epstriple_to_cone(epstriple):
    """
    Convert (xeps, yeps, zeps) to the corresponding cone.

    13 cases, see ``dic_eps_to_cone``.
    """
    try:
        return dic_eps_to_cone[epstriple]
    except KeyError:
        raise ValueError,"The limit epstriple %s does not exist." % epstriple

def plot_limit_cone_of_vertex(x, y, cone, color='red', r=0.03):
    """
    plot a cone or a ray or a point 
    """
    orig = vector(RDF, (x, y))
    if len(cone) == 0:
        p = point([orig], color=color, size=20, zorder=-1)
    elif len(cone) == 1:
        ray1 = vector(RDF, cone[0])
        p = line([orig, orig + ray1 * r / ray1.norm()], color=color, zorder=-3, thickness=3)
        p += point([orig], color='white', size=20, zorder=-2)
    elif len(cone) ==2:
        ray1 = vector(RDF, cone[0])
        ray2 = vector(RDF, cone[1])
        phi1 = arctan2(ray1[1], ray1[0])
        phi2 = arctan2(ray2[1], ray2[0])
        if phi1 > phi2:
            phi1, phi2 = phi2, phi1
        if phi2 - phi1 > pi:
            phi1, phi2 = phi2, phi1 + 2 * pi
        p = disk(orig, r, (phi1, phi2), color=color, zorder=-5)
        p += line([orig, orig + ray1 * r / ray1.norm()], color='white', zorder=-4, thickness=3)
        p += line([orig, orig + ray2 * r / ray2.norm()], color='white', zorder=-4, thickness=3)
    else:
        raise ValueError, "The cone %s is not well defined." % cone
    return p

def plot_2d_additive_limit_vertices(fn):
    """
    Temporary code. Show additivity information on a 2d-diagram
    """
    p = plot_2d_complex(fn)
    additive_vertices = generate_additive_vertices(fn, reduced=False)
    if fn.is_continuous():
        additive_vertices = {(x,y) for (x, y, z, xeps, yeps, zeps) in additive_vertices}
        p += point(list(additive_vertices), \
                            color = "mediumspringgreen", size = 50, legend_label="Additivity", zorder=-1)
        p += point([ (y,x) for (x,y) in additive_vertices ], color = "mediumspringgreen", size = 50, zorder=-1)
    elif additive_vertices != set([]):
        for (x, y, z, xeps, yeps, zeps) in additive_vertices:
            p += plot_limit_cone_of_vertex(x, y, epstriple_to_cone((xeps, yeps, zeps)), color="mediumspringgreen")
            if x != y:
                p += plot_limit_cone_of_vertex(y, x, epstriple_to_cone((yeps, xeps, zeps)), color="mediumspringgreen")
        # add legend_label
        p += point([(0,0)], color = "mediumspringgreen", size = 50, legend_label="Additivity", zorder=-10)
        p += point([(0,0)], color = "white", size = 50, zorder=-9)
    return p

def generate_symbolic_general(function, components, field=None, f=None):
    """
    Construct a vector-space-valued piecewise linear function
    compatible with the given function.  

    Each of the components of
    the function has a slope that is a basis vector of the vector
    space. Each discontinuous point has a left or/and right jump
    that is a basis vector of the vector space.
    """
    if f is None:
        f = find_f(function)
    n = len(components)
    intervals_and_slopes = []
    for component, slope in itertools.izip(components, range(n)):
        intervals_and_slopes.extend([ (interval, slope) for interval in component ])
    #intervals in components are coho intervals.
    #was: intervals_and_slopes.sort()
    intervals_and_slopes.sort(key=lambda (i, s): coho_interval_left_endpoint_with_epsilon(i))
    bkpt = [ field(interval[0]) for interval, slope in intervals_and_slopes ] + [field(1)]
    limits = [function.limits(x) for x in bkpt]
    if function.is_two_sided_discontinuous():
        num_jumps = len(bkpt) - 1
        num_left_jumps = bkpt.index(f)
    else:
        num_jumps = sum([x[0] != x[1] for x in limits[0:-1]])
        num_left_jumps = sum([(function.limit(x,-1) != function(x)) for x in bkpt if x > 0 and x <= f/2]) + \
                         sum([(function.limit(x,1) != function(x)) for x in bkpt if x < f/2])
    vector_space = VectorSpace(field, n + num_jumps)
    unit_vectors = vector_space.basis()
    slopes = [ unit_vectors[slope] for interval, slope in intervals_and_slopes ]
    jump_vectors = unit_vectors[n:n+num_left_jumps:] + unit_vectors[n+num_left_jumps-1:n-1:-1] + unit_vectors[n+num_left_jumps::] + unit_vectors[:n+num_left_jumps-1:-1]
    m = len(slopes)
    # Set up symbolic function
    current_value = zeros = vector_space.zero()
    pieces = []
    j = 0
    for i in range(m):
        pieces.append([singleton_interval(bkpt[i]), FastLinearFunction(zeros, current_value)])
        if function.is_two_sided_discontinuous() or limits[i][0] != limits[i][1]: # jump at bkpt[i]+
            current_value += jump_vectors[j]
            j += 1
        pieces.append([open_interval(bkpt[i], bkpt[i+1]), FastLinearFunction(slopes[i], current_value - slopes[i]*bkpt[i])])
        current_value += slopes[i] * (bkpt[i+1] - bkpt[i])
        if function.is_two_sided_discontinuous() or limits[i+1][-1] != limits[i+1][0]: # jump at bkpt[i+1]-
            current_value += jump_vectors[j]
            j += 1
    pieces.append([singleton_interval(bkpt[m]), FastLinearFunction(zeros, current_value)])
    symbolic_function = FastPiecewise(pieces, merge=True)
    if logging.getLogger().isEnabledFor(logging.DEBUG):
        logging.debug("Let v in R^%s.\nThe i-th entry of v represents the slope parameter on the i-th component of %s if i<=%s, or the function value jump parameter at breakpoint if i>%s. (The symmetry condition is considered so as to reduce the number of jump parameters).\nSet up the symbolic function sym: [0,1] -> R^%s, so that pert(x) = sym(x) * v.\nThe symbolic function sym is %s." % (n + num_jumps, components, n, n, n + num_jumps,  symbolic_function))
    return symbolic_function

def generate_additivity_equations_general(function, symbolic, field, f=None, bkpt=None):
    """
    Using additivity, set up a finite-dimensional system of linear equations
    that must be satisfied by any perturbation.
    """
    if f is None:
        f = find_f(function)
    vs = list(generate_additive_vertices(function, reduced = not function.is_two_sided_discontinuous(), bkpt=bkpt))
    equations = [symbolic(f), symbolic(field(1))]+[delta_pi_general(symbolic, x, y, (xeps, yeps, zeps)) for (x, y, z, xeps, yeps, zeps) in vs]
    M = matrix(field, equations)
    if not logging.getLogger().isEnabledFor(logging.DEBUG):
        return M
    pivot_r =  list(M.pivot_rows())
    for i in pivot_r:
        if i == 0:
            logging.debug("Condition pert(f) = 0 gives the equation\n%s * v = 0." % (symbolic(f)))
        elif i == 1:
            logging.debug("Condition pert(1) = 0 gives the equation\n%s * v = 0." % (symbolic(field(1))))
        else:
            (x, y, z, xeps, yeps, zeps) = vs[i-2]
            eqn = equations[i]
            logging.debug("Condition pert(%s%s) + pert(%s%s) = pert(%s%s) gives the equation\n%s * v = 0." % (x, print_sign(xeps),  y, print_sign(yeps), z, print_sign(zeps), eqn))
    return M[pivot_r]

def find_epsilon_interval_general(fn, perturb):
    """Compute the interval [minus_epsilon, plus_epsilon] such that
    (fn + epsilon * perturb) is subadditive for epsilon in this interval.
    Assumes that fn is subadditive.

    If one of the epsilons is 0, the function bails out early and returns 0, 0.
    """
    logging.info("Finding epsilon interval for perturbation...")
    fn_bkpt = fn.end_points()
    perturb_bkpt = perturb.end_points()
    bkpt = merge_bkpt(fn_bkpt,perturb_bkpt)
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]

    fn_limits = [fn.limits(x) for x in bkpt]
    perturb_limits = [perturb.limits(x) for x in bkpt]

    best_minus_epsilon_lower_bound = -10000
    best_plus_epsilon_upper_bound = +10000
    # type1check
    for i in range(len(bkpt)):
        perturb_x = perturb_limits[i]
        fn_x = fn_limits[i]
        for j in range(i,len(bkpt)):
            perturb_y = perturb_limits[j]
            fn_y = fn_limits[j]
            z = fractional(bkpt[i] + bkpt[j])
            perturb_z = perturb.limits(z)
            fn_z = fn.limits(z)
            if fn.is_continuous():
                eps_to_check = [(0, 0, 0)]
            elif fn_x[0] == fn_x[1] == fn_x[-1] and fn_y[0] == fn_y[1] == fn_y[-1] and \
                 perturb_x[0] == perturb_x[1] == perturb_x[-1] and perturb_y[0] == perturb_y[1] == perturb_y[-1]:
                ## if fn is continuous at x and y, then so is perturb.
                # both fn and perturb are continuous at x and y. ( needed if two-sided discontinuous at 0)
                eps_to_check = [(0, 0, 0)] + list(continuous_xy_eps)
            else:
                eps_to_check = [(0, 0, 0)] + list(nonzero_eps)
            for (xeps, yeps, zeps) in eps_to_check:
                delta_perturb = perturb_x[xeps] + perturb_y[yeps] - perturb_z[zeps]
                if delta_perturb != 0:
                    delta_fn = fn_x[xeps] + fn_y[yeps] - fn_z[zeps]
                    if delta_fn == 0:
                        logging.info("Zero epsilon encountered for x = %s%s, y = %s%s, z=%s%s" % (bkpt[i], print_sign(xeps), \
                                bkpt[j], print_sign(yeps), z, print_sign(zeps)) )
                        return 0, 0 # See docstring
                    epsilon_upper_bound = delta_fn / abs(delta_perturb)
                    if delta_perturb > 0:
                        if -epsilon_upper_bound > best_minus_epsilon_lower_bound:
                            best_minus_epsilon_lower_bound = -epsilon_upper_bound
                    else:
                        if epsilon_upper_bound < best_plus_epsilon_upper_bound:
                            best_plus_epsilon_upper_bound = epsilon_upper_bound
    # type2check
    for i in range(len(bkpt)):
        perturb_x = perturb_limits[i]
        fn_x = fn_limits[i]
        for k2 in range(i + 1, i + len(bkpt) - 1):
            if k2 < len(bkpt):
                k = k2
            else:
                k = k2 - len(bkpt) + 1
            perturb_z = perturb_limits[k]
            fn_z = fn_limits[k]
            y = bkpt2[k2] - bkpt[i]
            perturb_y = perturb.limits(y)
            fn_y = fn.limits(y)

            if fn.is_continuous():
                eps_to_check = [(0, 0, 0)]
            elif not (fn_y[0] == fn_y[1] == fn_y[-1]):
                # then y is a in bkpt. this is done in type1check.
                # for two_sided_discontinuous, could add (optional)
                # "or not (perturb_y[0] == perturb_y[1] == perturb_y[-1]):"
                eps_to_check = []
            else:
                # consider only y not being in bkpt.
                # so fn and perturb are both continuous at y. type2_reduced_eps works.
                eps_to_check = [(0, 0, 0)] + list(type2_reduced_eps)

            for (xeps, yeps, zeps) in eps_to_check:
                delta_perturb = perturb_x[xeps] + perturb_y[yeps] - perturb_z[zeps]
                if delta_perturb != 0:
                    delta_fn = fn_x[xeps] + fn_y[yeps] - fn_z[zeps]
                    if delta_fn == 0:
                        logging.info("Zero epsilon encountered for x = %s, y = %s" % (bkpt[i], y) )
                        return 0, 0 # See docstring
                    epsilon_upper_bound = delta_fn / abs(delta_perturb)
                    if delta_perturb > 0:
                        if -epsilon_upper_bound > best_minus_epsilon_lower_bound:
                            best_minus_epsilon_lower_bound = -epsilon_upper_bound
                    else:
                        if epsilon_upper_bound < best_plus_epsilon_upper_bound:
                            best_plus_epsilon_upper_bound = epsilon_upper_bound
    logging.info("Finding epsilon interval for perturbation... done.  Interval is %s", [best_minus_epsilon_lower_bound, best_plus_epsilon_upper_bound])
    return best_minus_epsilon_lower_bound, best_plus_epsilon_upper_bound

def delta_pi_general(fn, x, y, (xeps, yeps, zeps)=(0,0,0)):
    r"""
    return delta_pi = fn(x, xeps) + fn(y, yeps) - fn(z, zeps).

    EXAMPLES::

        sage: logging.disable(logging.INFO)
        sage: h = piecewise_function_from_breakpoints_and_limits(
        ....:       bkpt=[0, 1/5, 2/5, 3/5, 4/5, 1], 
        ....:       limits = [{-1:0, 0:0, 1:0},{-1:1, 0:1, 1:1}, 
        ....:                 {-1:0, 0:2/5, 1:2/5}, {-1:2/5, 0:1/2, 1:3/5},
        ....:                 {-1:3/5, 0:3/5, 1:1}, {-1:0, 0:0, 1:0}])
        sage: delta_pi_general(h, 2/5, 4/5, (-1, 1, -1))
        0
        sage: delta_pi_general(h, 2/5, 4/5, (-1, 0, -1))
        -2/5
    """
    return fn.limit(fractional(x), xeps) + fn.limit(fractional(y), yeps) - fn.limit(fractional(x + y), zeps)

def containing_eps_1d(x, interval):
    """
    Input:  

    - `x`: the projection of vertex `v` (of a face `F`), 

    - interval: the projection of the face `F`. 

    The projection direction is I/J/K. Note that `x` is in `interval`.

    Return: 

    - The projection of approching limits (`\subseteq \{x^-, x, x^+\}`) that need to be considered at `v` for testing the additivity of `F`.
    """
    if len(interval) == 1:
        return [0, 1, -1]
    elif x == interval[0]:
        return [1]
    elif x == interval[1]:
        return [-1]
    else:
        return [0]

def generate_containing_eps_triple(vertex, triple):
    """
    Given vertex `v` of face `F`, and the 3-projection-interval triple of `F`.
    Return the approching limits {(xeps, yeps, zeps)}
    pointing inwards at `v` from containning faces of `F`,
    that should be considered for testing the additivity of `F`.
    """
    xeps_list = containing_eps_1d(vertex[0], triple[0])
    yeps_list = containing_eps_1d(vertex[1], triple[1])
    zeps_list = containing_eps_1d(vertex[0] + vertex[1], triple[2])
    return [(xeps, yeps, zeps) for xeps in xeps_list for yeps in yeps_list for zeps in zeps_list]

def is_additive_face(fn, face):
    """
    Test whether the given face is additive 
    by taking the appropriate limits (pointing inwards) at the vertices.
    """
    if face.is_2D():
        for vertex in face.vertices:
            for eps_triple in generate_containing_eps_triple(vertex, face.minimal_triple):
                if delta_pi_general(fn, vertex[0], vertex[1], eps_triple) != 0:
                    return False
        return True
    elif face.is_1D():
        vertex_0 = face.vertices[0]
        vertex_1 = face.vertices[1]
        eps_triple_0 = generate_containing_eps_triple(vertex_0, face.minimal_triple)
        eps_triple_1 = generate_containing_eps_triple(vertex_1, face.minimal_triple)
        # FIXME: both eps_triple_0 and _1 have length 3? in compatible order? Yes.
        for i in range(3):
            if delta_pi_general(fn, vertex_0[0], vertex_0[1], eps_triple_0[i]) == 0 and \
               delta_pi_general(fn, vertex_1[0], vertex_1[1], eps_triple_1[i]) == 0:
                return True
        return False
    else:
        vertex = face.vertices[0]
        for eps_triple in nonzero_eps:
            if delta_pi_general(fn, vertex[0], vertex[1], eps_triple) == 0:
                return False
        return delta_pi_general(fn, vertex[0], vertex[1], (0,0,0)) == 0

def x_y_swapped_face(face):
    vert = face.vertices
    vert_sym = [(vertex[1], vertex[0]) for vertex in vert]
    trip = face.minimal_triple
    return Face( (trip[1], trip[0], trip[2]), vertices=vert_sym )

def generate_maximal_additive_faces_general(function):
    logging.info("Computing maximal additive faces...")
    bkpt = function.end_points()
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]
    n = len(bkpt) - 1
    I_list = J_list = [ (bkpt[i], bkpt[i+1]) for i in range(n) ]
    K_list = [ (bkpt2[i], bkpt2[i+1]) for i in range(2*n) ]

    faces = []
    # 2D faces
    for i in range(n):
        for j in range(i, n):
            IplusJ = interval_sum(I_list[i],J_list[j])
            for k in generate_overlapping_interval_indices(IplusJ, bkpt2):
                # Check if int(I+J) intersects int(K) is non-empty.
                if len(interval_intersection(IplusJ,K_list[k])) == 2:
                    face = Face( (I_list[i], J_list[j], K_list[k]) )
                    if is_additive_face(function, face): 
                        faces.append(face)
                        if i != j:
                            faces.append(x_y_swapped_face(face))
    # 1D horizontal and vertical faces
    for i in range(n):
        for j in range(n):
            I = [bkpt[i]]
            IplusJ = (bkpt[i] + bkpt[j], bkpt[i] + bkpt[j+1])
            for k in generate_overlapping_interval_indices(IplusJ, bkpt2):
                if len(interval_intersection(IplusJ, K_list[k])) == 2:
                    face = Face( (I, J_list[j], K_list[k]) )
                    if is_additive_face(function, face): 
                        faces.append(face)
                        faces.append(x_y_swapped_face(face))
    # 1D diagonal faces
    for i in range(n):
        for j in range(i, n):
            interval_K = interval_sum(I_list[i],J_list[j])
            for k in generate_overlapping_interval_indices(interval_K, bkpt2):
                if interval_K[0] < bkpt2[k] < interval_K[1]:
                    face = Face( (I_list[i], J_list[j], [bkpt2[k]]) )
                    if is_additive_face(function, face): 
                        faces.append(face)
                        if i != j:
                            faces.append(x_y_swapped_face(face))
    # 0D faces
    additive_vertices = {(x,y) for (x, y, z, xeps, yeps, zeps) in generate_additive_vertices(function) if x != 1 and y != 1}
    additive_vertices_seen = {vertex for face in faces for vertex in face.vertices}
    additive_vertices_new = additive_vertices.difference(additive_vertices_seen)
    for (x, y) in additive_vertices_new:
        face = Face(([x], [y], [x+y]))
        faces.append(face)
        if x != y:
            face = Face(([y], [x], [x+y]))
            faces.append(face)

    logging.info("Computing maximal additive faces... done")
    return faces
