from six.moves import range
from six.moves import zip
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
    r"""A generator of vertices.

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
            limits_x = limits[i]
            limits_y = limits[j]
            limits_z = fn.limits(z)
            def check(xeps, yeps, zeps):
                if all(l is not None for l in (limits_x[xeps], limits_y[yeps], limits_z[zeps])):
                    if comparison(limits_x[xeps] + limits_y[yeps] - limits_z[zeps], 0):
                        return True
                return False
            if check(0, 0, 0):
                yield (x, y, x+y, 0, 0, 0)
            if not fn.is_continuous():
                if reduced and limits_x[0] == limits_x[1] == limits_x[-1] and limits_y[0] == limits_y[1] == limits_y[-1]:
                    eps_to_check = continuous_xy_eps # continuous at x and y
                else:
                    eps_to_check = nonzero_eps
                for (xeps, yeps, zeps) in eps_to_check:
                    if check(xeps, yeps, zeps):
                        yield (x, y, x+y, xeps, yeps, zeps)

def generate_type_2_vertices_general(fn, comparison, reduced=True, bkpt=None):
    r"""
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
            limits_x = limits[i]
            limits_z = limits[k]
            limits_y = fn.limits(y)
            def check(xeps, yeps, zeps):
                if all(l is not None for l in (limits_x[xeps], limits_y[yeps], limits_z[zeps])):
                    if comparison(limits_x[xeps] + limits_y[yeps] - limits_z[zeps], 0):
                        return True
                return False
            if check(0, 0, 0):
                yield (x, y, z, 0, 0, 0)
            if not fn.is_continuous():
                # no trouble at 0- and 1+ since 0 < y < 1.
                if not (limits_y[0] == limits_y[1] == limits_y[-1]):
                    # then y is a in bkpt. this is done in type1check_general.
                    continue
                if reduced:
                    eps_to_check = type2_reduced_eps
                else:
                    eps_to_check = nonzero_eps
                for (xeps, yeps, zeps) in eps_to_check:
                    if check(xeps, yeps, zeps):
                       yield (x, y, x+y, xeps, yeps, zeps)

def generate_nonsymmetric_vertices_general(fn, f):
    r"""
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
    r"""
    Convert (xeps, yeps, zeps) to the corresponding cone.

    13 cases, see ``dic_eps_to_cone``.
    """
    try:
        return dic_eps_to_cone[epstriple]
    except KeyError:
        raise ValueError("The limit epstriple %s does not exist." % epstriple)

plot_limit_cone_style = 'sectors'
plot_limit_cone_wrap = True
plot_limit_cone_arrow_rel_distance = 0.008 / 0.055
plot_limit_cone_arrow_length = 0.055
plot_limit_cone_width = 1
plot_limit_cone_vertex_size = 20

from cutgeneratingfunctionology.igp.subadditivity_slack_diagrams.limit_arrow import limit_arrow

def plot_limit_cone_of_vertex_in_face(x, y, F, color='red', **kwds):
    v = vector([x, y])
    feasible_cone = Polyhedron(rays=[ vector(u) - v for u in F.vertices if vector(u) != v])
    rays = feasible_cone.rays_list()
    return plot_limit_cone_of_vertex(x, y, rays, color=color, **kwds)

def plot_limit_cone_of_vertex(x, y, cone, color='red', r=0.03, style=None, zorder=None,
                              width=None, vertex_size=None, **kwds):
    r"""
    plot a cone or a ray or a point 
    """
    orig = vector(RDF, (x, y))
    if style is None:
        style = plot_limit_cone_style
    if style == 'arrows':
        if color == 'white':
            return Graphics()
        if zorder is None:
            zorder = 10
        if width is None:
            width = plot_limit_cone_width
        if vertex_size is None:
            vertex_size = plot_limit_cone_vertex_size
        if len(cone) == 0:
            return point([orig], color=color, size=plot_limit_cone_vertex_size,
                         zorder=zorder, **kwds)  # on top of the complex
        else:
            uv = sum(vector(ray) for ray in cone)
            if plot_limit_cone_wrap:
                def wrap(coord, direction):
                    if coord == 0 and direction < 0:
                        return 1
                    if coord == 1 and direction > 0:
                        return 0
                    return coord
                orig = vector(RDF, (wrap(s, t) for s, t in zip(orig, uv)))
            uv = vector(RDF, uv)
            arrowsize = 3
            if uv.norm() > plot_limit_cone_arrow_length:
                uv = uv / uv.norm() * plot_limit_cone_arrow_length
            elif uv.norm() < 0.33 * plot_limit_cone_arrow_length:
                uv = uv / uv.norm() * 0.33 * plot_limit_cone_arrow_length
                arrowsize = 1.5
            return limit_arrow(orig + (1 + plot_limit_cone_arrow_rel_distance) * uv,
                               orig + plot_limit_cone_arrow_rel_distance * uv, color=color, arrowsize=arrowsize,
                               zorder=zorder, width=width)
    # default: 'sectors'.  FIXME: Ignores zorder.
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
        raise ValueError("The cone %s is not well defined." % cone)
    return p

def plot_2d_additive_limit_vertices(fn):
    r"""
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

generate_symbolic_two_sided_discontinuous_basis_functions = ('slopes', 'jumps')  # default is classic

def generate_symbolic_general(function, components, field=None, f=None, basis_functions=None):
    r"""
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
    for component, slope in zip(components, list(range(n))):
        intervals_and_slopes.extend([ (interval, slope) for interval in component ])
    #intervals in components are coho intervals.
    #was: intervals_and_slopes.sort()
    intervals_and_slopes.sort(key=lambda i_s: coho_interval_left_endpoint_with_epsilon(i_s[0]))
    basis = []

    if basis_functions is None:
        if function.is_two_sided_discontinuous():
            basis_functions = generate_symbolic_two_sided_discontinuous_basis_functions
        else:
            basis_functions = ('slopes', 'jumps')

    if basis_functions == ('midpoints', 'slopes'):
        # Use slopes and values at midpoints of intervals (including singletons).
        # Because in the two-sided discontinuous case we have jumps at every breakpoint,
        # there are no relations from continuity. 
        # In this basis, additivities are expressible with equations of small support.
        # Dictionary mapping an interval (a, b) to a list of pairs "(index, coeff)"
        # suitable for VectorSpace.sum_of_terms.
        midpoint_value_dict = { x: []                        # Already fixed by symmetry
                                for x in [0, f/2, f, (1+f)/2, 1] }
        def midpoint_value_lincomb(a, b):
            "Allocate a variable to be the midpoint value on the interval."
            interval = (a, b)
            midpoint = fractional((a + b) / 2)
            if midpoint not in midpoint_value_dict:
                v = [(len(basis), 1)]     # (index, coefficient)
                midpoint_value_dict[midpoint] = v
                midpoint_value_dict[fractional(f - midpoint)] = [(len(basis), -1)]
                basis.append(('function value at', midpoint))
            return midpoint_value_dict[midpoint]
        for interval, slope in intervals_and_slopes:
            a = interval[0]
            b = interval[1]
            assert a < b   # no singletons in components.
            midpoint_value_lincomb(a, a)   # allocates
            midpoint_value_lincomb(a, b)   # allocates
            midpoint_value_lincomb(b, b)   # allocates
        slope_dict = {}
        # Put slopes at the end.  They are dense columns.
        for slope_index in range(n):
            slope_dict[slope_index] = [(len(basis), 1)]
            basis.append(('slope of component', slope_index))
        dimension = len(basis)
    elif basis_functions == ('slopes', 'jumps'):
        bkpt = [ field(interval[0]) for interval, slope in intervals_and_slopes ] + [field(1)]
        limits = [function.limits(x) for x in bkpt]
        if function.is_two_sided_discontinuous():
            num_jumps = len(bkpt) - 1
            num_left_jumps = bkpt.index(f)   # this is assuming the interval right of f is in a component.
        else:
            num_jumps = sum([x[0] != x[1] for x in limits[0:-1]])
            num_left_jumps = sum([(function.limit(x,-1) != function(x)) for x in bkpt if x > 0 and x <= f/2]) + \
                             sum([(function.limit(x,1) != function(x)) for x in bkpt if x < f/2])
        dimension = n + num_jumps
    else:
        raise ValueError("this type of basis_functions is not supported")

    #import pdb; pdb.set_trace()
    vector_space = VectorSpace(field, dimension)
    pieces = []

    if basis_functions == ('midpoints', 'slopes'):
        # Now we have the space.
        def midpoint_value(a, b):
            return vector_space.sum_of_terms(midpoint_value_lincomb(a, b))
        def slope(slope_index):
            return vector_space.sum_of_terms(slope_dict[slope_index])
        last_x = None
        for interval, slope_index in intervals_and_slopes:
            a = interval[0]
            b = interval[1]
            assert a < b   # no singletons in components.
            if last_x is None or last_x < a:
                pieces.append((singleton_interval(a),
                               FastLinearFunction(vector_space.zero(), midpoint_value(a, a))))
            midpoint = (a + b) / 2
            # the function is: slope * (x - midpoint) + value_at_midpoint = slope * x + value_at_midpoint - slope * midpoint.
            pieces.append((open_interval(a, b),
                           FastLinearFunction(slope(slope_index),
                                              midpoint_value(a, b) - slope(slope_index) * midpoint)))
            pieces.append((singleton_interval(b),
                           FastLinearFunction(vector_space.zero(), midpoint_value(b, b))))
            last_x = b
        if last_x < 1:
            pieces.append((singleton_interval(1),
                           FastLinearFunction(vector_space.zero(), midpoint_value(1, 1))))
        # FIXME: Add logging.debug with details on the variables.
    else:
        # Use slopes and jumps from left to right.
        # This uses known continuity to reduce the number of variables.
        # Warning: These basis functions are monotonically increasing and therefore NOT periodic.
        unit_vectors = vector_space.basis()
        slopes = [ unit_vectors[slope] for interval, slope in intervals_and_slopes ]
        m = len(slopes)
        for slope_index in range(n):
            basis.append(('slope of component', slope_index))
        jumps = list(range(len(basis), len(basis) + num_jumps))
        basis.extend([None] * num_jumps)
        jump_coordinates = jumps[0:num_left_jumps:] + jumps[num_left_jumps-1::-1] + jumps[num_left_jumps::] + jumps[:num_left_jumps-1:-1]
        assert len(jump_coordinates) == 2 * num_jumps
        # Set up symbolic function
        current_value = zeros = vector_space.zero()
        j = 0
        for i in range(m):
            pieces.append([singleton_interval(bkpt[i]), FastLinearFunction(zeros, current_value)])
            if function.is_two_sided_discontinuous() or limits[i][0] != limits[i][1]: # jump at bkpt[i]+
                current_value += unit_vectors[jump_coordinates[j]]
                if basis[jump_coordinates[j]] is None:
                    basis[jump_coordinates[j]] = ('jump', bkpt[i], +1)
                j += 1
            pieces.append([open_interval(bkpt[i], bkpt[i+1]), FastLinearFunction(slopes[i], current_value - slopes[i]*bkpt[i])])
            current_value += slopes[i] * (bkpt[i+1] - bkpt[i])
            if function.is_two_sided_discontinuous() or limits[i+1][-1] != limits[i+1][0]: # jump at bkpt[i+1]-
                current_value += unit_vectors[jump_coordinates[j]]
                if basis[jump_coordinates[j]] is None:
                    basis[jump_coordinates[j]] = ('jump', bkpt[i+1], -1)
                j += 1
        assert j == 2 * num_jumps
        pieces.append([singleton_interval(bkpt[m]), FastLinearFunction(zeros, current_value)])
        if logging.getLogger().isEnabledFor(logging.DEBUG):
            logging.debug("Let v in R^%s.\nThe i-th entry of v represents the slope parameter on the i-th component of %s if i<=%s, or the function value jump parameter at breakpoint if i>%s. (The symmetry condition is considered so as to reduce the number of jump parameters).\n" % (n + num_jumps, components, n, n))

    symbolic_function = FastPiecewise(pieces, merge=True)
    symbolic_function.basis = basis
    if logging.getLogger().isEnabledFor(logging.DEBUG):
        logging.debug("Set up the symbolic function sym: [0,1] -> R^%s, so that pert(x) = sym(x) * v.\nThe symbolic function sym is %s." % (dimension, symbolic_function))
    return symbolic_function

def generate_additivity_equations_general(function, symbolic, field, f=None, bkpt=None,
                                          reduce_system=None, return_vertices=False,
                                          undefined_ok=False, vertices=None):
    r"""
    Using additivity, set up a finite-dimensional system of linear equations
    that must be satisfied by any perturbation.
    """
    if vertices is None:
        vertices = generate_additive_vertices(function, reduced = not function.is_two_sided_discontinuous(), bkpt=bkpt)
    vs = list(vertices)
    def generate_labeled_equations():
        yield 'f', symbolic(f)
        yield '1', symbolic(field(1))
        for (x, y, z, xeps, yeps, zeps) in vs:
            try:
                deltafn = delta_pi_general(symbolic, x, y, (xeps, yeps, zeps))
            except ValueError:
                if undefined_ok:
                    pass
                else:
                    raise
            yield (x, y, z, xeps, yeps, zeps), deltafn

    vs, equations = zip(*generate_labeled_equations())
    M = matrix(field, equations)
    if reduce_system is None:
        reduce_system = logging.getLogger().isEnabledFor(logging.DEBUG)
    if not reduce_system:
        if return_vertices:
            return M, vs
        else:
            return M
    pivot_r =  list(M.pivot_rows())
    for i in pivot_r:
        if vs[i] == 'f':
            logging.debug("Condition pert(f) = 0 gives the equation\n%s * v = 0." % (symbolic(f)))
        elif vs[i] == '1':
            logging.debug("Condition pert(1) = 0 gives the equation\n%s * v = 0." % (symbolic(field(1))))
        else:
            (x, y, z, xeps, yeps, zeps) = vs[i]
            if hasattr(function, "_vertices_used"):
                function._vertices_used.append(vs[i])
            eqn = equations[i]
            logging.debug("Condition pert(%s%s) + pert(%s%s) = pert(%s%s) gives the equation\n%s * v = 0." % (x, print_sign(xeps),  y, print_sign(yeps), z, print_sign(zeps), eqn))
    M = M[pivot_r]
    if return_vertices:
        vs = [ vs[i] for i in pivot_r ]
        return M, vs
    else:
        return M

def find_epsilon_interval_general(fn, perturb):
    r"""Compute the interval [minus_epsilon, plus_epsilon] such that
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

def delta_pi_general(fn, x, y, xyz_eps=(0, 0, 0)):
    r"""
    return delta_pi = fn(x, xeps) + fn(y, yeps) - fn(z, zeps).

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
    r"""
    xeps, yeps, zeps = xyz_eps
    return fn.limit(fractional(x), xeps) + fn.limit(fractional(y), yeps) - fn.limit(fractional(x + y), zeps)

def delta_pi_of_face(fn, x, y, F):
    def lim_mod_1(z, K):
        ## Annoying code because we still haven't made FastPiecewise periodic!
        fz = fractional(z)
        if fz != z:
            a, b = interval_to_endpoints(K)
            K = [a + (fz - z), b + (fz - z)]
        return fn.limit_within_relint(fz, K)
    return (fn.limit_within_relint(x, F.minimal_triple[0])
            + fn.limit_within_relint(y, F.minimal_triple[1])
            - lim_mod_1(x + y, F.minimal_triple[2]))

def delta_pi_of_face_symbolic(fn, x, y, F, simplify=False):
    def generic_point(I):
        a, b = interval_to_endpoints(I)
        return (a + b) / 2
    I, J, K = F.minimal_triple
    I = result_concrete_value(None, I)
    J = result_concrete_value(None, J)
    K = result_concrete_value(None, K)
    if simplify:
        if len(I) == 1:
            x = I[0]
        if len(J) == 1:
            y = J[0]
        if len(K) == 1:
            if len(I) > 1:
                y = K[0] - x
            elif len(J) > 1:
                x = K[0] - y
    if K[0] >= 1:
        t = -1
    else:
        t = 0
    return (fn.which_function(generic_point(I))(x)
            + fn.which_function(generic_point(J))(y)
            - fn.which_function(generic_point(K) + t)(x + y + t))

def containing_eps_1d(x, interval, with_limits=True):
    r"""
    Input:  

    - `x`: the projection of vertex `v` (of a face `F`), 

    - interval: the projection of the face `F`. 

    The projection direction is I/J/K. Note that `x` is in `interval`.

    Return: 

    - The projection of approaching limits (`\subseteq \{x^-, x, x^+\}`) that need to be considered at `v` for testing the additivity of `F`.

    If ``with_limits=False``, this computation is done for additivity of `F` sans limits.

    """
    if len(interval) == 1:
        if with_limits:
            return [0, 1, -1]
        else:
            return [0]
    elif x == interval[0]:
        return [1]
    elif x == interval[1]:
        return [-1]
    else:
        return [0]

def generate_containing_eps_triple(vertex, triple, with_limits=True):
    r"""
    Given vertex `v` of face `F`, and the 3-projection-interval triple (minimal triple) of `F`.
    Return the approaching limits {(xeps, yeps, zeps)}
    pointing inwards at `v` from containing faces of `F`,
    that should be considered for testing the additivity of `F`.

    If ``with_limits=False``, this computation is done for additivity of `F` sans limits.
    """
    xeps_list = containing_eps_1d(vertex[0], triple[0], with_limits=with_limits)
    yeps_list = containing_eps_1d(vertex[1], triple[1], with_limits=with_limits)
    zeps_list = containing_eps_1d(vertex[0] + vertex[1], triple[2], with_limits=with_limits)
    return [(xeps, yeps, zeps) for xeps in xeps_list for yeps in yeps_list for zeps in zeps_list]

# FIXME: Doctest this function
def is_additive_face(fn, face):
    r"""
    Test whether the given face is additive 
    by taking the appropriate limits (pointing inwards) at the vertices.
    """
    if face.is_2D():
        for vertex in face.vertices:
            for eps_triple in generate_containing_eps_triple(vertex, face.minimal_triple):
                try: 
                    if delta_pi_general(fn, vertex[0], vertex[1], eps_triple) != 0:
                        return False
                except ValueError: # undefined
                    return False
        return True
    elif face.is_1D():
        vertex_0 = face.vertices[0]
        vertex_1 = face.vertices[1]
        eps_triple_0 = generate_containing_eps_triple(vertex_0, face.minimal_triple)
        eps_triple_1 = generate_containing_eps_triple(vertex_1, face.minimal_triple)
        # FIXME: both eps_triple_0 and _1 have length 3? in compatible order? Yes.
        for i in range(3):
            try:
                if delta_pi_general(fn, vertex_0[0], vertex_0[1], eps_triple_0[i]) == 0 and \
                   delta_pi_general(fn, vertex_1[0], vertex_1[1], eps_triple_1[i]) == 0:
                    return True
            except ValueError: # undefined
                pass
        return False
    else:
        vertex = face.vertices[0]
        for eps_triple in itertools.chain([(0,0,0)], nonzero_eps):
            try:
                if delta_pi_general(fn, vertex[0], vertex[1], eps_triple) == 0:
                    return True
            except ValueError: # undefined
                pass
        return False

def x_y_swapped_face(face):
    vert = face.vertices
    vert_sym = [(vertex[1], vertex[0]) for vertex in vert]
    trip = face.minimal_triple
    return Face( (trip[1], trip[0], trip[2]), vertices=vert_sym )

def generate_additive_faces_general(function):
    """
    Implementation of ``generate_additive_faces`` for discontinuous piecewise linear functions.
    """
    logging.info("Computing additive faces...")
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
    #dditive_vertices_seen = {vertex for face in faces for vertex in face.vertices}
    #additive_vertices_new = additive_vertices.difference(additive_vertices_seen)
    for (x, y) in additive_vertices:
        face = Face(([x], [y], [x+y]))
        faces.append(face)
        if x != y:
            face = Face(([y], [x], [x+y]))
            faces.append(face)

    logging.info("Computing additive faces... done")
    return faces
