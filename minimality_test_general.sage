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

def limits_with_left0_right1(fn, fn_limits, x):
    if x == 0:  # take care of limit at 0-
        return fn_limits[0]
    elif x == 1: # take care of limit at 1+
        return fn_limits[-1]
    else:
        return fn.limits(x)

def generate_type_1_vertices_general(fn, comparison, continuity=True, reduced=True):
    """A generator...
    "...'general' refers to the fact that it outputs 6-tuples (x,y,z,xeps,yeps,zeps).
    When reduced=True:
        only outputs fewer triples satisfying `comparison' relation, for the purpose of minimality_check or setting up system of equations.
    When reduced=False:
        outputs all triples satisfying `comparison' relation, for the purpose of plotting nonsubadditive or additive_limit_vertices.
    """
    bkpt = fn.end_points()
    if not continuity:
        limits = fn.limits_at_end_points()
        limits[0][-1] = limits[-1][-1]
        limits[-1][1] = limits[0][1]
    for i in range(len(bkpt)):
        for j in range(i,len(bkpt)):
            x = bkpt[i]
            y = bkpt[j]
            z = fractional(x + y)
            if comparison(fn.values_at_end_points()[i] + fn.values_at_end_points()[j], fn(z)):
                yield (x, y, x+y, 0, 0, 0)
            if not continuity:
                limits_x = limits[i]
                limits_y = limits[j]
                limits_z = limits_with_left0_right1(fn, limits, z)
                if reduced and limits_x[0] == limits_x[1] == limits_x[-1] and limits_y[0] == limits_y[1] == limits_y[-1]:
                    eps_to_check = continuous_xy_eps # continuous at x and y
                else:
                    eps_to_check = nonzero_eps
                for (xeps, yeps, zeps) in eps_to_check:
                    if comparison(limits_x[xeps] + limits_y[yeps] - limits_z[zeps], 0):
                       yield (x, y, x+y, xeps, yeps, zeps)

def generate_type_2_vertices_general(fn, comparison, continuity=True, reduced=True):
    """
    When reduced=True:
        only outputs fewer triples satisfying `comparison' relation, for the purpose of minimality_check or setting up equations.
        Note: if fn is continuous at y, then fn(y-) = fn(y) = fn(y+)
    When reduced=False:
        outputs all triples satisfying `comparison' relation, for the purpose of plotting nonsubadditive or additive_limit_vertices.
    """
    bkpt = fn.end_points()
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]
    if not continuity:
        limits = fn.limits_at_end_points()
        limits[0][-1] = limits[-1][-1]
        limits[-1][1] = limits[0][1]
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
            if comparison(fn.values_at_end_points()[i] + fn(y), fn.values_at_end_points()[k]):
                yield (x, y, z, 0, 0, 0)
            if not continuity:
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

@cached_function
def generate_nonsubadditive_vertices_general(fn, continuity=True, reduced=True):
    """
    We are returning a set of 6-tuples (x, y, z, xeps, yeps, zeps),
    so that duplicates are removed, and so the result can be cached for later use.

    When reduced=True:
        only outputs fewer triples satisfying `comparison' relation, for the purpose of minimality_check.
    When reduced=False:
        outputs all triples satisfying `comparison' relation, for the purpose of plotting nonsubadditive_limit_vertices.
    """
    return set(itertools.chain( \
                generate_type_1_vertices_general(fn, operator.lt, continuity=continuity, reduced=reduced),\
                generate_type_2_vertices_general(fn, operator.lt, continuity=continuity, reduced=reduced))  )

@cached_function
def generate_additive_vertices_general(fn, continuity=True, reduced=True):
    """
    We are returning a set of 6-tuples (x, y, z, xeps, yeps, zeps),
    so that duplicates are removed, and so the result can be cached for later use.

    When reduced=True:
        only outputs fewer triples satisfying `comparison' relation, for the purpose of setting up the system of equations.
    When reduced=False:
        outputs all triples satisfying `comparison' relation, for the purpose of plotting additive_limit_vertices.
    """
    return set(itertools.chain( \
                generate_type_1_vertices_general(fn, operator.eq, continuity=continuity, reduced=reduced),\
                generate_type_2_vertices_general(fn, operator.eq, continuity=continuity, reduced=reduced)) )

def subadditivity_check_general(fn, continuity=None):
    """
    Check if fn is subadditive. Works for discontinuous functions as well.
    """
    bkpt = fn.end_points()
    if continuity == None:
        continuity = fn.is_continuous_defined()
    result = True
    for (x, y, z, xeps, yeps, zeps) in generate_nonsubadditive_vertices_general(fn, continuity=continuity, reduced=True):
        logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(xeps), y, print_sign(yeps), z, print_sign(zeps)))
        result = False
    if result:
        logging.info("pi is subadditive.")
    else:
        logging.info("Thus pi is not subadditive.")
    return result

def symmetric_check_general(fn, f, continuity=None):
    """
    Check if fn is symmetric. Works for discontinuous functions as well.
    """
    result = True
    if fn(f) != 1:
        logging.info('pi(f) is not equal to 1.')
        result = False
    bkpt = fn.end_points()
    if continuity == None:
        continuity = fn.is_continuous_defined()
    if not continuity:
        limits = fn.limits_at_end_points()
        limits[0][-1] = limits[-1][-1]
        limits[-1][1] = limits[0][1]
    for i in range(len(bkpt)):
        x = bkpt[i]
        if x == f:
            continue
        if x < f:
            y = f - x
        else:
            y = 1 + f - x
        if fn.values_at_end_points()[i] + fn(y) != 1:
            logging.info('pi(%s) + pi(%s) is not equal to 1' % (x, y))
            result = False
        if not continuity:
            limits_x = limits[i]
            limits_y = fn.limits(y)
            # no trouble at 0- and 1+ since 0 < y < 1.
            if limits_x[-1] + limits_y[1] != 1:
                logging.info('pi(%s-) + pi(%s+) is not equal to 1' % (x, y))
                result = False
            if limits_x[1] + limits_y[-1] != 1:
                logging.info('pi(%s+) + pi(%s-) is not equal to 1' % (x, y))
                result = False
    if result:
        logging.info("pi is symmetric.")
    else:
        logging.info("Thus pi is not symmetric.")
    return result

def minimality_test_general(fn, show_plots=False, f=None):
    """
    Check if fn is minimal with respect to f. Works for discontinuous functions as well.
    """
    for x in fn.values_at_end_points():
        if (x < 0) or (x > 1):
            logging.info('pi is not minimal because it does not stay in the range of [0, 1].')
            return False
    if f==None:
        f = find_f(fn, no_error_if_not_minimal_anyway=True)
        if f==None:
            return False
    if fn(0) != 0:
        logging.info('pi is NOT minimal because pi(0) is not equal to 0.')
        return False
    logging.info('pi(0) = 0')
    bkpt = fn.end_points()
    continuity = fn.is_continuous_defined()
    if not continuity:
        limits = fn.limits_at_end_points()
        limits[0][-1] = limits[-1][-1]
        limits[-1][1] = limits[0][1]
        for x in limits:
            if not ((0 <= x[-1] <=1) and (0 <= x[1] <=1)):
                logging.info('pi is not minimal because it does not stay in the range of [0, 1].')
                return False
    if subadditivity_check_general(fn, continuity) and symmetric_check_general(fn, f, continuity):
        logging.info('Thus pi is minimal.')
        is_minimal = True
    else:
        logging.info('Thus pi is NOT minimal.')
        is_minimal = False
    if show_plots:
        logging.info("Plotting 2d diagram...")
        show_plot( plot_2d_diagram_general(fn, show_function=False, continuity=continuity, known_minimal=is_minimal),\
                     show_plots, tag='2d_diagram' )
        logging.info("Plotting 2d diagram... done")
    return is_minimal

def plot_2d_diagram_general(fn, show_function=False, continuity=None, known_minimal=False):
    """
    To show only a part of it, use
    `show(plot_2d_diagram(h), xmin=0.25, xmax=0.35, ymin=0.25, ymax=0.35)`

    EXAMPLES::
    sage: h = FastPiecewise([[closed_interval(0,1/4), FastLinearFunction(4, 0)], \
    ...          [open_interval(1/4, 1), FastLinearFunction(4/3, -1/3)], \
    ...          [singleton_interval(1), FastLinearFunction(0,0)]], merge=False)
    sage: plot_2d_diagram_general(h)

    sage: h = FastPiecewise([[closed_interval(0,1/4), FastLinearFunction(4, 0)], \
    ...           [open_interval(1/4,1/2), FastLinearFunction(3, -3/4)], \
    ...           [closed_interval(1/2, 3/4), FastLinearFunction(-2, 7/4)], \
    ...           [open_interval(3/4,1), FastLinearFunction(3, -2)], \
    ...           [singleton_interval(1), FastLinearFunction(0,0)]], merge=False)
    """
    if continuity == None:
        continuity = fn.is_continuous_defined()
    p = plot_2d_complex(fn)

    ## FIXME: Need to develope code for Face of discontinuous functions
    #y = var('y')
    #faces = generate_maximal_additive_faces(function)
    #kwds = { 'legend_label': "Additive face" }
    #for face in faces:
    #    p += face.plot(**kwds)
    #    if 'legend_label' in kwds:
    #        del kwds['legend_label']

    ### For non-subadditive functions, show the points where delta_pi
    ### is negative.
    if not known_minimal:
        nonsubadditive_vertices = generate_nonsubadditive_vertices_general(fn, continuity=continuity, reduced=False)
        if continuity:
            nonsubadditive_vertices = {(x,y) for (x, y, z, xeps, yeps, zeps) in nonsubadditive_vertices}
            p += point(list(nonsubadditive_vertices), \
                                color = "red", size = 50, legend_label="Subadditivity violated", zorder=-1)
            p += point([ (y,x) for (x,y) in nonsubadditive_vertices ], color = "red", size = 50, zorder=-1)
        elif nonsubadditive_vertices:
            for (x, y, z, xeps, yeps, zeps) in nonsubadditive_vertices:
                p += plot_limit_cone_of_vertex(x, y, epstriple_to_cone((xeps, yeps, zeps)))
                if x != y:
                    p += plot_limit_cone_of_vertex(y, x, epstriple_to_cone((yeps, xeps, zeps)))
            # add legend_label
            p += point([(0,0)], color = "red", size = 50, legend_label="Subadditivity violated", zorder=-10)
            p += point([(0,0)], color = "white", size = 50, zorder=-9)
    if show_function:
        x = var('x')
        #FIXME parametric_plot doesn't work for discontinuous functions.
        p += parametric_plot((lambda x: x, lambda x: 0.3 * float(fn(x)) + 1), \
                                                (x, 0, 1), color='blue', legend_label="Function pi")
        p += parametric_plot((lambda x: - 0.3 * float(fn(x)), lambda x: x), \
                                                (x, 0, 1), color='blue')
    return p
def epstriple_to_cone(epstriple):
    try:
        return dic_eps_to_cone[epstriple]
    except KeyError:
        raise ValueError,"The limit epstriple %s does not exist." % epstriple

def plot_limit_cone_of_vertex(x, y, cone, color='red', r=0.03):
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

###### Temporary code. Look for covered intervals on a 2d-diagram ########
def plot_2d_additive_limit_vertices(fn, continuity=None):
    if continuity == None:
        continuity = fn.is_continuous_defined()
    p = plot_2d_complex(fn)
    additive_vertices = generate_additive_vertices_general(fn, continuity=continuity, reduced=False)
    if continuity:
        additive_vertices = {(x,y) for (x, y, z, xeps, yeps, zeps) in badditive_vertices}
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

def finite_dimensional_extremality_test_general(function, show_plots=False, f=None, components=None):
    """
    EXAMPLES::
    sage: logging.disable(logging.INFO)
    sage: h1 = drlm_not_extreme_2()
    sage: finite_dimensional_extremality_test_general(h1,show_plots=True, components = [ [[0, 1/4], [1/4, 1/2], [1/2, 3/4], [3/4, 1]] ])
    False
    sage: h2 = drlm_3_slope_limit()
    sage: finite_dimensional_extremality_test_general(h,show_plots=True, components = [ [[0, 1/5]], [[1/5, 1]] ])
    True
    """
    #FIXME: parameter `components' is temporary. should be removed once generate_covered_intervals_general has been developped.
    #if function is discontinuous, need to provide components, a list recording covered_intervals.
    # dg_2_step_mir_limit() example: components = [ [[0, 1/5], [1/5, 2/5], [2/5, 3/5], [3/5, 1]] ]
    # drlm_3_slope_limit() example: components = [ [[0, 1/5]], [[1/5, 1]] ]
    # rlm_dpl1_fig3_lowerleft() example: components = [ [[0, 1/4]], [[1/4, 5/8], [5/8, 1]] ]
    # drlm_not_extreme_2() example: components = [ [[0, 1/4], [1/4, 1/2], [1/2, 3/4], [3/4, 1]] ]

    continuity = function.is_continuous_defined()
    # temporary
    if continuity:
        components = generate_covered_intervals(function)
    elif components == None:
        raise ValueError,"Need to provide ``components'' for discontinuous funcitons."

    field = function(0).parent().fraction_field()
    symbolic = generate_symbolic_general(function, components, field=field)
    equation_matrix = generate_additivity_equations_general(function, symbolic, field, f=f, continuity=continuity)
    slope_jump_vects = equation_matrix.right_kernel().basis()
    logging.info("Solution space has dimension %s" % len(slope_jump_vects))
    if len(slope_jump_vects) == 0:
        logging.info("Thus the function is extreme.")
        return True
    else:
        for basis_index in range(len(slope_jump_vects)):
            slope_jump = slope_jump_vects[basis_index]
            perturbation = function._perturbation = slope_jump * symbolic
            check_perturbation_general(function, perturbation, continuity=continuity,
                                        show_plots=show_plots, show_plot_tag='perturbation-%s' % (basis_index + 1),
                                        legend_title="Basic perturbation %s" % (basis_index + 1))
        return False

def generate_symbolic_general(function, components, field=None):
    n = len(components)
    intervals_and_slopes = []
    for component, slope in itertools.izip(components, range(n)):
        intervals_and_slopes.extend([ (interval, slope) for interval in component ])
    intervals_and_slopes.sort()
    bkpt = [ field(interval[0]) for interval, slope in intervals_and_slopes ] + [field(1)]
    limits = [function.limits(x) for x in bkpt]
    num_jumps = sum([(x[-1] != x[0]) + (x[0] != x[1]) for x in limits]) - 2 # None at 0- and 1+, so sum-2
    vector_space = VectorSpace(field, n + num_jumps)
    unit_vectors = vector_space.basis()
    slopes = [ unit_vectors[slope] for interval, slope in intervals_and_slopes ]
    m = len(slopes)
    # Set up symbolic function
    current_value = zeros = vector_space.zero()
    pieces = []
    j = n
    for i in range(m):
        pieces.append([singleton_interval(bkpt[i]), FastLinearFunction(zeros, current_value)])
        if limits[i][0] != limits[i][1]: # jump at bkpt[i]+
            current_value += unit_vectors[j]
            j += 1
        pieces.append([open_interval(bkpt[i], bkpt[i+1]), FastLinearFunction(slopes[i], current_value - slopes[i]*bkpt[i])])
        current_value += slopes[i] * (bkpt[i+1] - bkpt[i])
        if limits[i+1][-1] != limits[i+1][0]: # jump at bkpt[i+1]-
            current_value += unit_vectors[j]
            j += 1
    pieces.append([singleton_interval(bkpt[m]), FastLinearFunction(zeros, current_value)])
    return FastPiecewise(pieces, merge=True)

def generate_additivity_equations_general(function, symbolic, field, f=None, continuity=True):
    equations = []
    if f == None:
        f = find_f(function)
    equations.append(symbolic(f))
    equations.append(symbolic(field(1)))
    for (x, y, z, xeps, yeps, zeps) in generate_additive_vertices_general(function, continuity=continuity, reduced=True):
        # handle 0- and 1+
        x, xeps = periodic_one_limit(x, xeps, field)
        y, yeps = periodic_one_limit(y, yeps, field)
        z, zeps = periodic_one_limit(z, zeps, field)
        new_equation = symbolic.limit(x, xeps) + symbolic.limit(y, yeps) - symbolic.limit(z, zeps)
        equations.append(new_equation)
    return  matrix(field, equations)

def find_epsilon_interval_general(fn, perturb, continuity=True):
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
    fn_limits[0][-1] = fn_limits[-1][-1]
    fn_limits[-1][1] = fn_limits[0][1]
    perturb_limits = [perturb.limits(x) for x in bkpt]
    perturb_limits[0][-1] = perturb_limits[-1][-1]
    perturb_limits[-1][1] = perturb_limits[0][1]

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
            perturb_z = limits_with_left0_right1(perturb, perturb_limits, z)
            fn_z = limits_with_left0_right1(fn, fn_limits, z)
            if continuity:
                eps_to_check = {(0, 0, 0)}
            elif fn_x[0] == fn_x[1] == fn_x[-1] and fn_y[0] == fn_y[1] == fn_y[-1]:
                # if fn is continuous at x and y, then so is perturb.
                eps_to_check = continuous_xy_eps
            else:
                eps_to_check = nonzero_eps
            for (xeps, yeps, zeps) in eps_to_check:
                delta_perturb = perturb_x[xeps] + perturb_y[yeps] - perturb_z[zeps]
                if delta_perturb != 0:
                    delta_fn = fn_x[xeps] + fn_y[yeps] - fn_z[zeps]
                    if delta_fn == 0:
                        logging.info("Zero epsilon encountered for x = %s, y = %s" % (bkpt[i], bkpt[j]) )
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
            perturb_y = limits_with_left0_right1(perturb, perturb_limits, y)
            fn_y = limits_with_left0_right1(fn, fn_limits, y)

            if continuity:
                eps_to_check = {(0, 0, 0)}
            elif not (fn_y[0] == fn_y[1] == fn_y[-1]):
                # then y is a in bkpt. this is done in type1check.
                eps_to_check = {}
            else:
                eps_to_check = type2_reduced_eps

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

def check_perturbation_general(fn, perturb, continuity=True, \
                                show_plots=False, show_plot_tag='perturbation', xmin=0, xmax=1, **show_kwds):
    epsilon_interval = fn._epsilon_interval = find_epsilon_interval_general(fn, perturb, continuity=continuity)
    epsilon = min(abs(epsilon_interval[0]), epsilon_interval[1])
    logging.info("Epsilon for constructed perturbation: %s" % epsilon)
    if show_plots:
        logging.info("Plotting perturbation...")
        p = plot(rescale_to_amplitude_general(perturb, 1/10), xmin=xmin, xmax=xmax, color='magenta', legend_label="perturbation (rescaled)")

        p += plot(fn + epsilon_interval[0] * perturb, xmin=xmin, xmax=xmax, color='red', legend_label="-perturbed (min)")
        p += plot(fn + epsilon_interval[1] * perturb, xmin=xmin, xmax=xmax, color='blue', legend_label="+perturbed (max)")
        if -epsilon != epsilon_interval[0]:
            p += plot(fn + (-epsilon) * perturb, xmin=xmin, xmax=xmax, color='orange', legend_label="-perturbed (matches max)")
        elif epsilon != epsilon_interval[1]:
            p += plot(fn + epsilon * perturb, xmin=xmin, xmax=xmax, color='cyan', legend_label="+perturbed (matches min)")
        p += plot(fn, xmin=xmin, xmax=xmax, color='black', thickness=2, legend_label="original function")
        show_plot(p, show_plots, tag=show_plot_tag, **show_kwds)
        logging.info("Plotting perturbation... done")
    assert epsilon > 0, "Epsilon should be positive, something is wrong"
    logging.info("Thus the function is not extreme.")

def rescale_to_amplitude_general(perturb, amplitude):
    """For plotting purposes, rescale the function `perturb` so that its
    maximum absolute function value is `amplitude`.
    """
    current_amplitude = max([ abs(x) for limits in perturb.limits_at_end_points() for x in limits ])
    if current_amplitude != 0:
        return perturb * (amplitude/current_amplitude)
    else:
        return perturb

def periodic_one_limit(x, xeps, field):
    x = fractional(x)
    if x == field(0) and xeps == -1:
        x = field(1)
    elif x == field(1) and xeps == 1:
        x = field(0)
    return x, xeps
