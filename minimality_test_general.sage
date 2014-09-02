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
                if z == 0:
                    # take care of limit at 0-
                    limits_z = limits[0]
                elif z == 1:
                    # take care of limit at 1+
                    limits_z = limits[-1]
                else:
                    limits_z = fn.limits(z)
                if reduced and limits_x[0] == limits_x[1] == limits_x[-1] and limits_y[0] == limits_y[1] == limits_y[-1]:
                    # continuous at x and y
                    for zeps in [1, -1]:
                        # note that (0 0 0) has alreadly been checked.
                        if comparison(limits_x[0] + limits_y[0] - limits_z[eps], 0):
                            yield (x, y, x+y, zeps, zeps, zeps)
                else:
                    for eps in [1, -1]:
                        # (+ + +), (- - -). note that (0 0 0) has alreadly been checked.
                        if comparison(limits_x[eps] + limits_y[eps] - limits_z[eps], 0):
                            yield (x, y, x+y, eps, eps, eps)
                    for eps in [1, -1]:
                        for zeps in [0, 1, -1]:
                            # (+ - 0), (+ - +), (+ - -), (- + 0), (- + +), (- + -)
                            if comparison(limits_x[eps] + limits_y[-eps] - limits_z[zeps], 0):
                                yield (x, y, x+y, eps, -eps, zeps)
                        if comparison(limits_x[0] + limits_y[eps] - limits_z[eps], 0):
                            # (0 + +), (0 - -)
                            yield (x, y, x+y, 0, eps, eps)
                        if comparison(limits_x[eps] + limits_y[0] - limits_z[eps], 0):
                            # (+ 0 +), (- 0 -)
                            yield (x, y, x+y, eps, 0, eps)

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
                    for zeps in [1, -1]:
                        for xeps in [0, 1, -1]:
                            if comparison(limits_x[xeps] + limits_y[0], limits_z[zeps]):
                                yield (x, y, z, xeps, zeps, zeps)
                    for xeps in [1, -1]:
                        if comparison(limits_x[xeps] + limits_y[0], limits_z[0]):
                            yield (x, y, z, xeps, -xeps, 0)
                else:
                    for eps in [1, -1]:
                        # (+ + +), (- - -). note that (0 0 0) has alreadly been checked.
                        if comparison(limits_x[eps] + limits_y[eps] - limits_z[eps], 0):
                            yield (x, y, x+y, eps, eps, eps)
                    for eps in [1, -1]:
                        for zeps in [0, 1, -1]:
                            # (+ - 0), (+ - +), (+ - -), (- + 0), (- + +), (- + -)
                            if comparison(limits_x[eps] + limits_y[-eps] - limits_z[zeps], 0):
                                yield (x, y, x+y, eps, -eps, zeps)
                        if comparison(limits_x[0] + limits_y[eps] - limits_z[eps], 0):
                            # (0 + +), (0 - -)
                            yield (x, y, x+y, 0, eps, eps)
                        if comparison(limits_x[eps] + limits_y[0] - limits_z[eps], 0):
                            # (+ 0 +), (- 0 -)
                            yield (x, y, x+y, eps, 0, eps)
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
    dic = { (-1,-1,-1): [(-1, 0), (0, -1)], \
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
    try:
        return dic[epstriple]
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
    #FIXME: parameter `components' is temporary. should be removed once generate_covered_intervals_general has been developped.
    continuity = function.is_continuous_defined()
    if continuity:    # Copy from functions.sage
        symbolic, components, field = symbolic_piecewise(function)
        equations = generate_additivity_equations(function, symbolic, field, f=f)
        slopes_vects = equations.right_kernel().basis()
        logging.info("Solution space has dimension %s" % len(slopes_vects))
        if len(slopes_vects) == 0:
            logging.info("Thus the function is extreme.")
            return True
        else:
            for basis_index in range(len(slopes_vects)):
                slopes = list(slopes_vects[basis_index])
                perturbation = function._perturbation = generate_compatible_piecewise_function(components, slopes)
                check_perturbation(function, perturbation,
                                   show_plots=show_plots, show_plot_tag='perturbation-%s' % (basis_index + 1),
                                   legend_title="Basic perturbation %s" % (basis_index + 1))
            return False
    else:
        #FIXME: if function is discontinuous, need to provide components, a list recording covered_intervals.
        # dg_2_step_mir_limit() example: components = [ [[0, 1/5], [1/5, 2/5], [2/5, 3/5], [3/5, 1]] ]
        # drlm_3_slope_limit() example: components = [ [[0, 1/5]], [[1/5, 1]] ]
        # rlm_dpl1_fig3_lowerleft_not_extreme() example: components = [ [[0, 1/4]], [[1/4, 5/8], [5/8, 1]] ]
        if components == None:
            raise ValueError,"Need to provide ``components'' for discontinuous funcitons."
        field = function(0).parent().fraction_field()
        bkpt_set = {interval[0] for component in components for interval in component}
        # bkpt_set.add(1) are the possible discontinuous points
        n = len(components)
        m = len(bkpt_set)
        vector_space = VectorSpace(field, n + 2*m)
        component_slopes = vector_space.basis()[0:n]
        intervals_and_slopes = []
        for component, slope in itertools.izip(components, component_slopes):
            intervals_and_slopes.extend([ (interval, slope) for interval in component ])
        intervals_and_slopes.sort()
        bkpt = [ field(interval[0]) for interval, slope in intervals_and_slopes ] + [field(1)]
        # Note: bkpt[0] == 0; bkpt[m] == 1
        slopes = [ slope for interval, slope in intervals_and_slopes ]
        jumps = vector_space.basis()[n: n + 2*m]
        # Initialize symbolic function
        current_value = zeros = vector_space.zero()
        pieces = []
        # Set up symbolic function
        for i in range(m):
            pieces.append([singleton_interval(bkpt[i]), FastLinearFunction(zeros, current_value)])
            current_value += jumps[2*i]
            pieces.append([open_interval(bkpt[i], bkpt[i+1]), FastLinearFunction(slopes[i], current_value - slopes[i]*bkpt[i])])
            current_value += slopes[i] * (bkpt[i+1] - bkpt[i])+ jumps[2*i + 1]
        pieces.append([singleton_interval(bkpt[m]), FastLinearFunction(zeros, current_value)])
        symbolic = FastPiecewise(pieces, merge=True)
        global symfn
        symfn = symbolic
        # Set up system of equations
        equations = []
        limit_values = function.limits(bkpt[0])
        for i in range(m):
            if limit_values[0] == limit_values[1]:
                # if \pi is right continuous at bkpt[i], so is any \phi
                equations.append(jumps[2*i])
            limit_values = function.limits(bkpt[i+1])
            if limit_values[-1] == limit_values[0]:
                # if \pi is left continuous at bkpt[i+1], so is any \phi
                equations.append(jumps[2*i + 1])
        if f == None:
            f = find_f(function)
        equations.append(symbolic(f))
        equations.append(symbolic(bkpt[m]))
        for (x, y, z, xeps, yeps, zeps) in generate_additive_vertices_general(function, continuity=continuity, reduced=True):
            # handle 0- and 1+
            x, xeps = periodic_one_limit(x, xeps, field)
            y, yeps = periodic_one_limit(y, yeps, field)
            z, zeps = periodic_one_limit(z, zeps, field)
            new_equation = symbolic.limit(x, xeps) + symbolic.limit(y, yeps) - symbolic.limit(z, zeps)
            equations.append(new_equation)
        equation_matrix = matrix(field, equations)
        # Solve system of equations.
        # A solution gives [slope_value[0],..., slope_value[n-1], jump_value[0], ..., jump_value[2*m -1]]
        slope_jump_vects = equation_matrix.right_kernel().basis()
        logging.info("Solution space has dimension %s" % len(slope_jump_vects))
        if len(slope_jump_vects) == 0:
            logging.info("Thus the function is extreme.")
            return True
        else:
            #FIXME: perturbations...
            #for basis_index in range(len(slope_jump_vects)):
            #    slope_jump = list(slope_jump_vects[basis_index])
            #    #perturbation = function._perturbation = generate_compatible_piecewise_function(components, slopes)
            #    check_perturbation(function, perturbation,
            #                        show_plots=show_plots, show_plot_tag='perturbation-%s' % (basis_index + 1),
            #                        legend_title="Basic perturbation %s" % (basis_index + 1))
            return False

#def generate_compatible_piecewise_function_general(components, slope_vars, field=None):

def periodic_one_limit(x, xeps, field):
    x = fractional(x)
    if x == field(0) and xeps == -1:
        x = field(1)
    elif x == field(1) and xeps == 1:
        x = field(0)
    return x, xeps
