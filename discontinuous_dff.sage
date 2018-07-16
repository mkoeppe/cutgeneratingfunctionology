# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

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

def generate_maximal_additive_faces_general_dff(function):
    logging.info("Computing maximal additive faces...")
    bkpt = function.end_points()
    n = len(bkpt) - 1
    I_list = J_list = K_list = [ (bkpt[i], bkpt[i+1]) for i in range(n) ]
    faces = []
    # 2D faces
    for i in range(n):
        for j in range(i, n):
            IplusJ = interval_sum(I_list[i],J_list[j])
            for k in generate_overlapping_interval_indices(IplusJ, bkpt):
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
            for k in generate_overlapping_interval_indices(IplusJ, bkpt):
                if len(interval_intersection(IplusJ, K_list[k])) == 2:
                    face = Face( ([bkpt[i]], J_list[j], K_list[k]) )
                    if is_additive_face(function, face): 
                        faces.append(face)
                        faces.append(x_y_swapped_face(face))
    # 1D diagonal faces
    for i in range(n):
        for j in range(i, n):
            interval_K = interval_sum(I_list[i],J_list[j])
            for k in generate_overlapping_interval_indices(interval_K, bkpt):
                if interval_K[0] < bkpt[k] < interval_K[1]:
                    face = Face( (I_list[i], J_list[j], [bkpt[k]]) )
                    if is_additive_face(function, face): 
                        faces.append(face)
                        if i != j:
                            faces.append(x_y_swapped_face(face))

    # 0D faces
    additive_vertices = {(x,y) for (x, y, z, xeps, yeps, zeps) in generate_additive_vertices_general_dff(function) if x != 1 and y != 1}
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

def generate_type_1_vertices_general_dff(fn, comparison):
    bkpt=fn.end_points()
    limits=[fn.limits(x) for x in bkpt]
    for i in range(len(bkpt)):
        for j in range(i,len(bkpt)):
            x=bkpt[i]
            y=bkpt[j]
            z=x+y
            if z<=1:
                if comparison(limits[i][0] + limits[j][0], fn(z)):
                    yield (x, y, x+y, 0, 0, 0)
                limits_x = limits[i]
                limits_y = limits[j]
                limits_z = fn.limits(z)
                for (xeps, yeps, zeps) in nonzero_eps:
                    if x==0 and xeps==-1:
                        continue
                    if y==0 and yeps==-1:
                        continue
                    if z==0 and zeps==-1:
                        continue
                    if x==1 and xeps==1:
                        continue
                    if y==1 and yeps==1:
                        continue
                    if z==1 and zeps==1:
                        continue
                    if comparison(limits_x[xeps] + limits_y[yeps] - limits_z[zeps], 0):
                        yield (x, y, x+y, xeps, yeps, zeps)

def generate_type_2_vertices_general_dff(fn, comparison):
    bkpt=fn.end_points()
    limits=[fn.limits(x) for x in bkpt]
    for i in range(len(bkpt)):
        for j in range(i,len(bkpt)):
            x=bkpt[i]
            z=bkpt[j]
            y=z-x
            if comparison(limits[i][0] + fn(y), limits[j][0]):
                yield (x, y, x+y, 0, 0, 0)
            limits_x = limits[i]
            limits_y = fn.limits(y)
            limits_z = limits[j]
            for (xeps, yeps, zeps) in nonzero_eps:
                if x==0 and xeps==-1:
                    continue
                if y==0 and yeps==-1:
                    continue
                if z==0 and zeps==-1:
                    continue
                if x==1 and xeps==1:
                    continue
                if y==1 and yeps==1:
                    continue
                if z==1 and zeps==1:
                    continue
                if comparison(limits_x[xeps] + limits_y[yeps] - limits_z[zeps], 0):
                    yield (x, y, x+y, xeps, yeps, zeps)

def generate_nonsymmetric_vertices_general_dff(fn):
    bkpt = fn.end_points()
    limits = fn.limits_at_end_points()
    for i in range(len(bkpt)):
        x=bkpt[i]
        y=1-x
        if limits[i][0] + fn(y) != 1:
            yield (x, y, 0, 0)
        if not fn.is_continuous():
            limits_x = limits[i]
            limits_y = fn.limits(y)
            if x==0:
                limits_x[-1]=0
                limits_y[1]=1
            if x==1:
                limits_x[1]=1
                limits_y[-1]=0
            if limits_x[-1] + limits_y[1] != 1:
                yield (x, y, -1, 1)
            if limits_x[1] + limits_y[-1] != 1:
                yield (x, y, 1, -1)

def generate_nonsuperadditive_vertices_general(fn):
    return set(itertools.chain( \
                generate_type_1_vertices_general_dff(fn, operator.gt),\
                generate_type_2_vertices_general_dff(fn, operator.gt)))

def generate_additive_vertices_general_dff(fn):
    return set(itertools.chain( \
                generate_type_1_vertices_general_dff(fn, operator.eq),\
                generate_type_2_vertices_general_dff(fn, operator.eq)))

def symmetry_test_general_dff(fn):
    result = True
    for (x, y, xeps, yeps) in generate_nonsymmetric_vertices_general_dff(fn):
        logging.info('f(%s%s) + f(%s%s) is not equal to 1' % (x, print_sign(xeps), y, print_sign(yeps)))
        result = False
    if result:
        logging.info('f is symmetric.')
    else:
        logging.info('Thus f is not symmetric.')
    return result

def superadditivity_test_general(fn):
    """
    Check if `fn` is superadditive.
    """
    result = True
    for (x, y, z, xeps, yeps, zeps) in generate_nonsuperadditive_vertices_general(fn):
        logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) > 0" % (x, print_sign(xeps), y, print_sign(yeps), z, print_sign(zeps)))
        result = False
    if result:
        logging.info("f is superadditive.")
    else:
        logging.info("Thus f is not superadditive.")
    return result

def maximality_test_general_dff(fn):
    bkpt = fn.end_points()
    limits=[fn.limits(x) for x in bkpt]
    for i in range(len(limits)):
        for j in range(3):
            if (limits[i][j] < 0) or (limits[i][j] > 1):
                logging.info('f is not maximal because it does not stay in the range of [0, 1].')
                return False
    if fn(0) != 0:
        logging.info('f is NOT maximal because f(0) is not equal to 0.')
        return False
    logging.info('f(0) = 0')
    if superadditivity_test_general(fn) and symmetry_test_general_dff(fn):
        logging.info('Thus f is maximal.')
        is_maximal = True
    else:
        logging.info('Thus f is NOT maximal.')
        is_maximal = False
    return is_maximal

def generate_symbolic_general_dff(function, components):
    n=len(components)
    intervals_and_slopes = []
    for component, slope in itertools.izip(components, range(n)):
        intervals_and_slopes.extend([ (interval, slope) for interval in component ])
    intervals_and_slopes.sort(key=lambda (i, s):coho_interval_left_endpoint_with_epsilon(i))
    field=function(0).parent().fraction_field()
    bkpt = [ field(interval[0]) for interval, slope in intervals_and_slopes ] + [field(1)]
    limits = [function.limits(x) for x in bkpt]
    num_jumps = sum([x[0] != x[1] for x in limits[0:-1]])
    vector_space = VectorSpace(field, n + num_jumps)
    unit_vectors = vector_space.basis()
    slopes = [ unit_vectors[slope] for interval, slope in intervals_and_slopes ]
    jump_vectors=unit_vectors[n:n+num_jumps]
    current_value = zeros = vector_space.zero()
    pieces = []
    j = 0
    k=0
    m = len(slopes)
    for i in range(m):
        pieces.append([singleton_interval(bkpt[i]), FastLinearFunction(zeros, current_value)])
        if limits[i][0] != limits[i][1]:
            current_value += jump_vectors[j] 
            j=j+1   
        pieces.append([open_interval(bkpt[i], bkpt[i+1]), FastLinearFunction(slopes[i], current_value - slopes[i]*bkpt[i])])
        current_value += slopes[i] * (bkpt[i+1] - bkpt[i])
        if limits[i+1][-1] != limits[i+1][0] and (1-bkpt[i+1]) in bkpt and function.limits(1-bkpt[i+1])[0] != function.limits(1-bkpt[i+1])[1]:
            current_value += jump_vectors[num_jumps-k-1]
            k=k+1
    pieces.append([singleton_interval(1), FastLinearFunction(zeros, current_value)])        
    symbolic_function = FastPiecewise(pieces, merge=True)
    return symbolic_function

def delta_pi_dff_general(fn,x,y,xeps,yeps,zeps):
    return fn.limit(x, xeps) + fn.limit(y, yeps) - fn.limit(x+y, zeps)

def generate_additivity_equations_general_dff(function, symbolic, field):
    bkpt=function.end_points()
    bkpt1=symbolic.end_points()
    vs = list(generate_additive_vertices_general_dff(function))
    if function.limit(bkpt[1],-1)==0:
        equations = [delta_pi_dff_general(symbolic, x, y,xeps, yeps, zeps) \
                     for (x, y, z, xeps, yeps, zeps) in generate_additive_vertices_general_dff(function) ] \
                   + [symbolic.limit(bkpt[1],-1)] \
                   + [symbolic(1)]
    else:
        equations = [delta_pi_dff_general(symbolic, x, y,xeps, yeps, zeps) \
                     for (x, y, z, xeps, yeps, zeps) in generate_additive_vertices_general_dff(function) ] \
                   + [symbolic(1)]
    return matrix(field, equations)
 
def find_epsilon_interval_general_dff(fn, perturb):   
    

    logging.info("Finding epsilon interval for perturbation...")
    fn_bkpt = fn.end_points()
    perturb_bkpt = perturb.end_points()
    bkpt = merge_bkpt(fn_bkpt,perturb_bkpt)
    

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
            if bkpt[i] + bkpt[j]<=1:
                z = bkpt[i] + bkpt[j]
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
                        if delta_fn <0:
                            epsilon_upper_bound = -delta_fn / abs(delta_perturb)
                            if delta_perturb > 0:
                                if -epsilon_upper_bound > best_minus_epsilon_lower_bound:
                                    best_minus_epsilon_lower_bound = -epsilon_upper_bound
                            else:
                                if epsilon_upper_bound < best_plus_epsilon_upper_bound:
                                    best_plus_epsilon_upper_bound = epsilon_upper_bound
    # type2check
    # if x==0, there is no left limit, need to be fixed
    for i in range(len(bkpt)):
        perturb_x = perturb_limits[i]
        fn_x = fn_limits[i]
        for k in range(i, len(bkpt)):
            perturb_z = perturb_limits[k]
            fn_z = fn_limits[k]
            y = bkpt[k] - bkpt[i]
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
                    if delta_fn<0:
                        epsilon_upper_bound = -delta_fn / abs(delta_perturb)
                        if delta_perturb > 0:
                            if -epsilon_upper_bound > best_minus_epsilon_lower_bound:
                                best_minus_epsilon_lower_bound = -epsilon_upper_bound
                        else:
                            if epsilon_upper_bound < best_plus_epsilon_upper_bound:
                                best_plus_epsilon_upper_bound = epsilon_upper_bound
    logging.info("Finding epsilon interval for perturbation... done.  Interval is %s", [best_minus_epsilon_lower_bound, best_plus_epsilon_upper_bound])
    return best_minus_epsilon_lower_bound, best_plus_epsilon_upper_bound

def generate_perturbations_finite_dimensional_general_dff(f):
    covered_intervals=generate_covered_intervals(f)
    uncovered_intervals=generate_uncovered_intervals(f)
    if uncovered_intervals:
        components = copy(covered_intervals)
        components.extend([int] for int in uncovered_intervals)
    else:
        components=covered_intervals
    field = f(0).parent().fraction_field()
    symbolic=generate_symbolic_general_dff(f,components)
    equation_matrix = generate_additivity_equations_general_dff(f, symbolic, field)
    slope_jump_vects = equation_matrix.right_kernel().basis()
    logging.info("Finite dimensional test: Solution space has dimension %s" % len(slope_jump_vects))
    for basis_index in range(len(slope_jump_vects)):
        slope_jump = slope_jump_vects[basis_index]
        perturbation = slope_jump * symbolic
        yield perturbation

def extremality_test_general_dff(fn):
    """Still in progress
    """
    covered_intervals=generate_covered_intervals(fn)
    uncovered_intervals=generate_uncovered_intervals(fn)
    if not maximality_test_general_dff(fn):
        logging.info("Not maximal, thus NOT extreme.")
        return False
    if uncovered_intervals:
        gen=generate_perturbations_equivariant(fn)
        for perturbation in gen:
            epsilon_interval = find_epsilon_interval_general_dff(fn, perturbation)
            epsilon = min(abs(epsilon_interval[0]), abs(epsilon_interval[1]))
            if epsilon>0:
                logging.info("Not extreme")
                return False
    gen=generate_perturbations_finite_dimensional_general_dff(fn)
    for perturbation in gen:
        epsilon_interval = find_epsilon_interval_general_dff(fn, perturbation)
        epsilon = min(abs(epsilon_interval[0]), abs(epsilon_interval[1]))
        if epsilon>0:
            logging.info("Not extreme")
            return False
    logging.info("extreme")
    return True

