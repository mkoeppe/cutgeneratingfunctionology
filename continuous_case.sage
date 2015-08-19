# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

########## Code for Continuous Case ###########

def generate_nonsymmetric_vertices_continuous(fn, f):
    bkpt = fn.end_points()
    for i in range(len(bkpt)):
        x = bkpt[i]
        if x == f:
            continue
        if x < f:
            y = f - x
        else:
            y = 1 + f - x
        if fn.values_at_end_points()[i] + fn(y) != 1:
            yield (x, y, 0, 0)

def generate_type_1_vertices_continuous(fn, comparison):
    """Output 6-tuples (x, y, z,xeps, yeps, zeps).
    """
    bkpt = fn.end_points()
    return ( (x, y, x+y, 0, 0, 0) for x in bkpt for y in bkpt if x <= y and comparison(delta_pi(fn,x,y), 0) ) # generator comprehension

def generate_type_2_vertices_continuous(fn, comparison):
    bkpt = fn.end_points()
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]
    return ( (x, z-x, z, 0, 0, 0) for x in bkpt for z in bkpt2 if x < z < 1+x and comparison(delta_pi(fn, x, z-x), 0) ) # generator comprehension

def modified_delta_pi(fn, fn_values, pts, i, j):
    return fn_values[i] + fn_values[j] - fn(fractional(pts[i]+pts[j])) 

def modified_delta_pi2(fn, fn_values2, pts, i, j):
    return fn_values2[i] + fn(fractional(pts[j] - pts[i])) - fn_values2[j]  

def generate_maximal_additive_faces_continuous(function):
    logging.info("Computing maximal additive faces...")
    bkpt = function.end_points()
    bkpt2 = []
    for i in range(len(bkpt)-1):
        bkpt2.append(bkpt[i])
    for i in range(len(bkpt)):
        bkpt2.append(bkpt[i]+1)      
   
    I_list = []
    J_list = []
    K_list = []
    
    intervals = []
    intervalsK = []
    
    for i in range(len(bkpt)-1):
        intervals.append([bkpt[i],bkpt[i+1]])
        
    for i in range(len(bkpt2)-1):
        intervalsK.append([bkpt2[i],bkpt2[i+1]])
        
    I_list = intervals
    J_list = intervals
    K_list = intervalsK
    
    additive_face = {}
    additive_vertices = {}
    faces = []
    for i in range(len(I_list)):
        for j in range(i, len(J_list)):
            for k in range(len(K_list)):
                # Check if int(I+J) intersects int(K) is non-empty.
                if len(interval_intersection(interval_sum(I_list[i],J_list[j]),K_list[k])) == 2:
                    temp_verts = verts(I_list[i],J_list[j],K_list[k])
                    temp = []
                    keep = False
                    if temp_verts != []:
                        for vertex in temp_verts:
                            if delta_pi(function, vertex[0],vertex[1]) == 0:
                                temp.append(vertex)
                                keep = True
                    if len(temp) == 2:
                        if temp[0][0] == temp[1][0]:
                            if temp[1][0] == I_list[i][0]:
                                if (i-1,j,k) in additive_face:
                                   keep = False
                            else:
                                keep = False
                        elif temp[0][1] == temp[1][1]: 
                            if temp[1][1] == J_list[j][0]:
                                if (i,j-1,k) in additive_face:
                                    keep = False
                            else:
                                keep = False
                        elif temp[0][0] + temp[0][1] == temp[1][0] + temp[1][1]: 
                            if temp[1][0] + temp[1][1] == K_list[k][0]:
                                if (i,j,k-1) in additive_face:
                                    keep = False
                            else:
                                keep = False
                    elif len(temp) == 1:
                        if temp[0][0] == I_list[i][0] and temp[0][1] == J_list[j][0] \
                            and temp[0][0] + temp[0][1] != K_list[k][1]:
                            keep = True
                        elif temp[0][0] == I_list[i][0] and temp[0][0] + temp[0][1] == K_list[k][0] \
                            and temp[0][1] != J_list[j][1]:
                            keep = True
                        elif temp[0][1] == J_list[j][0] and temp[0][0] + temp[0][1] == K_list[k][0] \
                            and temp[0][0] != I_list[i][1]:     
                            keep = True
                        else:
                            keep = False
                        if keep:
                            if (temp[0][0],temp[0][1]) in additive_vertices:
                                keep = False  
                            # if keep:
                            #     print I_list[i], J_list[j], K_list[k]      
                    if keep:
                        
                        additive_face[(i,j,k)] = temp
                        
                        for vert in temp:
                            additive_vertices[vert] = True
                        
                        trip = projections(temp)
                        faces.append(Face(trip, vertices=temp, is_known_to_be_minimal=True))

                        if i != j:
                            temp_swap = []
                            for vert in temp:
                                vert_new = [vert[1],vert[0]]
                                temp_swap.append(vert_new)
                            trip_swap = [trip[1], trip[0], trip[2]] #same as: trip_swap = projections(temp_swap)
                            faces.append(Face(trip_swap, vertices=temp_swap, is_known_to_be_minimal=True))
                        
    logging.info("Computing maximal additive faces... done")
    return faces

def find_epsilon_interval_continuous(fn, perturb):
    """Compute the interval [minus_epsilon, plus_epsilon] such that 
    (fn + epsilon * perturb) is subadditive for epsilon in this interval.
    Assumes that fn is subadditive.

    If one of the epsilons is 0, the function bails out early and returns 0, 0.
    """
    logging.info("Finding epsilon interval for perturbation...")
    fn_bkpt = fn.end_points()
    perturb_bkpt = perturb.end_points()
    bkpt_refinement = merge_bkpt(fn_bkpt,perturb_bkpt)
    bkpt_refinement2 = []
    length1 = len(bkpt_refinement)
    for i in range(length1 - 1):
        bkpt_refinement2.append(bkpt_refinement[i])
    for i in range(length1):
        bkpt_refinement2.append(bkpt_refinement[i]+1)
    length2 = length1 + length1 - 1  

    fn_values = []
    perturb_values = []
    for pt in bkpt_refinement2:
        fn_values.append(fn(fractional(pt)))
        perturb_values.append(perturb(fractional(pt)))
    
    best_minus_epsilon_lower_bound = -10000
    best_plus_epsilon_upper_bound = +10000
    # FIXME: We want to say infinity instead; but bool(SR(2) < infinity) ==> False
    for i in range(length1):
        for j in range(i,length1):
            a = modified_delta_pi(perturb, perturb_values, bkpt_refinement, i, j)
            if a != 0:
                b = modified_delta_pi(fn, fn_values, bkpt_refinement, i, j) 
                if b == 0:
                    logging.info("Zero epsilon encountered for x = %s, y = %s" % (bkpt_refinement[i], bkpt_refinement[j]))
                    return 0, 0 # See docstring
                epsilon_upper_bound = b/(abs(a))
                if a > 0:
                    if -epsilon_upper_bound > best_minus_epsilon_lower_bound:
                        best_minus_epsilon_lower_bound = -epsilon_upper_bound
                else:
                    if epsilon_upper_bound < best_plus_epsilon_upper_bound:
                        best_plus_epsilon_upper_bound = epsilon_upper_bound

    for i in range(length1):
        for j in range(length2):
            if bkpt_refinement2[j] - bkpt_refinement[i] > 0:
                a = modified_delta_pi2(perturb, perturb_values, bkpt_refinement2, i, j)
                if a != 0:
                    b = modified_delta_pi2(fn, fn_values, bkpt_refinement2, i, j) 

                    if b == 0:
                        logging.info("Zero epsilon encountered for x = %s, y = %s" % (bkpt_refinement2[i], bkpt_refinement2[j] - bkpt_refinement2[i]))
                        return 0, 0 # See docstring
                    epsilon_upper_bound = b/(abs(a)) 
                    if a > 0:
                        if -epsilon_upper_bound > best_minus_epsilon_lower_bound:
                            best_minus_epsilon_lower_bound = -epsilon_upper_bound
                    else:
                        if epsilon_upper_bound < best_plus_epsilon_upper_bound:
                            best_plus_epsilon_upper_bound = epsilon_upper_bound
    logging.info("Finding epsilon interval for perturbation... done.  Interval is %s", [best_minus_epsilon_lower_bound, best_plus_epsilon_upper_bound])
    return best_minus_epsilon_lower_bound, best_plus_epsilon_upper_bound

def generate_symbolic_continuous(function, components, field=None):
    """
    Construct a vector-space-valued piecewise linear function
    compatible with the given `function`.  Each of the components of
    the function has a slope that is a basis vector of the vector
    space. 
    """
    n = len(components)
    if field is None:
        vector_space = VectorSpace(QQ, n)
    else:
        vector_space = VectorSpace(field, n)
    unit_vectors = vector_space.basis()    
    intervals_and_slopes = []
    for component, slope in itertools.izip(components, unit_vectors):
        intervals_and_slopes.extend([ (interval, slope) for interval in component ])
    intervals_and_slopes.sort()
    if field is None:
        bkpt = [ interval[0] for interval, slope in intervals_and_slopes ] + [1]
        bkpt = nice_field_values(bkpt)
        field = bkpt[0].parent().fraction_field()
    else:
        bkpt = [ field(interval[0]) for interval, slope in intervals_and_slopes ] + [field(1)]
    slopes = [ slope for interval, slope in intervals_and_slopes ]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field)

def generate_additivity_equations_continuous(function, symbolic, field, f=None):
    if f is None:
        f = find_f(function)
    equations = [delta_pi(symbolic, x, y) \
                     for (x, y, z, xeps, yeps, zeps) in generate_additive_vertices(function) ] \
                   + [symbolic(f)] \
                   + [symbolic(1)]
    return matrix(field, equations)

