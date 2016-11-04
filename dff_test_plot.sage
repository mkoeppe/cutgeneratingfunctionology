def find_epsilon_interval_continuous_dff(fn, perturb):
    """Compute the interval [minus_epsilon, plus_epsilon] such that 
    (fn + epsilon * perturb) is superadditive for epsilon in this interval.
    Assumes that fn is superadditive.
    If one of the epsilons is 0, the function bails out early and returns 0, 0.
    """
    logging.info("Finding epsilon interval for perturbation...")
    fn_bkpt = fn.end_points()
    perturb_bkpt = perturb.end_points()
    bkpt_refinement = merge_bkpt(fn_bkpt,perturb_bkpt)
    length1 = len(bkpt_refinement)
    fn_values = []
    perturb_values = []
    for pt in bkpt_refinement:
        fn_values.append(fn(pt))
        perturb_values.append(perturb(pt))
    
    best_minus_epsilon_lower_bound = -10000
    best_plus_epsilon_upper_bound = +10000
    # FIXME: We want to say infinity instead; but bool(SR(2) < infinity) ==> False
    for i in range(length1):
        for j in range(i,length1):
            if bkpt_refinement[i]+ bkpt_refinement[j]<=1:
                a = modified_delta_pi_dff(perturb, perturb_values, bkpt_refinement, i, j)
                if a != 0:
                    b = modified_delta_pi_dff(fn, fn_values, bkpt_refinement, i, j) 
                    if b == 0:
                        logging.info("Zero epsilon encountered for x = %s, y = %s" % (bkpt_refinement[i], bkpt_refinement[j]))
                        return 0, 0 # See docstring
                    epsilon_upper_bound = -b/(abs(a))
                    if a > 0:
                        if -epsilon_upper_bound > best_minus_epsilon_lower_bound:
                            best_minus_epsilon_lower_bound = -epsilon_upper_bound
                    else:
                        if epsilon_upper_bound < best_plus_epsilon_upper_bound:
                            best_plus_epsilon_upper_bound = epsilon_upper_bound

    for i in range(length1):
        for j in range(length1):
            if bkpt_refinement[j] - bkpt_refinement[i] > 0:
                a = modified_delta_pi2_dff(perturb, perturb_values, bkpt_refinement, i, j)
                if a != 0:
                    b = modified_delta_pi2_dff(fn, fn_values, bkpt_refinement, i, j) 

                    if b == 0:
                        logging.info("Zero epsilon encountered for x = %s, y = %s" % (bkpt_refinement[i], bkpt_refinement[j] - bkpt_refinement[i]))
                        return 0, 0 # See docstring
                    epsilon_upper_bound = -b/(abs(a)) 
                    if a > 0:
                        if -epsilon_upper_bound > best_minus_epsilon_lower_bound:
                            best_minus_epsilon_lower_bound = -epsilon_upper_bound
                    else:
                        if epsilon_upper_bound < best_plus_epsilon_upper_bound:
                            best_plus_epsilon_upper_bound = epsilon_upper_bound
    logging.info("Finding epsilon interval for perturbation... done.  Interval is %s", [best_minus_epsilon_lower_bound, best_plus_epsilon_upper_bound])
    return best_minus_epsilon_lower_bound, best_plus_epsilon_upper_bound

def generate_perturbations_equivariant_continuous_dff(fn):
    generator = generate_generic_seeds_with_completion(fn) 
    for seed, stab_int, walk_list in generator:
        perturb = approx_discts_function(walk_list, stab_int, function=fn)
        perturb._seed = seed
        perturb._stab_int = stab_int
        perturb._walk_list = walk_list
        yield perturb

def generate_additivity_equations_dff_continuous(function, symbolic, field):
    bkpt=function.end_points()
    bkpt1=symbolic.end_points()
    #I need to merge these breakpoints, since perturbation function may have different breakpoints
    if function(bkpt[1])==0:
        equations = [delta_pi_for_nonmodule_one(symbolic, x, y) \
                     for (x, y, z, xeps, yeps, zeps) in generate_additive_vertices_dff(function) ] \
                   + [symbolic(bkpt[1])] \
                   + [symbolic(1)]
    else:
        equations = [delta_pi_for_nonmodule_one(symbolic, x, y) \
                     for (x, y, z, xeps, yeps, zeps) in generate_additive_vertices_dff(function) ] \
                   + [symbolic(1)]
    return matrix(field, equations)

def generate_perturbations_finite_dimensional_dff(f):
    covered_intervals=generate_covered_intervals(f)
    uncovered_intervals=generate_uncovered_intervals(f)
    if uncovered_intervals:
        components = copy(covered_intervals)
        components.extend([int] for int in uncovered_intervals)
    else:
        components=covered_intervals
    field = f(0).parent().fraction_field()
    symbolic=generate_symbolic_continuous(f,components,field)
    equation_matrix = generate_additivity_equations_dff_continuous(f, symbolic, field)
    slope_jump_vects = equation_matrix.right_kernel().basis()
    logging.info("Finite dimensional test: Solution space has dimension %s" % len(slope_jump_vects))
    for basis_index in range(len(slope_jump_vects)):
        slope_jump = slope_jump_vects[basis_index]
        perturbation = slope_jump * symbolic
        yield perturbation

def extremality_test_continuous_dff(fn):
    """Still in progress
    """
    covered_intervals=generate_covered_intervals(fn)
    uncovered_intervals=generate_uncovered_intervals(fn)
    if not maximality_test_continuous_dff(fn):
        logging.info("Not maximal, thus NOT extreme.")
        return False
    if uncovered_intervals:
        gen=generate_perturbations_equivariant_continuous_dff(fn)
        for perturbation in gen:
            epsilon_interval = find_epsilon_interval_continuous_dff(fn, perturbation)
            epsilon = min(abs(epsilon_interval[0]), abs(epsilon_interval[1]))
            if epsilon>0:
                logging.info("Not extreme")
                return False
    gen=generate_perturbations_finite_dimensional_dff(fn)
    for perturbation in gen:
        epsilon_interval = find_epsilon_interval_continuous_dff(fn, perturbation)
        epsilon = min(abs(epsilon_interval[0]), abs(epsilon_interval[1]))
        if epsilon>0:
            logging.info("Not extreme")
            return False
    logging.info("extreme")
    return True
    

def generate_additive_vertices_dff(fn):
    return set(itertools.chain( \
                generate_type_1_vertices_for_nonmodule_one(fn, operator.eq),\
                generate_type_2_vertices_for_nonmodule_one(fn, operator.eq)))


def maximality_test_continuous_dff(fn):
    for x in fn.values_at_end_points():
        if (x < 0) or (x > 1):
            logging.info('f is not maximal because it does not stay in the range of [0, 1].')
            return False
    if fn(0) != 0:
        logging.info('f is NOT maximal because f(0) is not equal to 0.')
        return False
    logging.info('f(0) = 0')
    bkpt = fn.end_points()
    if superadditivity_test_continuous(fn) and symmetry_test_dff(fn):
        logging.info('Thus f is maximal.')
        is_maximal = True
    else:
        logging.info('Thus f is NOT maximal.')
        is_maximal = False
    return is_maximal


def symmetry_test_dff(fn):
    bkpt = fn.end_points()
    result=True
    for x in bkpt:
        if fn(x)+fn(1-x)!=1:
            logging.info("pi(%s) + pi(%s) is not equal to 1" % (x, 1-x))
            result= False
    if result:
        logging.info("f is symmetric")
    else:
        logging.info("f is not symmetric")
    return result

def generate_nonsuperadditive_vertices(fn, reduced=True):
    return set(itertools.chain( \
                generate_type_1_vertices_for_nonmodule_one(fn, operator.gt),\
                generate_type_2_vertices_for_nonmodule_one(fn, operator.gt)))

def superadditivity_test_continuous(fn):
    """
    Check if `fn` is superadditive.
    """
    result = True
    for (x, y, z, xeps, yeps, zeps) in generate_nonsuperadditive_vertices(fn, reduced=True):
        logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) > 0" % (x, print_sign(xeps), y, print_sign(yeps), z, print_sign(zeps)))
        result = False
    if result:
        logging.info("f is superadditive.")
    else:
        logging.info("Thus f is not superadditive.")
    return result

def delta_pi_for_nonmodule_one(fn,x,y):
    return fn(x)+fn(y)-fn(x+y)

def modified_delta_pi_dff(fn, fn_values, pts, i, j):
    return fn_values[i] + fn_values[j] - fn(pts[i]+pts[j]) 

def modified_delta_pi2_dff(fn, fn_values2, pts, i, j):
    return fn_values2[i] + fn(pts[j] - pts[i]) - fn_values2[j]  


def generate_type_1_vertices_for_nonmodule_one(fn, comparison):
    """Output 6-tuples (x, y, z,xeps, yeps, zeps).
    """
    bkpt = fn.end_points()
    return ( (x, y, x+y, 0, 0, 0) for x in bkpt for y in bkpt if x+y<=1 and comparison(delta_pi_for_nonmodule_one(fn,x,y), 0) ) 

def generate_type_2_vertices_for_nonmodule_one(fn, comparison):
    bkpt = fn.end_points()
    return ( (x, z-x, z, 0, 0, 0) for x in bkpt for z in bkpt if x <= z and comparison(delta_pi_for_nonmodule_one(fn, x, z-x), 0) ) 

def plot_2d_complex_dff(function):
    """
    Return a plot of the horizonal lines, vertical lines, and diagonal lines of the complex.
    """
    bkpt = function.end_points()
    x = var('x')
    p = Graphics()
    ##kwd = ticks_keywords(function, True)
    ##kwd['legend_label'] = "Complex Delta pi"
    ##plot_kwds_hook(kwd)
    ## We now use lambda functions instead of Sage symbolics for plotting, 
    ## as those give strange errors when combined with our RealNumberFieldElement.
    for i in range(1,len(bkpt)):
        #p += plot(lambda x: bkpt[i]-x, (x, 0, bkpt[i]), color='grey', **kwd)
        p += line([(0,  bkpt[i]), (bkpt[i], 0)], color='grey')
        ##kwd = {}
    for i in range(len(bkpt)-1):
        p += plot(bkpt[i], (0, 1-bkpt[i]), color='grey')
    y=var('y')
    for i in range(len(bkpt)-1):
        p += parametric_plot((bkpt[i],y),(y,0,1-bkpt[i]), color='grey')
    return p

def plot_2d_diagram_dff(fn, show_function=True, known_minimal=False, colorful=False):
    faces = generate_maximal_additive_faces_continuous_dff(fn)
    p = Graphics()
    kwds = { 'legend_label': "Additive face" }
    plot_kwds_hook(kwds)
    if colorful:
        covered_intervals = generate_covered_intervals(fn)
        colors = rainbow(len(covered_intervals))
        interval_color = [(interval[0], i) \
                          for (i, component) in enumerate(covered_intervals) \
                          for interval in component]
        interval_color.sort()
    else:
        covered_intervals = None
    for face in faces:
        if not covered_intervals is None and face.is_2D():
            I = face.minimal_triple[0]
            x = (I[0] + I[1]) / 2
            j = bisect_left(interval_color, (x, len(covered_intervals) + 1)) # should be bisect
            i = interval_color[j-1][1]
            p += face.plot(fill_color = colors[i], **kwds)
        else:
            p += face.plot(**kwds)
        delete_one_time_plot_kwds(kwds)

    ### For non-subadditive functions, show the points where delta_pi is negative.
    if not known_minimal:
        nonsuperadditive_vertices = generate_nonsuperadditive_vertices(fn)
        kwds = { 'legend_label' : "Superadditivity violated" }
        plot_kwds_hook(kwds)
        if fn.is_continuous():
            nonsuperadditive_vertices = {(x,y) for (x, y, z, xeps, yeps, zeps) in nonsuperadditive_vertices}
            p += point(list(nonsuperadditive_vertices),
                       color = "red", size = 50, zorder=-1, **kwds)
            p += point([ (y,x) for (x,y) in nonsuperadditive_vertices ], color = "red", size = 50, zorder=-1)
        else:
            new_legend_label = False
            for (x, y, z, xeps, yeps, zeps) in nonsubadditive_vertices:
                new_legend_label = True
                p += plot_limit_cone_of_vertex(x, y, epstriple_to_cone((xeps, yeps, zeps)))
                if x != y:
                    p += plot_limit_cone_of_vertex(y, x, epstriple_to_cone((yeps, xeps, zeps)))
            if new_legend_label:
                # add legend_label
                p += point([(0,0)], color = "red", size = 50, zorder=-10, **kwds)
                p += point([(0,0)], color = "white", size = 50, zorder=-9)        
    p += plot_2d_complex_dff(fn)
    if show_function:
        p += plot_function_at_borders(fn, covered_components = covered_intervals)
    return p


def generate_maximal_additive_faces_continuous_dff(function):
    logging.info("Computing maximal additive faces...")
    bkpt = function.end_points()     
   
    I_list = []
    J_list = []
    K_list = []
    
    intervals = []
    
    for i in range(len(bkpt)-1):
        intervals.append([bkpt[i],bkpt[i+1]])
    I_list = intervals
    J_list = intervals
    K_list = intervals
    
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
                        else:
                            keep = False
                            logging.warn("Additivity appears only in the interior for some face. This is not shown on the diagram.")
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


