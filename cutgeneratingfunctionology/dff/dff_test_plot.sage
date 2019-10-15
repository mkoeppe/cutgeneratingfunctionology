"""
Classic Dual-Feasible Functions (cDFF):  Maximality and extremality test; plotting.
"""

from six.moves import range

def maximality_test_dff(f):
    """
    Test whether a piecewise linear cDFF is maximal.
    """
    if f.is_continuous():
        return maximality_test_continuous_dff(f)
    else:
        return maximality_test_general_dff(f)

def extremality_test_dff(f):
    """
    Test whether a piecewise linear cDFF is extreme.
    """
    if f.is_continuous():
        return extremality_test_continuous_dff(f)
    else:
        return extremality_test_general_dff(f)

def superadditive_test(f):
    """
    Test whether a piecewise linear cDFF is superadditive.
    """
    if f.is_continuous():
        return superadditivity_test_continuous(f)
    else:
        return superadditivity_test_general(f)



def is_equal(f,h):
    # FIXME: Should fix the problem in FastPiecewise.__eq__ instead.
    bkpt1=f.end_points()
    bkpt2=h.end_points()
    if bkpt1!=bkpt2 or f.limits_at_end_points()!=h.limits_at_end_points():
        return False
    return True

def is_bj(f):
    if not f.is_continuous():
        return False
    bkpt=f.end_points()
    if bkpt[0]==0 and bkpt[1]==1 and f(0)==0 and f(1)==1: 
        logging.info("The function is x")
        return True    
    c=1/bkpt[2]
    g=phi_bj_1(c)
    bkpt1=g.end_points()
    if bkpt!=bkpt1 or f.values_at_end_points()!=g.values_at_end_points():
        return False
    logging.info("The function is phi_bj_1(%s)", c)     
    return True



def pwl_compose(f,g):
    bkpt1=f.end_points()
    bkpt2=g.end_points()
    bkpt=[g.end_points()[0]]
    for i in range(len(bkpt2)-1):
        x=bkpt2[i]
        y=bkpt2[i+1]
        gx=g.limits(x)[1]
        gy=g.limits(y)[-1]
        if gx<gy:
            for j in range(len(bkpt1)):
                if gx<bkpt1[j]<gy:
                    bkpt.append(x+(y-x)*(bkpt1[j]-gx)/(gy-gx))
        if gx>gy:
            for j in range(len(bkpt1)):
                if gx>bkpt1[j]>gy:
                    bkpt.append(x+(y-x)*(bkpt1[j]-gx)/(gy-gx))
        bkpt.append(y)
    pieces=[]
    for i in range(len(bkpt)-1):
        a=bkpt[i]
        b=bkpt[i+1]
        pieces.append([singleton_interval(a),FastLinearFunction(0,f(g(a)))])
        ga=g.limits(a)[1]
        gb=g.limits(b)[-1]
        if ga<gb:
            fga=f.limits(ga)[1]
            fgb=f.limits(gb)[-1]
        if ga>gb:
            fga=f.limits(ga)[-1]
            fgb=f.limits(gb)[1]
        if ga==gb:
            fga=f(ga)
            fgb=fga
        slope=(fgb-fga)/(b-a)
        intercept=fga-slope*a
        pieces.append([open_interval(a,b), FastLinearFunction(slope,intercept)])
    pieces.append([singleton_interval(b),FastLinearFunction(0,f(g(b)))])    
    return FastPiecewise(pieces)    


def find_epsilon_interval_continuous_dff(fn, perturb):
    r"""Compute the interval [minus_epsilon, plus_epsilon] such that 
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

def generate_perturbations_finite_dimensional_continuous_dff(f):
    covered_intervals=generate_covered_intervals(f)
    uncovered_intervals=generate_uncovered_intervals(f)
    if uncovered_intervals:
        components = copy(covered_intervals)
        components.extend([int] for int in uncovered_intervals)
    else:
        components=covered_intervals
    field = f(0).parent().fraction_field()
    symbolic=generate_symbolic_continuous(f,components)
    equation_matrix = generate_additivity_equations_dff_continuous(f, symbolic, field)
    slope_jump_vects = equation_matrix.right_kernel().basis()
    logging.info("Finite dimensional test: Solution space has dimension %s" % len(slope_jump_vects))
    for basis_index in range(len(slope_jump_vects)):
        slope_jump = slope_jump_vects[basis_index]
        perturbation = slope_jump * symbolic
        yield perturbation


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



def extremality_test_continuous_dff(fn, show_plots=False):
    r"""
    DFF extremality test, continuous case.
    """
    if not maximality_test_continuous_dff(fn):
        logging.info("Not maximal, thus NOT extreme.")
        return False
    # Generate the maximal additive faces, which will set the attribute.
    # This is not strictly necessary because, due to the symmetry condition,
    # IGP faces (I, J, K) with K subset of [1,2]
    # map to (1-I, 1-J, 2-K), where K subset of [0,1]. 
    # Hence the covered components are exactly the same as those
    # that would be computed by igp.generate_maximal_additive_faces().
    generate_maximal_additive_faces_dff(fn)
    covered_intervals=generate_covered_intervals(fn)
    uncovered_intervals=generate_uncovered_intervals(fn)
    if uncovered_intervals:
        gen=generate_perturbations_equivariant_continuous_dff(fn)
        for perturbation in gen:
            epsilon_interval = find_epsilon_interval_continuous_dff(fn, perturbation)
            epsilon = min(abs(epsilon_interval[0]), abs(epsilon_interval[1]))
            assert epsilon > 0
            logging.info("Not extreme")
            return False
    gen = generate_perturbations_finite_dimensional_continuous_dff(fn)
    for perturbation in gen:
        epsilon_interval = find_epsilon_interval_continuous_dff(fn, perturbation)
        epsilon = min(abs(epsilon_interval[0]), abs(epsilon_interval[1]))
        assert epsilon > 0
        logging.info("Not extreme")
        return False
    logging.info("extreme")
    return True
    

def generate_additive_vertices_continuous_dff(fn):
    return set(itertools.chain( \
                generate_type_1_vertices_for_nonmodule_one(fn, operator.eq),\
                generate_type_2_vertices_for_nonmodule_one(fn, operator.eq)))

def generate_additive_vertices_dff(fn):
    if fn.is_continuous():
        return generate_additive_vertices_continuous_dff(fn)
    else:
        return generate_additive_vertices_general_dff(fn)


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
    if fn.is_continuous():
        return symmetry_test_continuous_dff(fn)
    else:
        return symmetry_test_general_dff(fn)

def symmetry_test_continuous_dff(fn):
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

def generate_nonsuperadditive_vertices(fn):
    if fn.is_continuous():
        return generate_nonsuperadditive_vertices_continuous(fn)
    else:
        return generate_nonsuperadditive_vertices_general(fn)


def generate_nonsuperadditive_vertices_continuous(fn, reduced=True):
    return set(itertools.chain( \
                generate_type_1_vertices_for_nonmodule_one(fn, operator.gt),\
                generate_type_2_vertices_for_nonmodule_one(fn, operator.gt)))

def superadditivity_test_continuous(fn):
    r"""
    Check if `fn` is superadditive.
    """
    result = True
    for (x, y, z, xeps, yeps, zeps) in generate_nonsuperadditive_vertices_continuous(fn, reduced=True):
        logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) > 0" % (x, print_sign(xeps), y, print_sign(yeps), z, print_sign(zeps)))
        result = False
    if result:
        logging.info("f is superadditive.")
    else:
        logging.info("Thus f is not superadditive.")
    return result

def delta_pi_for_nonmodule_one(fn,x,y):
    return fn(x)+fn(y)-fn(x+y)

def delta_pi_general_dff(fn, x, y, xxx_todo_changeme=(0,0,0)):
    (xeps, yeps, zeps) = xxx_todo_changeme
    return fn.limit(x, xeps) + fn.limit(y, yeps) - fn.limit(x + y, zeps)

def modified_delta_pi_dff(fn, fn_values, pts, i, j):
    return fn_values[i] + fn_values[j] - fn(pts[i]+pts[j]) 

def modified_delta_pi2_dff(fn, fn_values2, pts, i, j):
    return fn_values2[i] + fn(pts[j] - pts[i]) - fn_values2[j]  


def generate_type_1_vertices_for_nonmodule_one(fn, comparison):
    r"""Output 6-tuples (x, y, z,xeps, yeps, zeps).
    """
    bkpt = fn.end_points()
    return ( (x, y, x+y, 0, 0, 0) for x in bkpt for y in bkpt if x+y<=1 and comparison(delta_pi_for_nonmodule_one(fn,x,y), 0) ) 

def generate_type_2_vertices_for_nonmodule_one(fn, comparison):
    bkpt = fn.end_points()
    return ( (x, z-x, z, 0, 0, 0) for x in bkpt for z in bkpt if x <= z and comparison(delta_pi_for_nonmodule_one(fn, x, z-x), 0) ) 

def generate_nonsymmetric_vertices_dff(fn):
    if fn.is_continuous():
        return generate_nonsymmetric_vertices_continuous_dff(fn)
    else:
        return generate_nonsymmetric_vertices_general_dff(fn)

def generate_nonsymmetric_vertices_continuous_dff(fn):
    bkpt=fn.end_points()
    for i in range(len(bkpt)):
        x=bkpt[i]
        y=1-x
        if fn(x)+fn(y) !=1:
            yield(x,y,0,0)


def plot_2d_complex_dff(function,show_tag=True):
    r"""
    Return a plot of the horizonal lines, vertical lines, and diagonal lines of the complex.
    """
    bkpt = function.end_points()
    x = var('x')
    p = Graphics()
    kwd = ticks_keywords(function, True)  # FIXME: This is for IGP and a warning from find_f will be generated
    if show_tag:
        kwd['legend_label'] = "Complex Delta pi"
    plot_kwds_hook(kwd)
    for i in range(1,len(bkpt)):
        #p += plot(lambda x: bkpt[i]-x, (x, 0, bkpt[i]), color='grey', **kwd)
        p += line([(0,  bkpt[i]), (bkpt[i], 0)], color='grey', **kwd)
        kwd = {}
    for i in range(len(bkpt)-1):
        p += plot(bkpt[i], (0, 1-bkpt[i]), color='grey')
    y=var('y')
    for i in range(len(bkpt)-1):
        p += parametric_plot((bkpt[i],y),(y,0,1-bkpt[i]), color='grey')
    return p

def plot_2d_diagram_dff(fn, show_function=True, show_projections=True, known_maximal=False, colorful=False):
    faces = generate_maximal_additive_faces_dff(fn)
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

    ### For non-superadditive functions, show the points where nabla_fn is negative.
    if not known_maximal:
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
            for (x, y, z, xeps, yeps, zeps) in nonsuperadditive_vertices:
                new_legend_label = True
                p += plot_limit_cone_of_vertex(x, y, epstriple_to_cone((xeps, yeps, zeps)))
                if x != y:
                    p += plot_limit_cone_of_vertex(y, x, epstriple_to_cone((yeps, xeps, zeps)))
            if new_legend_label:
                # add legend_label
                p += point([(0,0)], color = "red", size = 50, zorder=-10, **kwds)
                p += point([(0,0)], color = "white", size = 50, zorder=-9)   
        nonsymmetric_vertices = generate_nonsymmetric_vertices_dff(fn)
        kwds = { 'legend_label' : "Symmetry violated" }
        plot_kwds_hook(kwds)
        if fn.is_continuous():
            nonsymmetric_vertices = unique_list((x,y) for (x, y, xeps, yeps) in nonsymmetric_vertices)
            p += point(nonsymmetric_vertices, color = "mediumvioletred", size = 50, zorder=5, **kwds)
            p += point([ (y,x) for (x,y) in nonsymmetric_vertices], color = "mediumvioletred", size = 50, zorder=5)
        else:
            new_legend_label = False
            for (x, y, xeps, yeps) in nonsymmetric_vertices:
                new_legend_label = True
                if (xeps, yeps) == (0, 0):
                    p += point([x, y], color="mediumvioletred", size=20, zorder=5)
                else:
                    p += disk([x, y], 3/100, (yeps* pi/2, (1 - xeps) * pi/2), color="mediumvioletred", zorder=5)
                if x != y:
                    if (xeps, yeps) == (0, 0):
                        p += point([y, x], color="mediumvioletred", size=20, zorder=5)
                    else:
                        p += disk([y, x], 3/100, (xeps* pi/2, (1 - yeps) * pi/2), color="mediumvioletred", zorder=5)
            if new_legend_label:
                # add legend_label
                p += point([(0,0)], color = "mediumvioletred", size = 50, zorder=-10, **kwds)
                p += point([(0,0)], color = "white", size = 50, zorder=-9)     
    p += plot_2d_complex_dff(fn)
    if show_projections:
        p += plot_projections_at_borders(fn)
    if show_function:
        p += plot_function_at_borders(fn, covered_components = covered_intervals)
    return p

def plot_2d_diagram_dff_no_lable(fn, show_function=True, show_projections=True, known_maximal=False, colorful=False):
    # Version of plot_2d_diagram_dff without legend_label.
    faces = generate_maximal_additive_faces_dff(fn)
    p = Graphics()
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
            p += face.plot(fill_color = colors[i])
        else:
            p += face.plot()

    ### For non-superadditive functions, show the points where nabla_fn is negative.
    if not known_maximal:
        nonsuperadditive_vertices = generate_nonsuperadditive_vertices(fn)
        if fn.is_continuous():
            nonsuperadditive_vertices = {(x,y) for (x, y, z, xeps, yeps, zeps) in nonsuperadditive_vertices}
            p += point(list(nonsuperadditive_vertices),
                       color = "red", size = 50, zorder=-1)
            p += point([ (y,x) for (x,y) in nonsuperadditive_vertices ], color = "red", size = 50, zorder=-1)
        else:
            new_legend_label = False
            for (x, y, z, xeps, yeps, zeps) in nonsuperadditive_vertices:
                new_legend_label = True
                p += plot_limit_cone_of_vertex(x, y, epstriple_to_cone((xeps, yeps, zeps)))
                if x != y:
                    p += plot_limit_cone_of_vertex(y, x, epstriple_to_cone((yeps, xeps, zeps)))
            if new_legend_label:
                # add legend_label
                p += point([(0,0)], color = "red", size = 50, zorder=-10)
                p += point([(0,0)], color = "white", size = 50, zorder=-9)   
        nonsymmetric_vertices = generate_nonsymmetric_vertices_dff(fn)
        if fn.is_continuous():
            nonsymmetric_vertices = unique_list((x,y) for (x, y, xeps, yeps) in nonsymmetric_vertices)
            p += point(nonsymmetric_vertices, color = "mediumvioletred", size = 50, zorder=5)
            p += point([ (y,x) for (x,y) in nonsymmetric_vertices], color = "mediumvioletred", size = 50, zorder=5)
        else:
            new_legend_label = False
            for (x, y, xeps, yeps) in nonsymmetric_vertices:
                new_legend_label = True
                if (xeps, yeps) == (0, 0):
                    p += point([x, y], color="mediumvioletred", size=20, zorder=5)
                else:
                    p += disk([x, y], 3/100, (yeps* pi/2, (1 - xeps) * pi/2), color="mediumvioletred", zorder=5)
                if x != y:
                    if (xeps, yeps) == (0, 0):
                        p += point([y, x], color="mediumvioletred", size=20, zorder=5)
                    else:
                        p += disk([y, x], 3/100, (xeps* pi/2, (1 - yeps) * pi/2), color="mediumvioletred", zorder=5)
            if new_legend_label:
                # add legend_label
                p += point([(0,0)], color = "mediumvioletred", size = 50, zorder=-10)
                p += point([(0,0)], color = "white", size = 50, zorder=-9)     
    p += plot_2d_complex_dff(fn,show_tag=False)
    if show_projections:
        p += plot_projections_at_borders_no_lable(fn)
    if show_function:
        p += plot_function_at_borders_no_lable(fn, covered_components = covered_intervals)
    return p


def plot_projections_at_borders_no_lable(fn):
    # Copy of plot_projections_of_faces from functions.sage that does not do legend_label
    r"""
    Plot the projections p1(F), p2(F), p3(F) of all full-dimensional
    additive faces F of `fn` as gray shadows: p1(F) at the top border,
    p2(F) at the left border, p3(F) at the bottom and the right
    borders.
    """
    g = Graphics()
    I_J_verts = set()
    K_verts = set()
    kwds = { 'alpha': proj_plot_alpha, 'zorder': -10 }
    if proj_plot_colors[0] == proj_plot_colors[1] == proj_plot_colors[2]:
        IJK_kwds = [ kwds for i in range(3) ]
    elif proj_plot_colors[0] == proj_plot_colors[1]:
        IJK_kwds = [ kwds, kwds, copy(kwds) ]
    else:
        IJK_kwds = [ copy(kwds) for i in range(3) ]
    for i in range(3):
        #IJK_kwds[i]['legend_color'] = proj_plot_colors[i] # does not work in Sage 5.11
        IJK_kwds[i]['color'] = proj_plot_colors[i]
        plot_kwds_hook(IJK_kwds[i])
    for face in generate_maximal_additive_faces(fn):
        I, J, K = face.minimal_triple
        I_J_verts.update(I) # no need for J because the x-y swapped face will also be processed
        K_verts.update(K)
        g += plot_projections_of_one_face(face, IJK_kwds)
    for (x, y, z, xeps, yeps, zeps) in generate_nonsubadditive_vertices(fn):
        I_J_verts.add(x)
        I_J_verts.add(y)
        K_verts.add(z)
    # plot dashed help lines corresponding to non-breakpoint projections. 
    # (plot_2d_complex already draws solid lines for the breakpoints.)
    I_J_verts.difference_update(fn.end_points())
    for x in I_J_verts:
        g += line([(x, 0), (x, 1)], linestyle=':', color='grey')
        g += line([(0, x), (1, x)], linestyle=':', color='grey')
    K_verts.difference_update(fn.end_points())
    K_verts.difference_update(1 + x for x in fn.end_points())
    for z in K_verts:
        if z <= 1:
            g += line([(0, z), (z, 0)], linestyle=':', color='grey')
        else:
            g += line([(1, z-1), (z-1, 1)], linestyle=':', color='grey')
    return g


def plot_function_at_borders_no_lable(fn, color='blue', covered_components=None):
    # Copy of plot_function_at_borders which does not add 'legend_label'.
    r"""
    Plot the function twice, on the upper and the left border, 
    to decorate 2d diagrams.
    """
    p = Graphics()
    bkpt = fn.end_points()
    limits = fn.limits_at_end_points()
    if not covered_components is None:
        color = 'black'
    if limits[0][0] is not None and limits[0][0] != limits[0][1]:
        p += point([(0,1), (0,0)], color=color, size = 23, zorder=-1)
    for i in range(len(bkpt) - 1):
        x1 = bkpt[i]
        y1 = limits[i][1]
        x2 = bkpt[i+1]
        y2 = limits[i+1][-1]
        y3 = limits[i+1][0]
        y4 = limits[i+1][1]
        if y1 is not None and y2 is not None:
            p += line([(x1, (3/10)*y1 + 1), (x2, (3/10)*y2 + 1)], color=color, zorder=-2)
           
            p += line([(-3/10*y1, x1), (-(3/10)*y2, x2)], color=color, zorder=-2)
        if y1 is not None and limits[i][0] != y1:
            p += point([(x1, (3/10)*y1 + 1), (-(3/10)*y1, x1)], color=color, pointsize=23, zorder=-1)
            p += point([(x1, (3/10)*y1 + 1), (-(3/10)*y1, x1)], color='white', pointsize=10, zorder=-1)
        if y2 is not None and y2 != y3:
            p += point([(x2, (3/10)*y2 + 1), (-(3/10)*y2, x2)], color=color, pointsize=23, zorder=-1)
            p += point([(x2, (3/10)*y2 + 1), (-(3/10)*y2, x2)], color='white', pointsize=10, zorder=-1)
        if y3 is not None and ((y2 != y3) or ((i < len(bkpt) - 2) and (y3 != y4))) and \
                              ((i == len(bkpt)-2) or not (y3 == y4 and y2 is None) and \
                                                     not (y2 == y3 and y4 is None)):
            p += point([(x2, (3/10)*y3 + 1), (-(3/10)*y3, x2)], color=color, pointsize=23, zorder=-1)
    # plot function at borders with different colors according to slope values.
    p += plot_covered_components_at_borders(fn, covered_components)
    # add legend_label
    if fn.is_discrete():
        p += point([(0,0)], color=color, pointsize=23, zorder=-10)
        p += point([(0,0)], color='white', pointsize=23, zorder=-9)
    else:
        p += line([(0,0), (0,1)], color=color, zorder=-10)
        p += line([(0,0), (0,1)], color='white', zorder=-9)
    return p

def plot_2d_diagram_with_cones_dff(fn, show_function=True):
    g=plot_2d_complex_dff(fn)
    if show_function:
        g += plot_function_at_borders(fn)
    bkpt = fn.end_points()
    kwds = { 'legend_label' : "Superadditivity violated" }
    plot_kwds_hook(kwds)
    type_1_vertices = [(x, y, x+y) for x in bkpt for y in bkpt if x <= y and x+y<=1]
    type_2_vertices = [(x, z-x, z) for x in bkpt for z in bkpt if x < z]
    vertices = set(type_1_vertices + type_2_vertices)
    if fn.is_continuous():
        for (x, y, z) in vertices:
            deltafn = delta_pi_for_nonmodule_one(fn, x, y)
            if deltafn < 0:
                color = "white"
            elif deltafn == 0:
                color = "mediumspringgreen"
            else:
                color = "red"
            g += point([(x, y), (y, x)], color=color, size = 200, zorder=-1, **kwds)
    else:
        for (x, y, z) in vertices:
            for (xeps, yeps, zeps) in [(0,0,0)]+list(nonzero_eps):
                if (x==0 and xeps==-1) or (y==0 and yeps==-1) or (z==0 and zeps==-1):
                    continue
                if (x==1 and xeps==1) or (y==1 and yeps==1) or (z==1 and zeps==1):
                    continue
                deltafn = delta_pi_general_dff(fn, x, y, (xeps, yeps, zeps))
                if deltafn < 0:
                    color = "white"
                elif deltafn == 0:
                    color = "mediumspringgreen"
                else:
                    color = "red"
                g += plot_limit_cone_of_vertex(x, y, epstriple_to_cone((xeps, yeps, zeps)), color=color, r=0.03)
                g += plot_limit_cone_of_vertex(y, x, epstriple_to_cone((yeps, xeps, zeps)), color=color, r=0.03)
    kwds = { 'legend_label' : "Symmetry violated" }
    plot_kwds_hook(kwds)
    for i in range(len(bkpt)):
        x=bkpt[i]
        y=1-x
        if fn(x)+fn(y) !=1:
            g += point([(x, y), (y, x)], color="mediumvioletred", size = 200, zorder=-1, **kwds) 
        if not fn.is_continuous():
            limits_x = fn.limits(x)
            limits_y = fn.limits(y)
            if limits_x[-1] + limits_y[1] != 1:
                g += plot_limit_cone_of_vertex(x, y, epstriple_to_cone((-1, 1, 0)), color="mediumvioletred", r=0.03)
                g += plot_limit_cone_of_vertex(y, x, epstriple_to_cone((1, -1, 0)), color="mediumvioletred", r=0.03)
            if limits_x[1] + limits_y[-1] != 1:
                g += plot_limit_cone_of_vertex(x, y, epstriple_to_cone((1, -1, 0)), color="mediumvioletred", r=0.03)
                g += plot_limit_cone_of_vertex(y, x, epstriple_to_cone((-1, 1, 0)), color="mediumvioletred", r=0.03)
    return g

def generate_maximal_additive_faces_dff(fn):
    if hasattr(fn, '_maximal_additive_faces'):
        return fn._maximal_additive_faces
    if fn.is_continuous():
        result = generate_maximal_additive_faces_continuous_dff(fn)
    else:
        result = generate_maximal_additive_faces_general_dff(fn)
    # FIXME: The following line is crucial.
    # Various functions from the IGP module are reused
    # and will only give correct results for DFF if this attribute is precomputed!
    fn._maximal_additive_faces = result
    return result


def generate_maximal_additive_faces_continuous_dff(function):
    """
    EXAMPLES::

        sage: from cutgeneratingfunctionology.dff import *
        sage: logging.disable(logging.INFO)   # Suppress output in automatic tests.
        sage: h=phi_bj_1(3/2)
        sage: generate_maximal_additive_faces_continuous_dff(h)
        [<Face ([0, 1/3], [0, 1/3], [0, 1/3])>,
         <Face ([0], [1/3, 2/3], [1/3, 2/3])>,
         <Face ([1/3, 2/3], [0], [1/3, 2/3])>,
         <Face ([0, 1/3], [2/3, 1], [2/3, 1])>,
         <Face ([2/3, 1], [0, 1/3], [2/3, 1])>,
         <Face ([1/3, 2/3], [1/3, 2/3], [1])>]
    """
    logging.info("Computing maximal additive faces...")
    bkpt = function.end_points()     
    intervals = []
    for i in range(len(bkpt)-1):
        intervals.append([bkpt[i],bkpt[i+1]])
    I_list = J_list = K_list = intervals
    
    additive_face = {}
    additive_vertices = {}
    faces = []
    for i in range(len(I_list)):
        for j in range(i, len(J_list)):
            IplusJ = interval_sum(I_list[i],J_list[j])
            for k in generate_overlapping_interval_indices(IplusJ, bkpt):
                # Check if int(I+J) intersects int(K) is non-empty.
                if len(interval_intersection(IplusJ, K_list[k])) == 2:
                    temp_verts = verts(I_list[i],J_list[j],K_list[k])
                    temp = []
                    keep = False
                    if temp_verts != []:
                        for vertex in temp_verts:
                            if delta_pi(function, vertex[0],vertex[1]) == 0:
                                temp.append(vertex)
                                keep = True
                    if len(temp) == 2:
                        # edge
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
                            # diagonal
                            if temp[1][0] + temp[1][1] == K_list[k][0]:
                                if (i,j,k-1) in additive_face:
                                    keep = False
                            elif temp[1][0] + temp[1][1] == 1:
                                # Part of the big diagonal edge!
                                pass
                            else:
                                keep = False
                        else:
                            keep = False
                            logging.warning("Additivity appears only in the interior for some face. This is not shown on the diagram.")
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


