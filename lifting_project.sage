if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *


def lift_on_uncovered_components(fn, show_plots=False, L=[]):
    """
    EXAMPLES::

        sage: logging.disable(logging.info)
        sage: h = drlm_not_extreme_1()
        sage: hh = lift_on_uncovered_components(h)
        sage: len(hh)
        4
        sage: h = example7slopecoarse2()
        sage: hh = lift_on_uncovered_components(h)
        sage: len(hh)
        12
        sage: bkpts = [0, 1/13, 3/13, 7/26, 4/13,5/13,21/52,23/52,6/13,8/13,33/52,35/52,9/13,10/13,21/26,11/13,1]
        sage: values = [0,1,3/14,5/7,3/4,5/14,55/112,33/112,3/7,4/7,79/112,57/112,9/14,1/4,2/7,11/14,0]
        sage: h = piecewise_function_from_breakpoints_and_values(bkpts, values)
        sage: hh = lift_on_uncovered_components(h)   # long time. 30 mins?
        sage: len(hh)                                # long time
        128
        sage: extremality_test(hh[0])
        True
    """
    return list(generate_extreme_lifted_function_equiv(fn, show_plots=show_plots, L=L, use_polyhedron=True))

def generate_extreme_lifted_function_equiv(fn, show_plots=False, L=[], use_polyhedron=True):
    find_decomposition_into_stability_intervals_with_completion(fn)
    if L:
        perturbation_components = [fn._stability_orbits[i] for i in L]
        gen_pert = generate_perturbation_on_perturbation_components(fn, perturbation_components=perturbation_components, use_polyhedron=use_polyhedron)
        for perturbation in gen_pert:
            lifted = fn + perturbation
            if extremality_test(lifted):
                if show_plots:
                    p = plot(perturbation, color='magenta')
                    p += plot(fn, color='black')
                    p += plot_covered_intervals(lifted)
                    p.show()
                yield(lifted)
    else:
        ### Loop over all subsets L.
        #stability_orbits_gen = list(powerset(fn._stability_orbits))
        #for perturbation_components in stability_orbits_gen[:0:-1]: # card descreasing order, discard empty set.
        for k in range(len(fn._stability_orbits) + 1,0,-1):
            for perturbation_components in itertools.combinations(fn._stability_orbits, k):
                gen_pert = generate_perturbation_on_perturbation_components(fn, perturbation_components=perturbation_components, use_polyhedron=use_polyhedron)
                for perturbation in gen_pert:
                    lifted = fn + perturbation
                    if extremality_test(lifted):
                        if show_plots:
                            p = plot(perturbation, color='magenta')
                            p += plot(fn, color='black')
                            p += plot_covered_intervals(lifted)
                            p.show()
                        yield(lifted)

# def check_lifted_function_is_extreme(fn, perturbation):
#     lifted = fn + perturbation
#     # epsilon_interval = find_epsilon_interval(fn, perturbation)
#     # if not epsilon_interval[1] == 1:
#     #     #print "epsilon interval = (%s, %s) is not 1." % epsilon_interval
#     #     #print sage_input(fn)
#     #     #print sage_input(perturbation)
#     #     return False
#     # if generate_uncovered_components(lifted):
#     #     #print "The lifted function has uncovered intervals"
#     #     #print sage_input(fn)
#     #     #print sage_input(perturbation)
#     #     return False
#     if not extremality_test(lifted):
#         # print "This perturbation does not give extreme lifted function"
#         # print sage_input(fn)
#         # print sage_input(perturbation)
#         return False
#     return True

def generate_perturbation_on_perturbation_components(fn, perturbation_components=None, use_polyhedron=True, repeat=None):
    if perturbation_components is None:
        find_decomposition_into_stability_intervals_with_completion(fn)
        if not fn._stability_orbits:
            return
        perturbation_components = fn._stability_orbits
    num_perturbation_components = len(perturbation_components)
    if repeat is None:
        repeat = [1] * num_perturbation_components
    gen_lim_slopes = generate_lim_slopes_on_perturbation_components(fn, perturbation_components=perturbation_components, solver=None, use_polyhedron=use_polyhedron)
    for lim_slope in gen_lim_slopes:
        perturbation = zero_perturbation_partial_function([[(0,1)]],[]) #piecewise_function_from_breakpoints_and_values([0,1],[0,0])
        for i in range(num_perturbation_components):
            orbit = perturbation_components[i][0]
            int = orbit[0]
            seed = (int.a + int.b) / 2
            stab_int = closed_or_open_or_halfopen_interval(int.a - seed, int.b - seed, int.left_closed, int.right_closed)
            walk_dict = perturbation_components[i][1]
            s_left= lim_slope[2*i]
            s_right = lim_slope[2*i+1]
            if s_left == s_right == 0:
                continue
            if s_left * s_right < 0:
                x = s_right/(s_right - s_left)
                perturbation_template_bkpts_once = [0, x]
                perturbation_template_values_once = [0, interval_length(int) * (s_left * s_right)/(s_right - s_left)]
            else:
                s_fn = (fn(int.b) - fn(int.a)) / (int.b - int.a)
                s_plus, s_minus = limiting_slopes(fn)
                if s_left < 0:
                    s_mid = s_plus - s_fn
                else:
                    s_mid = s_minus - s_fn
                x = s_mid / (2*s_mid - s_left - s_right)
                perturbation_template_bkpts_once = [0, x, 1-x]
                perturbation_template_values_once = [0, interval_length(int) * x * s_left, -interval_length(int) * x * s_right]
            global perturbation_template_bkpts
            global perturbation_template_values
            perturbation_template_bkpts = []
            perturbation_template_values = []
            for j in range(repeat[i]):
                perturbation_template_bkpts += [j/repeat[i]+x/repeat[i] for x in perturbation_template_bkpts_once]
                perturbation_template_values += [x/repeat[i] for x in perturbation_template_values_once]
            perturbation_template_bkpts += [1]
            perturbation_template_values += [0]
            #igp.perturbation_template_bkpts = perturbation_template_bkpts  # for debug
            #igp.perturbation_template_values = perturbation_template_values # for debug
            #print perturbation_template_bkpts, perturbation_template_values, walk_dict, sage_input(stab_int)
            pert_i = approx_discts_function(walk_dict, stab_int, function=fn)
            #print pert_i
            perturbation = pert_i + perturbation
        yield perturbation
    
def generate_lim_slopes_on_perturbation_components(fn, perturbation_components=None, solver=None, use_polyhedron=True):
    mip = perturbation_lim_slopes_mip(fn, perturbation_components=perturbation_components, solver=solver)
    if use_polyhedron:
        vertices = mip.polyhedron().vertices()
    else:
        vertices = generate_random_mip_sol(mip)
    for vertex in vertices:
        good_lim_slopes = True
        for i in range(0, len(vertex), 2):
            sa = vertex[i]
            sb = vertex[i+1]
            #if (sa == 0 and sb != 0) or (sa != 0 and sb == 0):
            #if not (sa * sb < 0 or sa == 0 and sb == 0):
            #if sa * sb >= 0:
            if sa * sb == 0:
                good_lim_slopes = False
                #print "This pert_lim_slopes = %s is ineligible." % vertex
                break
        if good_lim_slopes:
            yield tuple(vertex)
    
def perturbation_lim_slopes_mip(fn, perturbation_components=None, solver=None):
    if perturbation_components is None:
        find_decomposition_into_stability_intervals_with_completion(fn)
        perturbation_components = fn._stability_orbits
    num_perturbation_components = len(perturbation_components)
    pert_lim_slopes = {1:{}, -1:{}, 0:{}}
    mip = MixedIntegerLinearProgram(solver=solver, base_ring=QQ)
    for i in range(num_perturbation_components):
        orbit = perturbation_components[i][0]
        walk_dict = perturbation_components[i][1]
        sa = mip[2*i]
        sb = mip[2*i+1]
        for interval in orbit:
            seed = (interval.a + interval.b) / 2
            sign = walk_dict[seed][0]
            if sign == 1:
                pert_lim_slopes[1][interval.a] = sa
                pert_lim_slopes[-1][interval.b] = sb
            else:
                pert_lim_slopes[1][interval.a] = sb
                pert_lim_slopes[-1][interval.b] = sa
    bkpt = uniq(fn.end_points()+pert_lim_slopes[1].keys()+pert_lim_slopes[-1].keys())
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt]
    slope2 = [(fn(fractional(bkpt2[i+1])) - fn(fractional(bkpt2[i]))) / (bkpt2[i+1] - bkpt2[i]) for i in range(len(bkpt2)-1)]
    fn_lim_slopes = {1:{}, -1:{}, 0:{}}
    type_1_vertices = ((x, y, x+y) for x in bkpt[:-1] for y in bkpt[:-1] if x <= y)
    type_2_vertices = ((x, z-x, z) for x in bkpt[:-1] for z in bkpt2[:-1] if x < z < 1+x)
    additive_vertices = unique_list((x,y,fractional(z)) for (x,y,z) in itertools.chain(type_1_vertices,type_2_vertices) if delta_pi(fn,x,y)==0)
    for (x, y, z) in additive_vertices:
        for w in [x,y,z]:
            if w in fn_lim_slopes[-1]:
                continue
            i = bisect_left(bkpt, w)
            if bkpt[i] != w:
                fn_lim_slopes[-1][w] = fn_lim_slopes[1][w] = slope2[i-1]
            else: 
                fn_lim_slopes[-1][w] = slope2[i-1]
                fn_lim_slopes[1][w] = slope2[i]
    constraints_set = set([])
    for (x, y, z) in additive_vertices:
        for (xeps, yeps, zeps) in [(1,0,1), (-1,0,-1), (0,1,1), (0,-1,-1), (1,-1,0), (-1,1,0)]:
            deltafn = fn_lim_slopes[xeps].get(x, 0)*xeps + fn_lim_slopes[yeps].get(y, 0)*yeps - fn_lim_slopes[zeps].get(z, 0)*zeps
            deltap = pert_lim_slopes[xeps].get(x, 0)*xeps + pert_lim_slopes[yeps].get(y, 0)*yeps - pert_lim_slopes[zeps].get(z, 0)*zeps
            if deltap in QQ:
                continue
            deltap_coef = tuple(deltap.coefficient(mip[i]) for i in range(num_perturbation_components*2))
            if all(coef == 0 for coef in deltap_coef):
                continue
            constraint_coef = tuple([deltafn]) + deltap_coef
            if constraint_coef in constraints_set:
                # don't add duplicated constraints.
                continue
            else:
                constraints_set.add(constraint_coef)
            #print (x,y,z), (xeps,yeps,zeps), deltafn + deltap
            mip.add_constraint(deltafn + deltap >= 0)
    return mip
