def generate_animation_2d_diagram(fn, show_plots=True, f=None):
    """
    sage: fn = gmic() # not tested
    sage: fname="./animation/gmic-%s.pdf" # not tested
    sage: plots = generate_animation_2d_diagram(fn, show_plots=fname) # not tested

    sage: fn = drlm_backward_3_slope(1/8,2/8)
    sage: h = interpolate_to_infinite_group(restrict_to_finite_group(fn), merge=False)
    sage: fname="./animation/drlm_backward_3_slope_res-%s.pdf"
    sage: plots = generate_animation_2d_diagram(h, show_plots=fname)
    """
    if f is None:
        f = find_f(fn, no_error_if_not_minimal_anyway=True)
    num_of_colors = len(generate_covered_intervals(fn))
    colors = permute_rainbow(num_of_colors)
    faces = generate_maximal_additive_faces(fn)
    plots=[]
    ### plot the grey 2d-diagram
    g0 = plot_2d_complex(fn)
    kwds = { 'legend_label': "Additive face" }
    plot_kwds_hook(kwds)
    for face in faces:
        g0 += face.plot(rgbcolor="grey", fill_color="grey", **kwds)
        delete_one_time_plot_kwds(kwds)
    #g0 += plot_projections_at_borders(fn)
    g0 += plot_function_at_borders(fn, color='black')
    g0.SHOW_OPTIONS['show_legend']=False
    IJK_kwds =  [{ 'alpha': 0, 'zorder': -10 , 'color': 'white'} for i in range(3)]
    phantom_projection= plot_projections_of_one_face(faces[0], IJK_kwds)+plot_projections_of_one_face(faces[-1], IJK_kwds)
    g0 += phantom_projection
    g = g0
    step = 0
    show_plot(g0, show_plots, tag=step , object=fn, show_legend=False, xmin=-0.3, xmax=1.02, ymin=-0.02, ymax=1.3)
    plots.append(g0)
    ### direct covered intervals
    colored_faces={}
    covered_intervals = []
    for face in faces:
        if face.is_2D():
            component = []
            for int1 in face.minimal_triple:
                component.append(interval_mod_1(int1))
            component.sort()
            component = merge_within_comp(component)
            covered_intervals.append(component)
            j = len(covered_intervals)-1
            colored_faces[j]= [face]
            g += face.plot(fill_color=colors[j], **kwds)
            IJK_kwds =  [{ 'alpha': 0.35, 'zorder': -10 , 'color': 'grey'} for i in range(3)]
            p_face = plot_projections_of_one_face(face, IJK_kwds)
            g += p_face
            step += 1
            show_plot(g, show_plots, tag=step , object=fn, show_legend=False, xmin=-0.3, xmax=1.02, ymin=-0.02, ymax=1.3)
            plots.append(g)
            for i in range(j):
                if find_interior_intersection(covered_intervals[i], covered_intervals[j]):
                    covered_intervals[j] = merge_two_comp(covered_intervals[i],covered_intervals[j])
                    covered_intervals[i] = []
                    colored_faces[j].extend(colored_faces[i])
                    colored_faces[i] = []
            g0 += p_face
            g = g0
            g += plot_colored_faces(colored_faces, colors)
            g += plot_colored_segments_at_borders(fn, covered_intervals, colors, **kwds)
            beams = plot_beams_of_one_face(face)
            step += 1
            show_plot(g+beams, show_plots, tag=step, object=fn, show_legend=False, xmin=-0.3, xmax=1.02, ymin=-0.02, ymax=1.3)
            plots.append(g+beams)
    ### undirect covered intervals
    edges = [ face.minimal_triple for face in faces if face.is_1D()]
    any_change = True
    while any_change:
        any_change = False
        for edge in edges:
            intervals = []
            # 0 stands for I; 1 stands for J; 2 stands for K
            IJK = []
            for i in range(len(edge)):
                if len(edge[i]) == 2:
                    intervals.append(edge[i])
                    IJK.append(i)
            if edge_merge(covered_intervals,intervals,IJK):
                any_change = True
                face = Face(edge)
                colored_edge = face.plot(rgbcolor="cyan", **kwds)
                g += colored_edge
                step += 1
                show_plot(g, show_plots, tag=step, object=fn, show_legend=False, xmin=-0.3, xmax=1.02, ymin=-0.02, ymax=1.3)
                plots.append(g)
                g = g0
                colored_faces = generate_colored_faces(faces, covered_intervals)
                g += plot_colored_faces(colored_faces, colors)
                g += plot_colored_segments_at_borders(fn, covered_intervals, colors, **kwds)
                beams = plot_beams_of_one_face(face)
                step += 1
                show_plot(g+beams+colored_edge, show_plots, tag=step, object=fn, show_legend=False, xmin=-0.3, xmax=1.02, ymin=-0.02, ymax=1.3)
                plots.append(g+beams+colored_edge)

    step += 1
    show_plot(g, show_plots, tag=step, object=fn, show_legend=False, xmin=-0.3, xmax=1.02, ymin=-0.02, ymax=1.3)
    plots.append(g)
    #print len(plots)
    #if isinstance(show_plots, str):
    #    animate(plots).save(show_plots % "animation" + ".gif", delay=100)
    return plots

def permute_rainbow(num_of_colors):
    n = int(num_of_colors / 3)
    while gcd(n, num_of_colors) != 1:
        n -= 1
    colors = [rainbow(num_of_colors)[(i * n) % num_of_colors] for i in range(num_of_colors)]
    return colors

def plot_beams_of_one_face(face):
    g = Graphics()
    I, J, K = face.minimal_triple
    if len(I) == 2:
        g += polygon([(I[0], 1), (I[1], 1), (I[1], 0), (I[0], 0)], color='yellow', fill=True, alpha=0.35, zorder=-5)
    if len(J) == 2:
        g += polygon([(0, J[0]), (0, J[1]), (1, J[1]), (1, J[0])], color='yellow', fill=True, alpha=0.35, zorder=-5)
    if len(K) == 2:
        if coho_interval_contained_in_coho_interval(K, [0,1]):
            g += polygon([(K[0], 0), (K[1], 0), (0, K[1]), (0, K[0])], color='yellow', fill=True, alpha=0.35, zorder=-5)
        elif coho_interval_contained_in_coho_interval(K, [1,2]):
            g += polygon([(1, K[0]-1), (1, K[1]-1), (K[1] - 1, 1), (K[0] - 1, 1)], color='yellow', fill=True, alpha=0.35, zorder=-5)
        else:
            raise ValueError, "Bad face: %s" % face
    return g

def plot_colored_faces(colored_faces, colors):
    g = Graphics()
    for (i, faces) in colored_faces.items():
        for face in faces:
            g += face.plot(fill_color=colors[i])
    return g

def generate_colored_faces(faces, covered_intervals):
    colored_faces = {i:[] for i in range(len(covered_intervals))}
    for face in faces:
        if face.is_2D():
            I = face.minimal_triple[0]
            x = (I[0] + I[1]) / 2
            j = -1
            for i, component in enumerate(covered_intervals):
                for interval in component:
                    if interval[0] <= x <= interval[1]:
                        j = i
                        break
                if j != -1:
                    break
            colored_faces[j].append(face)
    return colored_faces

def plot_colored_segments_at_borders(fn, covered_intervals, colors, **kwds):
    g = Graphics()
    for i, component in enumerate(covered_intervals):
        for interval in component:
            x1 = interval[0]
            x2 = interval[1]
            y1 = fn.limit(x1, 1)
            y2 = fn.limit(x2, -1)
            g += line([(x1, 0.3*y1 + 1), (x2, 0.3*y2 + 1)], color=colors[i], zorder=-2, **kwds)
            g += line([(-0.3*y1, x1), (-0.3*y2, x2)], color=colors[i], zorder=-2, **kwds)
    return g
