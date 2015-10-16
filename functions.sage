# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *


import itertools

def fractional(num):
    """
    Reduce a number modulo 1.
    """
    parent = num.parent()
    one = parent._one_element
    zero = parent._zero_element
    while num > one:
        num = num - one
    while num < zero:
        num = num + one
    return num

def delta_pi(fn,x,y):
    """
    Compute the slack in subaddivity.
    """
    return fn(fractional(x))+fn(fractional(y))-fn(fractional(x+y))

def plot_2d_complex(function):
    """
    Return a plot of the horizonal lines, vertical lines, and diagonal lines of the complex.
    """
    bkpt = function.end_points()
    x = var('x')
    p = Graphics()
    kwd = ticks_keywords(function, True)
    kwd['legend_label'] = "Complex Delta pi"
    plot_kwds_hook(kwd)
    ## We now use lambda functions instead of Sage symbolics for plotting, 
    ## as those give strange errors when combined with our RealNumberFieldElement.
    for i in range(1,len(bkpt)):
        #p += plot(lambda x: bkpt[i]-x, (x, 0, bkpt[i]), color='grey', **kwd)
        p += line([(0,  bkpt[i]), (bkpt[i], 0)], color='grey', **kwd)
        kwd = {}
    for i in range(1,len(bkpt)-1):
        #p += plot(lambda x: (1+bkpt[i]-x), (x, bkpt[i], 1), color='grey')
        p += line([(bkpt[i], 1), (1, bkpt[i])], color='grey')
    for i in range(len(bkpt)):
        p += plot(bkpt[i], (0, 1), color='grey')
    y=var('y')
    for i in range(len(bkpt)):
        p += parametric_plot((bkpt[i],y),(y,0,1), color='grey')
    return p

##
##
##

def projection(vertices,linear_form):
    """
    Compute the projection of vertices based on the linear form.
    vertices is a list of vertices (2-tuples)
    linear_form is a 2-element list.
    Projection on x: [1,0]
    Projection on y: [0,1]
    Projection on x + y: [1,1]
    """
    temp = []
    for i in vertices:
        temp.append(i[0]*linear_form[0]+i[1]*linear_form[1])
    if max(temp) == min(temp):
        return [min(temp)]
    else:
        return [min(temp), max(temp)]

def projections(vertices):
    """
    Compute F(I,J,K)            
    """
    return [projection(vertices, [1,0]),projection(vertices, [0,1]),projection(vertices, [1,1])]    

def verts(I1, J1, K1):
    """
    Compute the vertices based on I, J, and K.        
    """
    temp = []
    for i in I1:
        for j in J1:
            if element_of_int(i+j,K1):
                temp.append((i,j))
    for i in I1:
        for k in K1:
            if element_of_int(k-i,J1) and (i,k-i) not in temp:
                temp.append((i,k-i))             
    for j in J1:
        for k in K1:
            if element_of_int(k-j,I1) and (k-j,j) not in temp:
                temp.append((k-j,j))
    
    if len(temp) > 0:
        return temp

def generate_maximal_additive_faces(fn):
    if hasattr(fn, '_maximal_additive_faces'):
        return fn._maximal_additive_faces
    if fn.is_discrete():
        result = generate_maximal_additive_faces_discrete(fn)
    elif fn.is_continuous():
        result = generate_maximal_additive_faces_continuous(fn)
    else:
        result = generate_maximal_additive_faces_general(fn)
    fn._maximal_additive_faces = result
    return result
            
### Create a new class representing a "face" (which knows its
### vertices, minimal triple, whether it's a translation/reflection,
### etc.; whether it's solid or dense).
class Face:
    def __init__(self, triple, vertices=None, is_known_to_be_minimal=False):
        """
        EXAMPLES::

            sage: logging.disable(logging.INFO)
            sage: f = generate_maximal_additive_faces(bhk_irrational(delta=(23/250,1/125)))
        """
        if not vertices:
            vertices = verts(triple[0], triple[1], triple[2])
            if not vertices:
                raise NotImplementedError, "An empty face. This could mean need to shift triple[2] by (1,1). Not implemented."
        self.vertices = vertices
        i, j, k = projections(vertices)
        self.minimal_triple = minimal_triple = (i, j, k)
        #self._warned_about_non_minimal_triple = False
        #if is_known_to_be_minimal and not triples_equal(minimal_triple, triple) and not self._warned_about_non_minimal_triple:
        #    logging.warn("Provided triple was not minimal: %s reduces to %s" % (triple, minimal_triple))
        #    self._warned_about_non_minimal_triple = True
            # FIXME: Check why (i,j,k) != (i,j,k+1) can happen.

    def __repr__(self):
        return '<Face ' + repr(self.minimal_triple) + '>'

    def plot(self, rgbcolor=(0.0 / 255.0, 250.0 / 255.0, 154.0 / 255.0), fill_color="mediumspringgreen", *args, **kwds):
        y = var('y')
        trip = self.minimal_triple
        vert = self.vertices
        if self.is_0D():
            return point((trip[0][0], \
                          trip[1][0]), rgbcolor = rgbcolor, size = 30, **kwds)
        elif self.is_horizontal():
            return line([(trip[0][0],trip[1][0]),(trip[0][1],trip[1][0])], rgbcolor = rgbcolor, thickness=2, **kwds)
        elif self.is_vertical():
            return line([(trip[0][0],trip[1][0]),(trip[0][0],trip[1][1])], rgbcolor = rgbcolor, thickness=2, **kwds)
        elif self.is_diagonal():
            return line([(trip[0][0],trip[2][0]-trip[0][0]),(trip[0][1],trip[2][0]-trip[0][1])], rgbcolor = rgbcolor, thickness=2, **kwds)
        elif self.is_2D():
            ## Sorting is necessary for this example:
            ## plot_2d_diagram(lift(piecewise_function_from_robert_txt_file("data/dey-richard-not-extreme.txt"))
            return polygon(convex_vert_list(vert), color=fill_color, **kwds)

    def is_directed_move(self):
        return self.is_1D()
        
    def directed_move_with_domain_and_codomain(self):
        """
        Maps a horizontal edge to a forward translation,
        a vertical edge to a backward translation, 
        a diagonal edge to a reflection.

        `domain` and `codomain` are lists of open intervals.
        Endpoints of `domain` and `codomain` will be taken care of by additive vertices.
        """
        (I, J, K) = self.minimal_triple
        if self.is_horizontal():
            K_mod_1 = interval_mod_1(K)
            t = K_mod_1[0] - I[0]
            return (1, t), [open_interval(* I)], [open_interval(* K_mod_1)]
        elif self.is_vertical():
            K_mod_1 = interval_mod_1(K)
            t = K_mod_1[0] - J[0]
            return (1, t), [open_interval(* J)], [open_interval(* K_mod_1)]
        elif self.is_diagonal():
            # attention: I, J are not sorted.
            return (-1, K[0]), [open_interval(* I), open_interval(* J)], [open_interval(* J), open_interval(* I)]
        else:
            raise ValueError, "Face does not correspond to a directed move: %s" % self

    def functional_directed_move(self):
        directed_move, domain, codomain = self.directed_move_with_domain_and_codomain()
        fdm = FunctionalDirectedMove(domain, directed_move)
        return fdm

    def is_0D(self):
        return len(self.vertices) == 1

    def is_1D(self):
        return len(self.vertices) == 2

    def is_2D(self):
        return len(self.vertices) > 2

    def is_horizontal(self):
        return self.is_1D() and self.vertices[0][1] == self.vertices[1][1]

    def is_vertical(self):
        return self.is_1D() and self.vertices[0][0] == self.vertices[1][0]

    def is_diagonal(self):
        return self.is_1D() and \
               self.vertices[0][0] + self.vertices[0][1] == self.vertices[1][0] + self.vertices[1][1]
    
def plot_faces(faces, **kwds):
    p = Graphics()
    for f in faces:
        if type(f) == list or type(f) == tuple: #legacy
            f = Face(f)
        p += f.plot(**kwds)
    return p

def plot_trivial_2d_diagram_with_grid(function, xgrid=None, ygrid=None): 
    """
    Return a plot of the 2d complex with vertices marked that 
    have delta_pi == 0.  Does not use any complicated code.
    Mainly used for visually double-checking the computation of 
    maximal additive faces.
    """
    if xgrid is None:
        xgrid = function.end_points()
    if ygrid is None:
        ygrid = function.end_points()
    return point([(x,y) for x in xgrid for y in ygrid \
                  if delta_pi(function, x, y) == 0],
                 color="cyan", size = 80)

def angle_cmp(a, b, center):
    # Adapted 
    # from http://stackoverflow.com/questions/6989100/sort-points-in-clockwise-order
    if a[0] - center[0] >= 0 and b[0] - center[0] < 0:
        return int(1)
    elif a[0] - center[0] < 0 and b[0] - center[0] >= 0:
        return int(-1)
    elif a[0] - center[0] == 0 and b[0] - center[0] == 0:
        return cmp(a[1], b[1])

    det = (a[0] - center[0]) * (b[1] - center[1]) - (b[0] - center[0]) * (a[1] - center[1])
    if det < 0:
        return int(1)
    elif det > 0:
        return int(-1)

    return int(0)

import operator

def convex_vert_list(vertices):
    if len(vertices) <= 3:
        return vertices
    else:
        center = reduce(operator.add, map(vector, vertices)) / len(vertices)
        return sorted(vertices, cmp = lambda a,b: angle_cmp(a, b, center))

def plot_kwds_hook(kwds):
    pass

def plot_2d_diagram(fn, show_function=True, show_projections=True, known_minimal=False, f=None, colorful=False):
    """
    Return a plot of the 2d complex (Delta P) of `fn` with shaded
    additive faces, i.e., faces where delta pi is 0.
    
    If `known_minimal` is False (the default), highlight
    non-subadditive or non-symmetric vertices of the 2d complex.

    If `show_function` is True (the default), plot the function at the left and top borders 
    of the diagram via `plot_function_at_borders`. 

    If `show_projections` is True (the default), plot the projections p1(F), p2(F), p3(F) of 
    all full-dimensional additive faces via `plot_projections_at_borders`.

    To show only a part of the diagram, use::

        sage: plot_2d_diagram(h).show(xmin=0.25, xmax=0.35, ymin=0.25, ymax=0.35)  # not tested

    EXAMPLES::

        sage: h = FastPiecewise([[closed_interval(0,1/4), FastLinearFunction(4, 0)],
        ...                      [open_interval(1/4, 1), FastLinearFunction(4/3, -1/3)],
        ...                      [singleton_interval(1), FastLinearFunction(0,0)]])
        sage: plot_2d_diagram(h)
        Graphics object ...
        sage: h = FastPiecewise([[closed_interval(0,1/4), FastLinearFunction(4, 0)],
        ...                      [open_interval(1/4,1/2), FastLinearFunction(3, -3/4)],
        ...                      [closed_interval(1/2, 3/4), FastLinearFunction(-2, 7/4)],
        ...                      [open_interval(3/4,1), FastLinearFunction(3, -2)],
        ...                      [singleton_interval(1), FastLinearFunction(0,0)]])
        sage: plot_2d_diagram(h)
        Graphics object ...
    """
    if f is None:
        f = find_f(fn, no_error_if_not_minimal_anyway=True)
    faces = generate_maximal_additive_faces(fn)
    p = plot_2d_complex(fn)
    kwds = { 'legend_label': "Additive face" }
    plot_kwds_hook(kwds)
    if colorful:
        covered_components = generate_covered_components(fn)
        colors = rainbow(len(covered_components))
        interval_color = [(interval[0], i) \
                          for (i, component) in enumerate(covered_components) \
                          for interval in component]
        interval_color.sort()
    else:
        covered_components = None
    for face in faces:
        if not covered_components is None and face.is_2D():
            I = face.minimal_triple[0]
            x = (I[0] + I[1]) / 2
            j = bisect_left(interval_color, (x, len(covered_components) + 1)) # should be bisect
            i = interval_color[j-1][1]
            p += face.plot(fill_color = colors[i], **kwds)
        else:
            p += face.plot(**kwds)
        delete_one_time_plot_kwds(kwds)

    ### For non-subadditive functions, show the points where delta_pi is negative.
    if not known_minimal:
        nonsubadditive_vertices = generate_nonsubadditive_vertices(fn, reduced=False)
        kwds = { 'legend_label' : "Subadditivity violated" }
        plot_kwds_hook(kwds)
        if fn.is_continuous():
            nonsubadditive_vertices = {(x,y) for (x, y, z, xeps, yeps, zeps) in nonsubadditive_vertices}
            p += point(list(nonsubadditive_vertices),
                       color = "red", size = 50, zorder=-1, **kwds)
            p += point([ (y,x) for (x,y) in nonsubadditive_vertices ], color = "red", size = 50, zorder=-1)
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
        if f is not None:
            nonsymmetric_vertices = generate_nonsymmetric_vertices(fn, f)
            kwds = { 'legend_label' : "Symmetry violated" }
            plot_kwds_hook(kwds)
            if fn.is_continuous():
                nonsymmetric_vertices = {(x,y) for (x, y, xeps, yeps) in nonsymmetric_vertices}
                p += point(list(nonsymmetric_vertices),
                           color = "mediumvioletred", size = 50, zorder=5, **kwds)
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
    if show_projections:
        p += plot_projections_at_borders(fn)
    if show_function:
        p += plot_function_at_borders(fn, covered_components = covered_components)
    return p

def plot_covered_components_at_borders(fn, covered_components=None, **kwds):
    """
    Colorful decoration.
    Plot function on covered intervals with different colors according to slope values,
    on the upper and the left border of the 2d diagrams.
    """
    p = Graphics()
    if not covered_components is None:
        bkpt = fn.end_points()
        limits = fn.limits_at_end_points()
        colors = rainbow(len(covered_components))
        delete_one_time_plot_kwds(kwds)
        for i, component in enumerate(covered_components):
            for interval in component:
                x1 = interval[0]
                x2 = interval[1]
                y1 = fn.limit(x1, 1)
                y2 = fn.limit(x2, -1)
                p += line([(x1, (3/10)*y1 + 1), (x2, (3/10)*y2 + 1)], color=colors[i], zorder=-2, **kwds)
                p += line([(-(3/10)*y1, x1), (-(3/10)*y2, x2)], color=colors[i], zorder=-2, **kwds)
    return p

def plot_function_at_borders(fn, color='blue', legend_label="Function pi", covered_components=None, **kwds):
    """
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
            p += line([(x1, (3/10)*y1 + 1), (x2, (3/10)*y2 + 1)], color=color, zorder=-2, **kwds)
            delete_one_time_plot_kwds(kwds)
            p += line([(-3/10*y1, x1), (-(3/10)*y2, x2)], color=color, zorder=-2, **kwds)
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
    p += plot_covered_components_at_borders(fn, covered_components, **kwds)
    # add legend_label
    kwds = { 'legend_label': legend_label }
    plot_kwds_hook(kwds)
    if fn.is_discrete():
        p += point([(0,0)], color=color, pointsize=23, zorder=-10, **kwds)
        p += point([(0,0)], color='white', pointsize=23, zorder=-9)
    else:
        p += line([(0,0), (0,1)], color=color, zorder=-10, **kwds)
        p += line([(0,0), (0,1)], color='white', zorder=-9)
    return p

proj_plot_width = 2/100
#proj_plot_colors = ['yellow', 'cyan', 'magenta']            # very clear but ugly
#proj_plot_colors = ['darkseagreen', 'darkseagreen', 'slategray']
proj_plot_colors = ['grey', 'grey', 'grey']
proj_plot_alpha = 0.35
#proj_plot_alpha = 1

def plot_projections_at_borders(fn):
    """
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
        kwds['legend_label'] = "projections p1(F), p2(F), p3(F)"
    elif proj_plot_colors[0] == proj_plot_colors[1]:
        IJK_kwds = [ kwds, kwds, copy(kwds) ]
        kwds['legend_label'] = "projections p1(F), p2(F)"
        IJK_kwds[2]['legend_label'] = "projections p3(F)"
    else:
        IJK_kwds = [ copy(kwds) for i in range(3) ]
        for i in range(3):
            IJK_kwds[i]['legend_label'] = "projections p_%s(F)" % (i+1)
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

def plot_projections_of_one_face(face, IJK_kwds):
    g = Graphics()
    I, J, K = face.minimal_triple
    if face.is_2D():
        # plot I at top border
        g += polygon([(I[0], 1), (I[1], 1), (I[1], 1 + proj_plot_width), (I[0], 1 + proj_plot_width)], **IJK_kwds[0])
        delete_one_time_plot_kwds(IJK_kwds[0])
        # plot J at left border
        g += polygon([(0, J[0]), (0, J[1]), (-proj_plot_width, J[1]), (-proj_plot_width, J[0])], **IJK_kwds[1])
        delete_one_time_plot_kwds(IJK_kwds[1])
        # plot K at right/bottom borders
        if coho_interval_contained_in_coho_interval(K, [0,1]):
            g += polygon([(K[0], 0), (K[1], 0), (K[1] + proj_plot_width, -proj_plot_width), (K[0] + proj_plot_width, -proj_plot_width)], **IJK_kwds[2])
        elif coho_interval_contained_in_coho_interval(K, [1,2]):
            g += polygon([(1, K[0]-1), (1, K[1]-1), (1 + proj_plot_width, K[1] - 1 - proj_plot_width), (1 + proj_plot_width, K[0] - 1 - proj_plot_width)], **IJK_kwds[2])
        else:
            raise ValueError, "Bad face: %s" % face
        delete_one_time_plot_kwds(IJK_kwds[2])
    return g
    
def interval_mod_1(interval):
    """
    Represent the given proper interval modulo 1
    as a subinterval of [0,1].

    EXAMPLES::

        sage: interval_mod_1([1,6/5])
        [0, 1/5]
        sage: interval_mod_1([1,2])
        [0, 1]
        sage: interval_mod_1([-3/10,-1/10])
        [7/10, 9/10]
        sage: interval_mod_1([-1/5,0])
        [4/5, 1]        
    """
    interval = copy(interval)
    if len(interval) == 0:
        return interval
    elif len(interval) == 1:
        while interval[0] >= 1:
            interval[0] -= 1
        while interval[0] < 0:
            interval[0] += 1
        return interval
    elif len(interval) == 2:
        assert interval[0] < interval[1]
        while interval[0] >= 1:
            interval[0] = interval[0] - 1
            interval[1] = interval[1] - 1
        while interval[1] <= 0:
            interval[0] = interval[0] + 1
            interval[1] = interval[1] + 1
        assert not(interval[0] < 1 and interval[1] > 1) 
        return interval
    else:
        raise ValueError, "Not an interval: %s" % interval

def generate_covered_components(function):
    fdms, covered_components = generate_directed_move_composition_completion(function)
    return covered_components

# alias
generate_covered_intervals = generate_covered_components

def generate_uncovered_intervals(function):
    covered_components = generate_covered_components(function)
    return uncovered_intervals_from_covered_components(covered_components)

def uncovered_intervals_from_covered_components(covered_components):
    uncovered_intervals = union_of_coho_intervals_minus_union_of_coho_intervals([[open_interval(0,1)]], covered_components, remove_closure=True)
    return uncovered_intervals

# alias
uncovered_intervals_from_covered_intervals = uncovered_intervals_from_covered_components
    
def ticks_keywords(function, y_ticks_for_breakpoints=False):
    """
    Compute `plot` keywords for displaying the ticks.
    """
    xticks = function.end_points()
    f = find_f(function, no_error_if_not_minimal_anyway=True)
    if f is not None and not f in xticks:
        xticks.append(f)
    xtick_formatter = [ "$%s$" % latex(x) for x in xticks ]
    #xtick_formatter = 'latex'  # would not show rationals as fractions
    ytick_formatter = None
    if y_ticks_for_breakpoints:
        yticks = xticks
        ytick_formatter = xtick_formatter
    else:
        #yticks = 1/5
        yticks = uniq([1] + [ y for limits in function.limits_at_end_points() for y in limits if y is not None ])
        ytick_formatter = [ "$%s$" % latex(y) for y in yticks ]
    ## FIXME: Can we influence ticks placement as well so that labels don't overlap?
    ## or maybe rotate labels 90 degrees?
    return {'ticks': [xticks, yticks], \

            'gridlines': True, \
            'tick_formatter': [xtick_formatter, ytick_formatter]}

def delete_one_time_plot_kwds(kwds):
    if 'legend_label' in kwds:
        del kwds['legend_label']
    if 'ticks' in kwds:
        del kwds['ticks']
    if 'tick_formatter' in kwds:
        del kwds['tick_formatter']

def plot_covered_intervals(function, covered_components=None, uncovered_color='black', labels=None, **plot_kwds):
    """
    Return a plot of the covered and uncovered intervals of `function`.
    """
    if covered_components is None:
        covered_components = generate_covered_components(function)
        uncovered_intervals = generate_uncovered_intervals(function)
    else:
        uncovered_intervals = uncovered_intervals_from_covered_components(covered_components)
    # Plot the function with different colors.
    # Each component has a unique color.
    # The uncovered intervals is by default plotted in black.
    colors = rainbow(len(covered_components))
    graph = Graphics()
    kwds = copy(plot_kwds)
    kwds.update(ticks_keywords(function))
    if uncovered_intervals:
        kwds.update({'legend_label': "not covered"})
        plot_kwds_hook(kwds)
        graph += plot(function, color = uncovered_color, **kwds)
        delete_one_time_plot_kwds(kwds)
    elif not function.is_continuous(): # to plot the discontinuity markers
        graph += plot(function, color = uncovered_color, **kwds)
        delete_one_time_plot_kwds(kwds)
    for i, component in enumerate(covered_components):
        if labels is None:
            label = "covered component %s" % (i+1)
        else:
            label = labels[i]
        kwds.update({'legend_label': label})
        plot_kwds_hook(kwds)
        for interval in component:
            # We do not plot anything if float(interval[0])==float(interval[1]) because
            # otherwise plot complains that
            # "start point and endpoint must be different"
            if float(interval[0])<float(interval[1]):
                graph += plot(function.which_function((interval[0] + interval[1])/2), (interval[0], interval[1]), color=colors[i], zorder=-1, **kwds)
                # zorder=-1 puts them below the discontinuity markers,
                # above the black function.
                delete_one_time_plot_kwds(kwds)
    return graph

def plot_directly_covered_intervals(function, uncovered_color='black', labels=None, **plot_kwds):
    components = generate_directly_covered_components(function)
    return plot_covered_intervals(function, components, uncovered_color=uncovered_color, labels=labels, **plot_kwds)

def number_of_components(fn):
    """
    Return the number of connected components of `fn`.

    This is an upper bound on `number_of_slopes`.

    EXAMPLES::

        sage: logging.disable(logging.WARN)
        sage: number_of_components(gmic())
        2
        sage: number_of_components(gomory_fractional())
        1
        sage: number_of_slopes(not_extreme_1())
        3
        sage: number_of_components(not_extreme_1())
        4
    """
    covered_components = generate_covered_components(fn)
    return len(covered_components)

def slopes_intervals_dict(fn):
    """
    Return a dictionary mapping slope values to a list of intervals of that slope.

    EXAMPLES::

        sage: logging.disable(logging.INFO)
        sage: slopes_intervals_dict(gmic(1/2))[2]
        [(0, 1/2)]
    """
    slopes_dict = dict()
    for i, f in fn.list():
        if interval_length(i) > 0:
            try: # Make sure we don't fail if passed a non-FastLinearFunction
                if f._slope not in slopes_dict:
                    slopes_dict[f._slope] = []
                slopes_dict[f._slope].append((i[0], i[1]))
            except AttributeError:
                pass
    return slopes_dict

def number_of_slopes(fn):
    """
    Return the number of different slopes of `fn`.

    If `fn` is discrete, this is defined as the number of different slopes
    of its piecewise linear continuous interpolation.

    EXAMPLES::

        sage: logging.disable(logging.WARN)
        sage: number_of_slopes(gmic())
        2
        sage: number_of_slopes(gomory_fractional())
        1
        sage: number_of_slopes(automorphism(restrict_to_finite_group(gmic(10/11)), 3))
        2
    """
    if fn.is_discrete():
        fn = interpolate_to_infinite_group(fn)
    return len(slopes_intervals_dict(fn))

def plot_with_colored_slopes(fn):
    """
    Return a plot of `fn`, with pieces of different slopes in different colors.
    """
    slopes_dict = slopes_intervals_dict(fn)
    return plot_covered_intervals(fn, slopes_dict.values(), labels=[ "Slope %s" % s for s in slopes_dict.keys() ])

### Minimality check.

def subadditivity_test(fn):
    """
    Check if `fn` is subadditive.
    """
    result = True
    for (x, y, z, xeps, yeps, zeps) in generate_nonsubadditive_vertices(fn, reduced=True):
        logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(xeps), y, print_sign(yeps), z, print_sign(zeps)))
        result = False
    if result:
        logging.info("pi is subadditive.")
    else:
        logging.info("Thus pi is not subadditive.")
    return result

def symmetric_test(fn, f):
    """
    Check if `fn` is symmetric.
    """
    result = True
    if fn(f) != 1:
        logging.info('pi(f) is not equal to 1.')
        result = False
    result = True
    for (x, y, xeps, yeps) in generate_nonsymmetric_vertices(fn, f):
        logging.info("pi(%s%s) + pi(%s%s) is not equal to 1" % (x, print_sign(xeps), y, print_sign(yeps)))
        result = False
    if result:
        logging.info('pi is symmetric.')
    else:
        logging.info('Thus pi is not symmetric.')
    return result

def find_f(fn, no_error_if_not_minimal_anyway=False):
    """
    Find the value of `f' for the given function `fn'.
    """
    if hasattr(fn, '_f'):
        return fn._f
    f = None
    for x in fn.end_points():
        if fn(x) > 1 or fn(x) < 0: 
            if no_error_if_not_minimal_anyway:
                logging.info('pi is not minimal because it does not stay in the range of [0, 1].')
                return None
            raise ValueError, "The given function does not stay in the range of [0, 1], so cannot determine f.  Provide parameter f to minimality_test or extremality_test."
    for x in fn.end_points():
        if fn(x) == 1:
            if not f is None:
                logging.warn("The given function has more than one breakpoint where the function takes the value 1; using f = %s.  Provide parameter f to minimality_test or extremality_test if you want a different f." % f)
                return f
            else:
                f = x
    if not f is None:
        fn._f = f
        return f
    if no_error_if_not_minimal_anyway:
        logging.info('pi is not minimal because it has no breakpoint where the function takes value 1.')
        return None
    raise ValueError, "The given function has no breakpoint where the function takes value 1, so cannot determine f.  Provide parameter f to minimality_test or extremality_test."

def minimality_test(fn, show_plots=False, f=None):
    """
    Check if `fn` is minimal with respect to the group relaxation with the given `f`. 

    If `f` is not provided, use the one found by `find_f`.

    If `show_plots` is True (default: False), show an illustrating diagram.

    This function verifies that function values stay between 0 and 1 and
    calls `subadditivity_test` and `symmetric_test`.

    EXAMPLES::

        sage: logging.disable(logging.INFO)
        sage: minimality_test(piecewise_function_from_breakpoints_and_values([0,1/5,4/5,1],[0,1/2,1,0]))
        False
        sage: minimality_test(piecewise_function_from_breakpoints_and_values([0,1/2,1], [0,2,0]))
        False
    """
    for x in fn.values_at_end_points():
        if (x < 0) or (x > 1):
            logging.info('pi is not minimal because it does not stay in the range of [0, 1].')
            return False
    if f is None:
        f = find_f(fn, no_error_if_not_minimal_anyway=True)
        if f is None:
            return False
    if fn(0) != 0:
        logging.info('pi is NOT minimal because pi(0) is not equal to 0.')
        return False
    logging.info('pi(0) = 0')
    bkpt = fn.end_points()
    if not fn.is_continuous():
        limits = fn.limits_at_end_points()
        for x in limits:
            if not ((x[-1] is None or 0 <= x[-1] <=1) and (x[1] is None or 0 <= x[1] <=1)):
                logging.info('pi is not minimal because it does not stay in the range of [0, 1].')
                return False
    if subadditivity_test(fn) and symmetric_test(fn, f):
        logging.info('Thus pi is minimal.')
        is_minimal = True
    else:
        logging.info('Thus pi is NOT minimal.')
        is_minimal = False
    if show_plots:
        logging.info("Plotting 2d diagram...")
        show_plot(plot_2d_diagram(fn, known_minimal=is_minimal, f=f),
                  show_plots, tag='2d_diagram', object=fn)
        logging.info("Plotting 2d diagram... done")
    return is_minimal

from sage.functions.piecewise import PiecewisePolynomial
from bisect import bisect_left

class FastPiecewise (PiecewisePolynomial):
    """
    Returns a piecewise function from a list of (interval, function)
    pairs.

    Uses binary search to allow for faster function evaluations
    than the standard class PiecewisePolynomial.
    """
    def __init__(self, list_of_pairs, var=None, periodic_extension=True, merge=True):
        """
        EXAMPLES::

            sage: h = FastPiecewise([[(3/10, 15/40), FastLinearFunction(1, 0)], [(13/40, 14/40), FastLinearFunction(1, 0)]], merge=True)
            sage: len(h.intervals())
            1
            sage: h.intervals()[0][0], h.intervals()[0][1]
            (3/10, 3/8)
            sage: h = FastPiecewise([[(3/10, 15/40), FastLinearFunction(1, 0)], [(13/40, 14/40), FastLinearFunction(1, 0)], [(17,18), FastLinearFunction(77,78)]], merge=True)
            sage: len(h.intervals())
            2
            sage: h.intervals()[0][0], h.intervals()[0][1]
            (3/10, 3/8)
        """
        # Sort intervals according to their left endpoints; In case of equality, place single point before interval. 
        list_of_pairs = sorted(list_of_pairs, key = lambda (i, f): coho_interval_left_endpoint_with_epsilon(i))
        if merge:
            merged_list_of_pairs = []
            intervals_to_scan = []
            singleton = None
            common_f = None
            for (i, f) in list_of_pairs:
                if len(i) == 1:
                    i = singleton_interval(i[0])            # upgrade to coho interval
                if common_f == f:
                    intervals_to_scan.append(i)
                    singleton = None
                elif common_f is not None and singleton is not None and common_f(singleton) == f(singleton):
                    intervals_to_scan.append(i)
                    singleton = None
                    common_f = f
                elif i[0] == i[1] and common_f is not None and common_f(i[0]) == f(i[0]):
                    intervals_to_scan.append(i)
                else:
                    merged_intervals = union_of_coho_intervals_minus_union_of_coho_intervals([[interval] for interval in intervals_to_scan], [],
                                                                                             old_fashioned_closed_intervals=True)
                    for merged_interval in merged_intervals:
                        merged_list_of_pairs.append((merged_interval, common_f))
                    intervals_to_scan = [i]
                    if i[0] == i[1]:
                        singleton = i[0]
                    else:
                        singleton = None
                    common_f = f
            merged_intervals = union_of_coho_intervals_minus_union_of_coho_intervals([[interval] for interval in intervals_to_scan], [],
                                                                                     old_fashioned_closed_intervals=True)
            for merged_interval in merged_intervals:
                merged_list_of_pairs.append((merged_interval, common_f))
            list_of_pairs = merged_list_of_pairs
            
        PiecewisePolynomial.__init__(self, list_of_pairs, var)

        intervals = self._intervals
        functions = self._functions
        # end_points are distinct.
        end_points = []
        # ith_at_end_points records in which interval the end_point first appears as a left_end or right_end.
        ith_at_end_points = []
        # record the value at each end_point, value=None if end_point is not in the domain.
        values_at_end_points = []
        # record function values at [x, x+, x-] for each endpoint x.
        limits_at_end_points = []
        left_limit = None
        for i in range(len(intervals)):
            left_value = None
            if len(intervals[i]) <= 2 or intervals[i].left_closed:
                left_value = functions[i](intervals[i][0])
            if intervals[i][0] != intervals[i][1]:
                right_limit = functions[i](intervals[i][0])
            else:
                right_limit = None
            if (end_points == []) or (end_points[-1] != intervals[i][0]):
                end_points.append(intervals[i][0])
                ith_at_end_points.append(i)
                values_at_end_points.append(left_value)
                if limits_at_end_points != []:
                    limits_at_end_points[-1][1]= None
                limits_at_end_points.append([left_value, right_limit, None])
            else:
                if left_value is not None:
                    values_at_end_points[-1] = left_value
                    limits_at_end_points[-1][0] = left_value
                limits_at_end_points[-1][1] = right_limit
            right_value = None
            if len(intervals[i]) <= 2 or intervals[i].right_closed:
                right_value = functions[i](intervals[i][1])
            if intervals[i][0] != intervals[i][1]:
                left_limit = functions[i](intervals[i][1])
                end_points.append(intervals[i][1])
                ith_at_end_points.append(i)
                values_at_end_points.append(right_value)
                limits_at_end_points.append([right_value, None, left_limit])
            elif right_value is not None:
                values_at_end_points[-1] = right_value        
        if periodic_extension and limits_at_end_points != []:
            #if values_at_end_points[0] != values_at_end_points[-1]:
            #    logging.warn("Function is actually not periodically extendable.")
            #    periodic_extension = False
            #else:
                limits_at_end_points[0][-1] = limits_at_end_points[-1][-1]
                limits_at_end_points[-1][1] = limits_at_end_points[0][1]
        self._end_points = end_points
        self._ith_at_end_points = ith_at_end_points
        self._values_at_end_points = values_at_end_points
        self._limits_at_end_points = limits_at_end_points
        self._periodic_extension = periodic_extension
        self._call_cache = dict()

        is_continuous = True
        if len(end_points) == 1 and end_points[0] is None:
            is_continuous = False
        elif len(end_points)>= 2:
            [m0, r0, l0] = limits_at_end_points[0]
            [m1, r1, l1] = limits_at_end_points[-1]
            if m0 is None or r0 is None or  m0 != r0 or l1 is None or m1 is None or l1 != m1:
                is_continuous = False
            else:
                for i in range(1, len(end_points)-1):
                    [m, r, l] = limits_at_end_points[i]
                    if l is None or m is None or r is None or not(l == m == r):
                        is_continuous = False
                        break
        self._is_continuous = is_continuous
        self._is_two_sided_discontinuous = not ( is_continuous or \
                                                 limits_at_end_points[0][0] == limits_at_end_points[0][1] or
                                                 limits_at_end_points[-1][-1] == limits_at_end_points[-1][0] )

    # The following makes this class hashable and thus enables caching
    # of the above functions; but we must promise not to modify the
    # contents of the instance.
    def __hash__(self):
        return id(self)

    def is_continuous(self):
        """
        return if function is continuous
        """
        return self._is_continuous

    def is_two_sided_discontinuous(self):
        """
        return if function is discontinuous at 0+ and at 1-.
        """
        return self._is_two_sided_discontinuous
        
    def is_discrete(self):
        """
        Return if the function is discrete, i.e., all pieces are singletons
        """
        return all(interval_length(interval) == 0 for interval in self.intervals())

    def end_points(self):
        """
        Returns a list of all interval endpoints for this function.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 2
            sage: f3(x) = 1-x
            sage: f4(x) = x^2-5
            sage: f = FastPiecewise([[open_interval(0,1),f1],[singleton_interval(1),f2],[open_interval(1,2),f3],[(2,3),f4]])
            sage: f.end_points()
            [0, 1, 2, 3]
            sage: f = FastPiecewise([[open_interval(0,1),f1],[open_interval(2,3),f3]])
            sage: f.end_points()
            [0, 1, 2, 3]
        """
        return self._end_points

    def values_at_end_points(self):
        """
        Returns a list of function values at all endpoints for this function.

        EXAMPLES::

            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = 4
            sage: f5(x) = sin(2*x)
            sage: f6(x) = x-3
            sage: f7(x) = 7
            sage: f = FastPiecewise([[right_open_interval(0,1),f1], \
            ...                      [right_open_interval(1,2),f2],\
            ...                      [open_interval(2,3),f3],\
            ...                      [singleton_interval(3),f4],\
            ...                      [left_open_interval(3,6),f5],\
            ...                      [open_interval(6,7),f6],\
            ...                      [(9,10),f7]])
            sage: f.values_at_end_points()
            [1, 0, None, 4, sin(12), None, 7, 7]
        """
        return self._values_at_end_points

    def limits_at_end_points(self):
        """
        Returns a list of 3-tuples [function value, right_limit, left_limit] at all endpoints for this function.

        EXAMPLES::

            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = 4
            sage: f5(x) = sin(2*x)
            sage: f6(x) = x-3
            sage: f7(x) = 7
            sage: f = FastPiecewise([[right_open_interval(0,1),f1], \
            ...                      [right_open_interval(1,2),f2],\
            ...                      [open_interval(2,3),f3],\
            ...                      [singleton_interval(3),f4],\
            ...                      [left_open_interval(3,6),f5],\
            ...                      [open_interval(6,7),f6],\
            ...                      [(9,10),f7]], periodic_extension= False)
            sage: f.limits_at_end_points()
            [[1, 1, None], [0, 0, 1], [None, e^2, -1], [4, sin(6), e^3], [sin(12), 3, sin(12)], [None, None, 4], [7, 7, None], [7, None, 7]]
        """
        return self._limits_at_end_points

    def which_function(self, x0):
        """
        Returns the function piece used to evaluate self at x0.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = FastPiecewise([[(0,1),f1],
            ...                      [(1,2),f2],
            ...                      [(2,3),f3],
            ...                      [(3,10),f4]])
            sage: f.which_function(0.5) is f1
            True
            sage: f.which_function(1) in [f1, f2]
            True
            sage: f.which_function(5/2) is f3
            True
            sage: f.which_function(3) in [f3, f4]
            True
            sage: f.which_function(-1)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point -1, outside of domain.
            sage: f.which_function(11)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 11, outside of domain.
            sage: f = FastPiecewise([[right_open_interval(0,1),f1],
            ...                      [right_open_interval(1,2),f2],
            ...                      [right_open_interval(2,3),f3],
            ...                      [closed_interval(3,10),f4]])
            sage: f.which_function(0.5) is f1
            True
            sage: f.which_function(1) is f2
            True
            sage: f.which_function(5/2) is f3
            True
            sage: f.which_function(3) is f4
            True
            sage: f = FastPiecewise([[open_interval(0,1),f1],
            ...                      [right_open_interval(2,3),f3]])
            sage: f.which_function(0)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 0, outside of domain.
            sage: f.which_function(0.5) is f1
            True
            sage: f.which_function(1)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 1, outside of domain.
            sage: f.which_function(3/2)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 3/2, outside of domain.
            sage: f.which_function(2) is f3
            True
            sage: f.which_function(5/2) is f3
            True
            sage: f.which_function(3)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 3, outside of domain.
        """
        endpts = self.end_points()
        ith = self._ith_at_end_points
        i = bisect_left(endpts, x0)
        if i >= len(endpts):
            raise ValueError,"Value not defined at point %s, outside of domain." % x0
        if x0 == endpts[i]:
            if self._values_at_end_points[i] is not None:
                if self.functions()[ith[i]](x0) == self._values_at_end_points[i]:
                    return self.functions()[ith[i]]
                else:
                    return self.functions()[ith[i]+1]
            else:
                raise ValueError,"Value not defined at point %s, outside of domain." % x0
        if i == 0:
            raise ValueError,"Value not defined at point %s, outside of domain." % x0
        if is_pt_in_interval(self._intervals[ith[i]],x0):
            return self.functions()[ith[i]]
        raise ValueError,"Value not defined at point %s, outside of domain." % x0

    def __call__(self,x0):
        """
        Evaluates self at x0. 
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = FastPiecewise([[(0,1),f1],
            ...                      [(1,2),f2],
            ...                      [(2,3),f3],
            ...                      [(3,10),f4]])
            sage: f(0.5)
            1
            sage: f(1)
            0
            sage: f(5/2)
            e^(5/2)
            sage: f(3)
            sin(6)
            sage: f(-1)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point -1, outside of domain.
            sage: f(11)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 11, outside of domain.
            sage: f = FastPiecewise([[right_open_interval(0,1),f1],
            ...                      [right_open_interval(1,2),f2],
            ...                      [right_open_interval(2,3),f3],
            ...                      [closed_interval(3,10),f4]])
            sage: f(0.5)
            1
            sage: f(1)
            0
            sage: f(5/2)
            e^(5/2)
            sage: f(3)
            sin(6)
            sage: f = FastPiecewise([[open_interval(0,1),f1],
            ...                      [right_open_interval(2,3),f3]])
            sage: f(0)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 0, outside of domain.
            sage: f(0.5)
            1
            sage: f(1)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 1, outside of domain.
            sage: f(3/2)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 3/2, outside of domain.
            sage: f(2)
            e^2
            sage: f(5/2)
            e^(5/2)
            sage: f(3)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 3, outside of domain.
        """
        # fast path 
        result = self._call_cache.get(x0)
        if result is not None:
            return result
        # Remember that intervals are sorted according to their left endpoints; singleton has priority.
        endpts = self.end_points()
        ith = self._ith_at_end_points
        i = bisect_left(endpts, x0)
        if i >= len(endpts):
            raise ValueError,"Value not defined at point %s, outside of domain." % x0
        if x0 == endpts[i]:
            if self._values_at_end_points[i] is not None:
                result = self._values_at_end_points[i]
                self._call_cache[x0] = result
                return result
            else:
                raise ValueError,"Value not defined at point %s, outside of domain." % x0
        if i == 0:
            raise ValueError,"Value not defined at point %s, outside of domain." % x0
        if is_pt_in_interval(self._intervals[ith[i]],x0):
            result = self.functions()[ith[i]](x0)
            self._call_cache[x0] = result
            return result
        raise ValueError,"Value not defined at point %s, outside of domain." % x0

    def limits(self, x0):
        """
        return [function value at x0, function value at x0+, function value at x0-].

        EXAMPLES::

            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = 4
            sage: f5(x) = sin(2*x)
            sage: f6(x) = x-3
            sage: f7(x) = 7
            sage: f = FastPiecewise([[right_open_interval(0,1),f1], \
            ...                      [right_open_interval(1,2),f2],\
            ...                      [open_interval(2,3),f3],\
            ...                      [singleton_interval(3),f4],\
            ...                      [left_open_interval(3,6),f5],\
            ...                      [open_interval(6,7),f6],\
            ...                      [(9,10),f7]], periodic_extension=False)
            sage: f.limits(1/2)
            [1, 1, 1]
            sage: f.limits(1)
            [0, 0, 1]
            sage: f.limits(2)
            [None, e^2, -1]
            sage: f.limits(3)
            [4, sin(6), e^3]
            sage: f.limits(6)
            [sin(12), 3, sin(12)]
            sage: f.limits(7)
            [None, None, 4]
            sage: f.limits(8)
            [None, None, None]
            sage: f.limits(9)
            [7, 7, None]
        """
        endpts = self.end_points()
        ith = self._ith_at_end_points
        i = bisect_left(endpts, x0)
        if i >= len(endpts):
            return [None, None, None]
        if x0 == endpts[i]:
            return self.limits_at_end_points()[i]
        if i == 0:
            return [None, None, None]
        if is_pt_in_interval(self._intervals[ith[i]],x0):
            result = self.functions()[ith[i]](x0)
            return [result, result, result]
        return [None, None, None]

    def limit(self, x0, epsilon):
        """
        return limit (from right if epsilon > 0, from left if epsilon < 0) value at x0;
        if epsilon == 0, return value at x0.

        EXAMPLES::

            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = 4
            sage: f5(x) = sin(2*x)
            sage: f6(x) = x-3
            sage: f7(x) = 7
            sage: f = FastPiecewise([[right_open_interval(0,1),f1], \
            ...                      [right_open_interval(1,2),f2],\
            ...                      [open_interval(2,3),f3],\
            ...                      [singleton_interval(3),f4],\
            ...                      [left_open_interval(3,6),f5],\
            ...                      [open_interval(6,7),f6],\
            ...                      [(9,10),f7]], periodic_extension=False)
            sage: f.limit(1,0)
            0
            sage: f.limit(1,1)
            0
            sage: f.limit(2,-1)
            -1
            sage: f.limit(2,0)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 2, outside of domain.
            sage: f.limit(7,1)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 7+, outside of domain.
            sage: f.limit(8,-1)
            Traceback (most recent call last):
            ...
            ValueError: Value not defined at point 8-, outside of domain.
        """
        result =self.limits(x0)[epsilon]
        if result is None:
            raise ValueError,"Value not defined at point %s%s, outside of domain." % (x0, print_sign(epsilon))
        return result

    def which_function_on_interval(self, interval):
        x = (interval[0] + interval[1]) / 2
        # FIXME: This should check that the given `interval` is contained in the defining interval!
        # This could be implemented by refactoring which_function using new function which_function_index.
        return self.which_function(x)

    def __add__(self,other):
        """
        Add `self` and another piecewise function.

        In contrast to PiecewisePolynomial.__add__, this does not do zero extension of domains.
        Rather, the result is only defined on the intersection of the domains.

        EXAMPLES::

            sage: f = FastPiecewise([[singleton_interval(1), FastLinearFunction(0,17)]])
            sage: g = FastPiecewise([[[0,2], FastLinearFunction(0,2)]])
            sage: (f+g).list()
            [[<Int{1}>, <FastLinearFunction 19>]]
            sage: h = FastPiecewise([[open_interval(1,3), FastLinearFunction(0,3)]])
            sage: (g+h).list()
            [[<Int(1, 2]>, <FastLinearFunction 5>]]
            sage: j = FastPiecewise([[open_interval(0,1), FastLinearFunction(0,1)], [[1, 3], FastLinearFunction(0, 5)]])
            sage: (g+j).list()
            [[<Int(0, 1)>, <FastLinearFunction 3>], [(1, 2), <FastLinearFunction 7>]]
        """
        intervals = intersection_of_coho_intervals([self.intervals(), other.intervals()])
        return FastPiecewise([ (interval, self.which_function_on_interval(interval) + other.which_function_on_interval(interval))
                               for interval in intervals ], merge=True)

    def __neg__(self):
        return FastPiecewise([[interval, -f] for interval,f in self.list()], merge=True)
        
    def __mul__(self,other):
        """
        Multiply `self` by a scalar or another piecewise function.

        In contrast to PiecewisePolynomial.__mul__, this does not do zero extension of domains.
        Rather, the result is only defined on the intersection of the domains.
        """
        if not isinstance(other, FastPiecewise):
            # assume scalar multiplication
            return FastPiecewise([[interval, other*f] for interval,f in self.list()])
        else:
            intervals = intersection_of_coho_intervals([self.intervals(), other.intervals()])
            return FastPiecewise([ (interval, self.which_function_on_interval(interval) * other.which_function_on_interval(interval))
                                   for interval in intervals ], merge=True)

    __rmul__ = __mul__

    def __div__(self, other):
        return self * (1 / other)

    def __sub__(self, other):
        return self + (-other)

    ## Following just fixes a bug in the plot method in piecewise.py
    ## (see doctests below).  Also adds plotting of single points
    ## and discontinuity markers.
    def plot(self, *args, **kwds):
        """
        Returns the plot of self.
        
        Keyword arguments are passed onto the plot command for each piece
        of the function. E.g., the plot_points keyword affects each
        segment of the plot.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = FastPiecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            sage: P = f.plot(rgbcolor=(0.7,0.1,0), plot_points=40)
            sage: P
            Graphics object...
        
        Remember: to view this, type show(P) or P.save("path/myplot.png")
        and then open it in a graphics viewer such as GIMP.

        TESTS:

        We should not add each piece to the legend individually, since
        this creates duplicates (:trac:`12651`). This tests that only
        one of the graphics objects in the plot has a non-``None``
        ``legend_label``::

            sage: f1(x) = sin(x)
            sage: f2(x) = cos(x)
            sage: f = FastPiecewise([[(-1,0), f1],[(0,1), f2]])
            sage: p = f.plot(legend_label='$f(x)$')
            sage: lines = [
            ...     line
            ...     for line in p._objects
            ...     if line.options()['legend_label'] is not None ]
            sage: len(lines)
            1

        The implementation of the plot method in Sage 5.11 piecewise.py
        is incompatible with the use of the xmin and xmax arguments.  Test that
        this has been fixed::

            sage: q = f.plot(xmin=0, xmax=3)
            sage: q = plot(f, xmin=0, xmax=3)
            sage: q = plot(f, 0, 3)
            sage: q = plot(f, 0, 3, color='red')
        
        The implementation should crop according to the given xmin, xmax::

            sage: q = plot(f, 1/2, 3)
            sage: q = plot(f, 1, 2)
            sage: q = plot(f, 2, 3)
        
        Also the following plot syntax should be accepted::

            sage: q = plot(f, [2, 3])

        """
        from sage.plot.all import plot, Graphics

        g = Graphics()
        if not 'rgbcolor' in kwds:
            color = 'blue'
        else:
            color = kwds['rgbcolor']
        if not 'plot_points' in kwds:
            plot_pts = 200
        else:
            plot_pts = kwds['plot_points']
        ### Code duplication with xmin/xmax code in plot.py.
        n = len(args)
        xmin = None
        xmax = None
        if n == 0:
            # if there are no extra args, try to get xmin,xmax from
            # keyword arguments
            xmin = kwds.pop('xmin', None)
            xmax = kwds.pop('xmax', None)
        elif n == 1:
            # if there is one extra arg, then it had better be a tuple
            xmin, xmax = args[0]
            args = []
            ## The case where the tuple is longer than 2 elements is for the 
            ## case of symbolic expressions; it does not apply here.
            ## FIXME: We should probably signal an error.
        elif n == 2:
            # if there are two extra args, they should be xmin and xmax
            xmin = args[0]
            xmax = args[1]
            args = []
        ## The case with three extra args is for the case of symbolic
        ## expressions; it does not apply here.  FIXME: We should
        ## probably signal an error.
        point_kwds = dict()
        if 'alpha' in kwds:
            point_kwds['alpha'] = kwds['alpha']
        if 'legend_label' in kwds and self.is_discrete():
            point_kwds['legend_label'] = kwds['legend_label']
        # Whether to plot discontinuity markers
        discontinuity_markers = kwds.pop('discontinuity_markers', True)
        # record last right endpoint, then compare with next left endpoint to decide whether it needs to be plotted.
        last_end_point = []
        last_closed = True
        for (i, f) in self.list():
            a = i[0]
            b = i[1]
            left_closed = True
            right_closed = True
            if len(i) > 2: # coho interval
                left_closed = i.left_closed
                right_closed = i.right_closed
            # using the above data.
            if (xmin is not None) and (a < xmin):
                a = xmin
                left_closed = True
            if (xmax is not None) and (b > xmax):
                b = xmax
                right_closed = True
            if discontinuity_markers:
                # Handle open/half-open intervals here
                if a < b or (a == b and left_closed and right_closed):
                    if not (last_closed or last_end_point == [a, f(a)] and left_closed):
                        # plot last open right endpoint
                        g += point(last_end_point, color=color, pointsize=23, **point_kwds)
                        delete_one_time_plot_kwds(point_kwds)
                        g += point(last_end_point, rgbcolor='white', pointsize=10, **point_kwds)
                    if last_closed and last_end_point != [] and last_end_point != [a, f(a)] and not left_closed:
                        # plot last closed right endpoint
                        g += point(last_end_point, color=color, pointsize=23, **point_kwds)
                        delete_one_time_plot_kwds(point_kwds)
                    if not (left_closed or last_end_point == [a, f(a)] and last_closed):
                        # plot current open left endpoint
                        g += point([a, f(a)], color=color, pointsize=23, **point_kwds)
                        delete_one_time_plot_kwds(point_kwds)
                        g += point([a, f(a)], rgbcolor='white', pointsize=10, **point_kwds)
                    if left_closed and last_end_point != [] and last_end_point != [a, f(a)] and not last_closed:
                        # plot current closed left endpoint
                        g += point([a, f(a)], color=color, pointsize=23, **point_kwds)
                        delete_one_time_plot_kwds(point_kwds)
                    last_closed = right_closed
                    last_end_point = [b, f(b)]
            if a < b and (float(b) - float(a))/(plot_pts-1) != float(0):
                # We do not plot anything if (float(b) - float(a))/(plot_pts-1) == float(0) because
                # otherwise the plot method in src/plot/misc.py complains that
                # "start point and endpoint must be different"
                g += plot(f, *args, xmin=a, xmax=b, zorder=-1, **kwds)
                # If it's the first piece, pass all arguments. Otherwise,
                # filter out 'legend_label' so that we don't add each
                # piece to the legend separately (trac #12651).
                delete_one_time_plot_kwds(kwds)
                #delete_one_time_plot_kwds(point_kwds)
            elif a == b and left_closed and right_closed:
                g += point([a, f(a)], color=color, pointsize=23, **point_kwds)
                delete_one_time_plot_kwds(point_kwds)
        # plot open rightmost endpoint. minimal functions don't need this.
        if discontinuity_markers and not last_closed:
            g += point(last_end_point, color=color,pointsize=23, **point_kwds)
            delete_one_time_plot_kwds(point_kwds)
            g += point(last_end_point, rgbcolor='white', pointsize=10, **point_kwds)
        return g

    def is_continuous_defined(self, xmin=0, xmax=1):
        """
        return True if self is defined on [xmin,xmax] and is continuous on [xmin,xmax]
        """
        bkpt = self._end_points
        if xmin < bkpt[0] or xmax > bkpt[-1]:
            return False
        if xmin == xmax:
            return (self(xmin) is not None)
        limits = self._limits_at_end_points
        i = 0
        while bkpt[i] < xmin:
            i += 1
        if bkpt[i] == xmin:
            if limits[i][0] is None or limits[i][1] is None or limits[i][0] != limits[i][1]:
                return False 
            i += 1
        while bkpt[i] < xmax:
            if limits[i][-1] is None or limits[i][0] is None or limits[i][1] is None or \
                                        not (limits[i][-1] == limits[i][0] == limits[i][1]):
                return False
            i += 1
        if bkpt[i] == xmax:
            if limits[i][0] is None or limits[i][-1] is None or limits[i][0] != limits[i][-1]:
                return False
        return True

    def __repr__(self):
        rep = "<FastPiecewise with %s parts, " % len(self._functions)
        for interval, function in itertools.izip(self._intervals, self._functions):
            rep += "\n " + repr(interval) + "\t" + repr(function) \
                   + "\t values: " + repr([function(interval[0]), function(interval[1])])
        rep += ">"
        return rep

    def _sage_input_(self, sib, coerced):
        """
        Produce an expression which will reproduce this value when evaluated.
        """
        # FIXME: Add keyword arguments
        # FIXME: "sage_input(..., verify=True)" does not work yet
        # because of module trouble?
        return sib.name('FastPiecewise')(sib(self.list()))

    def sha1(self):
        """
        Return a SHA-1 hash of the function.

        The hash is intended to stay stable even when the code is updated.

        Merged and unmerged versions have the same hash.

        TESTS::

            sage: h1 = piecewise_function_from_breakpoints_and_slopes([0, 1/4, 1/2, 1], [1, 1, -1], merge=False)
            sage: h1.sha1()
            'c562cf38581076609876b1c4fab604756690db7b'
            sage: h2 = piecewise_function_from_breakpoints_and_slopes([0, 1/4, 1/2, 1], [1, 1, -1], merge=True)
            sage: h2.sha1()
            'c562cf38581076609876b1c4fab604756690db7b'

        """
        from hashlib import sha1
        self_merged = self * 1            # in case we were constructed with merge=False!
        data = zip(self_merged.end_points(), self_merged.limits_at_end_points())
        is_rational, _ = is_all_QQ(flatten(data))
        if not is_rational:
            logging.warn("For functions with non-rational data, cannot guarantee a stable SHA-1 hash.")
        stable_str = str(data)
        return sha1(stable_str).hexdigest()

def singleton_piece(x, y):
    return (singleton_interval(x), FastLinearFunction(0, y))

def open_piece(p, q):
    return (open_interval(p[0], q[0]), linear_function_through_points(p, q))

def closed_piece(p, q):
    return (closed_interval(p[0], q[0]), linear_function_through_points(p, q))

def left_open_piece(p, q):
    return (left_open_interval(p[0], q[0]), linear_function_through_points(p, q))

def right_open_piece(p, q):
    return (right_open_interval(p[0], q[0]), linear_function_through_points(p, q))

        
def print_sign(epsilon):
    if epsilon > 0:
        return "+"
    elif epsilon < 0:
        return "-"
    else:
        return ""

default_precision = 53

default_field = RealNumberField   # can set to SR instead to keep fully symbolic

def can_coerce_to_QQ(x):
    try:
        QQ(x)
        return True
    except ValueError:
        pass
    except TypeError:
        pass
    return False

def is_all_QQ(values):
    is_rational = False
    try:
        values = [ QQ(x) for x in values ]
        is_rational = True
    except ValueError:
        pass
    except TypeError:
        pass
    return is_rational, values

def nice_field_values(symb_values, field=None):
    """
    Coerce the real numbers in the list `symb_values` into a convenient common field
    and return a list, parallel to `symb_values`, of the coerced values.

    If all given numbers are rational, the field will be the rational
    field (`QQ`).  

    Otherwise, if the numbers are algebraic, the field
    will be a suitable algebraic field extension of the rational
    numbers, embedded into the real numbers, in the form of a
    `RealNumberField`.  

    Otherwise, the given numbers are returned as is.
    """
    ### Add tests!
    if field is None:
        field = default_field
    is_rational, field_values = is_all_QQ(symb_values)
    if is_rational:
        logging.info("Rational case.")
        return field_values
    is_realnumberfield, field_values = is_all_the_same_real_number_field(symb_values)
    if is_realnumberfield:
        return field_values
    if field == RealNumberField and not is_rational and not is_realnumberfield:
        # Try to make it a RealNumberField:
        try:
            all_values = [ AA(x) for x in symb_values ]
            #global number_field, number_field_values, morphism, exact_generator, embedded_field, embedding_field, hom, embedded_field_values
            number_field, number_field_values, morphism = number_field_elements_from_algebraics(all_values)
            # Now upgrade to a RealNumberField
            exact_generator = morphism(number_field.gen(0))
            # Use our own RealNumberField.
            symbolic_generator = SR(exact_generator)  # does not quite work --> we won't recover our nice symbolic expressions that way
            if number_field.polynomial().degree() == 2:
                embedding_field = RR  # using a RIF leads to strange infinite recursion
            else:
                embedding_field = RealIntervalField(default_precision)
            embedded_generator = embedding_field(exact_generator)
            embedded_field = RealNumberField(number_field.polynomial(), number_field.variable_name(), \
                                             embedding=embedded_generator, exact_embedding=symbolic_generator)
            hom = number_field.hom([embedded_field.gen(0)])
            embedded_field_values = map(hom, number_field_values)
            # Store symbolic expression
            for emb, symb in itertools.izip(embedded_field_values, symb_values):
                if symb in SR and type(emb) == RealNumberFieldElement:
                    emb._symbolic = symb
            # Transform given data
            field_values = embedded_field_values
            logging.info("Coerced into real number field: %s" % embedded_field)
        except ValueError:
            logging.info("Coercion to a real number field failed, keeping it symbolic")
            pass
        except TypeError:
            logging.info("Coercion to a real number field failed, keeping it symbolic")
            pass
    return field_values

#@logger
def piecewise_function_from_breakpoints_slopes_and_values(bkpt, slopes, values, field=None, merge=True):
    """
    Create a continuous piecewise function from `bkpt`, `slopes`, and `values`.

    `bkpt` and `values` are two parallel lists; it is assumed that `bkpt` is 
    sorted in increasing order. 

    `slopes` is one element shorter and represents the slopes of the interpolation.

    The function is overdetermined by these data.  The consistency of the data is 
    currently not checked.

    The data are coerced into a common convenient field via `nice_field_values`.

    If `merge` is True (the default), adjacent pieces of equal slopes are merged into one.
    """
    if field is None:
        field = default_field
    # global symb_values
    symb_values = bkpt + slopes + values
    field_values = nice_field_values(symb_values, field)
    bkpt, slopes, values = field_values[0:len(bkpt)], field_values[len(bkpt):len(bkpt)+len(slopes)], field_values[-len(values):]
    intercepts = [ values[i] - slopes[i]*bkpt[i] for i in range(len(slopes)) ]
    # Make numbers nice
    ## slopes = [ canonicalize_number(slope) for slope in slopes ]
    ## intercepts = [ canonicalize_number(intercept) for intercept in intercepts ]
    #print slopes
    return FastPiecewise([ [(bkpt[i],bkpt[i+1]), 
                            fast_linear_function(slopes[i], intercepts[i])] for i in range(len(bkpt)-1) ],
                         merge=merge)

def piecewise_function_from_breakpoints_and_values(bkpt, values, field=None, merge=True):
    """
    Create a continuous piecewise function from `bkpt` and `values`.

    `bkpt` and `values` are two parallel lists; assuming `bpkt` is sorted (increasing).

    The data are coerced into a common convenient field via `nice_field_values`.

    If `merge` is True (the default), adjacent pieces of equal slopes are merged into one.
    """
    if len(bkpt)!=len(values):
        raise ValueError, "Need to have the same number of breakpoints and values."
    slopes = [ (values[i+1]-values[i])/(bkpt[i+1]-bkpt[i]) for i in range(len(bkpt)-1) ]
    return piecewise_function_from_breakpoints_slopes_and_values(bkpt, slopes, values, field, merge=merge)

def piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None, merge=True):
    """
    Create a continuous piecewise function from `bkpt` and `slopes`.

    `bkpt` and `slopes` are two parallel lists (except that `bkpt` is
    one element longer); assuming `bpkt` is sorted (increasing).  The
    function always has value 0 on bkpt[0].  

    The data are coerced into a common convenient field via `nice_field_values`.

    If `merge` is True (the default), adjacent pieces of equal slopes are merged into one.
    """
    if len(bkpt)!=len(slopes)+1:
        raise ValueError, "Need to have one breakpoint more than slopes."
    values = [0]
    for i in range(1,len(bkpt)-1):
        values.append(values[i-1] + slopes[i - 1] * (bkpt[i] - bkpt[i-1]))
    return piecewise_function_from_breakpoints_slopes_and_values(bkpt, slopes, values, field, merge=merge)

def piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes, field=None, merge=True):
    """
    Create a continuous piecewise function from `interval_lengths` and `slopes`.

    The function always has value 0 on 0. `interval_lengths` and
    `slopes` are two parallel lists that define the function values to
    the right of 0.

    The data are coerced into a common convenient field via `nice_field_values`.

    If `merge` is True (the default), adjacent pieces of equal slopes are merged into one.
    """
    if len(interval_lengths)!=len(slopes):
        raise ValueError, "Number of given interval_lengths and slopes needs to be equal."
    bkpt = []
    bkpt.append(0)
    for i in range(len(interval_lengths)):
        if interval_lengths[i] < 0:
            raise ValueError, "Interval lengths must be non-negative."
        bkpt.append(bkpt[i]+interval_lengths[i])
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field, merge=merge)

def discrete_function_from_points_and_values(points, values, field=None):
    """
    Create a function defined on a finite list of `points`. 

    `points` and `values` are two parallel lists.

    The data are coerced into a common convenient field via `nice_field_values`.
    """
    if field is None:
        field = default_field
    # global symb_values
    symb_values = points + values
    field_values = nice_field_values(symb_values, field)
    points, values = field_values[0:len(points)], field_values[-len(values):]
    pieces = [ (singleton_interval(point), FastLinearFunction(0, value))
               for point, value in itertools.izip(points, values) ]
    return FastPiecewise(pieces)

def limiting_slopes(fn):
    """
    Compute the limiting slopes on the right and the left side of the
    origin.
    
    The function `fn` is assumed minimal.

    EXAMPLES::

        sage: logging.disable(logging.WARN) # Suppress output in automatic tests.
        sage: limiting_slopes(gmic(f=4/5))
        (5/4, -5)
        sage: limiting_slopes(gmic_disjoint_with_singletons(f=4/5))
        (5/4, -5)
        sage: limiting_slopes(minimal_no_covered_interval())
        (+Infinity, -Infinity)
        sage: limiting_slopes(drlm_2_slope_limit_1_1())
        (2, -Infinity)
        sage: limiting_slopes(restrict_to_finite_group(gmic(f=4/5)))
        (5/4, -5)
    """
    breakpoints = fn.end_points()
    limits = fn.limits_at_end_points()
    assert breakpoints[0] == 0
    if limits[0][0] > 0 or limits[0][1] > 0:
        limit_plus = +Infinity
    elif limits[1][-1] is not None:
        limit_plus = limits[1][-1] / breakpoints[1]
    else:
        limit_plus = limits[1][0] / breakpoints[1]
    assert breakpoints[-1] == 1
    if limits[-1][0] > 0 or limits[-1][-1] > 0:
        limit_minus = -Infinity
    elif limits[-2][+1] is not None:
        limit_minus = -limits[-2][+1] / (1 - breakpoints[-2])
    else:
        limit_minus = -limits[-2][0] / (1 - breakpoints[-2])
    return limit_plus, limit_minus

maximal_asymmetric_peaks_around_orbit = 'maximal_asymmetric_peaks_around_orbit'
maximal_symmetric_peaks_around_orbit = 'maximal_symmetric_peaks_around_orbit'
narrow_symmetric_peaks_around_orbit = 'narrow_symmetric_peaks_around_orbit'
recentered_symmetric_peaks = 'recentered_symmetric_peaks'
recentered_peaks_with_slopes_proportional_to_limiting_slopes_for_positive_epsilon = 'recentered_peaks_with_slopes_proportional_to_limiting_slopes_for_positive_epsilon'
recentered_peaks_with_slopes_proportional_to_limiting_slopes_for_negative_epsilon = 'recentered_peaks_with_slopes_proportional_to_limiting_slopes_for_negative_epsilon'

default_perturbation_style = maximal_asymmetric_peaks_around_orbit

def approx_discts_function(perturbation_list, stability_interval, field=default_field, \
                           perturbation_style=default_perturbation_style, function=None):
    """
    Construct a function that has peaks of +/- 1 around the points of the orbit.
    perturbation_list actually is a dictionary.
    """
    perturb_points = sorted(perturbation_list.keys())
    fn_values = [0]
    fn_bkpt = [0]
    # This width is chosen so that the peaks are disjoint, and
    # so a nice continuous piecewise linear function is constructed.
    if perturbation_style==maximal_asymmetric_peaks_around_orbit or perturbation_style==maximal_symmetric_peaks_around_orbit or narrow_symmetric_peaks_around_orbit:
        width = min(abs(stability_interval.a),stability_interval.b)
        assert width > 0, "Width of stability interval should be positive"
        assert stability_interval.a < 0 < stability_interval.b, \
            "Stability interval should contain 0 in it s interior"
    for pt in perturb_points:
        sign = perturbation_list[pt][0] # the "walk_sign" (character) at the point
        if perturbation_style==maximal_asymmetric_peaks_around_orbit:
            if sign == 1:
                left = pt + stability_interval.a
                right = pt + stability_interval.b
            else:
                left = pt - stability_interval.b
                right = pt - stability_interval.a
        elif perturbation_style==maximal_symmetric_peaks_around_orbit:
            left = pt - width
            right = pt + width
        elif perturbation_style==narrow_symmetric_peaks_around_orbit:
            left = pt - width/1000
            right = pt + width/1000
        elif perturbation_style==recentered_symmetric_peaks:
            if sign == 1:
                left = pt + stability_interval.a
                right = pt + stability_interval.b
            else:
                left = pt - stability_interval.b
                right = pt - stability_interval.a
            pt = (left + right) /2
        elif perturbation_style==recentered_peaks_with_slopes_proportional_to_limiting_slopes_for_positive_epsilon:
            if function is None:
                raise ValueError, "This perturbation_style needs to know function"
            slope_plus, slope_minus = limiting_slopes(function)
            current_slope = function.which_function(pt + (stability_interval.b + stability_interval.a)/2)._slope 
            x = (stability_interval.b - stability_interval.a) * (slope_minus - current_slope)/(slope_minus-slope_plus)
            if sign == 1:
                left = pt + stability_interval.a
                right = pt + stability_interval.b
                pt = left + x
            else:
                left = pt - stability_interval.b
                right = pt - stability_interval.a
                pt = right - x
        elif perturbation_style==recentered_peaks_with_slopes_proportional_to_limiting_slopes_for_negative_epsilon:
            if function is None:
                raise ValueError, "This perturbation_style needs to know function"
            slope_plus, slope_minus = limiting_slopes(function)
            current_slope = function.which_function(pt + (stability_interval.b + stability_interval.a)/2)._slope 
            x = (stability_interval.b - stability_interval.a) * (slope_plus - current_slope)/(slope_plus-slope_minus)
            if sign == 1:
                left = pt + stability_interval.a
                right = pt + stability_interval.b
                pt = left + x
            else:
                left = pt - stability_interval.b
                right = pt - stability_interval.a
                pt = right - x
        else:
            raise ValueError, "Unknown perturbation_style: %s" % perturbation_style
        assert (left >= fn_bkpt[len(fn_bkpt)-1])
        if (left > fn_bkpt[len(fn_bkpt)-1]):
            fn_bkpt.append(left)
            fn_values.append(0)
        fn_bkpt.append(pt)
        fn_values.append(sign)
        fn_bkpt.append(right)
        fn_values.append(0)
    assert (1 >= fn_bkpt[len(fn_bkpt)-1])
    if (1 > fn_bkpt[len(fn_bkpt)-1]):
        fn_bkpt.append(1)
        fn_values.append(0)
    return piecewise_function_from_breakpoints_and_values(fn_bkpt, fn_values, field)

def merge_bkpt(bkpt1, bkpt2):
    i = 0
    j = 0
    bkpt_new = []
    while i < len(bkpt1) and j < len(bkpt2):
        if bkpt1[i] > bkpt2[j]:
            bkpt_new.append(bkpt2[j])
            j = j + 1
        elif bkpt1[i] < bkpt2[j]:
            bkpt_new.append(bkpt1[i])
            i = i + 1
        else:
            bkpt_new.append(bkpt1[i])
            i = i + 1
            j = j + 1
    if i == len(bkpt1) and j != len(bkpt2):
        bkpt_new = bkpt_new + bkpt2[j:len(bkpt2)]
    elif i != len(bkpt1) and j == len(bkpt2):
        bkpt_new = bkpt_new + bkpt1[i:len(bkpt1)]
    return bkpt_new

@cached_function
def find_epsilon_interval(fn, perturb):
    if fn.is_continuous() or fn.is_discrete():
        return find_epsilon_interval_continuous(fn, perturb)
    else:
        return find_epsilon_interval_general(fn, perturb)

def find_largest_epsilon(fn, perturb):
    """
    Compute the proper rescaling of a given perturbation function.
    If the largest epsilon is zero, we should try a different perturbation instead.
    """
    minus_epsilon, plus_epsilon = find_epsilon_interval(fn, perturb)
    return min(abs(minus_epsilon), plus_epsilon)

###
### Moves
###

class FunctionalDirectedMove (FastPiecewise):
    # FIXME: At the moment, does not reduce modulo 1, in contrast to old code!
    """
    Return a pieceweise function to represent a functional directed move
    from a list of domain intervals and the functional directed move. 
    """
    def __init__(self, domain_intervals, directed_move):
        function = fast_linear_function(directed_move[0], directed_move[1])
        pieces = [ (interval, function) for interval in domain_intervals ]
        FastPiecewise.__init__(self, pieces)
        self.directed_move = directed_move       # needed?

    def __repr__(self):
        return "<FunctionalDirectedMove %s with domain %s, range %s>" % (self.directed_move, self.intervals(), self.range_intervals())

    def sign(self):
        """
        Return the sign of the move

        EXAMPLES::

            sage: h = FunctionalDirectedMove([[0.3, 0.4]], (1,0))
            sage: h.sign()
            1
        """
        return self.directed_move[0]

    def is_functional(self):
        return True

    def __getitem__(self, item):
        return self.directed_move[item]

    def can_apply(self, x):
        """
        Determine if self can apply on x

        EXAMPLES::

            sage: h = FunctionalDirectedMove([[0.3, 0.4], [0.58, 0.68]], (1,0))
            sage: h.can_apply(0.3)
            True
            sage: h.can_apply(0.2)
            False            
        """
        try:
            self(x)
            return True
        except ValueError:
            return False

    def apply_ignoring_domain(self, x):
        """
        Apply self on x by ignoring the domain (use modulo 1)

        EXAMPLES::

            sage: h = FunctionalDirectedMove([[0.3, 0.4], [0.58, 0.68]], (1,0))
            sage: h.apply_ignoring_domain(1/10)
            1/10
            sage: h = FunctionalDirectedMove([[0.1, 0.6]], (-1,1))
            sage: h.apply_ignoring_domain(1/2)
            1/2

        """
        move_sign = self.sign()
        if move_sign == 1:
            next_x = fractional(x + self.directed_move[1])
        elif move_sign == -1:
            next_x = fractional(self.directed_move[1]-x)
        return next_x

    def apply_to_coho_interval(self, interval, inverse=False):
        # This does not do error checking.  Some code depends on this fact!
        # FIXME: This should be made clear in the name of this function.
        """
        Return a range inverval from a given interval by applying the move. 
        If the move sign is 1, the user can take the inverse of the operation,
        i.e y = x - t_1

        EXAMPLES::

            sage: h = FunctionalDirectedMove([[0.3, 0.4]], (-1, 1))
            sage: h.apply_to_coho_interval([1/10, 1/2])
            <Int[1/2, 9/10]>
            sage: h = FunctionalDirectedMove([[0.3, 0.4]], (1, 1/10))
            sage: h.apply_to_coho_interval([1/10, 1/2])
            <Int[1/5, 3/5]>
            sage: h.apply_to_coho_interval([1/10, 1/2], inverse=True)
            <Int[0, 2/5]>
        """
        if len(interval) <= 2:
            interval = coho_interval_from_interval(interval) # FIXME: Can be removed if FastPiecewise exclusively uses coho intervals.
        directed_move = self.directed_move
        move_sign = directed_move[0]
        if move_sign == 1:
            if inverse:
                result = closed_or_open_or_halfopen_interval(interval[0] - directed_move[1], interval[1] - directed_move[1], \
                                                             interval.left_closed, interval.right_closed)
            else:
                result = closed_or_open_or_halfopen_interval(interval[0] + directed_move[1], interval[1] + directed_move[1], \
                                                             interval.left_closed, interval.right_closed)
        elif move_sign == -1:
            result = closed_or_open_or_halfopen_interval(directed_move[1] - interval[1], directed_move[1] - interval[0], \
                                                         interval.right_closed, interval.left_closed)
        else:
            raise ValueError, "Move not valid: %s" % list(move)
        return result

    def range_intervals(self):
        """
        Return the range intervals of self.
        """
        return [ self.apply_to_coho_interval(interval) for interval in self.intervals() ] 

    def is_identity(self):
        """
        Determine if self is a identity function or not

        EXAMPLES::

            sage: h = FunctionalDirectedMove([[0.3, 0.4]], (1, 0))
            sage: h.is_identity()
            True
            sage: h = FunctionalDirectedMove([[0.3, 0.4]], (-1, 1))
            sage: h.is_identity()
            False
        """
        return self.directed_move[0] == 1 and self.directed_move[1] == 0

    def restricted(self, intervals):
        """
        Not used.
        Return a new move that is the restriction of domain and codomain of `self` to `intervals`.
        (The result may have the empty set as its domain.)
        """
        domain = self.intervals()                        # sorted.
        preimages = [ self.apply_to_coho_interval(interval, inverse=True) for interval in intervals ]
        preimages.sort(key=coho_interval_left_endpoint_with_epsilon)
        new_domain = list(intersection_of_coho_intervals([domain, intervals, preimages]))
        return FunctionalDirectedMove(new_domain, self.directed_move)

    def restricting(self, components):
        """ 
        Return a new move by removing self.restricted(component) for component in components.
        (The result may have the empty set as its domain.)
        """

        domain = self.intervals()                        # sorted.
        restricting_domain_list = []
        for component in components:
            preimages = [ self.apply_to_coho_interval(interval, inverse=True) for interval in component ]
            preimages.sort(key=coho_interval_left_endpoint_with_epsilon)
            restricting_domain_list.append(list(intersection_of_coho_intervals([component, preimages])))
        new_domain = union_of_coho_intervals_minus_union_of_coho_intervals([domain], restricting_domain_list, remove_closure=True )
        return FunctionalDirectedMove(new_domain, self.directed_move)

    def plot(self, *args, **kwds):
        kwds = copy(kwds)
        kwds['discontinuity_markers'] = False
        return FastPiecewise.plot(self, *args, **kwds)

@cached_function
def generate_functional_directed_moves(fn):
    """
    Compute the moves (translations and reflections).
    """
    moves = dict()
    for face in generate_maximal_additive_faces(fn):
        if face.is_directed_move():
            fdm = face.functional_directed_move()
            if fdm.intervals():
                if fdm.directed_move in moves:
                    moves[fdm.directed_move] = merge_functional_directed_moves(moves[fdm.directed_move], fdm)
                else:
                    moves[fdm.directed_move] = fdm
    return list(moves.values())

def plot_walk(walk_dict, color="black", ymin=0, ymax=1, **kwds):
    #return point([ (x,0) for x in walk_dict.keys()])
    g = Graphics()
    kwds['legend_label'] = "reachable orbit"
    for x in walk_dict.keys():
        g += line([(x,ymin), (x,ymax)], color=color, zorder=-4, **kwds)
        delete_one_time_plot_kwds(kwds)
    return g

class UnimplementedError (Exception):
    pass

def generate_symbolic(fn, components, field=None):
    if fn.is_continuous() or fn.is_discrete():
        return generate_symbolic_continuous(fn, components, field=field)
    else:
        return generate_symbolic_general(fn, components, field=field)

def generate_additivity_equations(fn, symbolic, field, f=None, bkpt=None):
    if fn.is_continuous() or fn.is_discrete():
        return generate_additivity_equations_continuous(fn, symbolic, field, f=f, bkpt=bkpt)
    else:
        return generate_additivity_equations_general(fn, symbolic, field, f=f, bkpt=bkpt)

def rescale_to_amplitude(perturb, amplitude):
    """For plotting purposes, rescale the function `perturb` so that its
    maximum (supremum) absolute function value is `amplitude`.
    """
    current_amplitude = max([ abs(x) for limits in perturb.limits_at_end_points() for x in limits if x is not None])
    if current_amplitude != 0:
        return perturb * (amplitude/current_amplitude)
    else:
        return perturb

# Global figsize for all plots made by show_plots.
show_plots_figsize = 10

def show_plot(graphics, show_plots, tag, object=None, **show_kwds):
    """
    Display or save `graphics`.

    `show_plots` can be one of: `False` (do nothing), 
    `True` (use `show` to display on screen),
    a string (file name format such as "FILENAME-%s.pdf", 
    where %s is replaced by `tag`.
    """
    plot_kwds_hook(show_kwds)
    if isinstance(show_plots, str):
        graphics.save(show_plots % tag, figsize=show_plots_figsize, **show_kwds)
    elif show_plots:
        graphics.show(figsize=show_plots_figsize, **show_kwds)

def plot_rescaled_perturbation(perturb, xmin=0, xmax=1, **kwds):
    return plot(rescale_to_amplitude(perturb, 1/10), xmin=xmin,
                xmax=xmax, color='magenta', legend_label="perturbation (rescaled)", **kwds)

check_perturbation_plot_three_perturbations = True

def basic_perturbation(fn, index):
    """
    Get a basic perturbation of `fn`.  `index` counts from 1 (to match the labels in the diagrams). 
    """
    if not hasattr(fn, '_perturbations'):
        extremality_test(fn, show_plots=False)
    if hasattr(fn, '_perturbations'):
        try: 
            return fn._perturbations[index-1]
        except IndexError:
            raise IndexError, "Bad perturbation index"
    raise ValueError, "No perturbations"

def plot_perturbation_diagram(fn, perturbation=None, xmin=0, xmax=1):
    """
    Plot a perturbation of `fn`.
    
    `perturbation` is either a perturbation function, or an integer
    (which designates a basic perturbation of `fn` via
    `basic_perturbation`).  If `perturbation` is not provided, it
    defaults to the perturbation indexed 1.

    To show only a part of the diagram, use::

        sage: plot_perturbation_diagram(h, 1).show(xmin=0.25, xmax=0.35, ymin=0.25, ymax=0.35)  # not tested
    """
    if perturbation is None:
       perturbation = 1
    if isinstance(perturbation, Integer):
        perturbation = basic_perturbation(fn, perturbation)
    epsilon_interval = find_epsilon_interval(fn, perturbation)
    epsilon = min(abs(epsilon_interval[0]), epsilon_interval[1])
    p = plot_rescaled_perturbation(perturbation, xmin=xmin, xmax=xmax)
    if check_perturbation_plot_three_perturbations:
        p += plot(fn + epsilon_interval[0] * perturbation, xmin=xmin, xmax=xmax, color='red', legend_label="-perturbed (min)")
        p += plot(fn + epsilon_interval[1] * perturbation, xmin=xmin, xmax=xmax, color='blue', legend_label="+perturbed (max)")
        if -epsilon != epsilon_interval[0]:
            p += plot(fn + (-epsilon) * perturbation, xmin=xmin, xmax=xmax, color='orange', legend_label="-perturbed (matches max)")
        elif epsilon != epsilon_interval[1]:
            p += plot(fn + epsilon * perturbation, xmin=xmin, xmax=xmax, color='cyan', legend_label="+perturbed (matches min)")
    else:
        label = "-perturbed"
        if -epsilon == epsilon_interval[0]:
            label += " (min)"
        else:
            label += " (matches max)"
        p += plot(fn - epsilon * perturbation, xmin=xmin, xmax=xmax, color='red', legend_label=label)
        label = "+perturbed"
        if epsilon == epsilon_interval[1]:
            label += " (max)"
        else:
            label += " (matches min)"
        p += plot(fn + epsilon * perturbation, xmin=xmin, xmax=xmax, color='blue', legend_label=label)
    p += plot(fn, xmin=xmin, xmax=xmax, color='black', thickness=2,
             legend_label="original function", **ticks_keywords(fn))
    return p

def check_perturbation(fn, perturb, show_plots=False, show_plot_tag='perturbation', xmin=0, xmax=1, **show_kwds):
    epsilon_interval = find_epsilon_interval(fn, perturb)
    epsilon = min(abs(epsilon_interval[0]), epsilon_interval[1])
    #logging.info("Epsilon for constructed perturbation: %s" % epsilon)
    if show_plots:
        logging.info("Plotting perturbation...")
        p = plot_perturbation_diagram(fn, perturb, xmin=xmin, xmax=xmax)
        show_plot(p, show_plots, tag=show_plot_tag, object=fn, **show_kwds)
        logging.info("Plotting perturbation... done")
    assert epsilon > 0, "Epsilon should be positive, something is wrong"
    #logging.info("Thus the function is not extreme.")  ## Now printed by caller.

def generate_perturbations_finite_dimensional(function, show_plots=False, f=None):
    ## FIXME: Perhaps we want an `oversampling` parameter as in generate_perturbations_simple??
    """
    Generate (with "yield") perturbations for `finite_dimensional_extremality_test`.
    """
    fdms, covered_components = generate_directed_move_composition_completion(function, show_plots=show_plots)
    uncovered_intervals = generate_uncovered_intervals(function)
    if uncovered_intervals:
        uncovered_components = generate_uncovered_components(function, show_plots=show_plots)
        components = covered_components + uncovered_components
    else:
        components = copy(covered_components)
    # FIXME: fraction_field() required because parent could be Integer
    # Ring.  This happens, for example, for three_slope_limit().  
    # We really should have a function to retrieve the field of
    # a FastPiecewise.  But now even .base_ring() fails because
    # FastLinearFunction does not have a .base_ring() method.
    field = function(0).parent().fraction_field()
    symbolic = generate_symbolic(function, components, field=field)
    bkpt = merge_bkpt(function.end_points(), symbolic.end_points())
    equation_matrix = generate_additivity_equations(function, symbolic, field, f=f, bkpt=bkpt)
    slope_jump_vects = equation_matrix.right_kernel().basis()
    logging.info("Finite dimensional test: Solution space has dimension %s" % len(slope_jump_vects))
    for basis_index in range(len(slope_jump_vects)):
        slope_jump = slope_jump_vects[basis_index]
        perturbation = slope_jump * symbolic
        yield perturbation

def finite_dimensional_extremality_test(function, show_plots=False, f=None, warn_about_uncovered_intervals=True, 
                                        show_all_perturbations=False):
    """
    Solve a homogeneous linear system of additivity equations with one
    slope variable for every component (including every non-covered
    interval) and one jump variable for each (left/right) discontinuity.

    Return a boolean that indicates whether the system has a nontrivial solution.

    EXAMPLES::

        sage: logging.disable(logging.WARN)
        sage: h1 = drlm_not_extreme_2()
        sage: finite_dimensional_extremality_test(h1, show_plots=True)
        False
        sage: h2 = drlm_3_slope_limit()
        sage: finite_dimensional_extremality_test(h2, show_plots=True)
        True
    """
    if show_all_perturbations is None:
        show_all_perturbations = show_plots
    if function.is_discrete():
        return simple_finite_dimensional_extremality_test(function, oversampling=1, show_all_perturbations=show_all_perturbations)
    seen_perturbation = False
    function._perturbations = []
    for index, perturbation in enumerate(generate_perturbations_finite_dimensional(function, show_plots=show_plots, f=f)):
        function._perturbations.append(perturbation)
        check_perturbation(function, perturbation,
                           show_plots=show_plots, show_plot_tag='perturbation-%s' % (index + 1),
                           legend_title="Basic perturbation %s" % (index + 1))
        if not seen_perturbation:
            seen_perturbation = True
            logging.info("Thus the function is NOT extreme.")
            if not show_all_perturbations:
                break
    if not seen_perturbation:
        logging.info("Finite dimensional extremality test did not find a perturbation.")
        uncovered_intervals = generate_uncovered_intervals(function)
        if uncovered_intervals:
            if warn_about_uncovered_intervals:
                logging.warn("There are non-covered intervals, so this does NOT prove extremality.")
        else:
            logging.info("Thus the function is extreme.")
    return not seen_perturbation

def generate_type_1_vertices(fn, comparison, reduced=True, bkpt=None):
    if fn.is_continuous() or fn.is_discrete():
        return generate_type_1_vertices_continuous(fn, comparison, bkpt=bkpt)
    else:
        return generate_type_1_vertices_general(fn, comparison, reduced=reduced, bkpt=bkpt)

def generate_type_2_vertices(fn, comparison, reduced=True, bkpt=None):
    if fn.is_continuous() or fn.is_discrete():
        return generate_type_2_vertices_continuous(fn, comparison, bkpt=bkpt)
    else:
        return generate_type_2_vertices_general(fn, comparison, reduced=reduced, bkpt=bkpt)

def generate_additive_vertices(fn, reduced=True, bkpt=None):
    """
    We are returning a set of 6-tuples (x, y, z, xeps, yeps, zeps),
    so that duplicates are removed, and so the result can be cached for later use.

    When reduced=True:
        only outputs fewer triples satisfying `comparison' relation, for the purpose of setting up the system of equations.

    When reduced=False:
        outputs all triples satisfying `comparison' relation, for the purpose of plotting additive_limit_vertices.
    """
    return set(itertools.chain( \
                generate_type_1_vertices(fn, operator.eq, reduced=reduced, bkpt=bkpt),\
                generate_type_2_vertices(fn, operator.eq, reduced=reduced, bkpt=bkpt)) )

@cached_function
def generate_nonsubadditive_vertices(fn, reduced=True):
    """
    We are returning a set of 6-tuples (x, y, z, xeps, yeps, zeps),
    so that duplicates are removed, and so the result can be cached for later use.

    When reduced=True:
        only outputs fewer triples satisfying `comparison' relation, for the purpose of minimality_test.

    When reduced=False:
        outputs all triples satisfying `comparison' relation, for the purpose of plotting nonsubadditive_limit_vertices.
    """
    return set(itertools.chain( \
                generate_type_1_vertices(fn, operator.lt, reduced=reduced),\
                generate_type_2_vertices(fn, operator.lt, reduced=reduced))  )

def generate_nonsymmetric_vertices(fn, f):
    if fn.is_continuous() or fn.is_discrete():
        return generate_nonsymmetric_vertices_continuous(fn, f)
    else:
        return generate_nonsymmetric_vertices_general(fn, f)

class MaximumNumberOfIterationsReached(Exception):
    pass

def extremality_test(fn, show_plots = False, f=None, max_num_it = 1000, perturbation_style=default_perturbation_style, phase_1 = False, finite_dimensional_test_first = False, show_all_perturbations=False):
    """Check if `fn` is extreme for the group relaxation with the given `f`. 

    If `fn` is discrete, it has to be defined on a cyclic subgroup of
    the reals containing 1, restricted to [0, 1].  The group
    relaxation is the corresponding cyclic group relaxation.

    Otherwise `fn` needs to be defined on the interval [0, 1], and the
    group relaxation is the infinite group relaxation.

    If `f` is not provided, uses the one found by `find_f()`.

    If `show_plots` is True (default: False), show many illustrating diagrams.

    The function first runs `minimality_test`.
    
    In the infinite group case, if `finite_dimensional_test_first` is
    True (default: False), after testing minimality of `fn`, we first
    check if the `finite_dimensional_extremality_test` finds a
    perturbation; otherwise (default) we first check for an
    equivariant perturbation.

    EXAMPLES::

        sage: logging.disable(logging.INFO) # to disable output in automatic tests.
        sage: h = piecewise_function_from_breakpoints_and_values([0, 1/2, 1], [0, 1, 0])
        sage: # This example has a unique candidate for "f", so we don't need to provide one.
        sage: extremality_test(h, False)
        True
        sage: # Same, with plotting:
        sage: extremality_test(h, True) # not tested
        ... lots of plots shown ...
        True
        sage: h = multiplicative_homomorphism(gmic(f=4/5), 3) 
        sage: # This example has several candidates for "f", so provide the one we mean:
        sage: extremality_test(h, True, f=4/15) # not tested
        ... lots of plots shown ...
        True
        sage: g = gj_2_slope()
        sage: gf = restrict_to_finite_group(g)
        sage: # This is now a finite (cyclic) group problem.
        sage: extremality_test(gf, True) # not tested
        ... lots of plots shown ...
        True
    """
    if show_all_perturbations is None:
        show_all_perturbations = show_plots
    do_phase_1_lifting = False
    if f is None:
        f = find_f(fn, no_error_if_not_minimal_anyway=True)
    if f is None or not minimality_test(fn, show_plots=show_plots, f=f):
        logging.info("Not minimal, thus NOT extreme.")
        if not phase_1:
            return False
        else:
            do_phase_1_lifting = True
    if do_phase_1_lifting:
        finite_dimensional_test_first = True
    seen_perturbation = False
    generator = generate_perturbations(fn, show_plots=show_plots, f=f, max_num_it=max_num_it, finite_dimensional_test_first=finite_dimensional_test_first, perturbation_style=perturbation_style)
    fn._perturbations = []
    for index, perturbation in enumerate(generator):
        fn._perturbations.append(perturbation)
        check_perturbation(fn, perturbation, show_plots=show_plots, 
                           show_plot_tag='perturbation-%s' % (index + 1), 
                           legend_title="Basic perturbation %s" % (index + 1))
        if not seen_perturbation:
            seen_perturbation = True
            logging.info("Thus the function is NOT extreme.")
            if not show_all_perturbations:
                break
    if not seen_perturbation:
        logging.info("Thus the function is extreme.")
    return not seen_perturbation

def generate_perturbations(fn, show_plots=False, f=None, max_num_it=1000, perturbation_style=default_perturbation_style, finite_dimensional_test_first = False):
    """
    Generate (with "yield") perturbations for `extremality_test`.
    """
    if fn.is_discrete():
        all = generate_perturbations_simple(fn, show_plots=show_plots, f=f, oversampling=None)
    else:
        fdms, covered_components= generate_directed_move_composition_completion(fn, show_plots=show_plots)
        finite = generate_perturbations_finite_dimensional(fn, show_plots=show_plots, f=f)
        uncovered_intervals = generate_uncovered_intervals(fn)
        if show_plots:
            logging.info("Plotting covered intervals...")
            show_plot(plot_covered_intervals(fn), show_plots, tag='covered_intervals', object=fn)
            logging.info("Plotting covered intervals... done")
        if not uncovered_intervals:
            logging.info("All intervals are covered (or connected-to-covered). %s components." % number_of_components(fn))
            all = finite
        else:
            logging.info("Uncovered intervals: %s", (uncovered_intervals,))
            equi = generate_perturbations_equivariant(fn, show_plots=show_plots, f=f, max_num_it=max_num_it, perturbation_style=perturbation_style)
            if finite_dimensional_test_first:
                all = itertools.chain(finite, equi)
            else:
                all = itertools.chain(equi, finite)
    for perturbation in all:
        yield perturbation

def generate_perturbations_equivariant(fn, show_plots=False, f=None, max_num_it=1000, perturbation_style=default_perturbation_style):
    if not fn.is_continuous():
        logging.warning("Code for detecting perturbations using moves is EXPERIMENTAL in the discontinuous case.")
    generator = generate_generic_seeds_with_completion(fn, show_plots=show_plots, max_num_it=max_num_it) # may raise MaximumNumberOfIterationsReached
    #seen_perturbation = False
    for seed, stab_int, walk_list in generator:
        # for debugging only:
        #global last_seed, last_stab_int, last_walk_list = seed, stab_int, walk_list
        perturb = approx_discts_function(walk_list, stab_int, perturbation_style=perturbation_style, function=fn)
        perturb._seed = seed
        perturb._stab_int = stab_int
        perturb._walk_list = walk_list
        if show_plots:
            logging.info("Plotting completion diagram with perturbation...")
            g = plot_completion_diagram(fn, perturb)        # at this point, the perturbation has not been stored yet
            show_plot(g, show_plots, tag='completion', object=fn._completion, legend_title="Completion of moves, perturbation", legend_loc="upper left")
            logging.info("Plotting completion diagram with perturbation... done")
        #seen_perturbation = True
        yield perturb
    #if not seen_perturbation:
    #    logging.info("Dense orbits in all non-covered intervals.")

def plot_completion_diagram(fn, perturbation=None):
    """
    Return a plot of the completion diagram.
    
    To view a part only, use::

        sage: plot_completion_diagram(h).show(xmin=0.3, xmax=0.55, ymin=0.3, ymax=0.55) # not tested
    """
    if not (hasattr(fn, '_completion') and fn._completion.is_complete):
        extremality_test(fn, show_plots=False)
    if fn._completion.plot_background is None:
        fn._completion.plot_background = plot_completion_diagram_background(fn)
    g = fn._completion.plot() 
    if perturbation is None:
        if hasattr(fn, '_perturbations') and fn._perturbations:
            perturbation = fn._perturbations[0]
    elif isinstance(perturbation, Integer):
        perturbation = basic_perturbation(fn, perturbation)
    if perturbation is not None:
        g += plot_function_at_borders(rescale_to_amplitude(perturbation, 1/10), color='magenta', legend_label='perturbation (rescaled)')
    if hasattr(perturbation, '_walk_list'):
        g += plot_walk_in_completion_diagram(perturbation._seed, perturbation._walk_list)
    return g

def lift(fn, show_plots = False, which_perturbation = 1, **kwds):
    # FIXME: Need better interface for perturbation selection.
    if not hasattr(fn, '_perturbations') and extremality_test(fn, show_plots=show_plots, **kwds):
        return fn
    else:
        perturbation = fn._perturbations[0]
        epsilon_interval = find_epsilon_interval(fn, perturbation)
        perturbed = fn._lifted = fn + epsilon_interval[which_perturbation] * perturbation
        ## Following is strictly experimental: It may change what "f" is.
        if 'phase_1' in kwds and kwds['phase_1']:
            perturbed = rescale_to_amplitude(perturbed, 1)
        return perturbed

def lift_until_extreme(fn, show_plots = False, pause = False, **kwds):
    next, fn = fn, None
    while next != fn:
        fn = next
        next = lift(fn, show_plots=show_plots, **kwds)
        if pause and next != fn:
            raw_input("Press enter to continue")
    return next

##############
def lift_new(fn, order, show_plots = False, which_perturbation = 1, **kwds):
    # FIXME: Need better interface for perturbation selection.
    if not hasattr(fn, '_perturbations') and simple_finite_dimensional_extremality_test(fn, show_plots=show_plots, order=order):
        return fn
    else:
        perturbation = fn._perturbations[0]
        epsilon_interval = find_epsilon_interval(fn, perturbation)
        perturbed = fn._lifted = fn + epsilon_interval[which_perturbation] * perturbation
        ## Following is strictly experimental: It may change what "f" is.
        if 'phase_1' in kwds and kwds['phase_1']:
            perturbed = rescale_to_amplitude(perturbed, 1)
        return perturbed

def lift_new_until_extreme(fn, show_plots = False, pause = False, first_oversampling = 4, **kwds):
    order = finite_group_order_from_function_f_oversampling_order(fn, oversampling=first_oversampling)
    next = lift_new(fn, order, show_plots, **kwds)
    while next != fn:
        fn = next
        next = lift_new(fn, order, show_plots=show_plots, **kwds)
        if pause and next != fn:
            raw_input("Press enter to continue")
    return next

##############
def last_lifted(fn):
    while hasattr(fn, '_lifted'):
        fn = fn._lifted
    return fn

def piecewise_function_from_robert_txt_file(filename):
    """The .txt files have 4 rows.  
    1st row = Y values
    2nd row = X values (I don't use these, but I included them in case you want them)
    3rd row = f   (the x coordinate for which I use as f)
    4th row = value at f  (I don't normalize this to 1.  This allows the Y values to range from 0 to this values)

    Also, I don't include the last value (pi(1)) ever because this is
    the same as pi(0) due to periodicity.  So, if you need this last
    value, please attach a 0 to the end of the Y values and an extra x
    value.
    """
    with open(filename) as f:
        yvalues = [QQ(x) for x in f.readline().split()]
        xvalues = [QQ(x) for x in f.readline().split()]
        if xvalues != range(len(yvalues)):
            raise ValueError, "Line 2 (xvalues) need to be consecutive integers"
        xscale = len(xvalues)
        xf = QQ(f.readline())
        yf = QQ(f.readline())
    if yvalues[xf] != yf:
        raise ValueError, "Lines 3/4 on f and value at f are not consistent with line 1."
    return piecewise_function_from_breakpoints_and_values([ x / xscale for x in xvalues ] + [1], [y / yf for y in yvalues] + [0])

def random_piecewise_function(xgrid=10, ygrid=10, continuous_proba=1, symmetry=True):
    """
    Return a random, continuous or discontinuous piecewise linear function defined on [0, 1]
    with breakpoints that are multiples of 1/`xgrid` and values that are multiples of 1/`ygrid`.

    `continuous_proba` (a real number in [0,1]) indicates the probability that the function is (left/right) continuous at a breakpoint. 
    Use continuous_proba = 1 (the default) to get a continuous piecewise linear function.

    Use symmetry=True (the default) to get a symmetric function. 

    EXAMPLES::

        sage: h = random_piecewise_function(10, 10)
        sage: h = random_piecewise_function(10, 10, continuous_proba=4/5, symmetry=True)
        sage: h = random_piecewise_function(10, 10, continuous_proba=4/5, symmetry=False)
    """
    xvalues = [0] + [ x/xgrid for x in range(1, xgrid) ] + [1]
    f = randint(1, xgrid - 1)
    left_midpoint = f / 2
    right_midpoint = (f+xgrid) / 2
    #yvalues = [0] + [ randint(0, ygrid) / ygrid for i in range(1, f) ] + [1] + [ randint(0, ygrid) / ygrid for i in range(f+1, xgrid) ]+ [0]
    yvalues = [0] + [ randint(1, ygrid-1) / ygrid for i in range(1, f) ] + [1] + [ randint(1, ygrid-1) / ygrid for i in range(f+1, xgrid) ]+ [0]
    if symmetry:
        for i in range(1, ceil(left_midpoint)):
            yvalues[f-i] = 1 - yvalues[i]
        if left_midpoint in ZZ:
            yvalues[left_midpoint] = 1/2
        for i in range(f+1, ceil(right_midpoint)):
            yvalues[f + xgrid - i] = 1 - yvalues[i]
        if right_midpoint in ZZ:
            yvalues[right_midpoint] = 1/2
    if continuous_proba == 1:
        return piecewise_function_from_breakpoints_and_values(xvalues, yvalues)
    else:
        piece1 = [ [singleton_interval(xvalues[i]), FastLinearFunction(0, yvalues[i])] for i in range(xgrid+1) ]
        leftlimits = [0]
        rightlimits = []
        for i in range(0, ygrid):
            p = random()
            if p > continuous_proba:
                rightlimits.append(randint(0, ygrid) / ygrid)
            else:
                rightlimits.append(yvalues[i])
            p = random()
            if p > continuous_proba:
                leftlimits.append(randint(0, ygrid) / ygrid)
            else:
                leftlimits.append(yvalues[i+1])
        rightlimits.append(0)
        if symmetry:
            for i in range(1, ceil(left_midpoint)):
                leftlimits[f-i] = 1 - rightlimits[i]
                rightlimits[f-i] = 1 - leftlimits[i]
            if left_midpoint in ZZ:
                rightlimits[left_midpoint] = 1 - leftlimits[left_midpoint]
            leftlimits[f] = 1 - rightlimits[0]
            for i in range(f+1, ceil(right_midpoint)):
                leftlimits[f + xgrid - i] = 1 - rightlimits[i]
                rightlimits[f + xgrid - i] = 1 - leftlimits[i]
            if right_midpoint in ZZ:
                rightlimits[right_midpoint] = 1 - leftlimits[right_midpoint]
            leftlimits[xgrid] = 1 - rightlimits[f]
        slopes = [ (leftlimits[i+1] - rightlimits[i]) * xgrid for i in range(0, xgrid) ]
        intercepts = [ rightlimits[i] - xvalues[i] * slopes[i] for i in range(0, xgrid) ]
        piece2 = [ [open_interval(xvalues[i], xvalues[i+1]), FastLinearFunction(slopes[i], intercepts[i])] for i in range(xgrid) ]
        pieces = [piece1[0]]
        for i in range(xgrid):
            pieces += [piece2[i], piece1[i+1]]
        return FastPiecewise(pieces, merge=True)

def is_all_QQ_fastpath(values):
    """
    This version does not do the full check whether it can be coerced to QQ,
    which is slow for RealNumberField.
    """
    for x in values:
        if not isinstance(x, (int, long, Rational, Integer)):
            return False
    return True

from sage.rings.number_field.number_field_element import is_NumberFieldElement

def is_all_the_same_number_field_fastpath(values):
    """
    This version does not try coercions and compares fields using 'is', rather than their comparison operator.
    """
    number_field_seen = None
    for x in values:
        if is_NumberFieldElement(x):
            if number_field_seen:
                if number_field_seen is not x.parent():
                    return False
            else:
                number_field_seen = x.parent()
        else:
            return False
    return True

def is_QQ_linearly_independent(*numbers):
    """
    Test if `numbers` are linearly independent over `QQ`.

    EXAMPLES::

        sage: logging.disable(logging.INFO)  # Suppress output in automatic tests.
        sage: is_QQ_linearly_independent()
        True
        sage: is_QQ_linearly_independent(1)
        True
        sage: is_QQ_linearly_independent(0)
        False
        sage: is_QQ_linearly_independent(1,2)
        False
        sage: is_QQ_linearly_independent(1,sqrt(2))
        True
        sage: is_QQ_linearly_independent(1+sqrt(2),sqrt(2),1)
        False
    """
    # trivial cases
    if len(numbers) == 0:
        return True
    elif len(numbers) == 1:
        return numbers[0] != 0
    # fast path for rationals
    if is_all_QQ_fastpath(numbers):
        return False
    if not is_all_the_same_number_field_fastpath(numbers):
        # try to coerce to common number field
        numbers = nice_field_values(numbers, RealNumberField)
        if not is_NumberFieldElement(numbers[0]):
            if is_all_QQ(numbers):
                return False
            raise ValueError, "Q-linear independence test only implemented for algebraic numbers"
    coordinate_matrix = matrix(QQ, [x.list() for x in numbers])
    return rank(coordinate_matrix) == len(numbers)

def compose_functional_directed_moves(A, B, show_plots=False):
    """
    Compute the directed move that corresponds to the directed move `A` after `B`.
    
    EXAMPLES::

        sage: compose_functional_directed_moves(FunctionalDirectedMove([(5/10,7/10)],(1, 2/10)),FunctionalDirectedMove([(2/10,4/10)],(1,2/10)))
        <FunctionalDirectedMove (1, 2/5) with domain [(3/10, 2/5)], range [<Int[7/10, 4/5]>]>
    """
    assert(A.is_functional() and B.is_functional())
    A_domain_preimages = [ B.apply_to_coho_interval(A_domain_interval, inverse=True) \
                           for A_domain_interval in A.intervals() ]
    A_domain_preimages.sort(key=coho_interval_left_endpoint_with_epsilon)
    result_domain_intervals = list(intersection_of_coho_intervals([A_domain_preimages, B.intervals()])) # generator!
    if result_domain_intervals:
        result = FunctionalDirectedMove(result_domain_intervals, (A[0] * B[0], A[0] * B[1] + A[1]))
    else:
        result = None
    if show_plots:
        p = plot(A, color="green", legend_label="A")
        p += plot(B, color="blue", legend_label="B")
        if result:
            p += plot(result, color="red", legend_label="C = A after B")
        show_plot(p, show_plots, tag='compose_directed_moves')
    return result

def merge_functional_directed_moves(A, B, show_plots=False):
    """
    EXAMPLES::

        sage: merge_functional_directed_moves(FunctionalDirectedMove([(3/10, 7/20), (9/20, 1/2)], (1,0)),FunctionalDirectedMove([(3/10, 13/40)], (1,0)))
        <FunctionalDirectedMove (1, 0) with domain [(3/10, 7/20), (9/20, 1/2)], range [<Int[3/10, 7/20]>, <Int[9/20, 1/2]>]>
    """
    if A.directed_move != B.directed_move:
        raise ValueError, "Cannot merge, moves have different operations"
    C = FunctionalDirectedMove(\
                               A.intervals() + B.intervals(),  # constructor takes care of merging
                               A.directed_move)
    if show_plots:
        p = plot(C, color="cyan", legend_label="C = A merge B", thickness=10)
        p += plot(A, color="green", legend_label="A = %s" % A )
        p += plot(B, color="blue", legend_label="B = %s" % B)
        show_plot(p, show_plots, tag='merge_functional_directed_moves')
    return C

def plot_directed_moves(dmoves, **kwds):
    g = Graphics()
    for dm in dmoves:
        g += plot(dm, **kwds)
        delete_one_time_plot_kwds(kwds)
    return g

def plot_dense_moves(component, **kwds):
    return sum([polygon(((domain[0], domain[0]), (domain[1], domain[0]), (domain[1], domain[1]), (domain[0], domain[1])), rgbcolor=kwds.get("rgbcolor", "cyan"), alpha=0.5) + polygon(((domain[0], domain[0]), (domain[1], domain[0]), (domain[1], domain[1]), (domain[0], domain[1])), color="red", fill=False) for domain in component])


def reduce_covered_components(covered_components):
    reduced_components = []
    remaining_components = covered_components
    while remaining_components:
        given_component = remaining_components[0]
        other_components = remaining_components[1::]
        given_component, remaining_components = merge_components_with_given_component(given_component, other_components)
        reduced_components.append(given_component)
    return reduced_components

def merge_components_with_given_component(given_component, other_components):
    remaining_components = []
    for component in other_components:
        if list(intersection_of_coho_intervals([given_component, component])):
            # if the two components intersect, the intersection is full-dimensional since intervals are all open intervals.
            # extend the given_component to the union.
            given_component = union_of_coho_intervals_minus_union_of_coho_intervals([given_component, component], [])
        else:
            # if component is not merged into the given_component, it remains in remaining_components.
            remaining_components.append(component)
    return given_component, remaining_components

def check_for_strip_lemma(m1, m2):
    dense_intervals = []
    if m1.sign() == 1 and m2.sign() == 1:
        t1 = m1[1]
        t2 = m2[1]
        if is_QQ_linearly_independent(t1, t2) and t1 >= 0 and t2 >= 0:
            for m1i in m1.intervals():
                for m2i in m2.intervals():
                    l1, u1 = m1i[0], m1i[1]
                    l2, u2 = m2i[0], m2i[1]
                    L = max(l1, l2)
                    U = min(u1 + t1, u2 + t2)
                    if t1 + t2 <= U - L:
                        dense_intervals.append(open_interval(L, U))
    if not dense_intervals:
        return None
    else:
        # sort the dense intervals and return the union
        d = union_of_coho_intervals_minus_union_of_coho_intervals([[interval] for interval in dense_intervals], [])
        return d

class DirectedMoveCompositionCompletion:

    def __init__(self, fdms, covered_components=[], proj_add_vert=set(), show_plots=False, plot_background=None, function_at_border=None):
        self.show_plots = show_plots
        self.plot_background = plot_background
        # To show colorful components at borders, need the function. Otherwise set it to default None
        self.function_at_border = function_at_border
        self.move_dict = dict()
        self.sym_init_moves = dict() # updated by self.add_backward_moves() in round 0.
        self.covered_components = covered_components
        if covered_components:
            self.any_change_components = True
        else:
            self.any_change_components = False
        self.any_change_moves = set()
        for fdm in fdms:
            self.add_move(fdm)
        self.num_rounds = -1
        self.is_complete = False
        self.proj_add_vert = proj_add_vert

    def add_move(self, fdm):
        """
        Add a functional_directed_move to self.
        Merge or restrict functional directed move if necessary.
        """
        reduced_fdm = fdm.restricting(self.covered_components)
        if not reduced_fdm.intervals():
            return
        directed_move = fdm.directed_move
        if directed_move in self.move_dict:
            merged = merge_functional_directed_moves(self.move_dict[directed_move], reduced_fdm, show_plots=False)
            if merged.intervals() != self.move_dict[directed_move].intervals():
                # Cannot compare the functions themselves because of the "hash" magic of FastPiecewise.
                self.move_dict[directed_move] = merged
                self.any_change_moves.add(directed_move)
        else:
            self.move_dict[directed_move] = reduced_fdm
            self.any_change_moves.add(directed_move)

    def generate_zero_perturbation_points(self):
        zero_perturbation_points = copy(self.proj_add_vert)
        for fdm in self.move_dict.values():
            moved_projections = [fdm(x) for x in self.proj_add_vert if fdm.can_apply(x)]
            zero_perturbation_points.update(moved_projections)
            # sign condiction point
            if fdm.sign() == -1:
                invariant_point = fdm[1] / 2
                if fdm.can_apply(invariant_point):
                    assert fdm(invariant_point) == invariant_point
                    zero_perturbation_points.add(invariant_point)
        #  return a set, unsorted
        return zero_perturbation_points

    def plot(self, *args, **kwds):
        g = plot_directed_moves(list(self.move_dict.values()), **kwds)
        zero_perturbation =  zero_perturbation_partial_function(self.covered_components, \
                                                                self.generate_zero_perturbation_points())
        if zero_perturbation:
            g += plot_function_at_borders(zero_perturbation, color='magenta', legend_label='fixed perturbation (mod interpol)', thickness=3)
        if self.plot_background:
            g += self.plot_background
        if self.function_at_border and self.covered_components:
            g +=  plot_covered_components_at_borders(self.function_at_border, self.covered_components, **kwds)
        return g

    def maybe_show_plot(self, current_dense_move_plot=None):
        if (self.any_change_components or self.any_change_moves) and self.show_plots:
            logging.info("Plotting...")
            if self.is_complete:
                tag = 'completion'
                title = "Completion of moves" 
            elif self.num_rounds == -1:
                tag = 'completion-initial'
                title = "Initial moves"
            else:
                tag = 'completion-%s' % self.num_rounds
                title = "Moves after %s completion round%s" % (self.num_rounds, "s" if self.num_rounds > 1 else "")
            g = self.plot(legend_label='moves')
            if current_dense_move_plot:
                g += current_dense_move_plot
            show_plot(g, self.show_plots, tag, legend_title=title, legend_loc="upper left", object=self)
            logging.info("Plotting... done")

    def extend_components_by_moves(self):
        """
        Compose covered_components with fdms of self.
        """
        # try extending each component by applying all the moves
        new_components = []
        self.any_change_components = False
        for component in self.covered_components:
            extended_component = component
            for dm in self.move_dict.keys():
                if self.sym_init_moves.has_key(dm) and list(intersection_of_coho_intervals([extended_component, self.move_dict[dm].intervals()])):
                    #can extend the compnent by self.sym_init_moves[dm]
                    self.any_change_components = True
                    fdm = self.sym_init_moves[dm]
                    overlapped_ints = intersection_of_coho_intervals([extended_component, fdm.intervals()]) # generator
                    moved_intervals = [[fdm.apply_to_coho_interval(overlapped_int)] for overlapped_int in overlapped_ints ]
                    extended_component = union_of_coho_intervals_minus_union_of_coho_intervals([extended_component] + moved_intervals, [])
            new_components.append(extended_component)
        # got a extended component list, now merge components
        self.covered_components = reduce_covered_components(new_components)
        # kill moves
        self.reduce_moves_by_components()

    def extend_components_by_strip_lemma(self):
        move_keys = self.any_change_moves
        current_dense_move = []
        for dm_a in move_keys:
            for dm_b in move_keys:
                a = self.move_dict.get(dm_a, None)
                b = self.move_dict.get(dm_b, None)
                if not a or not b:
                    continue                        # move has been killed
                # check if the strip lemma gives new dense intervals
                d = check_for_strip_lemma(a, b)
                # d is not dominated by covered_components because we start with reduced moves.
                if d:
                    self.any_change_components = True
                    current_dense_move += d
                    logging.info("New dense move from strip lemma: %s" % d)
                    # merge components with d
                    dense_intervals, self.covered_components = merge_components_with_given_component(d, self.covered_components)
                    self.covered_components.append(dense_intervals)
                    # kill moves
                    self.reduce_moves_by_components(given_components = [dense_intervals])
        if current_dense_move:
            return plot_dense_moves(current_dense_move)
        else:
            return None

    def extend_moves_by_composition_of_moves(self):
        move_keys = self.any_change_moves
        self.any_change_moves = set()
        current_dense_move = []
        for dm_a in move_keys:
            a = self.move_dict.get(dm_a, None)
            if not a:
                continue                            # move has been killed
            # compose move a with b from sym_init_moves
            for (dm_b, b) in self.sym_init_moves.items():
                if not self.move_dict.has_key(dm_b):
                    continue                        # move has been killed
                c = compose_functional_directed_moves(a, b)
                if c:
                    self.add_move(c)

    def reduce_moves_by_components(self, given_components=None):
        """
        Reduce moves with covered_components of self.
        """
        new_move_dict = dict()
        for (dm, fdm) in self.move_dict.items():
            if not given_components is None:
                reduced_fdm = fdm.restricting(given_components)
            else:
                reduced_fdm = fdm.restricting(self.covered_components)
            if reduced_fdm.intervals():
                new_move_dict[dm] = reduced_fdm
        self.move_dict = new_move_dict

    def complete_one_round(self):
        logging.info("Completing %d functional directed moves and %d covered components..." % (len(self.move_dict), len(self.covered_components)))
        if self.any_change_components:
            self.extend_components_by_moves()
        current_dense_move_plot = self.extend_components_by_strip_lemma()
        self.extend_moves_by_composition_of_moves()
        self.num_rounds += 1
        self.maybe_show_plot(current_dense_move_plot)

    def add_backward_moves(self):
        self.num_rounds = 0
        move_keys = self.any_change_moves
        self.any_change_moves = set()
        for dm in move_keys:
            if dm[0] == 1:
                forward_fdm = self.move_dict[dm]
                backward_fdm = FunctionalDirectedMove(forward_fdm.range_intervals(), (1, -dm[1]))
                self.add_move(backward_fdm)
        if self.any_change_moves: # has new backward moves
            self.maybe_show_plot()
        self.any_change_moves.update(move_keys)
        self.sym_init_moves = copy(self.move_dict)

    def complete(self, max_num_rounds=None, error_if_max_num_rounds_exceeded=True):
        if self.num_rounds == -1:
            if self.any_change_moves:
                # do not show move diagram if there is no moves.
                self.maybe_show_plot()
            self.add_backward_moves()

        while (self.any_change_components or self.any_change_moves) and (max_num_rounds is None or self.num_rounds < max_num_rounds):
            self.complete_one_round()
        if max_num_rounds is not None and self.num_rounds == max_num_rounds:
            if error_if_max_num_rounds_exceeded:
                raise MaximumNumberOfIterationsReached, "Reached %d rounds of the completion procedure, found %d directed moves and %d covered components, stopping." % (self.num_rounds, len(self.move_dict), len(self.covered_components))
            else:
                logging.info("Reached %d rounds of the completion procedure, found %d directed moves and %d covered components, stopping." % (self.num_rounds, len(self.move_dict), len(self.covered_components)))
        else:
            self.is_complete = True
            logging.info("Completion finished.  Found %d directed moves and %d covered components."
                         % (len(self.move_dict), len(self.covered_components)))


    def results(self):
        return self.move_dict.values(), self.covered_components


def directed_move_composition_completion(fdms, covered_components=[], proj_add_vert=set(), show_plots=False, plot_background=None, function_at_border=None, max_num_rounds=None, error_if_max_num_rounds_exceeded=True):
    completion = DirectedMoveCompositionCompletion(fdms, covered_components=covered_components, \
                                                   proj_add_vert = proj_add_vert, \
                                                   show_plots=show_plots, \
                                                   plot_background=plot_background, \
                                                   function_at_border=function_at_border)
    completion.complete(max_num_rounds=max_num_rounds, error_if_max_num_rounds_exceeded=error_if_max_num_rounds_exceeded)
    return completion.results()

def plot_completion_diagram_background(fn):
    plot_background = plot_function_at_borders(fn, color='black', **ticks_keywords(fn, y_ticks_for_breakpoints=True))
    plot_background += polygon2d([[0,0], [0,1], [1,1], [1,0]], fill=False, color='grey')
    return plot_background

def generate_directly_covered_components(fn):
    if hasattr(fn, '_directly_covered_components'):
        return fn._directly_covered_components
    covered_components = []
    for face in generate_maximal_additive_faces(fn):
        if face.is_2D():
            (I, J, K) = face.minimal_triple
            K_mod_1 = interval_mod_1(K)
            component = union_of_coho_intervals_minus_union_of_coho_intervals([[open_interval(* I)], [open_interval(* J)], [open_interval(* K_mod_1)]],[])
            covered_components.append(component)
    return reduce_covered_components(covered_components)

# alias
generate_directly_covered_intervals = generate_directly_covered_components

@cached_function
def generate_directed_move_composition_completion(fn, show_plots=False, max_num_rounds=None, error_if_max_num_rounds_exceeded=True):
    completion = getattr(fn, "_completion", None)
    if completion is None:
        functional_directed_moves = generate_functional_directed_moves(fn)
        covered_components = generate_directly_covered_components(fn)
        proj_add_vert = projections_of_additive_vertices(fn)
        if show_plots:
            plot_background = plot_completion_diagram_background(fn)
        else:
            plot_background = None
        completion = fn._completion = DirectedMoveCompositionCompletion(functional_directed_moves,
                                                                        covered_components = covered_components,
                                                                        proj_add_vert = proj_add_vert,
                                                                        show_plots=show_plots,
                                                                        plot_background=plot_background,
                                                                        function_at_border=fn)
        # To show colorful components at borders, need the function_at_border. Otherwise set it to default None
        completion.complete(max_num_rounds=max_num_rounds, error_if_max_num_rounds_exceeded=error_if_max_num_rounds_exceeded)
    return completion.results()

def plot_walk_in_completion_diagram(seed, walk_dict):
    g = line([(seed,0), (seed,1)], color="limegreen", legend_label="seed value", linestyle=':')
    kwds = { 'legend_label': "reachable orbit" }
    for x in walk_dict.keys():
        g += line([(0, x), (seed, x)], color="limegreen", linestyle=':', **kwds)
        delete_one_time_plot_kwds(kwds)
    return g

def scan_domains_of_moves(functional_directed_moves):
     scans = [ scan_coho_interval_list(fdm.intervals(), fdm) for fdm in functional_directed_moves ]
     return merge(*scans)

def projections_of_additive_vertices(function):
    proj_add_vert = set()
    for (x, y, z, xeps, yeps, zeps) in generate_additive_vertices(function, reduced=True):
        proj_add_vert.update([x, y, fractional(z)])
    # returns a set, not sorted
    return proj_add_vert

def scan_zero_perturbation_point(x):
    dummy_zero_fdm = FunctionalDirectedMove([[0, 1]], (1,0))
    yield ((x, 0), 0, dummy_zero_fdm)
    yield ((x, 1), 0, dummy_zero_fdm)

def scan_zero_perturbation_points(zero_perturbation_points):
    scans = [ scan_zero_perturbation_point(x) for x in zero_perturbation_points ]
    return  merge(*scans)

def find_decomposition_into_intervals_with_same_moves(functional_directed_moves, zero_perturbation_points=set()):
    scan_of_zero_perturbation_points = scan_zero_perturbation_points(zero_perturbation_points)
    scan = merge(scan_domains_of_moves(functional_directed_moves), \
                 scan_of_zero_perturbation_points)
    moves = set()
    (on_x, on_epsilon) = (None, None)
    for ((x, epsilon), delta, move) in scan:
        if on_x and on_x < x:
            if moves:
                int = closed_or_open_or_halfopen_interval(on_x, x,
                                                          on_epsilon == 0, epsilon > 0)
                yield (int, list(moves))
        (on_x, on_epsilon) = (x, epsilon)
        if delta == -1:                         # beginning of interval
            assert move not in moves
            moves.add(move)
        elif delta == +1:                      # end of interval
            moves.remove(move)
        elif delta == 0:                       # a zero perturbation point (including sign contradiction point)
            pass
        else:
            raise ValueError, "Bad scan item"

def find_decomposition_into_stability_intervals_with_completion(fn, show_plots=False, max_num_it=None):
    if hasattr(fn, '_stability_orbits'):
        return
    fn._stability_orbits = []
    fdms, covered_components= generate_directed_move_composition_completion(fn, show_plots=show_plots)

    zero_perturbation_points = fn._completion.generate_zero_perturbation_points()
    decomposition = find_decomposition_into_intervals_with_same_moves(fdms, zero_perturbation_points )
    done_intervals = set()

    for (interval, moves) in decomposition:
        if interval not in done_intervals:
            #print interval
            orbit = set()
            walk_dict = dict()
            seed = (interval.a + interval.b) / 2
            for move in moves:
                moved_interval = move.apply_to_coho_interval(interval)
                #print "Applying move %s to %s gives %s." % (move, interval, moved_interval)

                moved_seed = move(seed)
                walk_sign= move.sign()
                done_intervals.add(moved_interval)
                orbit.add(moved_interval)
                walk_dict[moved_seed] = [walk_sign, None, None] 
            stability_orbit = (list(orbit), walk_dict, None)
            fn._stability_orbits.append(stability_orbit)
    logging.info("Total: %s stability orbits, lengths: %s" % (len(fn._stability_orbits), \
                    [ ("%s+" if to_do else "%s") % len(shifted_stability_intervals) \
                      for (shifted_stability_intervals, walk_dict, to_do) in fn._stability_orbits ]))

def zero_perturbation_partial_function(components, zero_perturbation_points):
    """
    Compute the partial function for which the perturbation, modulo
    perturbations that are interpolations of values at breakpoints, is
    known to be zero.
    """
    zero_function = FastLinearFunction(0, 0)
    pieces = []
    if zero_perturbation_points:
        pieces += [ (singleton_interval(x), zero_function) for x in zero_perturbation_points]
    if components:
        pieces += [ (interval, zero_function) for component in components for interval in component ]
    if pieces:
        return FastPiecewise(pieces)
    else:
        return None

def generate_uncovered_components(fn, show_plots=False):
    if not hasattr(fn, '_stability_orbits'):
        find_decomposition_into_stability_intervals_with_completion(fn, show_plots=show_plots)
    uncovered_components = [ sorted(orbit, key=coho_interval_left_endpoint_with_epsilon) \
                             for (orbit, _, _) in fn._stability_orbits]
    return uncovered_components

def stab_int_length(x):
    (orbit, walk_dict, _) = x
    int = orbit[0]
    return interval_length(int)

def generate_generic_seeds_with_completion(fn, show_plots=False, max_num_it=None):
    # Ugly compatibility interface.
    find_decomposition_into_stability_intervals_with_completion(fn, show_plots=show_plots)
    for (orbit, walk_dict, _) in sorted(fn._stability_orbits, key=stab_int_length, reverse=True):
        int = orbit[0]
        if interval_length(int) > 0:
            seed = (int.a + int.b) / 2
            stab_int = closed_or_open_or_halfopen_interval(int.a - seed, int.b - seed, int.left_closed, int.right_closed)
            yield (seed, stab_int, walk_dict)

def stuff_with_random_irrational_function():
    while True:
        del1 = randint(1, 100) / 1000
        del2 = randint(1, 60) / 1000
        print "del1 = %s, del2 = %s" % (del1, del2)
        try:
            h = the_irrational_function_t1_t2(del1=del1, del2=del2)
            break
        except ValueError:
            print "... parameters do not describe a function, retrying."
    dmoves = generate_functional_directed_moves(h)
    fdms, _ = directed_move_composition_completion(dmoves, max_num_rounds=None)
    plot_directed_moves(fdms).show(figsize=40)
    
def merit_index(fn):
    r"""
    Compute Gomory--Johnson's merit index.

    It is defined as twice the area of additive points in the unit square
    (see [RD, Definition 19.28], which refers to [61]).

    EXAMPLES:

    The merit index for GMIC is 2f^2 - 2f + 1 [RD]::

        sage: def merit_index_gmic(f):
        ....:     return 2*f^2 - 2*f + 1
        sage: merit_index(gmic(1/4)) == merit_index_gmic(1/4)
        True
        sage: merit_index(gmic(1/2)) == merit_index_gmic(1/2)
        True

    The merit index for `drlm_2_slope_limit_1_1` is 4 f^2 - 6 f + 3 [40, p. 164]::

        sage: def merit_index_drlm_2_slope_limit_1_1(f):
        ....:     return 4 * f^2 - 6 * f + 3
        sage: merit_index(drlm_2_slope_limit_1_1(3/5)) == merit_index_drlm_2_slope_limit_1_1(3/5)
        True
        sage: merit_index(drlm_2_slope_limit_1_1(1/2)) == merit_index_drlm_2_slope_limit_1_1(1/2)
        True

    Reference:

     - [40] S. S. Dey, J.-P. P. Richard, Y. Li, and L. A. Miller, Extreme
       inequalities for infinite group problems, Mathematical Programming 121
       (2010), 145-170.

     - [61] R. E. Gomory and E. L. Johnson, T-space and cutting planes,
       Mathematical Programming 96 (2003), 341-375.

     - [RD] J.-P. P. Richard and S. S. Dey, The group-theoretic approach in mixed
       integer programming, 50 Years of Integer Programming 1958-2008
       (M. Juenger, T. M. Liebling, D. Naddef, G. L. Nemhauser, W. R. Pulleyblank,
       G. Reinelt, G. Rinaldi, and L. A. Wolsey, eds.), Springer Berlin Heidelberg,
       2010, pp. 727-801, doi:10.1007/978-3-540-68279-0_19, ISBN 978-3-540-68274-5.
    """
    return 2 * add(Polyhedron(face.vertices).volume()
                   for face in generate_maximal_additive_faces(fn)
                   if face.is_2D())

def arithmetic_complexity(fn, f=None, q=None):
    r"""
    Compute the arithmetic complexity

    It is defined as the least common denominator of the values fn(i/q) for i=0,1,...,q,
    where fn is a piecewise linear function with rational breakpoints in (1/q)Z,
    or a discrete function with its domain contained in (1/q)Z.

    EXAMPLES::

        sage: logging.disable(logging.INFO)
        sage: h = gmic(7/9)
        sage: arithmetic_complexity(h)
        14
        sage: arithmetic_complexity(restrict_to_finite_group(h))
        14
        sage: arithmetic_complexity(dg_2_step_mir_limit())
        24

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    if q is None:
        q = finite_group_order_from_function_f_oversampling_order(fn, f=None, oversampling=None, order=None)
    if fn.is_discrete():
        v = fn.values_at_end_points()
    elif fn.is_continuous():
        v = [fn(x/q) for x in range(q+1)]
    else:
        v = []
        for x in range(q+1):
            for epsilon in [0,1,-1]:
                try:
                    y = fn.limit(x/q, epsilon)
                    v.append(y)
                except:
                    pass
    is_rational_v, v = is_all_QQ(v)
    if is_rational_v:
        v_denominator = [denominator(y) for y in v]
        return lcm(v_denominator)
    else:
        raise ValueError, "This is a function has irrational value on (1/%s)Z, so the arithmetic_complexity is not available." % q
