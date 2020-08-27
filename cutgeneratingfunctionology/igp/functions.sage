"""
    sage: import warnings
    sage: warnings.filterwarnings('ignore', 'Matplotlib is building the font cache using fc-list. This may take a moment.')
    sage: from cutgeneratingfunctionology.igp import *
"""

import itertools
from six.moves import map
from six.moves import range
from six.moves import zip
from functools import reduce
from six.moves import input
from itertools import chain
from bisect import bisect_left


def unique_list(iterator):
    r"""
    Returns the list of the elements in the iterator without repetition.
    """
    l = []
    s = set()
    for i in iterator:
        if i not in s:
            s.add(i)
            l.append(i)
    return l

def fractional(num):
    r"""
    Reduces a number modulo `1`.
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
    r"""
    Computes the slack in subaddivity. See also ``delta_pi_general`` for discontinuous case.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: fn = not_minimal_2()
        sage: delta_pi(fn, 1/5, 3/5)
        0
    """
    return fn(fractional(x))+fn(fractional(y))-fn(fractional(x+y))

plot_2d_complex_discontinuous_kwds = plot_2d_complex_continuous_kwds = {'color': 'grey'}

def plot_2d_complex(function, continuous_plot_kwds=None, discontinuous_plot_kwds=None):
    r"""
    Returns a plot of the horizontal lines, vertical lines, and diagonal lines of the complex.
    """
    bkpt = function.end_points()
    x = var('x')
    p = Graphics()
    kwd = ticks_keywords(function, True)
    kwd['legend_label'] = "Complex Delta P"
    
    if continuous_plot_kwds is None:
        continuous_plot_kwds = plot_2d_complex_continuous_kwds
    if discontinuous_plot_kwds is None:
        discontinuous_plot_kwds = plot_2d_complex_discontinuous_kwds
    continuous_plot_kwds.update(kwd)
    discontinuous_plot_kwds.update(kwd)
    plot_kwds_hook(continuous_plot_kwds)
    plot_kwds_hook(discontinuous_plot_kwds)

    def plot_kwds_at(x):
        l = function.limits(x)
        left_continuous = l[-1] is None or l[-1] == l[0]
        right_continuous = l[1] is None or l[1] == l[0]
        if left_continuous and right_continuous:
            return continuous_plot_kwds
        else:
            return discontinuous_plot_kwds
    ## We now use lambda functions instead of Sage symbolics for plotting, 
    ## as those give strange errors when combined with our RealNumberFieldElement.
    for i in range(1,len(bkpt)):
        #p += plot(lambda x: bkpt[i]-x, (x, 0, bkpt[i]), color='grey', **kwd)
        p += line([(0,  bkpt[i]), (bkpt[i], 0)], **plot_kwds_at(bkpt[i]))
        delete_one_time_plot_kwds(continuous_plot_kwds)
        delete_one_time_plot_kwds(discontinuous_plot_kwds)
    for i in range(1,len(bkpt)-1):
        #p += plot(lambda x: (1+bkpt[i]-x), (x, bkpt[i], 1), color='grey')
        p += line([(bkpt[i], 1), (1, bkpt[i])], **plot_kwds_at(bkpt[i]))
        delete_one_time_plot_kwds(continuous_plot_kwds)
        delete_one_time_plot_kwds(discontinuous_plot_kwds)
    for i in range(len(bkpt)):
        p += plot(bkpt[i], (0, 1), **plot_kwds_at(bkpt[i]))
        delete_one_time_plot_kwds(continuous_plot_kwds)
        delete_one_time_plot_kwds(discontinuous_plot_kwds)
    y=var('y')
    for i in range(len(bkpt)):
        p += parametric_plot((bkpt[i],y), (y,0,1), **plot_kwds_at(bkpt[i]))
        delete_one_time_plot_kwds(continuous_plot_kwds)
        delete_one_time_plot_kwds(discontinuous_plot_kwds)
    return p

##
##
##

def projection(vertices,linear_form):
    r"""
    Computes the projection of vertices based on the linear form.

    vertices is a list of vertices (2-tuples).
    
    linear_form is a 2-element list:
    
    - Projection on `x`: `[1,0]`

    - Projection on `y`: `[0,1]`

    - Projection on `x + y`: `[1,1]`
    """
    temp = []
    for i in vertices:
        temp.append(i[0]*linear_form[0]+i[1]*linear_form[1])
    if max(temp) == min(temp):
        return [min(temp)]
    else:
        return [min(temp), max(temp)]

def projections(vertices):
    r"""
    Computes `F(I,J,K)`.
    """
    return [projection(vertices, [1,0]),projection(vertices, [0,1]),projection(vertices, [1,1])]    

def verts(I1, J1, K1):
    r"""
    Computes the vertices based on `I, J` and `K`.
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

equiv7_mode = False

## FIXME: Document what it does when input is not subadditive.
def generate_maximal_additive_faces(fn):
    r"""Compute the list of maximal additive faces of the subadditive function ``fn``.

    In the discontinuous case, this includes the faces determined by
    limit-additivities as well: Additive faces are defined in the
    discontinuous case according to
    :cite:`hildebrand-koeppe-zhou:algo-paper`, Definition 9.7
    ({def:set-of-additivities-and-limit-additivities-of-a-face-F}),
    following Lemma 9.6 (lemma:additive-face-discontinuous): The face `E
    \in \Delta\mathcal P` is additive if there exists an enclosing face
    `F \in \Delta\mathcal P` such that for some (equivalently, all) `(x,
    y) \in \mathop{\mathrm{rel int}}(E)`, the limit `\Delta``fn``_F(x,
    y) = 0``.

    The additive faces (together with the empty face) form a subcomplex
    of `Delta\mathcal P`.  Thus, they are determined by the maximal
    elements.

    (However, in the current implementation, in the discontinuous case,
    this function actually computes all additive faces, not just the
    maximal ones.)

    Note that `\Delta` ``fn`` is invariant under swapping `x` and `y`.
    This function does not break this symmetry.
"""
    if hasattr(fn, '_maximal_additive_faces'):
        return fn._maximal_additive_faces
    if fn.is_discrete():
        result = generate_additive_faces_discrete(fn)   # all 0-dimensional, so maximal
    elif fn.is_continuous():
        result = generate_maximal_additive_faces_continuous(fn)
    else:
        additive_faces = generate_additive_faces_general(fn)
        additive_faces_2D = [F for F in additive_faces if F.is_2D()]
        maximal_additive_face_1D = [F1 for F1 in additive_faces if F1.is_1D() and all(not set(F1.vertices).issubset(set(F2.vertices)) for F2 in additive_faces_2D)]
        result = additive_faces_2D + maximal_additive_face_1D
        for F0 in additive_faces:
            if F0.is_0D():
                v = F0.vertices[0]
                if all(not v in set(F.vertices) for F in additive_faces_2D) and all(not v in set(F.vertices) for F in maximal_additive_face_1D):
                    result.append(F0)
        ####result = generate_additive_faces_general(fn)
        ##### FIXME: Filter out non-maximal faces
    fn._maximal_additive_faces = result
    return result

def generate_additive_faces(fn):
    r"""
    Return the list of all additive faces of the subadditive function ``fn``.

    In the discontinuous case, this includes the faces determined by
    limit-additivities as well: Additive faces are defined in the
    discontinuous case according to
    :cite:`hildebrand-koeppe-zhou:algo-paper`, Definition 9.7
    ({def:set-of-additivities-and-limit-additivities-of-a-face-F}),
    following Lemma 9.6 (lemma:additive-face-discontinuous): The face `E
    \in \Delta\mathcal P` is additive if there exists an enclosing face
    `F \in \Delta\mathcal P` such that for some (equivalently, all) `(x,
    y) \in \mathop{\mathrm{rel int}}(E)`, the limit `\Delta``fn``_F(x,
    y) = 0``.

    The additive faces (together with the empty face) form a subcomplex
    of `Delta\mathcal P`.  Thus, they are determined by the maximal
    elements; these are provided by ``generate_maximal_additive_faces``.

    Note that `\Delta` ``fn`` is invariant under swapping `x` and `y`.
    This function does not break this symmetry.

    See also: ``generate_additive_faces_sans_limits``.
    """
    if hasattr(fn, '_additive_faces'):
        return fn._additive_faces
    if fn.is_discrete():
        result = generate_additive_faces_discrete(fn)
    else:
        result = generate_additive_faces_general(fn)
    fn._additive_faces = result
    return result

def is_additive_face_sans_limits(fn, face):
    r"""Test whether ``face`` is additive-sans-limits for ``fn``.

    See ``generate_additive_faces_sans_limits``.
    """
    ver = face.vertices
    n = len(ver)
    mx, my = sum([x for (x,y) in ver])/n, sum([y for (x,y) in ver])/n
    return delta_pi(fn, mx, my) == 0

def generate_additive_faces_sans_limits(fn):
    r"""Return a list of all additive faces of the subadditive function ``fn`` without taking limits into consideration.

    A face ``E`` is additive-sans-limits if for some (equivalently, all)
    `(x, y) \in \mathop{\mathrm{rel int}}(E)`, we have `\Delta``fn``(x,
    y) = 0``.

    Such a face ``E`` is the convex hull of vertices `(x, y)` that have
    the limits `\Delta``fn``_E(x, y) = 0`.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = kzh_minimal_has_only_crazy_perturbation_1()
        sage: x = h.end_points()
        sage: F = Face(([x[6]],[1+x[4]-x[6]],[1+x[4]]))
        sage: F in generate_additive_faces(h)
        True
        sage: F in generate_maximal_additive_faces(h)
        False
        sage: F in generate_additive_faces_sans_limits(h)
        True
    """
    if not hasattr(fn, '_additive_faces_sans_limits'):
        fn._additive_faces_sans_limits = [ face for face in generate_additive_faces(fn)
                                           if is_additive_face_sans_limits(fn, face) ]
    return fn._additive_faces_sans_limits

additive_fill_color = additive_color = "mediumspringgreen"

import six

### Create a new class representing a "face" (which knows its
### vertices, minimal triple, whether it's a translation/reflection,
### etc.; whether it's solid or dense).
class Face:
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: F = Face([[1/5, 3/10], [3/4, 17/20], [1, 6/5]])
        sage: F.vertices
        [(1/5, 17/20), (3/10, 3/4), (3/10, 17/20), (1/5, 4/5), (1/4, 3/4)]
        sage: F.minimal_triple
        ([1/5, 3/10], [3/4, 17/20], [1, 23/20])
    """
    def __init__(self, triple, vertices=None, is_known_to_be_minimal=False):
        if not vertices:
            vertices = verts(triple[0], triple[1], triple[2])
            if not vertices:
                raise NotImplementedError("An empty face. This could mean need to shift triple[2] by (1,1). Not implemented.")
        self.vertices = vertices
        i, j, k = projections(vertices)
        self.minimal_triple = minimal_triple = (i, j, k)
        #self._warned_about_non_minimal_triple = False
        #if is_known_to_be_minimal and not triples_equal(minimal_triple, triple) and not self._warned_about_non_minimal_triple:
        #    logging.warning("Provided triple was not minimal: %s reduces to %s" % (triple, minimal_triple))
        #    self._warned_about_non_minimal_triple = True
            # FIXME: Check why (i,j,k) != (i,j,k+1) can happen.

    def __repr__(self):
        return '<Face ' + repr(self.minimal_triple) + '>'

    def __contains__(self, point):
        xyz = (point[0], point[1], point[0] + point[1])
        return all(element_of_int(xyz[i], self.minimal_triple[i]) for i in range(3))

    def interior_contains(self, point):
        if not self.is_2D():
            return False
        xyz = (point[0], point[1], point[0] + point[1])
        return all(self.minimal_triple[i][0] < xyz[i] < self.minimal_triple[i][1] for i in range(3))

    def plot(self, rgbcolor=(0.0 / 255.0, 250.0 / 255.0, 154.0 / 255.0), fill_color=None, edge_thickness=2, vertex_size=30, *args, **kwds):
        y = var('y')
        trip = self.minimal_triple
        vert = self.vertices
        if fill_color is None:
            fill_color = additive_fill_color
        if self.is_0D():
            return point((trip[0][0], \
                          trip[1][0]), rgbcolor = rgbcolor, size=vertex_size, **kwds)
        elif self.is_horizontal():
            return line([(trip[0][0],trip[1][0]),(trip[0][1],trip[1][0])], rgbcolor = rgbcolor, thickness=edge_thickness, **kwds)
        elif self.is_vertical():
            return line([(trip[0][0],trip[1][0]),(trip[0][0],trip[1][1])], rgbcolor = rgbcolor, thickness=edge_thickness, **kwds)
        elif self.is_diagonal():
            return line([(trip[0][0],trip[2][0]-trip[0][0]),(trip[0][1],trip[2][0]-trip[0][1])], rgbcolor = rgbcolor, thickness=edge_thickness, **kwds)
        elif self.is_2D():
            ## Sorting is necessary for this example:
            ## plot_2d_diagram(lift(piecewise_function_from_robert_txt_file("data/dey-richard-not-extreme.txt"))
            return polygon(convex_vert_list(vert), color=fill_color, **kwds)

    def polyhedron(self):
        return Polyhedron(self.vertices)

    def is_directed_move(self):
        return self.is_1D()
        
    def directed_move_with_domain_and_codomain(self):
        r"""
        Maps a horizontal or vertical edge to a forward translation move.  Maps a diagonal edge to a reflection move.

        .. NOTE::

            Backward translation moves will be added to ``DirectedMoveCompositionCompletion`` by ``add_backward_moves`` in round 0.

            The domain and codomain are lists of open intervals. Their endpoints will be considered when treating the additive vertices.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: from cutgeneratingfunctionology.igp import *
            sage: face_hor = Face([[2/5, 3/5],[4/5],[6/5,7/5]])
            sage: face_hor.directed_move_with_domain_and_codomain()
            ((1, -1/5), [<Int(2/5, 3/5)>], [<Int(1/5, 2/5)>])
            sage: face_ver = Face([[4/5],[2/5, 3/5],[6/5,7/5]])
            sage: face_ver.directed_move_with_domain_and_codomain() == face_hor.directed_move_with_domain_and_codomain()
            True
            sage: face_dia = Face([[2/5, 3/5],[2/5, 1/5],[4/5]])
            sage: face_dia.directed_move_with_domain_and_codomain()
            ((-1, 4/5),
             [<Int(2/5, 3/5)>, <Int(1/5, 2/5)>],
             [<Int(1/5, 2/5)>, <Int(2/5, 3/5)>])
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
            raise ValueError("Face does not correspond to a directed move: %s" % self)

    def functional_directed_move(self, is_backward_translation=False):
        r"""
        Returns the (forward by default or backward if ``is_backward_translation=True``) translation move if given a horizontal/vertical edge, or the reflection move if given a diagonal edge.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: from cutgeneratingfunctionology.igp import *
            sage: face_hor = Face([[2/5, 3/5],[4/5],[6/5,7/5]])
            sage: face_hor.functional_directed_move()
            <FunctionalDirectedMove (1, -1/5) with domain [<Int(2/5, 3/5)>], range [<Int(1/5, 2/5)>]>
            sage: face_hor.functional_directed_move(is_backward_translation=True)
            <FunctionalDirectedMove (1, 1/5) with domain [<Int(1/5, 2/5)>], range [<Int(2/5, 3/5)>]>
            sage: face_ver = Face([[4/5],[2/5, 3/5],[6/5,7/5]])
            sage: face_ver.functional_directed_move() == face_hor.functional_directed_move()
            True
            sage: face_ver.functional_directed_move(is_backward_translation=True) == face_hor.functional_directed_move(is_backward_translation=True)
            True
            sage: face_dia = Face([[2/5, 3/5],[2/5, 1/5],[4/5]])
            sage: face_dia.functional_directed_move()
            <FunctionalDirectedMove (-1, 4/5) with domain [<Int(1/5, 2/5)>, <Int(2/5, 3/5)>], range [<Int(2/5, 3/5)>, <Int(1/5, 2/5)>]>
            sage: face_dia.functional_directed_move(is_backward_translation=True)
            <FunctionalDirectedMove (-1, 4/5) with domain [<Int(1/5, 2/5)>, <Int(2/5, 3/5)>], range [<Int(2/5, 3/5)>, <Int(1/5, 2/5)>]>

        """
        directed_move, domain, codomain = self.directed_move_with_domain_and_codomain()
        if not is_backward_translation:
            fdm = FunctionalDirectedMove(domain, directed_move)
        else:
            fdm = FunctionalDirectedMove(codomain, (directed_move[0], -directed_move[0]*directed_move[1]))
        return fdm

    def covered_component(self):
        if self.is_2D():
            (I, J, K) = self.minimal_triple
            K_mod_1 = interval_mod_1(K)
            return union_of_coho_intervals_minus_union_of_coho_intervals([[open_interval(*I)], [open_interval(*J)], [open_interval(*K_mod_1)]],[])
        raise ValueError("Face does not give a covered component")

    def dimension(self):
        if self.is_0D():
            return 0
        elif self.is_1D():
            return 1
        elif self.is_2D():
            return 2
        else:
            return -1

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
    def __hash__(self):
        return sum([hash(x) for i in self.minimal_triple for x in i])

    if six.PY2:
        def __cmp__(left, right):
            return cmp(left.minimal_triple, right.minimal_triple)
    else:
        def __eq__(left, right):
            return left.minimal_triple == right.minimal_triple

def plot_faces(faces, **kwds):
    p = Graphics()
    for f in faces:
        if type(f) == list or type(f) == tuple: #legacy
            f = Face(f)
        p += f.plot(**kwds)
    return p

def plot_trivial_2d_diagram_with_grid(function, xgrid=None, ygrid=None): 
    r"""
    Returns a plot of the 2d complex with vertices marked that 
    have `\Delta \pi = 0`.  

    Does not use any complicated code.
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

import operator

from six.moves import reduce
    
def convex_vert_list(vertices):
    if len(vertices) <= 3:
        return vertices
    else:
        center = reduce(operator.add, list(map(vector, vertices))) / len(vertices)
        return sorted(vertices, key=lambda x: atan2(x[1] - center[1], x[0] - center[0]))

def plot_kwds_hook(kwds):
    pass

def plot_2d_diagram(fn, show_function=True, show_projections=True, known_minimal=False, f=None, colorful=False, additive_color=None, function_color="blue"):
    r"""
    Returns a plot of the 2d complex (`\\Delta P`) of fn with shaded
    additive faces, i.e., faces where `\\Delta \\pi = 0`.
    
    - If known_minimal is ``False`` (the default), highlight
    non-subadditive or non-symmetric vertices of the 2d complex.

    - If show_function is ``True`` (the default), plot the function at the left and top borders of the diagram via ``plot_function_at_borders``. 

    - If show_projections is ``True`` (the default), plot the projections `p_1(F), p_2(F), p_3(F)` of 
    all full-dimensional additive faces via ``plot_projections_at_borders``.

    To show only a part of the diagram, use::

        sage: from cutgeneratingfunctionology.igp import *
        sage: plot_2d_diagram(h).show(xmin=0.25, xmax=0.35, ymin=0.25, ymax=0.35)  # not tested

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = FastPiecewise([[closed_interval(0,1/4), FastLinearFunction(4, 0)],
        ....:                    [open_interval(1/4, 1), FastLinearFunction(4/3, -1/3)],
        ....:                    [singleton_interval(1), FastLinearFunction(0,0)]])
        sage: plot_2d_diagram(h)
        Graphics object ...
        sage: h = FastPiecewise([[closed_interval(0,1/4), FastLinearFunction(4, 0)],
        ....:                    [open_interval(1/4,1/2), FastLinearFunction(3, -3/4)],
        ....:                    [closed_interval(1/2, 3/4), FastLinearFunction(-2, 7/4)],
        ....:                    [open_interval(3/4,1), FastLinearFunction(3, -2)],
        ....:                    [singleton_interval(1), FastLinearFunction(0,0)]])
        sage: plot_2d_diagram(h)
        Graphics object ...
    """
    if f is None:
        f = find_f(fn, no_error_if_not_minimal_anyway=True)
    if additive_color is None:
        import cutgeneratingfunctionology
        additive_color = cutgeneratingfunctionology.igp.additive_color
        additive_fill_color = cutgeneratingfunctionology.igp.additive_fill_color
    else:
        additive_fill_color = additive_color
    faces = generate_maximal_additive_faces(fn)
    p = Graphics()
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
            p += face.plot(rgbcolor = additive_color, fill_color = additive_fill_color, **kwds)
        delete_one_time_plot_kwds(kwds)

    ### For non-subadditive functions, show the points where delta_pi is negative.
    if not known_minimal:
        nonsubadditive_vertices = generate_nonsubadditive_vertices(fn, reduced=False)
        kwds = { 'legend_label' : "Subadditivity violated" }
        plot_kwds_hook(kwds)
        if fn.is_continuous():
            nonsubadditive_vertices = unique_list((x,y) for (x, y, z, xeps, yeps, zeps) in nonsubadditive_vertices)
            p += point(nonsubadditive_vertices,
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
                nonsymmetric_vertices = unique_list((x,y) for (x, y, xeps, yeps) in nonsymmetric_vertices)
                p += point(nonsymmetric_vertices,
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
    p += plot_2d_complex(fn)
    if show_projections:
        p += plot_projections_at_borders(fn)
    if show_function:
        p += plot_function_at_borders(fn, covered_components = covered_components, color=function_color)
    return p

def plot_covered_components_at_borders(fn, covered_components=None, **kwds):
    r"""
    Colorful decoration.

    Plot ``fn`` on covered intervals with different colors according to slope values,
    on the upper and the left border of the 2d diagrams.
    """
    p = Graphics()
    if not covered_components is None:
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

def color_for_delta(deltafn, additive_color=None):
    if additive_color is None:
        import cutgeneratingfunctionology
        additive_color = cutgeneratingfunctionology.igp.additive_color
    if deltafn > 0:
        return "white"
    elif deltafn == 0:
        return additive_color
    else:
        return "red"

def plot_2d_diagram_with_eps_cones(fn, show_function=True, f=None, conesize=200,
                                   additive_color=None, function_color="blue"):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = zhou_two_sided_discontinuous_cannot_assume_any_continuity()
        sage: g = plot_2d_diagram_with_eps_cones(h)
        sage: h = not_minimal_2()
        sage: g = plot_2d_diagram_with_eps_cones(h)
    """
    if f is None:
        f = find_f(fn, no_error_if_not_minimal_anyway=True)
    g = plot_2d_complex(fn)
    if show_function:
        g += plot_function_at_borders(fn, color=function_color)
    bkpt = fn.end_points()
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]
    type_1_vertices = [(x, y, x+y) for x in bkpt for y in bkpt if x <= y]
    type_2_vertices = [(x, z-x, z) for x in bkpt for z in bkpt2 if x < z < 1+x]
    vertices = set(type_1_vertices + type_2_vertices)
    if fn.is_continuous():
        for (x, y, z) in vertices:
            deltafn = delta_pi(fn, x, y)
            color = color_for_delta(deltafn, additive_color=additive_color)
            g += point([(x, y), (y, x)], color=color, size=conesize, zorder=-1)
    else:
        for (x, y, z) in vertices:
            eps_list = [(0,0,0)]+list(nonzero_eps)
            delta_list = [ delta_pi_general(fn, x, y, xyz_eps) for xyz_eps in eps_list ]
            if all(delta != 0 for delta in delta_list):
                # Don't plot anything if there's no additivity or limit-additivity at a vertex
                continue
            for deltafn, (xeps, yeps, zeps) in zip(delta_list, eps_list):
                color = color_for_delta(deltafn, additive_color=additive_color)
                g += plot_limit_cone_of_vertex(x, y, epstriple_to_cone((xeps, yeps, zeps)), color=color, r=0.03)
                g += plot_limit_cone_of_vertex(y, x, epstriple_to_cone((yeps, xeps, zeps)), color=color, r=0.03)
    return g

def plot_2d_diagram_with_face_cones(fn, show_function=True, f=None, conesize=200,
                                    additive_color=None, function_color="blue"):
    r"""
    This shows larger cones, corresponding to enclosing faces rather than
    cones corresponding to epsilons.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = zhou_two_sided_discontinuous_cannot_assume_any_continuity()
        sage: g = plot_2d_diagram_with_face_cones(h)
        sage: h = not_minimal_2()
        sage: g = plot_2d_diagram_with_face_cones(h)
    """
    if f is None:
        f = find_f(fn, no_error_if_not_minimal_anyway=True)
    g = plot_2d_complex(fn)
    if show_function:
        g += plot_function_at_borders(fn, color=function_color)

    from cutgeneratingfunctionology.spam.real_set import RealSet   # param-safe version of RealSet
    domain = RealSet([0, 1])
    faces = set(Face(triple)
                for triple in generate_triples_with_projections_intersecting(fn, domain, break_symmetry=False, halfopen_cover_only=False))
    for F in faces:
        for v in F.vertices:
            deltafn = delta_pi_of_face(fn, v[0], v[1], F)
            color = color_for_delta(deltafn, additive_color=additive_color)
            g += plot_limit_cone_of_vertex_in_face(v[0], v[1], F, color)

    return g

def plot_2d_diagram_with_cones(fn, show_function=True, f=None, conesize=200,
                               additive_color=None, function_color="blue"):
    if fn.is_continuous():
        plot_it = plot_2d_diagram_with_eps_cones
    else:
        plot_it = plot_2d_diagram_with_face_cones
    return plot_it(fn, show_function=show_function, f=f, conesize=conesize,
                   additive_color=additive_color, function_color=function_color)

def plot_2d_diagram_additive_domain_sans_limits(fn, show_function=True, f=None, additive_color=additive_color, function_color='blue', **kwds):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = hildebrand_discont_3_slope_1()
        sage: g = plot_2d_diagram_additive_domain_sans_limits(h)
    """
    if f is None:
        f = find_f(fn, no_error_if_not_minimal_anyway=True)
    g = Graphics()
    for face in generate_additive_faces_sans_limits(fn):
        g += face.plot(rgbcolor=additive_color, fill_color=additive_color, **kwds)
    g += plot_2d_complex(fn)
    if show_function:
        g += plot_function_at_borders(fn, color=function_color)
    return g

plot_function_at_borders_kwds = {}

def plot_function_at_borders(fn, color='blue', legend_label="Function pi", covered_components=None, **kwds):
    r"""
    Plot the function twice, on the upper and the left border, 
    to decorate 2d diagrams.
    """
    p = Graphics()
    bkpt = fn.end_points()
    limits = fn.limits_at_end_points()
    if not covered_components is None:
        color = 'black'
    kwds.update(plot_function_at_borders_kwds)
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
    r"""
    Plot the projections `p_1(F), p_2(F), p_3(F)` of all full-dimensional
    additive faces `F` of fn as gray shadows:

    - `p_1(F)` at the top border,
    - `p_2(F)` at the left border, 
    - `p_3(F)` at the bottom and the right borders.
    """
    return plot_projections_of_faces(end_points=fn.end_points(),
                                    nonsubadditive_vertices=generate_nonsubadditive_vertices(fn),
                                    additive_faces=generate_maximal_additive_faces(fn))

def plot_projections_of_faces(additive_faces, nonsubadditive_vertices=(), end_points=()):
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
    for face in additive_faces:
        I, J, K = face.minimal_triple
        I_J_verts.update(I) # no need for J because the x-y swapped face will also be processed
        K_verts.update(K)
        g += plot_projections_of_one_face(face, IJK_kwds)
    for (x, y, z, xeps, yeps, zeps) in nonsubadditive_vertices:
        I_J_verts.add(x)
        I_J_verts.add(y)
        K_verts.add(z)
    # plot dashed help lines corresponding to non-breakpoint projections. 
    # (plot_2d_complex already draws solid lines for the breakpoints.)
    I_J_verts.difference_update(end_points)
    for x in I_J_verts:
        g += line([(x, 0), (x, 1)], linestyle=':', color='grey')
        g += line([(0, x), (1, x)], linestyle=':', color='grey')
    K_verts.difference_update(end_points)
    K_verts.difference_update(1 + x for x in end_points)
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
            raise ValueError("Bad face: %s" % face)
        delete_one_time_plot_kwds(IJK_kwds[2])
    return g
    
def interval_mod_1(interval):
    r"""
    Represent the given proper interval modulo `1` as a subinterval of `[0,1]`.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: from cutgeneratingfunctionology.igp import *
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
        raise ValueError("Not an interval: %s" % interval)

def generate_covered_components(function):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = hildebrand_discont_3_slope_1()
        sage: generate_covered_components(h)
        [[<Int(0, 1/8)>, <Int(3/8, 1/2)>],
         [<Int(1/8, 3/8)>, <Int(1/2, 5/8)>, <Int(7/8, 1)>],
         [<Int(5/8, 7/8)>]]

    Also covered components that come from the strip lemma
    in the irrational case are included::

        sage: h = bhk_irrational()
        sage: len(generate_covered_components(h))
        3
    """
    if hasattr(function,  "_completion"):
        return function._completion.covered_components
    bkpts = function.end_points()
    field = bkpts[0].parent()
    if isinstance(field, ParametricRealField):
        bkpts_are_rational = is_all_QQ_fastpath(result_concrete_value(field, bkpts))
    else:
        bkpts_are_rational = is_all_QQ_fastpath(bkpts)
    if not bkpts_are_rational:
        fdms, covered_components = generate_directed_move_composition_completion(function)
        return covered_components
    global strategical_covered_components
    if strategical_covered_components:
        covered_components = generate_covered_components_strategically(function)
        return covered_components
    # Rational case.  We run a special version of the code in
    # generate_directed_move_composition_completion here,
    # which only runs extend_components_by_moves and does not
    # supply proj_add_vert.  (Why?)
    #proj_add_vert = projections_of_additive_vertices(function)
    functional_directed_moves = generate_functional_directed_moves(function)
    covered_components = generate_directly_covered_components(function)
    completion = DirectedMoveCompositionCompletion(functional_directed_moves,
                                                   covered_components=covered_components,
                                                   #proj_add_vert=proj_add_vert,
                                                   function=function)
    completion.add_backward_moves()
    while completion.any_change_components:
        completion.extend_components_by_moves()
        completion.num_rounds += 1
    completion.is_complete = True
    #function._completion = completion
    logging.info("Completion finished.  Found %d covered components."
                 % len(completion.covered_components))
    return completion.covered_components

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

y_ticks_use_solidus = False

def latex_tick_formatter_x(x):
    return "$%s$" % latex(x)

def latex_tick_formatter_y(y):
    if y_ticks_use_solidus and y in QQ and y not in ZZ:
        y = QQ(y)
        return "$%s/%s$" % (y.numerator(), y.denominator())
    else:
        return "$%s$" % latex(y)

def default_igp_ticks_keywords(function, y_ticks_for_breakpoints=False, extra_xticks=[], extra_yticks=[1]):
    r"""
    Compute ``plot`` keywords for displaying the ticks.
    """
    xticks = function.end_points()
    try:
        f = find_f(function, no_error_if_not_minimal_anyway=True)
        if f is not None and not f in xticks:
            xticks.append(f)
    except ValueError:
        pass
    xticks = sorted(set(list(extra_xticks) + xticks))
    xtick_formatter = [ latex_tick_formatter_x(x) for x in xticks ]
    #xtick_formatter = 'latex'  # would not show rationals as fractions
    ytick_formatter = None
    if y_ticks_for_breakpoints:
        yticks = xticks
    else:
        #yticks = 1/5
        yticks = sorted(set(list(extra_yticks) + [ y for limits in function.limits_at_end_points() for y in limits if y is not None ]))
    ytick_formatter = [ latex_tick_formatter_y(y) for y in yticks ]
    ## FIXME: Can we influence ticks placement as well so that labels don't overlap?
    ## or maybe rotate labels 90 degrees?
    return {'ticks': [xticks, yticks], \

            'gridlines': True, \
            'tick_formatter': [xtick_formatter, ytick_formatter]}

ticks_keywords = default_igp_ticks_keywords

def delete_one_time_plot_kwds(kwds):
    if 'legend_label' in kwds:
        del kwds['legend_label']
    if 'ticks' in kwds:
        del kwds['ticks']
    if 'tick_formatter' in kwds:
        del kwds['tick_formatter']

def plot_covered_intervals(function, covered_components=None, uncovered_color='black', labels=None,
                           show_one_point_overlap_markers=None, **plot_kwds):
    r"""
    Returns a plot of the covered and uncovered intervals of the function.
    """
    if covered_components is None:
        covered_components = generate_covered_components(function)
    uncovered_intervals = uncovered_intervals_from_covered_components(covered_components)
    # Plot the function with different colors.
    # Each component has a unique color.
    # The uncovered intervals is by default plotted in black.
    colors = rainbow(len(covered_components))
    graph = Graphics()
    kwds = copy(ticks_keywords(function))
    kwds.update(plot_kwds)
    if uncovered_intervals:
        kwds.update({'legend_label': "not covered"})
        plot_kwds_hook(kwds)
        graph += plot(function, color = uncovered_color, **kwds)
        delete_one_time_plot_kwds(kwds)
    elif not function.is_continuous(): # to plot the discontinuity markers
        graph += plot(function, color = uncovered_color, **kwds)
        delete_one_time_plot_kwds(kwds)
    if show_one_point_overlap_markers is None:
        show_one_point_overlap_markers = not function.is_continuous()
    for i, component in enumerate(covered_components):
        if labels is None:
            label = "covered component %s" % (i+1)
        else:
            label = labels[i]
        kwds.update({'legend_label': label})
        plot_kwds_hook(kwds)
        last_endpoint = None
        for interval in component:
            linear = function.which_function((interval[0] + interval[1])/2)
            # We do not plot anything if float(interval[0])==float(interval[1]) because
            # otherwise plot complains that
            # "start point and endpoint must be different"
            if float(interval[0])<float(interval[1]):
                graph += plot(linear, (interval[0], interval[1]), color=colors[i], zorder=-1, **kwds)
                # zorder=-1 puts them below the discontinuity markers,
                # above the black function.
                delete_one_time_plot_kwds(kwds)
                # Show a little marker where adjacent intervals of the same component end
                # if the function is continuous at that point.
                # For example, in zhou_two_sided_discontinuous_cannot_assume_any_continuity, or
                # hildebrand_discont_3_slope_1().
                if show_one_point_overlap_markers and interval[0] == last_endpoint:
                    limits = function.limits(last_endpoint)
                    if limits[0] == limits[1] == limits[2]:
                        slope = linear._slope
                        scale = 0.01
                        dx = scale * slope / sqrt(1 + slope**2) # FIXME: this tries to make it orthogonal to the slope
                        dy = -scale / sqrt(1 + slope**2)        # but fails to take the aspect ratio of the plot into account.
                        graph += line([(last_endpoint - dx, function(last_endpoint) - dy), (last_endpoint + dx, function(last_endpoint) + dy)],
                                      color=colors[i], zorder=-1)
            last_endpoint = interval[1]
    return graph

def plot_directly_covered_intervals(function, uncovered_color='black', labels=None, **plot_kwds):
    components = generate_directly_covered_components(function)
    return plot_covered_intervals(function, components, uncovered_color=uncovered_color, labels=labels, **plot_kwds)

def number_of_components(fn):
    r"""
    Returns the number of connected covered components of fn.

    This is an upper bound on ``number_of_slopes``.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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

def slopes_intervals_dict(fn, ignore_nonlinear=False):
    r"""
    Returns a dictionary that maps a slope value to a list of intervals with that slope.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
                if ignore_nonlinear:
                    pass
                else:
                    raise ValueError("Nonlinear piece in function")
    return slopes_dict

def number_of_slopes(fn):
    r"""
    Returns the number of different slopes of fn.

    If fn is discrete, this is defined as the number of different slopes
    of its piecewise linear continuous interpolation.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.WARN)
        sage: number_of_slopes(gmic())
        2
        sage: number_of_slopes(gomory_fractional())
        1
        sage: number_of_slopes(automorphism(restrict_to_finite_group(gmic(10/11)), 3))
        2

    For non-piecewise linear functions, this raises an exception::

        sage: number_of_slopes(california_ip())
        Traceback (most recent call last):
        ...
        ValueError: Nonlinear piece in function
        sage: number_of_slopes(bccz_counterexample())
        Traceback (most recent call last):
        ...
        AttributeError: ... has no attribute 'is_discrete'
        sage: number_of_slopes(bcds_discontinuous_everywhere())
        Traceback (most recent call last):
        ...
        AttributeError: ... has no attribute 'is_discrete'
        sage: number_of_slopes(kzh_minimal_has_only_crazy_perturbation_1_perturbation())
        Traceback (most recent call last):
        ...
        AttributeError: ... has no attribute 'is_discrete'
    """
    if fn.is_discrete():
        fn = interpolate_to_infinite_group(fn)
    return len(slopes_intervals_dict(fn))

def plot_with_colored_slopes(fn, **plot_kwds):
    r"""
    Returns a plot of fn, with pieces of different slopes in different colors.
    """
    slopes_dict = slopes_intervals_dict(fn, ignore_nonlinear=True)
    return plot_covered_intervals(fn, list(slopes_dict.values()), labels=[ "Slope %s" % s for s in slopes_dict.keys() ],
                                  **plot_kwds)

### Minimality check.

def subadditivity_test(fn, full_certificates=True):
    r"""
    Check if fn is subadditive.
    """
    result = True
    for (x, y, z, xeps, yeps, zeps) in generate_nonsubadditive_vertices(fn, reduced=True):
        logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(xeps), y, print_sign(yeps), z, print_sign(zeps)))
        if not full_certificates:
            logging.info("Thus pi is not subadditive.")
            return False
        else:
            result = False
    if result:
        logging.info("pi is subadditive.")
    else:
        logging.info("Thus pi is not subadditive.")
    return result

def symmetric_test(fn, f, full_certificates=True):
    r"""
    Check if fn is symmetric.
    """
    result = True
    if fn(f) != 1:
        logging.info('pi(f) is not equal to 1.')
        if not full_certificates:
            logging.info('pi is not symmetric.')
            return False
        else:
            result = False
    result = True
    for (x, y, xeps, yeps) in generate_nonsymmetric_vertices(fn, f):
        logging.info("pi(%s%s) + pi(%s%s) is not equal to 1" % (x, print_sign(xeps), y, print_sign(yeps)))
        if not full_certificates:
            logging.info('pi is not symmetric.')
            return False
        else:
            result = False
    if result:
        logging.info('pi is symmetric.')
    else:
        logging.info('Thus pi is not symmetric.')
    return result

def find_f(fn, no_error_if_not_minimal_anyway=False):
    r"""
    Find the value of `f` for the given function .
    """
    if hasattr(fn, '_f'):
        return fn._f
    f = None
    for x in fn.end_points():
        if fn(x) > 1 or fn(x) < 0: 
            if no_error_if_not_minimal_anyway:
                logging.info('pi is not minimal because it does not stay in the range of [0, 1].')
                return None
            raise ValueError("The given function does not stay in the range of [0, 1], so cannot determine f.  Provide parameter f to minimality_test or extremality_test.")
    for x in fn.end_points():
        if fn(x) == 1:
            if not f is None:
                logging.info("The given function has more than one breakpoint where the function takes the value 1; using f = %s.  Provide parameter f to minimality_test or extremality_test if you want a different f." % f)
                return f
            else:
                f = x
    if not f is None:
        fn._f = f
        return f
    if no_error_if_not_minimal_anyway:
        logging.info('pi is not minimal because it has no breakpoint where the function takes value 1.')
        return None
    raise ValueError("The given function has no breakpoint where the function takes value 1, so cannot determine f.  Provide parameter f to minimality_test or extremality_test.")

def minimality_test(fn, show_plots=False, f=None, full_certificates=True):
    r"""
    Checks if fn is minimal with respect to the group relaxation with the given `f`.  

    Assume that `0 \leq fn \leq 1`. This function calls ``subadditivity_test`` and ``symmetric_test``.

    If `f` is not provided, use the one found by ``find_f``.

    If show_plots is ``True`` (default: ``False``), show an illustrating diagram.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
    if subadditivity_test(fn, full_certificates=full_certificates) and \
       symmetric_test(fn, f, full_certificates=full_certificates):
        logging.info('Thus pi is minimal.')
        is_minimal = True
    else:
        logging.info('Thus pi is NOT minimal.')
        if not full_certificates:
            return False
        else:
            is_minimal = False
    if show_plots:
        logging.info("Plotting 2d diagram...")
        show_plot(plot_2d_diagram(fn, known_minimal=is_minimal, f=f),
                  show_plots, tag='2d_diagram', object=fn)
        logging.info("Plotting 2d diagram... done")
    return is_minimal

# Global variable to control repr of FastPiecewise.
show_values_of_fastpiecewise =  True

from .fast_linear import FastLinearFunction, fast_linear_function, linear_function_through_points
from .fast_piecewise import FastPiecewise

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
    r"""
    Check if all numbers in the list ``values`` are (or can be converted to) rationals.

    Returns a tuple of two values:

     - True if all rationals
     - a list of values (converted to elements of QQ if True, or the original elements otherwise).

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: is_QQ, QQ_values = is_all_QQ([1, 2/3])
        sage: is_QQ
        True
        sage: QQ_values
        [1, 2/3]
        sage: [ parent(x) for x in QQ_values ]
        [Rational Field, Rational Field]
        sage: is_all_QQ([1, pi])
        (False, [1, pi])
    """
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
    r"""
    Coerce the real numbers in the list symb_values into a convenient common field
    and return a list, parallel to symb_values, of the coerced values.

    If `field` is `False`, the given numbers are returned as is.

    If all given numbers are rational, the field will be the rational
    field (``QQ``).  

    Otherwise, if the numbers are algebraic, the field
    will be a suitable algebraic field extension of the rational
    numbers, embedded into the real numbers, in the form of a
    ``RealNumberField``.  

    Otherwise, the given numbers are returned as is.
    """
    ### Add tests!
    if field is False:
        # do nothing
        return symb_values
    if isinstance(field, ParametricRealField):
        syms = []
        vals = []
        for element in symb_values:
            if is_parametric_element(element):
                syms.append(element.sym())
                vals.append(element.val())
            else:
                syms.append(element)  # changed to not do SR. -mkoeppe
                vals.append(element)
        vals = nice_field_values(vals) #, field=RealNumberField)
        field_values = [ParametricRealFieldElement(field, vals[i],syms[i]) for i in range(len(symb_values))]
        return field_values

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
            embedded_field_values = list(map(hom, number_field_values))
            # Store symbolic expression
            for emb, symb in zip(embedded_field_values, symb_values):
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
    r"""
    Create a continuous piecewise function from bkpt, slopes, and values.

    - bkpt and values are two parallel lists; it is assumed that bkpt is sorted in increasing order. 

    - slopes is one element shorter and represents the slopes of the interpolation.

    - The function is overdetermined by these data.  The consistency of the data is currently not checked.

    - The data are coerced into a common convenient field via ``nice_field_values``.

    - If merge is ``True`` (the default), adjacent pieces of equal slopes are merged into one.
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
    pieces = []
    for i in range(len(bkpt)-1):
        if bkpt[i] > bkpt[i+1]:
            raise ValueError("Breakpoints are not sorted in increasing order.")
        elif bkpt[i] == bkpt[i+1]:  #hasattr(field, '_big_cells') and field._big_cells, should be off-record
            logging.warning("Degenerate interval occurs at breakpoint %s" % bkpt[i])
            if values[i] != values[i+1]:
                raise ValueError("Degeneration leads to a discontinuous function.")
        else:
            pieces.append( [(bkpt[i],bkpt[i+1]),
                            fast_linear_function(slopes[i], intercepts[i])] )
    return FastPiecewise(pieces, merge=merge)

def piecewise_function_from_breakpoints_and_values(bkpt, values, field=None, merge=True):
    r"""
    Create a continuous piecewise function from bkpt and values.

    - bkpt and values are two parallel lists; assuming bpkt is sorted (increasing).

    - The data are coerced into a common convenient field via ``nice_field_values``.

    - If merge is ``True`` (the default), adjacent pieces of equal slopes are merged into one.
    """
    if len(bkpt)!=len(values):
        raise ValueError("Need to have the same number of breakpoints and values.")
    slopes = [ (values[i+1]-values[i])/(bkpt[i+1]-bkpt[i]) if bkpt[i+1] != bkpt[i] else 0 for i in range(len(bkpt)-1) ]
    return piecewise_function_from_breakpoints_slopes_and_values(bkpt, slopes, values, field, merge=merge)

def piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None, merge=True):
    r"""
    Create a continuous piecewise function from bkpt and slopes.

    - bkpt and slopes are two parallel lists (except that bkpt is one element longer); assuming bpkt is sorted (increasing).  The function always has value `0` on ``bkpt[0]``.  

    - The data are coerced into a common convenient field via ``nice_field_values``.

    - If merge is ``True`` (the default), adjacent pieces of equal slopes are merged into one.
    """
    if len(bkpt)!=len(slopes)+1:
        raise ValueError("Need to have one breakpoint more than slopes.")
    values = [0]
    for i in range(1, len(bkpt)):
        values.append(values[i-1] + slopes[i - 1] * (bkpt[i] - bkpt[i-1]))
    return piecewise_function_from_breakpoints_slopes_and_values(bkpt, slopes, values, field, merge=merge)

def piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes, field=None, merge=True):
    r"""
    Create a continuous piecewise function from interval_lengths and slopes.

    - The function always has value 0 on 0. interval_lengths and slopes are two parallel lists that define the function values to the right of 0.

    - The data are coerced into a common convenient field via ``nice_field_values``.

    - If merge is ``True`` (the default), adjacent pieces of equal slopes are merged into one.
    """
    if len(interval_lengths)!=len(slopes):
        raise ValueError("Number of given interval_lengths and slopes needs to be equal.")
    bkpt = []
    bkpt.append(0)
    for i in range(len(interval_lengths)):
        if interval_lengths[i] < 0:
            raise ValueError("Interval lengths must be non-negative.")
        bkpt.append(bkpt[i]+interval_lengths[i])
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field, merge=merge)

def piecewise_function_from_breakpoints_and_limits(bkpt, limits, field=None, merge=True):
    r"""
    Create a continuous or discontinuous piecewise function from bkpt and limits.

    - bkpt and limits are two parallel lists.  Assume that bkpt is a sorted (increasing). limits is a list of tuple of 3 numbers (mid, right, left).

    - The data are coerced into a common convenient field via ``nice_field_values``.

    - If merge is True (the default), adjacent pieces of equal slopes are merged into one.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.WARN) # Suppress output in automatic tests.
        sage: bkpt = [0, 1/8, 3/8, 1/2, 5/8, 7/8, 1]
        sage: limits = [(0, 0, 1/2), (1/4, 1/4, 3/4), (3/4, 1/4, 3/4), (1, 1/2, 1), (3/4, 3/4, 3/4), (1/4, 1/4, 1/4), (0, 0, 1/2)]
        sage: h = piecewise_function_from_breakpoints_and_limits(bkpt, limits)
        sage: h  == hildebrand_discont_3_slope_1()
        True
        sage: h = piecewise_function_from_breakpoints_and_limits(bkpt=[0, 1/5, 2/5, 3/5, 4/5, 1], limits = [{-1:0, 0:0, 1:0},{-1:1, 0:1, 1:1}, {-1:0, 0:2/5, 1:2/5}, {-1:2/5, 0:1/2, 1:3/5}, {-1:3/5, 0:3/5, 1:1}, {-1:0, 0:0, 1:0}])
        sage: h.limit(3/5, 1)
        3/5
        sage: h.limit(3/5, 0)
        1/2
        sage: h.limit(3/5, -1)
        2/5
    """
    if len(bkpt)!=len(limits):
        raise ValueError("Need to have the same number of breakpoints and limits.")
    n = len(bkpt)
    mid, right, left = [limit[0] for limit in limits], [limit[1] for limit in limits], [limit[-1] for limit in limits]
    symb_values = bkpt + mid + right + left
    field_values = nice_field_values(symb_values, field)
    bkpt, mid, right, left = field_values[0:n], field_values[n:2*n], field_values[2*n:3*n], field_values[3*n:4*n]
    pieces = []
    for i in range(n-1):
        pieces += [ singleton_piece(bkpt[i], mid[i]), \
                    open_piece((bkpt[i], right[i]), (bkpt[i+1], left[i+1])) ]
    pieces += [ singleton_piece(bkpt[n-1], mid[n-1]) ]
    return FastPiecewise(pieces, merge=merge)

def piecewise_function_from_breakpoints_slopes_and_jumps(bkpt, slopes, jumps, field=None, merge=True):
    r"""
    Create a continuous or discontinuous piecewise function from bkpt, slopes and jumps.

    - The function always has value 0 on the first breakpoint `0`. The list jumps describes the function value jumps on the left and the right endpoints of each slope.

    - The data are coerced into a common convenient field via ``nice_field_values``.

    - If merge is ``True`` (the default), adjacent pieces of equal slopes are merged into one.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.WARN) # Suppress output in automatic tests.
        sage: bkpt = [0, 1/8, 3/8, 1/2, 5/8, 7/8, 1]
        sage: slopes = [6, 2, 6, 2, -2, 2]
        sage: jumps = [0, -1/2, 0, 0, -1/2, 0, -1/2, 0, 0, 0, 0, -1/2]
        sage: h = piecewise_function_from_breakpoints_slopes_and_jumps(bkpt, slopes, jumps)
    """
    n = len(bkpt)
    if n != len(slopes)+1:
        raise ValueError("Need to have one breakpoint more than slopes.")
    if 2*(n-1) != len(jumps):
        raise ValueError("Need to have number of jumps = 2 * number of slopes.")
    symb_values = bkpt + slopes + jumps
    field_values = nice_field_values(symb_values, field)
    bkpt, slopes, jumps = field_values[0:n], field_values[n:2*n-1], field_values[2*n-1:]
    current_value = 0
    pieces = []
    for i in range(n-1):
        pieces.append([singleton_interval(bkpt[i]), FastLinearFunction(0, current_value)])
        current_value += jumps[2*i]
        pieces.append([open_interval(bkpt[i], bkpt[i+1]), FastLinearFunction(slopes[i], current_value - slopes[i]*bkpt[i])])
        current_value += slopes[i] * (bkpt[i+1] - bkpt[i]) + jumps[2*i+1]
    pieces.append([singleton_interval(bkpt[n-1]), FastLinearFunction(0, current_value)])
    return FastPiecewise(pieces, merge=merge)

def discrete_function_from_points_and_values(points, values, field=None):
    r"""
    Create a function defined on a finite list of points. 

    points and values are two parallel lists.

    The data are coerced into a common convenient field via ``nice_field_values``.
    """
    if field is None:
        field = default_field
    # global symb_values
    symb_values = points + values
    field_values = nice_field_values(symb_values, field)
    points, values = field_values[0:len(points)], field_values[-len(values):]
    pieces = [ (singleton_interval(point), FastLinearFunction(0, value))
               for point, value in zip(points, values) ]
    return FastPiecewise(pieces)

def limiting_slopes(fn):
    r"""
    Computes the limiting slopes on the right and the left side of the
    origin.
    
    The function fn is assumed minimal.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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

use_pwl_template = 'use_pwl_template'
slopes_proportional_to_limiting_slopes_for_positive_epsilon = 'slopes_proportional_to_limiting_slopes_for_positive_epsilon'
slopes_proportional_to_limiting_slopes_for_negative_epsilon = 'slopes_proportional_to_limiting_slopes_for_negative_epsilon'

perturbation_style = use_pwl_template
perturbation_template_bkpts = [0, 1/2, 1]
perturbation_template_values = [0, 1, 0]

def approx_discts_function(perturbation_list, stability_interval, field=default_field, function=None):
    assert (stability_interval.a == - stability_interval.b)
    perturb_points = sorted(perturbation_list.keys())
    if perturbation_style == use_pwl_template:
        template_bkpts = copy(perturbation_template_bkpts)
        template_values = copy(perturbation_template_values)
    elif perturbation_style == slopes_proportional_to_limiting_slopes_for_positive_epsilon:
        if function is None:
            raise ValueError("This perturbation_style needs to know function")
        slope_plus, slope_minus = limiting_slopes(function)
        current_slope = function.which_function(perturb_points[0])._slope
        template_bkpts = [0, (current_slope - slope_minus)/(slope_plus - slope_minus), 1]
        #print slope_plus, slope_minus, current_slope, template_bkpts
        template_values = [0, 1, 0]
    elif perturbation_style == slopes_proportional_to_limiting_slopes_for_negative_epsilon:
        if function is None:
            raise ValueError("This perturbation_style needs to know function")
        slope_plus, slope_minus = limiting_slopes(function)
        current_slope = function.which_function(perturb_points[0])._slope
        template_bkpts = [0, (slope_plus - current_slope)/(slope_plus - slope_minus), 1]
        #print slope_plus, slope_minus, current_slope, template_bkpts
        template_values = [0, 1, 0]
    else:
        raise ValueError("Unknown perturbation_style: %s" % perturbation_style)
    pos_pert_bkpts = [(2 * x - 1) * stability_interval.b  for x in template_bkpts]
    pos_pert_values = [x for x in template_values]
    neg_pert_bkpts =  [(1 - 2 * x) * stability_interval.b  for x in template_bkpts[::-1]]
    neg_pert_values = [-x for x in template_values[::-1]]
    fn_values = [0]
    fn_bkpts = [0]
    for pt in perturb_points:
        sign = perturbation_list[pt][0]
        if sign == 1:
            pert_bkpts = pos_pert_bkpts
            pert_values = pos_pert_values
        else:
            pert_bkpts = neg_pert_bkpts
            pert_values = neg_pert_values
        assert ( pt + pert_bkpts[0] >= fn_bkpts[len(fn_bkpts)-1])
        for i in range(len(pert_bkpts)):
            if i == 0 and pt + pert_bkpts[0] == fn_bkpts[len(fn_bkpts)-1]:
                continue
            fn_bkpts.append(pt + pert_bkpts[i])
            fn_values.append(pert_values[i])
    assert (1 >= fn_bkpts[len(fn_bkpts)-1])
    if (1 > fn_bkpts[len(fn_bkpts)-1]):
        fn_bkpts.append(1)
        fn_values.append(0)
    return piecewise_function_from_breakpoints_and_values(fn_bkpts, fn_values, field)

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

def find_epsilon_interval(fn, perturb):
    r"""
    Compute the largest (closed) interval of values `\epsilon` such
    that ` ``fn`` + \epsilon * ``perturb`` ` is a minimal valid
    function.
    """
    if not hasattr(fn, '_epsilon_interval_dict'):
        fn._epsilon_interval_dict = dict()
    if perturb not in fn._epsilon_interval_dict:
        if fn.is_continuous() or fn.is_discrete():
            fn._epsilon_interval_dict[perturb] = find_epsilon_interval_continuous(fn, perturb)
        else:
            fn._epsilon_interval_dict[perturb] = find_epsilon_interval_general(fn, perturb)
    return fn._epsilon_interval_dict[perturb]

def find_largest_epsilon(fn, perturb):
    r"""
    Compute the proper rescaling of a given perturbation function.

    If the largest epsilon is zero, we should try a different perturbation instead.
    """
    minus_epsilon, plus_epsilon = find_epsilon_interval(fn, perturb)
    return min(abs(minus_epsilon), plus_epsilon)

###
### Moves
###

show_moves_with_discontinuity_markers = True

from .move_semigroup import FullMoveSemigroup, FunctionalDirectedMove

def generate_functional_directed_moves(fn):
    r"""
    Compute the (translations and reflections) moves.
    """
    moves = dict()
    for face in generate_additive_faces(fn) if equiv7_mode else generate_maximal_additive_faces(fn):
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

def generate_all_components(fn, show_plots=False):
    fdms, covered_components = generate_directed_move_composition_completion(fn, show_plots=show_plots)
    if logging.getLogger().isEnabledFor(logging.DEBUG):
        logging.debug("The covered components are %s." % (covered_components))
    uncovered_intervals = generate_uncovered_intervals(fn)
    if uncovered_intervals:
        uncovered_components = generate_uncovered_components(fn, show_plots=show_plots)
        if logging.getLogger().isEnabledFor(logging.DEBUG):
            logging.debug("The uncovered components are %s." % (uncovered_components))
        components = covered_components + uncovered_components
    else:
        components = copy(covered_components)
    return components

def generate_symbolic(fn, components=None, field=None, f=None, show_plots=False, basis_functions=None, **kwds):
    r"""
    Construct a vector-space-valued piecewise linear function
    compatible with the given function ``fn``.

    Two systems of ``basis_functions`` are available:
     - ``('slopes', 'jumps')``
     - ``('midpoints', 'slopes')

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = hildebrand_discont_3_slope_1()
        sage: components =  generate_covered_intervals(h)
        sage: g = generate_symbolic(h, components, field=QQ)
        sage: g.list()
        [[<Int[0, 1/8)>, <FastLinearFunction ((1,0,0,0,0))*x>],
         [(1/8, 3/8), <FastLinearFunction ((0,1,0,0,0))*x + ((1/8,-1/8,0,1,0))>],
         [<Int(3/8, 1/2]>, <FastLinearFunction ((1,0,0,0,0))*x - ((1/4,-1/4,0,-2,0))>],
         [<Int(1/2, 5/8]>, <FastLinearFunction ((0,1,0,0,0))*x + ((1/4,-1/4,0,2,1))>],
         [<Int(5/8, 7/8]>,
          <FastLinearFunction ((0,0,1,0,0))*x + ((1/4,3/8,-5/8,2,1))>],
         [<Int(7/8, 1)>, <FastLinearFunction ((0,1,0,0,0))*x + ((1/4,-1/2,1/4,2,1))>],
         [<Int{1}>, <FastLinearFunction ((1/4,1/2,1/4,2,2))>]]   

        sage: h = extreme_functions.drlm_backward_3_slope()
        sage: sym = generate_symbolic(h, basis_functions=('midpoints', 'slopes'))
        sage: sym.basis
        [('function value at', 1/8),
         ('function value at', 1/6),
         ('slope of component', 0),
         ('slope of component', 1),
         ('slope of component', 2)]

    """
    if field is None:
        field = fn(0).parent().fraction_field()
    if components is None:
        components = generate_all_components(fn)
    if (fn.is_continuous() or fn.is_discrete()) and basis_functions in (None, ('slopes', 'jumps')):
        return generate_symbolic_continuous(fn, components, field=field, f=f, **kwds)
    else:
        return generate_symbolic_general(fn, components, field=field, f=f, basis_functions=basis_functions, **kwds)

def plot_symbolic(sym, indices=None, as_array=True, color='magenta', **kwds):
    """
    EXAMPLES:

    Basis functions for the continuous case.  Note that they are monotonically increasing
    and therefore not periodic functions::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = not_extreme_1()
        sage: sym = generate_symbolic(h)
        sage: plot_symbolic(sym, thickness=2, **ticks_keywords(h)).show(figsize=(6, 6))  # not tested

    Basis functions for the one-sided discontinuous case.  The first three are slopes,
    the last two are jumps.  Again they are monotonically increasing::

        sage: h = hildebrand_discont_3_slope_1()
        sage: sym = generate_symbolic(h)
        sage: plot_symbolic(sym, thickness=2, **ticks_keywords(h)).show(figsize=(6, 6))  # not tested

    Basis functions for the two-sided discontinuous case.  Here we cannot assume continuity at any
    uncovered breakpoint, so nothing is gained by making them increasing functions.  Instead we choose
    them as local on the components (and therefore periodic)::

        sage: h = minimal_no_covered_interval()
        sage: sym = generate_symbolic(h)
        sage: plot_symbolic(sym, thickness=2, **ticks_keywords(h)).show(figsize=(6, 6))  # not tested

    """
    space = sym(0).parent()
    if indices is None:
        indices = range(space.dimension())
    legends = [ 'basis function {}'.format(i) for i in indices ]
    vectors = [ space.gen(i) for i in indices ]
    if as_array:
        plots = [ (v * sym).plot(color=color, legend_label=legend, **kwds)
                  for legend, v in zip(legends, vectors) ]
        return graphics_array(plots, ncols=1)
    else:
        plots = [ (v * sym).plot(color=color, legend_label=legend, **kwds)
                  for legend, v, color in zip(legends, vectors, rainbow(len(vectors))) ]
        return sum(plots)

def generate_additivity_equations(fn, symbolic, field=None, f=None, bkpt=None,
                                  reduce_system=None, return_vertices=False, vertices=None,
                                  undefined_ok=False, **args):
    r"""
    Using additivity, set up a finite-dimensional system of linear equations
    that must be satisfied by any perturbation.

    Use ``vertices`` if provided, otherwise call uses the vertices obtained by
    ``generate_additive_vertices``.

    If ``reduce_system`` is True (the default when ``logging.DEBUG`` is enabled),
    remove redundant equations to make the system minimal.

    If ``return_vertices`` is True, return a list of vertices in the form
    ``(x, y, z, xeps, yeps, zeps)`` corresponding to the rows of the matrix.
    There are two special labels: 'f' and '1' that can appear instead of a vertex.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = hildebrand_discont_3_slope_1()
        sage: components = generate_covered_intervals(h)
        sage: symbolic = generate_symbolic(h, components, field=QQ)
        sage: generate_additivity_equations(h, symbolic, field=QQ, reduce_system=False)
        181 x 5 dense matrix ...
        sage: generate_additivity_equations(h, symbolic, field=QQ, reduce_system=True)
        [ 1/4  1/4    0    2    0]
        [ 1/4  1/2  1/4    2    2]
        [ 1/4  1/2  1/4    3    1]
        [ 1/8 -1/8    0    1    0]
        [   0  1/8 -1/8    2   -1]
        sage: M, vs = generate_additivity_equations(h, symbolic, field=QQ, reduce_system=True, return_vertices=True)
        sage: vs
        ['f', '1', (0, 1/8, 1/8, -1, 0, -1), (1/8, 1/8, 1/4, 0, 0, 0), (3/8, 3/8, 3/4, 1, 1, 1)]

    For two-sided discontinuous functions we can use a different basis, in which we see
    the structure of the subadditivity constraints::

        sage: h = hildebrand_2_sided_discont_2_slope_1()
        sage: symbolic = generate_symbolic(h, basis_functions=('midpoints', 'slopes'))
        sage: M, vs = generate_additivity_equations(h, symbolic, reduce_system=True, return_vertices=True)
        sage: for row, v in zip(M, vs): print("Additive vertex {:30}  {}".format(str(v), row))
        Additive vertex (0, 1/8, 1/8, -1, 1, 0)         (0, -1, 1, -1, 0, 0, 0, 0)
        Additive vertex (0, 1/8, 1/8, 1, -1, 0)         (2, -1, 0, 0, 0, 0, 0, 0)
        Additive vertex (0, 5/8, 5/8, -1, 0, -1)        (0, 0, 0, -2, 1, 0, 0, 0)
        Additive vertex (0, 5/8, 5/8, 1, 0, 1)          (1, 0, 0, 0, 1, -1, 0, 0)
        Additive vertex (1/8, 1/8, 1/4, -1, 1, 0)       (1, 0, 1, 0, 0, 0, -1/16, 1/16)
        Additive vertex (1/8, 1/8, 1/4, -1, 1, 1)       (1, 0, 2, 0, 0, 0, 0, 1/16)
        Additive vertex (1/8, 1/8, 1/4, -1, -1, -1)     (2, 0, -1, 0, 0, 0, -1/16, 1/8)
        Additive vertex (1/8, 5/8, 3/4, 1, 0, 1)        (0, 0, 1, 0, 1, 1, -1/16, 1/16)
    """
    if field is None:
        field = fn(0).parent().fraction_field()
    if f is None:
        f = find_f(fn)
    if not undefined_ok and (fn.is_continuous() or fn.is_discrete()):
        return generate_additivity_equations_continuous(fn, symbolic, field, f=f, bkpt=bkpt,
                                                        reduce_system=reduce_system, return_vertices=return_vertices, vertices=vertices, **args)
    else:
        return generate_additivity_equations_general(fn, symbolic, field, f=f, bkpt=bkpt,
                                                     reduce_system=reduce_system, return_vertices=return_vertices,
                                                     vertices=vertices, undefined_ok=undefined_ok, **args)

def rescale_to_amplitude(perturb, amplitude):
    r"""For plotting purposes, rescale the function perturb so that its
    maximum (supremum) absolute function value is amplitude.
    """
    current_amplitude = max([ abs(x) for limits in perturb.limits_at_end_points() for x in limits if x is not None])
    if current_amplitude != 0:
        return perturb * (amplitude/current_amplitude)
    else:
        return perturb

# Global figsize for all plots made by show_plots.
show_plots_figsize = 10

try:
    from sage.plot.multigraphics import GraphicsArray
except ImportError:
    from sage.plot.graphics import GraphicsArray

def show_plot(graphics, show_plots, tag, object=None, **show_kwds):
    r"""
    Display or save graphics.

    show_plots can be one of: 

    - ``False`` (do nothing), 
    - ``True`` (use ``show`` to display on screen), 
    - a string (file name format such as "FILENAME-%s.pdf", where %s is replaced by tag.
    """
    show_kwds = copy(show_kwds)
    plot_kwds_hook(show_kwds)
    if isinstance(graphics, GraphicsArray):
        # GraphicsArrays can't have legends.
        plot_kwds_hook_no_legend(show_kwds)
    if show_plots_figsize is not None:
        show_kwds['figsize'] = show_plots_figsize
    if isinstance(show_plots, str):
        graphics.save(show_plots % tag, **show_kwds)
    elif show_plots:
        graphics.show(**show_kwds)

def plot_rescaled_perturbation(perturb, xmin=0, xmax=1, **kwds):
    return plot(rescale_to_amplitude(perturb, 1/10), xmin=xmin,
                xmax=xmax, color='magenta', legend_label="perturbation (rescaled)", **kwds)

check_perturbation_plot_three_perturbations = True

def basic_perturbation(fn, index):
    r"""
    Get a basic perturbation of fn.  index counts from 1 (to match the labels in the diagrams). 
    """
    if not hasattr(fn, '_perturbations'):
        extremality_test(fn, show_plots=False)
    if hasattr(fn, '_perturbations'):
        try: 
            return fn._perturbations[index-1]
        except IndexError:
            raise IndexError("Bad perturbation index")
    raise ValueError("No perturbations")

check_perturbation_plot_rescaled_perturbation = True

def plot_perturbation_diagram(fn, perturbation=None, xmin=0, xmax=1):
    r"""
    Plot a perturbation of fn.
    
    perturbation is either a perturbation function, or an integer (which designates a basic perturbation of fn via ``basic_perturbation``).  If perturbation is not provided, it defaults to the perturbation indexed 1.

    To show only a part of the diagram, use::

        sage: plot_perturbation_diagram(h, 1).show(xmin=0.25, xmax=0.35, ymin=0.25, ymax=0.35)  # not tested
    """
    if perturbation is None:
       perturbation = 1
    if isinstance(perturbation, Integer):
        perturbation = basic_perturbation(fn, perturbation)
    epsilon_interval = find_epsilon_interval(fn, perturbation)
    epsilon = min(abs(epsilon_interval[0]), abs(epsilon_interval[1]))
    if check_perturbation_plot_rescaled_perturbation:
        p = plot_rescaled_perturbation(perturbation, xmin=xmin, xmax=xmax)
    else:
        p = Graphics()
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
    epsilon = min(abs(epsilon_interval[0]), abs(epsilon_interval[1]))
    #logging.info("Epsilon for constructed perturbation: %s" % epsilon)
    if show_plots:
        logging.info("Plotting perturbation...")
        p = plot_perturbation_diagram(fn, perturb, xmin=xmin, xmax=xmax)
        show_plot(p, show_plots, tag=show_plot_tag, object=fn, **show_kwds)
        logging.info("Plotting perturbation... done")
    assert epsilon > 0, "Epsilon should be positive, something is wrong"
    #logging.info("Thus the function is not extreme.")  ## Now printed by caller.

def generate_perturbations_finite_dimensional(function, show_plots=False, f=None, full_certificates=True):
    ## FIXME: Perhaps we want an oversampling parameter as in generate_perturbations_simple??
    r"""
    Generate (with ``yield``) perturbations for ``finite_dimensional_extremality_test``.
    """
    # FIXME: fraction_field() required because parent could be Integer
    # Ring.  This happens, for example, for three_slope_limit().  
    # We really should have a function to retrieve the field of
    # a FastPiecewise.  But now even .base_ring() fails because
    # FastLinearFunction does not have a .base_ring() method.
    field = function(0).parent().fraction_field()
    symbolic = generate_symbolic(function, field=field, f=f, show_plots=show_plots)
    bkpt = merge_bkpt(function.end_points(), symbolic.end_points())
    equation_matrix = generate_additivity_equations(function, symbolic, field, f=f, bkpt=bkpt)
    if logging.getLogger().isEnabledFor(logging.DEBUG):
        logging.debug("Solve the linear system of equations:\n%s * v = 0." % (equation_matrix))
    slope_jump_vects = equation_matrix.right_kernel().basis()
    logging.info("Finite dimensional test: Solution space has dimension %s." % len(slope_jump_vects))
    for basis_index in range(len(slope_jump_vects)):
        if not full_certificates:
            yield None
            return
        slope_jump = slope_jump_vects[basis_index]
        logging.debug("The {}-th solution is\nv = {}.".format(basis_index+1, slope_jump))
        perturbation = slope_jump * symbolic
        if logging.getLogger().isEnabledFor(logging.DEBUG):
            logging.debug("The %s-th solution is\nv = %s." % (basis_index+1, slope_jump))
            logging.debug("The corresponding perturbation function pert(x) is:\n%s." % (perturbation))
        yield perturbation

def finite_dimensional_extremality_test(function, show_plots=False, f=None, warn_about_uncovered_intervals=True, 
                                        show_all_perturbations=False, full_certificates=True):
    r"""
    Solve a homogeneous linear system of additivity equations with one
    slope variable for every component (including every non-covered
    interval) and one jump variable for each (left/right) discontinuity.

    Return a boolean that indicates whether the system has a nontrivial solution.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.WARN)
        sage: h1 = drlm_not_extreme_2()
        sage: finite_dimensional_extremality_test(h1, show_plots=True)
        False
        sage: h2 = drlm_3_slope_limit()
        sage: finite_dimensional_extremality_test(h2, show_plots=True)
        True
        sage: h = equiv7_example_1()
        sage: finite_dimensional_extremality_test(h)
        False
        sage: h._perturbations[0].list()
        [[<Int{0}>, <FastLinearFunction 0>],
         [<Int(0, 1/2)>, <FastLinearFunction x - 1/4>],
         [(1/2, 1), <FastLinearFunction 0>]]
    """
    if show_all_perturbations is None:
        show_all_perturbations = show_plots
    if function.is_discrete():
        return simple_finite_dimensional_extremality_test(function, oversampling=1, show_all_perturbations=show_all_perturbations, full_certificates=full_certificates)
    seen_perturbation = False
    function._perturbations = []
    for index, perturbation in enumerate(generate_perturbations_finite_dimensional(function, show_plots=show_plots, f=f, full_certificates=full_certificates)):
        if not full_certificates:
            logging.info("The function is NOT extreme.")
            return False
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
                logging.warning("There are non-covered intervals, so this does NOT prove extremality.")
        else:
            logging.info("Thus the function is extreme.")
    return not seen_perturbation

def generate_facet_covered_components(fn, show_plots=False):
    if not hasattr(fn, "_facet_covered_components"):
        additive_faces_sans_limits = generate_additive_faces_sans_limits(fn)
        fn._facet_covered_components = generate_covered_components_strategically(fn,
                                                                                 show_plots=show_plots,
                                                                                 additive_faces=additive_faces_sans_limits)
    return fn._facet_covered_components

def facet_test(fn, show_plots=False, known_minimal=False, known_extreme=False):
    """
    Test whether `fn` is a facet.

    EXAMPLES:

    Minimal, and all intervals are covered (sans limits), unique solution (sans limits), so a facet::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = extreme_functions.drlm_backward_3_slope(conditioncheck=False)
        sage: facet_test(h)
        True

    Extreme function and continuous piecewise linear, so a facet::

        sage: h = extreme_functions.bhk_irrational()
        sage: facet_test(h)
        True

    The algorithm cannot handle all cases::

        sage: h = kzh_minimal_has_only_crazy_perturbation_1()
        sage: facet_test(h, known_extreme=True)                              # long time - 30s
        Traceback (most recent call last):
        ...
        NotImplementedError: facet_test does not know how to continue

    As a side effect, ``facet_test`` stores some data in the function::

        sage: h._facet_symbolic.basis                                        # long time
        [('function value at', 0.01010000000000000?),
         ('function value at', 0.02020000000000000?), ...,
         ('function value at', 0.8202000000000000?),
         ('slope of component', 0),
         ('slope of component', 1)]
        sage: [h._facet_symbolic.basis.index(('function value at', x)) for x in [h.a0, h.a1, h.a2]]   # long time
        [11, 19, 25]
        sage: h._facet_equation_matrix.column(11)       # random # long time
        (0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0)
        sage: [ (variable, value) for variable, value in zip(h._facet_symbolic.basis, h._facet_equation_matrix[4]) if value ]  # long time
        [(('function value at', 0.01010000000000000?), 1),
         (('function value at', 0.1900000000000000?), 1),
         (('function value at', 0.2071236354684422?), -1),
         (('slope of component', 1), 0.00702363546844223?)]

    """
    if known_extreme:
        known_minimal = True
    if not known_minimal:
        if not minimality_test(fn):
            logging.info("Not minimal, so not a facet")
            return False
    if known_extreme and fn.is_continuous():
        logging.info("Extreme function and continuous piecewise linear, so a facet.")
        return True
    additive_faces_sans_limits = generate_additive_faces_sans_limits(fn)
    covered_components = generate_facet_covered_components(fn)
    uncovered_intervals = uncovered_intervals_from_covered_components(covered_components)
    if uncovered_intervals:
        logging.info("There are non-covered intervals: {}".format(uncovered_intervals))
    # Finite dimensional test.  Default basis functions do not handle the case of uncovered intervals.
    symbolic = fn._facet_symbolic = generate_symbolic(fn, covered_components, basis_functions=('midpoints', 'slopes'))
    vertices = fn._facet_vertices = generate_vertices_of_additive_faces_sans_limits(fn)
    equation_matrix, used_vertices = generate_additivity_equations(fn, symbolic, reduce_system=True,
                                                                   return_vertices=True, undefined_ok=True, vertices=vertices)
    fn._facet_equation_matrix = equation_matrix
    fn._facet_used_vertices = used_vertices
    if logging.getLogger().isEnabledFor(logging.DEBUG):
        logging.debug("Solve the linear system of equations:\n%s * v = 0." % (equation_matrix))
    solution_basis = fn._facet_solution_basis = [ b * symbolic for b in equation_matrix.right_kernel().basis() ]
    logging.info("Finite dimensional test (sans limits): Solution space has dimension %s." % len(solution_basis))
    if not uncovered_intervals and not solution_basis:
        logging.info("Minimal, and all intervals are covered (sans limits), unique solution (sans limits), so a facet.")
        return True
    if not known_extreme and not extremality_test(fn):
        logging.info("Not extreme, so not a facet.")
        return False
    if fn.is_continuous():
        logging.info("Extreme function and continuous piecewise linear, so a facet.")
        return True
    raise NotImplementedError("facet_test does not know how to continue")

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
    r"""Return the list of additive vertices in the form of 6-tuples ``(x, y, z, xeps, yeps, zeps)``.

    When reduced is:

    - ``True``: only outputs fewer triples, for the purpose of setting up the system of equations.
    - ``False``: outputs all triples, for the purpose of plotting ``additive_limit_vertices``.
    """
    # Return a list instead of a generator so that the result can be cached.
    # Using unique_list we remove duplicates.
    return unique_list(itertools.chain( \
                generate_type_1_vertices(fn, operator.eq, reduced=reduced, bkpt=bkpt),\
                generate_type_2_vertices(fn, operator.eq, reduced=reduced, bkpt=bkpt)) )

def generate_vertices_of_additive_faces_sans_limits(fn):
    r"""
    Compute the list of vertices of faces that are additive-sans-limits in the form of 6-tuples ``(x, y, z, xeps, yeps, zeps)``.
    Here ``z`` is not reduced modulo 1.

    TESTS::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = kzh_minimal_has_only_crazy_perturbation_1()
        sage: x = h.end_points()
        sage: (x[6], 1+x[4]-x[6], 1+x[4], 0, 0, 0) in generate_vertices_of_additive_faces_sans_limits(h)
        True

    """
    return unique_list((v[0], v[1], v[0]+v[1], eps[0], eps[1], eps[2])
                       for F in generate_additive_faces_sans_limits(fn)
                       for v in F.vertices
                       for eps in generate_containing_eps_triple(v, F.minimal_triple, with_limits=False))

@cached_function
def generate_nonsubadditive_vertices(fn, reduced=True):
    r"""
    We are returning a set of 6-tuples `(x, y, z, xeps, yeps, zeps)`,
    so that duplicates are removed, and so the result can be cached for later use.

    When reduced is:

    - ``True``: only outputs fewer triples satisfying ``comparison`` relation, for the purpose of ``minimality_test``.
    - ``False``: outputs all triples satisfying ``comparison`` relation, for the purpose of plotting ``nonsubadditive_limit_vertices``.
    """
    return unique_list(itertools.chain( \
                generate_type_1_vertices(fn, operator.lt, reduced=reduced),\
                generate_type_2_vertices(fn, operator.lt, reduced=reduced))  )

def generate_nonsymmetric_vertices(fn, f):
    if fn.is_continuous() or fn.is_discrete():
        return generate_nonsymmetric_vertices_continuous(fn, f)
    else:
        return generate_nonsymmetric_vertices_general(fn, f)

class MaximumNumberOfIterationsReached(Exception):
    pass

crazy_perturbations_warning = False

def extremality_test(fn, show_plots = False, f=None, max_num_it=1000, phase_1=False, finite_dimensional_test_first=False, show_all_perturbations=False, crazy_perturbations=True, full_certificates=True):
    r"""Check if fn is extreme for the group relaxation with the given `f`. 

    If fn is discrete, it has to be defined on a cyclic subgroup of
    the reals containing `1`, restricted to `[0, 1]`.  The group
    relaxation is the corresponding cyclic group relaxation.

    Otherwise fn needs to be defined on the interval `[0, 1]`, and the
    group relaxation is the infinite group relaxation.

    If `f` is not provided, uses the one found by ``find_f()``.

    If show_plots is ``True`` (default: ``False``), show many illustrating diagrams.

    The function first runs ``minimality_test``.
    
    In the infinite group case, if ``finite_dimensional_test_first`` is
    ``True`` (default: ``False``), after testing minimality of fn, we first
    check if the ``finite_dimensional_extremality_test`` finds a
    perturbation; otherwise (default) we first check for an
    equivariant perturbation.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
    global crazy_perturbations_warning
    if crazy_perturbations and \
       limiting_slopes(fn)[0] is +Infinity and \
       limiting_slopes(fn)[1] is -Infinity:
        # this needs checking that it does the right thing for parametric functions.
        crazy_perturbations_warning = True
    else:
        crazy_perturbations_warning = False
    if f is None or not minimality_test(fn, show_plots=show_plots, f=f, full_certificates=full_certificates):
        logging.info("Not minimal, thus NOT extreme.")
        if not phase_1:
            return False
        else:
            do_phase_1_lifting = True
    if do_phase_1_lifting:
        finite_dimensional_test_first = True
    if fn.is_continuous() and number_of_slopes(fn) == 2:
        if not full_certificates:
            logging.info("Gomory-Johnson's 2-slope theorem applies. The function is extreme.")
            return True
        logging.info("Gomory-Johnson's 2-slope theorem applies. The function is extreme.  Continuing anyway because full_certificates=True.")
    seen_perturbation = False
    generator = generate_perturbations(fn, show_plots=show_plots, f=f, max_num_it=max_num_it, finite_dimensional_test_first=finite_dimensional_test_first, full_certificates=full_certificates)
    fn._perturbations = []
    for index, perturbation in enumerate(generator):
        if not full_certificates:
            logging.info("The function is NOT extreme.")
            return False
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

def generate_perturbations(fn, show_plots=False, f=None, max_num_it=1000, finite_dimensional_test_first=False, full_certificates=True):
    r"""
    Generate (with ``yield``) perturbations for ``extremality_test``.
    """
    if fn.is_discrete():
        all = generate_perturbations_simple(fn, show_plots=show_plots, f=f, oversampling=None, full_certificates=full_certificates)
    else:
        fdms, covered_components= generate_directed_move_composition_completion(fn, show_plots=show_plots)
        finite = generate_perturbations_finite_dimensional(fn, show_plots=show_plots, f=f, full_certificates=full_certificates)
        uncovered_intervals = generate_uncovered_intervals(fn)
        if not full_certificates and uncovered_intervals:
            yield None
            return
        if show_plots:
            logging.info("Plotting covered intervals...")
            show_plot(plot_covered_intervals(fn), show_plots, tag='covered_intervals', object=fn)
            logging.info("Plotting covered intervals... done")
        if not uncovered_intervals:
            logging.info("All intervals are covered (or connected-to-covered). %s components." % number_of_components(fn))
            all = finite
        else:
            logging.info("Uncovered intervals: %s", (uncovered_intervals,))
            equi = generate_perturbations_equivariant(fn, show_plots=show_plots, f=f, max_num_it=max_num_it)
            if finite_dimensional_test_first:
                all = itertools.chain(finite, equi)
            else:
                all = itertools.chain(equi, finite)
    for perturbation in all:
        yield perturbation

def generate_perturbations_equivariant(fn, show_plots=False, f=None, max_num_it=1000):
    if fn.is_two_sided_discontinuous():
        logging.warning("Code for detecting perturbations using moves is EXPERIMENTAL in the two-sided discontinuous case.")
    generator = generate_generic_seeds_with_completion(fn, show_plots=show_plots, max_num_it=max_num_it) # may raise MaximumNumberOfIterationsReached
    #seen_perturbation = False
    for seed, stab_int, walk_list in generator:
        # for debugging only:
        #global last_seed, last_stab_int, last_walk_list = seed, stab_int, walk_list
        perturb = approx_discts_function(walk_list, stab_int, function=fn)
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
    r"""
    Return a plot of the completion diagram.
    
    To view a part only, use::

        sage: plot_completion_diagram(h).show(xmin=0.3, xmax=0.55, ymin=0.3, ymax=0.55) # not tested
    """
    if not (hasattr(fn, '_completion') and fn._completion.is_complete):
        extremality_test(fn, show_plots=False)
    if fn._completion.plot_background is None:
        fn._completion.plot_background = plot_completion_diagram_background(fn)
    g = None
    if perturbation is None:
        if hasattr(fn, '_perturbations') and fn._perturbations:
            perturbation = fn._perturbations[0]
    elif isinstance(perturbation, Integer):
        perturbation = basic_perturbation(fn, perturbation)
    if perturbation is not None:
        g = plot_function_at_borders(rescale_to_amplitude(perturbation, 1/10), color='magenta', legend_label='perturbation (rescaled)')
    if hasattr(perturbation, '_walk_list'):
        g += plot_walk_in_completion_diagram(perturbation._seed, perturbation._walk_list)
    return fn._completion.plot(extra_graphics = g)

def perturbation_polyhedron(fn, perturbs, undefined_ok=False, **kwds):
    r"""
    Given fn  and a list of basic perturbations that are pwl, satisfing the symmetry condition and pert(0)=pert(f)=0. Set up a polyhedron, one dimension for each basic perturbation, with the subadditivities.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO) # to disable output in automatic tests.
        sage: h = not_extreme_1()
        sage: finite_dimensional_extremality_test(h, show_all_perturbations=True)
        False
        sage: perturbs = h._perturbations
        sage: len(perturbs)
        2
        sage: pert_polyhedron = perturbation_polyhedron(h, perturbs)
        sage: pert_polyhedron
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices
        sage: pert_polyhedron.Hrepresentation()
        (An inequality (-9, -6) x + 20 >= 0,
         An inequality (6, 9) x + 20 >= 0,
         An inequality (3, 0) x + 4 >= 0,
         An inequality (0, -3) x + 4 >= 0)
        sage: pert_polyhedron.Vrepresentation()
        (A vertex at (20/3, -20/3),
         A vertex at (4/3, 4/3),
         A vertex at (-4/3, 4/3),
         A vertex at (-4/3, -4/3))

    Lift function by adding a perturbtion that corresponds to the vertex (-4/3, 4/3), i.e., set ``h_lift = h - 4/3*h._perturbations[0] + 4/3*h._perturbations[1]``. The lifted function is extreme::

        sage: vertex = pert_polyhedron.vertices()[2]
        sage: perturbation = perturbation_corresponding_to_vertex(perturbs, vertex)
        sage: h_lift = h + perturbation
        sage: extremality_test(h_lift)
        True

    The following function has irrational data::

        sage: h = chen_tricky_uncovered_intervals()
        sage: finite_dimensional_extremality_test(h, show_all_perturbations=True)
        False
        sage: perturbs = h._perturbations
        sage: pert_polyhedron = perturbation_polyhedron(h, perturbs)
        sage: pert_polyhedron
        A 2-dimensional polyhedron in (Real Number Field in `a` as the root of the defining polynomial y^2 - 3 near 1.732050807568878?)^2 defined as the convex hull of 4 vertices
        sage: pert_polyhedron.Vrepresentation()
        (A vertex at (-2.36220486286011?, 0.967307929548895?),
         A vertex at (2.797434948471088?, 0.967307929548895?),
         A vertex at (-3.61183490350498?, -1.248914311409209?),
         A vertex at (1.79481389229748?, -3.45932770938056?))

    The following function is 2-sided discontinous at the origin::

        sage: h = zhou_two_sided_discontinuous_cannot_assume_any_continuity()
        sage: import cutgeneratingfunctionology.igp as igp
        sage: igp.generate_symbolic_two_sided_discontinuous_basis_functions = ('slopes', 'jumps')   # Restore classic behavior
        sage: finite_dimensional_extremality_test(h, show_all_perturbations=True)
        False
        sage: perturbs = h._perturbations
        sage: pert_polyhedron = perturbation_polyhedron(h, perturbs)
        sage: pert_polyhedron
        A 1-dimensional polyhedron in QQ^1 defined as the convex hull of 2 vertices
        sage: pert_polyhedron.Vrepresentation()
        (A vertex at (4/3), A vertex at (-4/9))
        sage: h_lift = h + perturbation_corresponding_to_vertex(perturbs, pert_polyhedron.vertices()[0])
        sage: extremality_test(h_lift)
        True
    """
    (ieqs, eqns) = perturbation_polyhedron_ieqs_eqns(fn, perturbs, undefined_ok=undefined_ok)
    pert_polyhedron = Polyhedron(ieqs = ieqs, eqns = eqns, **kwds)
    return pert_polyhedron

def perturbation_polyhedron_ieqs_eqns(fn, perturbs, undefined_ok=False):
    bkpt = copy(fn.end_points())
    for pert in perturbs:
        bkpt += pert.end_points()
    bkpt = unique_list(bkpt)
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]
    type_1_vertices = ((x, y, x+y) for x in bkpt for y in bkpt if x <= y)
    type_2_vertices = ((x, z-x, z) for x in bkpt for z in bkpt2 if x < z < 1+x)
    vertices = unique_list(itertools.chain(type_1_vertices,type_2_vertices))
    if fn.is_continuous() and all(pert.is_continuous() for pert in perturbs):
        limitingeps = []
        lim_xeps = [0]
    else:
        limitingeps = list(nonzero_eps) # nonzero_eps is defined in discontinuous_case.sage
        lim_xeps = [0, 1, -1]
    ieqset = []
    eqnset = []
    # assume that the basic perturbations come from finite_dimensional_extremality_test(), so that the symmetry constraints and the condition pert(0)=pert(f)=0 are always satisfied.
    # record the constraints 0 <= fn(x) + pert(x) <= 1 for any breakpoint x.
    # only need record the side >=0,  the other side is implied by symmetry.
    for x in bkpt:
        for xeps in lim_xeps:
            valuefn = fn.limit(x, xeps)
            try:
                valuep = [pert.limit(x, xeps) for pert in perturbs]
            except ValueError:
                if undefined_ok:
                    continue
                else:
                    raise
            constraint_coef = tuple([valuefn]) + tuple(valuep)
            ieqset.append(constraint_coef)
    # record the subadditivity constraints
    for (x, y, z) in vertices:
        for (xeps, yeps, zeps) in [(0,0,0)]+limitingeps:
            deltafn = delta_pi_general(fn, x, y, (xeps, yeps, zeps))
            try:
                deltap = [delta_pi_general(pert, x, y, (xeps, yeps, zeps)) for pert in perturbs]
            except ValueError:
                if undefined_ok:
                    continue
                else:
                    raise
            constraint_coef = tuple([deltafn]) + tuple(deltap)
            ieqset.append(constraint_coef)
            # if deltafn > 0:
            #     ieqset.append(constraint_coef)
            # else:
            #     eqnset.append(constraint_coef)
            #     # this is always true for basic perturbations coming from finite_dimensional_extremality_test().
    return unique_list(ieqset), unique_list(eqnset)

def perturbation_mip(fn, perturbs, solver=None, field=None, undefined_ok=False):
    r"""
    Given fn and a list of basic perturbations that are pwl, satisfing the symmetry condition and pert(0)=pert(f)=0. Set up a mip, one dimension for each basic perturbation, with the subadditivities.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO) # to disable output in automatic tests.
        sage: h = not_extreme_1()
        sage: finite_dimensional_extremality_test(h, show_all_perturbations=True)
        False
        sage: len(h._perturbations)
        2
        sage: pert_mip = perturbation_mip(h, h._perturbations,'ppl')

    We set ``solver=ppl`` here.  Note that we can also set ``solver=InteractiveLP``. The coefficients in the constraints are rational numbers, rather than ``float`` used by the default GLPK solver::

        sage: pert_mip.show()
        Maximization:
        <BLANKLINE>
        Constraints:
          constraint_0: -1/10 x_0 <= 1/3
          constraint_1: 1/20 x_0 <= 1/3
          constraint_2: -1/20 x_0 <= 2/3
          constraint_3: 1/10 x_0 <= 2/3
          constraint_4: -1/10 x_1 <= 2/3
          constraint_5: 1/20 x_1 <= 2/3
          constraint_6: -1/20 x_1 <= 1/3
          constraint_7: 1/10 x_1 <= 1/3
          constraint_8: -1/4 x_0 <= 1/3
          constraint_9: -1/10 x_0 + 1/10 x_1 <= 2/3
          constraint_10: -1/10 x_0 - 3/20 x_1 <= 1/3
          constraint_11: 3/20 x_0 + 1/10 x_1 <= 1/3
          constraint_12: 1/20 x_0 - 1/20 x_1 <= 2/3
          constraint_13: -1/20 x_0 + 1/20 x_1 <= 4/3
          constraint_14: -1/20 x_0 - 1/5 x_1 <= 1
          constraint_15: 1/5 x_0 + 1/20 x_1 <= 1
          constraint_16: 1/10 x_0 - 1/10 x_1 <= 4/3
          constraint_17: 1/4 x_1 <= 1/3
        Variables:
          x_0 is a continuous variable (min=-oo, max=+oo)
          x_1 is a continuous variable (min=-oo, max=+oo)

    Since rational coefficients are used in ``ppl`` solver, we can ask for the polyhedron defined by the Linear Program. This would fail if we set ``solver=GLPK`` and if coefficient are not integers, due to ``AttributeError: type object 'float' has no attribute 'fraction_field'``::

        sage: pert_poly = pert_mip.polyhedron()
        sage: pert_poly
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices
        sage: vertices = pert_poly.vertices()
        sage: vertices
        (A vertex at (20/3, -20/3),
         A vertex at (4/3, 4/3),
         A vertex at (-4/3, 4/3),
         A vertex at (-4/3, -4/3))

    Lifting the function by adding a perturbation that corresponds to a vertex, we obtain an extreme function::

        sage: h_lift = h + perturbation_corresponding_to_vertex(h._perturbations, vertices[2])
        sage: extremality_test(h_lift)
        True
    """
    bkpt = copy(fn.end_points())
    for pert in perturbs:
        bkpt += pert.end_points()
    bkpt = unique_list(bkpt)
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]
    type_1_vertices = ((x, y, x+y) for x in bkpt for y in bkpt if x <= y)
    type_2_vertices = ((x, z-x, z) for x in bkpt for z in bkpt2 if x < z < 1+x)
    vertices = unique_list(itertools.chain(type_1_vertices,type_2_vertices))
    if fn.is_continuous() and all(pert.is_continuous() for pert in perturbs):
        limitingeps = []
        lim_xeps = [0]
    else:
        limitingeps = list(nonzero_eps) # nonzero_eps is defined in discontinuous_case.sage
        lim_xeps = [0, 1, -1]
    if field is None:
        field = bkpt[0].parent().fraction_field()
    mip = MixedIntegerLinearProgram(solver=solver, base_ring=field)
    n = len(perturbs)
    constraints_set = set([])
    # assume that the basic perturbations come from finite_dimensional_extremality_test(), so that the symmetry constraints and the condition pert(0)=pert(f)=0 are always satisfied.
    # record the constraints 0 <= fn(x) + pert(x) <= 1 for any breakpoint x.
    # only need record the side >=0,  the other side is implied by symmetry.
    for x in bkpt:
        for xeps in lim_xeps:
            valuefn = fn.limit(x, xeps)
            try:
                valuep = [pert.limit(x, xeps) for pert in perturbs]
            except ValueError:
                if undefined_ok:
                    continue
                else:
                    raise
            if all(coef == 0 for coef in valuep):
                continue
            constraint_coef = tuple([valuefn]) + tuple(valuep)
            if not constraint_coef in constraints_set:
                constraints_set.add(constraint_coef)
                constraint_linear_func = sum(valuep[i] * mip[i] for i in range(n))
                mip.add_constraint(constraint_linear_func + valuefn >= 0)
    # record the subadditivity constraints
    for (x, y, z) in vertices:
        for (xeps, yeps, zeps) in [(0,0,0)]+limitingeps:
            deltafn = delta_pi_general(fn, x, y, (xeps, yeps, zeps))
            try:
                deltap = [delta_pi_general(pert, x, y, (xeps, yeps, zeps)) for pert in perturbs]
            except ValueError: # undefined
                if undefined_ok:
                    continue
                else:
                    raise
            if all(coef == 0 for coef in deltap):
                continue
            constraint_coef = tuple([deltafn]) + tuple(deltap)
            constraint_linear_func = sum(deltap[i] * mip[i] for i in range(n))
            if constraint_coef in constraints_set:
                # don't add duplicated constraints.
                # `if constraint_linear_func in mip.constraints:` doesn't work.
                continue
            else:
                constraints_set.add(constraint_coef)
            mip.add_constraint(constraint_linear_func + deltafn >= 0)
            # if deltafn > 0:
            #     mip.add_constraint(constraint_linear_func + deltafn >= 0)
            # else: # deltafn == 0
            #     mip.add_constraint(constraint_linear_func == 0)
            #     # this is always true for basic perturbations coming from finite_dimensional_extremality_test().
    return mip

def generate_lifted_functions(fn, perturbs=None, use_polyhedron=False, **kwds):
    r"""
    A generator of lifted functions.

    If ``use_polyhedron=False`` (the default), set up an LP (as a 
    ``MixedIntegerLinearProgram``) with one dimension for 
    each basic perturbation, with the subadditivities. Shoot random directions 
    as objective functions. Solve the MIP. Lift the function by adding the 
    perturbation that corresponds to the MIP solution. MIP solver keywords can be 
    passed in ``kwds``.

    If ``use_polyhedron=True``, set up a ``Polyhedron`` instead.  ``Polyhedron``
    constructor keywords can be passed in ``kwds``.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.WARN) # to disable output in automatic tests.
        sage: h = not_extreme_1()
        sage: h_lift = next(generate_lifted_functions(h, solver='ppl'))
        sage: extremality_test(h_lift)
        True

    The above mip problem use ``ppl`` backend. We can use other backends,
    such as ``solver=InteractiveLP``::

        sage: h = not_extreme_1()
        sage: gen = generate_lifted_functions(h, solver='InteractiveLP')
        sage: h_lift = next(gen)
        sage: extremality_test(h_lift)
        True

    If the solver argument is not specified, the code can figure it out,
    using field (``base_ring``).

    The option ``solver=InteractiveLP`` is able to deal with irrational numbers::

        sage: h = chen_tricky_uncovered_intervals()
        sage: gen = generate_lifted_functions(h, perturbs=None, solver='InteractiveLP', field=None)
        sage: h_lift = next(gen)
        sage: extremality_test(h_lift)
        True

    By setting ``use_polyhedron=True``, we use ``perturbation_polyhedron()`` rather than ``perturbation_mip()`` to generate lifted functions::

        sage: h = not_extreme_1()
        sage: gen = generate_lifted_functions(h, use_polyhedron=True)
        sage: len([h_lift for h_lift in gen])
        4
    """
    if perturbs is None:
        #if not hasattr(fn, '_perturbations'):
        # do extremeality test anyway, because if previously show_all_perturbations=False, only one perturbation is stored in fn._perturbations.
        #finite_dimensional_extremality_test or extremality_test?
        #extremality_test(fn, show_all_perturbations=True)
        finite_dimensional_extremality_test(fn, show_all_perturbations=True)
        perturbs = fn._perturbations
    if use_polyhedron:
        pert_polyhedron = perturbation_polyhedron(fn, perturbs, **kwds)
        vertices = pert_polyhedron.vertices()
    else:
        pert_mip = perturbation_mip(fn, perturbs, **kwds)
        vertices = generate_random_mip_sol(pert_mip)
    for vertex in vertices:
        logging.info("vertex = %s" % str(vertex))
        perturb = perturbation_corresponding_to_vertex(perturbs, vertex)
        yield fn + perturb

def perturbation_corresponding_to_vertex(perturbs, vertex):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO) # to disable output in automatic tests.
        sage: h = not_extreme_1()
        sage: finite_dimensional_extremality_test(h, show_all_perturbations=True)
        False
        sage: perturbs = h._perturbations
        sage: pert_mip = perturbation_mip(h, perturbs, 'ppl')
        sage: pert_mip.solve()
        0
        sage: mip_sol = pert_mip.get_values([pert_mip[0], pert_mip[1]])
        sage: mip_sol
        [-4/3, -4/3]
        sage: perturbation = perturbation_corresponding_to_vertex(perturbs, mip_sol)
        sage: h_lift = h + perturbation
        sage: extremality_test(h_lift)
        True
    """
    n = len(perturbs)
    perturb = vertex[0] * perturbs[0]
    for i in range(1, n):
        perturb += vertex[i] * perturbs[i]
    return perturb

def solve_mip_with_random_objective_function(mip):
    n = mip.number_of_variables()
    obj_fun = 0
    for i in range(n):
        random_coeff = QQ(sage.misc.prandom.random() - 1/2)
        # integer coefficient ZZ.random_element() was used, due to a printing error, but this has been solved.
        # When objective function has zero coefficient, solver='InteractiveLP'
        # sometimes gives non-vertex optimial solution, which comes from the standard-form back transformation of a vertex.
        # To avoid such case, we generate another random coefficient.
        while random_coeff == 0:
            random_coeff = QQ(sage.misc.prandom.random() - 1/2) #ZZ.random_element()
        obj_fun += random_coeff * mip[i]
    mip.set_objective(obj_fun)
    opt_val = mip.solve()
    opt_sol = mip.get_values([mip[i] for i in range(n)])
    return tuple(opt_sol)

def generate_random_mip_sol(mip):
    seen_solutions = set([])
    while True:
        mip_sol = solve_mip_with_random_objective_function(mip)
        if not mip_sol in seen_solutions:
            seen_solutions.add(mip_sol)
            yield(mip_sol)

def lift(fn, show_plots = False, use_all_perturbations=True, use_largest_absolute_epsilon=True, **kwds):
    if use_all_perturbations:
        kwds['show_all_perturbations'] = True
    if not 'finite_dimensional_test_first' in kwds:
        kwds['finite_dimensional_test_first'] = True
    if not hasattr(fn, '_perturbations') and extremality_test(fn, show_plots=show_plots, crazy_perturbations=False, **kwds):
        return fn
    else:
        perturbed = fn
        for perturbation in fn._perturbations:
            epsilon_interval = find_epsilon_interval(perturbed, perturbation)
            if epsilon_interval[0] > epsilon_interval[1]:
                # TODO: Find an epsilon that decrease the subadditivity violation and lift?
                continue
            if abs(epsilon_interval[0]) > abs(epsilon_interval[1]):
                which_perturbation = 0
            else:
                which_perturbation = 1
            if not use_largest_absolute_epsilon:
                which_perturbation = 1 - which_perturbation
            perturbed = fn._lifted = perturbed + epsilon_interval[which_perturbation] * perturbation
            ## Following is strictly experimental: It may change what "f" is.
            if 'phase_1' in kwds and kwds['phase_1']:
                perturbed = rescale_to_amplitude(perturbed, 1)
            if not use_all_perturbations:
                break
        if perturbed == fn:
            logging.info("Lifting fails. Try generate_lifted_functions() via polyhedron of perturbation space.")
        return perturbed

def lift_extreme_function_for_finite_group_to_infinite_group(fn, show_plots = False, show_all_lifting=True):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = drlm_not_extreme_1()
        sage: hh = lift_extreme_function_for_finite_group_to_infinite_group(h)
        sage: len(hh)
        2
        sage: h = example7slopecoarse2()
        sage: hh = lift_extreme_function_for_finite_group_to_infinite_group(h)
        sage: len(hh)
        12
        sage: bkpts = [0, 1/13, 3/13, 7/26, 4/13,5/13,21/52,23/52,6/13,8/13,33/52,35/52,9/13,10/13,21/26,11/13,1]
        sage: values = [0,1,3/14,5/7,3/4,5/14,55/112,33/112,3/7,4/7,79/112,57/112,9/14,1/4,2/7,11/14,0]
        sage: h = piecewise_function_from_breakpoints_and_values(bkpts, values)
        sage: hh = lift_extreme_function_for_finite_group_to_infinite_group(h)
        sage: len(hh)
        2
        sage: extremality_test(hh[0])
        True

    BUG EXAMPLES::

        sage: h = FastPiecewise([[(QQ(0), 1/18), FastLinearFunction(QQ(18), QQ(0))], [(1/18, 1/9), FastLinearFunction(-126/11, 18/11)], [(1/9, 1/6), FastLinearFunction(18/55, 18/55)], [(1/6, 2/9), FastLinearFunction(342/55, -36/55)], [(2/9, 5/18), FastLinearFunction(-126/11, 36/11)], [(5/18, 1/3), FastLinearFunction(666/55, -36/11)], [(1/3, 7/18), FastLinearFunction(-306/55, 144/55)], [(7/18, 4/9), FastLinearFunction(18/55, 18/55)], [(4/9, 1/2), FastLinearFunction(342/55, -126/55)], [(1/2, 5/9), FastLinearFunction(-126/11, 72/11)], [(5/9, 11/18), FastLinearFunction(342/55, -36/11)], [(11/18, 2/3), FastLinearFunction(18/55, 18/55)], [(2/3, 13/18), FastLinearFunction(-306/55, 234/55)], [(13/18, 7/9), FastLinearFunction(666/55, -468/55)], [(7/9, 5/6), FastLinearFunction(-126/11, 108/11)], [(5/6, 8/9), FastLinearFunction(342/55, -54/11)], [(8/9, 17/18), FastLinearFunction(18/55, 18/55)], [(17/18, QQ(1)), FastLinearFunction(-126/11, 126/11)]])
        sage: lift_extreme_function_for_finite_group_to_infinite_group(h)
        []
    """
    bkpts = fn.end_points()
    values = fn.values_at_end_points()
    assert is_all_QQ_fastpath(bkpts + values) == True
    q = lcm([denominator(x) for x in bkpts])
    v = lcm([denominator(x) for x in values])
    f = find_f(fn, no_error_if_not_minimal_anyway=True)
    fdms, covered_components= generate_directed_move_composition_completion(fn, show_plots=show_plots)
    for pert_fin in  generate_perturbations_finite_dimensional(fn, show_plots=show_plots, f=f):
        raise ValueError("The given function is not extreme for the finite group (1/%s)Z" % q)
    uncovered_intervals = generate_uncovered_intervals(fn)
    if show_plots:
        logging.info("Plotting covered intervals...")
        show_plot(plot_covered_intervals(fn), show_plots, tag='covered_intervals', object=fn)
        logging.info("Plotting covered intervals... done")
    if not uncovered_intervals:
        logging.info("All intervals are covered (or connected-to-covered). %s components." % number_of_components(fn))
        return [fn]
    logging.info("Uncovered intervals: %s", (uncovered_intervals,))
    all_lifting = [fn]
    global perturbation_template_bkpts
    global perturbation_template_values
    generator = generate_generic_seeds_with_completion(fn)
    for seed, stab_int, walk_list in generator:
        lifted_functions = []
        for perturbed in all_lifting:
            perturbation_template_values = [0, 1, 0]
            perturbation_template_bkpts = [0, 1/(4*q*v), 1]
            perturb_1 = approx_discts_function(walk_list, stab_int, function=fn)
            epsilon_interval_1 = find_epsilon_interval(perturbed, perturb_1)
            perturbation_template_bkpts = [0, 1-1/(4*q*v), 1]
            perturb_2 = approx_discts_function(walk_list, stab_int, function=fn)
            epsilon_interval_2 = find_epsilon_interval(perturbed, perturb_2)
            if epsilon_interval_1[1] != 10000 and epsilon_interval_2[1] != 10000:
                s_left = 4*q*v * epsilon_interval_1[1]
                s_right = - 4*q*v * epsilon_interval_2[1]
                if s_left == s_right == 0:
                    lifted_functions.append(perturbed)
                else:
                    perturbation_template_bkpts = [0, s_right/(s_right - s_left), 1]
                    perturbation_template_values = [0, s_left*s_right/(s_right - s_left), 0]
                    perturbation = approx_discts_function(walk_list, stab_int, function=fn)
                    epsilon_interval = find_epsilon_interval(perturbed, perturbation)
                    #assert epsilon_interval[1] == 1
                    if epsilon_interval[1] == 1:
                        lifted_functions.append(perturbed + perturbation)
            if (show_all_lifting or not lifted_functions) and (epsilon_interval_1[0] != -10000 and epsilon_interval_2[0] != -10000):
                s_left = 4*q*v * epsilon_interval_1[0]
                s_right = - 4*q*v * epsilon_interval_2[0]
                if s_left == s_right == 0:
                    if not perturbed == lifted_functions[-1]:
                        lifted_functions.append(perturbed)
                else:
                    perturbation_template_bkpts = [0, s_right/(s_right - s_left), 1]
                    perturbation_template_values = [0, s_left*s_right/(s_right - s_left), 0]
                    perturbation = approx_discts_function(walk_list, stab_int, function=fn)
                    epsilon_interval = find_epsilon_interval(perturbed, perturbation)
                    #assert epsilon_interval[1] == 1
                    if epsilon_interval[1] == 1:
                        lifted_functions.append(perturbed + perturbation)
        if not lifted_functions:
                raise ValueError("can not find continuous perturbation")
        all_lifting = lifted_functions
    lifted_functions = []
    for h in all_lifting:
        if not generate_uncovered_components(h): #extremality_test(h):
            lifted_functions.append(h)
            if show_plots:
                perturbation = h - fn
                p = plot_rescaled_perturbation(perturbation)
                p += plot(perturbed, color='black')
                p += plot_with_colored_slopes(h)
                p.show()
    return lifted_functions

def lift_until_extreme(fn, show_plots = False, pause = False, covered_length=True, use_all_perturbations=True, **kwds):
    next, fn = fn, None
    while next != fn:
        if covered_length:
            covered_components = generate_covered_components(next)
            covered_intervals = union_of_coho_intervals_minus_union_of_coho_intervals(covered_components, [])
            covered_length = sum(interval_length(i) for i in covered_intervals)
            if show_plots:
                plot_covered_intervals(next).show(title = 'covered length = %s' % covered_length)
            logging.info("Covered length = %s" % covered_length)
        fn = next
        next = lift(fn, show_plots=False , use_all_perturbations=use_all_perturbations, **kwds)  #show_plots=show_plots
        if pause and next != fn:
            input("Press enter to continue")
    return next

##############
def last_lifted(fn):
    while hasattr(fn, '_lifted'):
        fn = fn._lifted
    return fn

def piecewise_function_from_robert_txt_file(filename):
    r"""The .txt files have 4 rows.  
    
    - 1st row = `Y` values

    - 2nd row = `X` values (I do not use these, but I included them in case you want them)

    - 3rd row = `f` (the `x` coordinate for which I use as `f`)
    
    - 4th row = value at `f`  (I do not normalize this to 1.  This allows the `Y` values to range from 0 to this values)

    Also, I do not include the last value (`\pi(1)`) ever because this is
    the same as `\pi(0)` due to periodicity.  So, if you need this last
    value, please attach a 0 to the end of the `Y` values and an extra `x`
    value.
    """
    with open(filename) as f:
        yvalues = [QQ(x) for x in f.readline().split()]
        xvalues = [QQ(x) for x in f.readline().split()]
        if xvalues != list(range(len(yvalues))):
            raise ValueError("Line 2 (xvalues) need to be consecutive integers")
        xscale = len(xvalues)
        xf = QQ(f.readline())
        yf = QQ(f.readline())
    if yvalues[xf] != yf:
        raise ValueError("Lines 3/4 on f and value at f are not consistent with line 1.")
    return piecewise_function_from_breakpoints_and_values([ x / xscale for x in xvalues ] + [1], [y / yf for y in yvalues] + [0])

def random_piecewise_function(xgrid=10, ygrid=10, continuous_proba=1, symmetry=True):
    r"""
    Return a random, continuous or discontinuous piecewise linear function defined on `[0, 1]` with breakpoints that are multiples of `\\frac{1}{xgrid}` and values that are multiples of `\\frac{1}{ygrid}`.

    - continuous_proba (a real number in `[0,1]`) indicates the probability that the function is (left/right) continuous at a breakpoint. 
    
    - Use ``continuous_proba = 1`` (the default) to get a continuous piecewise linear function.

    - Use ``symmetry = True`` (the default) to get a symmetric function. 

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
        for i in range(0, xgrid):
            p = sage.misc.prandom.random()
            if p > continuous_proba:
                rightlimits.append(randint(0, ygrid) / ygrid)
            else:
                rightlimits.append(yvalues[i])
            p = sage.misc.prandom.random()
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

import six

def is_all_QQ_fastpath(values):
    r"""
    This version does not do the full check whether it can be coerced to ``QQ``,
    which is slow for ``RealNumberField``.
    """
    for x in values:
        if not (isinstance(x, (Rational, Integer))
                or type(x) in six.integer_types):
            return False
    return True

from sage.rings.number_field.number_field_element import is_NumberFieldElement

def is_all_the_same_number_field_fastpath(values):
    r"""
    This version does not try coercions and compares fields using ``is``, rather than their comparison operator.
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
    r"""
    Test if numbers are linearly independent over ``QQ``.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
        sage: is_QQ_linearly_independent(pi)
        True
        sage: is_QQ_linearly_independent(1,pi)
        Traceback (most recent call last):
        ...
        ValueError: Q-linear independence test only implemented for algebraic numbers
    """
    # trivial cases
    if len(numbers) == 0:
        return True
    elif len(numbers) == 1:
        return bool(numbers[0] != 0)
    # fast path for rationals
    if is_all_QQ_fastpath(numbers):
        return False
    # param
    if any(is_parametric_element(x) for x in numbers):
        field = None
        is_independent = is_QQ_linearly_independent(*(result_concrete_value(field, numbers)))
        #numbers[0].parent().record_independence_of_pair(numbers, is_independent)
        return is_independent
    if not is_all_the_same_number_field_fastpath(numbers):
        # try to coerce to common number field
        numbers = nice_field_values(numbers, RealNumberField)
        if not is_NumberFieldElement(numbers[0]):
            is_QQ, QQ_numbers = is_all_QQ(numbers)
            if is_QQ:
                return False
            raise ValueError("Q-linear independence test only implemented for algebraic numbers")
    coordinate_matrix = matrix(QQ, [x.list() for x in numbers])
    return rank(coordinate_matrix) == len(numbers)

def compose_functional_directed_moves(A, B, show_plots=False):
    r"""
    Compute the directed move that corresponds to the directed move `A` after `B`.
    
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: merge_functional_directed_moves(
        ....:    FunctionalDirectedMove([(3/10, 7/20), (9/20, 1/2)], (1,0)),
        ....:    FunctionalDirectedMove([(3/10, 13/40)], (1,0)))
        <FunctionalDirectedMove (1, 0) with domain [(3/10, 7/20), (9/20, 1/2)], range [<Int[3/10, 7/20]>, <Int[9/20, 1/2]>]>
    """
    if A.directed_move != B.directed_move:
        raise ValueError("Cannot merge, moves have different operations")
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

def _plot_component_rectangle(domain, range, **kwds):
    "Default colors are for strip lemma rectangles (densely covered)"
    corners = ((domain[0], range[0]), (domain[1], range[0]), (domain[1], range[1]), (domain[0], range[1]))
    rgbcolor = kwds.pop("rgbcolor", "cyan")
    frame_color = kwds.pop("frame_color", "red")
    kwds.pop("edgecolor", None)
    g = polygon(corners, rgbcolor=rgbcolor, edgecolor=rgbcolor, alpha=0.5, **kwds)
    delete_one_time_plot_kwds(kwds)
    g += polygon(corners, color=frame_color, fill=False, **kwds)
    return g

def plot_covered_component_as_rectangles(component, **kwds):
    g = Graphics()
    for domain in component:
        for range in component:
            g += _plot_component_rectangle(domain, range, **kwds)
            delete_one_time_plot_kwds(kwds)
    return g

def plot_covered_components_as_rectangles(components, **kwds):
    g = Graphics()
    for i, (component, color) in enumerate(zip(components, rainbow(len(components)))):
        g += plot_covered_component_as_rectangles(component,
                                                  rgbcolor=color, frame_color='grey',
                                                  legend_label='Covered component {}'.format(i + 1),
                                                  **kwds)
    return g

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

def partition_overlapping_components(given_component, components):
    overlapping_components = []
    remaining_components = []
    for component in components:
        if list(intersection_of_coho_intervals([given_component, component])):
            # if the two components intersect, the intersection is full-dimensional since intervals are all open intervals.
            overlapping_components.append(component)
        else:
            remaining_components.append(component)
    return overlapping_components, remaining_components

def check_for_strip_lemma_fastpath(m1, m2):
    r"""
    Given two moves `m_1` and `m_2`, return the dense interval by trying to apply the strip lemma.
    Not used, because we want to add covered intervals back to the domains of moves so that L-U is large.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: t1 = 2/15; t2 = sqrt(2)/15; 
        sage: l1 = 4/15; u1 = 3/4 - t1; l2 = 1/3; u2 = 14/15-t2
        sage: m1 = FunctionalDirectedMove([open_interval(l1, u1)], (1, t1))
        sage: m2 = FunctionalDirectedMove([open_interval(l2, u2)], (1, t2))
        sage: check_for_strip_lemma_fastpath(m1, m2)
        [<Int(4/15, 14/15)>]
    """
    if not check_for_strip_lemma_linear_independence(m1, m2):
        return []
    return check_for_strip_lemma_small_translations(m1.intervals(), m2.intervals(), m1[1], m2[1])

def check_for_strip_lemma_linear_independence(m1, m2):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: t1 = 2/15; t2 = sqrt(2)/15; t3 = 3/15;
        sage: l1 = 4/15; u1 = 3/4 - t1; l2 = 1/3; u2 = 14/15-t2
        sage: m1 = FunctionalDirectedMove([open_interval(l1, u1)], (1, t1))
        sage: m2 = FunctionalDirectedMove([open_interval(l2, u2)], (1, t2))
        sage: m3 = FunctionalDirectedMove([open_interval(l2, u2)], (1, t3))
        sage: check_for_strip_lemma_linear_independence(m1, m2)
        True
        sage: check_for_strip_lemma_linear_independence(m1, m3)
        False
    """
    return (m1.sign() == 1 and m2.sign() == 1 and \
            m1[1] >= 0 and m2[1] >= 0 and \
            is_QQ_linearly_independent(m1[1], m2[1]))

def check_for_strip_lemma_small_translations(domain1, domain2, t1, t2):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: t1 = 2/15; t2 = sqrt(2)/15;
        sage: l1 = 4/15; u1 = 3/4 - t1; l2 = 1/3; u2 = 14/15-t2
        sage: m1 = FunctionalDirectedMove([open_interval(l1, u1)], (1, t1))
        sage: m2 = FunctionalDirectedMove([open_interval(l2, u2)], (1, t2))
        sage: check_for_strip_lemma_small_translations(m1.intervals(), m2.intervals(), m1[1], m2[1])
        [<Int(4/15, 14/15)>]
    """
    dense_intervals = []
    for i1 in domain1:
        for i2 in domain2:
            l1, u1 = i1[0], i1[1]
            l2, u2 = i2[0], i2[1]
            L = max(l1, l2)
            LL = min(l1, l2)
            U = min(u1 + t1, u2 + t2)
            UU = max(u1 + t1, u2 + t2)
            if t1 + t2 <= U - L:
                dense_intervals.append(open_interval(LL, UU))
    d = union_of_coho_intervals_minus_union_of_coho_intervals([[interval] for interval in dense_intervals], [])
    return d

# following is used in move_semigroup.py
def interval_including_endpoints_if_continuous(interval, pts_of_discontinuity=True, fdm=None):
    assert pts_of_discontinuity is not None
    if pts_of_discontinuity is True:
        # everywhere discontinuous
        return interval
    # pts of discontinuity is a list
    left_closed = (interval[0] not in pts_of_discontinuity)
    right_closed = (interval[1] not in pts_of_discontinuity)
    if fdm:
        left_closed= left_closed and (fdm.apply_ignoring_domain(interval[0]) not in pts_of_discontinuity)
        right_closed = right_closed and (fdm.apply_ignoring_domain(interval[1]) not in pts_of_discontinuity)
    return closed_or_open_or_halfopen_interval(interval[0], interval[1], left_closed, right_closed)

def plot_points_of_discontinuity_at_borders(pts_of_discontinuity):
    if pts_of_discontinuity is True:
        # Everywhere discontinuous.
        return line(((0, 0), (0, 1), (1, 1)), color='white', linestyle=':', thickness=2, zorder=10)
    if not pts_of_discontinuity:
        return Graphics()
    cts_intervals = union_of_coho_intervals_minus_union_of_coho_intervals([[[0, 1]]],
                                                                          [ [[x]] for x in pts_of_discontinuity ])
    cts_function = FastPiecewise([ (interval, FastLinearFunction(0,0)) for interval in cts_intervals ])
    g = plot_function_at_borders(cts_function, color='black')
    for x in pts_of_discontinuity:
        g += line(((0, x), (1, x)), color='gray', linestyle=':') + line(((x, 0), (x, 1)), color='gray', linestyle=':')
    return g

show_translations_and_reflections_separately = False
show_translations_and_reflections_by_color = True
show_covered_components_as_rectangles = True

class DirectedMoveCompositionCompletion:

    def __init__(self, fdms=[], covered_components=[], proj_add_vert=set(), show_plots=False, plot_background=None, function=None, pts_of_discontinuity=None, show_zero_perturbation=True, show_kwds={}):
        r"""
        Computation of the continuous closure of (the inverse semigroup generated by) a set of moves.

        Input for abstract move ensembles:
         - fdms: list of ``FunctionalDirectedMove`` elements (each is a move with a list of domains)
         - covered_components: list of lists of intervals
         - pts_of_discontinuity: Either ``True'' (everywhere discontinuous), or a list of points
           where discontinuity is allowed.  If not provided (or ``None''), uses the following defaults:
           If ``function`` is provided, uses the allowed discontinuities of its effective perturbations.
           Otherwise, defaults to ``True`` (everywhere discontinuous).

        Additional input for move ensembles coming from a function:
         - function: a ``FastPiecewise`` function
        """
        self.show_plots = show_plots
        self._show_kwds = show_kwds
        self.plot_background = plot_background
        # To show colorful components at borders, and in extend_domain_of_move_by_adding_covered_intervals, need the function. Otherwise set it to default None
        self.function = function
        if function is not None:
            if pts_of_discontinuity is not None:
                raise ValueError("only one of function and pts_of_discontinuity can be provided")
            if function.is_continuous():
                pts_of_discontinuity = []
            elif function.is_two_sided_discontinuous():
                pts_of_discontinuity = True
            elif pts_of_discontinuity is None:
                # function is discontinuous but not two-sided discontinuous at the origin.
                # pts_of_discontinuity is not given. We compute itl
                pts_of_discontinuity = []
                for x in function.end_points():
                    limits = function.limits(x)
                    if not limits[-1]==limits[0]==limits[1]:
                        pts_of_discontinuity.append(x)
        else:
            if pts_of_discontinuity is None:
                pts_of_discontinuity = True
        self.pts_of_discontinuity = pts_of_discontinuity
        self.move_dict = dict()
        self.sym_init_moves = dict() # updated by self.add_backward_moves() in round 0.
        self.covered_components = covered_components
        self._max_component_intervals = 0
        self._max_directed_move_intervals = 0
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
        self._show_zero_perturbation = show_zero_perturbation

    def update_stats(self):
        r"""
        Update ``max_component_intervals`` and ``max_directed_move_intervals``.
        """
        self._max_component_intervals = max(self._max_component_intervals,
                                            self.num_component_intervals())
        self._max_directed_move_intervals = max(self._max_directed_move_intervals,
                                                self.num_directed_move_intervals())

    def num_component_intervals(self):
        return sum(len(c) for c in self.covered_components)

    def num_directed_move_intervals(self):
        return sum(len(m.intervals()) for m in six.itervalues(self.move_dict))

    def max_component_intervals(self):
        return self._max_component_intervals

    def max_directed_move_intervals(self):
        return self._max_directed_move_intervals

    def add_move(self, fdm):
        r"""
        Add a functional directed move to self.
        Merge or restrict functional directed move if necessary.
        """
        reduced_fdm = fdm.reduced_by_components(self.covered_components, self.pts_of_discontinuity)
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
            # sign contradiction point
            if fdm.sign() == -1:
                invariant_point = fdm[1] / 2
                if fdm.can_apply(invariant_point):
                    assert fdm(invariant_point) == invariant_point
                    zero_perturbation_points.add(invariant_point)
        #  return a set, unsorted
        return zero_perturbation_points

    def _plot_directed_moves(self, moves, extra_graphics=None, **kwds):
        g = Graphics()
        if show_covered_components_as_rectangles:
            g += plot_covered_components_as_rectangles(self.covered_components)
        g += plot_directed_moves(moves, **kwds)
        g += plot_points_of_discontinuity_at_borders(self.pts_of_discontinuity)
        if self._show_zero_perturbation:
            zero_perturbation = zero_perturbation_partial_function(self.covered_components,
                                                                   self.generate_zero_perturbation_points())
            if zero_perturbation:
                g += plot_function_at_borders(zero_perturbation, color='magenta', legend_label='fixed perturbation (mod interpol)', thickness=3)
        if self.plot_background:
            g += self.plot_background
        if self.function and self.covered_components:
            g +=  plot_covered_components_at_borders(self.function, self.covered_components, **kwds)
        if extra_graphics:
            g += extra_graphics
        return g

    def plot(self, legend_label='moves', *args, **kwds):
        if show_translations_and_reflections_separately:
            gt = self._plot_directed_moves([ fdm for fdm in self.move_dict.values() if fdm.sign() == +1 ], legend_label='translation moves', **kwds)
            gr = self._plot_directed_moves([ fdm for fdm in self.move_dict.values() if fdm.sign() == -1 ], legend_label='reflection moves', **kwds)
            return graphics_array([gt, gr])
        else:
            return self._plot_directed_moves(list(self.move_dict.values()),
                                             legend_label=legend_label, **kwds)

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
            show_plot(g, self.show_plots, tag, legend_title=title, legend_loc="upper left", object=self, **self._show_kwds)
            logging.info("Plotting... done")

    def extend_components_by_moves(self):
        r"""
        Compose ``self.covered_components`` with the functional directed moves.
        """
        # try extending each component by applying all the moves
        new_components = []
        self.any_change_components = False
        for component in self.covered_components:
            extended_component = component
            for (dm, fdm) in self.move_dict.items():
                overlapped_ints = intersection_of_coho_intervals([extended_component, fdm.intervals()]) # generator
                moved_ints = [[fdm.apply_to_coho_interval(overlapped_int)] for overlapped_int in overlapped_ints]
                new_ints = union_of_coho_intervals_minus_union_of_coho_intervals(moved_ints, [extended_component])
                if new_ints:
                    #can extend the compnent by applying fdm
                    self.any_change_components = True
                    extended_component = union_of_coho_intervals_minus_union_of_coho_intervals([extended_component] + moved_ints, [])
                    extended_component = self.extend_component_by_continuity(extended_component)
            new_components.append(extended_component)
        # got a list of extended components, now merge components
        self.covered_components = reduce_covered_components(new_components)
        # kill moves
        self.reduce_moves_by_components()

    def extend_components_by_strip_lemma(self):
        global crazy_perturbations_warning
        all_move_keys = list(self.move_dict.keys())
        strip_lemma_intervals = []
        for dm_a in all_move_keys:
            for dm_b in all_move_keys:
                a = self.move_dict.get(dm_a, None)
                b = self.move_dict.get(dm_b, None)
                if not a or not b:
                    continue                        # move has been killed
                # check if the strip lemma gives new dense intervals
                d = check_for_strip_lemma_fastpath(a, b)
                if d:
                    if crazy_perturbations_warning:
                        logging.warning("This function is two-sided discontinuous at the origin. Crazy perturbations might exist.")
                    self.any_change_components = True
                    strip_lemma_intervals += d
                    logging.info("New dense move from strip lemma: %s" % d)
                    # merge components with d
                    dense_intervals, self.covered_components = merge_components_with_given_component(d, self.covered_components)
                    dense_intervals = self.extend_component_by_continuity(dense_intervals)
                    self.covered_components.append(dense_intervals)
                    # kill moves
                    self.reduce_moves_by_components()
        if strip_lemma_intervals:
            return sum(plot_covered_component_as_rectangles([interval])
                       for interval in strip_lemma_intervals)
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
                if dm_b not in self.move_dict:
                    continue                        # move has been killed
                c = compose_functional_directed_moves(a, b)
                if c:
                    self.add_move(c)

    def reduce_moves_by_components(self):
        r"""
        Reduce moves with ``self.covered_components``.
        """
        # NOTE: it does not modify any_change_moves even if longer moves are created.
        new_move_dict = dict()
        for (dm, fdm) in self.move_dict.items():
            reduced_fdm = fdm.reduced_by_components(self.covered_components, self.pts_of_discontinuity)
            if reduced_fdm.intervals():
                new_move_dict[dm] = reduced_fdm
        self.move_dict = new_move_dict

    def complete_one_round(self):
        logging.info("Completing %d functional directed moves and %d covered components..." % (len(self.move_dict), len(self.covered_components)))
        if self.any_change_components:
            self.extend_components_by_moves()
        current_dense_move_plot = self.extend_components_by_strip_lemma()
        self.extend_moves_by_composition_of_moves()
        self.update_stats()
        self.num_rounds += 1
        self.maybe_show_plot(current_dense_move_plot)

    def extend_component_by_continuity(self, component):
        assert self.pts_of_discontinuity is not None
        if self.pts_of_discontinuity is not True:
            intervals = [interval_including_endpoints_if_continuous(interval, self.pts_of_discontinuity, None) for interval in component]
            union_intervals = union_of_coho_intervals_minus_union_of_coho_intervals([intervals],[])
            new_component = [open_interval(i[0],i[1]) for i in union_intervals]
            return new_component
        else:
            # discontinuous everywhere
            return component

    def extend_components_by_continuity(self):
        self.covered_components = [self.extend_component_by_continuity(component) for component in self.covered_components]

    def add_backward_moves(self):
        self.num_rounds = 0
        self.extend_components_by_continuity()
        for dm in list(self.any_change_moves):
            # For reflections coming from faces, this is a no-op.
            # But we do this anyway, to be consistent for abstract moves diagram.
            #if dm[0] == -1: continue
            # extend initial moves by continuity
            forward_fdm = self.move_dict[dm].reduced_by_components(self.covered_components, self.pts_of_discontinuity)
            backward_fdm = ~forward_fdm
            self.add_move(backward_fdm)
        # to see that unnecessary discontinuity marks have disappeared,
        # need to set kwds['discontinuity_markers'] = True in FunctionalDirectedMove.plot()
        # and call the following no matter self.any_change_moves is True or False
        self.maybe_show_plot()
        self.any_change_moves = set(self.move_dict.keys())
        self.sym_init_moves = copy(self.move_dict)

    def complete(self, max_num_rounds=None, error_if_max_num_rounds_exceeded=True):
        if self.num_rounds == -1:
            if self.any_change_moves or show_covered_components_as_rectangles:
                # do not show move diagram if there is no moves.
                self.maybe_show_plot()
            self.add_backward_moves()
            self.reduce_moves_by_components()
        self.update_stats()

        while (self.any_change_components or self.any_change_moves) and (max_num_rounds is None or self.num_rounds < max_num_rounds):
            self.complete_one_round()
        if max_num_rounds is not None and self.num_rounds == max_num_rounds:
            if error_if_max_num_rounds_exceeded:
                raise MaximumNumberOfIterationsReached("Reached %d rounds of the completion procedure, found %d directed moves and %d covered components, stopping." % (self.num_rounds, len(self.move_dict), len(self.covered_components)))
            else:
                logging.info("Reached %d rounds of the completion procedure, found %d directed moves and %d covered components, stopping." % (self.num_rounds, len(self.move_dict), len(self.covered_components)))
        else:
            self.is_complete = True
            logging.info("Completion finished.  Found %d directed moves and %d covered components."
                         % (len(self.move_dict), len(self.covered_components)))


    def results(self):
        return list(self.move_dict.values()), self.covered_components
    
    def __repr__(self):
        if self.is_complete:
            return "<DirectedMoveCompositionCompletion (complete) with {} directed moves and {} covered components>".format(len(self.move_dict), len(self.covered_components))
        else:
            return "<DirectedMoveCompositionCompletion (incomplete, {} rounds) with {} directed moves and {} covered components>".format(self.num_rounds, len(self.move_dict), len(self.covered_components))

def directed_move_composition_completion(fdms, covered_components=[], proj_add_vert=set(), show_plots=False, plot_background=None, function=None, pts_of_discontinuity=None, max_num_rounds=None, error_if_max_num_rounds_exceeded=True):
    r"""
    Only used in ``stuff_with_random_irrational_function()``.
    """
    completion = DirectedMoveCompositionCompletion(fdms, covered_components=covered_components, \
                                                   proj_add_vert = proj_add_vert, \
                                                   show_plots=show_plots, \
                                                   plot_background=plot_background, \
                                                   function=function, \
                                                   pts_of_discontinuity=pts_of_discontinuity)
    completion.complete(max_num_rounds=max_num_rounds, error_if_max_num_rounds_exceeded=error_if_max_num_rounds_exceeded)
    return completion.results()

def plot_completion_diagram_background(fn):
    plot_background = plot_function_at_borders(fn, color='black', **ticks_keywords(fn, y_ticks_for_breakpoints=True))
    plot_background += polygon2d([[0,0], [0,1], [1,1], [1,0]], fill=False, color='grey')
    return plot_background

# Global variable ``strategical_covered_components`` to control whether generate_covered_components_strategically() is used in place of generate_covered_components.
strategical_covered_components = False

def generate_covered_components_strategically(fn, show_plots=False, additive_faces=None):
    r"""Return both directly and indirectly covered components of ``fn``.

    Directly covered components are obtained by using the interval lemma
    (convex additivity domain lemma) on two-dimensional maximal additive
    faces of ``fn``.

    Indirectly covered components are obtained using the one-dimensional
    maximal additive faces of ``fn``; this includes using
    limit-additivities.  This is justified by Theorem 3.3 in
    :cite:`koeppe-zhou:crazy-perturbation` (which covers the two-sided
    discontinuous case).

    If ``additive_faces`` is provided, use exactly these faces, assuming that
    they are additive.

    Set ``logging.getLogger().setLevel(logging.DEBUG)`` to see a human-readable 
    proof of covered components.

    Set ``show_plots=True`` to visualize the proof.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: fn = bcdsp_arbitrary_slope(k=4)
        sage: two_slope_fn = two_slope_fill_in_extreme(fn, 1)
        sage: minimal_fn = 9/10 * two_slope_fn + 1/10 * gmic(f=find_f(two_slope_fn))
        sage: minimality_test(minimal_fn)
        True
        sage: len(generate_covered_components(minimal_fn)) # long time
        4
        sage: len(generate_covered_components_strategically(minimal_fn)) # long time
        4
    """
    if additive_faces is None and hasattr(fn, '_strategical_covered_components'):
        return fn._strategical_covered_components
    step = 0
    if show_plots:
        g = plot_2d_diagram(fn, function_color='black', additive_color="grey")
        show_plot(g, show_plots, tag=step , object=fn, show_legend=False, xmin=-0.3, xmax=1.02, ymin=-0.02, ymax=1.3)
    if additive_faces is None:
        use_faces = list(generate_maximal_additive_faces(fn))
        edges = [ face for face in use_faces if face.is_horizontal() or face.is_diagonal() ] # break symmetry ]
    else:
        use_faces = list(additive_faces)
        edges = [ face for face in use_faces if face.is_1D() ] # don't assume both face and x-y-swapped face are present in caller-provided list
    faces = [ face for face in use_faces if face.is_2D() ]
    covered_components = []
    max_size = 1
    while max_size > 0:
        max_size = -1
        max_face = None
        def face_size(face):
            (I, J, K) = face.minimal_triple
            K_mod_1 = interval_mod_1(K)
            newly_covered = union_of_coho_intervals_minus_union_of_coho_intervals([[I], [J], [K_mod_1]], covered_components)
            return sum(interval_length(i) for i in newly_covered)
        def edge_size(edge):
            fdm = edge.functional_directed_move()
            sym_fdm = [fdm]
            if fdm.sign() == 1:
                backward_fdm = FunctionalDirectedMove(fdm.range_intervals(), (1, -fdm[1]))
                sym_fdm.append(backward_fdm)
            moved_intervals = []
            for covered_component in covered_components:
                for fdm in sym_fdm:
                    overlapped_ints = intersection_of_coho_intervals([covered_component, fdm.intervals()]) # generator
                    moved_intervals += [ [fdm.apply_to_coho_interval(overlapped_int)] for overlapped_int in overlapped_ints ]
            newly_covered = union_of_coho_intervals_minus_union_of_coho_intervals(moved_intervals, covered_components)
            return sum(interval_length(i) for i in newly_covered)

        for face in faces:
            size = face_size(face)
            if size > max_size:
                max_size = size
                max_face = face
        for edge in edges:
            size = edge_size(edge)
            if size > max_size:
                max_size = size
                max_face = edge
        if max_size <= 0:
            break

        if max_face.is_2D():
            face = max_face
            step += 1
            faces.remove(face)
            (I, J, K) = face.minimal_triple
            K_mod_1 = interval_mod_1(K)
            component = union_of_coho_intervals_minus_union_of_coho_intervals([[open_interval(* I)], [open_interval(* J)], [open_interval(* K_mod_1)]],[])
            if logging.getLogger().isEnabledFor(logging.DEBUG):
                if hasattr(fn, "_faces_used"):
                    fn._faces_used.append(face)
                logging.debug("Step %s: Consider the 2d additive %s.\n%s is directly covered." % (step, face, component))
            if show_plots:
                if fn.is_continuous():
                    g += face.plot(rgbcolor='red', fill_color='red')
                else:
                    g += face.plot(fill_color='red')
                g += plot_covered_components_at_borders(fn, covered_components=[component])
                show_plot(g + plot_beams_of_one_face(face), show_plots, tag=step, object=fn, show_legend=False, xmin=-0.3, xmax=1.02, ymin=-0.02, ymax=1.3)
            new_component, remaining_components = merge_components_with_given_component(component, covered_components)
            if new_component != component:
                if logging.getLogger().isEnabledFor(logging.DEBUG):
                    logging.debug("We obtain a new covered component %s, with overlapping components merged in." % (new_component))
                if show_plots:
                    show_plot(g, show_plots, tag="{}a".format(step), object=fn, show_legend=False, xmin=-0.3, xmax=1.02, ymin=-0.02, ymax=1.3)
            covered_components = remaining_components + [new_component]

        elif max_face.is_1D():
            edge = max_face
            step += 1
            if hasattr(fn, "_faces_used"):
                fn._faces_used.append(edge)
            if logging.getLogger().isEnabledFor(logging.DEBUG):
                logging.debug("Step %s: Consider the 1d additive %s." % (step, edge))
            fdm = edge.functional_directed_move()
            sym_fdm = [fdm]
            if fdm.sign() == 1:
                backward_fdm = FunctionalDirectedMove(fdm.range_intervals(), (1, -fdm[1]))
                sym_fdm.append(backward_fdm)
            for covered_component in covered_components:
                component = []
                for fdm in sym_fdm:
                    overlapped_ints = list(intersection_of_coho_intervals([covered_component, fdm.intervals()]))
                    moved_intervals = [[fdm.apply_to_coho_interval(overlapped_int)] for overlapped_int in overlapped_ints ]
                    newly_covered = union_of_coho_intervals_minus_union_of_coho_intervals(moved_intervals, covered_components)
                    if newly_covered:
                        if logging.getLogger().isEnabledFor(logging.DEBUG):
                            logging.debug("%s is indirectly covered." % (newly_covered))
                        if show_plots:
                            g += plot_covered_components_at_borders(fn, covered_components=[newly_covered])
                        component = union_of_coho_intervals_minus_union_of_coho_intervals(moved_intervals + [overlapped_ints] + [component], [])
                if component:
                    new_component, remaining_components = merge_components_with_given_component(component, covered_components)
                    if new_component != component and logging.getLogger().isEnabledFor(logging.DEBUG):
                       logging.debug("We obtain a new covered component %s, with overlapping components merged in." % (new_component))
                    covered_components = remaining_components + [new_component]
            if show_plots:
                g += edge.plot(rgbcolor='red')
                show_plot(g + plot_beams_of_one_face(edge), show_plots, tag=step , object=fn, show_legend=False, xmin=-0.3, xmax=1.02, ymin=-0.02, ymax=1.3)

    # There will be no more new covered intervals.
    # covered_components are mutually exclusive.
    # But perhaps merging of components will happen.
    for face in faces:
        (I, J, K) = face.minimal_triple
        K_mod_1 = interval_mod_1(K)
        projections_component = union_of_coho_intervals_minus_union_of_coho_intervals([[open_interval(* I)], [open_interval(* J)], [open_interval(* K_mod_1)]],[])
        overlapping_components, remaining_components = partition_overlapping_components(projections_component, covered_components)
        if len(overlapping_components) > 1:
            new_component = union_of_coho_intervals_minus_union_of_coho_intervals(overlapping_components,[])
            step += 1
            if hasattr(fn, "_faces_used"):
                fn._faces_used.add(face)
            if logging.getLogger().isEnabledFor(logging.DEBUG):
                logging.debug("Step %s: By merging components that overlap with projections of the 2d additive %s, we obtain a larger covered component %s" % (step, face, new_component))
            covered_components = remaining_components + [new_component]

    for edge in edges:
        fdm = edge.functional_directed_move()
        sym_fdm = [fdm]
        if fdm.sign() == 1:
            backward_fdm = FunctionalDirectedMove(fdm.range_intervals(), (1, -fdm[1]))
            sym_fdm.append(backward_fdm)
        for covered_component in covered_components:
            component = []
            for fdm in sym_fdm:
                overlapped_ints = list(intersection_of_coho_intervals([covered_component, fdm.intervals()]))
                moved_intervals = [[fdm.apply_to_coho_interval(overlapped_int)] for overlapped_int in overlapped_ints ]
                component = union_of_coho_intervals_minus_union_of_coho_intervals(moved_intervals + [overlapped_ints] + [component], [])
            if component:
            # component will be in the same (final) covered component as the current covered_component.
            # component is not necessarily new, and it can be a subset of covered_component.
                overlapping_components, remaining_components = partition_overlapping_components(component, covered_components)
                if len(overlapping_components) > 1:
                    new_component = union_of_coho_intervals_minus_union_of_coho_intervals(overlapping_components,[])
                    step += 1
                    if hasattr(fn, "_faces_used"):
                        fn._faces_used.add(edge)
                    if logging.getLogger().isEnabledFor(logging.DEBUG):
                        logging.debug("Step %s: By merging components that are connected by the 1d additive %s, we obtain a larger covered component %s." % (step, edge, new_component))
                    covered_components = remaining_components + [new_component]
    if additive_faces is None: # Don't cache if computation is done with provided list of faces

        fn._strategical_covered_components = covered_components
    return covered_components

def generate_directly_covered_components(fn):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = gj_2_slope(3/5,1/3)
        sage: generate_directly_covered_components(h)
        [[<Int(0, 7/30)>, <Int(11/30, 3/5)>], [<Int(7/30, 11/30)>, <Int(3/5, 1)>]]
    """
    if hasattr(fn, '_directly_covered_components'):
        return fn._directly_covered_components
    covered_components = []
    for face in generate_maximal_additive_faces(fn):
        if face.is_2D():
            covered_components.append(face.covered_component())
    return reduce_covered_components(covered_components)

# alias
generate_directly_covered_intervals = generate_directly_covered_components

def generate_directed_move_composition_completion(fn, show_plots=False, max_num_rounds=None, error_if_max_num_rounds_exceeded=True, plot_background=None):
    completion = getattr(fn, "_completion", None)
    if completion is None:
        functional_directed_moves = generate_functional_directed_moves(fn)
        global strategical_covered_components
        if strategical_covered_components:
            # compute both directly and indirectly covered components.
            covered_components = generate_covered_components_strategically(fn, show_plots=show_plots)
        else:
            covered_components = generate_directly_covered_components(fn)
        proj_add_vert = projections_of_additive_vertices(fn)
        if plot_background is None:
            if show_plots:
                plot_background = plot_completion_diagram_background(fn)
        completion = fn._completion = DirectedMoveCompositionCompletion(functional_directed_moves,
                                                                        covered_components = covered_components,
                                                                        proj_add_vert = proj_add_vert,
                                                                        show_plots=show_plots,
                                                                        plot_background=plot_background,
                                                                        function=fn)
        # To show colorful components at borders and in extend_domain_of_move_by_adding_covered_intervals, need the function. Otherwise set it to default None
        completion.complete(max_num_rounds=max_num_rounds, error_if_max_num_rounds_exceeded=error_if_max_num_rounds_exceeded)
    return completion.results()

def plot_walk_in_completion_diagram(seed, walk_dict):
    g = line([(seed,0), (seed,1)], color="limegreen", legend_label="seed value", linestyle=':')
    kwds = { 'legend_label': "reachable orbit" }
    for x in walk_dict.keys():
        g += line([(0, x), (seed, x)], color="limegreen", linestyle=':', **kwds)
        delete_one_time_plot_kwds(kwds)
    return g

def merge_with_key_xeps_delta(*iterables):
    def xeps_delta(xeps_delta_index):
        xeps, delta, tag = xeps_delta_index
        return xeps, delta
    if six.PY2:
        # merge does not support 'key'
        for i in merge(*iterables):
            yield i
    else:
        # Python >= 3.5 has 'key'
        for i in merge(*iterables, key = xeps_delta):
            yield i

def scan_domains_of_moves(functional_directed_moves):
     scans = [ scan_coho_interval_list(fdm.intervals(), fdm) for fdm in functional_directed_moves ]
     return merge_with_key_xeps_delta(*scans)

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
    return merge_with_key_xeps_delta(*scans)

def find_decomposition_into_intervals_with_same_moves(functional_directed_moves, zero_perturbation_points=set()):
    scan_of_zero_perturbation_points = scan_zero_perturbation_points(zero_perturbation_points)
    scan = merge_with_key_xeps_delta(scan_domains_of_moves(functional_directed_moves),
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
            raise ValueError("Bad scan item")

def find_decomposition_into_stability_intervals_with_completion(fn, show_plots=False, max_num_it=None):
    if hasattr(fn, '_stability_orbits'):
        return
    fn._stability_orbits = []
    long_fdms, covered_components= generate_directed_move_composition_completion(fn, show_plots=show_plots)
    fdms = [fdm.restricting(covered_components) for fdm in long_fdms]
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
    r"""
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
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = equiv7_example_1()
        sage: generate_uncovered_components(h)
        [[<Int(0, 1/4)>, <Int(1/4, 1/2)>]]
        sage: h = minimal_no_covered_interval()
        sage: generate_uncovered_components(h)
        [[<Int(0, 1/4)>, <Int(1/4, 1/2)>], [<Int(1/2, 3/4)>, <Int(3/4, 1)>]]
    """
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
        print("del1 = %s, del2 = %s" % (del1, del2))
        try:
            h = the_irrational_function_t1_t2(del1=del1, del2=del2)
            break
        except ValueError:
            print("... parameters do not describe a function, retrying.")
    dmoves = generate_functional_directed_moves(h)
    fdms, _ = directed_move_composition_completion(dmoves, max_num_rounds=None)
    plot_directed_moves(fdms).show(figsize=40)
    
def merit_index(fn):
    r"""
    Compute Gomory--Johnson's merit index.

    It is defined as twice the area of additive points in the unit square
    (see [RD, Definition 19.28], which refers to [61]).

    EXAMPLES:

    The merit index for GMIC is `2f^2 - 2f + 1` [RD]::

        sage: from cutgeneratingfunctionology.igp import *
        sage: def merit_index_gmic(f):
        ....:     return 2*f^2 - 2*f + 1
        sage: merit_index(gmic(1/4)) == merit_index_gmic(1/4)
        True
        sage: merit_index(gmic(1/2)) == merit_index_gmic(1/2)
        True

    The merit index for ``drlm_2_slope_limit_1_1`` is `4 f^2 - 6 f + 3` [40, p. 164]::

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
    Compute the arithmetic complexity.

    It is defined as the least common denominator of the values fn(i/q) for i=0,1, ... ,q,
    where fn is a piecewise linear function with rational breakpoints in (1/q)Z,
    or a discrete function with its domain contained in (1/q)Z.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
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
        raise ValueError("This is a function has irrational value on (1/%s)Z, so the arithmetic_complexity is not available." % q)
