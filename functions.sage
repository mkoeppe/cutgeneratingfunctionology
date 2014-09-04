import logging

# this is a version of Python's StreamHandler which prints log
# messages to the stream *currently* pointed to by sys.stderr (not the
# one when StreamHandler is set up).  This is useful in a Sage notebook, 
# where every cell has its own set of streams.

class DynamicStdErrStreamHandler(logging.StreamHandler):
    """
    A handler class which writes logging records, appropriately formatted,
    to a stream. Note that this class does not close the stream, as
    sys.stdout or sys.stderr may be used.
    """

    def __init__(self):
        logging.StreamHandler.__init__(self, sys.stderr)

    def flush(self):
        self.stream = sys.stderr
        logging.StreamHandler.flush(self)

    def emit(self, record):
        self.stream = sys.stderr
        logging.StreamHandler.emit(self, record)

#logging.basicConfig(format='%(levelname)s: %(asctime)s %(message)s', level=logging.INFO)

if not logging.getLogger().handlers:
    fmt = logging.Formatter('%(levelname)s: %(asctime)s %(message)s', None)
    hdlr = DynamicStdErrStreamHandler()
    hdlr.setFormatter(fmt)
    logging.getLogger().addHandler(hdlr)
    logging.getLogger().setLevel(logging.INFO)


def logger(func):
    def inner(*args, **kwargs): #1
        print "Arguments to %s were: %s, %s" % (func, args, kwargs)
        result = func(*args, **kwargs) #2
        print "Result is: %s" % (result)
        return result
        
    return inner


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
    ## We now use lambda functions instead of Sage symbolics for plotting, 
    ## as those give strange errors when combined with our RealNumberFieldElement.
    for i in range(1,len(bkpt)):
        p += plot(lambda x: bkpt[i]-x, (x, 0, bkpt[i]), color='grey', **kwd)
        kwd = {}
    for i in range(1,len(bkpt)-1):
        p += plot(lambda x: (1+bkpt[i]-x), (x, bkpt[i], 1), color='grey')
    for i in range(len(bkpt)):
        p += plot(bkpt[i], (0, 1), color='grey')
    y=var('y')
    for i in range(len(bkpt)):
        p += parametric_plot((bkpt[i],y),(y,0,1), color='grey')
    return p

## 
## A lightweight representation of closed bounded intervals, possibly empty or degenerate.
##

def interval_sum(int1, int2):
    """
    Return the sum of two intervals.
    """
    if len(int1) == 1 and len(int2) == 1:
        return [int1[0]+int2[0]]
    elif len(int1) == 2 and len(int2) == 1:
        return [int1[0]+int2[0],int1[1]+int2[0]]
    elif len(int1) == 1 and len(int2) == 2:
        return [int1[0]+int2[0],int1[0]+int2[1]]
    else:    
        return [int1[0]+int2[0],int1[1]+int2[1]]

def interval_intersection(int1, int2):
    """
    Return the intersection of two intervals.

    EXAMPLES::
        sage: interval_intersection([1], [2])
        []
        sage: interval_intersection([1,3], [2,4])
        [2, 3]
        sage: interval_intersection([1,3], [2])
        [2]
        sage: interval_intersection([2], [2,3])
        [2]
        sage: interval_intersection([1,3], [3, 4])
        [3]
        sage: interval_intersection([1,3], [4, 5])
        []
        sage: interval_intersection([], [4, 5])
        []
    """
    if len(int1) == 0 or len(int2) == 0:
        return []
    if len(int1) == 1 and len(int2) == 1:
        if int1[0] == int2[0]:
            return [int1[0]]
        else:
            return []
    elif len(int1) == 2 and len(int2) == 1:
        if int1[0] <= int2[0] <= int1[1]:
            return [int2[0]]
        else:
            return []
    elif len(int1) == 1 and len(int2) == 2:
        if int2[0] <= int1[0] <= int2[1]:
            return [int1[0]]
        else:
            return []    
    else:        
        max0 = max(int1[0],int2[0])
        min1 = min(int1[1],int2[1])
        if max0 > min1:
            return []
        elif max0 == min1:
            return [max0]
        else:
            return [max0,min1]

def interval_empty(interval):
    """
    Determine whether an interval is empty.            
    """
    if len(interval) == 0:
        return True
    else:
        return False

def interval_equal(int1, int2):
    """
    Determine whether two intervals are equal.
    (This ignores whether the intervals are represented as tuples or lists.)
    """
    return tuple(int1) == tuple(int2)

def element_of_int(x,int):
    """
    Determine whether value `x` is inside the interval `int`.

    EXAMPLES::
    sage: element_of_int(1, [])
    False
    sage: element_of_int(1, [1])
    True
    sage: element_of_int(1, [2])
    False
    sage: element_of_int(1, [0,2])
    True
    sage: element_of_int(1, [1,2])
    True
    sage: element_of_int(2, [3,4])
    False
    """
    if len(int) == 0:
        return False
    elif len(int) == 1:
        if x == int[0]:
            return True
        else:
            return False
    elif int[0] <= x <= int[1]:
        return True
    else:
        return False

def interval_to_endpoints(int):
    """
    Convert a possibly degenerate interval to a pair (a,b)
    of its endpoints, suitable for specifying pieces of a `FastPiecewise`.

    EXAMPLES::
        sage: interval_to_endpoints([1])
        (1, 1)
        sage: interval_to_endpoints([1,3])
        (1, 3)
    """
    if len(int) == 0:
        raise ValueError, "An empty interval does not have a pair representation"
    elif len(int) == 1:
        return (int[0], int[0])
    elif len(int) == 2:
        return (int[0], int[1])
    else:
        raise ValueError, "Not an interval: %s" % int

def interval_contained_in_interval(I, J):
    I = interval_to_endpoints(I)
    J = interval_to_endpoints(J)
    return J[0] <= I[0] and I[1] <= J[1]

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

# Remove duplicates in a list.
# FIXME: use some builtin instead.--Matthias
def remove_duplicate(myList):
    if myList:
        myList.sort()
        last = myList[-1]
        for i in range(len(myList)-2, -1, -1):
            if last == myList[i]:
                del myList[i]
            else:
                last = myList[i]

def triples_equal(a, b):
    return interval_equal(a[0], b[0]) and interval_equal(a[1], b[1]) and interval_equal(a[2], b[2])

@cached_function
def generate_maximal_additive_faces(function):
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

            
### Create a new class representing a "face" (which knows its
### vertices, minimal triple, whether it's a translation/reflection,
### etc.; whether it's solid or dense).
class Face:
    def __init__(self, triple, vertices=None, is_known_to_be_minimal=False):
        """
        EXAMPLES::
        sage: load("extreme_functions_in_literature.sage")
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
        trip = self.minimal_triple
        vert = self.vertices
        if self.is_0D():
            return point((trip[0][0], \
                          trip[1][0]), rgbcolor = rgbcolor, size = 30, **kwds)
        elif self.is_horizontal():
            return parametric_plot((y,trip[1][0]),\
                                   (y,trip[0][0], trip[0][1]), rgbcolor = rgbcolor, thickness=2, **kwds)
        elif self.is_vertical():
            return parametric_plot((trip[0][0],y),\
                                   (y,trip[1][0], trip[1][1]), rgbcolor = rgbcolor, thickness=2, **kwds)
        elif self.is_diagonal():
            return parametric_plot((lambda y: y, lambda y: trip[2][0]-y),\
                                   (y,trip[0][0],trip[0][1]), rgbcolor = rgbcolor, thickness=2, **kwds)
        elif self.is_2D():
            ## Sorting is necessary for this example:
            ## plot_2d_diagram(lift(piecewise_function_from_robert_txt_file("/Users/mkoeppe/w/papers/basu-hildebrand-koeppe-papers/reu-2013/Sage_Function/dey-richard-not-extreme.txt"))
            return polygon(convex_vert_list(vert), color=fill_color, **kwds)

    def is_directed_move(self):
        return self.is_1D() #or self.is_0D()
        # FIXME: Do we need additive vertices?
        
    def directed_move_with_domain_and_codomain(self):
        """
        Maps a horizontal edge to a forward translation,
        a vertical edge to a backward translation, 
        a diagonal edge to a reflection.
        """
        trip = self.minimal_triple
        if self.is_horizontal():
            return (1, trip[1][0]), [trip[0]], [trip[2]]
        elif self.is_vertical():
            return (1, -trip[0][0]), [trip[2]], [trip[1]]
        elif self.is_diagonal():
            return (-1, trip[2][0]), [trip[0], trip[1]], [trip[1], trip[0]]
        else:
            raise ValueError, "Face does not correspond to a directed move: %s" % self

    def functional_directed_move(self, intervals=None):
        directed_move, domain, codomain = self.directed_move_with_domain_and_codomain()
        if not intervals == None:
            domain = interval_list_intersection(domain, intervals)
        return FunctionalDirectedMove(domain, directed_move)

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

def plot_2d_diagram(function, show_function = True):
    """
    Return a plot of the 2d complex with shaded faces where delta_pi is 0.        
    To show only a part of it, use 
    `show(plot_2d_diagram(h), xmin=0.25, xmax=0.35, ymin=0.25, ymax=0.35)`
    """
    bkpt = function.end_points()
    faces = generate_maximal_additive_faces(function)

    y = var('y')

    p = plot_2d_complex(function)
    kwds = { 'legend_label': "Additive face" }
    for face in faces:
        p += face.plot(**kwds)
        if 'legend_label' in kwds:
            del kwds['legend_label']
    ### For non-subadditive functions, show the points where delta_pi
    ### is negative.
    nonsubadditive_vertices = generate_nonsubadditive_vertices(function)
    p += point(list(nonsubadditive_vertices), \
                                  color = "red", size = 50, legend_label="Subadditivity violated")
    p += point([ (y,x) for (x,y) in nonsubadditive_vertices ], \
                                  color = "red", size = 50)
    if show_function:
        x = var('x')
        p += parametric_plot((lambda x: x, lambda x: 0.3 * float(function(x)) + 1), \
                                                (x, 0, 1), color='blue', legend_label="Function pi")
        p += parametric_plot((lambda x: - 0.3 * float(function(x)), lambda x: x), \
                                                (x, 0, 1), color='blue')
    return p

# Assume component is sorted.
def merge_within_comp(component, one_point_overlap_suffices=False):   
    for i in range(len(component)-1):
        if component[i][1] > component[i+1][0]  \
           or (one_point_overlap_suffices and component[i][1] == component[i+1][0]):
            component[i+1] = [component[i][0],max(component[i][1],component[i+1][1])]
            component[i] = []
    component_new = []
    for int in component:
        if len(int) == 2 and max(int) <= 1:
            component_new.append(int)
    return component_new


# Assume comp1 and comp2 are sorted.    
def merge_two_comp(comp1,comp2, one_point_overlap_suffices=False):
    temp = []
    i = 0
    j = 0
    while i < len(comp1) and j < len(comp2):
        if comp1[i][0] < comp2[j][0]:
            temp.append(comp1[i])
            i = i+1
        else:
            temp.append(comp2[j])
            j = j+1
    if i == len(comp1):
        temp = temp + comp2[j:len(comp2)]
    else:
        temp = temp + comp1[i:len(comp1)]
    temp = merge_within_comp(temp, one_point_overlap_suffices=one_point_overlap_suffices)
    return temp
            

def partial_overlap(interval,component):
    """
    Return a list of the intersections of the interiors 
    of `interval` and the intervals in `component`.

    EXAMPLES::
    sage: partial_overlap([2,3], [[1,2], [3,5]])
    []
    sage: partial_overlap([2,6], [[1,3], [5,7], [7,9]])
    [[2, 3], [5, 6]]
    """
    overlap = []
    for int1 in component:
        overlapped_int = interval_intersection(interval,int1)
        if len(overlapped_int) == 2:
            overlap.append(overlapped_int)
    return overlap


def remove_empty_comp(comps):
    """
    Return a new list that includes all non-empty lists of `comps`.

    EXAMPLES::
    sage: remove_empty_comp([[[1,2]], [], [[3,4],[5,6]]])
    [[[1, 2]], [[3, 4], [5, 6]]]
    """
    temp = []
    for int in comps:
        if len(int) > 0:
            temp.append(int)
    return temp
    

def partial_edge_merge(comps, partial_overlap_intervals, ijk, ijk2, intervals, i, IJK):
    """
    Modifies the list `comps`.
    Returns whether any change occurred.
    """
    any_change = False
    for int1 in partial_overlap_intervals:
        front = int1[0] - intervals[ijk][0]
        back = intervals[ijk][1] - int1[1]
        
        # If it is not the pair I and J, then the action is a translation.
        if IJK != [0,1]:
            other = [intervals[ijk2][0]+front, intervals[ijk2][1]-back]
        # I and J implies a reflection
        else:
            other = [intervals[ijk2][0]+back, intervals[ijk2][1]-front]
        other = interval_mod_1(other)
        #print "other: ", other
            
        overlapped_component_indices = []
        i_included = False
        all_other_overlaps = []
        for k in range(len(comps)):
            other_overlap = partial_overlap(other,comps[k])
            #print "other_overlap:", other_overlap
            if other_overlap:
                #print "overlap with component", k, "is: ", other_overlap
                all_other_overlaps = merge_two_comp(all_other_overlaps, other_overlap)
                if k < i:
                    overlapped_component_indices.append(k)
                elif k > i and i_included == False:
                    overlapped_component_indices.append(i)
                    overlapped_component_indices.append(k)
                    i_included = True
                else:
                    overlapped_component_indices.append(k)
        if overlapped_component_indices == [i] :
            ## Only overlap within component i.
            # print "Self-overlap only"
            if (partial_overlap(other, comps[i]) == [other]):
                pass
            else:
                comps[overlapped_component_indices[-1]] = merge_two_comp(comps[overlapped_component_indices[-1]], [other])
                any_change = True
        elif len(overlapped_component_indices) > 0:
            ## Overlap with some other components; this will cause some merging.
            #print "Have overlapped components: ", overlapped_component_indices, "with ", i
            comps[overlapped_component_indices[-1]] = merge_two_comp(comps[overlapped_component_indices[-1]], [other])
            for j in range(len(overlapped_component_indices)-1):
                comps[overlapped_component_indices[j+1]] =  merge_two_comp(comps[overlapped_component_indices[j]],\
                     comps[overlapped_component_indices[j+1]])
                comps[overlapped_component_indices[j]] = []
            any_change = True

        # previous non-covered:
        #print "other: ", other, "all_other_overlaps: ", all_other_overlaps
        noncovered_overlap = interval_minus_union_of_intervals(other, all_other_overlaps)
        if noncovered_overlap:
            # print "Previously non-covered: ", uncovered_intervals_from_covered_intervals(comps)
            # print "Newly covered: ", noncovered_overlap
            any_change = True
            comps[i] = merge_two_comp(comps[i], noncovered_overlap)
            # print "Now non-covered: ", uncovered_intervals_from_covered_intervals(comps)
    return any_change
                  

def edge_merge(comps,intervals,IJK):
    #print "edge_merge(%s,%s,%s)" % (comps, intervals, IJK)
    any_change = False
    for i in range(len(comps)): 
        partial_overlap_intervals = partial_overlap(intervals[0],comps[i])
        # If there is overlapping...
        if len(partial_overlap_intervals) > 0:
            if partial_edge_merge(comps, partial_overlap_intervals, 0, 1, intervals, i, IJK):
                any_change = True
        # Repeat the same procedure for the other interval.
        partial_overlap_intervals = partial_overlap(intervals[1],comps[i])
        if len(partial_overlap_intervals) > 0:
            if partial_edge_merge(comps, partial_overlap_intervals, 1, 0, intervals, i, IJK):
                any_change = True
    return any_change
    
# Assume the lists of intervals are sorted.                
def find_interior_intersection(list1, list2):
    """
    Tests whether `list1` and `list2` contain a pair of intervals
    whose interiors intersect.

    Assumes both lists are sorted.
    
    EXAMPLES:
    sage: find_interior_intersection([[1, 2], [3, 4]], [[2, 3], [4, 5]])
    False
    sage: find_interior_intersection([[1, 2], [3, 5]], [[2, 4]])
    True
    """
    i=0
    j=0
    while i < len(list1) and j < len(list2):
        if len(interval_intersection(list1[i], list2[j])) == 2:
            return True
        else:
            if list1[i][0] < list2[j][0]:
                i = i + 1
            else:
                j = j + 1
    return False

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

@cached_function
def generate_covered_intervals(function):
    logging.info("Computing covered intervals...")
    faces = generate_maximal_additive_faces(function)

    covered_intervals = []      
    for face in faces:
        if face.is_2D():
            component = []
            for int1 in face.minimal_triple:
                component.append(interval_mod_1(int1))
            component.sort()
            component = merge_within_comp(component)
            covered_intervals.append(component)
            
    remove_duplicate(covered_intervals)
    
    #show(plot_covered_intervals(function, covered_intervals), xmax=1.5)

    for i in range(len(covered_intervals)):
        for j in range(i+1, len(covered_intervals)):
            if find_interior_intersection(covered_intervals[i], covered_intervals[j]):
                covered_intervals[j] = merge_two_comp(covered_intervals[i],covered_intervals[j])
                covered_intervals[i] = []
                    
    covered_intervals = remove_empty_comp(covered_intervals)
    # debugging plot:
    # show(plot_covered_intervals(function, covered_intervals), \
    #      legend_fancybox=True, \
    #      legend_title="Directly covered, merged", \
    #      legend_loc=2) # legend in upper left

    edges = [ face.minimal_triple for face in faces if face.is_1D()]

    any_change = True
    ## FIXME: Here we saturate the covered interval components
    ## with the edge relations.  There should be a smarter way
    ## to avoid this while loop.  Probably by keeping track 
    ## of a set of non-covered components (connected by edges).
    ## --Matthias
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

    covered_intervals = remove_empty_comp(covered_intervals)
    logging.info("Computing covered intervals... done")
    return covered_intervals

def interval_minus_union_of_intervals(interval, remove_list):
    """Compute a list of intervals that represent the
    set difference of `interval` and the union of the 
    intervals in `remove_list`.

    Assumes `remove_list` is sorted (and pairwise essentially
    disjoint), and returns a sorted list.

    EXAMPLES::
    sage: interval_minus_union_of_intervals([0, 10], [[-1, 0], [2, 3], [9,11]]) 
    [[0, 2], [3, 9]]
    sage: interval_minus_union_of_intervals([0, 10], [[-1, 0], [2, 3]]) 
    [[0, 2], [3, 10]]
    sage: interval_minus_union_of_intervals([0, 10], [[-1, 0], [2, 3], [9,11], [13, 17]])
    [[0, 2], [3, 9]]
    """
    scan = scan_union_of_coho_intervals_minus_union_of_coho_intervals([[interval]], [remove_list])
    return list(proper_interval_list_from_scan(scan))

def uncovered_intervals_from_covered_intervals(covered_intervals):
    """Compute a list of uncovered intervals, given the list of components
    of covered intervals.

    EXAMPLES::
    sage: uncovered_intervals_from_covered_intervals([[[10/17, 11/17]], [[5/17, 6/17], [7/17, 8/17]]])
    [[0, 5/17], [6/17, 7/17], [8/17, 10/17], [11/17, 1]]
    """

    covered = reduce(merge_two_comp, covered_intervals)
    return interval_minus_union_of_intervals([0,1], covered)

@cached_function
def generate_uncovered_intervals(function):
    """
    Compute a sorted list of uncovered intervals.
    """
    covered_intervals = generate_covered_intervals(function)
    return uncovered_intervals_from_covered_intervals(covered_intervals)

def ticks_keywords(function, y_ticks_for_breakpoints=False):
    """
    Compute `plot` keywords for displaying the ticks.
    """
    xticks = function.end_points()
    xtick_formatter = [ "$%s$" % latex(x) for x in xticks ]
    #xtick_formatter = 'latex'  # would not show rationals as fractions
    ytick_formatter = None
    if y_ticks_for_breakpoints:
        yticks = xticks
        ytick_formatter = xtick_formatter
    else:
        #yticks = 1/5
        yticks = uniq([ function(x) for x in function.end_points() ])
        ytick_formatter = [ "$%s$" % latex(y) for y in yticks ]
    ## FIXME: Can we influence ticks placement as well so that labels don't overlap?
    ## or maybe rotate labels 90 degrees?
    return {'ticks': [xticks, yticks], \

            'gridlines': True, \
            'tick_formatter': [xtick_formatter, ytick_formatter]}

def plot_covered_intervals(function, covered_intervals=None, **plot_kwds):
    """
    Return a plot of the covered and uncovered intervals of `function`.
    """
    if covered_intervals is None:
        covered_intervals = generate_covered_intervals(function)
        uncovered_intervals = generate_uncovered_intervals(function)
    else:
        uncovered_intervals = uncovered_intervals_from_covered_intervals(covered_intervals)
    # Plot the function with different colors.
    # Each component has a unique color.
    # The uncovered intervals is by default plotted in black.
    colors = rainbow(len(covered_intervals))
    graph = Graphics()
    kwds = copy(plot_kwds)
    kwds.update(ticks_keywords(function))
    if uncovered_intervals:
        graph += plot(function, [0,1],
                      color = "black", legend_label="not covered", **kwds)
        kwds = {}
    elif not function.is_continuous_defined(): # to plot the discontinuity markers
        graph += plot(function, [0,1], color = "black", **kwds)
        kwds = {}
    for i, component in enumerate(covered_intervals):
        kwds.update({'legend_label': "covered component %s" % (i+1)})
        for interval in component:
            graph += plot(function.which_function((interval[0] + interval[1])/2), interval, color=colors[i], **kwds)
            if 'legend_label' in kwds:
                del kwds['legend_label']
            if 'ticks' in kwds:
                del kwds['ticks']
            if 'tick_formatter' in kwds:
                del kwds['tick_formatter']
    return graph

### Minimality check.

def subadditivity_check(fn):
    """
    Check if a function is subadditive.
    Could take quite a while. (O(n^2))
    """
    result = True
    for (x, y) in generate_nonsubadditive_vertices(fn):
        logging.info("For x = %s, y = %s, x+y = %s" % (x, y, x+y))
        logging.info("    Delta pi(x,y) = %s + %s - %s = %s < 0" % 
                     (fn(fractional(x)), fn(fractional(y)), fn(fractional(x+y)), delta_pi(fn,x,y)))
        result = False
    if result:
        logging.info("pi is subadditive.")
    else:
        logging.info("Thus pi is not subadditive.")
    return result

def symmetric_test(fn, f):
    result = True
    if fn(f) != 1:
        logging.info('pi(f) is not equal to 1.')
        result = False
    else:
        bkpt = fn.end_points()
        for i in bkpt:
            if i < f:
                j = f - i
            else:
                j = 1 + f - i
            if delta_pi(fn, i, j) != 0:
                logging.info('For x = %s, y = %s' % (i,j))
                logging.info('    Delta pi is equal to %s, not equal to 0' % delta_pi(fn, i, j))
                result = False
    if result:
        logging.info('pi is symmetric.')
    else:
        logging.info('Thus pi is not symmetric.')
    return result

@cached_function
def find_f(fn, no_error_if_not_minimal_anyway=False):
    """
    Find the value of `f' for the given function `fn'.
    """
    f = None
    for x in fn.end_points():
        if fn(x) > 1 or fn(x) < 0: 
            if no_error_if_not_minimal_anyway:
                logging.info('pi is not minimal because it does not stay in the range of [0, 1].')
                return None
            raise ValueError, "The given function does not stay in the range of [0, 1], so cannot determine f.  Provide parameter f to minimality_test or extremality_test."
    for x in fn.end_points():
        if fn(x) == 1:
            if f:
                logging.warn("The given function has more than one breakpoint where the function takes the value 1; using f = %s.  Provide parameter f to minimality_test or extremality_test if you want a different f." % f)
                return f
            else:
                f = x
    if f:
        return f
    if no_error_if_not_minimal_anyway:
        logging.info('pi is not minimal because it has no breakpoint where the function takes value 1.')
        return None
    raise ValueError, "The given function has no breakpoint where the function takes value 1, so cannot determine f.  Provide parameter f to minimality_test or extremality_test."

def minimality_test(fn, show_plots=False, f=None):
    """
    Check if function `fn` is minimal with respect to the given `f`.
    If `f` is not provided, uses the one found by `find_f`.
    Return a boolean.

    EXAMPLES::
        sage: logging.disable(logging.INFO)
        sage: minimality_test(piecewise_function_from_breakpoints_and_values([0,1/5,4/5,1],[0,1/2,1,0]))
        False
        sage: minimality_test(piecewise_function_from_breakpoints_and_values([0,1/2,1], [0,2,0]))
        False
    """
    if f==None:
        f = find_f(fn, no_error_if_not_minimal_anyway=True)
        if f==None:
            return False
    if fn(0) != 0:
        logging.info('pi is not minimal because pi(0) is not equal to 0.')
        return False
    else:
        logging.info('pi(0) = 0')
        if show_plots:
            logging.info("Plotting 2d diagram...")
            show_plot(plot_2d_diagram(fn), show_plots, tag='2d_diagram')
            logging.info("Plotting 2d diagram... done")
        if subadditivity_check(fn) and symmetric_test(fn, f):
            logging.info('Thus pi is minimal.')
            return True
        else:
            logging.info('Thus pi is not minimal.')
            return False

global default_field 
default_field = QQ

## Lambda functions are the fastest to evaluate, but we cannot add them,
## which may be inconvenient.
def very_fast_linear_function(slope, intercept, field=default_field):
    """
    Return a linear function.
    """
    # Ignore field in this implementation.
    return lambda x: slope * x + intercept

# Linear univariate polynomials are 10 times slower than lambda functions to evaluate,
# but still 10 times faster to evaluate than symbolic expressions.
# Note that this implementation is NOT compatible with symbolic expressions.
# For example slope=0, intercept=sqrt(2) leads to the result
# being just a symbolic expression, which is not callable.
def fast_addable_linear_function(slope, intercept, field=default_field):
     RK = PolynomialRing(field, 'x')
     x = RK.0
     return slope * x + intercept

from sage.functions.piecewise import PiecewisePolynomial
from bisect import bisect_left

## FIXME: Its __name__ is "Fast..." but nobody so far has timed
## its performance against the other options. --Matthias
class FastLinearFunction :

    def __init__(self, slope, intercept, field=default_field):
        self._slope = slope
        self._intercept = intercept

    def __call__(self, x):
        if type(x) == float:
            # FIXME: There must be a better way.
            return float(self._slope) * x + float(self._intercept)
        else:
            return self._slope * x + self._intercept

    def __float__(self):
        return self

    def __add__(self, other):
        return FastLinearFunction(self._slope + other._slope,
                                  self._intercept + other._intercept)

    def __mul__(self, other):
        # scalar multiplication
        return FastLinearFunction(self._slope * other,
                                  self._intercept * other)


    def __neg__(self):
        return FastLinearFunction(-self._slope,
                                  -self._intercept)

    __rmul__ = __mul__

    def __eq__(self, other):
        if not isinstance(other, FastLinearFunction):
            return False
        return self._slope == other._slope and self._intercept == other._intercept

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        try:
            return '<FastLinearFunction ' + sage.misc.misc.repr_lincomb([('x', self._slope), (1, self._intercept)], strip_one = True) + '>'
        except TypeError:
            return '<FastLinearFunction (%s)*x + (%s)>' % (self._slope, self._intercept)

    ## FIXME: To be continued.

fast_linear_function = FastLinearFunction

class FastPiecewise (PiecewisePolynomial):
    """
    Returns a piecewise function from a list of (interval, function)
    pairs.

    Uses binary search to allow for faster function evaluations
    than the standard class PiecewisePolynomial.
    """
    def __init__(self, list_of_pairs, var=None, periodic_extendable=True, merge=True):
        """
        EXAMPLES:
        sage: h = FastPiecewise([[(3/10, 15/40), FastLinearFunction(1, 0)], [(13/40, 14/40), FastLinearFunction(1, 0)]], merge=True)
        sage: len(h.intervals())
        1
        sage: h.intervals()[0][0], h.intervals()[0][1]
        (3/10, 15/40)
        """
        # Sort intervals according to their left endpoints; In case of equality, place single point before interval. 
        # This setting would be helpful in plotting discontinuous functions   
        list_of_pairs = sorted(list_of_pairs, key = coho_interval_left_endpoint_with_epsilon)
        if merge:
            merged_list_of_pairs = []
            intervals_to_scan = []
            singleton = None
            common_f = None
            for (i, f) in list_of_pairs:
                if common_f == f:
                    intervals_to_scan.append(i)
                    singleton = None
                elif common_f != None and singleton != None and common_f(singleton) == f(singleton):
                    intervals_to_scan.append(i)
                    singleton = None
                    common_f = f
                elif i[0] == i[1] and common_f != None and common_f(i[0]) == f(i[0]):
                    intervals_to_scan.append(i)
                else:
                    merged_intervals = coho_interval_list_from_scan(scan_coho_interval_list(intervals_to_scan))
                    for merged_i in merged_intervals:
                        if merged_i.left_closed and merged_i.right_closed and merged_i[0] != merged_i[1]:
                            merged_interval = (merged_i[0], merged_i[1])
                        else:
                            merged_interval = merged_i
                        merged_list_of_pairs.append((merged_interval, common_f))
                    intervals_to_scan = [i]
                    if i[0] == i[1]:
                        singleton = i[0]
                    else:
                        singleton = None
                    common_f = f
            merged_intervals = union_of_coho_intervals_minus_union_of_coho_intervals([[interval] for interval in intervals_to_scan], [])
            for merged_i in merged_intervals:
                if merged_i.left_closed and merged_i.right_closed and merged_i[0] != merged_i[1]:
                    merged_interval = (merged_i[0], merged_i[1])
                else:
                    merged_interval = merged_i
                merged_list_of_pairs.append((merged_interval, common_f))
            list_of_pairs = merged_list_of_pairs
            
        PiecewisePolynomial.__init__(self, list_of_pairs, var)

        intervals = self._intervals
        functions = self._functions
        # end_points are distinct.
        end_points = []
        # ith_at_end_points records in which interval the end_point first appears as a left_end or right_end.
        ith_at_end_points = []
        # record the value at each end_point, value==None if end_point is not in the domain.
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
                if left_value != None:
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
            elif right_value != None:
                values_at_end_points[-1] = right_value
        if periodic_extendable:
            if values_at_end_points[0] != values_at_end_points[-1]:
                logging.warn("Function is actually not periodically extendable.")
                periodic_extendable = False
            else:
                limits_at_end_points[0][-1] = limits_at_end_points[-1][-1]
                limits_at_end_points[-1][1] = limits_at_end_points[0][1]
        self._end_points = end_points
        self._ith_at_end_points = ith_at_end_points
        self._values_at_end_points = values_at_end_points
        self._limits_at_end_points = limits_at_end_points
        self._periodic_extendable = periodic_extendable

    # The following makes this class hashable and thus enables caching
    # of the above functions; but we must promise not to modify the
    # contents of the instance.
    def __hash__(self):
        return id(self)

    def end_points(self):
        """
        Returns a list of all interval endpoints for this function.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 2
            sage: f3(x) = 1-x
            sage: f4(x) = x^2-5
            sage: f = FastPiecewise([[open_interval(0,1),f1],[singleton_interval(1),f2],[open_interval(1,2),f3],[(2,3),f4]], merge=False)
            sage: f.end_points()
            [0, 1, 2, 3]
            sage: f = FastPiecewise([[open_interval(0,1),f1],[open_interval(2,3),f3]], merge=False)
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
            ...                      [(9,10),f7]], merge=False)
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
            ...                      [(9,10),f7]], periodic_extendable= False, merge=False)
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
            ...                      [closed_interval(3,10),f4]], merge=False)
            sage: f.which_function(0.5) is f1
            True
            sage: f.which_function(1) is f2
            True
            sage: f.which_function(5/2) is f3
            True
            sage: f.which_function(3) is f4
            True
            sage: f = FastPiecewise([[open_interval(0,1),f1],
            ...                      [right_open_interval(2,3),f3]], merge=False)
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
            if not self._values_at_end_points[i] == None:
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

    @cached_method
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
            ...                      [closed_interval(3,10),f4]], merge=False)
            sage: f(0.5)
            1
            sage: f(1)
            0
            sage: f(5/2)
            e^(5/2)
            sage: f(3)
            sin(6)
            sage: f = FastPiecewise([[open_interval(0,1),f1],
            ...                      [right_open_interval(2,3),f3]], merge=False)
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
        # Remember that intervals are sorted according to their left endpoints; singleton has priority.
        endpts = self.end_points()
        ith = self._ith_at_end_points
        i = bisect_left(endpts, x0)
        if i >= len(endpts):
            raise ValueError,"Value not defined at point %s, outside of domain." % x0
        if x0 == endpts[i]:
            if not self._values_at_end_points[i] == None:
                return self._values_at_end_points[i]
            else:
                raise ValueError,"Value not defined at point %s, outside of domain." % x0
        if i == 0:
            raise ValueError,"Value not defined at point %s, outside of domain." % x0
        if is_pt_in_interval(self._intervals[ith[i]],x0):
            return self.functions()[ith[i]](x0)
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
            ...                      [(9,10),f7]], periodic_extendable=False, merge=False)
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
            ...                      [(9,10),f7]], periodic_extendable=False, merge=False)
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
        if result == None:
            raise ValueError,"Value not defined at point %s%s, outside of domain." % (x0, print_sign(epsilon))
        return result

    def which_function_on_interval(self, interval):
        x = (interval[0] + interval[1]) / 2
        # FIXME: This should check that the given `interval` is contained in the defining interval!
        # This could be implemented by refactoring which_function using new function which_function_index.
        return self.which_function(x)

    def __add__(self,other):
        """
        In contrast to PiecewisePolynomial.__add__, this does not do zero extension of domains.
        Rather, the result is only defined on the intersection of the domains.

        EXAMPLES::
        sage: f = FastPiecewise([[singleton_interval(1), FastLinearFunction(0,17)]], merge=False)
        sage: g = FastPiecewise([[[0,2], FastLinearFunction(0,2)]], merge=False)
        sage: (f+g).list()
        [[<Int{1}>, <FastLinearFunction 19>]]
        sage: h = FastPiecewise([[open_interval(1,3), FastLinearFunction(0,3)]], merge=False)
        sage: (g+h).list()
        [[<Int(1, 2]>, <FastLinearFunction 5>]]
        sage: j = FastPiecewise([[open_interval(0,1), FastLinearFunction(0,1)], [[1, 3], FastLinearFunction(0, 5)]], merge=False)
        sage: (g+j).list()
        [[<Int(0, 1)>, <FastLinearFunction 3>], [(1, 2), <FastLinearFunction 7>]]
        """
        intervals = intersection_of_coho_intervals([self.intervals(), other.intervals()])
        return FastPiecewise([ (interval, self.which_function_on_interval(interval) + other.which_function_on_interval(interval))
                               for interval in intervals ], merge=True)

    def __neg__(self):
        return FastPiecewise([[interval, -f] for interval,f in self.list()], merge=True)
        
    def __mul__(self,other):
        """In contrast to PiecewisePolynomial.__mul__, this does not do zero extension of domains.
        Rather, the result is only defined on the intersection of the domains."""
        if not isinstance(other, FastPiecewise):
            # assume scalar multiplication
            return FastPiecewise([[interval, other*f] for interval,f in self.list()])
        else:
            intervals = intersection_of_coho_intervals([self.intervals(), other.intervals()])
            return FastPiecewise([ (interval, self.which_function_on_interval(interval) * other.which_function_on_interval(interval))
                                   for interval in intervals ], merge=True)

    __rmul__ = __mul__

    def __sub__(self, other):
        return self + (-other)

    ## Following just fixes a bug in the plot method in piecewise.py
    ## (see doctests below).  Also adds plotting of single points.
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
        this has been fixed:

            sage: q = f.plot(xmin=0, xmax=3)
            sage: q = plot(f, xmin=0, xmax=3)
            sage: q = plot(f, 0, 3)
            sage: q = plot(f, 0, 3, color='red')
        
        The implementation should crop according to the given xmin, xmax.

            sage: q = plot(f, 1/2, 3)
            sage: q = plot(f, 1, 2)
            sage: q = plot(f, 2, 3)
        
        Also the following plot syntax should be accepted.

            sage: q = plot(f, [2, 3])

        """
        from sage.plot.all import plot, Graphics

        g = Graphics()
        if not 'rgbcolor' in kwds:
            color = 'blue'
        else:
            color = kwds['rgbcolor']
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
            # Handle open/half-open intervals here
            if (a < b) or (a == b) and (left_closed) and (right_closed):
                if not (last_closed or last_end_point == [a, f(a)] and left_closed):
                    # plot last open right endpoint
                    g += point(last_end_point, color=color, pointsize=23)
                    g += point(last_end_point, rgbcolor='white', pointsize=10)
                if last_closed and last_end_point != [] and last_end_point != [a, f(a)] and not left_closed:
                    # plot last closed right endpoint
                    g += point(last_end_point, color=color, pointsize=23)
                if not (left_closed or last_end_point == [a, f(a)] and last_closed):
                    # plot current open left endpoint
                    g += point([a, f(a)], color=color, pointsize=23)
                    g += point([a, f(a)], rgbcolor='white', pointsize=10)
                if left_closed and last_end_point != [] and last_end_point != [a, f(a)] and not last_closed:
                    # plot current closed left endpoint
                    g += point([a, f(a)], color=color, pointsize=23)
                last_closed = right_closed
                last_end_point = [b, f(b)]
            if a < b:
                # We do not plot anything if a==b because
                # otherwise plot complains that
                # "start point and endpoint must be different"
                g += plot(f, *args, xmin=a, xmax=b, zorder=-1, **kwds)
                # If it's the first piece, pass all arguments. Otherwise,
                # filter out 'legend_label' so that we don't add each
                # piece to the legend separately (trac #12651).
                if 'legend_label' in kwds:
                    del kwds['legend_label']
            elif (a == b) and (left_closed) and (right_closed):
                g += point([a, f(a)], color=color, pointsize=23)
        # plot open rightmost endpoint. minimal functions don't need this.
        if not last_closed:
            g += point(last_end_point, color=color,pointsize=23)
            g += point(last_end_point, rgbcolor='white', pointsize=10)
        return g

    def is_continuous_defined(self, xmin=0, xmax=1):
        """
        return True if self is defined on [xmin,xmax] and is continuous on [xmin,xmax]
        """
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
            if a < xmin:
                a = xmin
                left_closed = True
            if b > xmax:
                b = xmax
                right_closed = True
            if (a < b) or (a == b) and (left_closed) and (right_closed):
                if last_end_point == [] and not a == xmin:
                    return False
                if not (last_closed or last_end_point == [a, f(a)] and left_closed):
                    return False
                if not (left_closed or last_end_point == [a, f(a)] and last_closed):
                    return False
                last_closed = right_closed
                last_end_point = [b, f(b)]
        if not (last_closed and last_end_point[0] == xmax):
            return False
        return True

    def __repr__(self):
        rep = "<FastPiecewise with %s parts, " % len(self._functions)
        for interval, function in itertools.izip(self._intervals, self._functions):
            rep += "\n " + repr(interval) + "\t" + repr(function) \
                   + "\t values: " + repr([function(interval[0]), function(interval[1])])
        rep += ">"
        return rep
        
def print_sign(epsilon):
    if epsilon > 0:
        return "+"
    elif epsilon < 0:
        return "-"
    else:
        return ""

def is_pt_in_interval(i, x0):
    """
    retrun whether the point x0 is contained in the (ordinary or coho) interval i.
    """
    if len(i) == 2:
        return bool(i[0] <= x0 <= i[1])
    else:  
        if i.left_closed and i.right_closed:
            return bool(i.a <= x0 <= i.b)
        if i.left_closed and not i.right_closed:
            return bool(i.a <= x0 < i.b)
        if not i.left_closed and i.right_closed:
            return bool(i.a < x0 <= i.b)
        if not i.left_closed and not i.right_closed:
            return bool(i.a < x0 < i.b)

from sage.rings.number_field.number_field_element_quadratic import NumberFieldElement_quadratic

### def can_coerce(coercer, values):
###     try:
###         coerced_values = [ coercer(value) for value in values ]
###         return True
###     except ValueError:
###         return False
###     except TypeError:
###         return False

### def fast_field(values):
###     "Return a field and a coercing function."
###     ## If everything is rational, we are happy.
###     field, coercer = QQ, QQ
###     if can_coerce(coercer, values):
###         return field, coercer
###     ## ## See if we can make do with a quadratic number field
###     ## radicant_wildcard = SR.wild(0)
###     ## radicals = set()
###     ## for value in values:
###     ##     if parent(value) == SR:
###     ##         radicals |= set(value.find(sqrt(radicant_wildcard)))
###     ## logging.info("Set of radicals: %s" % (radicals,))
###     ## if len(radicals) == 1:
###     ##     ## 
###     ##     radicant = radicals[0].op[0]
###     ##     field.<root> = QuadraticField(radicant, name='sqrt%s' % radicant)
###     ##     def coercer(x, field=field, root=root, radicant=radicant):
###     ##         if parent(x) == SR:
###     ##             return x.subs(sqrt(radicant)==root) 
###     ##             ## Does not work because we get immediate
###     ##             ## back-substitution, when sqrt2 is inserted in the
###     ##             ## symbolic ring.
###     ##         else
###     ##     if can_coerce(coercer, values):
###     ##         return field.values
###     try:
###         field, field_values, morphism = number_field_elements_from_algebraics(values)
###     except ValueError:
###         pass
###     except TypeError:
###         pass
    
import sage.rings.number_field.number_field
from sage.rings.number_field.number_field import NumberField_absolute, NumberField_quadratic

import sage.rings.number_field.number_field_element
from sage.rings.number_field.number_field_element import NumberFieldElement_absolute

class RealNumberFieldElement(NumberFieldElement_absolute):

    ## def __init__(self, parent, f):
    ##     NumberFieldElement_absolute.__init__(self, parent, f)
    ##     self._embedding = parent.coerce_embedding()

    #@logger

    def embedded(self):
        e = getattr(self, '_embedded', None)
        if e is None:
            parent = self.parent()
            embedding = parent.coerce_embedding()
            self._embedded = e = embedding(self)
        return e

    def __cmp__(left, right):
    #     return (left)._cmp(right)
    # def _cmp(left, right):
        # print "cmp", left, "and", right
        #return sign(left - right) ## Whoa, sign is defined in "generalized.py"

        if NumberFieldElement_absolute.__cmp__(left, right) == 0:
            return 0
        result = cmp(left.embedded(), right.embedded())
        if result == 0:
            raise UnimplementedError, "Precision of real interval field not sufficient to continue"
        return result

    def __abs__(self):
        if self.sign() >= 0:
            return self
        else:
            return -self

    def sign(self):
        parent = self.parent()
        return cmp(self, parent._zero_element)

    def __float__(self):
        embedded = self.embedded()
        #print "__float..."
        #try:
        result = (float(embedded.lower()) + float(embedded.upper())) / 2
        # except Exception:
        #     print "ERROR:", self
        # print "... end, result = %s" %result
        return result

    ## def __richcmp__(left, right, op):
    ##     print "rich"
    ##     return left._richcmp(right, op)

    ## def _richcmp(left, right, op):
    ##     return sign(left - right)

    def __repr__(self):
        embedded = self.embedded()
        symbolic = getattr(self, '_symbolic', None)
        if symbolic is None:
            return 'RNF%s' % embedded
        else:
            return '<RNF%s=%s>' % (symbolic, embedded)

    def _latex_(self):
        symbolic = getattr(self, '_symbolic', None)
        if symbolic is None:
            parent = self.parent()
            embedding = parent._exact_embedding
            symbolic = embedding(self)
        #print embedding(self), "%s" % (latex(embedding(self)))
        return "%s" % (latex(symbolic))

    def _maxima_(self, session=None):
        symbolic = getattr(self, '_symbolic', None)
        if symbolic is None:
            raise UnimplementedError, "Cannot make %s a Maxima number" % self
        else:
            return symbolic._maxima_()

    def _add_(self, other):
        result = NumberFieldElement_absolute._add_(self, other)
        # self_e = getattr(self, '_embedded', None)
        # if self_e is not None:
        #     other_e = getattr(other, '_embedded', None)
        #     if other_e is not None:
        #         result._embedded = self_e + other_e
        self_e = self.embedded()
        other_e = other.embedded()
        result._embedded = self_e + other_e
        return result

    def _sub_(self, other):
        result = NumberFieldElement_absolute._sub_(self, other)
        # self_e = getattr(self, '_embedded', None)
        # if self_e is not None:
        #     other_e = getattr(other, '_embedded', None)
        #     if other_e is not None:
        #         result._embedded = self_e - other_e
        self_e = self.embedded()
        other_e = other.embedded()
        result._embedded = self_e - other_e
        return result

    def __neg__(self):
        result = NumberFieldElement_absolute._neg_(self)
        # self_e = getattr(self, '_embedded', None)
        # if self_e is not None:
        #     result._embedded = -self_e
        self_e = self.embedded()
        result._embedded = -self_e
        return result

    def _mul_(self, other):
        result = NumberFieldElement_absolute._mul_(self, other)
        self_e = self.embedded()
        other_e = other.embedded()
        result._embedded = self_e * other_e
        return result
        
    def _div_(self, other):
        result = NumberFieldElement_absolute._div_(self, other)
        self_e = self.embedded()
        other_e = other.embedded()
        result._embedded = self_e / other_e
        return result

class RealNumberField_absolute(NumberField_absolute):
    """
    A RealNumberField knows its embedding into a RealIntervalField
    and use that for <, > comparisons.
    == comparison is exact, using the underlying numberfield.
    
    A RealIntervalField also knows an embedding into an exact field (SR)
    for the purpose of latexing.

    EXAMPLES::
    sage: field, field_values, morphism = number_field_elements_from_algebraics((sqrt(2), sqrt(3)))
    sage: emb_field = RealNumberField(field.polynomial(), 'a', embedding=morphism(field.gen(0)), exact_embedding=SR(morphism(field.gen(0))))
    sage: hom = field.hom([emb_field.gen(0)])
    sage: Integer(7)/5 < hom(field_values[0])
    True
    sage: hom(field_values[0]) < Integer(3)/2
    True
    sage: hom(field_values[0]) < hom(field_values[1])
    True
    sage: Integer(3)/2 < hom(field_values[1])
    True
    sage: hom(field_values[1]) < 2
    True
    """

    def __init__(self, polynomial, name=None, latex_name=None, check=True, embedding=None,
                 assume_disc_small=False, maximize_at_primes=None, exact_embedding=None):
        """
        Create a real number field.
        """
        NumberField_absolute.__init__(self, polynomial, name=name, check=check,
                                      embedding=embedding, latex_name=latex_name,
                                      assume_disc_small=assume_disc_small, maximize_at_primes=maximize_at_primes)
        self._standard_embedding = True
        self._element_class = RealNumberFieldElement
        self._zero_element = self(0)
        self._one_element =  self(1)
        self._exact_embedding = self.hom([exact_embedding])

    ## def specified_complex_embedding(self):
    ##     ## This is so that _symbolic_ works.
    ##     return self._exact_embedding
    ### FIXME: _symbolic_ leads to infinite recursion of LazyWrappers etc.


#### Repeat the whole exercise for quadratics... 

from sage.rings.number_field.number_field_element_quadratic import NumberFieldElement_quadratic

class RealNumberFieldElement_quadratic(NumberFieldElement_quadratic):
    
    def embedded(self):
        parent = self.parent()
        embedding = parent.coerce_embedding()
        e = embedding(self)
        return e
    
    def __repr__(self):
        embedded = self.embedded()
        return 'RNF%s' % embedded

class RealNumberField_quadratic(NumberField_quadratic):
    def __init__(self, polynomial, name=None, latex_name=None, check=True, embedding=None,
                 assume_disc_small=False, maximize_at_primes=None, exact_embedding=None):
        """
        Create a real number field.
        """
        #### This is copy/paste from number_field.py, Sage 5.11,
        #### modified to change the element type.
        NumberField_quadratic.__init__(self, polynomial, name=name, check=check,
                                      embedding=embedding, latex_name=latex_name,
                                      assume_disc_small=assume_disc_small, maximize_at_primes=maximize_at_primes)
        #### Note: Modify "NumberField_absolute.__init__()" to "NumberField_quadratic.__init__()" as the former causes TypeError in Sage 5.12. --Yuan
        self._standard_embedding = True
        self._element_class = RealNumberFieldElement_quadratic
        c, b, a = [Rational(t) for t in self.defining_polynomial().list()]
        # set the generator
        Dpoly = b*b - 4*a*c
        D = (Dpoly.numer() * Dpoly.denom()).squarefree_part(bound=10000)
        self._D = D
        parts = -b/(2*a), (Dpoly/D).sqrt()/(2*a)
        self._NumberField_generic__gen = self._element_class(self, parts)

        # we must set the flag _standard_embedding *before* any element creation
        # Note that in the following code, no element is built.
        emb = self.coerce_embedding()
        if emb is not None:
            rootD = RealNumberFieldElement_quadratic(self, (QQ(0),QQ(1)))
            if D > 0:
                from sage.rings.real_double import RDF
                self._standard_embedding = RDF.has_coerce_map_from(self) and RDF(rootD) > 0
            else:
                raise DataError, "RealNumberField_quadratic needs a positive discriminant"

        # we reset _NumberField_generic__gen has the flag standard_embedding
        # might be modified
        self._NumberField_generic__gen = self._element_class(self, parts)

        # NumberField_absolute.__init__(...) set _zero_element and
        # _one_element to NumberFieldElement_absolute values, which is
        # wrong (and dangerous; such elements can actually be used to
        # crash Sage: see #5316).  Overwrite them with correct values.
        #### Note: In my Sage 5.12, self(0) and self(1) return NumberFieldElement_quadratic. Should be RNFE_quadratic --Yuan
        #### self._zero_element = self(0) 
        #### self._one_element =  self(1)
        self._zero_element = RealNumberFieldElement_quadratic(self, QQ(0)) 
        self._one_element =  RealNumberFieldElement_quadratic(self, QQ(1))       
    
# The factory.

def RealNumberField(polynomial, name=None, latex_name=None, check=True, embedding=None,
                    assume_disc_small=False, maximize_at_primes=None, exact_embedding=None):
    if polynomial.degree() == 2:
        K = RealNumberField_quadratic(polynomial, name, latex_name, check, embedding,
                                      assume_disc_small=assume_disc_small, 
                                      maximize_at_primes=maximize_at_primes, 
                                      exact_embedding=exact_embedding)
    else:
        K = RealNumberField_absolute(polynomial, name, latex_name, check, embedding,
                                     assume_disc_small=assume_disc_small, 
                                     maximize_at_primes=maximize_at_primes, 
                                     exact_embedding=exact_embedding)
    return K

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

def is_real_number_field_element(x):
    try:
        x.embedded() # FIXME: this is a hack
        return True
    except AttributeError:
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

def is_all_the_same_real_number_field(values):
    real_number_field_seen = None
    for x in values:
        if is_real_number_field_element(x):
            if real_number_field_seen:
                if real_number_field_seen != x.parent():
                    return False, values
            else:
                real_number_field_seen = x.parent()
        elif not can_coerce_to_QQ(x):
            return False, values
    if real_number_field_seen:
        # Coerce any rationals to the number field.
        return True, [ real_number_field_seen(x) for x in values ]
    else:
        return False, values

def nice_field_values(symb_values, field=None):
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
def piecewise_function_from_breakpoints_slopes_and_values(bkpt, slopes, values, field=None):
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
                            fast_linear_function(slopes[i], intercepts[i])] for i in range(len(bkpt)-1) ])

def piecewise_function_from_breakpoints_and_values(bkpt, values, field=None):
    """
    bkpt and values are two parallel lists; assuming bpkt is sorted (increasing).
    Return a function.
    """
    if len(bkpt)!=len(values):
        raise ValueError, "Need to have the same number of breakpoints and values."
    slopes = [ (values[i+1]-values[i])/(bkpt[i+1]-bkpt[i]) for i in range(len(bkpt)-1) ]
    return piecewise_function_from_breakpoints_slopes_and_values(bkpt, slopes, values, field)

def piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None):
    """
    bkpt and slopes are two parallel lists (except that bkpt is one element longer).
    The function always has value 0 on bkpt[0].
    Return a function.
    """
    if len(bkpt)!=len(slopes)+1:
        raise ValueError, "Need to have one breakpoint more than slopes."
    values = [0]
    for i in range(1,len(bkpt)-1):
        values.append(values[i-1] + slopes[i - 1] * (bkpt[i] - bkpt[i-1]))
    return piecewise_function_from_breakpoints_slopes_and_values(bkpt, slopes, values, field)

def piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes, field=None):
    """
    interval_lengths and slopes are two parallel lists.
    The function always has value 0 on 0.
    Return a function.
    """
    if len(interval_lengths)!=len(slopes):
        raise ValueError, "Number of given interval_lengths and slopes needs to be equal."
    bkpt = []
    bkpt.append(0)
    for i in range(len(interval_lengths)):
        if interval_lengths[i] < 0:
            raise ValueError, "Interval lengths must be non-negative."
        bkpt.append(bkpt[i]+interval_lengths[i])
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field)

def discrete_function_from_points_and_values(points, values, field=None):
    if field is None:
        field = default_field
    # global symb_values
    symb_values = points + values
    field_values = nice_field_values(symb_values, field)
    points, values = field_values[0:len(points)], field_values[-len(values):]
    pieces = [ (singleton_interval(point), FastLinearFunction(0, value))
               for point, value in itertools.izip(points, values) ]
    return FastPiecewise(pieces, merge=False)

def limiting_slopes(fn):
    functions = fn.functions()
    return functions[0]._slope, functions[-1]._slope

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

def modified_delta_pi(fn, fn_values, pts, i, j):
    return fn_values[i] + fn_values[j] - fn(fractional(pts[i]+pts[j])) 

def modified_delta_pi2(fn, fn_values2, pts, i, j):
    return fn_values2[i] + fn(fractional(pts[j] - pts[i])) - fn_values2[j]  

def find_epsilon_interval(fn, perturb):
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

    def __init__(self, domain_intervals, directed_move):
        function = fast_linear_function(directed_move[0], directed_move[1])
        pieces = [ [interval_to_endpoints(interval), function] for interval in domain_intervals ]
        FastPiecewise.__init__(self, pieces)
        self.directed_move = directed_move       # needed?

    def __repr__(self):
        return "<FunctionalDirectedMove %s with domain %s>" % (self.directed_move, self.intervals())

    def sign(self):
        return self.directed_move[0]

    def is_functional(self):
        return True

    def __getitem__(self, item):
        return self.directed_move[item]

    def can_apply(self, x):
        try:
            self(x)
            return True
        except ValueError:
            return False

    def apply_ignoring_domain(self, x):
        move_sign = self.sign()
        if move_sign == 1:
            next_x = fractional(x + self.directed_move[1])
        elif move_sign == -1:
            next_x = fractional(self.directed_move[1]-x)
        return next_x

    def apply_to_interval(self, interval, inverse=False):
        # FIXME: Interval has to be a pair (a, b). This is only
        # compatible with intervals as used in FastPiecewise, but
        # nothing else in our code. 
        #
        # This does not do error checking.  Some code depends on this fact!
        # FIXME: This should be made clear in the name of this function.
        directed_move = self.directed_move
        move_sign = directed_move[0]
        if move_sign == 1:
            if inverse:
                result = (interval[0] - directed_move[1], interval[1] - directed_move[1])
            else:
                result = (interval[0] + directed_move[1], interval[1] + directed_move[1])
        elif move_sign == -1:
            result = (directed_move[1] - interval[1], directed_move[1] - interval[0])
        else:
            raise ValueError, "Move not valid: %s" % list(move)
        return result

    def apply_to_coho_interval(self, interval, inverse=False):
        # This does not do complete error checking.
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
        return [ self.apply_to_interval(interval) for interval in self.intervals() ] 

    def is_identity(self):
        return self.directed_move[0] == 1 and self.directed_move[1] == 0

    def minimal_triples(self): # unused
        """
        Does not output symmetric pairs!  Rather, maps positive translations to horizontal faces
        and negative translations to vertical faces.
        """
        if self.directed_move[0] == 1:                      # translation
            t = self.directed_move[1]
            if t >= 0:
                return [ (interval, [t], (interval[0] + t, interval[1] + t)) for interval in self.intervals() ]
            else:
                return [ ([-t], (interval[0] + t, interval[1] + t), interval) for interval in self.intervals() ]
        elif self.directed_move[0] == -1: 
            r = self.directed_move[1]
            return [ (interval, (r - interval[0], r - interval[1]), r) for interval in self.intervals() ]
        else:
            raise ValueError, "Move not valid: %s" % list(move)

@cached_function
def generate_functional_directed_moves(fn, intervals=None):
    """
    Compute the moves (translations and reflections) relevant for the given intervals
    (default: all uncovered intervals).
    """
    ### FIXME: Do we also have to take care of edges of some
    ### full-dimensional additive faces sometimes?
    if intervals==None:
        # Default is to generate moves for ALL uncovered intervals
        intervals = generate_uncovered_intervals(fn)
    moves = set()
    for face in generate_maximal_additive_faces(fn):
        if face.is_directed_move():
            fdm = face.functional_directed_move(intervals)
            if not fdm.is_identity() and find_interior_intersection(fdm.intervals(), intervals): #FIXME: why interior?
                moves.add(fdm)
    return list(moves)

def is_directed_move_possible(x, move):
    return move.can_apply(x)

def plot_moves(seed, moves, colors=None, ymin=0, ymax=1):
    if colors == None:
        colors = rainbow(len(moves))
    g = Graphics()
    g += line([(seed,ymin), (seed,ymax)], color="magenta")
    y = 0
    covered_interval = [0,1]
    for move, color in itertools.izip(moves, colors):
        next_x = move.apply_ignoring_domain(seed)
        arrow_interval = [min(seed, next_x), max(seed, next_x)]
        if (len(interval_intersection(covered_interval, arrow_interval)) == 2):
            y += 0.04
            covered_interval = arrow_interval
        else:
            y += 0.002
            covered_interval[0] = min(covered_interval[0], arrow_interval[0])
            covered_interval[1] = max(covered_interval[1], arrow_interval[1])
        midpoint_x = (seed + next_x) / 2
        if move[0] == -1:
            # Reflection
            bezier_y = y + min(0.03, 0.3 * float(abs(move[1]/2 - seed)))
            g += arrow(path = [[(seed, y), (midpoint_x, bezier_y), \
                                (next_x, y)]], \
                       head = 2, \
                       # legend_label = "xyzz"
                       color = color, \
                       zorder = 7
            )
            ## Plot the invariant point somehow?
            #g += point((move[1]/2, y + 0.03), color=color)
        elif move[0] == 1:
            # Translation
            g += arrow((seed, y), (next_x, y), color=color, zorder = 7)
        else:
            raise ValueError, "Bad move: %s" % list(move)
        g += text("%s" % list(move), (midpoint_x, y), \
                  vertical_alignment="bottom", \
                  horizontal_alignment="center", \
                  color=color, zorder = 7)
    return g

def plot_possible_and_impossible_directed_moves(seed, directed_moves, fn):
    colors = [ "blue" if is_directed_move_possible(seed, directed_move) \
               else "red" for directed_move in directed_moves ]
    return plot_moves(seed, directed_moves, colors)


def plot_walk(walk_dict, color="black", ymin=0, ymax=1, **options):
    #return point([ (x,0) for x in walk_dict.keys()])
    g = Graphics()
    for x in walk_dict.keys():
        g += line([(x,ymin), (x,ymax)], color=color, zorder = -4, **options)
    return g

def plot_intervals(intervals, ymin=0, ymax=1):
    g = Graphics()
    for interval in intervals:
        g += polygon([(interval[0], ymin), (interval[1], ymin), \
                      (interval[1], ymax), (interval[0], ymax)], 
                     color="yellow", zorder = -8)
    return g

import collections
_closed_or_open_or_halfopen_interval = collections.namedtuple('Interval', ['a', 'b', 'left_closed', 'right_closed'])

class closed_or_open_or_halfopen_interval (_closed_or_open_or_halfopen_interval):
    def __repr__(self):
        if self.a == self.b and self.left_closed and self.right_closed:
            r = "{" + repr(self.a) + "}"
        else:
            r = ("[" if self.left_closed else "(") \
                + repr(self.a) + ", " + repr(self.b) \
                + ("]" if self.right_closed else ")")
        return "<Int" + r + ">"

def closed_interval(a, b):
    return closed_or_open_or_halfopen_interval(a, b, True, True)

def open_interval(a, b):
    return closed_or_open_or_halfopen_interval(a, b, False, False)

def singleton_interval(a):
    return closed_or_open_or_halfopen_interval(a, a, True, True)

def left_open_interval(a, b):
    return closed_or_open_or_halfopen_interval(a, b, False, True)

def right_open_interval(a, b):
    return closed_or_open_or_halfopen_interval(a, b, True, False)

def interval_length(interval):
    if interval[1] >= interval[0]:
        return interval[1] - interval[0]
    return 0

def coho_interval_left_endpoint_with_epsilon(i):
    """Return (x, epsilon)
    where x is the left endpoint
    and epsilon is 0 if the interval is left closed and 1 otherwise.
    """
    if len(i) == 2:
        # old-fashioned closed interval
        return i[0], 0 # Scanning from the left, turn on at left endpoint.
    else:
        # coho interval
        return i.a, 0 if i.left_closed else 1

def coho_interval_right_endpoint_with_epsilon(i):
    """Return (x, epsilon)
    where x is the right endpoint
    and epsilon is 1 if the interval is right closed and 0 otherwise.
    """
    if len(i) == 2:
        # old-fashioned closed interval
        return i[1], 1 # Scanning from the left, turn off at right endpoint plus epsilon
    else:
        # coho interval
        return i.b, 1 if i.right_closed else 0

def scan_coho_interval_left_endpoints(interval_list, tag=None, delta=-1):
    """Generate events of the form `(x, epsilon), delta, tag.`

    This assumes that `interval_list` is sorted from left to right,
    and that the intervals are pairwise disjoint.
    """
    for i in interval_list:
        yield coho_interval_left_endpoint_with_epsilon(i), delta, tag

def scan_coho_interval_right_endpoints(interval_list, tag=None, delta=+1):
    """Generate events of the form `(x, epsilon), delta, tag.`

    This assumes that `interval_list` is sorted from left to right,
    and that the intervals are pairwise disjoint.
    """
    for i in interval_list:
        yield coho_interval_right_endpoint_with_epsilon(i), delta, tag

def scan_coho_interval_list(interval_list, tag=None, on_delta=-1, off_delta=+1):
    """Generate events of the form `(x, epsilon), delta, tag.`

    This assumes that `interval_list` is sorted, and 
    that the intervals are pairwise disjoint.

    delta is -1 for the beginning of an interval ('on').
    delta is +1 for the end of an interval ('off'). 

    This is so that the events sort lexicographically in a way that if
    we have intervals whose closures intersect in one point, such as
    [a, b) and [b, c], we see first the 'on' event and then the 'off'
    event.  In this way consumers of the scan can easily implement merging 
    of such intervals. 

    If merging is not desired, set on_delta=+1, off_delta=-1. 

    sage: list(scan_coho_interval_list([closed_or_open_or_halfopen_interval(1, 2, True, False), closed_or_open_or_halfopen_interval(2, 3, True, True)]))
    [((1, 0), -1, None), ((2, 0), -1, None), ((2, 0), 1, None), ((3, 1), 1, None)]
    """
    return merge(scan_coho_interval_left_endpoints(interval_list, tag, on_delta), 
                 scan_coho_interval_right_endpoints(interval_list, tag, off_delta))

## def scan_set_difference(a, b):
##     """`a` and `b` should be event generators."""

from heapq import *

def scan_union_of_coho_intervals_minus_union_of_coho_intervals(interval_lists, remove_lists):
    # Following uses the lexicographic comparison of the tuples.
    scan = merge(merge(*[scan_coho_interval_list(interval_list, True) for interval_list in interval_lists]),
                 merge(*[scan_coho_interval_list(remove_list, False) for remove_list in remove_lists]))
    interval_indicator = 0
    remove_indicator = 0
    on = False
    for ((x, epsilon), delta, tag) in scan:
        was_on = on
        if tag:                                       # interval event
            interval_indicator -= delta
            assert(interval_indicator) >= 0
        else:                                           # remove event
            remove_indicator -= delta
            assert(remove_indicator) >= 0
        now_on = interval_indicator > 0 and remove_indicator == 0
        if not was_on and now_on: # switched on
            yield (x, epsilon), -1, None
        elif was_on and not now_on: # switched off
            yield (x, epsilon), +1, None
        on = now_on
    # No unbounded intervals:
    assert interval_indicator == 0
    assert remove_indicator == 0

def intersection_of_coho_intervals(interval_lists):
    """Compute the intersection of the union of intervals. 
    
    Each interval_list must be sorted, but intervals may overlap.  In
    this case, the output is broken into non-overlapping intervals at
    the points where the overlap multiplicity changes.
    
    EXAMPLES:
    sage: list(intersection_of_coho_intervals([[[1,2]], [[2,3]]]))
    [<Int{2}>]
    sage: list(intersection_of_coho_intervals([[[1,2], [2,3]], [[0,4]]]))
    [<Int[1, 2)>, <Int{2}>, <Int(2, 3]>]
    sage: list(intersection_of_coho_intervals([[[1,3], [2,4]], [[0,5]]]))
    [<Int[1, 2)>, <Int[2, 3]>, <Int(3, 4]>]
    sage: list(intersection_of_coho_intervals([[[1,2], left_open_interval(2,3)], [[0,4]]]))
    [<Int[1, 2]>, <Int(2, 3]>]
    sage: list(intersection_of_coho_intervals([[[1,3]], [[2,4]]]))
    [<Int[2, 3]>]
    """
    scan = merge(*[scan_coho_interval_list(interval_list, tag=index) for index, interval_list in enumerate(interval_lists)])
    interval_indicators = [ 0 for interval_list in interval_lists ]
    (on_x, on_epsilon) = (None, None)
    for ((x, epsilon), delta, index) in scan:
        was_on = all(on > 0 for on in interval_indicators)
        interval_indicators[index] -= delta
        assert interval_indicators[index] >= 0
        now_on = all(on > 0 for on in interval_indicators)
        if was_on: 
            assert on_x is not None
            assert on_epsilon >= 0
            assert epsilon >= 0
            if (on_x, on_epsilon) < (x, epsilon):
                yield closed_or_open_or_halfopen_interval(on_x, x,
                                                          on_epsilon == 0, epsilon > 0)
        if now_on:
            (on_x, on_epsilon) = (x, epsilon)
        else:
            (on_x, on_epsilon) = (None, None)
    assert all(on == 0 for on in interval_indicators) # no unbounded intervals

def coho_interval_list_from_scan(scan):
    """Actually returns a generator."""
    indicator = 0
    (on_x, on_epsilon) = (None, None)
    for ((x, epsilon), delta, tag) in scan:
        was_on = indicator > 0
        indicator -= delta
        assert indicator >= 0
        now_on = indicator > 0
        if not was_on and now_on:                        # switched on
            (on_x, on_epsilon) = (x, epsilon)
        elif was_on and not now_on:                     # switched off
            assert on_x is not None
            assert on_epsilon >= 0
            assert epsilon >= 0
            if (on_x, on_epsilon) < (x, epsilon):
                yield closed_or_open_or_halfopen_interval(on_x, x,
                                                          on_epsilon == 0, epsilon > 0)
            (on_x, on_epsilon) = (None, None)
    assert indicator == 0

def union_of_coho_intervals_minus_union_of_coho_intervals(interval_lists, remove_lists):
    """Compute a list of closed/open/half-open intervals that represent
    the set difference of `interval` and the union of the intervals in
    `remove_list`.

    Assume each of the lists in `interval_lists' and `remove_lists` are sorted (and
    each pairwise disjoint).  Returns a sorted list.

    EXAMPLES::
    sage: union_of_coho_intervals_minus_union_of_coho_intervals([[[0,10]]], [[[2,2], [3,4]]])
    [<Int[0, 2)>, <Int(2, 3)>, <Int(4, 10]>]
    sage: union_of_coho_intervals_minus_union_of_coho_intervals([[[0, 10]]], [[[1, 7]], [[2, 5]]])
    [<Int[0, 1)>, <Int(7, 10]>]
    sage: union_of_coho_intervals_minus_union_of_coho_intervals([[[0,10], closed_or_open_or_halfopen_interval(10, 20, False, True)]], [])
    [<Int[0, 20]>]
    """
    gen = coho_interval_list_from_scan(scan_union_of_coho_intervals_minus_union_of_coho_intervals(interval_lists, remove_lists))
    return list(gen)

def proper_interval_list_from_scan(scan):
    """Return a generator of the proper intervals [a, b], a<b, in the `scan`.

    This ignores whether intervals are open/closed/half-open.
    """
    indicator = 0
    (on_x, on_epsilon) = (None, None)
    for ((x, epsilon), delta, tag) in scan:
        was_on = indicator > 0
        indicator -= delta
        assert indicator >= 0
        now_on = indicator > 0
        if not was_on and now_on:                        # switched on
            (on_x, on_epsilon) = (x, epsilon)
        elif was_on and not now_on:                     # switched off
            assert on_x is not None
            assert on_epsilon >= 0
            assert epsilon >= 0
            if on_x < x:
                yield [on_x, x]
            (on_x, on_epsilon) = (None, None)
    assert indicator == 0


# size has to be a positive integer
def lattice_plot(A, A0, t1, t2, size):
    size = size + 1
    x0 = A + (A0-A)/2
    p1 = points((x,y) for x in range(size) for y in range(size)) + points((-x,y) for x in range(size) for y in range(size))
    p2 = points((-x,-y) for x in range(size) for y in range(size)) + points((x,-y) for x in range(size) for y in range(size))
    p3 = plot((A-x0-x*t1)/t2, (x,-size + 1, size - 1), color = "red")
    p4 = plot((A0-x0-x*t1)/t2, (x,-size + 1,size - 1), color = "red")
    return p1+p2+p3+p4

# 

class UnimplementedError (Exception):
    pass

def generate_compatible_piecewise_function(components, component_slopes, field=None):
    intervals_and_slopes = []
    for component, slope in itertools.izip(components, component_slopes):
        intervals_and_slopes.extend([ (interval, slope) for interval in component ])
    intervals_and_slopes.sort()
    bkpt = [ int[0] for int, slope in intervals_and_slopes ] + [1]
    slopes = [ slope for int, slope in intervals_and_slopes ]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field)

def symbolic_piecewise(function):
    """
    Construct a vector-space-valued piecewise linear function
    compatible with the given `function`.  Each of the components of
    the function has a slope that is a basis vector of the vector
    space.  Returns three values: The constructed function, the list
    of components (including one component for every non-covered
    interval), and the field over which the original `function` is
    defined.
    """
    covered_intervals = generate_covered_intervals(function)
    uncovered_intervals = generate_uncovered_intervals(function)
    if uncovered_intervals:
        logging.warn(\
                     """There are non-covered intervals, so (1) the symbolic piecewise is
                     not suitable for proving extremality; and (2) in the current
                     implementation, there may be too many slope variables, since the 
                     relations between non-covered intervals are not taken into account.""")
        components = copy(covered_intervals)
        components.extend([int] for int in uncovered_intervals)
    else:
        components = covered_intervals
    # FIXME: fraction_field() required because parent could be Integer
    # Ring.  This happens, for example, for three_slope_limit().  
    # We really should have a function to retrieve the field of
    # a FastPiecewise.  But now even .base_ring() fails because
    # FastLinearFunction does not have a .base_ring() method.
    field = function(0).parent().fraction_field()
    vector_space = VectorSpace(field,len(components))
    slope_vars = vector_space.basis()
    symbolic = generate_compatible_piecewise_function(components, slope_vars, field)
    return symbolic, components, field

def generate_additivity_equations(function, symbolic, field, f=None):
    if f==None:
        f = find_f(function)
    equations = matrix(field, ([delta_pi(symbolic, x, y) 
                                for (x, y) in generate_additive_vertices(function) ]
                               + [symbolic(f)]
                               + [symbolic(1)]))
    return equations

def rescale_to_amplitude(perturb, amplitude):
    """For plotting purposes, rescale the function `perturb` so that its
    maximum absolute function value is `amplitude`.
    """
    current_amplitude = max([ abs(perturb(x)) for x in perturb.end_points() ])
    if current_amplitude != 0:
        return perturb * (amplitude/current_amplitude)
    else:
        return perturb

# Global figsize for all plots made by show_plots.
show_plots_figsize = 10

def show_plot(graphics, show_plots, tag, **show_kwds):
    """
    Display or save `graphics`.

    `show_plots` should be `False` (do nothing), 
    `True` (use `show` to display on screen),
    a string (file name format such as "FILENAME-%s.pdf", 
    where %s is replaced by `tag`.
    """
    if isinstance(show_plots, str):
        graphics.save(show_plots % tag, figsize=show_plots_figsize, **show_kwds)
    elif show_plots:
        graphics.show(figsize=show_plots_figsize, **show_kwds)

def check_perturbation(fn, perturb, show_plots=False, show_plot_tag='perturbation', xmin=0, xmax=1, **show_kwds):
    epsilon_interval = fn._epsilon_interval = find_epsilon_interval(fn, perturb)
    epsilon = min(abs(epsilon_interval[0]), epsilon_interval[1])
    logging.info("Epsilon for constructed perturbation: %s" % epsilon)
    if show_plots:
        logging.info("Plotting perturbation...")
        p = plot(fn, xmin=xmin, xmax=xmax, color='black', thickness=2, legend_label="original function")
        p += plot(fn + epsilon_interval[0] * perturb, xmin=xmin, xmax=xmax, color='red', legend_label="-perturbed (min)")
        p += plot(fn + epsilon_interval[1] * perturb, xmin=xmin, xmax=xmax, color='blue', legend_label="+perturbed (max)")
        if -epsilon != epsilon_interval[0]:
            p += plot(fn + (-epsilon) * perturb, xmin=xmin, xmax=xmax, color='orange', legend_label="-perturbed (matches max)")
        elif epsilon != epsilon_interval[1]:
            p += plot(fn + epsilon * perturb, xmin=xmin, xmax=xmax, color='cyan', legend_label="+perturbed (matches min)")
        p += plot(rescale_to_amplitude(perturb, 1/10), xmin=xmin, xmax=xmax, color='magenta', legend_label="perturbation (rescaled)")
        show_plot(p, show_plots, tag=show_plot_tag, **show_kwds)
        logging.info("Plotting perturbation... done")
    assert epsilon > 0, "Epsilon should be positive, something is wrong"
    logging.info("Thus the function is not extreme.")

def finite_dimensional_extremality_test(function, show_plots=False, f=None):
    """
    Solve a homogeneous linear system of additivity equations with one
    slope variable for every component (including every non-covered
    interval).  Return a boolean that indicates whether the system has
    a nontrivial solution.
    """
    if not function.is_continuous_defined():
        logging.warn("This is a discontinuous function; finite dimensional extremality test does not handle it yet.")
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

def generate_type_1_vertices(fn, comparison):
    """A generator...
    "...'general' refers to the fact that it outputs 6-tuples (x,xeps,y,yeps,z,zeps).
    FIXME: it currently does not take care of any discontinuities at all.
    """
    bkpt = fn.end_points()
    return ( (x,0,y,0,x+y,0) for x in bkpt for y in bkpt if x <= y and comparison(delta_pi(fn,x,y), 0) ) # generator comprehension
    ## Equivalent to previous line:
    # for x in bkpt:
    #     for y in bkpt:
    #         if x <= y and comparison(delta_pi(fn,x,y), 0):
    #             yield (x,0,y,0,x+y,0)

def generate_type_2_vertices(fn, comparison):
    bkpt = fn.end_points()
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]
    return ( (x,0,z-x,0,z,0) for x in bkpt for z in bkpt2 if x < z < 1+x and comparison(delta_pi(fn, x, z-x), 0) ) # generator comprehension

@cached_function
def generate_additive_vertices(fn):
    """
    This function does not output limit vertices for discontinuous function in any way.
    We are returning a set, so that duplicates are removed, and so the result can be cached for later use.
    """
    return { (x, y) 
             for (x, xeps, y, yeps, z, zeps) in itertools.chain(generate_type_1_vertices(fn, operator.eq),
                                                                generate_type_2_vertices(fn, operator.eq))
             if xeps==yeps==zeps==0 }

@cached_function
def generate_nonsubadditive_vertices(fn):
    """
    This function does not output limit vertices for discontinuous function in any way.
    We are returning a set, so that duplicates are removed, and so the result can be cached for later use.
    """
    return { (x, y) 
             for (x, xeps, y, yeps, z, zeps) in itertools.chain(generate_type_1_vertices(fn, operator.lt),
                                                                generate_type_2_vertices(fn, operator.lt))
             if xeps==yeps==zeps==0 }

class MaximumNumberOfIterationsReached(Exception):
    pass

def extremality_test(fn, show_plots = False, f=None, max_num_it = 1000, perturbation_style=default_perturbation_style, phase_1 = False, finite_dimensional_test_first = False, use_new_code=True):
    do_phase_1_lifting = False
    if f == None:
        f = find_f(fn, no_error_if_not_minimal_anyway=True)
    if f == None or not minimality_test(fn, show_plots=show_plots, f=f):
        logging.info("Not minimal, thus not extreme.")
        if not phase_1:
            return False
        else:
            do_phase_1_lifting = True
    covered_intervals = generate_covered_intervals(fn)
    uncovered_intervals = generate_uncovered_intervals(fn)
    if show_plots:
        logging.info("Plotting covered intervals...")
        show_plot(plot_covered_intervals(fn), show_plots, tag='covered_intervals')
        logging.info("Plotting covered intervals... done")
    if not uncovered_intervals:
        logging.info("All intervals are covered (or connected-to-covered). %s components." % len(covered_intervals))
        return finite_dimensional_extremality_test(fn, show_plots, f=f)
    else:
        logging.info("Uncovered intervals: %s", (uncovered_intervals,))
        if do_phase_1_lifting or finite_dimensional_test_first:
            # First try the finite dimensional one.
            if not finite_dimensional_extremality_test(fn, show_plots):
                return False
        # Now do the magic.
        moves = generate_functional_directed_moves(fn)
        logging.debug("Moves relevant for these intervals: %s" % (moves,))
        if use_new_code:
            seed, stab_int, walk_list = find_generic_seed_with_completion(fn, show_plots=show_plots, max_num_it=max_num_it) # may raise MaximumNumberOfIterationsReached
            if not seed:
                logging.info("Dense orbits in all non-covered intervals.  According to conjectures, this means that the function is extreme.")
                return True
        else:
            seed, stab_int, walk_list = find_generic_seed(fn, max_num_it=max_num_it) # may raise MaximumNumberOfIterationsReached
        fn._seed = seed
        fn._stab_int = stab_int
        fn._walk_list = walk_list
        if show_plots:
            logging.info("Plotting moves and reachable orbit...")
            # FIXME: Visualize stability intervals?
            g = (plot_walk(walk_list,thickness=0.7) + 
                 plot_possible_and_impossible_directed_moves(seed, moves, fn) + 
                 plot_intervals(uncovered_intervals) + plot_covered_intervals(fn))
            show_plot(g, show_plots, tag='moves')
            logging.info("Plotting moves and reachable orbit... done")
        perturb = fn._perturbation = approx_discts_function(walk_list, stab_int, perturbation_style=perturbation_style, function=fn)
        check_perturbation(fn, perturb, show_plots=show_plots)
        return False

def lift(fn, show_plots = False, **kwds):
    if extremality_test(fn, show_plots, **kwds):
        return fn
    else:
        perturbed = fn._lifted = fn + fn._epsilon_interval[1] * fn._perturbation
        ## Following is strictly experimental: It may change what "f" is.
        if 'phase_1' in kwds and kwds['phase_1']:
            perturbed = rescale_to_amplitude(perturbed, 1)
        return perturbed

def lift_until_extreme(fn, show_plots = False, pause = False, **kwds):
    next, fn = fn, None
    while next != fn:
        fn = next
        next = lift(fn, show_plots, **kwds)
        if pause and next != fn:
            raw_input("Press enter to continue")
    return next

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

def random_piecewise_function(xgrid, ygrid, continuity=True):
    xvalues = [0] + [ x/xgrid for x in range(1, xgrid) ] + [1]
    f = randint(1, xgrid - 1)
    yvalues = [0] + [ randint(0, ygrid) / ygrid for i in range(1, f) ] + [1] + [ randint(0, ygrid) / ygrid for i in range(f+1, xgrid) ]+ [0]
    if continuity:        return piecewise_function_from_breakpoints_and_values(xvalues, yvalues)
    else:
        piece1 = [ [singleton_interval(xvalues[i]), FastLinearFunction(0, yvalues[i])] for i in range(xgrid+1) ]
        leftlimits = [0] + [ randint(0, ygrid) / ygrid for i in range(1, xgrid+1) ]
        rightlimits = [ randint(0, ygrid) / ygrid for i in range(1, xgrid+1) ] + [0]
        slopes = [ (leftlimits[i+1] - rightlimits[i]) * xgrid for i in range(0, xgrid) ]
        intercepts = [ rightlimits[i] - xvalues[i] * slopes[i] for i in range(0, xgrid) ]
        piece2 = [ [open_interval(xvalues[i], xvalues[i+1]), FastLinearFunction(slopes[i], intercepts[i])] for i in range(xgrid) ]
        pieces = [piece1[0]]
        for i in range(xgrid):
            pieces += [piece2[i], piece1[i+1]]
        return FastPiecewise(pieces, merge=True)

def is_QQ_linearly_independent(*numbers):
    """
    Test if `numbers` are linearly independent over `QQ`.

    EXAMPLES::
        sage: logging.disable(logging.INFO); is_QQ_linearly_independent()
        True
        sage: logging.disable(logging.INFO); is_QQ_linearly_independent(1)
        True
        sage: logging.disable(logging.INFO); is_QQ_linearly_independent(0)
        False
        sage: logging.disable(logging.INFO); is_QQ_linearly_independent(1,2)
        False
        sage: logging.disable(logging.INFO); is_QQ_linearly_independent(1,sqrt(2))
        True
        sage: logging.disable(logging.INFO); is_QQ_linearly_independent(1+sqrt(2),sqrt(2),1)
        False
    """
    # trivial cases
    if len(numbers) == 0:
        return True
    elif len(numbers) == 1:
        return numbers[0] != 0
    # fast path for rationals
    all_QQ, numbers = is_all_QQ(numbers)
    if all_QQ:
        return False
    # otherwise try to coerce to common number field
    numbers = nice_field_values(numbers, RealNumberField)
    if not is_real_number_field_element(numbers[0]):
        raise ValueError, "Q-linear independence test only implemented for algebraic numbers"
    coordinate_matrix = matrix(QQ, [x.parent().0.coordinates_in_terms_of_powers()(x) for x in numbers])
    return rank(coordinate_matrix) == len(numbers)

def interval_list_intersection(interval_list_1, interval_list_2, interiors=False):
    """
    Return a list of the intersections of the intervals
    in `interval_list_1` and the intervals in `interval_list_2`.

    Assumes the input lists are sorted.
    The output is sorted.

    EXAMPLES::
    sage: interval_list_intersection([[2,3]], [[1,2], [3,5]])
    [[2], [3]]
    sage: interval_list_intersection([[2,6]], [[1,3], [5,7], [7,9]])
    [[2, 3], [5, 6]]
    """
    overlap = []
    # FIXME: Should be able to speed up by merging algorithm.
    for int1 in interval_list_1:
        for int2 in interval_list_2:
            overlapped_int = interval_intersection(int1,int2)
            if len(overlapped_int) >= 2 or (not interiors and len(overlapped_int) >= 1):
                overlap.append(overlapped_int)
    return sorted(overlap)

def compose_directed_moves(A, B, interiors=False, show_plots=False):
    """
    Compute the directed move that corresponds to the directed move `A` after `B`.
    
    EXAMPLES::
        sage: compose_directed_moves(FunctionalDirectedMove([(5/10,7/10)],(1, 2/10)),FunctionalDirectedMove([(2/10,4/10)],(1,2/10)))
        <FunctionalDirectedMove (1, 2/5) with domain [(3/10, 2/5)]>
    """
    #print result_domain_intervals
    if A.is_functional() and B.is_functional():
        A_domain_preimages = [ B.apply_to_interval(A_domain_interval, inverse=True) \
                               for A_domain_interval in A.intervals() ]
        result_domain_intervals = interval_list_intersection(A_domain_preimages, B.intervals(), interiors=interiors)
        if len(result_domain_intervals) > 0:
            result = FunctionalDirectedMove([ interval_to_endpoints(I) for I in result_domain_intervals ], (A[0] * B[0], A[0] * B[1] + A[1]))
        else:
            result = None
    elif not A.is_functional() and B.is_functional():
        A_domain_preimages = [ B.apply_to_interval(A_domain_interval, inverse=True) \
                               for A_domain_interval in A.intervals() ]
        # FIXME: This is a version of interval_list_intersection.  Should be able to speed up by merging algorithm. 
        interval_pairs = []
        for A_domain_preimage, A_range in itertools.izip(A_domain_preimages, A.range_intervals()):
            for B_domain in B.intervals():
                overlapped_int = interval_intersection(A_domain_preimage, B_domain)
                if len(overlapped_int) >= 2 or (not interiors and len(overlapped_int) >= 1):
                    interval_pairs.append((interval_to_endpoints(overlapped_int), A_range))
        if interval_pairs:
            result = DenseDirectedMove(interval_pairs)
        else:
            result = None
    elif A.is_functional() and not B.is_functional():
        interval_pairs = []
        for A_domain_interval in A.intervals():
            for B_domain_interval, B_range_interval in B.interval_pairs():
                overlapped_int = interval_intersection(A_domain_interval, B_range_interval)
                if len(overlapped_int) >= 2 or (not interiors and len(overlapped_int) >= 1):
                    interval_pairs.append((B_domain_interval, A.apply_to_interval(interval_to_endpoints(overlapped_int))))
        if interval_pairs:
            result = DenseDirectedMove(interval_pairs)
        else:
            result = None
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
        <FunctionalDirectedMove (1, 0) with domain [(3/10, 7/20), (9/20, 1/2)]>
    """
    if A.directed_move != B.directed_move:
        raise ValueError, "Cannot merge, moves have different operations"
    #merge_two_comp(A.intervals(), B.intervals(), one_point_overlap_suffices=True), 
    C = FunctionalDirectedMove(\
                               A.intervals() + B.intervals(),  # constructor takes care of merging
                               A.directed_move)
    if show_plots:
        p = plot(C, color="cyan", legend_label="C = A merge B", thickness=10)
        p += plot(A, color="green", legend_label="A = %s" % A )
        p += plot(B, color="blue", legend_label="B = %s" % B)
        show_plot(p, show_plots, tag='merge_functional_directed_moves')
    return C

def plot_directed_moves(dmoves):
    return sum(plot(dm) for dm in dmoves)

def reduce_with_dense_moves(functional_directed_move, dense_moves, show_plots=False):
    """
    EXAMPLES::
        sage: reduce_with_dense_moves(FunctionalDirectedMove([[3/10,7/10]],(1, 1/10)), [DenseDirectedMove([[[2/10,6/10],[2/10,6/10]]])])
        <FunctionalDirectedMove (1, 1/10) with domain [(1/2, 7/10)]>
        sage: reduce_with_dense_moves(FunctionalDirectedMove([[1/10,7/10]],(1, 1/10)), [DenseDirectedMove([[[7/20,5/10],[3/10,5/10]]]), DenseDirectedMove([[[6/20,6/10],[4/10,6/10]]])])
        <FunctionalDirectedMove (1, 1/10) with domain [(1/10, 3/10), (1/2, 7/10)]>
    """
    remove_lists = []
    for domain, codomain in itertools.chain(*[ dense_move.interval_pairs() for dense_move in dense_moves ]):
        remove_list = []
        int = interval_intersection(functional_directed_move.apply_to_interval(codomain, inverse=True) , domain)
        if len(int) == 2:
            remove_list.append(closed_or_open_or_halfopen_interval(int[0], int[1], True, True))
        remove_lists.append(remove_list)  # Each remove_list is sorted because each interval is a subinterval of the domain interval.
    #print remove_lists
    scan_difference = scan_union_of_coho_intervals_minus_union_of_coho_intervals([functional_directed_move.intervals()], remove_lists)
    proper_difference = list(proper_interval_list_from_scan(scan_difference))
    if proper_difference:
        result = FunctionalDirectedMove(proper_difference, functional_directed_move.directed_move)
    else:
        result = None
    if show_plots:
        p = plot(functional_directed_move, color="yellow", thickness=8)
        p += plot_directed_moves(dense_moves)
        if result:
            p += plot(result)
        show_plot(p, show_plots, tag='reduce_with_dense_moves')
    return result

class DirectedMoveCompositionCompletion:

    def __init__(self, directed_moves, show_plots=False):
        self.show_plots = show_plots
        self.move_dict = dict()
        self.dense_moves = set()
        self.any_change = False
        for move in directed_moves:
            self.add_move(move)

    def reduce_move_dict_with_dense_moves(self, dense_moves):
        new_move_dict = dict()
        for key, move in self.move_dict.items():
            new_move = reduce_with_dense_moves(move, dense_moves)
            if new_move:
                new_move_dict[key] = new_move
        self.move_dict = new_move_dict

    def upgrade_or_reduce_dense_interval_pair(self, a_domain, a_codomain):
        another_pass = True
        while another_pass:
            another_pass = False
            for b in self.dense_moves:
                for (b_domain, b_codomain) in b.interval_pairs():
                    if (b_domain[0] <= a_domain[0] and a_domain[1] <= b_domain[1]
                        and b_codomain[0] <= a_codomain[0] and a_codomain[1] <= b_codomain[1]):
                        # is dominated by existing rectangle, exit.
                        return None, None
                    elif (a_domain[0] <= b_domain[0] and b_domain[1] <= a_domain[1]
                        and a_codomain[0] <= b_codomain[0] and b_codomain[1] <= a_codomain[1]):
                        # dominates existing rectangle, do nothing (we take care of that later).
                        pass
                    elif (a_domain[0] == b_domain[0] and a_domain[1] == b_domain[1]
                          and len(interval_intersection(a_codomain, b_codomain)) >= 1):
                          # simple vertical merge
                        a_codomain = ((min(a_codomain[0], b_codomain[0]), max(a_codomain[1], b_codomain[1])))
                        another_pass = True
                    elif (a_codomain[0] == b_codomain[0] and a_codomain[1] == b_codomain[1]
                          and len(interval_intersection(a_domain, b_domain)) >= 1):
                          # simple horizontal merge
                        a_domain = ((min(a_domain[0], b_domain[0]), max(a_domain[1], b_domain[1])))
                        another_pass = True
                    elif len(interval_intersection(a_domain, b_domain)) == 2 \
                       and len(interval_intersection(a_codomain, b_codomain)) == 2:
                        # full-dimensional intersection, extend to big rectangle.
                        logging.info("Applying rectangle lemma")
                        a_domain = ((min(a_domain[0], b_domain[0]), max(a_domain[1], b_domain[1])))
                        a_codomain = ((min(a_codomain[0], b_codomain[0]), max(a_codomain[1], b_codomain[1])))
                        another_pass = True
        return a_domain, a_codomain

    def add_move(self, c):
        if c.is_functional():
            reduced = reduce_with_dense_moves(c, self.dense_moves)
            if reduced is None:
                return
            cdm = c.directed_move
            if cdm in self.move_dict:
                merged = merge_functional_directed_moves(self.move_dict[cdm], reduced, show_plots=False)

                if merged.intervals() != self.move_dict[cdm].intervals():
                    # Cannot compare the functions themselves because of the "hash" magic of FastPiecewise.
                    #print "merge: changed from %s to %s" % (self.move_dict[cdm], merged)
                    self.move_dict[cdm] = merged
                    self.any_change = True
                else:
                    #print "merge: same"
                    pass
            # elif is_move_dominated_by_dense_moves(c, self.dense_moves):
            #     pass
            else:
                self.move_dict[cdm] = reduced
                self.any_change = True
        else:
            # dense move.
            new_dense_moves = []
            for (c_domain, c_codomain) in c.interval_pairs():
                c_domain, c_codomain = self.upgrade_or_reduce_dense_interval_pair(c_domain, c_codomain)
                if c_domain:
                    new_dense_moves.append(DenseDirectedMove([(c_domain, c_codomain)]))
            if not new_dense_moves:
                return 
            dominated_dense_list = [ move for move in self.dense_moves if is_move_dominated_by_dense_moves(move, new_dense_moves) ]
            for move in dominated_dense_list:
                self.dense_moves.remove(move)
            for m in new_dense_moves:
                self.dense_moves.add(m)
            # dominated_functional_key_list = [ key for key, move in self.move_dict.items() if is_move_dominated_by_dense_moves(move, self.dense_moves) ]
            # for key in dominated_functional_key_list:
            #     self.move_dict.pop(key)
            self.reduce_move_dict_with_dense_moves(new_dense_moves)
            self.any_change = True

    def plot(self, *args, **kwargs):
        return plot_directed_moves(list(self.dense_moves) + list(self.move_dict.values()))

    def maybe_show_plot(self):
        if self.show_plots:
            logging.info("Plotting...")
            show_plot(self.plot(), self.show_plots, 'completion')
            logging.info("Plotting... done")

    def complete_one_round(self):
        self.maybe_show_plot()
        logging.info("Completing %d functional directed moves and %d dense directed moves..." % (len(self.move_dict), len(self.dense_moves)))
        self.any_change = False
        critical_pairs = ([ (a, b) for a in list(self.dense_moves) for b in list(self.dense_moves) ] 
                          + [ (a, b) for a in self.move_dict.keys() for b in list(self.dense_moves) + self.move_dict.keys() ] 
                          + [ (a, b) for a in list(self.dense_moves) for b in  self.move_dict.keys() ])
        for (a, b) in critical_pairs:
            # Get the most current versions of the directed moves.
            # FIXME: We should rather implement a better completion
            # algorithm.
            if type(a) == tuple or type(a) == list:
                a = self.move_dict.get(a, None)
            if type(b) == tuple or type(b) == list:
                b = self.move_dict.get(b, None)

            if not a or not b: 
                continue                                    # critical pair has been killed

            if not a.is_functional() and not b.is_functional():
                new_pairs = []
                for (a_domain, a_codomain) in a.interval_pairs():
                    for (b_domain, b_codomain) in b.interval_pairs():
                        if len(interval_intersection(a_domain, b_domain)) == 2 \
                           and len(interval_intersection(a_codomain, b_codomain)) == 2:
                            # full-dimensional intersection, extend to big rectangle;
                            # but this is taken care of in add_move.
                            pass
                        elif len(interval_intersection(a_codomain, b_domain)) == 2:
                            # composition of dense moves
                            new_pairs.append((a_domain, b_codomain))
                if new_pairs:
                    d = DenseDirectedMove(new_pairs)
                    if not is_move_dominated_by_dense_moves(d, self.dense_moves):
                        logging.info("New dense move from dense-dense composition: %s" % d)
                        self.add_move(d)
                        # self.maybe_show_plot()
            else:
                if a.is_functional() and b.is_functional():
                    d = check_for_dense_move(a, b)
                    if d and not is_move_dominated_by_dense_moves(d, self.dense_moves):
                        logging.info("New dense move from strip lemma: %s" % d)
                        self.add_move(d)
                        # self.maybe_show_plot()
                c = compose_directed_moves(a, b, interiors=True)   ## experimental - only interiors
                if c:
                    self.add_move(c)

    def complete(self, max_num_rounds=8, error_if_max_num_rounds_exceeded=True):
        num_rounds = 0
        while self.any_change and (not max_num_rounds or num_rounds < max_num_rounds):
            self.complete_one_round()
            num_rounds += 1
        if max_num_rounds and num_rounds == max_num_rounds:
            if error_if_max_num_rounds_exceeded:
                raise MaximumNumberOfIterationsReached, "Reached %d rounds of the completion procedure, found %d directed moves and %d dense directed moves, stopping." % (num_rounds, len(self.move_dict), len(self.dense_moves))
            else:
                logging.info("Reached %d rounds of the completion procedure, found %d directed moves and %d dense directed moves, stopping." % (num_rounds, len(self.move_dict), len(self.dense_moves)))
        else:
            logging.info("Completion finished.  Found %d directed moves and %d dense directed moves." 
                         % (len(self.move_dict), len(self.dense_moves)))


    def results(self):
        ## FIXME: Should return the dense moves somehow as well.
        # if self.dense_moves:
        #     raise UnimplementedError, "Dense moves found, handling them in the following code is not implemented yet."
        return list(self.move_dict.values())


def directed_move_composition_completion(directed_moves, show_plots=False, max_num_rounds=8, error_if_max_num_rounds_exceeded=True):
    completion = DirectedMoveCompositionCompletion(directed_moves, show_plots=show_plots)
    completion.complete(max_num_rounds=max_num_rounds, error_if_max_num_rounds_exceeded=error_if_max_num_rounds_exceeded)
    return completion.results()

@cached_function
def generate_directed_move_composition_completion(fn, show_plots=False, max_num_rounds=8, error_if_max_num_rounds_exceeded=True):
    completion = getattr(fn, "_completion", None)
    if not completion:
        functional_directed_moves = generate_functional_directed_moves(fn)
        completion = fn._completion = DirectedMoveCompositionCompletion(functional_directed_moves, show_plots=show_plots)
    completion.complete(max_num_rounds=max_num_rounds, error_if_max_num_rounds_exceeded=error_if_max_num_rounds_exceeded)
    return completion.results()

def apply_functional_directed_moves(functional_directed_moves, seed):
    """
    Return a dictionary whose keys are the reachable orbit of `seed`.

    If `functional_directed_moves` is complete under compositions,
    then this computes the reachable orbit of `seed`, just like
    `deterministic_walk` would.
    """
    #orbit = set()
    orbit = dict()
    for fdm in functional_directed_moves:
        try:
            #orbit.add(fdm(seed))
            element = fdm(seed)
            if element in orbit:
                orbit[element].append(fdm)
            else:
                orbit[element] = [fdm]
        except ValueError:
            pass
    #return sorted(orbit)
    return orbit

def scan_sign_contradiction_point(fdm):
    if fdm.sign() == -1:
        invariant_point = fdm[1] / 2
        if fdm.can_apply(invariant_point):
            assert fdm(invariant_point) == invariant_point
        yield ((invariant_point, 0), 0, fdm)
        yield ((invariant_point, 1), 0, fdm)

def scan_sign_contradiction_points(functional_directed_moves):
     scans = [ scan_sign_contradiction_point(fdm) for fdm in functional_directed_moves ]
     return merge(*scans)

def scan_domains_of_moves(functional_directed_moves):
     scans = [ scan_coho_interval_list(fdm.intervals(), fdm) for fdm in functional_directed_moves ]
     return merge(*scans)

def find_decomposition_into_intervals_with_same_moves(functional_directed_moves, separate_by_sign_contradiction=True):
    scan = scan_domains_of_moves(functional_directed_moves)
    if separate_by_sign_contradiction:
        scan = merge(scan, \
                     scan_sign_contradiction_points(functional_directed_moves))
    moves = set()
    (on_x, on_epsilon) = (None, None)
    for ((x, epsilon), delta, move) in scan:
        if on_x and (on_x, on_epsilon) < (x, epsilon):
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
        elif delta == 0:                       # an invariant point
            pass
        else:
            raise ValueError, "Bad scan item"

def find_decomposition_into_stability_intervals_with_completion(fn, show_plots=False, max_num_it=None):
    fn._stability_orbits = []
    completion = generate_directed_move_composition_completion(fn, show_plots=show_plots)

    decomposition = find_decomposition_into_intervals_with_same_moves(completion)
    done_intervals = set()

    for (interval, moves) in decomposition:
        if interval not in done_intervals:
            #print interval
            orbit = set()
            walk_dict = dict()
            seed = (interval.a + interval.b) / 2
            sign_contradiction = False
            for move in moves:
                moved_interval = move.apply_to_coho_interval(interval)
                #print "Applying move %s to %s gives %s." % (move, interval, moved_interval)
                moved_seed = move(seed)
                walk_sign = move.sign()
                done_intervals.add(moved_interval)
                orbit.add(moved_interval)
                if moved_seed in walk_dict and walk_dict[moved_seed][0] != walk_sign:
                    sign_contradiction = True
                walk_dict[moved_seed] = [walk_sign, None, None] 
            if sign_contradiction:
                for y in walk_dict.values():
                    y[0] = 0
            stability_orbit = (list(orbit), walk_dict, None)
            fn._stability_orbits.append(stability_orbit)
    logging.info("Total: %s stability orbits, lengths: %s" % (len(fn._stability_orbits), \
                    [ ("%s+" if to_do else "%s") % len(shifted_stability_intervals) \
                      for (shifted_stability_intervals, walk_dict, to_do) in fn._stability_orbits ]))

def find_generic_seed_with_completion(fn, show_plots=False, max_num_it=None):
    # Ugly compatibility interface.
    find_decomposition_into_stability_intervals_with_completion(fn, show_plots=show_plots)
    for (orbit, walk_dict, _) in fn._stability_orbits:
        int = orbit[0]
        if interval_length(int) > 0:
            seed = (int.a + int.b) / 2
            stab_int = closed_or_open_or_halfopen_interval(int.a - seed, int.b - seed, int.left_closed, int.right_closed)
            return (seed, stab_int, walk_dict)
    return None, None, None

class DenseDirectedMove ():

    def __init__(self, interval_pairs):
        self._interval_pairs = interval_pairs

    def __repr__(self):
        return "<DenseDirectedMove %s>" % self._interval_pairs

    def is_functional(self):
        return False

    def plot(self, *args, **kwds):
        return sum([polygon(((domain[0], codomain[0]), (domain[1], codomain[0]), (domain[1], codomain[1]), (domain[0], codomain[1])), rgbcolor=kwds.get("rgbcolor", "cyan"), alpha=0.5) + polygon(((domain[0], codomain[0]), (domain[1], codomain[0]), (domain[1], codomain[1]), (domain[0], codomain[1])), color="red", fill=False) for (domain, codomain) in self._interval_pairs])

    def intervals(self):
        return [ domain_interval for (domain_interval, range_interval) in self._interval_pairs ]

    def range_intervals(self):
        return [ range_interval for (domain_interval, range_interval) in self._interval_pairs ]

    def interval_pairs(self):
        return self._interval_pairs
    
def check_for_dense_move(m1, m2):
    # the strip lemma.
    if m1.sign() == 1 and m2.sign() == 1: 
        t1 = m1[1]
        t2 = m2[1]
        if is_QQ_linearly_independent(t1, t2) and t1 >= 0 and t2 >= 0:
            dense_intervals = []
            for (l1, u1) in m1.intervals():
                for (l2, u2) in m2.intervals():
                    L = max(l1, l2)
                    U = min(u1 + t1, u2 + t2)
                    if t1 + t2 <= U - L:
                        dense_intervals.append((L, U))
            if dense_intervals:
                return DenseDirectedMove([(I, I) for I in dense_intervals])
    return None

def is_interval_pair_dominated_by_dense_move(domain_interval, range_interval, dense_move):
    for (dense_domain_interval, dense_range_interval) in dense_move.interval_pairs():
        if interval_contained_in_interval(domain_interval, dense_domain_interval) \
           and interval_contained_in_interval(range_interval, dense_range_interval):
            return True
    return False

def is_interval_pair_dominated_by_dense_moves(domain_interval, range_interval, dense_moves):
    for dense_move in dense_moves:
        if is_interval_pair_dominated_by_dense_move(domain_interval, range_interval, dense_move):
            return True
    return False

def is_move_dominated_by_dense_moves(move, dense_moves):
    for (domain_interval, range_interval) in itertools.izip(move.intervals(), move.range_intervals()):
        if is_interval_pair_dominated_by_dense_moves(domain_interval, range_interval, dense_moves):
            pass
        else:
            return False
    #print "Dominated: %s" % move
    return True

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
    completion = directed_move_composition_completion(dmoves, max_num_rounds=None)
    plot_directed_moves(completion).show(figsize=40)
    
