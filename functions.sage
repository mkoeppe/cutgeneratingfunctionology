import logging

logging.basicConfig(format='%(levelname)s: %(asctime)s %(message)s', level=logging.INFO)

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
    ## We now use lambda functions instead of Sage symbolics for plotting, 
    ## as those give strange errors when combined with our RealNumberFieldElement.
    for i in range(1,len(bkpt)):
        p += plot(lambda x: bkpt[i]-x, (x, 0, bkpt[i]), color='grey', \
                  **ticks_keywords(function, True))
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
                        
def face_0D(face):
    if len(face[0]) == 1 and len(face[1]) == 1:
        return True
    else:
        return False
        
def face_2D(face):
    if len(face[0]) == 2 and len(face[1]) == 2 and len(face[2]) == 2:
        return True
    else:
        return False

def face_horizontal(face):
    if len(face[0]) == 2 and len(face[1]) == 1:
        return True
    else:
        return False

def face_vertical(face):
    if len(face[0]) == 1 and len(face[1]) == 2:
        return True
    else:
        return False                

def face_diagonal(face):
    if len(face[0]) == 2 and len(face[1]) == 2 and len(face[2]) == 1:
        return True
    else:
        return False 

def face_1D(face):
    return face_horizontal(face) or face_vertical(face) or face_diagonal(face)

@cached_function
def generate_vert_face_additive(function):
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
    vert_face_additive = []
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
                        
                        vert_face_additive.append(temp)
                        if i != j:
                            temp_swap = []
                            for vert in temp:
                                vert_new = [vert[1],vert[0]]
                                temp_swap.append(vert_new)
                            vert_face_additive.append(temp_swap)
                        
    logging.info("Computing maximal additive faces... done")
    return vert_face_additive  

@cached_function
def generate_minimal_triples(function):
    """
    Compute the minimal triples (projections) of the 
    (maximal) additive faces of the given function.
    """
    logging.info("Computing minimal triples representing maximal additive faces...")
    vert_face_additive = generate_vert_face_additive(function)
    minimal_triples = []
    for i in vert_face_additive:
        minimal_triples.append(projections(i))
    logging.info("Computing minimal triples representing maximal additive faces... done")
    return minimal_triples    

def triples_equal(a, b):
    return interval_equal(a[0], b[0]) and interval_equal(a[1], b[1]) and interval_equal(a[2], b[2])

            
### Create a new class representing a "face" (which knows its
### vertices, minimal triple, whether it's a translation/reflection,
### etc.; whether it's solid or dense).

### FIXME: Refactor some old code using it. 
        
class Face:
    def __init__(self, triple, vertices=None, is_known_to_be_minimal=False):
        if not vertices:
            vertices = verts(triple[0], triple[1], triple[2])
            if not vertices:
                raise NotImplementedError, "An empty face. This could mean need to shift triple[2] by (1,1). Not implemented."
        self.vertices = vertices
        i, j, k = projections(vertices)
        self.minimal_triple = minimal_triple = (i, j, interval_mod_1(k))
        if is_known_to_be_minimal and not triples_equal(minimal_triple, triple):
            raise ValueError, "Provided triple was not minimal: %s reduces to %s" % (triple, minimal_triple)

    def __repr__(self):
        return '<Face ' + repr(self.minimal_triple) + '>'

    def plot(self, *args, **kwds):
        return plot_face(self.minimal_triple, self.vertices, **kwds)

    def is_directed_move(self):
        return face_1D(self.minimal_triple) or face_0D(self.minimal_triple)
        
    def directed_move_with_domain_and_codomain(self):
        """
        Maps a horizontal edge to a forward translation,
        a vertical edge to a backward translation, 
        a diagonal edge to a reflection.
        """
        trip = self.minimal_triple
        if face_horizontal(trip):
            return (1, trip[1][0]), [trip[0]], [trip[2]]
        elif face_vertical(trip):
            return (1, -trip[0][0]), [trip[2]], [trip[1]]
        elif face_diagonal(trip):
            return (-1, trip[2][0]), [trip[0], trip[1]], [trip[1], trip[0]]
        else:
            raise ValueError, "Face does not correspond to a directed move: %s" % self

    def functional_directed_move(self):
        directed_move, domain, codomain = self.directed_move_with_domain_and_codomain()
        return FunctionalDirectedMove(domain, directed_move)

    
def plot_faces(faces, **kwds):
    p = Graphics()
    for f in faces:
        if type(f) == list or type(f) == tuple: #legacy
            f = Face(f)
        p += f.plot(**kwds)
    return p

@cached_function
def generate_maximal_additive_faces(function):
    return [ Face(trip, vertices=vertices, is_known_to_be_minimal=True) for trip, vertices in itertools.izip(generate_minimal_triples(function),generate_vert_face_additive(function)) ]

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

def plot_face(trip, vert, rgbcolor=(0,1,0), fill_color="mediumspringgreen", **kwds):
    if face_0D(trip):
        return point((trip[0][0], \
            trip[1][0]), rgbcolor = rgbcolor, size = 30)
    elif face_horizontal(trip):
        return parametric_plot((y,trip[1][0]),\
                             (y,trip[0][0], trip[0][1]), rgbcolor = rgbcolor, thickness=2)
    elif face_vertical(trip):
        return parametric_plot((trip[0][0],y),\
                             (y,trip[1][0], trip[1][1]), rgbcolor = rgbcolor, thickness=2)
    elif face_diagonal(trip):
        return parametric_plot((lambda y: y, lambda y: trip[2][0]-y),\
                             (y,trip[0][0],trip[0][1]), rgbcolor = rgbcolor, thickness=2)
    elif face_2D(trip):
        ## Sorting is necessary for this example:
        ## plot_2d_diagram(lift(piecewise_function_from_robert_txt_file("/Users/mkoeppe/w/papers/basu-hildebrand-koeppe-papers/reu-2013/Sage_Function/dey-richard-not-extreme.txt"))
        return polygon(convex_vert_list(vert), color=fill_color) 

def plot_2d_diagram(function, show_function = False):
    """
    Return a plot of the 2d complex with shaded faces where delta_pi is 0.        
    To show only a part of it, use 
    `show(plot_2d_diagram(h), xmin=0.25, xmax=0.35, ymin=0.25, ymax=0.35)`
    """
    bkpt = function.end_points()
    vert_face_additive = generate_vert_face_additive(function)
    minimal_triples = generate_minimal_triples(function)
        
    y = var('y')

    p = plot_2d_complex(function)
    for trip, vert in itertools.izip(minimal_triples, vert_face_additive):
        p += plot_face(trip, vert)
    ### For non-subadditive functions, show the points where delta_pi
    ### is negative.
    nonsubadditive_vertices = generate_nonsubadditive_vertices(function)
    p += point(nonsubadditive_vertices, \
                                  color = "red", size = 50)
    p += point([ (y,x) for (x,y) in nonsubadditive_vertices ], \
                                  color = "red", size = 50)
    if show_function:
        x = var('x')
        p += parametric_plot((lambda x: x, lambda x: 0.3 * float(function(x)) + 1), \
                                                (x, 0, 1), color='black')
        p += parametric_plot((lambda x: - 0.3 * float(function(x)), lambda x: x), \
                                                (x, 0, 1), color='black')
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
    minimal_triples = generate_minimal_triples(function)
            
    covered_intervals = []      
    # face = (I,J,K)
    for face in minimal_triples:
        if face_2D(face):
            component = []
            for int1 in face:
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

    edges = [ face for face in minimal_triples if face_1D(face) ]

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
    # the last test currently fails!
    bracketed_list = [[interval[0],interval[0]]] + remove_list + [[interval[1],interval[1]]]
    difference = []
    for i in range(len(bracketed_list) - 1):
        if (bracketed_list[i][1] < bracketed_list[i+1][0]): 
            difference.append([bracketed_list[i][1], bracketed_list[i+1][0]])
    return difference

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
        graph += plot(lambda x: function(x), \
                      [0,1], \
                      color = "black", legend_label="not covered", \
                      **kwds)
        kwds = {}
    for i, component in enumerate(covered_intervals):
        kwds.update({'legend_label': "covered component %s" % (i+1)})
        for interval in component:
            graph += plot(lambda x: function(x), interval, color=colors[i], **kwds)
            if 'legend_label' in kwds:
                del kwds['legend_label']
            if 'ticks' in kwds:
                del kwds['ticks']
            if 'tick_formatter' in kwds:
                del kwds['tick_formatter']
    return graph

### Minimality check.



# Fix x and y.
def type1check(fn):
    k=1
    bkpt = fn.end_points()
    for i in bkpt:
        for j in bkpt:
            if delta_pi(fn,i,j)<0:
                print "For x = ", i, ", y = ", j, ", x+y = ", i+j
                print "    Delta pi(x,y) = ", fn(fractional(i)), "+", fn(fractional(j)), "-", fn(fractional(i+j)), " = ", \
                    delta_pi(fn,i,j), " < 0"
                k=0
    if k==1:
        return True
    return False

# Fix x and x+y.
# By symmetry, the case in which y and x+y are fixed is also done. 
def type2check(fn):
    bkpt = fn.end_points()
    bkpt2 = []
    for i in range(len(bkpt)-1):
        bkpt2.append(bkpt[i])
    for i in range(len(bkpt)):
        bkpt2.append(bkpt[i]+1)     

    k=1
    
    for i in bkpt:
        for j in bkpt2:
            if j - i > 0 and delta_pi(fn,i,j-i)<0:
                print "For x = ", i, ", y = ", j-i, ", x+y = ",j
                print "    Delta pi(x,y) = ", fn(fractional(i)), "+", fn(fractional(j-i)), "-", fn(fractional(j)), " = ", \
                    delta_pi(fn,i,j-i), " < 0"
                k=0
    if k==1:
        return True
    return False     


def subadditivity_check(fn):
    """
    Check if a function is subadditive.
    Could take quite a while. (O(n^2))
    """
    if type1check(fn) and type2check(fn):
        print "pi is subadditive!"
        return True
    else:
        print "pi is not subadditive!"
        return False

def symmetric_test(fn, f):
    k = 1
    if fn(f) != 1:
        print 'pi(f) is not equal to 1. pi is not symmetric.'
        return False
    else:
        bkpt = fn.end_points()
        for i in bkpt:
            for j in bkpt:
                if i + j == f or i + j == 1 + f:
                    if delta_pi(fn, i, j) != 0:
                        print 'For x = ',i,'; ','y = ',j
                        print '    Delta pi is equal to ',delta_pi(fn, i, j),',not equal to 1'
                        k = 0
    if k == 1:
        print 'pi is symmetric.'
        return True
    return False

@cached_function
def find_f(fn):
    """
    Find the value of `f' for the given function `fn'.
    """
    for x in fn.end_points():
        if fn(x) > 1 or fn(x) < 0: 
            raise ValueError, "The given function does not stay in the range of [0, 1]"
    for x in fn.end_points():
        if fn(x) == 1:
            return x
    raise ValueError, "The given function has no breakpoint where the function takes value 1."

def minimality_test(fn, f=None):
    """
    Check if function `fn' is minimal with respect to the given `f'.
    Print diagnostics and return a boolean.
    """
    if f==None:
        f = find_f(fn)
    if fn(0) != 0:
        print 'pi is not minimal because pi(0) is not equal to 0.'
        return False
    else:
        print 'pi(0) = 0'
        if subadditivity_check(fn) and symmetric_test(fn, f):
            print 'Thus pi is minimal.'
            return True
        else:
            print 'Thus pi is not minimal.'
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
    def __init__(self, list_of_pairs, var=None):
        # Ensure sorted
        list_of_pairs = sorted(list_of_pairs, key = lambda ((a,b), f): a)
        # If adjacent functions are the same, just merge the pieces
        merged_list_of_pairs = []
        merged_interval_a = None
        merged_interval_b = None
        last_f = None
        for (a,b), f in list_of_pairs:
            #print f, last_f, f != last_f
            if f != last_f or (last_f != None and merged_interval_b < a):
                # Different function or a gap in the domain,
                # so push out the accumulated merged interval
                if last_f != None:
                    merged_list_of_pairs.append(((merged_interval_a, merged_interval_b), last_f))
                last_f = f
                merged_interval_a = a
                merged_interval_b = b
            else:
                if last_f != None:
                    merged_interval_b = max(b, merged_interval_b)
                else:
                    merged_interval_b = b
        if last_f != None:
            merged_list_of_pairs.append(((merged_interval_a, merged_interval_b), last_f))
            
        PiecewisePolynomial.__init__(self, merged_list_of_pairs, var)

        intervals = self._intervals
        functions = self._functions
        end_points = [ intervals[0][0] ] + [b for a,b in intervals]
        self._end_points = end_points
        values_at_end_points = [ functions[0](end_points[0]) ]
        for i in range(len(functions)):
            value = functions[i](intervals[i][1])
            values_at_end_points.append(value)
        self._values_at_end_points = values_at_end_points


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
            sage: f2(x) = 1-x
            sage: f3(x) = x^2-5
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3]])
            sage: f.end_points()
            [0, 1, 2, 3]
        """
        return self._end_points

    @cached_method
    def __call__(self,x0):
        """
        Evaluates self at x0. Returns the average value of the jump if x0
        is an interior endpoint of one of the intervals of self and the
        usual value otherwise.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            sage: f(0.5)
            1
            sage: f(5/2)
            e^(5/2)
            sage: f(5/2).n()
            12.1824939607035
            sage: f(1)
            1/2
        """
        #print "EVAL at ", x0
        endpts = self.end_points()
        i = bisect_left(endpts, x0)
        if i >= len(endpts):
            raise ValueError,"Value not defined at point %s, outside of domain." % x0
        if x0 == endpts[i]:
            return self._values_at_end_points[i]
        if i == 0:
            raise ValueError,"Value not defined at point %s, outside of domain." % x0
        if self._intervals[i-1][0] <= x0 < self._intervals[i-1][1]:
            return self.functions()[i-1](x0)
        raise ValueError,"Value not defined at point %s, outside of domain." % x0

    def __add__(self,other):
        F, G, intervals = self._make_compatible(other)
        fcn = []
        for a,b in intervals:
            fcn.append([(a,b), F.which_function(b)+G.which_function(b)])        
        return FastPiecewise(fcn)
        
    def __mul__(self,other):
        if not isinstance(other, FastPiecewise):
            # assume scalar multiplication
            return FastPiecewise([[(a,b), other*f] for (a,b),f in self.list()])
        else:
            F, G, intervals = self._make_compatible(other)
            fcn = []
            for a,b in intervals:
                fcn.append([(a,b),F.which_function(b)*G.which_function(b)])     
            return FastPiecewise(fcn)

    __rmul__ = __mul__

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
        for ((a,b), f) in self.list():
            if xmin is not None:
                a = max(a, xmin)
            if xmax is not None:
                b = min(b, xmax)
            if a < b:
                # We do not plot anything if a==b because
                # otherwise plot complains that
                # "start point and end point must be different"
                g += plot(f, *args, xmin=a, xmax=b, **kwds)
                # If it's the first piece, pass all arguments. Otherwise,
                # filter out 'legend_label' so that we don't add each
                # piece to the legend separately (trac #12651).
                if 'legend_label' in kwds:
                    del kwds['legend_label']
            elif a == b:
                g += point([a, f(a)])
        return g

    def __repr__(self):
        rep = "<FastPiecewise with %s parts, " % len(self._functions)
        for interval, function in itertools.izip(self._intervals, self._functions):
            rep += "\n " + repr(interval) + "\t" + repr(function) \
                   + "\t values: " + repr([function(interval[0]), function(interval[1])])
        rep += ">"
        return rep

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
    sage: emb_field = RealNumberField(field.polynomial(), 'a', embedding=morphism(field.gen(0)))
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
        NumberField_absolute.__init__(self, polynomial, name=name, check=check,
                                      embedding=embedding, latex_name=latex_name,
                                      assume_disc_small=assume_disc_small, maximize_at_primes=maximize_at_primes)
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
        self._zero_element = self(0)
        self._one_element =  self(1)
    
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
            global number_field, number_field_values, morphism, exact_generator, embedded_field, embedding_field, hom, embedded_field_values
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
    global symb_values
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
        bkpt.append(bkpt[i]+interval_lengths[i])
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field)

def limiting_slopes(fn):
    functions = fn.functions()
    return functions[0]._slope, functions[-1]._slope

maximal_asymmetric_peaks_around_orbit = 'maximal_asymmetric_peaks_around_orbit'
maximal_symmetric_peaks_around_orbit = 'maximal_symmetric_peaks_around_orbit'
recentered_symmetric_peaks = 'recentered_symmetric_peaks'
recentered_peaks_with_slopes_proportional_to_limiting_slopes_for_positive_epsilon = 'recentered_peaks_with_slopes_proportional_to_limiting_slopes_for_positive_epsilon'

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
    if perturbation_style==maximal_asymmetric_peaks_around_orbit or perturbation_style==maximal_symmetric_peaks_around_orbit:
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

# Tree of closed intervals [a, b] to speed up Delta pi calculations

# Type1: We check x breakpoint, y breakpoint against x + y arbitrary point.
# So use discrete branching on x, y.

# Persistent tree: In one variable.
# Saves maximum and minimum information of the function over nodes of tree.

tree_threshold = 0  # Means allow branch down to singleton intervals

class IntervalNode:
    # An interval can be broken into three subintervals:
    # [a, v] (v, w) [w, b] 
    # where a, v, w, b are breakpoints 
    # and (v, w) is a subinterval that contains no breakpoint.
    # The indices a_index and b_index are the indices of the endpoints of the function
    # that are a and b
    
    def __init__(self, fn, a, b, a_index, b_index):
        self.fn = fn
        self.a = a
        self.b = b
        self.a_index = a_index
        self.b_index = b_index
        self.end_node = True
        self.min_max()

    def __repr__(self):
        return self._repr_with_level_(0)

    def _repr_with_level_(self, level, recurse = False):
        rep = "IntervalNode a: %s, b: %s, a_index: %s, b_index: %s, min: %s, max: %s" \
              % (self.a, self.b, self.a_index, self.b_index, self._min, self._max)
        if self.end_node:
            return "  " * level + "<" + rep + ">"
        else:
            rep += ", v: %s, w: %s " % (self.v, self.w)
            rep = "  " * level + "<" + rep 
            if recurse:
                rep += "\n" \
                       + self.left._repr_with_level_(level + 1, recurse=recurse) + "\n" \
                       + self.right._repr_with_level_(level + 1, recurse=recurse)
            rep += ">"
            return rep

    def print_tree(self):
        print self._repr_with_level_(0, True)

    def children(self):
        return self.left, self.right

    def min_max(self):
        mini = getattr(self, '_min', None)
        maxi = getattr(self, '_max', None)
        if mini is None:
            fn = self.fn
            end_points = fn.end_points()
            a_index = self.a_index
            b_index = self.b_index
            a = self.a
            b = self.b
            if b_index - a_index > tree_threshold:
                # Branch
                v_index = (b_index + a_index) // 2
                w_index = v_index + 1
                assert a_index <= v_index < w_index <= b_index
                v = self.v = end_points[v_index]
                w = self.w = end_points[w_index]
                self.left = IntervalNode(fn, a, v, a_index, v_index)
                self.right = IntervalNode(fn, w, b, w_index, b_index)
                left_min, left_max = self.left.min_max()
                right_min, right_max = self.right.min_max()
                mini = min(left_min, right_min)
                maxi = max(left_max, right_max)
                self.end_node = False
            else:
                mini = maxi = fn(end_points[self.b_index])
                for i in range(self.a_index, self.b_index):
                    value = fn(end_points[i])
                    if value < mini:
                        mini = value
                    elif value > maxi:
                        maxi = value
            self._min = mini
            self._max = maxi
        return mini, maxi
    
    def min_max_on_interval(self, c, d):
        ## later we want to bail out earlier based on a given range as well.
        #print "Enter node ", self
        c = max(self.a, c)
        d = min(self.b, d)

        if c > d:
            return None, None
        if c == self.a and d == self.b:
            return self._min, self._max

        minmax = [None, None]

        def update_minmax(small, big):
            #print "update: from ", minmax, "with ", small, big,  
            if small is not None and (minmax[0] is None or small < minmax[0]):
                minmax[0] = small
            if big is not None and (minmax[1] is None or big > minmax[1]):
                minmax[1] = big
            #print "to ", minmax 

        fn = self.fn
        if self.end_node:
            # Test at c and d, which may or may not be endpoints of the function.
            fnc = fn(c)
            update_minmax(fnc, fnc)
            fnd = fn(d)
            update_minmax(fnd, fnd)
            end_points = fn.end_points()
            for i in range(self.a_index, self.b_index):
                if c <= end_points[i] <= d:
                    fnx = fn(end_points[i])
                    update_minmax(fnx, fnx)
        else:    
            left, right = self.children()
            if c <= self.v:
                left_min, left_max = left.min_max_on_interval(c, d)
                update_minmax(left_min, left_max)
            if self.w <= d:
                right_min, right_max = right.min_max_on_interval(c, d)
                update_minmax(right_min, right_max)
            c = max(self.v, c)
            d = min(self.w, d)
            if c <= d:
                if self.v < c:
                    fnc = fn(c)
                    update_minmax(fnc, fnc)
                if c < d < self.w:
                    fnd = fn(d)
                    update_minmax(fnd, fnd)
        #print "min_max: ", minmax
        return minmax[0], minmax[1]

def generate_min_max_tree(fn):
    return IntervalNode(fn, 0, 1, 0, len(fn.end_points()) - 1)

# Branching tree: In two variables.


## class BranchNode:
    
##     def variable_intervals(self):
##         return self._variable_intervals

    

## def branch



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
        # FIXME: This does not do complete
        # error checking.
        directed_move = self.directed_move
        move_sign = directed_move[0]
        if move_sign == 1:
            if inverse:
                result = (interval[0] - directed_move[1], interval[1] - directed_move[1])
            else:
                result = (interval[0] + directed_move[1], interval[1] + directed_move[1])
        elif move_sign == -1:
            result = (directed_move[1] - interval[0], directed_move[1] - interval[1])
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
            result = closed_or_open_or_halfopen_interval(directed_move[1] - interval[0], directed_move[1] - interval[1], \
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
            fdm = face.functional_directed_move()
            if not fdm.is_identity() and find_interior_intersection(fdm.intervals(), intervals): #FIXME: why interior?
                moves.add(fdm)
    return list(moves)

def is_directed_move_possible(x, move):
    return move.can_apply(x)

def find_possible_directed_moves(x,directed_moves):
    """Find the directed moves applicable at x.
    """
    possible_directed_moves = []
    for directed_move in directed_moves:
        if is_directed_move_possible(x, directed_move):
            possible_directed_moves.append(directed_move)
    return possible_directed_moves

def find_impossible_directed_moves(x,directed_moves):
    """
    Find the directed moves NOT applicable at x.
    """
    impossible_directed_moves = []
    for directed_move in directed_moves:
        if not is_directed_move_possible(x, directed_move):
            impossible_directed_moves.append(directed_move)
    return impossible_directed_moves

def apply_directed_move(x, directed_move):
    return directed_move(x)

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



def random_walk(seed, moves, fn, num_it):
    """
    Find the orbit of a seed through a random walk.  It may not find
    all the possible elements in the orbit (if unlucky).  It may not
    be very efficient if the size of orbit is small because it will
    continue to run until the predetermined number of iterations is
    reached.

    deterministic_walk is preferable in most cases.
    """
    import random
    x = seed
    xlist = {x:[1,None,None]}
    walk_sign = 1
    # points_plot = point((x,0))
    for i in range(num_it):
        possible_directed_moves = find_possible_directed_moves(x,moves)
        move = possible_directed_moves[random.randint(0,len(possible_directed_moves)-1)]
        move_sign = move.sign()
        directed_move = move
        
        next_x = move(x)

        next_walk_sign = move_sign * walk_sign
        if next_x not in xlist:
            pass
            # points_plot += point((next_x,0))
        else:
            if xlist[next_x][2] != None:
                x = xlist[next_x][1]
                directed_move = xlist[next_x][2]
            if xlist[next_x][0] != next_walk_sign:
                next_walk_sign = 0
                print next_x, ","
                previous = xlist[next_x][1]
                while previous != next_x and previous != None:
                    print previous, ","
                    previous = xlist[previous][1]
                print next_x, "."
        # Prep for next
        xlist[next_x] = [next_walk_sign, x, directed_move] 
        walk_sign = next_walk_sign
        x = next_x
    # show(points_plot + plot(-0.1+x-x, [A,A0]) + plot(-0.1+x-x, [f-A0,f-A])\
    #     + plot(-0.2+x-x, [A,A+a0], color = "green") + plot(-0.3+x-x, [A,A+a1], color = "green") + \
    #     plot(-0.4+x-x, [A,A+a2], color = "green"), ymin = -1, ymax = 0)
    return xlist

class MaximumNumberOfIterationsReached(Exception):
    pass

class SignContradiction(Exception):
    pass

def deterministic_walk(seed, moves, fn=None, max_num_it = 1000, error_if_sign_contradiction=False, error_if_max_num_it_exceeded=True):
    """
    Compute the orbit of a given seed. (Done by a breadth-first search.)
    To avoid infinite computations in the case of a dense orbit,
    there is a maximum number of iterations (by default, it is 1000).

    The moves know the domains on which they are allowed.

    Returns a dictionary:
    - keys are elements of the orbit
    - values are lists of the form [walk_sign, predecessor, directed_move_from_predecessor].

    If `error_if_max_num_it_exceeded` is `False`, then 
    a secondary result is the to_do list.
    """
    ## FIXME: If `fn` is provided, store the dictionary in `fn.walk_dict`
    ## and the to_do list in `fn.walk_to_do`, to allow us to resume
    ## an interrupted BFS.  OR: Turn it into an iterator/generator (yield).

    logging.info("Breadth-first search to discover the reachable orbit...")
    to_do = [seed]
    # xlist is actually a dictionary.
    
    xlist = {seed:[1,None,None]}
    walk_sign = 1
    # points_plot = point((seed,0))
    contradiction_reached = False
    num_it = 0
    while to_do and num_it < max_num_it:
        if (num_it > 0 and num_it % 100 == 0):
            logging.info("(Iteration %d, to do list has %d items)" % (num_it, len(to_do)))
        x = to_do.pop(0)
        possible_directed_moves = find_possible_directed_moves(x, moves)
        for directed_move in possible_directed_moves:
            walk_sign = xlist[x][0]
            move_sign = directed_move.sign()
            move_for_next_x = directed_move
            next_x = directed_move(x)
            
            next_walk_sign = move_sign * walk_sign
            
            if next_x not in xlist:
                # points_plot += point((next_x,0))
                xlist[next_x] = [next_walk_sign, x, move_for_next_x]
                #if next_x not in to_do:
                to_do.append(next_x)
                
            else:
                # If next_x is already in xlist, we do not need to modify anything, except probably the sign.
                if xlist[next_x][0] != next_walk_sign:
                    xlist[next_x][0] = 0
                    if contradiction_reached == False:
                        logging.info('A contradiction of signs was reached. All the elements in the reachable orbit take the value 0.')
                        if error_if_sign_contradiction:
                            raise SignContradiction
                        else:
                            contradiction_reached = True
                            for pt in xlist.keys():
                                xlist[pt][0] = 0
        num_it = num_it + 1

    if error_if_max_num_it_exceeded:
        if num_it == max_num_it:
            raise MaximumNumberOfIterationsReached, "Reached %d iterations, to do list has still %d items" % (num_it, len(to_do))
        logging.info("Breadth-first search to discover the reachable orbit... done")
        return xlist
    else:
        if num_it == max_num_it:
            logging.info("Breadth-first search to discover the reachable orbit reached %d iterations, to do list has still %d items" % (num_it, len(to_do)))
        return xlist, to_do

def plot_walk(walk_dict, color="black", ymin=0, ymax=1, **options):
    #return point([ (x,0) for x in walk_dict.keys()])
    g = Graphics()
    for x in walk_dict.keys():
        g += line([(x,ymin), (x,ymax)], color=color, zorder = -4, **options)
    return g

def plot_shifted_stability_intervals(shifted_stability_intervals, color="cyan", ymin=0, ymax=1, **options):
    g = Graphics()
    for i in shifted_stability_intervals:
        g += polygon2d([(i.a, ymin), (i.b, ymin), (i.b, ymax), (i.a, ymax)], \
                       color=color, zorder = -5, **options)
    return g

def plot_walk_with_shifted_stability_intervals(stab_int, walk_dict, color="black", stab_color="cyan", **options):
    shifted_stability_intervals = generate_shifted_stability_intervals(stab_int, walk_dict)
    return plot_shifted_stability_intervals(shifted_stability_intervals, color=stab_color, **options) \
        + plot_walk(walk_dict, color=color, **options)

def plot_intervals(intervals, ymin=0, ymax=1):
    g = Graphics()
    for interval in intervals:
        g += polygon([(interval[0], ymin), (interval[1], ymin), \
                      (interval[1], ymax), (interval[0], ymax)], 
                     color="yellow", zorder = -8)
    return g

def plot_orbit_comparison(fn, orbits):
    """
    `orbits` can be a list of:
    tuples (stab_int, walk_dict)
    or numbers (seeds)
    """
    g = plot_intervals(generate_uncovered_intervals(fn))
    ymin = 0.5
    ymax = 1
    for orbit, color in itertools.izip(orbits, rainbow(len(orbits), 'rgbtuple')):
        ymax -= 0.2
        ymin -= 0.2
        stab_color = Color(color).lighter() 
        try:
            stab_int, walk_dict = orbit
        except TypeError:
            seed = orbit
            stab_int, walk_dict = find_stability_interval_with_deterministic_walk_list(seed, generate_uncovered_intervals(fn), generate_functional_directed_moves(fn), fn)
        g += plot_walk_with_shifted_stability_intervals(stab_int, walk_dict, stab_color=stab_color, color=color, ymin=ymin, ymax=ymax)
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

def one_step_stability_interval(x, intervals, moves):
    """Returns the stability interval, i.e., an open, half-open, or closed
    interval I), such that for all moves, `move`(`x`) lies in
    `intervals` if and only if `move`(`x` + t) lies in `intervals` for
    t in I."""
    a = -10
    b = 10
    left_closed = True
    right_closed = True
    for move in moves:
        next_x = directed_move(x)
        for interval in intervals:
            if (interval[0] <= next_x <= interval[1]):
                # Move was possible, so:
                if move[0] == 1:
                    if (interval[0] - next_x) > a:
                        a = interval[0] - next_x
                        left_closed = True
                    if (interval[1] - next_x) < b:
                        b = interval[1] - next_x
                        right_closed = True
                elif move[0] == -1:
                    if (next_x - interval[0]) < b:
                        b = next_x - interval[0]
                        right_closed = True
                    if (interval[1] - next_x) < -a:
                        a = next_x - interval[1]
                        left_closed = True      
                else:
                    raise ValueError, "Move not valid: %s" % list(move)
                break # next_x can only lie in one interval 
        else:
            # Move was not possible.
            if move[0] == 1:
                for interval in intervals:
                    temp = interval[0] - next_x
                    if 0 < temp <= b:
                        b = temp
                        right_closed = False
                    temp2 = interval[1] - next_x
                    if a <= temp2 < 0:
                        a = temp2
                        left_closed = False
            elif move[0] == -1:
                for interval in intervals:
                    temp = interval[0] - next_x
                    if 0 < temp <= -a:
                        a = -temp
                        left_closed = False
                    temp2 = next_x - interval[1]
                    if 0 < temp2 <= b:
                        b = temp2
                        right_closed = False
            else:
                raise ValueError, "Move not valid: %s" % list(move)
    return closed_or_open_or_halfopen_interval(a, b, left_closed, right_closed)

# def closed_or_open_or_halfopen_interval_intersection(int1,int2):
#     if max(int1[0],int2[0]) > min(int1[1],int2[1]):
#         return []
#     if max(int1[0],int2[0]) == min(int1[1],int2[1]):
#         return [max(int1[0],int2[0])]
#     else:
#         return [max(int1[0],int2[0]),min(int1[1],int2[1])]

def one_step_stability_refinement(interval, intervals, moves):
    if len(interval) == 2:
        interval = closed_or_open_or_halfopen_interval(interval[0], interval[1], True, True)
    seed = (interval[0] + interval[1]) / 2
    (stab_a, stab_b, stab_left_closed, stab_right_closed) = one_step_stability_interval(seed, intervals, moves)
    stab_a += seed
    stab_b += seed
    left = None
    right = None
    #print "  Shifted stability interval: ", stab_a, stab_b, stab_left_closed, stab_right_closed
    #print interval.a < stab_a, interval.left_closed
    if interval.a < stab_a \
       or (interval.a == stab_a and interval.left_closed and not stab_left_closed):
        left = closed_or_open_or_halfopen_interval(interval.a, stab_a, interval.left_closed, not stab_left_closed)
        #print "  Left: ", left
    else:
        stab_a = interval.a
        stab_left_closed = interval.left_closed
    if stab_b < interval.b \
       or (stab_b == interval.b and interval.right_closed and not stab_right_closed):
        right = closed_or_open_or_halfopen_interval(stab_b, interval.b, not stab_right_closed, interval.right_closed)
        #print "  Right: ", right
    else:
        stab_b = interval.b
        stab_right_closed = interval.right_closed
    #print "  Modified stability interval: ", stab_a, stab_b, stab_left_closed, stab_right_closed
    refinement = []
    if left != None:
        refinement.append(left)
    refinement.append(closed_or_open_or_halfopen_interval(stab_a, stab_b, stab_left_closed, stab_right_closed))
    if right != None:
        refinement.append(right)
    return refinement

def interval_length(interval):
    if interval[1] >= interval[0]:
        return interval[1] - interval[0]
    return 0

from heapq import *

def iterative_stability_refinement(intervals, moves):
    interval_pq = [ (-interval_length(interval), interval) \
                    for interval in intervals ]
    heapify(interval_pq)
    finished_list = []
    while (interval_pq):
        priority, int = heappop(interval_pq)
        length = -priority
        logging.info("%s of length %s " % (int, float(length)))
        refinement = one_step_stability_refinement(int, intervals, moves)
        logging.debug("  Refinement: %s" % refinement)
        if len(refinement) == 1:
            finished_list.append(refinement[0])
        else:
            for new_int in refinement:
                heappush(interval_pq, (-interval_length(new_int), new_int))
    return sorted(finished_list + \
                  [ interval for (priority, interval) in interval_pq ])

def find_decomposition_into_stability_intervals(fn):
    ## experimental.
    intervals = generate_uncovered_intervals(fn)
    moves = generate_functional_directed_moves(fn)
    return iterative_stability_refinement(intervals, moves)

def scan_coho_interval_list(interval_list, tag=None):
    """Generate events of the form `(x, epsilon), delta, tag.`"""
    for i in interval_list:
        if len(i) == 2:
            # old-fashioned closed interval
            yield (i[0], 0), +1, tag                             # Turn on at left endpoint
            yield (i[1], 1), -1, tag                             # Turn off at right endpoint plus epsilon
        else:
            # coho interval
            yield (i.a, 0 if i.left_closed else 1), +1, tag
            yield (i.b, 1 if i.right_closed else 0), -1, tag

## def scan_set_difference(a, b):
##     """`a` and `b` should be event generators."""

def scan_union_of_coho_interval_minus_union_of_coho_intervals(interval_list, remove_list):
    # Following uses the lexicographic comparison of the tuples.
    scan = merge(scan_coho_interval_list(interval_list, True),
                 scan_coho_interval_list(remove_list, False))
    interval_indicator = 0
    remove_indicator = 0
    on = False
    for ((x, epsilon), delta, tag) in scan:
        was_on = on
        if tag:                                       # interval event
            interval_indicator += delta
            assert(interval_indicator) >= 0
        else:                                           # remove event
            remove_indicator += delta
            assert(remove_indicator) >= 0
        now_on = interval_indicator > 0 and remove_indicator == 0
        if not was_on and now_on: # switched on
            yield (x, epsilon), +1, None
        elif was_on and not now_on: # switched off
            yield (x, epsilon), -1, None
        on = now_on
    # No unbounded intervals:
    assert interval_indicator == 0
    assert remove_indicator == 0

def coho_interval_list_from_scan(scan):
    """Actually returns a generator."""
    indicator = 0
    (on_x, on_epsilon) = (None, None)
    for ((x, epsilon), delta, tag) in scan:
        was_on = indicator > 0
        indicator += delta
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

def union_of_coho_interval_minus_union_of_coho_intervals(interval_list, remove_list):
    """Compute a list of closed/open/half-open intervals that represent
    the set difference of `interval` and the union of the intervals in
    `remove_list`.

    Assumes `interval_list' and `remove_list` are both sorted (and
    each pairwise disjoint), and returns a sorted list.

    EXAMPLES::
    sage: union_of_coho_interval_minus_union_of_coho_intervals([[0,10]], [[2,2], [3,4]])
    [<Int[0, 2)>, <Int(2, 3)>, <Int(4, 10]>]
    """
    gen = coho_interval_list_from_scan(scan_union_of_coho_interval_minus_union_of_coho_intervals(interval_list, remove_list))
    return [ int for int in gen ]

def find_decomposition_into_stability_intervals(fn, show_plots=False, max_num_it=1000):
    fn._stability_orbits = []
    uncovered_intervals = generate_uncovered_intervals(fn)
    intervals = uncovered_intervals
    moves = generate_functional_directed_moves(fn)
    orbits = []
    while intervals:
        #print "Intervals: ", intervals
        seed = intervals[0][0] + (intervals[0][1] - intervals[0][0]) / 3
        print "Seed: ", seed
        walk_dict, to_do = deterministic_walk(seed, moves, fn, max_num_it = max_num_it, \
                                              error_if_sign_contradiction=False, 
                                              error_if_max_num_it_exceeded=False)
        # When to_do is nonempty, the BFS did not finish.
        # We compute the stability interval anyway -- it gives an interval of points that ALL have such a long orbit!
        # But we don't try to make the stability orbits disjoint in this case (this is not well-defined);
        # instead we count on our scan functions to deal with the "multiset" case properly.
        make_disjoint = not to_do
        stab_int = compute_stability_interval(walk_dict, uncovered_intervals, moves, fn, make_disjoint = make_disjoint)
        orbits.append((stab_int, walk_dict))
        if show_plots:
            plot_orbit_comparison(fn, orbits).show(dpi=500)
        #print "Stability interval: ", stab_int
        shifted_stability_intervals = generate_shifted_stability_intervals(stab_int, walk_dict)
        print "Stability orbit: ", shifted_stability_intervals[0], ", ... (length ", len(walk_dict), ")"
        fn._stability_orbits.append((shifted_stability_intervals, walk_dict, to_do))
        remaining = union_of_coho_interval_minus_union_of_coho_intervals(intervals, shifted_stability_intervals)
        intervals = remaining
        
    logging.info("Total: %s stability orbits, lengths: %s" \
                 % (len(fn._stability_orbits), \
                    [ ("%s+" if to_do else "%s") % len(shifted_stability_intervals) \
                      for (shifted_stability_intervals, walk_dict, to_do) in fn._stability_orbits ]))

def compute_stability_interval(deterministic_walk_list, intervals, moves, fn, make_disjoint = True):
    ## FIXME: Refactor using above.
    a = -10
    b = 10
    left_closed = True
    right_closed = True
    for pt in deterministic_walk_list.keys():
        if deterministic_walk_list[pt][0] == 0:
            return closed_or_open_or_halfopen_interval(0, 0, True, True)
        for interval in intervals:
            if element_of_int(pt, interval):
                if deterministic_walk_list[pt][0] == 1:
                    if (interval[0] - pt) > a:
                        a = interval[0] - pt
                        left_closed = True
                    if (interval[1] - pt) < b:
                        b = interval[1] - pt
                        right_closed = True
                elif deterministic_walk_list[pt][0] == -1:
                    if (pt - interval[0]) < b:
                        b = pt - interval[0]
                        right_closed = True
                    if (interval[1] - pt) < -a:
                        a = pt - interval[1]
                        left_closed = True      
        impossible_directed_moves = find_impossible_directed_moves(pt, moves)
        ### We now take the set difference of
        ### the __directed__ moves with the possible directed moves.
        ### This was not done in the old code: --Matthias
        # impossible_directed_moves = []
        # for move in moves:
        #     if move not in possible_directed_moves:
        #         impossible_directed_moves.append(move)
        
        for move in impossible_directed_moves:
            if move.sign() == 1:
                impossible_next_x = fractional(pt + move[1])
                for interval in intervals:
                    temp = interval[0] - impossible_next_x
                    if 0 < temp <= b:
                        b = temp
                        right_closed = False
                    temp2 = interval[1] - impossible_next_x
                    if a <= temp2 < 0:
                        a = temp2
                        left_closed = False
            elif move.sign() == -1:
                impossible_next_x = fractional(move[1] - pt)
                for interval in intervals:
                    temp = interval[0] - impossible_next_x
                    if 0 < temp <= -a:
                        a = -temp
                        left_closed = False
                    temp2 = impossible_next_x - interval[1]
                    if 0 < temp2 <= b:
                        b = temp2
                        right_closed = False
    if make_disjoint:
        ### Now we make the stability orbit disjoint by looking at adjacent intervals
        ### and shrinking if necessary.
        orbit = sorted(deterministic_walk_list.keys())
        min_midpt_dist = 1
        for i in range(len(orbit)-1):
            sign_1 = deterministic_walk_list[orbit[i]][0]
            sign_2 = deterministic_walk_list[orbit[i+1]][0]
            half_distance = (orbit[i+1]-orbit[i])/2
            ## Two intervals of different signs might overlap.
            if sign_1 == 1 and sign_2 == -1:
                if b > half_distance:
                    b = half_distance
                    logging.info("Making stability intervals disjoint: Reducing right to %s to separate %s and %s" % (b, orbit[i], orbit[i+1]))
                    right_closed = False
            elif sign_1 == -1 and sign_2 == 1:
                if -a > half_distance:
                    a = -half_distance
                    logging.info("Making stability intervals disjoint: Reducing left to %s to separate %s and %s" % (a, orbit[i], orbit[i+1]))
                    left_closed = False
            else:
                assert(-a + b <= 2 * half_distance)
    return closed_or_open_or_halfopen_interval(a, b, left_closed, right_closed)

def find_stability_interval_with_deterministic_walk_list(seed, intervals, moves, fn, max_num_it = 1000, error_if_sign_contradiction=False):
    """
    Returns the stability interval (an open, half-open, or closed interval)
    and the deterministic_walk_list.
    """
    walk_dict = deterministic_walk(seed, moves, fn, max_num_it, \
                                   error_if_sign_contradiction=error_if_sign_contradiction)
    return (compute_stability_interval(walk_dict, intervals, moves, fn), walk_dict)

def generate_shifted_stability_intervals(stab_int, walk_list):
    """The result is sorted."""
    orbit = sorted(walk_list.keys())
    intervals = []
    for i in orbit:
        sign = walk_list[i][0]
        if sign == 1:
            intervals.append(closed_or_open_or_halfopen_interval(stab_int.a + i, stab_int.b + i, stab_int.left_closed, stab_int.right_closed))
        elif sign == -1:
            intervals.append(closed_or_open_or_halfopen_interval(i - stab_int.b, i - stab_int.a, stab_int.right_closed, stab_int.left_closed))
        elif sign == 0:
            assert(stab_int.a == -stab_int.b)
            assert(stab_int.left_closed == stab_int.right_closed)
            intervals.append(closed_or_open_or_halfopen_interval(stab_int.a + i, stab_int.b + i, stab_int.left_closed, stab_int.right_closed))
    return intervals

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

def find_generic_seed(fn, max_num_it = 1000):
    intervals = generate_uncovered_intervals(fn)
    if not intervals:
        raise ValueError, "Need an uncovered interval"
    moves = generate_functional_directed_moves(fn)
    seed = intervals[0][1]
    while True:
        seed = 2/3 * intervals[0][0] + 1/3 * seed
        logging.info("Trying seed %s" % seed)
        try:
            stab_int, walk_list = find_stability_interval_with_deterministic_walk_list \
                                  (seed, intervals, moves, fn, max_num_it = max_num_it, \
                                   error_if_sign_contradiction = True)
            if not stab_int[0] < 0 < stab_int[1]:
                logging.info("Stability interval does not contain 0 in its interior, continuing search.")
                continue
            logging.info("Seed %s has a proper stability interval %s, reachable orbit has %s elements" % (seed, stab_int, len(walk_list)))
            return (seed, stab_int, walk_list)
        except SignContradiction:
            continue

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
    """Return a piecewise function compatible with the given `function`, 
    using symbolic slope coefficients, one for each component.
    As a second result, return a list of the symbolic slope coefficients."""
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
    slope_vars = list(var([ "slope%s" % (i+1) for i in range(len(components))]))
    symbolic = generate_compatible_piecewise_function(components, slope_vars, field=SR)
    return symbolic, slope_vars, components

def additivity_equation(symbolic, x, y):
    return delta_pi(symbolic, x, y) == 0

def generate_additivity_equations(function, symbolic):
    f = find_f(function)
    return uniq([additivity_equation(symbolic, x, y) \
                 for (x, y) in generate_additive_vertices(function) ] \
                + [symbolic(f) == 0] \
                + [symbolic(1) == 0])

def rescale_to_amplitude(perturb, amplitude):
    """For plotting purposes, rescale the function `perturb` so that its
    maximum absolute function value is `amplitude`.
    """
    current_amplitude = max([ abs(perturb(x)) for x in perturb.end_points() ])
    if current_amplitude != 0:
        return perturb * (amplitude/current_amplitude)
    else:
        return perturb

def check_perturbation(fn, perturb, show_plots=False, **show_kwds):
    epsilon_interval = fn._epsilon_interval = find_epsilon_interval(fn, perturb)
    epsilon = min(abs(epsilon_interval[0]), epsilon_interval[1])
    print "Epsilon for constructed perturbation: ", epsilon
    assert epsilon > 0, "Epsilon should be positive, something is wrong"
    print "Thus the function is not extreme."
    if show_plots:
        logging.info("Plotting perturbation...")
        (plot(fn, xmin=0, xmax=1, color='black', thickness=2, legend_label="original function") \
         + plot(fn + epsilon * perturb, xmin=0, xmax=1, color='blue', legend_label="+perturbed") \
         + plot(fn + (-epsilon) * perturb, xmin=0, xmax=1, color='red', legend_label="-perturbed") \
         + plot(rescale_to_amplitude(perturb, 1/10), xmin=0, xmax=1, color='magenta', legend_label="perturbation")) \
        .show(figsize=50, **show_kwds)
        logging.info("Plotting perturbation... done")

def finite_dimensional_extremality_test(function, show_plots=False):
    symbolic, slope_vars, components = symbolic_piecewise(function)
    ## FIXME: Get rid of maxima-based "solve" in favor of linear algebra code.
    equations = generate_additivity_equations(function, symbolic)
    logging.info("Equations: %s" % equations)
    solutions = solve(equations, slope_vars, solution_dict=True)
    logging.info("Solutions: %s" % solutions)
    assert len(solutions) == 1, "We need to see exactly one general solution."
    solution = solutions[0]
    rhs_variables = []
    for sol in solution.values():
        rhs_variables.extend(sol.variables())
    rhs_variables = uniq(rhs_variables)
    logging.info("Solution space has dimension %s" % len(rhs_variables))
    if len(rhs_variables) == 0:
        logging.info("Thus the function is extreme.")
        return True
    else:
        for basis_index in range(len(rhs_variables)):
            substitutions = {var: 1 if i == basis_index else 0 for i, var in enumerate(rhs_variables)}
            slopes = [ solution[slope_var].subs(substitutions) for slope_var in slope_vars ]
            perturbation = function._perturbation = generate_compatible_piecewise_function(components, slopes)
            check_perturbation(function, perturbation, show_plots=show_plots, legend_title="Basic perturbation %s" % (basis_index + 1))
    return False

@cached_function
def generate_additive_vertices(fn):
    bkpt = fn.end_points()
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]
    type1 = [ (x,y) for x in bkpt for y in bkpt if x <= y and delta_pi(fn,x,y) == 0 ]
    type2 = [ (x,z-x) for x in bkpt for z in bkpt2 if x < z < 1+x and delta_pi(fn, x, z-x) == 0]
    return uniq(type1 + type2)

def generate_nonsubadditive_vertices(fn):
    bkpt = fn.end_points()
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]
    type1 = [ (x,y) for x in bkpt for y in bkpt if x <= y and delta_pi(fn,x,y) < 0 ]
    type2 = [ (x,z-x) for x in bkpt for z in bkpt2 if x < z < 1+x and delta_pi(fn, x, z-x) < 0]
    return uniq(type1 + type2)

def extremality_test(fn, show_plots = False, max_num_it = 1000, perturbation_style=default_perturbation_style, phase_1 = False, finite_dimensional_test_first = False):
    do_phase_1_lifting = False
    if not minimality_test(fn):
        print "Not minimal, thus not extreme."
        if not phase_1:
            return False
        else:
            do_phase_1_lifting = True
    generate_minimal_triples(fn)
    if show_plots:
        logging.info("Plotting 2d diagram...")
        show(plot_2d_diagram(fn))
        logging.info("Plotting 2d diagram... done")
    covered_intervals = generate_covered_intervals(fn)
    uncovered_intervals = generate_uncovered_intervals(fn)
    if show_plots:
        logging.info("Plotting covered intervals...")
        show(plot_covered_intervals(fn))
        logging.info("Plotting covered intervals... done")
    if not uncovered_intervals:
        print "All intervals are covered (or connected-to-covered).", len(covered_intervals), "components."
        return finite_dimensional_extremality_test(fn, show_plots)
    else:
        print "Uncovered intervals: ", uncovered_intervals
        if do_phase_1_lifting or finite_dimensional_test_first:
            # First try the finite dimensional one.
            if not finite_dimensional_extremality_test(fn, show_plots):
                return False
        # Now do the magic.
        moves = generate_functional_directed_moves(fn)
        print "Moves relevant for these intervals: ", moves
        seed, stab_int, walk_list = find_generic_seed(fn, max_num_it=max_num_it) # may raise MaximumNumberOfIterationsReached
        fn._seed = seed
        fn._stab_int = stab_int
        fn._walk_list = walk_list
        if show_plots:
            logging.info("Plotting moves and reachable orbit...")
            # FIXME: Visualize stability intervals?
            (plot_walk(walk_list,thickness=0.7) + \
             plot_possible_and_impossible_directed_moves(seed, moves, fn) + \
             plot_intervals(uncovered_intervals) + plot_covered_intervals(fn)).show(figsize=50)
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

def random_piecewise_function(xgrid, ygrid):
    xvalues = [0] + [ x/xgrid for x in range(1, xgrid) ] + [1]
    f = randint(1, xgrid - 1)
    yvalues = [0] + [ randint(0, ygrid) / ygrid for i in range(1, f) ] + [1] + [ randint(0, ygrid) / ygrid for i in range(f+1, xgrid) ]+ [0]
    return piecewise_function_from_breakpoints_and_values(xvalues, yvalues)


def the_irrational_function_t1_t2_t3(d1 = 3/5, d3 = 1/10, f = 4/5, p = 15/100, \
                                     del1 = 1/200, del2 = 6*sqrt(2)/200, del3 = 1/500):

    d2 = f - d1 - d3

    c2 = 0
    c3 = -1/(1-f)
    c1 = (1-d2*c2-d3*c3)/d1

    d12 = d2 / 2
    d22 = d2 / 2

    d13 = c1 / (c1 - c3) * d12
    d43 = d13
    d11 = p - d13
    d31 = p - d12
    d51 = d11
    d41 = (d1 - d11 - d31 - d51)/2
    d21 = d41
    d23 = (d3 - d13 - d43)/2
    d33 = d23

    del13 = c1 * del1 / (c1 - c3)
    del11 = del1 - del13

    del23 = c1 * del2 / (c1 - c3)
    del21 = del2 - del23

    del33 = c1 * del3 / (c1 - c3)
    del31 = del3 - del33

    d21new = d21 - del11 - del21 - del31
    d41new = d41 - del11 - del21 - del31
    d23new = d23 - del13 - del23 - del33
    d33new = d33 - del13 - del23 - del33

    t1 = del1
    t2 = del1 + del2
    t3 = del1 + del2 + del3

    a0 = d11+d13
    a1 = a0 + t1
    a2 = a0 + t2
    a3 = a0 + t3

    A = a0+d21+d23
    A0 = A + d12

    slopes = [c1,c3,c1,c3,c1,c3,c1,c3,c1,c3,c2,c1,c2,c3,c1,c3,c1,c3,c1,c3,c1,c3,c1,c3]

    interval_lengths = [d11,d13,del11,del13,del21,del23,del31, del33, d21new,d23new,d12,d31,d22,d33new,d41new,del33,del31,del23,del21,del13,del11,d43,d51,1-f]

    return piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes) 

def is_QQ_linearly_independent(*numbers):
    """
    Test if `numbers` are linearly independent over `QQ`.

    EXAMPLES::
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
    all_QQ, numbers = is_all_QQ(numbers)
    if all_QQ:
        return False
    # otherwise try to coerce to common number field
    numbers = nice_field_values(numbers, RealNumberField)
    if not is_real_number_field_element(numbers[0]):
        raise ValueError, "Q-linear independence test only implemented for algebraic numbers"
    coordinate_matrix = matrix(QQ, [x.parent().0.coordinates_in_terms_of_powers()(x) for x in numbers])
    return rank(coordinate_matrix) == len(numbers)

def compose_directed_moves(A_move, B_move):
    "`A` after `B`."
    return (A_move[0] * B_move[0], A_move[0] * B_move[1] + A_move[1])

def interval_list_intersection(interval_list_1, interval_list_2):
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
            if len(overlapped_int) >= 1:
                overlap.append(overlapped_int)
    return sorted(overlap)

def compose_functional_directed_moves(A, B):
    """
    Compute the directed move that corresponds to the directed move `A` after `B`.
    
    EXAMPLES::
        sage: compose_functional_directed_moves(FunctionalDirectedMove([(5/10,7/10)],(1, 2/10)),FunctionalDirectedMove([(2/10,4/10)],(1,2/10)))
        <FastPiecewise with 1 parts, 
         (3/10, 2/5)\t<FastLinearFunction x + 2/5>\t values: [7/10, 4/5]>
    """
    result_directed_move = compose_directed_moves(A.directed_move, B.directed_move)
    A_domain_preimages = [ B.apply_to_interval(A_domain_interval, inverse=True) \
                           for A_domain_interval in A.intervals() ]
    result_domain_intervals = interval_list_intersection(A_domain_preimages, B.intervals())

    #print result_domain_intervals
    if len(result_domain_intervals) > 0:
        return FunctionalDirectedMove(result_domain_intervals, result_directed_move)
    else:
        return None

def plot_compose_functional_directed_moves(A, B):
    C = compose_functional_directed_moves(A, B)
    p = plot(A, color="green", legend_label="A")
    p += plot(B, color="blue", legend_label="B")
    p += plot(C, color="red", legend_label="C = A after B")
    return p

def merge_functional_directed_moves(A, B, show_plots=False):
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
        show(p)
    return C

def plot_directed_moves(dmoves):
    return sum(plot(dm) for dm in dmoves)

def functional_directed_move_composition_completion(functional_directed_moves, max_num_rounds=None, show_plots=False):
    # Only takes care of the functional case.
    move_dict = { fdm.directed_move : fdm for fdm in functional_directed_moves }
    any_change = True
    num_rounds = 0

    while any_change and (not max_num_rounds or num_rounds < max_num_rounds):
        if show_plots:
            show(plot_directed_moves(list(move_dict.values())))
        any_change = False
        critical_pairs = [ (a, b) for a in move_dict.values() for b in move_dict.values() ]
        for (a, b) in critical_pairs:
            c = compose_functional_directed_moves(a, b)
            if c:
                cdm = c.directed_move
                if cdm in move_dict:
                    merged = merge_functional_directed_moves(move_dict[cdm], c, show_plots=False)
                    if merged.end_points() != move_dict[cdm].end_points():
                        # Cannot compare the functions themselves because of the "hash" magic of FastPiecewise.
                        #print "merge: changed from %s to %s" % (move_dict[cdm], merged)
                        move_dict[cdm] = merged
                        any_change = True
                    else:
                        #print "merge: same"
                        pass
                else:
                    move_dict[cdm] = c
                    any_change = True
        num_rounds += 1

    return list(move_dict.values())

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
        if delta == 1:  # beginning of interval
            assert move not in moves
            moves.add(move)
        elif delta == -1:
            moves.remove(move)
        elif delta == 0:
            pass
        else:
            raise ValueError, "Bad scan item"

def find_decomposition_into_stability_intervals_with_completion(fn, show_plots=False, max_num_it=None):
    fn._stability_orbits = []
    uncovered_intervals = generate_uncovered_intervals(fn)
    intervals = uncovered_intervals
    functional_directed_moves = generate_functional_directed_moves(fn)
    completion = functional_directed_move_composition_completion(functional_directed_moves, max_num_rounds=None)

    decomposition = find_decomposition_into_intervals_with_same_moves(completion)
     
    done_intervals = set()
    
    for (interval, moves) in decomposition:
        if interval not in done_intervals:
            print interval
            orbit = set()
            walk_dict = dict()
            seed = (interval.a + interval.b) / 2
            sign_contradiction = False
            for move in moves:
                moved_interval = move.apply_to_coho_interval(interval)
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
    logging.info("Total: %s stability orbits, lengths: %s" \
                 % (len(fn._stability_orbits), \
                    [ ("%s+" if to_do else "%s") % len(shifted_stability_intervals) \
                      for (shifted_stability_intervals, walk_dict, to_do) in fn._stability_orbits ]))



