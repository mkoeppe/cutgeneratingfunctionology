import logging

logging.basicConfig(format='%(levelname)s: %(asctime)s %(message)s', level=logging.INFO)

import itertools

def fractional(num):
    """
    Reduce a number modulo 1.
    """
    while num > 1:
        num = num - 1
    while num < 0:
        num = num + 1
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
    p1 = plot([(bkpt[i]-x).plot((x, 0, bkpt[i])) for i in range(1,len(bkpt))], \
              **ticks_keywords(function, True))
    p2 = plot([(1+bkpt[i]-x).plot((x, bkpt[i], 1)) for i in range(1,len(bkpt)-1)])
    p3 = plot([(bkpt[i]+x-x).plot((x, 0, 1)) for i in range(len(bkpt))])
    y=var('y')
    for i in range(len(bkpt)):
        p3 += parametric_plot((bkpt[i],y),(y,0,1))
    return p1+p2+p3

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

def plot_2d_diagram(function, show_function = False):
    """
    Return a plot of the 2d complex with shaded faces where delta_pi is 0.        
    To show only a part of it, use 
    `show(plot_2d_diagram(h), xmin=0.25, xmax=0.35, ymin=0.25, ymax=0.35)`
    """
    ### FIXME: For non-subaddititive functions, the points where delta_pi
    ### is negative should be indicated somehow!!
    bkpt = function.end_points()
    vert_face_additive = generate_vert_face_additive(function)
    minimal_triples = generate_minimal_triples(function)
        
    y = var('y')

    plot_minimal_triples = plot_2d_complex(function)
    for trip, vert in itertools.izip(minimal_triples, vert_face_additive):
        if face_0D(trip):
            plot_minimal_triples += point((trip[0][0], \
                trip[1][0]), color = "magenta", size = 30)
        elif face_horizontal(trip):
            plot_minimal_triples += parametric_plot((y,trip[1][0]),\
                (y,trip[0][0], trip[0][1]), rgbcolor=(0, 1, 0))
        elif face_vertical(trip):
            plot_minimal_triples += parametric_plot((trip[0][0],y),\
                (y,trip[1][0], trip[1][1]), rgbcolor=(0, 1, 0))
        elif face_diagonal(trip):
            plot_minimal_triples += parametric_plot((y,trip[2][0]-y),\
                (y,trip[0][0],trip[0][1]), rgbcolor=(0, 1, 0))
        elif face_2D(trip):
            ## Sorting is necessary for this example:
            ## plot_2d_diagram(lift(piecewise_function_from_robert_txt_file("/Users/mkoeppe/w/papers/basu-hildebrand-koeppe-papers/reu-2013/Sage_Function/dey-richard-not-extreme.txt"))
            plot_minimal_triples += polygon(convex_vert_list(vert), rgbcolor=(1, 0, 0)) 
    if show_function:
        x = var('x')
        plot_minimal_triples += parametric_plot((lambda x: x, lambda x: 0.3 * function(x) + 1), \
                                                (x, 0, 1), color='black')
        plot_minimal_triples += parametric_plot((lambda x: - 0.3 * function(x), lambda x: x), \
                                                (x, 0, 1), color='black')
    return plot_minimal_triples

# Assume component is sorted.
def merge_within_comp(component):   
    for i in range(len(component)-1):
        if component[i][1] > component[i+1][0]:
            component[i+1] = [component[i][0],max(component[i][1],component[i+1][1])]
            component[i] = []
    component_new = []
    for int in component:
        if len(int) == 2 and max(int) <= 1:
            component_new.append(int)
    return component_new


# Assume comp1 and comp2 are sorted.    
def merge_two_comp(comp1,comp2):
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
    temp = merge_within_comp(temp)
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
    assert interval[0] < interval[1]
    while interval[0] >= 1:
        interval[0] = interval[0] - 1
        interval[1] = interval[1] - 1
    while interval[1] <= 0:
        interval[0] = interval[0] + 1
        interval[1] = interval[1] + 1
    assert not(interval[0] < 1 and interval[1] > 1) 
    return interval

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
    """
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

def directed_moves_from_moves(moves):
    directed_moves = []
    for move in moves:
        move_sign = move[0]
        if move_sign == 1:
            directed_moves.append(move)
            directed_moves.append([1,-move[1]])
        elif move_sign == -1:
            directed_moves.append(move)
        else:
            raise ValueError, "Move not valid: %s" % list(move)
    return directed_moves

def is_directed_move_possible(x, move, fn=None, intervals=None):
    if fn:
        move_sign = move[0]
        if move_sign == 1:
            if move[1] >= 0:
                return fn(x) + fn(move[1]) == fn(fractional(x + move[1]))
            else:
                return fn(fractional(x + move[1])) + fn(-move[1]) == fn(x)
        else:
            return fn(x) + fn(fractional(move[1]-x)) == fn(move[1])
    elif intervals:
        next_x = apply_directed_move(x, move)
        for interval in intervals:
            if (interval[0] <= next_x <= interval[1]):
                return True
        return False
    else:
        raise ValueError, "Need either fn or intervals to check if move is possible"

def find_possible_directed_moves(x,moves,fn=None, intervals=None):
    """Find the directed moves applicable at x.
    """
    directed_moves = directed_moves_from_moves(moves)
    possible_directed_moves = []
    for directed_move in directed_moves:
        if is_directed_move_possible(x, directed_move, fn, intervals):
            possible_directed_moves.append(directed_move)
    return possible_directed_moves

def find_impossible_directed_moves(x,moves,fn=None, intervals=None):
    """
    Find the directed moves NOT applicable at x.
    """
    directed_moves = directed_moves_from_moves(moves)
    impossible_directed_moves = []
    for directed_move in directed_moves:
        if not is_directed_move_possible(x, directed_move, fn, intervals):
            impossible_directed_moves.append(directed_move)
    return impossible_directed_moves

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
        possible_directed_moves = find_possible_directed_moves(x,moves,fn)
        move = possible_directed_moves[random.randint(0,len(possible_directed_moves)-1)]
        move_sign = move[0]
        directed_move = move
        
        if move_sign == 1:
            next_x = fractional(x + move[1])
        elif move_sign == -1:
            next_x = fractional(move[1]-x)
            
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

    def __eq__(self, other):
        if not isinstance(other, FastLinearFunction):
            return False
        return self._slope == other._slope and self._intercept == other._intercept

    def __ne__(self, other):
        return not (self == other)

    __rmul__ = __mul__

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
        if last_f != None:
            merged_list_of_pairs.append(((merged_interval_a, merged_interval_b), last_f))
            
        PiecewisePolynomial.__init__(self, merged_list_of_pairs, var)
        self.update_cache()

    # The following makes this class hashable and thus enables caching
    # of the above functions; but we must promise not to modify the
    # contents of the instance.
    def __hash__(self):
        return id(self)

    def update_cache(self):
        intervals = self._intervals
        functions = self._functions
        end_points = [ intervals[0][0] ] + [b for a,b in intervals]
        self._end_points = end_points
        values_at_end_points = [ functions[0](end_points[0]) ]
        for i in range(len(functions)):
            value = functions[i](intervals[i][1])
            values_at_end_points.append(value)
        self._values_at_end_points = values_at_end_points

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
        if endpts[i-1] < x0 < endpts[i]:
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
    ## (see doctests below).
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
        return g

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
    
        
        

def piecewise_function_from_breakpoints_and_values(bkpt, values, field=default_field):
    """
    bkpt and values are two parallel lists; assuming bpkt is sorted (increasing).
    Return a function.
    """
    if len(bkpt)!=len(values):
        raise DataError, "Need to have the same number of breakpoints and values."
    slopes = [ (values[i+1]-values[i])/(bkpt[i+1]-bkpt[i]) for i in range(len(bkpt)-1) ]
    return FastPiecewise([ [(bkpt[i],bkpt[i+1]), 
                        # lambda x,i=i: values[i] + slopes[i] * (x - bkpt[i])
                        fast_linear_function(slopes[i], values[i] - slopes[i]*bkpt[i])] for i in range(len(bkpt)-1) ])

def piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=default_field):
    """
    bkpt and slopes are two parallel lists (except that bkpt is one element longer).
    The function always has value 0 on bkpt[0].
    Return a function.
    """
    if len(bkpt)!=len(slopes)+1:
        raise DataError, "Need to have one breakpoint more than slopes."
    function_values = [0]
    for i in range(1,len(bkpt)-1):
        function_values.append(function_values[i-1] + slopes[i - 1] * (bkpt[i] - bkpt[i-1]))
    pieces = [[(bkpt[i],bkpt[i+1]),
               #lambda x,i=i: function_values[i]+slopes[i]*(x - bkpt[i])
               fast_linear_function(slopes[i], function_values[i] - slopes[i]*bkpt[i])] for i in range(len(bkpt)-1)] 
    return FastPiecewise(pieces)        

def piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes, field=default_field):
    """
    interval_lengths and slopes are two parallel lists.
    The function always has value 0 on 0.
    Return a function.
    """
    if len(interval_lengths)!=len(slopes):
        raise DataError, "Number of given interval_lengths and slopes needs to be equal."
    bkpt = []
    bkpt.append(0)
    for i in range(len(interval_lengths)):
        bkpt.append(bkpt[i]+interval_lengths[i])
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field)

def approx_discts_function(perturbation_list, stability_interval, field=default_field, \
                           symmetric=False):
    """
    Construct a function that has peaks of +/- 1 around the points of the orbit.
    perturbation_list actually is a dictionary.
    """
    perturb_points = sorted(perturbation_list.keys())
    fn_values = [0]
    fn_bkpt = [0]
    # This width is chosen so that the peaks are disjoint, and
    # so a nice continuous piecewise linear function is constructed.
    width = min(abs(stability_interval.a),stability_interval.b)
    assert width > 0, "Width of stability interval should be positive"
    assert stability_interval.a < 0 < stability_interval.b, \
        "Stability interval should contain 0 in it s interior"
    for pt in perturb_points:
        if symmetric:
            left = pt - width
            right = pt + width
        else:
            if perturbation_list[pt][0] == 1:
                left = pt + stability_interval.a
                right = pt + stability_interval.b
            else:
                left = pt - stability_interval.b
                right = pt - stability_interval.a
        assert (left >= fn_bkpt[len(fn_bkpt)-1])
        if (left > fn_bkpt[len(fn_bkpt)-1]):
            fn_bkpt.append(left)
            fn_values.append(0)
        fn_bkpt.append(pt)
        fn_values.append(perturbation_list[pt][0]) # the "walk_sign" (character) at the point
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

def canonicalize_number(number):
    """Make sure that if `number` is a rational number, then it is
    represented as an element of `Rational` (rather than an element of
    `QQbar`, for example).  This will make sure that we do not have
    two mathematically equal numbers with different `hash` values."""
    try:
        return QQ(number)
    except ValueError:
        return number
    except TypeError:
        return number

def apply_directed_move(x, directed_move):
    move_sign = directed_move[0]
    if move_sign == 1:
        next_x = fractional(x + directed_move[1])
    elif move_sign == -1:
        next_x = fractional(directed_move[1]-x)
    next_x = canonicalize_number(next_x)
    return next_x

class MaximumNumberOfIterationsReached(Exception):
    pass

class SignContradiction(Exception):
    pass

def deterministic_walk(seed, moves, fn=None, max_num_it = 1000, intervals=None, error_if_sign_contradiction=False):
    """
    Compute the orbit of a given seed. (Done by a breadth-first search.)
    To avoid infinite computations in the case of a dense orbit,
    there is a maximum number of iterations (by default, it is 1000).
    Which moves are allowed is decided by testing the Delta of `fn',
    or, if that is not provided, by testing whether we stay in `intervals'.    

    Returns a dictionary:
    - keys are elements of the orbit
    - values are lists of the form [walk_sign, predecessor, directed_move_from_predecessor].

    """
    ## FIXME: If `fn` is provided, store the dictionary in `fn.walk_dict`
    ## and the to_do list in `fn.walk_to_do`, to allow us to resume
    ## an interrupted BFS.  OR: Turn it into an iterator/generator (yield).

    logging.info("Breadth-first search to discover the reachable orbit...")
    seed = canonicalize_number(seed)
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
        possible_directed_moves = find_possible_directed_moves(x, moves, fn, intervals)
        for directed_move in possible_directed_moves:
            walk_sign = xlist[x][0]
            move_sign = directed_move[0]
            move_for_next_x = directed_move
            next_x = apply_directed_move(x, directed_move)
            
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
        
    if num_it == max_num_it:
        raise MaximumNumberOfIterationsReached, "Reached %d iterations, to do list has still %d items" % (num_it, len(to_do))
 
    logging.info("Breadth-first search to discover the reachable orbit... done")
    return xlist

def plot_walk(walk_dict, **options):
    #return point([ (x,0) for x in walk_dict.keys()])
    g = Graphics()
    for x in walk_dict.keys():
        g += line([(x,0), (x,1)], color="black", zorder = -4, **options)
    return g

def plot_intervals(intervals):
    g = Graphics()
    for interval in intervals:
        g += polygon([(interval[0], 0), (interval[1], 0), \
                      (interval[1], 1.0), (interval[0], 1.0)], 
                     color="yellow", zorder = -8)
    return g

def plot_moves(seed, moves, colors=None):
    if colors == None:
        colors = rainbow(len(moves))
    g = Graphics()
    g += line([(seed,0), (seed,1)], color="magenta")
    y = 0
    covered_interval = [0,1]
    for move, color in itertools.izip(moves, colors):
        next_x = apply_directed_move(seed, move)
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
            bezier_y = y + min(0.03, 0.3 * abs(move[1]/2 - seed))
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

def plot_possible_and_impossible_directed_moves(seed, moves, fn):
    directed_moves = directed_moves_from_moves(moves)
    colors = [ "blue" if is_directed_move_possible(seed, directed_move, fn) \
               else "red" for directed_move in directed_moves ]
    return plot_moves(seed, directed_moves, colors)

import collections
_closed_or_open_or_halfopen_interval = collections.namedtuple('Interval', ['a', 'b', 'left_closed', 'right_closed'])

class closed_or_open_or_halfopen_interval (_closed_or_open_or_halfopen_interval):
    def __repr__(self):
        return "<Int" \
            + ("[" if self.left_closed else "(") \
            + repr(self.a) + ", " + repr(self.b) \
            + ("]" if self.right_closed else ")") \
            + ">"

def one_step_stability_interval(x, intervals, moves):
    """Returns the stability interval, i.e., an open, half-open, or closed
    interval I), such that for all moves, `move`(`x`) lies in
    `intervals` if and only if `move`(`x` + t) lies in `intervals` for
    t in I."""
    a = -10
    b = 10
    left_closed = True
    right_closed = True
    for move in directed_moves_from_moves(moves):
        next_x = apply_directed_move(x, move)
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
                    if (interval[1] - next_x) < -1 * a:
                        a = next_x - interval[1]
                        left_closed = True      
                else:
                    raise ValueError, "Move not valid: %s" % list(move)
        # Move was not possible.
        if move[0] == 1:
            for interval in intervals:
                temp = interval[0] - next_x
                if 0 < temp <= b:
                    b = temp
                    right_closed = False
                temp2 = interval[1] - next_x
                if a <= temp2 < 0:
                    a = temp
                    left_closed = False
        elif move[0] == -1:
            for interval in intervals:
                temp = interval[0] - next_x
                if 0 < temp <= -1 * a:
                    a = -1 * temp
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
        logging.info("%s of length %s " % (int, RR(length)))
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
    moves = generate_moves(fn)
    return iterative_stability_refinement(intervals, moves)

def find_stability_interval_with_deterministic_walk_list(seed, intervals, moves, fn, max_num_it = 1000, error_if_sign_contradiction=False):
    """
    Returns the stability interval (an open, half-open, or closed interval)
    and the deterministic_walk_list.
    """
    ## FIXME: Refactor using above.
    a = -10
    b = 10
    left_closed = True
    right_closed = True
    deterministic_walk_list = deterministic_walk(seed,moves,fn, max_num_it, \
                                                 error_if_sign_contradiction=error_if_sign_contradiction)
    if deterministic_walk_list[seed][0] == 0:
        return (closed_or_open_or_halfopen_interval(0, 0, True, True), deterministic_walk_list)
    for pt in deterministic_walk_list.keys():
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
                    if (interval[1] - pt) < -1 * a:
                        a = pt - interval[1]
                        left_closed = True      
        impossible_directed_moves = find_impossible_directed_moves(pt, moves, fn)
        ### We now take the set difference of
        ### the __directed__ moves with the possible directed moves.
        ### This was not done in the old code: --Matthias
        # impossible_directed_moves = []
        # for move in moves:
        #     if move not in possible_directed_moves:
        #         impossible_directed_moves.append(move)
        
        for move in impossible_directed_moves:
            if move[0] == 1:
                impossible_next_x = fractional(pt + move[1])
                for interval in intervals:
                    temp = interval[0] - impossible_next_x
                    if 0 < temp <= b:
                        b = temp
                        right_closed = False
                    temp2 = interval[1] - impossible_next_x
                    if a <= temp2 < 0:
                        a = temp
                        left_closed = False
            elif move[0] == -1:
                impossible_next_x = fractional(move[1] - pt)
                for interval in intervals:
                    temp = interval[0] - impossible_next_x
                    if 0 < temp <= -1 * a:
                        a = -1 * temp
                        left_closed = False
                    temp2 = impossible_next_x - interval[1]
                    if 0 < temp2 <= b:
                        b = temp2
                        right_closed = False
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
                #print "half_distance wins:", half_distance 
                b = half_distance
                right_closed = False
        elif sign_1 == -1 and sign_2 == 1:
            if -a > half_distance:
                #print "half_distance wins:", half_distance 
                a = -half_distance
                left_closed = False
        else:
            assert(-a + b <= 2 * half_distance)
    return (closed_or_open_or_halfopen_interval(a, b, left_closed, right_closed), deterministic_walk_list)



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

@cached_function
def generate_moves(fn, intervals=None):
    """
    Compute the moves (translations and reflections) relevant for the given intervals
    (default: all uncovered intervals).
    """
    if intervals==None:
        # Default is to generate moves for ALL uncovered intervals
        intervals = generate_uncovered_intervals(fn)
    moves = set()
    for trip in generate_minimal_triples(fn):
        intersects = [find_interior_intersection([trip[i]], intervals) for i in range(3)]
        if face_horizontal(trip):
            if intersects[0] or intersects[2]:
                #assert intersects[0] and intersects[2]
                if not (intersects[0] and intersects[2]):
                    logging.warn("generate_moves: Tricky case hit, think harder!")
                moves.add((1, trip[1][0]))               # Translation
        elif face_diagonal(trip):
            if intersects[0] or intersects[1]:
                #assert intersects[0] and intersects[1]
                if not (intersects[0] and intersects[1]):
                    logging.warn("generate_moves: Tricky case hit, think harder!")
                moves.add((-1, fractional(trip[2][0])))  # Reflection
    moves.remove((1,0))                                  # Remove the trivial translation
    return list(moves)

def find_generic_seed(fn, max_num_it = 1000):
    intervals = generate_uncovered_intervals(fn)
    if not intervals:
        raise ValueError, "Need an uncovered interval"
    moves = generate_moves(fn)
    seed = intervals[0][1]
    while True:
        seed = 2/3 * intervals[0][0] + 1/3 * seed
        logging.info("Trying seed %s" % seed)
        try:
            stab_int, walk_list = find_stability_interval_with_deterministic_walk_list \
                                  (seed, intervals, moves, fn, max_num_it = max_num_it, \
                                   error_if_sign_contradiction = True)
            if stab_int[1] == stab_int[0]:
                logging.info("Stability interval is not proper, continuing search.")
                continue
            logging.info("Seed %s has a proper stability interval %s, reachable orbit has %s elements" % (seed, stab_int, len(walk_list)))
            return (seed, stab_int, walk_list)
        except SignContradiction:
            continue

class UnimplementedError (Exception):
    pass

def generate_compatible_piecewise_function(components, component_slopes):
    intervals_and_slopes = []
    for component, slope in itertools.izip(components, component_slopes):
        intervals_and_slopes.extend([ (interval, slope) for interval in component ])
    intervals_and_slopes.sort()
    bkpt = [ int[0] for int, slope in intervals_and_slopes ] + [1]
    slopes = [ slope for int, slope in intervals_and_slopes ]
    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes)

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
    symbolic = generate_compatible_piecewise_function(components, slope_vars)
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
    maximum absolute function value is 1.
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


def extremality_test(fn, show_plots = False, max_num_it = 1000):
    if not minimality_test(fn):
        print "Not minimal, thus not extreme."
        return False
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
        moves = generate_moves(fn)
        print "Moves relevant for these intervals: ", moves
        seed, stab_int, walk_list = find_generic_seed(fn, max_num_it=max_num_it) # may raise MaximumNumberOfIterationsReached
        if show_plots:
            logging.info("Plotting moves and reachable orbit...")
            # FIXME: Visualize stability intervals?
            (plot_walk(walk_list,thickness=0.7) + \
             plot_possible_and_impossible_directed_moves(seed, moves, fn) + \
             plot_intervals(uncovered_intervals) + plot_covered_intervals(fn)).show(figsize=50)
            logging.info("Plotting moves and reachable orbit... done")
        perturb = fn._perturbation = approx_discts_function(walk_list, stab_int)
        check_perturbation(fn, perturb, show_plots=show_plots)
        return False

def lift(fn, show_plots = False):
    if extremality_test(fn, show_plots=show_plots):
        return fn
    else:
        perturbed = fn._lifted = fn + fn._epsilon_interval[1] * fn._perturbation
        return perturbed

def lift_until_extreme(fn, show_plots = False):
    next, fn = fn, None
    while next != fn:
        fn = next
        next = lift(fn, show_plots=show_plots)
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

