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
    """A generator...
    "...'general' refers to the fact that it outputs 6-tuples (x,xeps,y,yeps,z,zeps).
    FIXME: it currently does not take care of any discontinuities at all.
    """
    bkpt = fn.end_points()
    return ( (x, y, x+y, 0, 0, 0) for x in bkpt for y in bkpt if x <= y and comparison(delta_pi(fn,x,y), 0) ) # generator comprehension

def generate_type_2_vertices_continuous(fn, comparison):
    bkpt = fn.end_points()
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]
    return ( (x, z-x, z, 0, 0, 0) for x in bkpt for z in bkpt2 if x < z < 1+x and comparison(delta_pi(fn, x, z-x), 0) ) # generator comprehension


@cached_function
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
