from six.moves import range
def search_example_continuous_dff(q):
    """
    Generator for the extreme continuous piecewise linear cDFFs on the 1/q grid.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.dff import *
        sage: logging.disable(logging.WARN)   # Suppress output in automatic tests.
        sage: max_q = 13
        sage: dffs_by_q = [ list(search_example_continuous_dff(q)) for q in range(max_q+1) ]   # optional - pynormaliz
        sage: [ len(dffs) for dffs in dffs_by_q ]    # optional - pynormaliz
        [0, 1, 1, 1, 1, 2, 1, 3, 3, 3, 3, 7, 6, 8]
        sage: [ len([ phi for phi in dffs if not is_bj(phi) ]) for dffs in dffs_by_q ]    # optional - pynormaliz
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 4, 2]

    """
    lp = MixedIntegerLinearProgram(base_ring=QQ)
    x = lp.new_variable()
    for i in range(q+1):
        lp.add_constraint(x[i]>=0)
    lp.add_constraint(x[0]==0)
    lp.add_constraint(x[q]==1)
    for i in range(1,q):
        for j in range(i,q):
            if (i+j)<q:
                lp.add_constraint(x[i]+x[j]-x[i+j]<=0)
            if (i+j)==q:
                lp.add_constraint(x[i]+x[j]-1==0)    
    p = lp.polyhedron(backend='normaliz')
    v=p.vertices()
    for k in range(len(v)):
        h=h_from_vertex_values(v[k])
        generate_maximal_additive_faces_dff(h)  # set attribute
        uncovered_intervals=generate_uncovered_intervals(h)
        if not uncovered_intervals:
            yield h

def search_example_general_dff(q): 
    """
    Generator for the extreme, possibly discontinuous piecewise linear cDFFs on the 1/q grid.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.dff import *
        sage: logging.disable(logging.INFO)   # Suppress output in automatic tests.
        sage: max_q = 7
        sage: dffs_by_q = [ list(search_example_general_dff(q)) for q in range(max_q+1) ]   # optional - pynormaliz
        sage: [ len(dffs) for dffs in dffs_by_q ]   # optional - pynormaliz
        [1, 1, 2, 3, 7, 11, 14, 44]
    """

    lp = MixedIntegerLinearProgram(base_ring=QQ)
    fn = lp.new_variable()
    for i in range(3*q-1):
        lp.add_constraint(fn[i]>=0)
    lp.add_constraint(fn[0] == 0)
    lp.add_constraint(fn[3*q-2] == 1)
    for i in range(1,q):
        for j in range(i,q):
            k=i+j
            x=3*i-1
            y=3*j-1
            z=3*k-1
            if k<q:
                lp.add_constraint(fn[x]+fn[y]-fn[z]<=0)
            if k==q:
                lp.add_constraint(fn[x]+fn[y]-1==0)
            for (xeps, yeps, zeps) in nonzero_eps:
                if k<q:
                    lp.add_constraint(fn[x+xeps]+fn[y+yeps]-fn[z+zeps]<=0)              
                if k==q:
                    if zeps==-1:
                        lp.add_constraint(fn[x+xeps]+fn[y+yeps]-1<=0)
                    if zeps==0:
                        lp.add_constraint(fn[x+xeps]+fn[y+yeps]-1==0)
    p = lp.polyhedron(backend='normaliz')
    v=p.vertices()
    for k in range(len(v)):
        h=h_from_vertex_values_general(v[k])
        generate_maximal_additive_faces_dff(h)  # set attribute
        uncovered_intervals=generate_uncovered_intervals(h)
        if not uncovered_intervals:
            yield h


## Another implmementation, using PPL.

def initial_cs_dff(q):
    cs=Constraint_System()
    fn=[Variable(i) for i in range(q+1)]
    cs.insert(fn[0] == 0)
    cs.insert(fn[q] == 1)
    cs.insert(fn[1] >= 0)
    for x in range(1,q):
        for y in range(x,q):
            if (x+y)<q:
                cs.insert(fn[x]+fn[y]-fn[x+y]<=0)
            if (x+y)==q:
                cs.insert(fn[x]+fn[y]-1==0)
    return cs

def search_kslope_example_continuous_dff(k_slopes,q):
    cs=initial_cs_dff(q)
    polytope=C_Polyhedron(cs)
    l=list(generate_vertex_values(k_slopes,q,polytope,set([])))
    for i in range(len(l)):
        h=h_from_vertex_values(l[i])
        uncovered_intervals=generate_uncovered_intervals(h)
        if not uncovered_intervals:
            yield h

def initial_cs_general_dff(q):
    cs=Constraint_System()
    fn=[Variable(i) for i in range(3*q-1)]
    cs.insert(fn[0] == 0)
    cs.insert(fn[3*q-2] == 1)
    cs.insert(fn[1] >= 0)
    for i in range(1,q):
        for j in range(i,q):
            k=i+j
            x=3*i-1
            y=3*j-1
            z=3*k-1
            if k<q:
                cs.insert(fn[x]+fn[y]-fn[z]<=0)
            if k==q:
                cs.insert(fn[x]+fn[y]-1==0)
            for (xeps, yeps, zeps) in nonzero_eps:
                if k<q:
                    cs.insert(fn[x+xeps]+fn[y+yeps]-fn[z+zeps]<=0)              
                if k==q:
                    if zeps==-1:
                        cs.insert(fn[x+xeps]+fn[y+yeps]-1<=0)
                    if zeps==0:
                        cs.insert(fn[x+xeps]+fn[y+yeps]-1==0)
    return cs

def search_kslope_example_general_dff(k_slopes,q):
    cs=initial_cs_general_dff(q)
    polytope=C_Polyhedron(cs)
    l=list(generate_vertex_values(k_slopes,q,polytope,set([])))
    for i in range(len(l)):
        h=h_from_vertex_values_general(l[i])
        uncovered_intervals=generate_uncovered_intervals(h)
        if not uncovered_intervals:
            yield h   

##

def h_from_vertex_values_general(l):
    q=(len(l)+1)/3
    m=max(l)
    value=[QQ(y) / m for y in l]
    pieces=[]
    pieces.append([singleton_interval(0),FastLinearFunction(0,0)])         
    for i in range(q-1):
        slope=(value[3*i+1]-value[3*i])*q
        int=value[3*i]-slope*i/q
        pieces.append([open_interval(i/q,(i+1)/q),FastLinearFunction(slope,int)])
        pieces.append([singleton_interval((i+1)/q),FastLinearFunction(0,value[3*i+2])]) 
    slope=(value[3*q-2]-value[3*q-3])*q
    int=1-slope
    pieces.append([left_open_interval(1-1/q,1),FastLinearFunction(slope,int)])  
    return FastPiecewise(pieces)
