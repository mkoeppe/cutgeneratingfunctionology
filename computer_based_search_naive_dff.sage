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
   
    