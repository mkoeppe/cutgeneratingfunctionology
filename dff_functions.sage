def phi_bj_1(c):
    """
    Summary:
        - Name: FBJ;
        - Dim= 1; Slopes = 2; Continuous;
        - Proven maximal [1] p.26, theorem 2.1。
    
    Parameters:
        c (real) \in [1,positive infinity).

    Examples:
        [1] p.25, Fig 2.4 ::
            sage: h=phi_bj_1(3/2))
            sage: superadditive_test(h)
            True

    References:
   
    - [1]: Alves, Clautiaux, Carvalho, Rietz, Dual-feasible functions for integer programming and combinatorial optimization, 2016.
    """

    n=floor(c)
    beta=c-n
    a=beta/c
    b=1/c-beta/c
    interval=[a,b]*n+[a]
    slope=[0,1/n/b]*n+[0]
    return piecewise_function_from_interval_lengths_and_slopes(interval, slope)

def phi_simple(c):
    n=floor(c)
    l=1/c
    pieces=[]
    for k in range(n):
        pieces = pieces + \
                 [[singleton_interval(l*k), FastLinearFunction(0,k/n)], \
                  [open_interval(l*k,l*k+l), FastLinearFunction(0,k/n)]]
    pieces.append([closed_interval(n/c,1), FastLinearFunction(0,1)])
    return FastPiecewise(pieces)

def phi_ccm_1(c):
    n=floor(c)
    l=1/c
    m=floor(1/2/l)
    pieces=[]
    for k in range(m):
        pieces = pieces + \
                 [[singleton_interval(l*k), FastLinearFunction(0,k/n)], \
                  [open_interval(l*k,l*k+l), FastLinearFunction(0,k/n)]]
    pieces.append([right_open_interval(m/c,1/2), FastLinearFunction(0,m/n)])
    pieces.append([singleton_interval(1/2), FastLinearFunction(0,1/2)])
    pieces.append([left_open_interval(1/2,1-m/c), FastLinearFunction(0,1-m/n)])
    for k in range(m):
        left_x=1-(m-k)*l
        pieces = pieces + \
                 [[left_open_interval(left_x,left_x+l), FastLinearFunction(0,1-(m-k-1)/n)]]
    return FastPiecewise(pieces)

def phi_fs_1(k):
    pieces=[]
    pieces.append([singleton_interval(0),FastLinearFunction(0,0)])
    for i in range(k+1):
        pieces=pieces + \
               [[open_interval(i/(k+1),(i+1)/(k+1)),FastLinearFunction(0,i/k)],\
                 [singleton_interval((i+1)/(k+1)),FastLinearFunction(0,(i+1)/(k+1))]]
    return FastPiecewise(pieces)

def phi_vb_2(k):
    l=1/k
    if mod(k,2)==1:
        pieces=[]
        pieces.append([singleton_interval(0),FastLinearFunction(0,0)])
        m=(k-1)/2
        for i in range(m):
            pieces.append([left_open_interval(i*l,i*l+l),FastLinearFunction(0,i/(k-1))])
        pieces.append([open_interval(1/2-l/2,1/2+l/2),FastLinearFunction(0,1/2)])
        for i in range(m):
            pieces.append([right_open_interval((m+1+i)*l,(m+1+i)*l+l),FastLinearFunction(0,(m+1+i)/(k-1))])
        pieces.append([singleton_interval(1),FastLinearFunction(0,1)])
        return FastPiecewise(pieces)       
    else:
        pieces=[]
        pieces.append([singleton_interval(0),FastLinearFunction(0,0)])
        m=k/2
        for i in range(m-1):
            pieces.append([left_open_interval(i*l,i*l+l),FastLinearFunction(0,i/(k-1))])
        pieces.append([open_interval(1/2-l,1/2),FastLinearFunction(0,(m-1)/(k-1))])
        pieces.append([singleton_interval(1/2),FastLinearFunction(0,1/2)])
        pieces.append([open_interval(1/2,1/2+l),FastLinearFunction(0,m/(k-1))])
        for i in range(m-1):
            pieces.append([right_open_interval((m+1+i)*l,(m+1+i)*l+l),FastLinearFunction(0,(m+1+i)/(k-1))])
        pieces.append([singleton_interval(1),FastLinearFunction(0,1)])
        return FastPiecewise(pieces)

def phi_ll_1(c,k):
    n=floor(c)
    beta=c-n
    l1=beta/c
    l2=(1-beta)/c/(k-1)
    pieces=[]
    for i in range(n):
        pieces.append([closed_interval(i/c,i/c+l1), FastLinearFunction(0,i/n)])    
        for j in range(k-1):
            if j==k-2:
               pieces.append([open_interval(i/c+l1+j*l2,i/c+l1+(j+1)*l2), FastLinearFunction(0,(i+(j+1)/k)/n)])
            else:
               pieces.append([left_open_interval(i/c+l1+j*l2,i/c+l1+(j+1)*l2), FastLinearFunction(0,(i+(j+1)/k)/n)])
    pieces.append([closed_interval(n/c,1), FastLinearFunction(0,1)])
    return FastPiecewise(pieces)

def phi_ll_2(c,k):
    n=floor(c)
    beta=c-n
    l1=beta/c
    l2=(1-beta)/c/(k-1)
    pieces=[]
    m=floor(n/2)
    for i in range(m):
        pieces.append([closed_interval(i/c,i/c+l1), FastLinearFunction(0,i/n)])    
        for j in range(k-1):
            if j==k-2:
               pieces.append([open_interval(i/c+l1+j*l2,i/c+l1+(j+1)*l2), FastLinearFunction(0,(i+(j+1)/k)/n)])
            else:
               pieces.append([left_open_interval(i/c+l1+j*l2,i/c+l1+(j+1)*l2), FastLinearFunction(0,(i+(j+1)/k)/n)])
    if (1/2-m/c)<=l1:
        pieces.append([right_open_interval(m/c,1/2), FastLinearFunction(0,m/n)])
        pieces.append([singleton_interval(1/2),FastLinearFunction(0,1/2)])
        pieces.append([left_open_interval(1/2,1-m/c), FastLinearFunction(0,1-m/n)])
    else:
        d=floor((1/2-m/c-l1)/l2)
        pieces.append([closed_interval(m/c,m/c+l1), FastLinearFunction(0,m/n)])
        if d<(1/2-m/c-l1)/l2:
            for t in range(d):
                pieces.append([left_open_interval(m/c+l1+t*l2,m/c+l1+(t+1)*l2), FastLinearFunction(0,(m+(t+1)/k)/n)])
            pieces.append([open_interval(m/c+l1+d*l2,1/2), FastLinearFunction(0,(m+(d+1)/k)/n)])
            pieces.append([singleton_interval(1/2),FastLinearFunction(0,1/2)])
            pieces.append([open_interval(1/2,1-(m/c+l1+d*l2)), FastLinearFunction(0,1-(m+(d+1)/k)/n)])
            for t in range(d):
                pieces.append([right_open_interval(1-m/c-l1-(d-t)*l2,1-m/c-l1-(d-t)*l2+l2), FastLinearFunction(0,1-(m+(d-t)/k)/n)])
        else:
            for t in range(d-1):
                pieces.append([left_open_interval(m/c+l1+t*l2,m/c+l1+(t+1)*l2), FastLinearFunction(0,(m+(t+1)/k)/n)])
            pieces.append([open_interval(1/2-l2,1/2), FastLinearFunction(0,(m+d/k)/n)])
            pieces.append([singleton_interval(1/2),FastLinearFunction(0,1/2)])
            pieces.append([open_interval(1/2,1/2+l2), FastLinearFunction(0,1-(m+d/k)/n)])
            for t in range(d-1):
                pieces.append([right_open_interval(1-m/c-l1-(d-t)*l2,1-m/c-l1-(d-t)*l2+l2), FastLinearFunction(0,1-(m+(d-t)/k)/n)])
        pieces.append([closed_interval(1-m/c-l1,1-m/c), FastLinearFunction(0,1-m/n)])
    for i in range(m):        
        for j in range(k-1):
            if j==0:
               pieces.append([open_interval(1-(m-i)/c+j*l2,1-(m-i)/c+(j+1)*l2), FastLinearFunction(0,1-(m-i-(j+1)/k)/n)])
            else:
               pieces.append([right_open_interval(1-(m-i)/c+j*l2,1-(m-i)/c+(j+1)*l2), FastLinearFunction(0,1-(m-i-(j+1)/k)/n)])            
        pieces.append([closed_interval(1-(m-i-1)/c-l1,1-(m-i-1)/c), FastLinearFunction(0,1-(m-1-i)/n)])    
    return FastPiecewise(pieces)

def phi_dg_1(c,k):
    n=floor(c)
    beta=c-n
    l1=beta/c
    l2=(1-beta)/c/(k-1)
    pieces=[]
    for i in range(n):
        pieces.append([closed_interval(i/c,i/c+l1), FastLinearFunction(0,i/n)])    
        for j in range(k-1):
            if j==k-2:
               pieces.append([open_interval(i/c+l1+j*l2,i/c+l1+(j+1)*l2), FastLinearFunction(0,(i+(j+1)/k)/n)])
            else:
               pieces.append([open_interval(i/c+l1+j*l2,i/c+l1+(j+1)*l2), FastLinearFunction(0,(i+(j+1)/k)/n)])
               pieces.append([singleton_interval(i/c+l1+(j+1)*l2), FastLinearFunction(0,(i+(j+1)/(k-1))/n)])
    pieces.append([closed_interval(n/c,1), FastLinearFunction(0,1)])
    return FastPiecewise(pieces)



