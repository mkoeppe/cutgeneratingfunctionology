# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *


def phi_bj_1_quasi(c=3/2):
    n=floor(c)
    beta=c-n
    a=beta/c
    t=1/c
    list_of_pairs=[[(0,a), FastLinearFunction(0,0)],[(a,t),linear_function_through_points([a, 0], [t, 1/n])]]
    return PiecewiseQuasiPeriodic(list_of_pairs)

def restricted_maximality_test_quasi(phi):
    if not maximality_test_quasi(phi):
        return False
    if phi(1)==1:
        return True
    else:
        return False

def strongly_maximality_test_quasi(phi):
    if not restricted_maximality_test_quasi(phi):
        return False
    limit_values=phi.limits_at_end_points()
    if limit_values[1][-1]>0:
        return False
    else:
        return True

    

def maximality_test_quasi(phi):
    bkpt=phi.end_points()
    #check whether it is linear function
    if len(bkpt)==2:
        if phi(1)<=1 and phi(1)>=0:
            return True
        else:
            return False

    if phi(0)!=0 or phi(bkpt[1])<0:
        return False
    
    #symmetry condition
    for i in range(len(bkpt)):
        x=bkpt[i]
        y=1-x
        if phi(x)+phi(y)!=1:
            return False
    
    #maximality test
    for i in range(len(bkpt)):
        for j in range(len(bkpt)):
            x=bkpt[i]
            y=bkpt[j]
            if phi(x)+phi(y)>phi(x+y):
                return False

    for i in range(len(bkpt)):
        for k in range(len(bkpt)):
            x=bkpt[i]
            z=bkpt[k]
            y=z-x
            if phi(x)+phi(y)>phi(z):
                return False

    return True
         


def phi_bj_1_gdff(c=3/2,periods=3):
    """
    Summary:
        - Name: f_BJ_1;
        - Dim= 1; Slopes = 2; Continuous;
        - Proven maximal for c>=1;
        - Proven extreme for c>=1.
        - quasi periodic function with period 1/c
        - plot periods number of pieces for x<0 and periods number of pieces for x>0

    Parameters:
        c (real) \in [1,positive infinity).

    Examples:
        [1] p.53, Fig 3.1 ::

            sage: logging.disable(logging.INFO)   # Suppress output in automatic tests.
            sage: h=phi_bj_1_gdff(c=3/2)


    References:
   
    - [1]: Alves, Clautiaux, Carvalho, Rietz, Dual-feasible functions for integer programming and combinatorial optimization, 2016.
    """


    n=floor(c)
    beta=c-n
    a=beta/c
    b=1/c-beta/c
    t=1/c
    #positive_slope=1/n/b
    bkpt_start=-periods*t
    bkpt=[bkpt_start]
    value=[-periods/n]
    for i in range(2*periods):
        last=bkpt[-1]
        bkpt.append(last+a)
        bkpt.append(last+t)
        value.append(value[-1])
        value.append(value[-1]+1/n)
    return piecewise_function_from_breakpoints_and_values(bkpt, value)
    



def phi_s_delta(delta=1/10, s=3/2,inf=5):
    """
    Create a family of piecewise linear general DFFs, which have breakpoints: [-inf,-delta, 0, delta, 1-delta,1,1+delta,inf] and slopes: [s,2s,0,1/(1-s*delta),0,2s,s] in each affine piece. The function is maximal if delta is small and s is large.

    """
    s_m=1/(1-2*delta)
    s1=2*s
    bkpt = [-inf,-delta, 0, delta, 1-delta,1,1+delta,inf]
    v1=-delta*s1-s*(inf-delta)
    v2=1+delta*s1+s*(inf-1-delta)
    values =[v1, -delta*s1, 0, 0, 1, 1, 1+delta*s1, v2]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)    

def is_superadditive_almost_strict(delta=1/10,s=3/2,inf=5):
    """
    Check if the special general DFF is almost strictly superadditive, which means that \phi(x+y)-\phi(x)-\phi(y)>t for (x,y) not near three lines: x=0, y=0, x+y=1.
    """
    phi=phi_s_delta(delta, s,inf)
    # Here we choose t=delta as the threshold of strict superadditivity.
    t=delta

    bkpt = [-inf,-delta, 0, delta, 1-delta,1,1+delta,inf]
    for i in range(len(bkpt)):
        for j in range(i,len(bkpt)):
            x=bkpt[i]
            y=bkpt[j]
            z=bkpt[i]+bkpt[j]
            if z<=inf and z>=-inf and (not -delta<x<delta) and (not -delta<y<delta) and (not 1-delta<z<1+delta):
                if phi(x)+phi(y)>phi(z)-t:
                    return False
    
    for i in range(len(bkpt)):
        for k in range(len(bkpt)):
            x=bkpt[i]
            z=bkpt[k]
            y=z-x
            if y<=inf and y>=-inf and (not -delta<x<delta) and (not -delta<y<delta) and (not 1-delta<z<1+delta):
                if phi(x)+phi(y)>phi(z)-t:
                    return False
    return True

def maximality_test_gdff_linear(phi):
    """
    Check whether a given piecewise linear function is a maximal general DFF.
    """ 
    bkpt=phi.end_points()
    inf=-bkpt[0]
    #check if phi(x)>=0 for x>=0
    for i in range(len(bkpt)):
        if bkpt[i]>0:
            if phi.limits_at_end_points()[i][2]<0:
                return False
            else:
                break

    if phi(0)==0 and symmetry_test_gdff_linear(phi) and superadditivity_test_gdff_linear(phi):
        return True
    return False

def symmetry_test_gdff_linear(phi):
    """
    Check whether a given piecewise linear function satisfies \phi(x)+\phi(1-x)=1.
    """
    bkpt=phi.end_points()
    inf=-bkpt[0]
    for i in range(len(bkpt)):
        x=bkpt[i]
        y=1-x
        if y<=inf and y>=-inf:
            if phi(x)+phi(y)!=1:
                return False
    return True


def superadditivity_test_gdff_linear(phi):
    """
    Check whether a given piecewise linear function is superadditive.
    """
    bkpt=phi.end_points()
    inf=-bkpt[0]
    for i in range(len(bkpt)):
        for j in range(i,len(bkpt)):
            x=bkpt[i]
            y=bkpt[j]
            z=bkpt[i]+bkpt[j]
            if z<=inf and z>=-inf:
                if phi(x)+phi(y)>phi(z):
                    return False
    
    for i in range(len(bkpt)):
        for k in range(len(bkpt)):
            x=bkpt[i]
            z=bkpt[k]
            y=z-x
            if y<=inf and y>=-inf:
                if phi(x)+phi(y)>phi(z):
                    return False
    return True

def linear_extension_from_cdff(phi,inf=5):
    """
    Extend a maximal classical DFF to a maximal general DFF using linear extension. Proposition 3.12.
    """
    if not maximality_test_dff(phi):
        logging.info("The function is not maximal")
        return False
    bkpt=phi.end_points()
    limiting_values=phi.limits_at_end_points()
    limiting_slope=1
    for i in range(1,len(bkpt)):
        max_value=max(limiting_values[i])
        t=max_value/bkpt[i]
        if t>limiting_slope:
            limiting_slope=t
    pieces=[[singleton_interval(-inf), FastLinearFunction(0,-inf*limiting_slope+1-limiting_slope)]]
    pieces.append([open_interval(-inf,0), FastLinearFunction(limiting_slope, 1-limiting_slope)])
    for k in range(len(bkpt)-1):
         pieces.append([singleton_interval(bkpt[k]), FastLinearFunction(0,limiting_values[k][0])])
         slope=(limiting_values[k+1][2]-limiting_values[k][1])/(bkpt[k+1]-bkpt[k])
         intercept=limiting_values[k][1]-slope*bkpt[k]
         pieces.append([open_interval(bkpt[k], bkpt[k+1]), FastLinearFunction(slope, intercept)])
    pieces.append([singleton_interval(1), FastLinearFunction(0,1)])
    pieces.append([open_interval(1,inf), FastLinearFunction(limiting_slope, 0)])
    pieces.append([singleton_interval(inf), FastLinearFunction(0,inf*limiting_slope)])
    return FastPiecewise(pieces)
    
def find_0index(phi):
    """
    find the index of 0 in the breakpoints.
    """
    bkpt=phi.end_points()  
    for i in range(len(bkpt)):
        if bkpt[i]==0:
            return i


def two_slope_approximation_gdff_linear(phi,epsilon):
    bkpt=phi.end_points()
    limiting_values=phi.limits_at_end_points()
    inf=-bkpt[0]
    if not maximality_test_gdff_linear(phi):
        logging.info("The function is not restricted maximal")

    bkpt_denominator = [denominator(x) for x in bkpt]
    q=lcm(bkpt_denominator)
    #find limiting slope
    l=len(bkpt)
    s=(limiting_values[l-1][0]-limiting_values[l-2][1])/(bkpt[l-1]-bkpt[l-2])
    delta=min((s-1)/(2*s),1/3,1/q)
    phi_perturb=phi_s_delta(delta,s,inf)
    phi_loose=(1-(epsilon)/(3*(s-1)))*phi+(epsilon)/(3*(s-1)) *phi_perturb

    index0=find_0index(phi_loose)
    bkpt_new=phi_loose.end_points()
    sm=(phi_loose(bkpt_new[index0-1]))/(bkpt_new[index0-1])

    #2slope fill in 
    bkpt_new_denominator = [denominator(x) for x in bkpt_new]
    q1=lcm(bkpt_new_denominator)
    q_new=q1
    while True:
        if q_new%2==0 and q_new>sm/min(epsilon/3,epsilon*delta/9/(s-1)):
            break
        else:
            q_new=q_new+q1
    bkpt_sampling=[-inf]
    for i in range(2*inf*q_new):
        bkpt_sampling.append(bkpt_sampling[-1]+1/q_new)
    pieces=[]
    for i in range(len(bkpt_sampling)-1):
        x=bkpt_sampling[i]
        next_x=bkpt_sampling[i+1]
        y=phi_loose(x)
        next_y=phi_loose(next_x)
        if x<1/2:
            if y==next_y or next_y==y+sm/q_new:
                pieces.append(closed_piece((x, y), (next_x, next_y)))  
            else:    
                mx=next_x-(next_y-y)/sm
                my=y
                pieces.append(closed_piece((x, y), (mx, my))) 
                pieces.append(closed_piece((mx, my), (next_x, next_y))) 
        else:
            if y==next_y or next_y==y+sm/q_new:
                pieces.append(closed_piece((x, y), (next_x, next_y)))  
            else:
                mx=x+(next_y-y)/sm
                my=next_y
                pieces.append(closed_piece((x, y), (mx, my))) 
                pieces.append(closed_piece((mx, my), (next_x, next_y)))     
    return FastPiecewise(pieces)
    
    
        







