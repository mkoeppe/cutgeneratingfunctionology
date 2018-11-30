# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

def phi_forward_3_slope(b=4/5, lambda_1=4/9, lambda_2=2/3, field=None):
    """
    Summary:
        - Dim = 1; Slopes = 3; Continuous. 

    Parameters:
        b (real) \in (0,positive infinity) and not an integer;
        lambda_1, lambda_2 (real) in (0,1].

    Function is known to be extreme under the conditions:
        0 <= lambda_1 <= 1/2, 0 <= lambda_2 <= 1, b>3, 0 < lambda_1 * f + lambda_2 * (f - 1) < lambda_1 * f, f=frac(b).

    Examples:
        sage: logging.disable(logging.INFO)   # Suppress output in automatic tests.
        sage: h=phi_forward_3_slope()
        sage: extremality_test_dff(h)
        False

    """
    n=floor(b)
    f=b-n
    a = lambda_1 * f / 2
    a1 = a + lambda_2 * (f - 1) / 2
    if not bool(0 < f < 1) & bool(0 < a1 < a < f / 2):
        raise ValueError, "Bad parameters. Unable to construct the function."    
    bkpts = [0, a1, a, f - a, f - a1, f]
    values = [0, (lambda_1 + lambda_2)/2, lambda_1 / 2, 1 - lambda_1 / 2, 1 - (lambda_1 + lambda_2)/2, 1]
    bkpt=[]
    for i in range(n+1):
        for j in range(len(bkpts)):
            bkpt.append((i+bkpts[j])/b)
    value=[]
    la=(2*a1)/(lambda_1 + lambda_2)
    for i in range(len(bkpt)):
        value.append((b*bkpt[i]-la*values[i%len(bkpts)])/(b-la))
    return piecewise_function_from_breakpoints_and_values(bkpt, value, field=field)


def phi_2_slope(b=3/5, lambda_1=1/6, field=None):
    """
    Summary:
        - Dim = 1; Slopes = 2; Continuous. 

    Parameters:
        b (real) \in (0,positive infinity) and not an integer;
        lambda_1 (real) in (0,1].

    Function is known to be extreme under the conditions:
        0 < lambda_1 <=1, lambda_1 < f/(1 - f), b>3, f=frac(b).

    Examples:
        sage: logging.disable(logging.INFO)   # Suppress output in automatic tests.
        sage: h=phi_2_slope()
        sage: extremality_test_dff(h)
        False

    """
    n=floor(b)
    f=b-n
    if not (bool(0 < f < 1) & bool(0 < lambda_1 < f/(1 - f))) & bool(lambda_1<=1):
        raise ValueError, "Bad parameters. Unable to construct the function."
    la=f-lambda_1/(1+lambda_1)
    bkpt=[]
    bkpts = [0, (f - lambda_1*(1 - f))/2, (f + lambda_1*(1 - f))/2, f]
    for i in range(n+1):
        for j in range(len(bkpts)):
            bkpt.append((i+bkpts[j])/b)
    values=[]
    value_1 = [0, (1 + lambda_1)/2, (1 - lambda_1)/2, 1]
    for i in range(len(bkpt)):
        values.append((b*bkpt[i]-la*value_1[i%len(bkpts)])/(b-la))
    return piecewise_function_from_breakpoints_and_values(bkpt, values, field=field)

def w_2slope_3covered_nonextreme():
    """
    A continuous 2-slope nonextreme function with 3 covered components. 
    """
    return FastPiecewise([[(QQ(0), 1/10), FastLinearFunction(QQ(0), QQ(0))], [(1/10, 3/20), FastLinearFunction(5/2, -1/4)], [(3/20, 1/4), FastLinearFunction(QQ(0), 1/8)], [(1/4, 7/20), FastLinearFunction(5/2, -1/2)], [(7/20, 2/5), FastLinearFunction(QQ(0), 3/8)], [(2/5, 9/20), FastLinearFunction(5/2, -5/8)], [(9/20, 11/20), FastLinearFunction(QQ(0), 1/2)], [(11/20, 3/5), FastLinearFunction(5/2, -7/8)], [(3/5, 13/20), FastLinearFunction(QQ(0), 5/8)], [(13/20, 3/4), FastLinearFunction(5/2, -QQ(1))], [(3/4, 17/20), FastLinearFunction(QQ(0), 7/8)], [(17/20, 9/10), FastLinearFunction(5/2, -5/4)], [(9/10, QQ(1)), FastLinearFunction(QQ(0), QQ(1))]])  

def w_2slope_3covered():
    """
    A continuous 2-slope extreme function with 3 covered components.
    """
    return FastPiecewise([[(QQ(0), 1/14), FastLinearFunction(QQ(0), QQ(0))], [(1/14, 3/28), FastLinearFunction(7/3, -1/6)], [(3/28, 5/28), FastLinearFunction(QQ(0), 1/12)], [(5/28, 1/4), FastLinearFunction(7/3, -1/3)], [(1/4, 2/7), FastLinearFunction(QQ(0), 1/4)], [(2/7, 9/28), FastLinearFunction(7/3, -5/12)], [(9/28, 11/28), FastLinearFunction(QQ(0), 1/3)], [(11/28, 3/7), FastLinearFunction(7/3, -7/12)], [(3/7, 13/28), FastLinearFunction(QQ(0), 5/12)], [(13/28, 15/28), FastLinearFunction(7/3, -2/3)], [(15/28, 4/7), FastLinearFunction(QQ(0), 7/12)], [(4/7, 17/28), FastLinearFunction(7/3, -3/4)], [(17/28, 19/28), FastLinearFunction(QQ(0), 2/3)], [(19/28, 5/7), FastLinearFunction(7/3, -11/12)], [(5/7, 3/4), FastLinearFunction(QQ(0), 3/4)], [(3/4, 23/28), FastLinearFunction(7/3, -QQ(1))], [(23/28, 25/28), FastLinearFunction(QQ(0), 11/12)], [(25/28, 13/14), FastLinearFunction(7/3, -7/6)], [(13/14, QQ(1)), FastLinearFunction(QQ(0), QQ(1))]])


def phi_bj_1(c=3/2):
    """
    Summary:
        - Name: f_BJ_1;
        - Dim= 1; Slopes = 2; Continuous;
        - Proven maximal [1] p.26, theorem 2.1;
        - Proven extreme for c>=2 and non-extreme for 1<c<2.
    
    Parameters:
        c (real) \in [1,positive infinity).

    Examples:
        [1] p.25, Fig 2.4 ::

            sage: logging.disable(logging.INFO)   # Suppress output in automatic tests.
            sage: h=phi_bj_1(3/2)
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

def phi_simple(c=3/2):
    """
    Summary:
        - Dim= 1; Slopes = 1; Discontinuous;
        - Not maximal.     
    
    Parameters:
        c (real) \in [1,positive infinity).

    Examples:
        [1] p.22, Fig 2.1 ::

            sage: logging.disable(logging.INFO)   # Suppress output in automatic tests.
            sage: h=phi_simple(3/2)
            sage: superadditive_test(h)
            True

    References:
   
    - [1]: Alves, Clautiaux, Carvalho, Rietz, Dual-feasible functions for integer programming and combinatorial optimization, 2016.
    """

    n=floor(c)
    l=1/c
    pieces=[]
    for k in range(n):
        pieces = pieces + \
                 [[singleton_interval(l*k), FastLinearFunction(0,k/n)], \
                  [open_interval(l*k,l*k+l), FastLinearFunction(0,k/n)]]
    pieces.append([closed_interval(n/c,1), FastLinearFunction(0,1)])
    return FastPiecewise(pieces)

def phi_ccm_1(c=3/2):
    """
    Summary:
        - Name: f_CCM_1;
        - Dim= 1; Slopes = 1; Discontinuous;
        - Always extreme.     
    
    Parameters:
        c (real) \in [1,positive infinity).

    Examples:
        [1] p.29, Fig 2.5 ::
            
            sage: logging.disable(logging.INFO)   # Suppress output in automatic tests.
            sage: h=phi_ccm_1(3/2)
            sage: extremality_test_dff(h)
            True

    References:
   
    - [1]: Alves, Clautiaux, Carvalho, Rietz, Dual-feasible functions for integer programming and combinatorial optimization, 2016.
    """

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

def phi_fs_1(k=3):
    """
    Summary:
        - Name: f_FS_1;
        - Dim= 1; Slopes = 1; Discontinuous.    
    
    Parameters:
        k (positive integer).

    Examples:
        [1] p.35, Fig 2.6 ::
            
            sage: logging.disable(logging.INFO)   # Suppress output in automatic tests.
            sage: h=phi_fs_1(3)
            sage: extremality_test_dff(h)
            True

    References:
   
    - [1]: Alves, Clautiaux, Carvalho, Rietz, Dual-feasible functions for integer programming and combinatorial optimization, 2016.
    """

    pieces=[]
    pieces.append([singleton_interval(0),FastLinearFunction(0,0)])
    for i in range(k+1):
        pieces=pieces + \
               [[open_interval(i/(k+1),(i+1)/(k+1)),FastLinearFunction(0,i/k)],\
                 [singleton_interval((i+1)/(k+1)),FastLinearFunction(0,(i+1)/(k+1))]]
    return FastPiecewise(pieces)

def phi_vb_2(k=3):
    """
    Summary:
        - Name: f_VB_2;
        - Dim= 1; Slopes = 1; Discontinuous;  
        - Always maximal.  
    
    Parameters:
        k (positive integer).

    Examples:
        [1] p.40, Fig 2.8 ::
            
            sage: logging.disable(logging.INFO)   # Suppress output in automatic tests.
            sage: h=phi_vb_2(3)
            sage: extremality_test_dff(h)
            True

    References:
   
    - [1]: Alves, Clautiaux, Carvalho, Rietz, Dual-feasible functions for integer programming and combinatorial optimization, 2016.
    """

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

def phi_ll_1(c=3/2,k=5):
    """
    Summary:
        - Name: f_LL_1;
        - Dim= 1; Slopes = 1; Discontinuous.    
    
    Parameters:
        k (positive integer);
        c (real) \in [1,positive infinity).

    Examples:
        [1] p.42, Fig 2.9 ::
            
            sage: logging.disable(logging.INFO)   # Suppress output in automatic tests.
            sage: h=phi_ll_1(3/2, 5)
            sage: maximality_test_dff(h)
            False

    References:
   
    - [1]: Alves, Clautiaux, Carvalho, Rietz, Dual-feasible functions for integer programming and combinatorial optimization, 2016.
    """


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

def phi_ll_2(c=3/2,k=5):
    """
    Summary:
        - Name: f_LL_2;
        - Dim= 1; Slopes = 1; Discontinuous.    
    
    Parameters:
        k (positive integer);
        c (real) \in [1,positive infinity).

    Function is know to be maximal under the conditions:
        c is not an integer and k >= ceil(1/frac(c))

    Examples:
        [1] p.42, Fig 2.9 ::
            
            sage: logging.disable(logging.INFO)   # Suppress output in automatic tests.
            sage: h=phi_ll_2(3/2, 5)
            sage: maximality_test_dff(h)
            True

    References:
   
    - [1]: Alves, Clautiaux, Carvalho, Rietz, Dual-feasible functions for integer programming and combinatorial optimization, 2016.
    """

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
                pieces.append([right_open_interval(1-m/c-l1-(d-t-1)*l2,1-m/c-l1-(d-t-1)*l2+l2), FastLinearFunction(0,1-(m+(d-t-1)/k)/n)])
        pieces.append([closed_interval(1-m/c-l1,1-m/c), FastLinearFunction(0,1-m/n)])
    for i in range(m):        
        for j in range(k-1):
            if j==0:
               pieces.append([open_interval(1-(m-i)/c+j*l2,1-(m-i)/c+(j+1)*l2), FastLinearFunction(0,1-(m-i-(j+1)/k)/n)])
            else:
               pieces.append([right_open_interval(1-(m-i)/c+j*l2,1-(m-i)/c+(j+1)*l2), FastLinearFunction(0,1-(m-i-(j+1)/k)/n)])            
        pieces.append([closed_interval(1-(m-i-1)/c-l1,1-(m-i-1)/c), FastLinearFunction(0,1-(m-1-i)/n)])    
    return FastPiecewise(pieces)

def phi_dg_1(c=3/2,k=5):
    """
    Summary:
        - Name: f_DG_1;
        - Dim= 1; Slopes = 1; Discontinuous.    
    
    Parameters:
        k (positive integer);
        c (real) \in [1,positive infinity).

    Function is know to be maximal under the conditions:
        c is not an integer and k >= ceil(1/frac(c))

    Examples:
        [1] p.44, Fig 2.11 ::
            
            sage: logging.disable(logging.INFO)   # Suppress output in automatic tests.
            sage: h=phi_dg_1(3/2, 5)
            sage: extremality_test_dff(h)
            False

    References:
   
    - [1]: Alves, Clautiaux, Carvalho, Rietz, Dual-feasible functions for integer programming and combinatorial optimization, 2016.
    """

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



