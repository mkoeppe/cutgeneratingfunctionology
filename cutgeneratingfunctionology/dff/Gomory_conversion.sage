from six.moves import range
def conversion_from_DFF_to_GJ(phi):
    bkpt=phi.end_points()
    value=phi.values_at_end_points()
    slope=[0]
    largest_slope_index=[]
    for i in range(len(bkpt)):
        if i>0:
            slope.append(value[i]/bkpt[i])
    largest_slope=max(slope)
    for i in range(len(slope)):
        if slope[i]==largest_slope:
            largest_slope_index.append(i)
    n=len(largest_slope_index)
    i0=largest_slope_index[0]
    b=1/bkpt[i0]
    f=b-n
    bkpt_pi=[]
    for i in range(i0+1):
        bkpt_pi.append(bkpt[i]*b)
    x=phi(f/b)
    la=(f-b*x)/(1-x)
    value_pi=[]
    for i in range(i0+1):
        value_pi.append((bkpt_pi[i]-(b-la)*value[i])/(la))
    pi=piecewise_function_from_breakpoints_and_values(bkpt_pi, value_pi)
    return pi
    
            


def possible_extreme_DFF(pi,n=1,left=0,right=1):
    bkpt=pi.end_points()
    f=find_f(pi)
    if pi.limits(0)[1]!= 0:
        logging.info("Not possible since discontinuous at 0.")
        return False
    s=pi.limits(bkpt[1])[-1]/bkpt[1]
    if n==0 and f*s==1:
        logging.info("The function is x on [0,1], and not extreme as general DFF.")
        return False
    la=1/s
    return conversion_from_GJ_to_DFF(pi,la,left,right,n)

def conversion_from_GJ_to_DFF(pi,la=0,left=0,right=1,n=1):
    f=find_f(pi)
    limit=pi.limits_at_end_points()
    bkpt=pi.end_points()
    if pi.limits(0)[1]!= 0:
        logging.info("Not possible since discontinuous at 0.")
        return False
    bkpt1=[]
    bkpt2=[]
    b=n+f
    l=1/b
    k1=floor(right/l)
    k2=ceil(left/l)
    pieces=[]
    if left != k2*l:
        left1=left/l-k2+1
        for i in range (len(bkpt)-1):
            if bkpt[i]<= left1 <bkpt[i+1]:
                break
        bkpt1.append(left)
        for j in range(len(bkpt)-i-1):    
            bkpt1.append((bkpt[i+j+1]-1+k2)*l)
        for s in range(len(bkpt1)-1):
            x=bkpt1[s]
            y=bkpt1[s+1]
            vx=(b*x-la*pi.limits(b*x-floor(b*x))[1])/(b-la)
            vy=(b*y-la*pi.limits(b*y-floor(b*y))[-1])/(b-la)
            slope=(vy-vx)/(y-x)
            intercept=vy-slope*y
            pieces = pieces + \
                 [[singleton_interval(x), FastLinearFunction(0,compute_value(pi,b,la,x))], \
[open_interval(x,y), FastLinearFunction(slope,intercept)]]
    for t in range(k1-k2):
        k=k2+t
        for i in range(len(bkpt)-1):
            x=(k+bkpt[i])/b
            y=(k+bkpt[i+1])/b
            vx=(b*x-la*limit[i][1])/(b-la)
            vy=(b*y-la*limit[i+1][-1])/(b-la)
            slope=(vy-vx)/(y-x)
            intercept=vy-slope*y
            pieces = pieces + \
                 [[singleton_interval(x), FastLinearFunction(0,compute_value(pi,b,la,x))], \
                  [open_interval(x,y), FastLinearFunction(slope,intercept)]]
    if right != k1*l:
        right1=right/l-k1
        for i in range (len(bkpt)-1):
            if bkpt[i]< right1 <=bkpt[i+1]:
                break
        for j in range (i+1):
            bkpt2.append((bkpt[j]+k1)*l)
        bkpt2.append(right)
        for s in range(len(bkpt2)-1):
            x=bkpt2[s]
            y=bkpt2[s+1]
            vx=(b*x-la*pi.limits(b*x-floor(b*x))[1])/(b-la)
            vy=(b*y-la*pi.limits(b*y-floor(b*y))[-1])/(b-la)
            slope=(vy-vx)/(y-x)
            intercept=vy-slope*y
            pieces = pieces + \
                 [[singleton_interval(x), FastLinearFunction(0,compute_value(pi,b,la,x))], \
                  [open_interval(x,y), FastLinearFunction(slope,intercept)]]
    pieces.append([singleton_interval(right), FastLinearFunction(0,compute_value(pi,b,la,right))])
    h=FastPiecewise(pieces)
    bkpt0=h.end_points()
    for i in range(len(bkpt0)):
        if bkpt0[i]==0:
            break
    if h.limits(bkpt0[i+1])[-1]<0:
        logging.info("Lambda is too large.")
        return False    
    return h


def compute_value(pi,b,la,x):
    return (b*x-la*pi(b*x-floor(b*x)))/(b-la)

    