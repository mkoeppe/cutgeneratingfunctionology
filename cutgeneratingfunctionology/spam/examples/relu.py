r"""
FM elimination examples: relu and clipped relu using BasicSemialgebraicSet_polyhedral_linear_system
"""

from __future__ import division, print_function, absolute_import

from cutgeneratingfunctionology.spam.basic_semialgebraic_linear_system import BasicSemialgebraicSet_polyhedral_linear_system

import itertools

def FM_relu_1d(base_ring, poly_ring):
    """
    One dimensional relu y=max(0,x). Original variables (x0,x1,x,y,z).
    
    EXAMPLE::
    
        sage: import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *
        sage: from cutgeneratingfunctionology.spam.examples.relu import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: K.<L,U,W,b>=ParametricRealField([QQ(-2),QQ(2),QQ(2),QQ(1/2)])
        sage: Q.<x0,x1,x,y,z>=K[]
        sage: bsa = FM_relu_1d(K,Q)
        sage: bsa.le_poly()
        {-y,
        (((-1)/W)~)*y + (((L*W + b)/W)~)*z,
        (1/W)~*y + (((-U*W - b)/W)~)*z,
        -x + (1/W)~*y + (((-L*W - b)/W)~)*z + L~,
        x + (((-1)/W)~)*y + (((U*W + b)/W)~)*z + (-U)~,
        W~*x - y + b~}
    """
    x0,x1,x,y,z = iter(poly_ring.gens())
    L,U,W,b = iter(base_ring.gens())
    le = [-y, W*x0-b*z+b, x0+U*z-U, -x0-L*z+L, x1-U*z, -x1+L*z]
    eq = [x0+x1-x, W*x1-y+b*z]
    bsa = BasicSemialgebraicSet_polyhedral_linear_system(poly_ring=poly_ring, le=le, eq=eq)
    bsa_eliminated = bsa.coordinate_projection([x0,x1])
    return bsa_eliminated

'''
def initialize_1d_clipped_relu_parameters():
    K.<L,U,W,b,C>=ParametricRealField([QQ(-2),QQ(2),QQ(1),QQ(1/2),QQ(2)])
    param={}
    param['L']=[L]
    param['U']=[U]
    param['W']=[W]
    param['b']=b
    param['C']=C
    return K, param

def initialize_2d_clipped_relu_parameters():
    K.<L1,U1,L2,U2,W1,W2,b,C>=ParametricRealField([QQ(-2),QQ(2),QQ(-1),QQ(3),QQ(1),QQ(2),QQ(1/2),QQ(2)])
    param={}
    param['L']=[L1,L2]
    param['U']=[U1,U2]
    param['W']=[W1,W2]
    param['b']=b
    param['C']=C
    return K, param

def find_new_inequalities_clipped_relu(param):
    L,U,W,b,C=param['L'],param['U'],param['W'],param['b'],param['C']
    n=len(W)

    # assumptions
    for i in range(n):
        assert L[i]<0
        assert U[i]>0
        assert W[i]>0

    assert C>0
    assert b+sum(W[i]*U[i] for i in range(n))>C
    assert b+sum(W[i]*L[i] for i in range(n))<0

    K.freeze()

    # know exponential family valid inequalities
    exponential_valid_inequalities=clipped_relu_known_valid_inqualities(param)
    exponential_valid_inequalities=naive_simplify(exponential_valid_inequalities,binary_variables=3)
    exponential_valid_inequalities=[normalize(row) for row in exponential_valid_inequalities]


    # inequalities after FM elimination
    A=clipped_relu_extended_matrix(param)
    S=initialize_fourier_system(A)
    for i in range(3*n):
        S.pivot()
    FM_matrix=naive_simplify(S.matrix,binary_variables=3)
    FM_matrix=[normalize(row) for row in FM_matrix]

    res=[]
    for i1 in FM_matrix:
        if any([check_redundancy_one_to_one(i1,i2,3) for i2 in exponential_valid_inequalities]):
            continue
        else:
            res.append(i1)
    return res

def initialize_fourier_system(A):
    M=[]
    for i in range(len(A)):
        M.append(InequalityRow(A[i],{i}))
    return FourierSystem(M)

def check_redundancy_one_to_one(i1,i2,binary_variables=3):
    """
    Given two inequalities i1,i2 (<=0), check whether i1 is redundant compared to i2.
    We assume all entries of i1 and i2 are the same except the binary variables, otherwise raise exception.
    """
    for i in range(len(i1)-4):
        try:
            if not i1[i]==i2[i]:
                return False
        except ParametricRealFieldFrozenError:
            return False
    for i in range(binary_variables+1):
        try:
            if i1[-i-1]>i2[-i-1]:
                return False
        except ParametricRealFieldFrozenError:
            return False
    return True

def normalize(row):
    for i in range(len(row)):
        if row[i]!=0:
            break
    if row[i]==0:
        return row
    pivot=row[i]
    for j in range(i,len(row)):
        if pivot>0:
            row[j]=row[j]/pivot
        else:
            row[j]=-row[j]/pivot
    return row

def naive_simplify(A,binary_variables=3):
    """
    The last binary_variables+1 columns represent all binary variables z_i, and the constant term.
    If all previous entries of a row in A is 0, drop the row.
    """
    res=[]
    for row in A:
        if row[:-binary_variables-1].any():
            res.append(row)
    return np.asarray(res)

def FM_relu_2d(**kwds):
    """
    Two dimensional relu. (x0_1,x1_1,x0_2,x1_2,x_1,x_2,y,z,1).
    """
    K.<L1,U1,L2,U2,W1,W2,b>=ParametricRealField([QQ(-2),QQ(2),QQ(-1),QQ(3),QQ(1),QQ(2),QQ(1/2)])
    A=[]
    A.append([1,1,0,0,-1,0,0,0,0])
    A.append([-1,-1,0,0,1,0,0,0,0])
    A.append([0,0,1,1,0,-1,0,0,0])
    A.append([0,0,-1,-1,0,-1,0,0,0])
    A.append([W1,0,W2,0,0,0,0,-b,b])
    A.append([0,0,0,0,0,0,-1,0,0])
    A.append([0,W1,0,W2,0,0,-1,b,0])
    A.append([0,-W1,0,-W2,0,0,1,-b,0])
    A.append([1,0,0,0,0,0,0,U1,-U1])
    A.append([0,0,1,0,0,0,0,U2,-U2])
    A.append([-1,0,0,0,0,0,0,-L1,L1])
    A.append([0,0,-1,0,0,0,0,-L2,L2])
    A.append([0,1,0,0,0,0,0,-U1,0])
    A.append([0,0,0,1,0,0,0,-U2,0])
    A.append([0,-1,0,0,0,0,0,L1,0])
    A.append([0,0,0,-1,0,0,0,L2,0])

    A=np.asarray(A)
    A1=FM_symbolic(A,**kwds)
    A2=FM_symbolic(A1,**kwds)
    A3=FM_symbolic(A2,**kwds)
    A4=FM_symbolic(A3,**kwds)
    return naive_simplify(A4,binary_variables=1)

def clipped_relu_extended_matrix(param):
    """
    Return Balas's formulation matrix for n dimensional clipped relu. (x1_1,...x1_n,x2_1,...x2_n,x3_1,...x3_n,x1,...,xn,y,z1,z2,z3,1)
    Given parameters L,U,W,b,C.
    """
    L,U,W,b,C=param['L'],param['U'],param['W'],param['b'],param['C']
    n=len(W)
    A=[]

    # z1+z2+z3=1
    A.append([0]*4*n+[0,1,1,1,-1])
    A.append([0]*4*n+[0,-1,-1,-1,1])

    # x1+x2+x3=x
    for i in range(n):
        A.append(unit_vector(n,i,False)*3+unit_vector(n,i,True)+[0]*5)
        A.append(unit_vector(n,i,True)*3+unit_vector(n,i,False)+[0]*5)

    # Wx1+bz1<=0
    A.append(W+[0]*3*n+[0,b,0,0,0])

    # Wx3+bz3>=Cz3
    A.append([0]*2*n+negation(W)+[0]*n+[0,0,0,C-b,0])

    # 0<=y-Cz3=Wx2+bz2<=Cz2
    A.append([0]*4*n+[-1,0,0,C,0])
    A.append([0]*4*n+[1,0,-C,-C,0])
    A.append([0]*n+W+[0]*2*n+[-1,0,b,C,0])
    A.append([0]*n+negation(W)+[0]*2*n+[1,0,-b,-C,0])

    # L<=x<=U
    for i in range(n):
        A.append(unit_vector(4*n,i,True)+[0,-U[i],0,0,0])
        A.append(unit_vector(4*n,i,False)+[0,L[i],0,0,0])
        A.append(unit_vector(4*n,n+i,True)+[0,0,-U[i],0,0])
        A.append(unit_vector(4*n,n+i,False)+[0,0,L[i],0,0])
        A.append(unit_vector(4*n,2*n+i,True)+[0,0,0,-U[i],0])
        A.append(unit_vector(4*n,2*n+i,False)+[0,0,0,L[i],0])
    return np.asarray(A)

def clipped_relu_known_valid_inqualities(param):
    """
    Return known valid inequalities for n dimensional clipped relu. (x1,...,xn,y,z1,z2,z3,1)
    Given parameters L,U,W,b,C.
    """
    L,U,W,b,C=param['L'],param['U'],param['W'],param['b'],param['C']
    n=len(W)
    A=[]

    # z1+z2+z3=1
    A.append([0]*n+[0,1,1,1,-1])
    A.append([0]*n+[0,-1,-1,-1,1])

    # Cz3<=y<=Cz2+Cz3
    A.append([0]*n+[1,0,-C,-C,0])
    A.append([0]*n+[-1,0,0,C,0])

    # L<=x<=U
    for i in range(n):
        A.append(unit_vector(n,i,True)+[0,-U[i],-U[i],-U[i],0])
        A.append(unit_vector(n,i,False)+[0,L[i],L[i],L[i],0])

    # exponential family of strengthening inequalities.
    for I in generate_all_binary_vectors(n):
        I_bar=complement_index(I)
        x_coefs=[I[i]*W[i] for i in range(n)]
        nx_coefs=[-t for t in x_coefs]
        xbar_coefs=[I_bar[i]*W[i] for i in range(n)]
        nxbar_coefs=[-t for t in xbar_coefs]
        IU=sum(x_coefs[i]*U[i] for i in range(n))
        IL=sum(x_coefs[i]*L[i] for i in range(n))
        IbarU=sum(xbar_coefs[i]*U[i] for i in range(n))
        IbarL=sum(xbar_coefs[i]*L[i] for i in range(n))

        # type 1 facet-defining inequalities
        A.append(nx_coefs+[1,IL,-b-IbarU,-b-IbarU,0])
        A.append(x_coefs+[-1,b+IbarL,b+IbarL,C-IU,0])

        # type 2 facet-defining inequalities
        A.append(x_coefs+[0,b+IbarL,-IU,-IU,0])
        A.append(nx_coefs+[0,IL,IL,C-b-IbarU,0])

        # type 3 redundant inequalities
        A.append(nx_coefs+[1,IL,-b-IbarU,IL-C,0])
        A.append(x_coefs+[-1,-IU,b+IbarL,C-IU,0])

    return np.asarray(A)

def FM_clipped_relu_1d(**kwds):
    """
    One dimensional relu. (x1,x2,x3,x,y,z1,z2,z3,1).
    """
    K.<L,U,W,b,C>=ParametricRealField([QQ(-2),QQ(2),QQ(1),QQ(1/2),QQ(2)])
    A=[]
    A.append([0,0,0,0,0,1,1,1,-1])
    A.append([0,0,0,0,0,-1,-1,-1,1])
    A.append([1,0,0,0,0,-U,0,0,0])
    A.append([-1,0,0,0,0,L,0,0,0])
    A.append([0,1,0,0,0,0,-U,0,0])
    A.append([0,-1,0,0,0,0,L,0,0])
    A.append([0,0,1,0,0,0,0,-U,0])
    A.append([0,0,-1,0,0,0,0,L,0])
    A.append([1,1,1,-1,0,0,0,0,0])
    A.append([-1,-1,-1,1,0,0,0,0,0])
    A.append([W,0,0,0,0,b,0,0,0])
    A.append([0,0,-W,0,0,0,0,C-b,0])
    A.append([0,0,0,0,-1,0,0,C,0])
    A.append([0,0,0,0,1,0,-C,-C,0])
    A.append([0,W,0,0,-1,0,b,C,0])
    A.append([0,-W,0,0,1,0,-b,-C,0])

    A=np.asarray(A)
    A1=FM_symbolic(A,**kwds)
    A2=FM_symbolic(A1,**kwds)
    A3=FM_symbolic(A2,**kwds)
    return naive_simplify(A3,binary_variables=3)

def FM_clipped_relu_2d(**kwds):
    """
    Two dimensional relu. (x1_1,x1_2,x2_1,x2_2,x3_1,x3_2,x1,x2,y,z1,z2,z3,1).
    """
    K.<L1,U1,L2,U2,W1,W2,b,C>=ParametricRealField([QQ(-2),QQ(2),QQ(-1),QQ(3),QQ(1),QQ(2),QQ(1/2),QQ(2)])
    A=[]
    A.append([0,0,0,0,0,0,0,0,0,1,1,1,-1])
    A.append([0,0,0,0,0,0,0,0,0,-1,-1,-1,1])
    A.append([1,0,0,0,0,0,0,0,0,-U1,0,0,0])
    A.append([-1,0,0,0,0,0,0,0,0,L1,0,0,0])
    A.append([0,1,0,0,0,0,0,0,0,-U2,0,0,0])
    A.append([0,-1,0,0,0,0,0,0,0,L2,0,0,0])
    A.append([0,0,1,0,0,0,0,0,0,0,-U1,0,0])
    A.append([0,0,-1,0,0,0,0,0,0,0,L1,0,0])
    A.append([0,0,0,1,0,0,0,0,0,0,-U2,0,0])
    A.append([0,0,0,-1,0,0,0,0,0,0,L2,0,0])
    A.append([0,0,0,0,1,0,0,0,0,0,0,-U1,0])
    A.append([0,0,0,0,-1,0,0,0,0,0,0,L1,0])
    A.append([0,0,0,0,0,1,0,0,0,0,0,-U2,0])
    A.append([0,0,0,0,0,-1,0,0,0,0,0,L2,0])

    A.append([1,0,1,0,1,0,-1,0,0,0,0,0,0])
    A.append([-1,0,-1,0,-1,0,1,0,0,0,0,0,0])
    A.append([0,1,0,1,0,1,0,-1,0,0,0,0,0])
    A.append([0,-1,0,-1,0,-1,0,1,0,0,0,0,0])

    A.append([W1,W2,0,0,0,0,0,0,0,b,0,0,0])
    A.append([0,0,0,0,-W1,-W2,0,0,0,0,0,C-b,0])
    A.append([0,0,0,0,0,0,0,0,-1,0,0,C,0])
    A.append([0,0,0,0,0,0,0,0,1,0,-C,-C,0])
    A.append([0,0,W1,W2,0,0,0,0,-1,0,b,C,0])
    A.append([0,0,-W1,-W2,0,0,0,0,1,0,-b,-C,0])

    A=np.asarray(A)
    A1=FM_symbolic(A,**kwds)
    A2=FM_symbolic(A1,**kwds)
    A3=FM_symbolic(A2,**kwds)
    A4=FM_symbolic(A3,**kwds)
    A5=FM_symbolic(A4,**kwds)
    A6=FM_symbolic(A5,**kwds)
    return naive_simplify(A6,binary_variables=3)

def FM_symbolic(A,check_feasible_redundancy=False):
    """
    Eliminate the first variable of the system defined by Ax<=0.
    The last column of A corresponds to constant term.
    """
    equ=[]
    low=[]
    upp=[]
    A1=[]
    for i in range(len(A)):
        row=A[i]
        if row[0]>0:
            upp.append(i)
        elif row[0]<0:
            low.append(i)
        else:
            equ.append(i)
    for e in equ:
        A1.append(A[e][1:])
    for l in low:
        for u in upp:
            positive_row=A[u]
            negative_row=A[l]
            new_row=positive_row[1:]*(-negative_row[0])+negative_row[1:]*positive_row[0]
            if check_feasible_redundancy:
                if not FM_check_feasible(new_row):
                    print('The system is not feasible.')
                    return False
                if FM_check_redundancy(new_row):
                    continue
                else:
                    A1.append(new_row)
            else:
                A1.append(new_row)
    A1=np.asarray(A1)
    return A1

def FM_check_redundancy(row):
    if np.all(row[:-1]==0) and row[-1]<=0:
        return True
    else:
        return False

def FM_check_feasible(row):
    if np.all(row[:-1]==0) and row[-1]>0:
        return False
    else:
        return True

def generate_all_binary_vectors(n):
    return list(itertools.product([0, 1], repeat=n))

def complement_index(I):
    return [1-i for i in I]

def unit_vector(dim,non_zero,positive=True):
    res=[0]*dim
    if positive:
        res[non_zero]=1
    else:
        res[non_zero]=-1
    return res

def negation(a):
    return [-v for v in a]

'''
