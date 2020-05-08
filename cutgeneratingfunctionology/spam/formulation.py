import numpy as np
import itertools
import random

from cutgeneratingfunctionology.igp import ParametricRealFieldFrozenError

class FourierSystem :
    
    r"""
    Class for FM elimination system using matrix A representing the set Ax<=0.
    The last column respresent the constant term.
    """

    def __init__(self, field, matrix, history_set, binary_variables = 0):
        r"""
        Initialize the system for FM elimination.
        The parameter remove_binary_only is specific for relu/clipped relu formulation.
        If remove_binary_only is True, then those rows with only binary variables will be removed.
        The binary variables always appear at the end.
        
        Examples::
        
            sage: import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *
            sage: from cutgeneratingfunctionology.spam.formulation import *
            sage: logging.disable(logging.INFO)
            sage: K.<L1,U1,L2,U2,W1,W2,b>=ParametricRealField([QQ(-2),QQ(2),QQ(-1),QQ(3),QQ(1),QQ(2),QQ(1/2)], allow_coercion_to_float=False)
            sage: FS=FM_relu_2d(K)
        """
        if len(matrix) == 0:
            raise ValueError("empty matrix")
        if len(matrix) != len(history_set):
            raise ValueError("Size of matrix and history_set doesn't match.")
        self.field=field
        self.matrix=matrix
        self.history_set=history_set
        self.binary_variables = binary_variables
        # Only normalize rows and remove binary only rows if no future elimination is needed.
        #if remove_binary_only:
            # with self.field.off_the_record():
                #self.remove_binary_variable_only_rows(binary_variables)
        #with self.field.off_the_record():
            #self.normalize()

    def normalize(self):
        r"""
        Normalize every row based on the first non-zero entry.
        """
        matrix=self.matrix.copy()
        for i in range(len(matrix)):
            row=matrix[i]
            for j in range(len(row)):
                if row[j]!=0:
                    break
            if row[j]==0:
                continue
            pivot=row[j]
            for k in range(j,len(row)):
                if pivot>0:
                    row[k]=row[k]/pivot
                else:
                    row[k]=-row[k]/pivot
            self.matrix[i]=row

    def remove_binary_variable_only_rows(self, binary_variables):
        r"""
        Remove those rows with only binary variables.
        """
        nonredundancy_list=[]
        for i in range(len(self.matrix)):
            row = self.matrix[i]
            if row[:-binary_variables-1].any():
                nonredundancy_list.append(i)
        self.matrix=np.asarray([self.matrix[i] for i in range(len(self.matrix)) if i in nonredundancy_list])
        self.history_set=[self.history_set[i] for i in range(len(self.history_set)) if i in nonredundancy_list]

    def remove_non_minimal_history_set(self):
        r"""
        Remove those rows with non minimal history_set.
        """
        redundancy_list=[]
        for i in range(len(self.history_set)):
            for j in range(len(self.history_set)):
                if i==j:
                    continue
                if self.history_set[i].issubset(self.history_set[j]) and self.history_set[i] != self.history_set[j]:
                    redundancy_list.append(j)
        self.matrix=np.asarray([self.matrix[i] for i in range(len(self.matrix)) if i not in redundancy_list])
        self.history_set=[self.history_set[i] for i in range(len(self.history_set)) if i not in redundancy_list]
    
    def remove_constant_only_rows(self):
        r"""
        Remove the constant row and check naive feasibility of constant row.
        """
        redundancy_list=[]
        for i in range(len(self.matrix)):
            row = self.matrix[i]
            if np.all(row[:-1]==0) and row[-1]>0:
                print('The system is not feasible.')
                return
            if np.all(row[:-1]==0) and row[-1]<=0:
                redundancy_list.append(i)
        self.matrix=np.asarray([self.matrix[i] for i in range(len(self.matrix)) if i not in redundancy_list])
        self.history_set=[self.history_set[i] for i in range(len(self.history_set)) if i not in redundancy_list]

    def one_step_elimination(self, minimal_history_set = True):
        r"""
        Perform one step FM elimination to eliminate the first column.
        """
        equ=[]
        low=[]
        upp=[]
        new_matrix=[]
        new_history_set=[]
        for i in range(len(self.matrix)):
            row=self.matrix[i]
            if row[0]>0:
                upp.append(i)
            elif row[0]<0:
                low.append(i)
            else:
                equ.append(i)
        for e in equ:
            new_matrix.append(self.matrix[e][1:])
            new_history_set.append(self.history_set[e])
        for l in low:
            for u in upp:
                positive_row=self.matrix[u]
                negative_row=self.matrix[l]
                new_row=positive_row[1:]*(-negative_row[0])+negative_row[1:]*positive_row[0]
                new_h=self.history_set[u].union(self.history_set[l])
                new_matrix.append(new_row)
                new_history_set.append(new_h)
        new_FourierSystem = FourierSystem(field = self.field, matrix = np.asarray(new_matrix), history_set = new_history_set, binary_variables = self.binary_variables)
        if minimal_history_set:
            new_FourierSystem.remove_non_minimal_history_set()
        #with new_FourierSystem.field.off_the_record():
        #new_FourierSystem.remove_constant_only_rows()
        return new_FourierSystem

    def is_binary_only_row(self, row):
        r"""
        Return True if the row only contains binary variables.
        """
        with self.field.off_the_record():
            if row[:-self.binary_variables-1].any():
                return False
            else:
                return True

    def is_constant_only_row(self, row):
        r"""
        Return True if the row only contains the constant term.
        """
        with self.field.off_the_record():
            if np.all(row[:-1]==0):
                return True
            else:
                return False

    def is_binary_only_row_provable(self, row):
        r"""
        Return True if the row only contains binary variables and the row is redundant.
        """
        if self.is_binary_only_row(row):
            constant_term=row[-1]
            binary_terms=row[-self.binary_variables-1:-1]
            all_binary_vectors=generate_all_binary_vectors(self.binary_variables)
            for binary_vector in all_binary_vectors:
                test_term=constant_term+sum(binary_terms[i]*binary_vector[i] for i in range(self.binary_variables))
                assert test_term-1000000<=0
                self.field._bsa.tighten_upstairs_by_mccormick(max_iter=5)
                self.field.freeze()
                try:
                    assert test_term<=0
                except (ParametricRealFieldFrozenError, AssertionError) as e:
                    self.field.unfreeze()
                    return False
                self.field.unfreeze()
            return True
        else:
            return False


def find_new_inequalities_clipped_relu(K, param):
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
        S.one_step_elimination()
    FM_matrix=naive_simplify(S.matrix,binary_variables=3)
    FM_matrix=[normalize(row) for row in FM_matrix]

    res=[]
    for i1 in FM_matrix:
        if any([check_redundancy_one_to_one(i1,i2,3) for i2 in exponential_valid_inequalities]):
            continue
        else:
            res.append(i1)
    return res

def initialize_history_set(matrix):
    return [{i} for i in range(len(matrix))]

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

def FM_relu_1d(K):
    """
    One dimensional relu. (x0,x1,x,y,z,1).
    
    EXAMPLES::
    
        sage: import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *
        sage: from cutgeneratingfunctionology.spam.formulation import *
        sage: logging.disable(logging.INFO)
        sage: K.<L,U,W,b> = ParametricRealField([QQ(-2),QQ(2),QQ(1),QQ(1/2)], allow_coercion_to_float=False)
        sage: FS = FM_relu_1d(K)
        sage: FS.matrix
        array([[0, -1, 0, 0],
        [0, 0, 0, 0],
        [0, 0, (-L*W - b)~, (L*W + b)~],
        [0, 0, (-L + U)~, (L - U)~],
        [0, 0, 0, 0],
        [0, 1, (-U*W - b)~, 0],
        [(-W)~, 1, (-L*W - b)~, (L*W)~],
        [0, -1, (L*W + b)~, 0],
        [0, 0, (L - U)~, 0],
        [-1, 0, 0, L~],
        [(W^2)~, (-W)~, 0, (W*b)~],
        [W~, 0, (-U*W - b)~, b~],
        [W~, -1, (U*W + b)~, (-U*W)~],
        [1, 0, 0, (-U)~]], dtype=object)
    """
    L,U,W,b = K.gens()
    A=[]
    A.append([1,1,-1,0,0,0])
    A.append([-1,-1,1,0,0,0])
    A.append([W,0,0,0,-b,b])
    A.append([0,0,0,-1,0,0])
    A.append([0,W,0,-1,b,0])
    A.append([0,-W,0,1,-b,0])
    A.append([1,0,0,0,U,-U])
    A.append([-1,0,0,0,-L,L])
    A.append([0,1,0,0,-U,0])
    A.append([0,-1,0,0,L,0])

    A=np.asarray(A)
    FS=FourierSystem(K,A,initialize_history_set(A), binary_variables = 1)
    FS1=FS.one_step_elimination()
    FS2=FS1.one_step_elimination()
    return FS2

def FM_relu_2d(K):
    """
    Two dimensional relu. (x0_1,x1_1,x0_2,x1_2,x_1,x_2,y,z,1).
    """
    L1,U1,L2,U2,W1,W2,b = K.gens()
    #K.<L1,U1,L2,U2,W1,W2,b>=ParametricRealField([QQ(-2),QQ(2),QQ(-1),QQ(3),QQ(1),QQ(2),QQ(1/2)])
    A=[]
    A.append([1,1,0,0,-1,0,0,0,0])
    A.append([-1,-1,0,0,1,0,0,0,0])
    A.append([0,0,1,1,0,-1,0,0,0])
    A.append([0,0,-1,-1,0,1,0,0,0])
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
    FS=FourierSystem(K,A,initialize_history_set(A),binary_variables=1)
    FS1=FS.one_step_elimination()
    FS2=FS1.one_step_elimination()
    FS3=FS2.one_step_elimination()
    FS4=FS3.one_step_elimination()
    return FS4

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

def FM_clipped_relu_1d(K, **kwds):
    """
    One dimensional relu. (x1,x2,x3,x,y,z1,z2,z3,1).
    """
    L,U,W,b,C = K.gens()
    #K.<L,U,W,b,C>=ParametricRealField([QQ(-2),QQ(2),QQ(1),QQ(1/2),QQ(2)])
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

def FM_clipped_relu_2d(K, **kwds):
    """
    Two dimensional relu. (x1_1,x1_2,x2_1,x2_2,x3_1,x3_2,x1,x2,y,z1,z2,z3,1).
    """
    L1,U1,L2,U2,W1,W2,b,C = K.gens()
    #K.<L1,U1,L2,U2,W1,W2,b,C>=ParametricRealField([QQ(-2),QQ(2),QQ(-1),QQ(3),QQ(1),QQ(2),QQ(1/2),QQ(2)])
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

def verifying_one_dim_case(number_of_cases,backend='normaliz'):
    i=0
    while i<number_of_cases:
        l=np.random.randint(-100,-1)
        u=np.random.randint(1,100)
        w=np.random.randint(-100,100)
        b=np.random.randint(-10,10)
        c=np.random.randint(1,1000)
        m_plus,m_minus=affine_function_range(l,u,w,b)
        if m_minus>=0 or m_plus<=c:
            continue
        p_relaxed=clipped_relu_polyhedron_relaxed(l,u,w,b,c,backend=backend)
        p_exact=clipped_relu_polyhedron_exact(l,u,w,b,c,backend=backend)
        if not p_relaxed==p_exact:
            return p_relaxed,p_exact
        i+=1
    return True

def k_dim_test_example(dim,max_number_cases=10,backend='normaliz'):
    i=0
    while i<max_number_cases:
        L=np.random.randint(-1000,-1,size=dim)
        U=np.random.randint(1,1000,size=dim)
        W=np.random.randint(-1000,1000,size=dim)
        b=np.random.randint(-100,100)
        C=np.random.randint(1,100)
        m_plus,m_minus=affine_function_range(L,U,W,b)
        if m_minus>=0 or m_plus<=C:
            continue
        p_exact=clipped_relu_polyhedron_exact(L,U,W,b,C,backend=backend)
        p_strengthened=clipped_relu_polyhedron_strengthened(L,U,W,b,C,backend=backend)
        p_more_strengthened=clipped_relu_polyhedron_more_strengthened(L,U,W,b,C,backend=backend)
        p_relaxed=clipped_relu_polyhedron_relaxed(L,U,W,b,C,backend=backend)
        assert(is_subpolyhedron(p_exact,p_strengthened))
        assert(is_subpolyhedron(p_strengthened,p_more_strengthened))
        assert(is_subpolyhedron(p_more_strengthened,p_relaxed))
        if not p_more_strengthened==p_exact:
            return L,U,W,b,C
        i+=1
    return True

def is_subpolyhedron(P,Q):
    """
    Return whether the bounded polyhedron P is contained in polyhedron Q.
    """
    return all(Q.contains(v) for v in P.Vrepresentation())

def affine_function_range(L,U,W,b):
    m_plus=b
    m_minus=b
    for i in range(len(W)):
        w=W[i]
        if w>=0:
            m_plus+=w*U[i]
            m_minus+=w*L[i]
        else:
            m_plus+=w*L[i]
            m_minus+=w*U[i]
    return m_plus,m_minus

def clipped_relu_polyhedron_relaxed(L,U,W,b,C,backend='normaliz'):
    m_plus,m_minus=affine_function_range(L,U,W,b)
    ieqs=generate_clipped_relu_ieqs(L,U,W,b,C,m_plus,m_minus)
    return Polyhedron(ieqs=ieqs,backend=backend)

def clipped_relu_polyhedron_strengthened(L,U,W,b,C,backend='normaliz'):
    m_plus,m_minus=affine_function_range(L,U,W,b)
    ieqs=generate_clipped_relu_ieqs_strengthened(L,U,W,b,C,m_plus,m_minus)
    return Polyhedron(ieqs=ieqs,backend=backend)

def clipped_relu_polyhedron_more_strengthened(L,U,W,b,C,backend='normaliz'):
    m_plus,m_minus=affine_function_range(L,U,W,b)
    ieqs=generate_clipped_relu_ieqs_more_strengthened(L,U,W,b,C,m_plus,m_minus)
    return Polyhedron(ieqs=ieqs,backend=backend)

def clipped_relu_polyhedron_exact(L,U,W,b,C,backend='normaliz'):
    p_left=clipped_relu_polyhedron_left(L,U,W,b,C,backend=backend)
    p_mid=clipped_relu_polyhedron_mid(L,U,W,b,C,backend=backend)
    p_right=clipped_relu_polyhedron_right(L,U,W,b,C,backend=backend)
    vertices=p_left.Vrepresentation()+p_mid.Vrepresentation()+p_right.Vrepresentation()
    return Polyhedron(vertices=vertices,backend=backend)

def clipped_relu_polyhedron_left(L,U,W,b,C,backend='normaliz'):
    m_plus,m_minus=affine_function_range(L,U,W,b)
    ieqs=generate_clipped_relu_ieqs_left(L,U,W,b,C,m_plus,m_minus)
    return Polyhedron(ieqs=ieqs,backend=backend)

def clipped_relu_polyhedron_mid(L,U,W,b,C,backend='normaliz'):
    m_plus,m_minus=affine_function_range(L,U,W,b)
    ieqs=generate_clipped_relu_ieqs_mid(L,U,W,b,C,m_plus,m_minus)
    return Polyhedron(ieqs=ieqs,backend=backend)

def clipped_relu_polyhedron_right(L,U,W,b,C,backend='normaliz'):
    m_plus,m_minus=affine_function_range(L,U,W,b)
    ieqs=generate_clipped_relu_ieqs_right(L,U,W,b,C,m_plus,m_minus)
    return Polyhedron(ieqs=ieqs,backend=backend)

def generate_clipped_relu_ieqs(L,U,W,b,C,m_plus,m_minus):
    """
    Generate coefficients (constant,x,y,z1,z2,z3) of inequalities of clipped relu formulation.
    """
    res=[]
    zero_coefs=tuple([0]*len(W))
    W_coefs=tuple(W)
    nW_coefs=tuple([-w for w in W])
    # c*z3<=y<=c*(z2+z3)
    res.append((0,)+zero_coefs+(1,0,0,-C))
    res.append((0,)+zero_coefs+(-1,0,C,C))
    # y<=w*x+b-m_minus*z1
    res.append((b,)+W_coefs+(-1,-m_minus,0,0))
    # y>=w*x+b+(C-m_plus)*z3
    res.append((-b,)+nW_coefs+(1,0,0,m_plus-C))
    # z1+z2+z3=1
    res.append((-1,)+zero_coefs+(0,1,1,1))
    res.append((1,)+zero_coefs+(0,-1,-1,-1))
    # L<=x<=U
    for i in range(len(W)):
        ei=[0]*len(W)
        nei=[0]*len(W)
        ei[i]=1
        nei[i]=-1
        res.append((-L[i],)+tuple(ei)+(0,0,0,0))
        res.append((U[i],)+tuple(nei)+(0,0,0,0))
    # 0<=zi<=1
    res.append((0,)+zero_coefs+(0,1,0,0))
    res.append((0,)+zero_coefs+(0,0,1,0))
    res.append((0,)+zero_coefs+(0,0,0,1))
    res.append((1,)+zero_coefs+(0,-1,0,0))
    res.append((1,)+zero_coefs+(0,0,-1,0))
    res.append((1,)+zero_coefs+(0,0,0,-1))
    return res

def generate_clipped_relu_ieqs_left(L,U,W,b,C,m_plus,m_minus):
    res=generate_clipped_relu_ieqs(L,U,W,b,C,m_plus,m_minus)
    zero_coefs=tuple([0]*len(W))
    # z1=1,z2=0,z3=0
    res.append((-1,)+zero_coefs+(0,1,0,0))
    res.append((0,)+zero_coefs+(0,0,-1,0))
    res.append((0,)+zero_coefs+(0,0,0,-1))
    return res

def generate_clipped_relu_ieqs_mid(L,U,W,b,C,m_plus,m_minus):
    res=generate_clipped_relu_ieqs(L,U,W,b,C,m_plus,m_minus)
    zero_coefs=tuple([0]*len(W))
    # z1=0,z2=1,z3=0
    res.append((0,)+zero_coefs+(0,-1,0,0))
    res.append((-1,)+zero_coefs+(0,0,1,0))
    res.append((0,)+zero_coefs+(0,0,0,-1))
    return res

def generate_clipped_relu_ieqs_right(L,U,W,b,C,m_plus,m_minus):
    res=generate_clipped_relu_ieqs(L,U,W,b,C,m_plus,m_minus)
    zero_coefs=tuple([0]*len(W))
    # z1=0,z2=0,z3=1
    res.append((0,)+zero_coefs+(0,-1,0,0))
    res.append((0,)+zero_coefs+(0,0,-1,0))
    res.append((-1,)+zero_coefs+(0,0,0,1))
    return res

def h_L(L,U,W,I):
    res=0
    for i in range(len(I)):
        if I[i]==0:
            continue
        if W[i]>=0:
            res+=W[i]*L[i]
        else:
            res+=W[i]*U[i]
    return res

def h_U(L,U,W,I):
    res=0
    for i in range(len(I)):
        if I[i]==0:
            continue
        if W[i]>=0:
            res+=W[i]*U[i]
        else:
            res+=W[i]*L[i]
    return res

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

def generate_clipped_relu_ieqs_strengthened(L,U,W,b,C,m_plus,m_minus):
    res=generate_clipped_relu_ieqs(L,U,W,b,C,m_plus,m_minus)
    for I in generate_all_binary_vectors(len(W)):
        I_bar=complement_index(I)
        x_coefs=[I[i]*W[i] for i in range(len(W))]
        nx_coefs=[-t for t in x_coefs]
        # strengthening inequality (12)
        res.append((0,)+tuple(x_coefs)+(-1,-h_L(L,U,W,I),b+h_U(L,U,W,I_bar),b+h_U(L,U,W,I_bar)))
        # strengthening inequality (13)
        res.append((0,)+tuple(nx_coefs)+(1,-b-h_L(L,U,W,I_bar),-b-h_L(L,U,W,I_bar),h_U(L,U,W,I)-C))
    return res

def generate_clipped_relu_ieqs_more_strengthened(L,U,W,b,C,m_plus,m_minus):
    res=generate_clipped_relu_ieqs_strengthened(L,U,W,b,C,m_plus,m_minus)
    for I in generate_all_binary_vectors(len(W)):
        I_bar=complement_index(I)
        x_coefs=[I[i]*W[i] for i in range(len(W))]
        nx_coefs=[-t for t in x_coefs]
        res.append((0,)+tuple(x_coefs)+(0,-h_L(L,U,W,I),-h_L(L,U,W,I),b-C+h_U(L,U,W,I_bar)))
        res.append((0,)+tuple(nx_coefs)+(0,-b-h_L(L,U,W,I_bar),h_U(L,U,W,I),h_U(L,U,W,I)))
    return res
