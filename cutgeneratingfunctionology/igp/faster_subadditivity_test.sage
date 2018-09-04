# Make sure current directory is in path.
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *
from itertools import chain

import itertools
import numpy as np


def find_best_slope_intercept(X,Y,lower_bound=True,solver='GLPK'):
    """
    Find the slope anf intercept of the affine lower/upper estimator.
    """
    if lower_bound:
        p = MixedIntegerLinearProgram(maximization=True, solver=solver)
        v = p.new_variable()
        m, b= v['m'], v['b']
        p.set_objective(m*sum(X)+b*len(X))
        for i in range(len(X)):
            p.add_constraint(m*X[i]+b<=Y[i])
    else:
        p = MixedIntegerLinearProgram(maximization=False, solver=solver)
        v = p.new_variable()
        m, b= v['m'], v['b']
        p.set_objective(m*sum(X)+b*len(X))
        for i in range(len(X)):
            p.add_constraint(m*X[i]+b>=Y[i])
    p.solve()
    return [p.get_values(m),p.get_values(b)]

def find_affine_bounds(fn,I,lower_bound=True,extended=False):
    """
    Find the affine lower/upper estimator of the function fn over the closed interval I.
    """
    if not extended:
        bkpts=fn.end_points()
    else:
        bkpts=fn.end_points()[:-1]+[1+bkpt for bkpt in fn.end_points()]
    I_index=find_possible_branching_bkpts_index(bkpts,I)
    X,Y=[I[0],I[1]],[fn(fractional(I[0])),fn(fractional(I[1]))]
    if I_index:
        for i in range(I_index[0], I_index[1]+1):
            X.append(bkpts[i])
            Y.append(fn(fractional(bkpts[i])))
    slope,intercept=find_best_slope_intercept(X,Y,lower_bound=lower_bound)
    return slope, intercept

def find_possible_branching_bkpts_index(bkpts,I):
    """
    Return [i,j] such that I[0]<bkpts[k]<I[1] for any i<=k<=j. Return None if empty.
    """
    i=bisect_left(bkpts, I[0])
    j=bisect_left(bkpts, I[1])-1
    if I[0]==bkpts[i]:
        i=i+1
    if i>j:
        return None
    else:
        return [i,j]

def fn_constant_bounds(fn,I,lower_bound=True,extended=False):
    """
    Find the constant lower/upper bound of the function fn over the closed interval I.
    """
    if not extended:
        bkpts=fn.end_points()
    else:
        bkpts=fn.end_points()[:-1]+[1+bkpt for bkpt in fn.end_points()]
    I_index=find_possible_branching_bkpts_index(bkpts,I)
    if lower_bound:
        if not I_index:
            return min(fn(fractional(I[0])), fn(fractional(I[1])))
        else:
            return min(fn(fractional(I[0])), fn(fractional(I[1])), min(fn(fractional(bkpts[i])) for i in range(I_index[0], I_index[1]+1)))
    else:
        if not I_index:
            return max(fn(fractional(I[0])), fn(fractional(I[1])))
        else:
            return max(fn(fractional(I[0])), fn(fractional(I[1])), min(fn(fractional(bkpts[i])) for i in range(I_index[0], I_index[1]+1)))

def find_branching_bkpt(fn,I_index,extended=False):
    bkpts1=fn.end_points()
    bkpts2=fn.end_points()[:-1]+[1+bkpt for bkpt in fn.end_points()]
    if extended:
        return bkpts2[floor((I_index[0]+I_index[1])/2)]
    else:
        return bkpts1[floor((I_index[0]+I_index[1])/2)]

def find_branching_direction(fn,I,J,K,I_index,J_index,K_index):
    lenI=I[1]-I[0]
    lenJ=J[1]-J[0]
    lenK=K[1]-K[0]
    candidates=[]
    if I_index:
        candidates.append(['I',lenI])
    if J_index:
        candidates.append(['J',lenJ])
    if K_index:
        candidates.append(['K',lenK])
    if not candidates:
        return None
    else:
        ind=np.argmax([candidate[1] for candidate in candidates])
        return candidates[ind][0]

def delta_pi_min(fn,I,J,K,approximation='constant'):
    """
    Return the minimum of the function delta fn in the region whose projections are I,J,K.
    """
    bkpts1=fn.end_points()
    bkpts2=fn.end_points()[:-1]+[1+bkpt for bkpt in fn.end_points()]
    vertices=verts(I,J,K)
    if len(vertices)==0:
        return 100
    elif len(vertices)==1:
        return delta_pi(fn,vertices[0][0],vertices[0][1])
    new_I,new_J,new_K=projections(vertices)
    I_index=find_possible_branching_bkpts_index(bkpts1,new_I)
    J_index=find_possible_branching_bkpts_index(bkpts1,new_J)
    K_index=find_possible_branching_bkpts_index(bkpts2,new_K)
    if approximation=='constant':
        alpha_I=fn_constant_bounds(fn,new_I,lower_bound=True,extended=False)
        alpha_J=fn_constant_bounds(fn,new_J,lower_bound=True,extended=False)
        beta_K=fn_constant_bounds(fn,new_K,lower_bound=False,extended=True)
        if alpha_I+alpha_J-beta_K>0:
            return alpha_I+alpha_J-beta_K
    elif approximation=='affine':
        m_I, b_I=find_affine_bounds(fn,new_I,lower_bound=True,extended=False)
        m_J, b_J=find_affine_bounds(fn,new_J,lower_bound=True,extended=False)
        m_K, b_K=find_affine_bounds(fn,new_K,lower_bound=False,extended=True)
        delta_lower_bound=min((m_I*vertex[0]+b_I)+(m_J*vertex[1]+b_J)-(m_K*fractional(vertex[0]+vertex[1])+b_K) for vertex in vertices)
        if delta_lower_bound>0:
            return delta_lower_bound
    upper_bound=min(delta_pi(fn,v[0],v[1]) for v in vertices)
    if upper_bound<0:
        return upper_bound
    branch_direction=find_branching_direction(fn,new_I,new_J,new_K,I_index,J_index,K_index)
    if not branch_direction:
        return upper_bound
    elif branch_direction=='I':
        bkpt=find_branching_bkpt(fn,I_index,extended=False)
        return min(delta_pi_min(fn,[I[0],bkpt],J,K,approximation=approximation),delta_pi_min(fn,[bkpt,I[1]],J,K,approximation=approximation))
    elif branch_direction=='J':
        bkpt=find_branching_bkpt(fn,J_index,extended=False)
        return min(delta_pi_min(fn,I,[J[0],bkpt],K,approximation=approximation),delta_pi_min(fn,I,[bkpt,J[1]],K,approximation=approximation))
    elif branch_direction=='K':
        bkpt=find_branching_bkpt(fn,K_index,extended=True)
        return min(delta_pi_min(fn,I,J,[K[0],bkpt],approximation=approximation),delta_pi_min(fn,I,J,[bkpt,K[1]],approximation=approximation))
    
