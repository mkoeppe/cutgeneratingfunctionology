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
        sage: bsa.le_poly() # The second and fifth inequalities are redundant.
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

