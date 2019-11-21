r"""
FM elimination examples: relu and clipped relu using BasicSemialgebraicSet_polyhedral_linear_system
"""

from __future__ import division, print_function, absolute_import

from cutgeneratingfunctionology.spam.basic_semialgebraic_linear_system import BasicSemialgebraicSet_polyhedral_linear_system

import itertools

def FM_relu_1d(base_ring, poly_ring):
    """
    One dimensional relu y=max(0, W*x+b) with x in [L,U].
    Original variables (x0,x1,x,y,z) in extended formulation.
    We have x0+x1=x and the goal is to eliminate x0,x1 and get a non-extended formulation.
    
    EXAMPLE::
    
        sage: import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *
        sage: from cutgeneratingfunctionology.spam.examples.relu import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: K.<L,U,W,b> = ParametricRealField([QQ(-2),QQ(2),QQ(2),QQ(1/2)])
        sage: Q.<x0,x1,x,y,z> = K[]
        sage: bsa = FM_relu_1d(K,Q)
        sage: with K.off_the_record():
        ....:     bsa.le_poly() # The fourth and last inequalities are redundant.
        {-z,
        z - 1,
        -y,
        -y + ((L*W + b)~)*z,
        y + ((-U*W - b)~)*z,
        ((-W)~)*x + y + ((-L*W - b)~)*z + (L*W)~,
        W~*x - y + b~,
        W~*x - y + ((U*W + b)~)*z + (-U*W)~}
    """
    x0,x1,x,y,z = iter(poly_ring.gens())
    L,U,W,b = iter(base_ring.gens())
    # assumptions
    assert L<0
    assert U>0
    assert W>0
    assert W*L+b<0
    assert W*U+b>0
    assert b>0
    
    base_ring.freeze()
    
    le = [-y, W*x0-b*z+b, x0+U*z-U, -x0-L*z+L, x1-U*z, -x1+L*z, -z, z-1]
    eq = [x0+x1-x, W*x1-y+b*z]
    bsa = BasicSemialgebraicSet_polyhedral_linear_system(poly_ring=poly_ring, le=le, eq=eq)
    bsa_eliminated = bsa.coordinate_projection([x0,x1])
    # normalize base on y, the second last generator.
    normalized_le = set([normalize(le, -2) for le in bsa_eliminated.le_poly()])
    return BasicSemialgebraicSet_polyhedral_linear_system(poly_ring=bsa_eliminated.poly_ring(), le=normalized_le)

def FM_relu_2d(base_ring,poly_ring):
    """
    Two dimensional relu y = max(0, W1*x_1 + W2*x_2 + b) with x_1 in [L1,U1] and x_2 in [L2,U2].
    Original variables (x0_1,x1_1,x0_2,x1_2,x_1,x_2,y,z) in extended formulation.
    We have x0+x1=x, and x0,x1,x are all two dimensional vectors.
    The goal is to eliminate x0,x1 and get a non-extended formulation.
    
    EXAMPLE::
    
        sage: import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *
        sage: from cutgeneratingfunctionology.spam.examples.relu import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: K.<L1,U1,L2,U2,W1,W2,b>=ParametricRealField([QQ(-2),QQ(2),QQ(-3),QQ(3),QQ(5),QQ(2),QQ(1/2)])
        sage: Q.<x0_1,x1_1,x0_2,x1_2,x_1,x_2,y,z> = K[]
        sage: bsa = FM_relu_2d(K,Q)
        Traceback (most recent call last):
        ...
        ParametricRealFieldFrozenError: Cannot prove that constraint is implied: -U1*W1 - b < 0
        sage: K._bsa.tighten_upstairs_by_mccormick()  # Actually the inequality -U1*W1 - b < 0 is provable, so we add mccormick inequality first.
        sage: bsa = FM_relu_2d(K,Q) # This time the inequality L1*W1 + b < 0 is not provable.
        Traceback (most recent call last):
        ...
        ParametricRealFieldFrozenError: Cannot prove that constraint is implied: L1*W1 + b < 0
            
    """
    x0_1,x1_1,x0_2,x1_2,x_1,x_2,y,z = iter(poly_ring.gens())
    L1,U1,L2,U2,W1,W2,b = iter(base_ring.gens())
    # assumptions
    assert L1<0
    assert L2<0
    assert U1>0
    assert U2>0
    assert W1>0
    assert W2>0
    assert W1*L1+W2*L2+b<0
    assert W1*U1+W2*U2+b>0
    assert b>0
    
    base_ring.freeze()
    
    le = [-y, W1*x0_1+W2*x0_2-b*z+b, x0_1+U1*z-U1, -x0_1-L1*z+L1, x0_2+U2*z-U2, -x0_2-L2*z+L2,\
          x1_1-U1*z, -x1_1+L1*z, x1_2-U2*z, -x1_2+L2*z, -z, z-1]
    eq = [x0_1+x1_1-x_1, x0_2+x1_2-x_2, W1*x1_1+W2*x1_2-y+b*z]
    bsa = BasicSemialgebraicSet_polyhedral_linear_system(poly_ring=poly_ring, le=le, eq=eq)
    bsa_eliminated = bsa.coordinate_projection([x0_1,x1_1,x0_2,x1_2])
    # normalize base on y, the second last generator.
    normalized_le = set([normalize(le, -2) for le in bsa_eliminated.le_poly()])
    return BasicSemialgebraicSet_polyhedral_linear_system(poly_ring=bsa_eliminated.poly_ring(), le=normalized_le)

def normalize(polynomial, index):
    poly_ring = polynomial.parent()
    base_ring = poly_ring.base_ring()
    v = poly_ring.gens()[index]
    c = base_ring(polynomial.monomial_coefficient(v))
    if c>0:
        return polynomial/c
    elif c<0:
        return polynomial/-c
    else:
        return polynomial

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

