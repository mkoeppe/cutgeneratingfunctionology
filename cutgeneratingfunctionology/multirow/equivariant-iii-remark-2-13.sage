## Example from Equivariant III, Remark 2.13

def remark_2_13_vertices():
    """
    The additivity domain from Remark 2.13.
    
    The following doctests verify some claims of this remark.

    TESTS::

        sage: vertices = remark_2_13_vertices()
        sage: F = Polyhedron(vertices)
        sage: len(vertices) == len(F.vertices())
        True
        sage: p1F, p2F, p3F = projections(F)
        sage: sorted(p1F.vertices_list())
        [[0, 3], [0, 4], [3, 4]]
        sage: sorted(p2F.vertices_list())
        [[4, 1], [5, 0], [5, 2]]
        sage: sorted(p3F.vertices_list())
        [[5, 5], [5, 6], [8, 4]]
        sage: F == FIJK(p1F, p2F, p3F)
        True

    EXAMPLES::

        Reproduce a part of Figure 4.

        sage: plot(p1F,polygon='red') + plot(p2F,polygon='blue') + plot(p3F,polygon='green') #not tested

    """
    vertices = [ vector (v)
                 for v in [ [ 0, 4, 5, 1 ], [ 0, 3, 5, 2 ], [ 2, 11/3, 4, 1],
                            [ 1, 4, 4, 1 ], [ 8/3, 35/9, 4, 1 ], [ 5/2, 4, 4, 1 ],
                            [ 1, 10/3, 5, 2 ], [ 0, 4, 5, 2 ], [ 3, 4, 5, 0 ] ] ]
    return vertices

def projections(F):
    """
    Return the projections p1, p2, p3 of a face F.

    Currently only implemented for the two-dimensional case.
    """
    vertices = F.vertices()
    p1 = matrix(QQ, [[1, 0, 0, 0], [0, 1, 0, 0]])
    p1F = Polyhedron([ p1 * v.vector() for v in vertices ])
    p2 = matrix(QQ, [[0, 0, 1, 0], [0, 0, 0, 1]])
    p2F = Polyhedron([ p2 * v.vector() for v in vertices ])
    p3 = p1 + p2
    p3F = Polyhedron([ p3 * v.vector() for v in vertices ])
    return p1F, p2F, p3F

def FIJK(I, J, K):
    """
    The face (polyhedron) F(I, J, K), where I, J, K are polyhedra.
    """
    zero = [0] * I.ambient_dim()
    ieqs = []
    eqns = []
    ieqs.extend( [ [i[0]] + list(i[1:]) + zero
                   for i in I.inequality_generator() ] )
    eqns.extend( [ [i[0]] + list(i[1:]) + zero
                   for i in I.equation_generator() ] )
    ieqs.extend( [ [i[0]] + zero + list(i[1:])
                   for i in J.inequality_generator() ] )
    eqns.extend( [ [i[0]] + zero + list(i[1:])
                   for i in J.equation_generator() ] )
    ieqs.extend( [ [i[0]] + list(i[1:]) + list(i[1:])
                   for i in K.inequality_generator() ] )
    eqns.extend( [ [i[0]] + list(i[1:]) + list(i[1:])
                   for i in K.equation_generator() ] )
    return Polyhedron(ieqs=ieqs, eqns=eqns)
