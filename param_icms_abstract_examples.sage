def test_posdef(a, b):
    return is_positive_definite(matrix([[a, b], [b, 1/4]]))

def is_positive_definite(self):
    # copyed from "sage/matrix/matrix2.pyx":12383, modified L12579.
    imaginary = True
    a = self.base_ring().an_element()
    try:
        a.conjugate()
    except AttributeError:
        imaginary = False
    if imaginary:
        if not self.is_hermitian():
            return False
    else:
        if not self.is_symmetric():
            return False
    try:
        if imaginary:
            _, d = self._indefinite_factorization('hermitian', check=False)
        else:
            _, d = self._indefinite_factorization('symmetric', check=False)
    except ValueError as e:
        # a zero singular leading principal submatrix is one
        # indicator that the matrix is not positive definite
        if str(e).find('leading principal submatrix is singular') != -1:
            return False
        else:
            raise ValueError(e)
    # Now have diagonal entries (hopefully real) and so can
    # test with a generator (which will short-circuit)
    # positive definite iff all entries of d are positive
    try:
        posdef = all( x > 0 for x in d )  # change RR(x)>0 to x>0, x is always real.
    except TypeError:
        universe = Sequence(d).universe()
        msg = "cannot convert computations from {0} into real numbers"
        raise TypeError(msg.format(universe))
    return posdef
    

def find_region_type(field, result):
    # Return these label because that's what the plotting code expects.
    # TODO: Make mapping customizable
    if result:
        return 'is_extreme'
    else:
        return 'not_extreme'
    
complex = SemialgebraicComplex(test_posdef, ['a', 'b'], find_region_type=find_region_type, default_var_bound=(-2,2))

complex.shoot_random_points(10)
# or sage: complex.bfs_completion(flip_ineq_step=1/100)
plot(complex)


# TODO:
# customize find_region_type mapping, more examples with bfs,
# flip_ineq_step should always >0, but choose shooting, mathematic, gradient, linear...
# it should be user-friendly; not specific to gomory-johnson; have cool examples; including examples where the output are different formulas (rather than just truth values) in different cells; ...
# For example a formula for the square of the smallest eigenvalue.
# example, number of positive eigenvalues?
