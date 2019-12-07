r"""
Utilities for programming with big cells -- implementation.
"""

from __future__ import division, print_function, absolute_import

trivial_parametric_real_field = None

# should export from a different module
def _common_parametric_real_field(iterable, key=None):
    """
    EXAMPLES::

        sage: from cutgeneratingfunctionology.spam.big_cells_impl import _common_parametric_real_field
        sage: from cutgeneratingfunctionology.igp import ParametricRealField
        sage: def generate_incompatible_real_field_elements():
        ....:     yield ParametricRealField([1], names=['f']).gen()
        ....:     yield ParametricRealField([1], names=['x']).gen()
        sage: _common_parametric_real_field(generate_incompatible_real_field_elements()) #should give error
        Traceback (most recent call last):
        ...
        RuntimeError: ...
    """
    from cutgeneratingfunctionology.igp import ParametricRealField
    if key is None:
        key = lambda i: i
    field = None
    for x in iterable:
        try:
            parent = key(x).parent()
        except AttributeError:
            parent = None
        if isinstance(parent, ParametricRealField):
            if field is None:
                field = parent
            elif field is not parent:
                raise TypeError("arguments from different parametric fields")
    if field:
        return field
    else:
        global trivial_parametric_real_field
        if trivial_parametric_real_field is None:
            trivial_parametric_real_field = ParametricRealField(names=(), big_cells=False)
        return trivial_parametric_real_field

def big_cells_min(iterable, *args, **kwds):
    """
    Compute the minimum of the values of the function ``key``
    on the elements of ``iterable``.

    EXAMPLES::

        sage: import logging
        sage: logging.disable(logging.WARNING)
        sage: from cutgeneratingfunctionology import igp
        sage: from cutgeneratingfunctionology.igp import ParametricRealField, result_symbolic_expression, SemialgebraicComplex
        sage: from cutgeneratingfunctionology.spam import big_cells
        sage: K.<a,b> = ParametricRealField([2, 1], big_cells=True, allow_refinement=True)
        sage: big_cells.min([a, b, 2])
        b~
        sage: sorted(K.get_le_factor(), key=str)
        [-a + b, b - 2]

     If ``allow_refinement`` is False, the same example yields an error.

        sage: K.<a,b> = ParametricRealField([2, 1], big_cells=True, allow_refinement=False)
        sage: big_cells.min([a, b, 2])
        Traceback (most recent call last):
        ...
        ParametricRealFieldRefinementError: ...
        sage: assert a <= 2
        sage: assert b <= 2
        sage: big_cells.min([a, b, 2])
        Traceback (most recent call last):
        ...
        ParametricRealFieldRefinementError: ...
        sage: assert b <= a
        sage: big_cells.min([a, b, 2])
        b~

    Compared to ``min``, this function aims to create bigger
    parameter cells::

        sage: igp.big_cells_default = True
        sage: igp.allow_refinement_default = True
        sage: def f(a=2, b=1): return min([a, b, 3])
        sage: def big_cells_f(a=2, b=1): return big_cells.min([a, b, 3])
        sage: complex = SemialgebraicComplex(f, ['a', 'b'],
        ....:                                find_region_type=result_symbolic_expression,
        ....:                                default_var_bound=(0,4))
        sage: complex.shoot_random_points(10)
        sage: len(complex.components)
        4
        sage: big_cells_complex = SemialgebraicComplex(big_cells_f, ['a', 'b'],
        ....:                                          find_region_type=result_symbolic_expression,
        ....:                                          default_var_bound=(0,4))
        sage: big_cells_complex.shoot_random_points(10)
        sage: len(big_cells_complex.components)
        3
        sage: graphics_array([complex.plot(), big_cells_complex.plot()])    # not tested

    TESTS:

        The alternative interface of the built-in function min is suppported::

            sage: big_cells.min(5, 2, 3)
            2

        An error is raised when there are no args::

            sage: big_cells.min()
            Traceback (most recent call last):
            ...
            TypeError: ...argument...

        Default::

            sage: big_cells.min([], default=77)
            77

    """
    if args:
        iterable = [iterable] + list(args)
    key = kwds.pop('key', None)
    field = kwds.pop('field', None)
    if key is None:
        key = lambda i: i
    iv_list = [ (i, key(i)) for i in iterable ]
    if not iv_list:
        if 'default' in kwds:
            return kwds['default']
        raise TypeError("min expected 1 arguments, got 0")
    if field is None:
        field = _common_parametric_real_field(iv_list, key=lambda iv: iv[1])
    from cutgeneratingfunctionology.igp import ParametricRealField
    if not isinstance(field, ParametricRealField):
        field = trivial_parametric_real_field
    if not field._big_cells:
        min_i, min_v = min(iv_list, key=lambda iv: iv[1])
        return min_i
    with field.off_the_record():
        min_i, min_v = min(iv_list, key=lambda iv: iv[1])
    if field._allow_refinement:
        assert all(min_v <= iv[1] for iv in iv_list)   # records
    else:
        from cutgeneratingfunctionology.igp import ParametricRealFieldFrozenError, ParametricRealFieldRefinementError
        try:
            with field.frozen():
                assert all(min_v <= iv[1] for iv in iv_list)
        except ParametricRealFieldFrozenError:
            raise ParametricRealFieldRefinementError("big_cells_min")
        # cannot have AssertionError because we compare to min_v 
    return min_i

def is_min_le(iterable, value, key=None, field=None):
    """
    Decide whether the minimum of the values of the function ``key``
    on the elements of ``iterable`` is <= ``value``.

    EXAMPLES::

        sage: import logging
        sage: logging.disable(logging.WARNING)
        sage: from cutgeneratingfunctionology import igp
        sage: from cutgeneratingfunctionology.igp import ParametricRealField, result_symbolic_expression, SemialgebraicComplex
        sage: from cutgeneratingfunctionology.spam import big_cells

    Trivial case::

        sage: K.<a,b> = ParametricRealField([2, 1], big_cells=True, allow_refinement=True)
        sage: big_cells.is_min_le([], 3)
        True

    A non-convex region::

        sage: K.<a,b> = ParametricRealField([2, 1], big_cells=True, allow_refinement=True)
        sage: big_cells.is_min_le([a, b], 3)
        True
        sage: K.get_le_factor(), K.get_lt_factor()
        ({b - 3}, set())

    This is (necessarily) a refinement of the actual non-convex
    region where min(a, b) <= 3, but a bigger cell than what is
    obtained by the following::

        sage: K.<a,b> = ParametricRealField([2, 1], big_cells=True, allow_refinement=True)
        sage: min([a, b]) <= 3
        True
        sage: K.get_le_factor(), K.get_lt_factor()
        ({b - 3}, {-a + b})

    If ``allow_refinement`` is False, the same example yields an error.

        sage: K.<a,b> = ParametricRealField([2, 1], big_cells=True, allow_refinement=False)
        sage: big_cells.is_min_le([a, b], 3)
        Traceback (most recent call last):
        ...
        ParametricRealFieldRefinementError: ...

    The same example with the addtional condition a <= 3.

        sage: K.<a,b> = ParametricRealField([2, 1], big_cells=True, allow_refinement=False)
        sage: assert a <= 3
        sage: big_cells.is_min_le([a, b], 3)
        True

    Another example with ``allow_refinement`` False::

        sage: K.<a,b> = ParametricRealField([2, 1], big_cells=True, allow_refinement=False)
        sage: assert b >= 0
        sage: big_cells.is_min_le([2*b, b], 3)
        True


    Example where the big cell is the entire space. 
    This requires a ParametricRealField set up with ``mutable_values=True``::

        sage: K.<a,b> = ParametricRealField([4, 1], big_cells=True, mutable_values=True, allow_refinement=False)
        sage: big_cells.is_min_le([a+2, 2-a, b+2, 2-b], 2)
        True

    Previous bug examples. Now fixed if ParametricRealField is set up with ``mutable_values=True``::

        sage: K.<a,b> = ParametricRealField([4, 1], big_cells=True, mutable_values=True, allow_refinement=False)
        sage: big_cells.is_min_le([3/4*a, 1/4*a], 2)
        True

        sage: K.<a,b> = ParametricRealField([4, 1], big_cells=True, mutable_values=True, allow_refinement=False)
        sage: big_cells.is_min_le([a, b], 2) # not convex
        Traceback (most recent call last):
        ...
        ParametricRealFieldRefinementError: ...

    But we should avoid ``mutable_values=True`` because it is extremely slow::

        sage: K.<a,b> = ParametricRealField([4, 1], big_cells=True, mutable_values=True, allow_refinement=False)
        sage: big_cells.is_min_le([a] * 10000, 5)      # not tested
        True

    Test that many repeated values do not require quadratic running time::

        sage: K.<a,b> = ParametricRealField([4, 1], big_cells=True, allow_refinement=False)
        sage: big_cells.is_min_le([a] * 10000, 5)      # long time - 40s
        True

    In fact, the big cells form a cover (arrangement), not a complex;
    there is a full-dimensional intersection::

        sage: igp.big_cells_default = False
        sage: igp.allow_refinement_default = True
        sage: def f(a=2, b=1): return min([a, b]) <= 3
        sage: complex = SemialgebraicComplex(f, ['a', 'b'],
        ....:                                find_region_type=result_symbolic_expression,
        ....:                                default_var_bound=(0,4))
        sage: complex.shoot_random_points(1)
        sage: complex.bfs_completion(wall_crossing_method='heuristic', goto_lower_dim=True)
        sage: len(complex.components)
        9
        sage: igp.big_cells_default = True
        sage: igp.allow_refinement_default = True
        sage: def big_cells_f(a=2, b=1): return big_cells.is_min_le([a, b], 3)
        sage: big_cells_complex = SemialgebraicComplex(big_cells_f, ['a', 'b'],
        ....:                                          find_region_type=result_symbolic_expression,
        ....:                                          default_var_bound=(0,4))
        sage: big_cells_complex.shoot_random_points(10)
        sage: len(big_cells_complex.components)
        3
        sage: graphics_array([complex.plot(), big_cells_complex.plot()])    # not tested

    """
    if key is None:
        key = lambda i: i
    iv_list = [ (i, key(i)) for i in iterable ]
    # Fast path for trivial cases
    if not iv_list:
        return True
    if len(iv_list) == 1:
        return iv_list[1] <= value    # records
    # Determine whether it's True or False.
    if field is None:
        field = _common_parametric_real_field(iv_list, key=lambda iv: iv[1])
    if not field._big_cells:
        min_i, min_v = min(iv_list, key=lambda iv: iv[1])
        return min_v <= value
    with field.off_the_record():
        min_i, min_v = min(iv_list, key=lambda iv: iv[1])
        is_le = min_v <= value
    # is_min_le is a comparison, not a computation.
    # when allow_refinement=False, comparisons are allowed to record;
    # this one should record the convex region; but if the region was non-convex, it would raise an error.
    if not is_le:
        # The "False" case is always convex and can just be recorded.
        assert all(iv[1] > value for iv in iv_list)   # records
    else:
        # In the "True" case, we try hard to keep a large region.
        if field._allow_refinement:
            assert min_v <= value    # records
        else:
            from cutgeneratingfunctionology.igp import ParametricRealFieldFrozenError, ParametricRealFieldRefinementError, ParametricRealFieldInconsistencyError
            import operator
            def is_known_op(lhs, op, rhs):
                try:
                    with field.frozen():
                        if op(lhs, value):
                            return True
                except ParametricRealFieldFrozenError:
                    return False
            # (1) if one element is known to be <= value, then nothing to record.
            if any(is_known_op(v, operator.le, value) for i, v in iv_list):
                return True
            # (2) in one linear pass, filter out elements known to be > value.
            iv_list = [ iv for iv in iv_list
                        if not is_known_op(iv[1], operator.gt, value) ]
            # (3) fast path for trivial cases
            if not iv_list:
                return True
            if len(iv_list) == 1:
                return iv_list[0][1] <= value    # records
            # (4) Find all candidates for min according to test point, and
            # filter out all elements that are known >= some of these.
            with field.off_the_record():
                small = []
                large = []
                for iv in iv_list:
                    if iv[1] == min_v:
                        small.append(iv)
                    else:
                        large.append(iv)
            # (5) Attempt to self-reduce the list of candidates.
            small_iv = small[0]
            small = [ small_iv ] + [ iv for iv in small
                                     if not is_known_op(small_iv[1], operator.le, iv[1]) ]
            ##import pdb; pdb.set_trace()
            # (6) Reduce large by small
            iv_list = small + [ iv for iv in large
                                if not any(is_known_op(small_iv[1], operator.le, iv[1]) for small_iv in small) ]
            # (7) test fast path again
            if len(iv_list) == 1:
                return iv_list[0][1] <= value    # records
            # (8)
            # new region = current region \cap ({iv_list[0][1] <= value} \cup ... \cup {iv_list[-1][1] <= value})
            # when not is_le_satisfied, asserting each iv_list[0][1] <= value would cut the current region of the field.
            # Let X = current region \cap {iv_list[0][1] > value} \cap ... \cap {iv_list[-1][1]>value is empty,
            # don't raise error if X is empty or if X = current region \cap {iv_list[0][i] > value} for some i.
            if field._mutable_values:
                # This requires a ParametricRealField set up with ``mutable_values=True``.
                with field.removed_test_point():
                    with field.temporary_assumptions():
                        for iv in iv_list:
                            try:
                                field.assume_comparison(iv[1].sym(), operator.gt, value)
                            except ParametricRealFieldInconsistencyError:
                                ### FIXME: Isn't this the same that is tested in (1) above?
                                return True
            # (9) don't raise error if X = current region \cap {iv_list[0][i] > value} for some i.
            # that is, if one iv is enough.
            for iv in iv_list:
                with field.off_the_record():
                    if iv[1] > value:
                        continue
                # at least one of the iv is <= value. 
                adding_iv_is_enough = True
                for iv_other in iv_list:
                    if iv is iv_other:
                        continue
                    with field.off_the_record():
                        if iv_other[1] <= value:
                            iv_other_le_value = True
                        else:
                            iv_other_le_value = False
                    if iv_other_le_value:
                        with field.temporary_assumptions():
                            assert iv_other[1] <= value # record temporary assumption
                            try:
                                with field.frozen():
                                    assert iv[1] <= value # always true becaus of the above check.
                            except ParametricRealFieldFrozenError:
                                adding_iv_is_enough = False
                                break
                    else:
                        #This requires a ParametricRealField set up with ``mutable_values=True``.
                        if not field._mutable_values:
                            raise ParametricRealFieldRefinementError("is_min_le - this case would require a ParametricRealField set up with mutable_values=True")
                        with field.removed_test_point():
                            with field.temporary_assumptions():
                                try:
                                    field.assume_comparison(iv_other[1].sym(), operator.le, value)
                                    with field.frozen():
                                        field.assume_comparison(iv[1].sym(), operator.le, value)
                                except (ParametricRealFieldInconsistencyError, ParametricRealFieldFrozenError):
                                    adding_iv_is_enough = False
                                    break
                if adding_iv_is_enough:
                    break # with adding_iv_is_enough = True, and iv
            if adding_iv_is_enough:
                assert iv[1] <= value # record
            else:
                raise ParametricRealFieldRefinementError("is_min_le")
    return is_le

def big_cells_sorted(iterable, key=None, reverse=False, field=None):
    """
    EXAMPLES::

        sage: import logging
        sage: logging.disable(logging.WARNING)
        sage: from cutgeneratingfunctionology import igp
        sage: from cutgeneratingfunctionology.igp import ParametricRealField
        sage: from cutgeneratingfunctionology.spam import big_cells
        sage: K.<x,y> = ParametricRealField([3/2, 1/2], big_cells=False, allow_refinement=True)
        sage: a = 2*x; b=4*y; c=x*x+y*y-4
        sage: big_cells.sorted([0,a,b,c])
        [(x^2 + y^2 - 4)~, 0, (4*y)~, (2*x)~]
        sage: sorted(K.get_lt_factor())
        [-4*y, -2*x, -2*x + 4*y, x^2 + y^2 - 4, x^2 + y^2 - 4*y - 4]
        sage: sorted(K._bsa.eq_poly()), sorted(K._bsa.lt_poly()), sorted(K._bsa.le_poly())
        ([], [-y, -x + 2*y, x^2 + y^2 - 4], [])
        sage: K.<x,y> = ParametricRealField([3/2, 1/2], big_cells=True, allow_refinement=True)
        sage: a = 2*x; b=4*y; c=x*x+y*y-4
        sage: big_cells.sorted([0,a,b,c])
        [(x^2 + y^2 - 4)~, 0, (4*y)~, (2*x)~]
        sage: sorted(K.get_le_factor())
        [-4*y, -2*x + 4*y, x^2 + y^2 - 4]
        sage: sorted(K._bsa.eq_poly()), sorted(K._bsa.lt_poly()), sorted(K._bsa.le_poly())
        ([], [], [-y, -x + 2*y, x^2 + y^2 - 4])

    TESTS:

        Key is supported::

            sage: big_cells.sorted([1, 0, -2], key=abs)
            [0, 1, -2]
            sage: K.<x,y> = ParametricRealField([3/2, 1/2], big_cells=False, allow_refinement=True)
            sage: big_cells.sorted([K(1), K(0), K(-2)], key=abs)
            [0, 1, -2]
    """
    iterable = list(iterable)
    if len(iterable) <= 1:                 # fast path
        return iterable
    if key is None:
        key = lambda i: i
    iv_list = [ (i, key(i)) for i in iterable ]
    if field is None:
        field = _common_parametric_real_field(iv_list, key=lambda iv: iv[1])
    from cutgeneratingfunctionology.igp import ParametricRealField
    if not isinstance(field, ParametricRealField):
        field = trivial_parametric_real_field
    if not field._big_cells:
        return sorted(iterable, key=key, reverse=reverse)
    with field.off_the_record():
        sorted_list = sorted(iterable, key=key)
    if field._allow_refinement:
        assert all(key(sorted_list[i-1]) <= key(sorted_list[i]) for i in range(1,len(sorted_list)))   # records
    else:
        from cutgeneratingfunctionology.igp import ParametricRealFieldFrozenError, ParametricRealFieldRefinementError
        try:
            with field.frozen():
                assert all(key(sorted_list[i-1]) <= key(sorted_list[i]) for i in range(1,len(sorted_list)))
        except ParametricRealFieldFrozenError:
            raise ParametricRealFieldRefinementError("big_cells_sorted")
        # cannot have AssertionError because we sorted the list
    if reverse:
        return list(reversed(sorted_list))
    else:
        return sorted_list
