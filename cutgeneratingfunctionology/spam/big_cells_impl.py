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
            trivial_parametric_real_field = ParametricRealField()
        return trivial_parametric_real_field

def big_cells_min(iterable, key=None, field=None):
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
        sage: sorted(K._le_factor, key=str)
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

    """
    if key is None:
        key = lambda i: i
    iv_list = [ (i, key(i)) for i in iterable ]
    if field is None:
        field = _common_parametric_real_field(iv_list, key=lambda iv: iv[1])
    from cutgeneratingfunctionology.igp import ParametricRealField
    if not isinstance(field, ParametricRealField):
        field = trivial_parametric_real_field
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
        sage: K.<a,b> = ParametricRealField([2, 1], big_cells=True, allow_refinement=True)
        sage: big_cells.is_min_le([a, b], 3)
        True
        sage: K._le_factor, K._lt_factor
        ({b - 3}, set())

    This is (necessarily) a refinement of the actual non-convex
    region where min(a, b) <= 3, but a bigger cell than what is
    obtained by the following::

        sage: K.<a,b> = ParametricRealField([2, 1], big_cells=True, allow_refinement=True)
        sage: min([a, b]) <= 3
        True
        sage: K._le_factor, K._lt_factor
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
    if field is None:
        field = _common_parametric_real_field(iv_list, key=lambda iv: iv[1])
    with field.off_the_record():
        min_i, min_v = min(iv_list, key=lambda iv: iv[1])
        is_le = min_v <= value
    if field._allow_refinement:
        if is_le:
            assert min_v <= value    # records
        else:
            assert all(iv[1] > value for iv in iv_list)   # records
    else:
        from cutgeneratingfunctionology.igp import ParametricRealFieldFrozenError, ParametricRealFieldRefinementError
        if is_le:
            is_le_satisfied = False
            for iv in iv_list:
                try:
                    with field.frozen():
                        assert iv[1] <= value
                    is_le_satisfied = True
                    break
                except ParametricRealFieldFrozenError:
                    pass
            if not is_le_satisfied:
                raise ParametricRealFieldRefinementError("is_min_le")
        else:
            try:
                with field.frozen():
                    assert all(iv[1] > value for iv in iv_list)
            except ParametricRealFieldFrozenError:
                raise ParametricRealFieldRefinementError("is_min_le")
    return is_le
