from six.moves import range
def generate_perturbations_simple(fn, show_plots=False, f=None, oversampling=3, order=None, full_certificates=True):
    r"""
    Generate (with "yield") perturbations for ``simple_finite_dimensional_extremality_test``. 
    """
    if f is None:
        f = find_f(fn, no_error_if_not_minimal_anyway=True)
    is_discontinuous = not fn.is_continuous() and not fn.is_discrete()
    field = fn(0).parent().fraction_field()
    grid_nb = finite_group_order_from_function_f_oversampling_order(fn, f=f, oversampling=oversampling, order=order)
    grid = field(1 / grid_nb)
    pts = [grid * i for i in range(grid_nb)]
    values = [fn(pt) for pt in pts]
    f_grid_index = int(f * grid_nb)
    if fn.is_continuous() or fn.is_discrete():
        basis_vects = VectorSpace(field,grid_nb).basis()
        # equation pi(0) = 0 and equation pi(f) = 1.
        equation_list = [list(basis_vects[0]), list(basis_vects[f_grid_index])]
        for i in range(grid_nb):
            for j in range(i, grid_nb, 1):
                k = (i + j) % grid_nb
                if values[i] + values[j] == values[k]:
                    equation = list(basis_vects[i] + basis_vects[j] - basis_vects[k])
                    equation_list.append(equation)
    else:
        if not order is None:
            try:
                q = finite_group_order_from_function_f_oversampling_order(fn)
            except ValueError:
                q = order
            oversampling = int(grid_nb / q)
        q = int(grid_nb / oversampling)
        basis_vects = VectorSpace(field, grid_nb + 2*q).basis()
        # notice that discontinuity can only appear at k/q, we record left and right limits at these points.
        # suppose grid_nb = q * oversampling
        # the first grid_nb components in basis_vects correspond to values at [0, 1/grid_nb, 2/grid_nb, ..., (grid_nb-1)/grid_nb]]
        # the next q components correspond to right-limits at 0, 1/q, 2/q, ..., 1- 1/q
        # the last q components correspond to left-limits at 0 (which is also 1), 1/q, 2/q, .. 1 - 1/q
        mid = basis_vects[0: grid_nb]
        right = basis_vects[grid_nb: grid_nb + q]
        left = basis_vects[grid_nb + q: grid_nb + 2*q]
        equation_list = [list(mid[0]), list(mid[f_grid_index])]
        for i in range(grid_nb):
            for j in range(i, grid_nb, 1):
                k = (i + j) % grid_nb
                if values[i] + values[j] == values[k]:
                    equation_list.append(list(mid[i] + mid[j] - mid[k]))
        limits = [fn.limits(pt) for pt in pts[0::oversampling]]
        for i in range(q):
            if limits[i][0] == limits[i][1]:
                equation_list.append(list(mid[i * oversampling] - right[i]))
            else:
                for k in range(grid_nb):
                    if k % oversampling != 0:
                        l = (i*oversampling + k) % grid_nb
                        if limits[i][1] + values[k] == values[l]:
                            equation_list.append(list(right[i] + mid[k] - mid[l]))
                for k in range(q):
                    l = (i + k) % q
                    if limits[i][1] + limits[k][0] == limits[l][1]:
                        equation_list.append(list(right[i] + mid[k * oversampling] - right[l]))
                for k in range(i + 1):
                    l = (i + k) % q
                    if limits[i][1] + limits[k][1] == limits[l][1]:
                        equation_list.append(list(right[i] + right[k] - right[l]))
                    if limits[i][1] + limits[k][-1] == limits[l][1]:
                        equation_list.append(list(right[i] + left[k] - right[l]))
                    if limits[i][1] + limits[k][-1] == limits[l][-1]:
                        equation_list.append(list(right[i] + left[k] - left[l]))
                    if limits[i][1] + limits[k][-1] == limits[l][0]:
                        equation_list.append(list(right[i] + left[k] - mid[l * oversampling]))
            if limits[i][0] == limits[i][-1]:
                equation_list.append(list(mid[i * oversampling] - left[i]))
            else:
                for k in range(grid_nb):
                    if k % oversampling != 0:
                        l = (i*oversampling + k) % grid_nb
                        if limits[i][-1] + values[k] == values[l]:
                            equation_list.append(list(left[i] + mid[k] - mid[l]))
                for k in range(q):
                    l = (i + k) % q
                    if limits[i][-1] + limits[k][0] == limits[l][-1]:
                        equation_list.append(list(left[i] + mid[k * oversampling] - left[l]))
                for k in range(i + 1):
                    l = (i + k) % q
                    if limits[i][-1] + limits[k][-1] == limits[l][-1]:
                        equation_list.append(list(left[i] + left[k] - left[l]))
                    if limits[i][-1] + limits[k][1] == limits[l][1]:
                        equation_list.append(list(left[i] + right[k] - right[l]))
                    if limits[i][-1] + limits[k][1] == limits[l][-1]:
                        equation_list.append(list(left[i] + right[k] - left[l]))
                    if limits[i][-1] + limits[k][1] == limits[l][0]:
                        equation_list.append(list(left[i] + right[k] - mid[l * oversampling]))

    equations = matrix(field, equation_list, sparse=True)
    # FIXME: When the grid is fine or there are many additive pairs, the number of equations is enormous. Very slow.
    # Maybe, should use the fact that equations is a sparse matrix.
    # Or compute rank instead?  extreme iff equations.rank() == grid_nb
    solutions = equations.right_kernel().basis()
    logging.info("Solution space has dimension %s" % len(solutions))
    # FIXME: cache data up to here.
    for sol_index in range(len(solutions)):
        if not full_certificates:
            yield None
            return
        solution = list(solutions[sol_index])
        # print solution
        if fn.is_discrete():
            perturbation = discrete_function_from_points_and_values(pts+[1], solution+[0], field)
        elif fn.is_continuous():
            perturbation = piecewise_function_from_breakpoints_and_values(pts+[1], solution+[0], field)
        else:
            mid = solution[0: grid_nb]
            right = solution[0: grid_nb]
            left = solution[0: grid_nb]
            for i in range(q):
                right[i * oversampling] = solution[grid_nb + i]
                left[i * oversampling] = solution[grid_nb + q + i]
            perturbation = discontinuous_interpolation(pts, mid, right, left)
        yield perturbation

def discontinuous_interpolation(pts, mid, right, left, merge=True):
    grid_nb = len(pts)
    pieces = []
    for i in range(grid_nb - 1):
        pieces += [ singleton_piece(pts[i], mid[i]), \
                    open_piece((pts[i], right[i]), (pts[i+1], left[i+1])) ]
    pieces += [ singleton_piece(pts[-1], mid[-1]), \
                open_piece((pts[-1], right[-1]), (1, left[0])), \
                singleton_piece(1, 0) ]
    return FastPiecewise(pieces, merge=merge)

def simple_finite_dimensional_extremality_test(fn, show_plots=False, f=None, oversampling=3, order=None, show_all_perturbations=False, full_certificates=True):
    r"""
    Simple finite dimensional extremality test for fn that does not go
    through the whole machinery of covered intervals etc., but rather
    sets up a grid 1/mq, where q is the lcm of the breakpoint
    denominator and m is the oversampling factor. We prove in
    Equivariant I that m = 4 always is enough for continuous functions; 
    in the survey we prove that, in fact, m = 3 (the default) suffices.

    Instead of the oversampling factor, the code also allows to
    specify a different grid using parameter order; this is particularly 
    useful for the irrational case.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO) # to disable output in automatic tests.
        sage: h = piecewise_function_from_breakpoints_and_values([0, 1/2, 1], [0, 1, 0])
        sage: simple_finite_dimensional_extremality_test(h, False)
        True
        sage: simple_finite_dimensional_extremality_test(drlm_not_extreme_1(), False, oversampling=1)
        True
        sage: simple_finite_dimensional_extremality_test(drlm_not_extreme_1(), False, oversampling=2)
        False
        sage: simple_finite_dimensional_extremality_test(minimal_has_uncovered_interval(), False)
        False
        sage: simple_finite_dimensional_extremality_test(drlm_not_extreme_2(), False)
        False
    """
    if show_all_perturbations is None:
        show_all_perturbations = show_plots
    if f is None:
        f = find_f(fn, no_error_if_not_minimal_anyway=True)
    if not minimality_test(fn, show_plots=show_plots, f=f):
        logging.info("Not minimal, thus NOT extreme.")
        return False
    if not full_certificates and fn.is_continuous() and number_of_slopes(fn) == 2:
        logging.info("Gomory-Johnson's 2-slope theorem applies. The function is extreme.")
        return True
    seen_perturbation = False
    fn._perturbations = []
    for index, perturbation in enumerate(generate_perturbations_simple(fn, show_plots=show_plots, f=f, oversampling=oversampling, order=order, full_certificates=full_certificates)):
        if not full_certificates:
            logging.info("The function is NOT extreme.")
            return False
        fn._perturbations.append(perturbation)
        check_perturbation(fn, perturbation,
                           show_plots=show_plots, show_plot_tag='perturbation-%s' % (index + 1),
                           legend_title="Basic perturbation %s" % (index + 1))
        if not seen_perturbation:
            seen_perturbation = True
            logging.info("Thus the function is NOT extreme.")
            if not show_all_perturbations:
                break

    if not seen_perturbation:
        if fn.is_discrete():
            logging.info("Thus the function is extreme.")
        elif oversampling is not None and oversampling >= 3 and order is None and fn.is_continuous():
            logging.info("Because oversampling >= 3, this means that the function is extreme.")
    return not seen_perturbation
