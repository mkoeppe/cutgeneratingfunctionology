def simple_finite_dimensional_extremality_test(fn, show_plots=False, m=3, irr_grid_nb=None):
    """
    simple finite dimensional extremality test that does not go through the whole machinery of covered intervals etc., but rather sets up a grid 1/mq, where q is the lcm of the breakpoint denominator and m is the "oversampling factor" (we prove in Equivariant I that m=4 always is enough; in fact m = 3 suffices; one of these two values should be the default). The code should allow, for the irrational case, to specify a different grid (as the breakpoints don't have a gcd in this case). This code can then be easily extended to the discontinuous case! (I would be quite interested to know if this code would be able to detect any perturbation for some non-extreme "irrational" functions.)...
    """
    f = find_f(fn)
    field = fn(0).parent()
    bkpt = fn.end_points()
    is_rational_bkpt, bkpt = is_all_QQ(bkpt)
    if is_rational_bkpt:
        bkpt_denominator = [denominator(x) for x in bkpt]
        q = lcm(bkpt_denominator)
        grid_nb = m * q
        logging.info("Rational breakpoints case; use grid = 1 / %s * %s" % (m, q))
    else:
        # FIXME: In irrational case, how can we ensure that f is on the grid?
        # We need that because there is the equation pi(f) = 1.
        # Maybe, assume that f is rational and f * irr_grid_nb is integer?
        grid_nb = irr_grid_nb
        logging.info("Irrational breakpoints case; use grid = 1 / %s" % (irr_grid_nb))
    grid = field(1 / grid_nb)
    pts = [grid * i for i in range(grid_nb)]
    values = [fn(pt) for pt in pts]
    # values = [fn(i / grid_nb) for i in range(grid_nb)]
    f_grid_index = int(f * grid_nb)
    if not minimality_test_general(fn, f):
        logging.info("The function is NOT extreme.")
        return False
    basis_vects = VectorSpace(field,grid_nb).basis()
    # equation pi(0) = 0 and equation pi(f) = 1.
    equation_list = [list(basis_vects[0]), list(basis_vects[f_grid_index])]
    for i in range(grid_nb):
        for j in range(i, grid_nb, 1):
            k = (i + j) % grid_nb
            if values[i] + values[j] == values[k]:
                equation = list(basis_vects[i] + basis_vects[j] - basis_vects[k])
                equation_list.append(equation)
    equations = matrix(field, equation_list, sparse=True)
    # FIXME: When the grid is fine or there are many additive pairs, the number of equations is enormous. Very slow.
    # Maybe, should use the fact that equations is a sparse matrix.
    # Or compute rank instead?  extreme iff equations.rank() == grid_nb
    solutions = equations.right_kernel().basis()
    logging.info("Solution space has dimension %s" % len(solutions))
    if len(solutions) == 0:
        logging.info("Thus the function is extreme.")
        return True
    else:
        # FIXME: How to generate pertubation from solutions? Does the following always work (at least for continuous case)?
        for sol_index in range(len(solutions)):
            solution = list(solutions[sol_index])
            perturbation = fn._perturbation = piecewise_function_from_breakpoints_and_values(pts+[1], solution+[0], field)
            check_perturbation(fn, perturbation, show_plots=show_plots, legend_title="Basic perturbation %s" % (sol_index + 1))
        return False


