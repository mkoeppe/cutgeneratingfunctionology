# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

def simple_finite_dimensional_extremality_test(fn, show_plots=False, f=None, oversampling=3, order=None):
    """
    simple finite dimensional extremality test that does not go through the whole machinery of covered intervals etc., but rather sets up a grid 1/mq, where q is the lcm of the breakpoint denominator and m is the "oversampling factor" (we prove in Equivariant I that m=4 always is enough; in fact m = 3 suffices; one of these two values should be the default). The code should allow, for the irrational case, to specify a different grid (as the breakpoints don't have a gcd in this case). This code can then be easily extended to the discontinuous case! (I would be quite interested to know if this code would be able to detect any perturbation for some non-extreme "irrational" functions.)...
    """
    if f is None:
        f = find_f(fn, no_error_if_not_minimal_anyway=True)
    field = fn(0).parent()
    grid_nb = finite_group_order_from_function_f_oversampling_order(fn, f=f, oversampling=oversampling, order=order)
    grid = field(1 / grid_nb)
    pts = [grid * i for i in range(grid_nb)]
    # or range(grid_nb+1)? It doesn't matter for continuous function, since if minimality_test is pass, fn(1) must equal to 0
    values = [fn(pt) for pt in pts]
    # values = [fn(i / grid_nb) for i in range(grid_nb)]
    f_grid_index = int(f * grid_nb)
    if not minimality_test(fn, show_plots=show_plots, f=f):
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
        
        for sol_index in range(len(solutions)):
            solution = list(solutions[sol_index])
            if fn.is_discrete():
                perturbation = discrete_function_from_points_and_values(pts+[1], solution+[0], field)
            elif fn.is_continuous():
                perturbation = piecewise_function_from_breakpoints_and_values(pts+[1], solution+[0], field)
            else:
                 # FIXME: Handle the discontinuous case properly !
                perturbation = piecewise_function_from_breakpoints_and_values(pts+[1], solution+[0], field)
            fn._perturbation = perturbation
            check_perturbation(fn, perturbation, show_plots=show_plots, legend_title="Basic perturbation %s" % (sol_index + 1))
        return False


