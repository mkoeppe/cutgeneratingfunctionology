from cutgeneratingfunctionology.igp import *
from cutgeneratingfunctionology.spam.basic_semialgebraic import EmptyBSA
from cutgeneratingfunctionology.spam.basic_semialgebraic_pplite import BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron
import csv
import logging
import os

minimal_funciton_cell_description_logger = logging.getLogger("cutgeneratingfunctionology.igp.minimal_funciton_cell_description")
minimal_funciton_cell_description_logger.setLevel(logging.INFO)

### Note to future reader: ###
### bkpt is assumed to be a breakpoint sequence of length n>= 2.
### Breakpoint sequence are sorted lists of real numbers in [0,1).
### A breakpoint sequences should always have 0 as an element.
### This is never strictly enforced in this file and it is assumed that
### the user is always provided a breakpoint sequence.

class RepElemGenFailure(Exception):
    pass


def nnc_poly_from_bkpt_sequence(bkpt, backend=None, log_paramateric_real_field = False,log_pw_functions = False):
    r"""
    Defines an NNC polyhedron P such that for all b in P the delta complex of b is isomoprhic to the delta complex of bkpt.

    INPUT:
    - ``bkpt`` - sorted list of length 1 or more of sage type in the interval [0,1).
    - ```backend`` - ``None``, ``str(pplite)``

    OUTPUT: class::``BasicSemialgebraicSet_veronese``

    EXAMPLES::
    sage: from cutgeneratingfunctionology.igp import *
    sage: logging.disable(logging.INFO) # suppress logging for tests
    sage: nnc_poly_from_bkpt_sequence([0, 4/5])
    BasicSemialgebraicSet_veronese(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(Constraint_System {x0==0, -x1+1>0, 2*x1-1>0}, names=[x0, x1]), polynomial_map=[lambda0, lambda1])
    """
    n = len(bkpt)
    coord_names = []
    bkpt_vals = bkpt
    vals = bkpt_vals[0:n]
    bkpt_extd = list(bkpt)+[1]
    if not log_paramateric_real_field:
        parametric_logging_level = logging.getLogger("cutgeneratingfunctionology.igp.parametric").getEffectiveLevel()
        logging.getLogger("cutgeneratingfunctionology.igp.parametric").setLevel(logging.ERROR)
    for i in range(0,n):
        coord_names.append('lambda'+str(i))
    K = ParametricRealField(names=coord_names, values=vals, mutable_values=True, big_cells=True, default_backend=backend)
    K.gens()[0] == 0
    for i in range(n-1):
        K.gens()[i] < K.gens()[i+1]
    K.gens()[n-1] < 1
    for i in range(n):
        for j in range(n):
            if bkpt[i]+bkpt[j] >= 1:
                w = 1
            else:
                w = 0
            for k in range(n):
                if bkpt_extd[k] < bkpt[i]+bkpt[j] - w and bkpt[i]+bkpt[j] - w < bkpt_extd[k+1]:
                    if k != n-1:
                        K.gens()[k] < K.gens()[i] + K.gens()[j] - w
                        K.gens()[i] + K.gens()[j] - w < K.gens()[k+1]
                    else:
                        K.gens()[k] < K.gens()[i] + K.gens()[j] - w
                        K.gens()[i] + K.gens()[j] - w < 1
                elif bkpt_extd[k] == bkpt[i]+bkpt[j] - w:
                    if k != n-1:
                        K.gens()[k] == K.gens()[i] + K.gens()[j] - w
                        K.gens()[i] + K.gens()[j] - w < K.gens()[k+1]
                    else:
                        K.gens()[k] == K.gens()[i] + K.gens()[j] - w
                        K.gens()[i] + K.gens()[j] - w < 1
    if not log_paramateric_real_field:
        logging.getLogger("cutgeneratingfunctionology.igp.parametric").setLevel(parametric_logging_level)
    return K._bsa

def add_breakpoints_and_find_equiv_classes(bkpt_poly):
    r"""
    Takes dim k-1 breakpoint NNC polyhedron (as a :class:`BasicSemialgebraicSet_base`), adds a dimension,
    and finds all possible representative elements for equivalence classes of polyhedral complexes.

    INPUT: :class:`BasicSemialgebraicSet_Polyhedral`

    OUTPUT: unique list of breakpoint sequnes of length k (as tuples).

    EXAMPLES::

    sage: from cutgeneratingfunctionology.igp import *
    sage: logging.disable(logging.INFO) # suppress logging for tests
    sage: add_breakpoints_and_find_equiv_classes(nnc_poly_from_bkpt_sequence([0,4/5]).upstairs())
    [(0, 9/14, 83/112), (0, 13/20, 33/40), (0, 9/14, 101/112), (0, 7/10, 17/20)]
    """
    # BSAs are highly mutable, work only with copies.
    B_cap_N_b = copy(bkpt_poly)
    B_cap_N_b.add_space_dimensions_and_embed(1)
    # get new number of breakpoints.
    k = B_cap_N_b.ambient_dim()
    if k < 1:
        raise ValueError("bkpt_poly should have space dim at least 1.")
    model_bound_bkpts = [0]*k
    model_bound_bkpts[k-1] = 1
    # 0 < lambda_k < 1
    # we assume the bsa is "polyhedra type"
    B_cap_N_b.add_linear_constraint(model_bound_bkpts, QQ(-1), operator.lt) # model bounds
    B_cap_N_b.add_linear_constraint(model_bound_bkpts, QQ(0), operator.gt) # model bounds
    bkpt_order = [0]*k
    bkpt_order[k-2] = 1
    bkpt_order[k-1] = -1
    B_cap_N_b.add_linear_constraint(bkpt_order, QQ(0), operator.lt) # order on bkpts
    rep_elems = []
    for j in range(k-1):
        for i in range(k):
            for interval_w in [0,1]:
                for line_w in [0,1]:
                    # which interval is (lambda_k,lambda_k) located in?
                    # modeled lambda_i op 2lambda_k - w < lambda_{i+1}
                    for interval_op in [operator.lt, operator.eq, operator.gt]:
                        for line_op in [operator.lt, operator.eq, operator.gt]:
                            # highly mutable objects, operate on the copy.
                            B_cap_N_b_copy = copy(B_cap_N_b)
                            lhs_i = [0]*k
                            lhs_i[k-1] = -2
                            lhs_i[i] = 1
                            B_cap_N_b_copy.add_linear_constraint(lhs_i, QQ(interval_w), interval_op)
                            lhs_i_plus_1 = [0]*k
                            lhs_i_plus_1[k-1] = -2
                            if i < k-1:
                                lhs_i_plus_1[i+1] = 1
                                B_cap_N_b_copy.add_linear_constraint(lhs_i_plus_1, QQ(interval_w), operator.gt)
                            else:
                                B_cap_N_b_copy.add_linear_constraint(lhs_i_plus_1, QQ(interval_w + 1), operator.gt)
                            if not B_cap_N_b_copy.is_empty():
                                # does the line x+y equiv lambda_k mod 1 lie on/above/below (lambda_j,lambda_j)?
                                # modeled by  2lambda_j op lambda_k + w
                                lhs_j = [0]*k
                                lhs_j[j] = 2
                                lhs_j[k-1] = -1
                                B_cap_N_b_copy.add_linear_constraint(lhs_j, QQ(-line_w), line_op)
                                try:
                                    rep_elem = B_cap_N_b_copy.find_point()
                                    rep_elems.append(tuple(rep_elem))
                                except EmptyBSA:
                                    pass
    return unique_list(rep_elems)

def make_rep_bkpts_with_len_n(n, k=1, bkpts=None, backend=None):
    r"""
    Produce representative elements of every isomorphism class of breakpoints complexes for breakpoint sequences of length n.

    INPUT:
    - n, integer, maximum length of breakpoint sequence.
    - k, assumed length of breakpoint sequences in ``bkpts``.
    - bkpts, list of breakpoint sequences of length k.

    OUTPUT: A list of representative elements of every isomorphism class of breakpoints complexes for breakpoint sequences of length n extrapolated from bkpts.

    EXAMPLES::

    sage: from cutgeneratingfunctionology.igp import *
    sage: logging.disable(logging.INFO) # suppress logging for tests
    sage: make_rep_bkpts_with_len_n(2)
    [(0, 1/2), (0, 13/18), (0, 5/18)]

    The number of representative elements grows quickly::

    sage: bkpts_rep_with_len_3 = make_rep_bkpts_with_len_n(3)
    sage: len(bkpts_rep_with_len_3)
    34

    Previous computations can be reused::

    sage: bkpts_rep_with_len_4 = make_rep_bkpts_with_len_n(4, 3, bkpts_rep_with_len_3)
    sage: len(bkpts_rep_with_len_4)
    329
    """
    # Matthias has suggested looking at a directed tree.
    # An alternative approach would be to look into using a (graded) lattice as a data structure.
    # We have bkpt \leq bkpt' if and only if dim(NNC(bkpt)) \leq dim(NNC(bkpt'))
    # and embed(NNC(bkpt), dim(NNC(bkpt')) \leq NNC(bkpt') in the poset of NNC polyhedra.
    # This might speed up/less the load of verifying unquiness of cells which is the time bounding task here.
    new_bkpts = []
    if n < 2:
        raise ValueError("n>=2")
    if k == n and bkpts is not None:
        minimal_funciton_cell_description_logger.warning(f"Initial inputs suggest that the bkpts provided are already correct. Returning the initial bkpts.")
        return bkpts
    if k == n and bkpts is None:
        raise ValueError("k<n")
    if k == 1 and bkpts is None:
        bkpts = [[0]]
    for bkpt in bkpts:
        new_bkpts += add_breakpoints_and_find_equiv_classes(nnc_poly_from_bkpt_sequence(bkpt, backend=backend).upstairs())
    new_bkpts = unique_list(new_bkpts)
    k += 1
    if k == n:
        minimal_funciton_cell_description_logger.info(f"Breakpoints of length {n} have been generated.")
        return new_bkpts
    else:
        minimal_funciton_cell_description_logger.info(f"Breakpoints of length {k} have been generated. Now generating breakpoints of length {k+1}.")
        return make_rep_bkpts_with_len_n(n, k, new_bkpts, backend)

def generate_assumed_symmetric_vertices_continuous(fn, f, bkpt):
    """
    Silently assumes the symmetry condition holds for all vertices (x,y) in bkpt's breakpoints complex
    such that x+y equiv f.

    See function::``generate_symmetric_vertices_continuous``.
    """
    for i in range(len(bkpt)):
        x = bkpt[i]
        if x == f:
            continue
        if x < f:
            y = f - x
        else:
            y = 1 + f - x
        fn(x) + fn(y) == 1
        yield (x, y, 0, 0)

def value_nnc_polyhedron_value_cords(bkpt, f_index, backend=None, log_paramateric_real_field = False, log_pw_functions = False):
    r"""
    For a given breakpoints sequence and f index, write the value polyhedron as a basic semialgebraic set in only the value parameters.

    INPUT:
    - ``bkpt`` - a breakpoinnt sequence of length 2 or more.
    - ``f_index`` - integer between 1 and length of ``len(bkpt) -1``.
    - ```backend`` - ``None``, ``str(pplite)``

    OUTPUT:
    - class::``BasicSemialgebraicSet_veronese`` - The value polyehdron in only the value parameters.

    EXAMPLES::

    sage: from cutgeneratingfunctionology.igp import *
    sage: logging.disable(logging.INFO) # suppress logging for tests
    sage: value_nnc_polyhedron_value_cords([0,4/5], 1) # gmic with f=4/5
    BasicSemialgebraicSet_veronese(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(Constraint_System {x1-1==0, x0==0}, names=[x0, x1]), polynomial_map=[gamma0, gamma1])
    """
    # this saves a slight amount of overhead when detemrining points for the value polyhedron since the assumed.
    # minimality test does not have to entirely go though parametric real field.
    # This is useful in application.
    n = len(bkpt)
    if not log_paramateric_real_field:
        parametric_logging_level = logging.getLogger("cutgeneratingfunctionology.igp.parametric").getEffectiveLevel()
        logging.getLogger("cutgeneratingfunctionology.igp.parametric").setLevel(logging.ERROR)
    if not log_pw_functions:
        pw_logging_level = logging.getLogger("cutgeneratingfunctionology.igp.functions").getEffectiveLevel()
        logging.getLogger("cutgeneratingfunctionology.igp.functions").setLevel(logging.ERROR)
    assert(n >= 2)
    assert(f_index >= 1)
    assert(f_index <= n - 1)
    if not isinstance(bkpt, list):
        bkpt = list(bkpt)
    coord_names = []
    val = [None]*(n)
    for i in range(n):
        coord_names.append('gamma'+str(i))
    K = ParametricRealField(names=coord_names, values=val, mutable_values=True, big_cells=True, allow_refinement=False, default_backend=backend)
    K.gens()[0] == 0
    for i in range(1, n):
        K.gens()[i] <= 1
        K.gens()[i] > 0
    h = piecewise_function_from_breakpoints_and_values(bkpt + [1], K.gens() + [0], merge=False)
    # Assumes minimality for the partially defined function.
    for vert in generate_type_1_vertices_continuous(h, operator.ge, bkpt + [1]):
        vert
    for vert in generate_type_2_vertices_continuous(h, operator.ge, bkpt + [1]):
        vert
    for vert in generate_assumed_symmetric_vertices_continuous(h, bkpt[f_index], bkpt + [1]):
        vert
    if not log_paramateric_real_field:
        logging.getLogger("cutgeneratingfunctionology.igp.parametric").setLevel(parametric_logging_level)
    if not log_pw_functions:
        logging.getLogger("cutgeneratingfunctionology.igp.functions").setLevel(pw_logging_level)
    return K._bsa

def value_nnc_polyhedron(bkpt, f_index, backend=None, log_paramateric_real_field = False,log_pw_functions = False):
    r"""
    For a given breakpoints sequence and f index, write the value polyhedron as a basic semialgebraic set in the full space of parameters.

    INPUTS:
    - ``bkpt`` - sorted list of length 2 or more of sage type in the interval [0,1).
    - ``f_index`` - integer between 1 and length of ``len(bkpt) -1``.
    - ``backend`` - ``None``, ``str(pplite)``

    OUTPUT:
    - class::``BasicSemialgebraicSet_veronese`` - The value polyehdron in only the breakpoint and value parameter space.

    EXAMPLES::

    sage: from cutgeneratingfunctionology.igp import *
    sage: logging.disable(logging.INFO) # suppress logging for tests
    sage: value_nnc_polyhedron([0,4/5], 1) # gmic with f=4/5, a trivial value polyhedron
    BasicSemialgebraicSet_veronese(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(Constraint_System {x3-1==0, x2==0, 5*x1-4==0, x0==0}, names=[x0, x1, x2, x3]), polynomial_map=[lambda0, lambda1, gamma0, gamma1])
    sage: P =  value_nnc_polyhedron([0,2/15,2/3, 4/5], 3) #  A non trivial value polyhedron.
    sage: P.upstairs()._polyhedron
    A 1-dimensional polyhedron in QQ^8 defined as the convex hull of 2 points
    """
    n = len(bkpt)
    if not log_paramateric_real_field:
        parametric_logging_level = logging.getLogger("cutgeneratingfunctionology.igp.parametric").getEffectiveLevel()
        logging.getLogger("cutgeneratingfunctionology.igp.parametric").setLevel(logging.ERROR)
    if not log_pw_functions:
        pw_logging_level = logging.getLogger("cutgeneratingfunctionology.igp.functions").getEffectiveLevel()
        logging.getLogger("cutgeneratingfunctionology.igp.functions").setLevel(logging.ERROR)
    assert(n >= 2)
    assert(f_index >= 1)
    assert(f_index <= n)
    coord_names = []
    bkpt_vals = list(bkpt)
    vals = bkpt_vals + [None]*(n)
    for i in range(n):
        coord_names.append('lambda'+str(i))
    for i in range(n):
        coord_names.append('gamma'+str(i))
    K = ParametricRealField(names=coord_names, values=vals, mutable_values=True, big_cells=True, allow_refinement=False, default_backend=backend)
    # breakpoint parameters are the measured breakpoint values.
    for i in range(n):
        K.gens()[i] == bkpt[i]
    # Necessary conditions on value parameters.
    K.gens()[n] == 0
    for i in range(1, n):
        K.gens()[i+n] <= 1
        K.gens()[i+n] > 0
    h = piecewise_function_from_breakpoints_and_values(K.gens()[0:n] + [1], K.gens()[n:2*n] + [0], merge=False)
    # Assumes minimality  for the partially defined function.
    for vert in generate_type_1_vertices_continuous(h, operator.ge, K.gens()[0:n] + [1]):
        vert
    for vert in generate_type_2_vertices_continuous(h, operator.ge, K.gens()[0:n] + [1]):
        vert
    for vert in generate_assumed_symmetric_vertices_continuous(h, K.gens()[f_index], [0] + K.gens()[0:n] + [1]):
        vert
    if not log_paramateric_real_field:
        logging.getLogger("cutgeneratingfunctionology.igp.parametric").setLevel(parametric_logging_level)
    if not log_pw_functions:
        logging.getLogger("cutgeneratingfunctionology.igp.functions").setLevel(pw_logging_level)
    return K._bsa

def bsa_of_rep_element(bkpt, vals, backend=None, log_paramateric_real_field = False,log_pw_functions = False):
    r"""
    Given pi_(bkpt, vals) is {minimal, not minimal}, find BSA subset of R^(2n) such that (bkpt, vals) in BSA and for all p
    in BSA, pi_p is {minimal, not minimal}.

    INPUT:
    - ``bkpt`` - a breakpoint sequence
    - ``vals`` - list like of sage numerical types corresponding values for the breakpoint sequence.
    - ``backend`` - None, ``str(pplite)``

    OUTPUT: A basic semialgebraic set.

    EXAMPLES::

    sage: from cutgeneratingfunctionology.igp import *
    sage: logging.disable(logging.INFO) # suppress logging for tests
    sage: bsa_of_rep_element([0,4/5], [0,1]) # bsa for GMIC
    BasicSemialgebraicSet_veronese(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(Constraint_System {x5-1==0, x4==0, x0==0, 2*x3-1>0, -2*x1+x2-x3+1>=0, -x3+1>0}, names=[x0, x1, x2, x3, x4, x5]), polynomial_map=[gamma0, lambda1^2*gamma0, lambda1*gamma0, lambda1, lambda0, gamma1])
    """
    n = len(bkpt)
    if not log_paramateric_real_field:
        parametric_logging_level = logging.getLogger("cutgeneratingfunctionology.igp.parametric").getEffectiveLevel()
        logging.getLogger("cutgeneratingfunctionology.igp.parametric").setLevel(logging.ERROR)
    if not log_pw_functions:
        pw_logging_level = logging.getLogger("cutgeneratingfunctionology.igp.functions").getEffectiveLevel()
        logging.getLogger("cutgeneratingfunctionology.igp.functions").setLevel(logging.ERROR)
    assert(n >= 2)
    coord_names = []
    for i in range(n):
        coord_names.append('lambda'+str(i))
    for i in range(n):
        coord_names.append('gamma'+str(i))
    K = ParametricRealField(names=coord_names, values=bkpt+vals, big_cells=True, default_backend=backend)
    h = piecewise_function_from_breakpoints_and_values(K.gens()[0:n] + [1], K.gens()[n:2*n] + [0], merge=False)
    minimality_test(h)
    if not log_paramateric_real_field:
        logging.getLogger("cutgeneratingfunctionology.igp.parametric").setLevel(parametric_logging_level)
    if not log_pw_functions:
        logging.getLogger("cutgeneratingfunctionology.igp.functions").setLevel(pw_logging_level)
    return K.make_proof_cell().bsa

def find_minimal_function_reps_from_bkpts(bkpts, prove_minimality=True, backend=None):
    r"""
    Finds representative elements of minimal functions from a given breakpoint sequence.

    INPUT:
    - ``bkpts`` - an iterable of breakpoint sequences.
    - ``prove_minimality`` - bool, proves minimality of parameterized function.
    - ``backend`` - None, ``str(pplite)``

    OUTPUT:
    - List of (breakpoint, value) pairs.

    EXAMPLES::

    sage: from cutgeneratingfunctionology.igp import *
    sage: logging.disable(logging.INFO) # suppress logging for tests
    sage: bkpts = make_rep_bkpts_with_len_n(2)
    sage: find_minimal_function_reps_from_bkpts(bkpts)
    [([0, 1/2], [0, 1]), ([0, 13/18], [0, 1]), ([0, 5/18], [0, 1])]

    """
    rep_elems = []
    if not log_pw_functions:
        pw_logging_level = logging.getLogger("cutgeneratingfunctionology.igp.functions").getEffectiveLevel()
        logging.getLogger("cutgeneratingfunctionology.igp.functions").setLevel(logging.ERROR)
    for bkpt in bkpts:
        n = len(bkpt)
        for f_index in range(1, n):
            poly_bsa = value_nnc_polyhedron_value_cords(list(bkpt), f_index, backend)
            gammas = poly_bsa.polynomial_map()[0].parent().gens()
            try:
                test_point = poly_bsa.upstairs().find_point()
            except EmptyBSA:
                raise RepElemGenFailure("The value polyhedron {} is empty. This should not be empty. Double check inputs".format(poly_bsa))
            test_val = []
            for gamma_i in gammas:
                test_val.append(test_point[poly_bsa.v_dict()[gamma_i]])
            if prove_minimality:
                h = piecewise_function_from_breakpoints_and_values(list(bkpt)+[1], test_val+[0])
                if not minimality_test(h): # The following error should never be raised when this function is used as intended.
                    raise ValueError(f"({bkpt}, {test_val}) parameterized by breakpoints and values is not a minimal function but assuming a breakpoint sequence is input, this should be minimal.")
            rep_elems.append((list(bkpt), test_val))
    if not log_pw_functions:
        logging.getLogger("cutgeneratingfunctionology.igp.functions").setLevel(pw_logging_level)
    return rep_elems

class BreakpointComplexClassContainer:
    """
    A container for the family of breakpoint complexes for piecewise linear functions
    with at most n breakpoints.

    The container will attempt to load from the minimalFunctionCache module or will generate
    data if the minimalFunctionCache is not available.

    Loading files generated from the container's `write_data()` method is supported.
    When manually loading data, it is assume that all loaded data is correct
    for the Initialized number of breakpoints.

    Initialization:
    n - integer 2 or larger.
    backend - (optional) `None` or `str(pplite)`
    manually_load_breakpoint_cache - (optional) bool
    loading_kwrds - file_or_folder, path_to_file_or_folder

    EXAMPLES::

    sage: from cutgeneratingfunctionology.igp import *
    sage: logging.disable(logging.INFO) # suppress logging for tests
    sage: bkpt_of_len_2 = BreakpointComplexClassContainer(2)
    sage: bkpt_of_len_2.num_rep_elems()
    3
    sage: [elem for elem in bkpt_of_len_2.get_rep_elems()]
    [(0, 1/2), (0, 13/18), (0, 5/18)]
    sage: make_rep_bkpts_with_len_n(2)
    [(0, 1/2), (0, 13/18), (0, 5/18)]

    A cell description of semialgebraic sets can be accessed::

    sage: all([isinstance(cell, BasicSemialgebraicSet_base) for cell in bkpt_of_len_2.get_nnc_poly_from_bkpt()])
    True

    (Advanced use) Data can be written and reused::

    sage: bkpt_of_len_2.write_data() # not tested
    sage: bkpt_of_len_2_loaded = BreakpointComplexClassContainer(2, manually_load_breakpoint_cache=True, file_or_folder="file", path_to_file_or_folder="bkpts_of_len_2.csv") # not tested
    """
    def __init__(self, n, backend=None, manually_load_breakpoint_cache=False, **loading_kwrds):
        self._n = n
        assert(self._n >= 2)
        self._backend = backend
        if not manually_load_breakpoint_cache:
            try:
                # Load the minimal function cache.
                from minimalFunctionCache.utils import minimal_function_cache_info, minimal_function_cache_loader
                cache_info = minimal_function_cache_info()
                # minimal_function_cache_loader throws a value error if a cache for n is not found.
                try:
                    self._data = minimal_function_cache_loader(self._n, "breakpoints")
                except ValueError:
                    minimal_funciton_cell_description_logger.info(f"The cache for {n} breakpoints has not been generated or could not be found.")
                    minimal_funciton_cell_description_logger.info("Generating representative breakpoints.")
                    avail_bkpts = [ i for i in cache_info["avail_bkpts"] if i < self._n]
                    if len(avail_bkpts) > 0:
                        k = max(avail_bkpts)
                        prev_gen_bkpts = self._data = minimal_function_cache_loader(k, "breakpoints")
                        if self._n > 4:
                            minimal_funciton_cell_description_logger.warning(f"This may take a while.")
                        self._data = make_rep_bkpts_with_len_n(self._n, k, prev_gen_bkpts, backend=self._backend)
                    else:
                        if self._n > 4:
                            minimal_funciton_cell_description_logger.warning(f"This may take a while.")
                        self._data = make_rep_bkpts_with_len_n(self._n, backend=self._backend)
            except ImportError:
                if self._n > 4:
                    minimal_funciton_cell_description_logger.warning(f"This may take a while. Try installing the minimalfunctioncache.")
                self._data = make_rep_bkpts_with_len_n(self._n, backend=self._backend)
                minimal_funciton_cell_description_logger.info("Finished generating breakpoints.")
        # For generating the cache and advanced use.
        elif manually_load_breakpoint_cache:
            try:
                if loading_kwrds["file_or_folder"].strip(" ").lower() == "folder":
                    loaded_data = []
                    path = loading_kwrds["path_to_file_or_folder"]
                    import os
                    files_in_folder = os.listdir(path)
                    for file in files_in_folder:
                        with open(os.path.join(path, file)) as csvfile:
                            file_reader = csv.reader(csvfile)
                            for row in file_reader:
                                loaded_data.append([QQ(data) for data in row])
                    k = len(loaded_data[0]) # assume all the data is correct and of the same length.
                    self._data = make_rep_bkpts_with_len_n(self._n, k, loaded_data, backend=self._backend)
                elif loading_kwrds["file_or_folder"].strip(" ").lower() == "file":
                    loaded_data = []
                    with open(loading_kwrds["path_to_file_or_folder"]) as csvfile:
                        file_reader = csv.reader(csvfile)
                        for row in file_reader:
                            loaded_data.append([QQ(data) for data in row])
                    k = len(loaded_data[0]) # assume all the data is correct and of the same length.
                    self._data = make_rep_bkpts_with_len_n(self._n, k, loaded_data, backend=self._backend)
                else:
                    raise ValueError("Check spelling of folder or file.")
            except KeyError:
                raise ValueError("No valid inputs provided. Maybe check your spelling?")
        else:
            raise ValueError("No elements have been loaded. Check inputs")

    def __repr__(self):
        return f"Container for the space breakpoint sequences of length {self._n} under equivalence of polyhedral complexes."

    def get_rep_elems(self):
        """
        Iterator yielding a breakpint sequence.
        """
        for bkpt in self._data:
            yield bkpt

    def get_nnc_poly_from_bkpt(self):
        """
        Iterator yielding a cell such that all elements in the cell correspond to 2d polyhedral complexes of the same combinatorial type.
        """
        for bkpt in self._data:
            yield nnc_poly_from_bkpt_sequence(bkpt, backend=self._backend)

    def num_rep_elems(self):
        """
        Number of repersentative elements.
        """
        return len(self._data)

    def n(self):
        """
        The length of breakpoint sequences.
        """
        return self._n

    def covers_space(self):
        """
        Not Implemented. Future use is intended to be a proof of correctness that all breakpoint sequences are covered.
        """        ### TODO: Impl.
        raise NotImplementedError

    def refine_space(self):
        """
        Not Implemented. Future use is to prove that the current cell description for breakpoints is, with respect to inclusion, minimal or in the case that it is not minimal, find a minimal description.
        """
        raise NotImplementedError
        ### TODO: Add a contaimenet checking method to BSA classes and finish impl.
        self._data = unique_list(self._data)
        cells_found = []
        for bkpt in self._data:
            bkpt_contained_in_cell = False
            for found_cell in cells_found:
                if found_cell.contains(bkpt):
                    bkpt_contained_in_cell = True
                    break
            if not contained_in_cell:
                cells_found.append(nnc_poly_from_bkpt_sequence(bkpt, backend=self._backend))
        new_data = []
        for cell in cells_found:
            new_data.append(cell.find_point())
        self._data = new_data

    def write_data(self, output_file_name_style=None, max_rows=None):
        """
        Writes representative element data to a `.csv` file with one column and rows of representative elements.
        Optionally, write many `.csv` files with at `most max_rows` rows per file.
        Files are named output_file_file_name_style_filenumber.csv.
        The default output_file_name_style="bkpts_of_len_n".
        """
        # TODO: Future, support writing different types of data such as polyhedra data.
        if output_file_name_style is None:
            file_name_base = "bkpts_of_len_{}".format(self._n)
        else:
            file_name_base = output_file_name_style
        if max_rows is not None:
            assert(max_rows >= 1)
            num_files = len(self._data)//max_rows + 1
            file_name_base = file_name_base+"_part_0"
        if max_rows is None:
            max_rows = len(self._data)
            num_files = 1
        output_file = file_name_base+".csv"
        for file_number in range(num_files):
            out_file = open(output_file, "w")
            data_writer = csv.writer(out_file, csv.QUOTE_NONE)
            for row in range(max_rows):
                try:
                    data_writer.writerow(self._data[max_rows * file_number + row])
                except IndexError:
                    break
            out_file.close()
            output_file = file_name_base[:-1]+"{}".format(file_number+1)+".csv"


class PiMinContContainer:
    """
    A container for the space of continuous piecewise linear minimal functions with at
    most n breakpoints parameterized by breakpoints and values using semialgebraic sets.

    The container will attempt to load from the `minimalFunctionCache module` or will generate
    data if the minimalFunctionCache is not available.

    Loading files generated from the container's `write_data()` method is supported.
    When manually loading data, it is assume that all loaded data is correct
    for the Initialized number of breakpoints.

    INPUTS:
    - n, an integer
    - backend, None or str(pplite)
    - loading_kwrds

    EXAMPLES::

    sage: from cutgeneratingfunctionology.igp import *
    sage: logging.disable(logging.INFO) # suppress logging for tests
    sage: PiMin_2 = PiMinContContainer(2)
    sage: all([minimality_test(pi) for pi in PiMin_2.get_rep_functions()])
    True
    sage: PiMin_2
    Space of minimal functions with at most 2 breakpoints parameterized by breakpoints and values using semialgebraic sets.

    A cell description of semialgebraic sets can be accessed::

    sage: all([isinstance(cell, BasicSemialgebraicSet_base) for cell in PiMin_2.get_semialgebraic_sets()])
    True

    Data is stored as repersenative elements. The number of repersentative elements grows quickly. ::

    sage: PiMin_4 = PiMinContContainer(4)  # not tested
    sage: len([rep_elem for rep_elem in PiMin_4.get_rep_elems()]) # not tested
    987
    sage: len([rep_elem for rep_elem in PiMin_2.get_rep_elems()])
    3

    The container provides methods of writing data.::

    sage: PiMin_2.write_data() # not tested

    Written data can be reused.::

    sage: PiMin_2_loaded_data = PiMinContContainer(2, manually_load_function_cache=True, file_or_folder="file", path_to_file_or_folder="Pi_Min_2.csv", breakpoints_or_rep_elems="rep_elems") # not tested
    sage: len([rep_elem for rep_elem in PiMin_2_loaded_data.get_rep_elems()])  # not tested
    3
    """
    def __init__(self, n, backend=None, manually_load_function_cache=False, **loading_kwrds):
        self._n = n
        assert(self._n >= 2)
        self._backend = backend
        if not manually_load_function_cache:
            try:
                # Load the minimal function cache.
                from minimalFunctionCache.utils import minimal_function_cache_loader
                try:
                    self._data = minimal_function_cache_loader(self._n, "rep_elems")
                # cache loader throws a value error if a cache for n is not found.
                except ValueError:
                    minimal_funciton_cell_description_logger.info(f"The cache for {n} breakpoints has not been generated or could not be found.")
                    minimal_funciton_cell_description_logger.info("Finding or generating representative breakpoints.")
                    bkpts = BreakpointComplexClassContainer(self._n, backend=self._backend).get_rep_elems()
                    if self._n > 4:
                        minimal_funciton_cell_description_logger.warning(f"This may take a while.")
                    minimal_funciton_cell_description_logger.info("Finished generating breakpoints.")
                    minimal_funciton_cell_description_logger.info("Computing repersentative elements.")
                    self._data = find_minimal_function_reps_from_bkpts(bkpts, backend=self._backend)
            except ImportError:
                if self._n > 4:
                    minimal_funciton_cell_description_logger.warning(f"This may take a while. Try installing the minimalFunctionCache.")
                bkpts = make_rep_bkpts_with_len_n(self._n, backend=self._backend)
                minimal_funciton_cell_description_logger.info("Finished generating breakpoints.")
                minimal_funciton_cell_description_logger.info("Computing repersentative elements.")
                self._data = find_minimal_function_reps_from_bkpts(bkpts, backend=self._backend)
            minimal_funciton_cell_description_logger.info("PiMin container, Reportin' for duty.")
        elif manually_load_function_cache:
            # this is for generating the cache and advanced use.
            minimal_funciton_cell_description_logger.info("loading files...")
            try:
                if loading_kwrds["breakpoints_or_rep_elems"].strip(" ").lower() == "breakpoints":
                    bkpts = BreakpointComplexClassContainer(self._n, backend=self._backend, manually_load_breakpoint_cache=True, file_or_folder=loading_kwrds["file_or_folder"], path_to_file_or_folder=loading_kwrds["path_to_file_or_folder"]).get_rep_elems()
                    if self._n > 4:
                        minimal_funciton_cell_description_logger.warning("This may take a while. Try installing the minimalFunctionCache.")
                    self._data = find_minimal_function_reps_from_bkpts(bkpts, backend=self._backend)
                    minimal_funciton_cell_description_logger.info("PiMin container, Reportin' for duty.")
                elif loading_kwrds["breakpoints_or_rep_elems"].strip(" ").lower() == "rep_elems":
                    if loading_kwrds["file_or_folder"].strip(" ").lower() == "folder":
                        self._data = []
                        path = loading_kwrds["path_to_file_or_folder"]
                        import os
                        files_in_folder = os.listdir(path)
                        for file in files_in_folder:
                            with open(os.path.join(path, file)) as csvfile:
                                file_reader = csv.reader(csvfile)
                                for row in file_reader:
                                    bkpt = [QQ(data) for data in row[0].strip("[]").split(",")]
                                    val = [QQ(data) for data in row[1].strip("[]").split(",")]
                                    self._data.append((bkpt, val))
                        minimal_funciton_cell_description_logger.info("PiMin container, Reportin' for duty.")
                    elif loading_kwrds["file_or_folder"].strip(" ").lower() == "file":
                        self._data = []
                        with open(loading_kwrds["path_to_file_or_folder"]) as csvfile:
                            file_reader = csv.reader(csvfile)
                            for row in file_reader:
                                bkpt = [QQ(data) for data in row[0].strip("[]").split(",")]
                                val = [QQ(data) for data in row[1].strip("[]").split(",")]
                                self._data.append((bkpt, val))
                        minimal_funciton_cell_description_logger.info("PiMin container, Reportin' for duty.")
                else:
                    raise ValueError("check spelling of folder or file.")
            except KeyError:
                raise ValueError("No valid inputs provided. Maybe check your spelling?")
        else:
            raise ValueError("No elements have been loaded. Check inputs")

    def __repr__(self):
        return "Space of minimal functions with at most {} breakpoints parameterized by breakpoints and values using semialgebraic sets.".format(self._n)

    def get_semialgebraic_sets(self):
        """Iterator for semialgebraic set description."""
        for b, v in self._data:
            yield bsa_of_rep_element(list(b), list(v), backend=self._backend)

    def get_rep_elems(self):
        """Iterator for representative elements."""
        for b, v in self._data:
            yield (list(b), list(v))

    def get_rep_functions(self):
        """Iterator for representative functions."""
        for b, v in self._data:
            yield piecewise_function_from_breakpoints_and_values(list(b)+[1], list(v)+[0])

    def n(self):
        """The maximum number of proper breakpoints of parameterized functions in the space."""
        return self._n

    def covers_space(self):
        """
        Not Implemented. Future use is intended to be a proof of correctness that all breakpoint sequences are covered.
        """
        ### TODO: Impl.
        raise NotImplementedError

    def refine_space(self):
        """
        Not Implemented. Future use is to prove that the current cell description for breakpoints is, with respect to inclusion, minimal or in the case that it is not minimal, find a minimal description.
        """
        ### TODO: Impl.
        raise NotImplementedError

    def write_data(self, output_file_name_style=None, max_rows=None):
        """
        Writes representative element data to a `.csv` file with one column and rows of representative elements.
        Optionally, write many `.csv` files with at `most max_rows` rows per file.
        Files are named output_file_file_name_style_filenumber.csv.
        The default output_file_name_style=f"Pi_Min_{n}".
        """
        # TODO: Future, support writing different types of data such as polyhedra data.
        if output_file_name_style is None:
            file_name_base = "Pi_Min_{}".format(self._n)
        else:
            file_name_base = output_file_name_style
        if max_rows is not None:
            assert(max_rows >= 1)
            num_files = len(self._data)//max_rows + 1
            file_name_base = file_name_base + "_part_0"
        if max_rows is None:
            max_rows = len(self._data)
            num_files = 1
        output_file = file_name_base + ".csv"
        for file_number in range(num_files):
            out_file = open(output_file, "w")
            data_writer = csv.writer(out_file, csv.QUOTE_NONE)
            for row in range(max_rows):
                try:
                    data_writer.writerow(self._data[max_rows * file_number + row])
                except IndexError:
                    break
            out_file.close()
            output_file = file_name_base[:-1]+"{}".format(file_number+1)+".csv"


def plot_2_d_polyhedral_complex(bkpt, highlight_index=-1, highlight_color='blue', **kwds):
    r"""
    Returns a plot of 2-d polyehdral complex of the breakpoint sequence of length n, Complex Delta P_bkpt.
    Lines generated from a partiuclar breakpoint can be highlighted.
    Assumes that the breakpoint sequence is correct and ``highlight_index`` is between 0 and ``len(bkpt)-1``, inclusive.

    EXAMPLES::
    sage: from cutgeneratingfunctiology.igp import *
    sage: plot_2_d_polyhedral_complex([0, 4/5]) # not tested

    A breakpoint can be highlighted::
    sage: plot_2_d_polyhedral_complex([0, 1/2, 4/5], highlight_index=1) # not tested

    An the color changed::
    sage: plot_2_d_polyhedral_complex([0, 4/5], highlight_index=1, highlight_color='green') # not tested
    """
    bkpt_ext = bkpt + [1]
    x = var('x')
    p = Graphics()
    n = len(bkpt)
    if highlight_index is None:
        highlight_index = -1
    for i in range(1,n+1):
        if i == highlight_index:
            p += line([(0,  bkpt_ext[i]), (bkpt_ext[i], 0)], color=highlight_color, linestyle='dashed')
        else:
            p += line([(0,  bkpt_ext[i]), (bkpt_ext[i], 0)], color='grey')
    for i in range(1,n):
        if i == highlight_index:
            p += line([(bkpt_ext[i], 1), (1, bkpt_ext[i])], color=highlight_color, linestyle='dashed')
        else:
            p += line([(bkpt_ext[i], 1), (1, bkpt_ext[i])], color='grey')
    for i in range(n+1):
        if i == highlight_index:
            p += plot(bkpt_ext[i], (0, 1),  color=highlight_color, linestyle='dashed')
        else:
            p += plot(bkpt_ext[i], (0, 1), color='grey')
    y = var('y')
    for i in range(n):
        if i == highlight_index:
            p += parametric_plot((bkpt_ext[i],y), (y,0,1), color=highlight_color, linestyle='dashed')
        else:
            p += parametric_plot((bkpt_ext[i],y), (y,0,1), color='grey')
    p += line([(1,0), (1,1)], color='grey')
    return p

def plot_2_d_polyhedral_complex_and_descendants(bkpt, max_row_length=5, max_number_additional_diagrams=15, **kwds):
    r"""
    Returns a plot of 2-d polyehdral complex of the breakpoint sequence of length n, Complex Delta P_bkpt
    together with plots of  Complex :math:`Delta P_{bkpt \cup \lambda^*}` where Complex Delta P_bkpt is a
    subcomplex of :math:`Delta P_{bkpt \cup \lambda^*}`.

    EXAMPLES::
    sage: from cutgeneratingfunctiology.igp import *
    sage: plot_2_d_polyhedral_complex_and_descendants([0, 1/2]) # not tested
    """
    n = len(bkpt)
    n_plus_one_breakpoints = unique_list(make_rep_bkpts_with_len_n(n+1, k=n, bkpts=[bkpt]))
    number_of_new_diagrams = len(n_plus_one_breakpoints)
    m = max_row_length
    i = 0
    plots = []
    while i < max_number_additional_diagrams and i < number_of_new_diagrams:
        plots.append((plot_2_d_polyhedral_complex(list(n_plus_one_breakpoints[i]), n, **kwds), (i % m , int(i/m) ,.75,.75)))
        i += 1
    if i < m:
        plots.append((plot_2_d_polyhedral_complex(bkpt), (int(i/2), -1,.75,.75)))
    else:
        plots.append((plot_2_d_polyhedral_complex(bkpt), (int(m/2), -1 ,.75,.75)))
    return multi_graphics(plots)

