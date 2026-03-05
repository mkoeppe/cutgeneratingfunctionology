from copy import deepcopy
from itertools import pairwise
from cutgeneratingfunctionology.igp import *
import csv
import os
from cutgeneratingfunctionology.spam.basic_semialgebraic import EmptyBSA
from cutgeneratingfunctionology.spam.basic_semialgebraic_pplite import BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron
import logging

minimal_function_cache_logger = logging.getLogger("cutgeneratingfunctionology.igp.minimal_function_cache")
minimal_function_cache_logger.setLevel(logging.INFO)

### Note to future reader: ###
### bkpt is assumed to be a breakpoint sequence of length n>= 2. 
### Breakpoint sequence are sorted lists of real numbers in [0,1). 
### A breakpoint sequences should always have 0 as an element.
### This is never strictly enforced in this file and it is assumed that
### the user is always provided a breakpoint sequence. 

# global defaults for logging from portions of igp; change to log different parts of igp.
log_paramateric_real_field = False
log_pw_functions = False


class RepElemGenFailure(Exception):
    pass


def nnc_poly_from_bkpt_sequence(bkpt, backend=None):
    r"""
    Defines an NNC polyhedron P such that for all b in P the delta complex of b is isomoprhic to the delta complex of bkpt.
    
    INPUT:
    - ``bkpt`` - sorted list of length 2 or more of sage type in the interval [0,1). 
    - ```backend`` - ``None``, ``str(pplite)`` 
    
    OUTPUT: class::``BasicSemialgebraicSet_veronese``
    
    EXAMPLES::
    sage: from cutgeneratingfunctionology.igp import * 
    sage: logging.disable(logging.INFO) # suppress logging for tests
    sage: nnc_poly_from_bkpt_sequence([0, 4/5])
    BasicSemialgebraicSet_veronese(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(Constraint_System {x0==0, -x1+1>0, 2*x1-1>0}, names=[x0, x1]), polynomial_map=[lambda0, lambda1])
    """
    n = len(bkpt)
    # assert(n >= 2)
    coord_names = []
    bkpt_vals = bkpt
    vals = bkpt_vals[0:n]
    bkpt_extd = list(bkpt)+[1]
    if not log_paramateric_real_field:
        parametric_logging_level = logging.getLogger("cutgeneratingfunctionology.igp.parametric").getEffectiveLevel()
        logging.getLogger("cutgeneratingfunctionology.igp.parametric").setLevel(logging.ERROR)
    for i in range(0,n):
        coord_names.append('lambda'+str(i))
    K = ParametricRealField(names=coord_names, values = vals, mutable_values=True, big_cells=True, default_backend=backend)
    K.gens()[0] == 0
    for i in range(n-1):
        K.gens()[i] < K.gens()[i+1]
    K.gens()[n-1] < 1
    for i in range(n):
        for j in range(n):
            if bkpt[i]+bkpt[j]>= 1:
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
    """
    Takes dim k-1 breakpoint NNC polyhedron (as a :class:`BasicSemialgebraicSet_base`), adds a dimension,
    and finds all possible representive elements for equivlance classes of polyhedral complexes. 
    
    INPUT: :class:`BasicSemialgebraicSet_Polyhedral`
    
    OUTPUT: unique list of breakpoint sequnes of lenght k (as tuples). 
    
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
    if k< 1:
        raise ValueError("bkpt_poly should have space dim at least 1.")
    model_bound_bkpts = [0]*k
    model_bound_bkpts[k-1] = 1
    # 0 < lambda_k <1
    # Using ppl and pplite type bsas have the same signature, so this should work regardless of backend.
    # this will fail if the BSA attached does not have these methods. 
    B_cap_N_b.add_linear_constraint(model_bound_bkpts, -1, operator.lt) # model bounds
    B_cap_N_b.add_linear_constraint(model_bound_bkpts, 0, operator.gt) # model bounds 
    bkpt_order = [0]*k
    bkpt_order[k-2] = 1
    bkpt_order[k-1] = -1
    B_cap_N_b.add_linear_constraint(bkpt_order, 0, operator.lt) # order on bkpts
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
                            B_cap_N_b_copy.add_linear_constraint(lhs_i, interval_w, interval_op)
                            lhs_i_plus_1 = [0]*k
                            lhs_i_plus_1[k-1] = -2
                            if i < k-1:
                                lhs_i_plus_1[i+1] = 1            
                                B_cap_N_b_copy.add_linear_constraint(lhs_i_plus_1, interval_w, operator.gt)
                            else:
                                B_cap_N_b_copy.add_linear_constraint(lhs_i_plus_1, interval_w + 1, operator.gt)
                            if not B_cap_N_b_copy.is_empty():
                                # does the line x+y equiv lambda_k mod 1 lie on/above/below (lambda_j,lambda_j)?
                                # modeled by  2lambda_j op lambda_k + w
                                lhs_j = [0]*k
                                lhs_j[j] = 2
                                lhs_j[k-1] = -1
                                B_cap_N_b_copy.add_linear_constraint(lhs_j, -line_w, line_op)
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
    - n,  integer, maximum length of breakpoint sequence.
    - k, assumed length of breakpoint sequences in ``bkpts``
    - bkpts, list of breakpoint sequenes (sorted    , length of every element, bkpts, an iterable of breakpoints all of length k.

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
    # Look into using a directed tree as an underlying data structure for generating elements.
    new_bkpts = []
    if n < 2:
        raise ValueError("n>=2")
    if k == n:
        raise ValueError("k<n")
    if k == 1 and bkpts is None:
        bkpts=[[0]]
    for bkpt in bkpts:
        new_bkpts += add_breakpoints_and_find_equiv_classes(nnc_poly_from_bkpt_sequence(bkpt).upstairs())
    new_bkpts = unique_list(new_bkpts)
    k += 1

    if k == n:
        minimal_function_cache_logger.info(f"Breakpoints of length {n} have been generated. ")
        return new_bkpts
    else:
        minimal_function_cache_logger.info(f"Breakpoints of lenght {k} have been generated. Now generating breakpoints of length {k+1}.")
        return make_rep_bkpts_with_len_n(n, k, new_bkpts, backend)


def generate_assumed_symmetric_vertices_continuous(fn, f, bkpt):
    """
    Silently assumes the symmetry condition holds for all vertices (x,y) in bkpt's breakpoints complex
    such that x+y equiv f.
    
    See fun:``generate_symmetric_vertices_continuous``.
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


def value_nnc_polyhedron_value_cords(bkpt, f_index, backend =None):
    """
    For a given breakpoints seqeunce and f index, write the value polyhedron as a basic semialgebraic set in only the value parameters.
    
    INPUT:
    - ``bkpt`` - sorted list of length 2 or more of sage type in the interval [0,1). 
    - ``f_index`` - integer between 1 and length of ``len(bkpt) -1``. 
    - ```backend`` - ``None``, ``str(pplite)`` 
    
    OUTPUT: 
    - class::``BasicSemialgebraicSet_veronese``
    
    EXAMPLES::
    
    sage: from cutgeneratingfunctionology.igp import * 
    sage: logging.disable(logging.INFO) # suppress logging for tests
    sage: value_nnc_polyhedron_value_cords([0,4/5], 1) # gmic with f=4/5
    BasicSemialgebraicSet_veronese(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(Constraint_System {x1-1==0, x0==0}, names=[x0, x1]), polynomial_map=[gamma0, gamma1])
    """
    # this saves a slight amount of overhead when detemrining points for the value polyhedron since the assumed
    # minimality test does not have to entierly go though parametric real field.
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
    K = ParametricRealField(names=coord_names, values = val, mutable_values=True, big_cells=True, allow_refinement=False, default_backend=backend)
    K.gens()[0] == 0
    for i in range(1, n):
        K.gens()[i] <=1
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

def value_nnc_polyhedron(bkpt, f_index, backend=None):
    """
    For a given breakpoints seqeunce and f index, write the value polyhedron as a basic semialgebraic set in the full space of parameters.
    
    INPUTS:
    - ``bkpt`` - sorted list of length 2 or more of sage type in the interval [0,1). 
    - ``f_index`` - integer between 1 and length of ``len(bkpt) -1``. 
    - ``backend`` - ``None``, ``str(pplite)`` 
    
    OUTPUT: 
    - class::``BasicSemialgebraicSet_veronese`` 
    
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
    K = ParametricRealField(names=coord_names, values = vals, mutable_values=True, big_cells=True, allow_refinement=False, default_backend=backend)
    # breakpoint parameters are the mesured breakpoint values. 
    for i in range(n):
        K.gens()[i] == bkpt[i]
    # necessary conditions on value parameters
    K.gens()[n] == 0 
    for i in range(1, n):
        K.gens()[i+n] <=1
        K.gens()[i+n] > 0
    h = piecewise_function_from_breakpoints_and_values(K.gens()[0:n]  + [1], K.gens()[n:2*n] + [0], merge=False)
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


def bsa_of_rep_element(bkpt, vals, backend=None):
    """
    Given pi_(bkpt, vals) is {minimal, not minimal}, find BSA subset of R^(2n) such that (bkpt, vals) in BSA and for all p
    in BSA, pi_p is {minimal, not minimal}.

    INPUT: 
    - ``bkpt`` - a breakpoint seqeunce
    - ``vals`` - list like of sage numerical types corrosponding values for the breakpoint sequence.
    - ``backend`` - None, ``str(pplite)``

    OUTPUT: A basic semialgebraic set.
    
    EXAMPLES::
    
    
    """
    n = len(bkpt)
    if not log_paramateric_real_field:
        parametric_logging_level = logging.getLogger("cutgeneratingfunctionology.igp.parametric").getEffectiveLevel()
        logging.getLogger("cutgeneratingfunctionology.igp.parametric").setLevel(logging.ERROR)
    if not log_pw_functions:
        pw_logging_level = logging.getLogger("cutgeneratingfunctionology.igp.functions").getEffectiveLevel()
        logging.getLogger("cutgeneratingfunctionology.igp.functions").setLevel(logging.ERROR)
    assert(n>=2)
    coord_names = []
    for i in range(n):
        coord_names.append('lambda'+str(i))
    for i in range(n):
        coord_names.append('gamma'+str(i))
    K = ParametricRealField(names=coord_names, values = bkpt+vals, big_cells=True, default_backend=backend)
    h = piecewise_function_from_breakpoints_and_values(K.gens()[0:n] + [1], K.gens()[n:2*n] + [0], merge=False)
    minimality_test(h)
    if not log_paramateric_real_field:
        logging.getLogger("cutgeneratingfunctionology.igp.parametric").setLevel(parametric_logging_level)
    if not log_pw_functions:
        logging.getLogger("cutgeneratingfunctionology.igp.functions").setLevel(pw_logging_level)
    return K.make_proof_cell().bsa


def find_minimal_function_reps_from_bkpts(bkpts, prove_minimality=True, backend=None):
    """
    Finds representative elements of minimal functions from a given breakpoint sequence.
    
    INPUT:
    - ``bkpts`` - an interable of breakpoint seqeunces.
    - ``prove_minimality`` - bool, proves minimality of paramaterized function
    - ``backend`` - None, ``str(pplite)``
    
    OUTPUT:
    - List of tuples of lists 
    
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
                    raise ValueError(f"({bkpt}, {test_val}) paramaterized by breakpoints and values is not a minimal function but assuming a breakpoint sequence is input, this should be minimal.")
            rep_elems.append((list(bkpt), test_val))
    if not log_pw_functions:
        logging.getLogger("cutgeneratingfunctionology.igp.functions").setLevel(pw_logging_level)
    return rep_elems


class BreakpointComplexClassContainer:
    """
    A container for the family of breakpoint complexes for peicewise linear functions
    with at most n breakpoints.
    
    The container assumes that loaded data is correct and performs no checking
    that the loaded data represnts the full space. 
    
    This class contains ways to read/write data for use with minimal function generation.
    
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
    
    A cell descrption of semialgebraic sets can be accessed::
    
    sage: all([isinstance(cell, BasicSemialgebraicSet_base) for cell in bkpt_of_len_2.get_nnc_poly_from_bkpt()])
    True
    
    
    """
    def __init__(self, n, **kwrds):
        self._n = n
        assert(self._n >= 2)
        if "backend" in kwrds.keys():
            self._backend = kwrds["backend"]
        else:
            self._backend = None
        if "load_rep_elem_data" in kwrds.keys():
            if kwrds[load_rep_elem_data] is None:
                minimal_function_cache_logger.info("Generating representative elements. This might take a while.")
                self._data = make_rep_bkpts_with_len_n(self._n, backend = self._backend)
            else:
                file_names = kwrds["load_bkpt_data"].split(",")
                self._data = []
                for file_name in file_names:
                    file = open(file_name, "r")
                    self._data += [eval(preparse(data)) for data in list(csv.reader(file))]
                    file.close()
                if "gen_elems_from_data" in kwrds.keys():
                    if kwrds[gen_elems_from_data] == True:
                        k = len(self._data[0])
                        if k < n:
                            self._data = make_rep_bkpts_with_len_n(n, k, self._data)
        else:
            minimal_function_cache_logger.info("Generating representative elements. This might take a while.")
            self._data = make_rep_bkpts_with_len_n(self._n,  backend=self._backend)

    def __repr__(self):
        return f"Container for the space breakpoint sequences of length {self._n} under equivlance of polyhedral complexes."

    def get_rep_elems(self):
        for bkpt in self._data:
            yield bkpt

    def get_nnc_poly_from_bkpt(self):
        for bkpt in self._data:
            yield nnc_poly_from_bkpt_sequence(bkpt, backend=self._backend)

    def num_rep_elems(self):
        return len(self._data)

    # def add_one_bkpt_to_all(self):
        # minimal_function_cache_logger.info("Generating representative elements. This might take a while.")
        # self._n = n+1
        # self._data = make_bkpts_with_len_n(self._n, self._n-1, self._data)
        
    def covers_space(self):
        ### TODO: This method should prove that container is correct.
        raise NotImplementedError

    def refine_space(self):
        """Ensures that repersentative elements are unique and contained in a single cell. """
        raise NotImplementedError
        ### TODO: Test this. I think I need to write  a containment method for the BSA 
        ### or pick the up the underlying polyhedron in the BSA. 
        self._data =  unique_list(self._data)
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
            file_name_base =output_file_name_style
        if max_rows is not None:
            assert(max_rows >= 1)
            num_files = len(self._data)//max_rows + 1
            file_name_base = file_name_base + "_part_0"
        if max_rows is None:
            max_rows = 0
        output_file = file_name_base +".csv"
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
    most n breakpoints paramaterized by breakpoints and values using semialgebraic sets.
    
    The container assumes that loaded data is correct and performs no checking.
    
    INPUTS:
    - n, an integer
    keywords:
    - backend, None or str(pplite)
    - load_bkpt_data, .csv file(s) of list of tuples of breakpoitns
    - load_rep_elem_data, .csv file(s) of list of tuples of representative elements of the space of minimal functions

    EXAMPLES::

    sage: from cutgeneratingfunctionology.igp import * 
    sage: logging.disable(logging.INFO) # suppress logging for tests
    sage: PiMin_2 = PiMinContContainer(2)
    sage: all([minimality_test(pi) for pi in PiMin_2.get_rep_functions()])
    True
    sage: PiMin_2
    Space of minimal functions with at most 2 breakpoints parameterized by breakpoints and values using semialgebraic sets.
    
    A cell descrption of semialgebraic sets can be accessed::
    
    sage: all([isinstance(cell, BasicSemialgebraicSet_base) for cell in PiMin_2.get_semialgebraic_sets()])
    True
    
    Data is stored as repersenative elements. The number of repersentative elements grows quickly. ::
    
    sage: PiMin_4 = PiMinContContainer(4)  # not tested
    sage: len([rep_elem for rep_elem in PiMin_4.get_rep_elems()]) # not tested
    987
    sage: len([rep_elem for rep_elem in PiMin_2.get_rep_elems()])
    3
    
    The container provides methods of writing data.::
    
    sage: PiMin_2.write_data() #
    
    Written data can be reused ::
    
    sage: PiMin_2_loaded_data = PiMinContContainer(2,  load_rep_elem_data="Pi_Min_2.csv")
    sage: len([rep_elem for rep_elem in PiMin_2_loaded_data.get_rep_elems()])
    3
    """
    def __init__(self, n, **kwrds):
        self._n = n
        assert(self._n >= 2)
        if "backend" in kwrds.keys():
            self._backend = kwrds["backend"]
        else:
            self._backend = None
        if "load_bkpt_data" in kwrds.keys() and "load_rep_elem_data" not in kwrds.keys():
            file_names = kwrds["load_bkpt_data"].split(",")
            bkpts = []
            for file_name in file_names:
                with open(file_name, newline='' ) as csvfile:
                    file_reader = csv.reader(csvfile)
                    for row in file_reader:
                        bkpts.append([eval(preparse(data)) for data in row])
            self._data = find_minimal_function_reps_from_bkpts(bkpts, backend = self._backend)
        elif "load_bkpt_data" not in kwrds.keys() and "load_rep_elem_data" in kwrds.keys():
            file_names = kwrds["load_rep_elem_data"].strip(" ").split(",")
            self._data = []
            for file_name in file_names:
                with open(file_name, newline='' ) as csvfile:
                    file_reader = csv.reader(csvfile)
                    for row in file_reader:
                        self._data.append([eval(preparse(data)) for data in row])
        else:
            minimal_function_cache_logger.info("Generating representative elements. This might take a while.")
            bkpts = make_rep_bkpts_with_len_n(self._n, backend=self._backend)
            minimal_function_cache_logger.info("Finished generating elements.")
            self._data = find_minimal_function_reps_from_bkpts(bkpts, backend=self._backend)

    def __repr__(self):
        return "Space of minimal functions with at most {} breakpoints parameterized by breakpoints and values using semialgebraic sets.".format(self._n)

    def get_semialgebraic_sets(self):
        """Iterator for semialgebraic set description"""
        for b, v in self._data:
            yield bsa_of_rep_element(list(b), list(v), backend=self._backend)

    def get_rep_elems(self):
        """Iterator for representative elements."""
        for b, v in self._data:
            yield (list(b), list(v))

    def get_rep_functions(self):
        """Iterator for representative functinos."""
        for b, v in self._data:
            yield piecewise_function_from_breakpoints_and_values(list(b)+[1], list(v)+[0])

    def n(self):
        """The maximum number of proper breakpoints of paramaterized functions in the space."""
        return self._n
        
    def covers_space(self):
        ### TODO: This method should prove that container is correct.
        raise NotImplementedError

    def refine_space(self):
        ### TODO: This method should be called when loading multiple file to ensure cells have not been duplicated.
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
            file_name_base =output_file_name_style
        if max_rows is not None:
            assert(max_rows >= 1)
            num_files = len(self._data)//max_rows + 1
            file_name_base = file_name_base + "_part_0"
        if max_rows is None:
            max_rows = len(self._data)
            num_files = 1
        output_file = file_name_base +".csv"
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

    