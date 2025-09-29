from copy import deepcopy
from itertools import pairwise
from cutgeneratingfunctionology.igp import *
import csv
import os


ppl_bsa = BasicSemialgebraicSet_veronese(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(0), polynomial_map=[], poly_ring=sym_ring, v_dict={})
pplite_bsa = BasicSemialgebraicSet_veronese(BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron(0), polynomial_map=[], poly_ring=sym_ring, v_dict={})

def mod_one(x):
    if x >= 0:
        return x - int(x)
    return x - int(x) + 1

def add_breakpoint(bkpt):
    """
    Given a breakpoint sequence, creates a list of breakpoint sequences of length one more such that
    each the breakpoint complex of each breakpoint sequence is not guaranteed to be isomorphic to any
    other breakpoint sequence in the list. 


    INPUT: Assume vector or list of n breakpoints with lambda_0=0 and lambda_i<lambda_i+1<1 (i.e. a breakpoint sequence)
    OUTPUT: List of vectors each of which is a breakpoint sequence with $n+1$ elements and no two polyhedral complexes based on the data from each breakpoint sequence is in the same isomoprhism class.

    TESTS::
    >>> add_breakpoint([0])
    [[0, 3/4], [0, 1/2], [0, 1/4]]
    >>> add_breakpoint([0,1/3])
    [[0, 1/12, 1/3], [0, 1/6, 1/3], [0, 1/4, 1/3], [0, 1/3, 2/3], [0, 1/3, 5/12], [0, 1/3, 1/2], [0, 1/3, 7/12], [0, 1/3, 5/6]]
    """
    possible_new_seqs = []
    possible_eq_lambda_star = [b/2 for b in bkpt+[1]] + [(1+b)/2 for b in bkpt]
    intervals_along_y_eq_x = sorted(list(set(tuple(bkpt + [1] + [b/2 for b in bkpt+[1]] + [(1+b)/2 for b in bkpt]))))
    possible_ne_lambda_star = []
    for lower_bound, upper_bound in pairwise(intervals_along_y_eq_x):
        possible_ne_lambda_star.append(1/2 * (lower_bound + upper_bound))
    possible_lambda_star = possible_eq_lambda_star + possible_ne_lambda_star
    for lambda_star in possible_eq_lambda_star + possible_ne_lambda_star:
        temp_bkpt = deepcopy(bkpt)
        temp_bkpt.append(mod_one(lambda_star))
        possible_new_seqs.append(sorted(temp_bkpt))
    possible_new_seqs = [list(y) for y in set([tuple(x) for x in possible_new_seqs]) if len(set(y)) == len(y)]
    return possible_new_seqs



def make_bkpts_with_len_n(n, k=1, bkpts=None):
    """
    Produce representative elements of every isomorphism class of breakpoints complexes for breakpoint sequences of length n.

    Note, this function does not check that the data input is correct and assume it is being used correctly. 

    INPUT: n, length of breakpoint sequence, k, length of every element, bkpts, an iterable of breakpoints all of length k.

    OUTPUT: A list of representative elements of every isomorphism class of breakpoints complexes for breakpoint sequences of length n extrapolated from bkpts.
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
        new_bkpts += add_breakpoint(bkpt)
    new_bkpts = [list(y) for y in set([tuple(x) for x in new_bkpts])]
    k += 1
    if k == n:
        return new_bkpts
    else:
        return make_bkpts_with_len_n(n, k, new_bkpts)


def nnc_poly_from_bkpt(bkpt, backend=None):
    n = len(bkpt)
    if backend is None:
        bsa = ppl_bsa.copy()
    if backend == "pplite":
        bsa = pplite_bsa.copy()
    assert(n >= 2)
    coord_names = []
    bkpt_vals = bkpt
    vals = bkpt_vals[0:n]
    for i in range(0,n):
        coord_names.append('lambda'+str(i))
    K = ParametricRealField(names=coord_names, values = vals, mutable_values=True, big_cells=True, partial_test_point_mode=True, bsa=bsa)    
    logging.disable(logging.INFO)
    K.gens()[0] == 0
    for i in range(n-1):
        K.gens()[i] < K.gens()[i+1]
    K.gens()[n-1] < 1    
    h = piecewise_function_from_breakpoints_and_values([0] + K.gens()[0:n-1] + [1], [0]*(n+1), merge=False)
    for vert in generate_type_1_vertices_continuous(h, operator.ge, [0] + K.gens()[0:n-1] + [1]):
        vert
    for vert in generate_type_2_vertices_continuous(h, operator.ge, [0] + K.gens()[0:n-1] + [1]):
        vert
    return K.make_proof_cell().bsa


class BreakpointComplexClassContainer:
    """
    A container for the family of breakpoint complexes for peicewise linear functions
    with at most n breakpoints.
    """
    def __init__(self, n, **kwrds):
        self._n = n
        assert(self._n >= 2)
        if "backend" in kwrds.keys():
            if kwrds[backend] == "pplite":
                self._backend = "pplite"
            else:
                self._backend = None
        if "load_rep_elem_data" in kwrds.keys():
            if kwrds[load_rep_elem_data] is None:
                logging.warning("Generating representative elements. This might take a while.")
                self._data = make_bkpts_with_len_n(self._n)
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
                            self._data = make_bkpts_with_len_n(n, k, self._data)
        else:
            logging.warning("Generating representative elements. This might take a while.")
            self._data = make_bkpts_with_len_n(self._n)

    def __repr__(self):
        return "Container of a family of breakpoint"

    def get_rep_elems(self):
        for bkpt in self._data:
            yield bkpt

    def get_nnc_poly_from_bkpt(self):
        for bkpt in self._data:
            yield nnc_poly_from_bkpt(bkpt, self._backend)

    def num_rep_elems(self):
        return len(self._data)

    def add_one_bkpt_to_all(self):
        logging.warning("Generating representative elements. This might take a while.")
        self._data = make_bkpts_with_len_n(self._n+1, self._n, self._data)

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
            num_files = len(self._data)//max_rows + 1
            file_name_base = file_name_base + "_part_0"
        if max_rows is None:
            max_rows = 0
        output_file = file_name_base +".csv"
        for file_number in range(num_files):
            out_file = open(output_file, "w")
            data_writer = csv.writer(out_file, csv.QUOTE_NONE)
            for row in range(max_row):
                data_writer.writerow(self._data[max_row * file_number + row])
            out_file.close()
            output_file = file_name_base[:-1]+"{}".format(file_number)


def generate_assumed_symmetric_vertices_continuous(fn, f, bkpt):
    """
    Assumes the symmetry condition holds for all vertices (x,y) in bkpt's breakpoints complex
    such that x+y equiv f.
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


def assume_minimality(bkpt, f_index, backend=None):
    """
    Given a breakpoint sequence, bkpt, and an index for f, f_index, determine if there is a (rep_bkpt, v)
    such that pi_(rep_bkpt,v) is minimal, pi_(bkpt,v)(lambda_f_index)=1, and rep_bkpt's breakpoint complex
    is isomorphic to bkpt's breakpoint complex.

    INPUT: bkpt a list or vector of length n. bkpt is assumed to be breakpoint sequence, f_index an integer.

    OUTPUT: (rep_bkpt, v), a pair of lists of length n with the described property.
    """
    n = len(bkpt)
    if backend is None:
        bsa = ppl_bsa.copy()
    if backend == "pplite":
        bsa = pplite_bsa.copy()
    assert(n >= 2)
    assert(f_index >= 1)
    assert(f_index <= n - 1)
    coord_names = []
    bkpt_vals = bkpt
    vals = bkpt_vals[1:n]+ [None]*(n-1)
    for i in range(1,n):
        coord_names.append('lambda'+str(i))
    for i in range(1,n):
        coord_names.append('gamma'+str(i))
    logging.disable(logging.INFO)
    K = ParametricRealField(names=coord_names, values = vals, mutable_values=True, big_cells=True, partial_test_point_mode=True, bsa=bsa)
    for i in range(n-1):
        K.gens()[i+n-1] <=1
        K.gens()[i+n-1] > 0 
    h = piecewise_function_from_breakpoints_and_values([0] + K.gens()[0:n-1] + [1], [0] + K.gens()[n-1:2*n-2] + [0], merge=False)
    for vert in generate_type_1_vertices_continuous(h, operator.ge, [0] + K.gens()[0:n-1] + [1]):
        vert
    for vert in generate_type_2_vertices_continuous(h, operator.ge, [0] + K.gens()[0:n-1] + [1]):
        vert
    for vert in generate_assumed_symmetric_vertices_continuous(h, K.gens()[f_index-1], [0] + K.gens()[0:n-1] + [1]):
        vert
    for i in range(n):
        K.gens()[i] == bkpt[i]
    K.find_test_point()
    h_2 = piecewise_function_from_breakpoints_and_values([0]+list(K._values[0:n-1])+[1], [0] + list(K._values[n-1:2*n-2])+[0])
    is_minimal = minimality_test(h_2)
    if is_minimal:
        rep_bkpt = [0] + list(K._values[0:n-1])
        v = [0] + list(K._values[n-1:2*n-2])
        return (rep_bkpt, v)



def bsa_of_rep_element(bkpt, vals):
    """ 
    Given pi_(bkpt, vals) is {minimal, not minimal}, find BSA subset of R^(2n) such that (bkpt, vals) in BSA and for all p
    in BSA, pi_p is {minimal, not minimal}. 

    INPUT: (bkpt, vals) are lists or vectors of length n and bkpt is a proper breakpoints sequence and vals
    is the corresponding value parameters.

    OUTPUT: A basic semialgebraic set.
    """
    n = len(bkpt)
    assert(n>=2)
    coord_names = []
    for i in range(0,n):
        coord_names.append('lambda'+str(i))
    for i in range(0,n):
        coord_names.append('gamma'+str(i))
    logging.disable(logging.INFO)
    K = ParametricRealField(names=coord_names, values = bkpt+vals, big_cells=True)
    h = piecewise_function_from_breakpoints_and_values([0] + K.gens()[1:n] + [1], [0] + K.gens()[n+1:2*n] + [0], merge=False)
    minimality_test(h)
    return K.make_proof_cell().bsa


def find_minimal_function_reps_from_bkpts(bkpts, backend=None):
    """
    Finds representative elements of minimal functions if they exist from the given breakpoint sequence.  
    """
    data = []
    for bkpt in bkpts:
        for i in range(1, len(bkpt)):
            result = assume_minimality(bkpt, i, backend)
            if result is not None:
                data.append(result)
    return data

class PiMinContContainer:
    """
    A container for the space of continuous piecewise linear minimal functions with at
    most n breakpoints paramaterized by breakpoints and values using semialgebraic sets.

    TESTS::

    >>> PiMin_at_most_4_breakpoints = PiMinContContainer(4)
    >>> all([minimality_test(pi) for PiMin_at_most_4_breakpoints.get_rep_functions()])
    True
    >>> 
    """
    def __init__(self, n, **kwrds):
        self._n = n
        assert(self._n >= 2)
        if "backend" in kwrds.keys():
            if kwrds[backend] == "pplite":
                self._backend = "pplite"
            else:
                self._backend = None
        if "load_bkpt_data" in kwrds.keys() and "load_rep_elem_data" not in kwrds.keys():
            file_names = kwrds["load_bkpt_data"].split(",")
            bkpts = []
            for file_name in file_names:
                file = open(file_name, "r")
                bkpts += [eval(preparse(data)) for data in list(csv.reader(file))]
                close(file)
            self._data = find_minimal_function_reps_from_bkpts(bkpts)
        elif "load_bkpt_data" not in kwrds.keys() and "load_rep_elem_data" in kwrds.keys():
            file_names = kwrds["load_rep_elem_data"].split(",")
            self._data = []
            for file_name in file_names:
                file = open(file_name, "r")
                self._data += [(eval(preparse(data[0])), eval(preparse(data[1]))) for data in list(csv.reader(file))]
                file.close()
        else:
            logging.warning("Generating representative elements. This might take a while.")
            bkpts = make_bkpts_with_len_n(self._n)
            self._data = find_minimal_function_reps_from_bkpts(bkpts, self._backend)

    def __repr__(self): 
        return "Space of minimal functions with at most {} breakpoints parameterized by breakpoints and values using semialgebraic sets.".format(self._n)

    def get_semialgebraic_sets(self):
        for b, v in self._data:
            yield bsa_of_rep_element(b, v) 

    def get_rep_elems(self):
        for b, v in self._data:
            yield (b, v)

    def get_rep_functions(self):
        for b, v in self._data:
            yield piecewise_function_from_breakpoints_and_values(list(b)+[1], list(v)+[0])

    def n(self):
        return self._n

    def covers_space(self):
        raise NotImplementedError

    def refine_space(self):
        raise NotImplementedError

    def write_data(self, output_file_name_style=None, max_rows=None):
        """
        Writes representative element data to a `.csv` file with one column and rows of representative elements.
        Optionally, write many `.csv` files with at `most max_rows` rows per file. 
        Files are named output_file_file_name_style_filenumber.csv.
        The default output_file_name_style="Pi_Min_n".
        """
        # TODO: Future, support writing different types of data such as polyhedra data. 
        if output_file_name_style is None:
            file_name_base = "Pi_Min_{}".format(self._n)
        else:
            file_name_base =output_file_name_style
        if max_rows is not None:
            num_files = len(self._data)//max_rows + 1
            file_name_base = file_name_base + "_part_0"
        if max_rows is None:
            max_rows = 0
        output_file = file_name_base +".csv"
        for file_number in range(num_files):
            out_file = open(output_file, "w")
            data_writer = csv.writer(out_file, csv.QUOTE_NONE)
            for row in range(max_row):
                data_writer.writerow(self._data[max_row * file_number + row])
            out_file.close()
            output_file = file_name_base[:-1]+"{}".format(file_number)


