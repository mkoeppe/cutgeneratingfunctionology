if '' not in sys.path:
    sys.path = [''] + sys.path
from igp import *
from sage.numerical.mip import MIPSolverException
from sage.misc.all import cached_method
from collections import defaultdict
from itertools import product
import sage.numerical.backends.glpk_backend as backend

current_grid_setting = "default"

available_settings = ["default", "one_bucket"]

def construct_2d_strip_by_FastPiecewise(f):
    """
    （Credit by Shuidie Yao)
    Take a FastPiecewise function, convert it to straight 2d.
    Points are extended to lines, intervals are converted to rectangles
    
    EXAMPLES::
    
        sage: logging.disable(logging.INFO)
        sage: h = gmic()
        sage: for p in construct_2d_strip_by_FastPiecewise(h):
        ....:     p.vertices()
        (A vertex at (4/5, 0), A vertex at (4/5, 1))
        (A vertex at (1, 0),
         A vertex at (1, 1),
         A vertex at (4/5, 0),
         A vertex at (4/5, 1))
        (A vertex at (0, 0),
         A vertex at (0, 1),
         A vertex at (4/5, 0),
         A vertex at (4/5, 1))
        (A vertex at (0, 0), A vertex at (0, 1))
        (A vertex at (1, 0), A vertex at (1, 1))
    """
    polyhedra = set([Polyhedron(vertices = [[intv[0],0],[intv[1],0],[intv[0],1],[intv[1],1]]) for intv in f.intervals()])
    for intv in f.intervals():
        if len(intv) == 2:
            polyhedra.add(Polyhedron(vertices = [[intv[0],0],[intv[0],1]]))
            polyhedra.add(Polyhedron(vertices = [[intv[1],0],[intv[1],1]]))
        else:
            if intv[0] != intv[1] and intv.left_closed:
                polyhedra.add(Polyhedron(vertices = [[intv[0],0],[intv[0],1]]))
            if intv[0] != intv[1] and intv.right_closed:
                polyhedra.add(Polyhedron(vertices = [[intv[1],0],[intv[1],1]]))
    return polyhedra

# This method is designed to change the grid setting globally.
# Locally, users can specify the grid setting when
# initiate a :class:`Grid`.
def grid_setting(new_setting=None):
    r"""
    Set or get the current setting of grid.
    """
    global current_grid_setting
    if new_setting is not None:
        if new_setting not in available_settings:
            raise ValueError("Style must be one of: {}".format(
                             ", ".join(available_settings)))
        current_grid_setting = new_setting
    return current_grid_setting

def ideal_width(number_of_polyhedra, ambient_dim, side_length=1, constant=1, convergent_index=20):
    r"""
    This method is used to determine the "ideal" width of each bucket in a grid by using
    the formula constant / (number_of_polyhedra ^ (1 / ambient_dim)).

    FIXME: need experiments to determine ``constant``
    Koeppe's comments on width:
    1) In the end, only computational experiment can determine what is the
    most efficient choice of width.
    2) For small objects (like points), one would expect something like
    constant / (number_of_polyhedra ^ (1 / ambient_dim)) (currently commented-out) 
    to be best. Again, only computational experiments can determine what
    constant to use.

    INPUT::

    - ``number_of_polyhedra`` -- an integer indicating the numbers of polyhedra in the grid
    
    - ``ambient_dim`` -- an integer indicating the ambient dimension of the polyhedra

    - ``side_length`` -- (default: 1) the length of the hypercube as ``self``
    
    - ``constant`` -- (default: 1) the constant in the formula of calculating the width of each bucket

    - ``convergent_index`` -- (default: 20) the parameter to :func:`convergent`, a method to approximate
      the width by continued fraction

    OUTPUT::

    - a rational number indicating the width of each bucket of the grid
    """
    width = constant / (number_of_polyhedra ^ (1 / ambient_dim))
    cf = continued_fraction(width)
    cf_conv = cf.convergent(convergent_index)
    d = cf_conv.denominator()
    n = cf_conv.numerator()
    return side_length / ceil(d / n)

###############################
# Polyhedral Arrangement
###############################

class PolyhedralArrangement(object):
    """
    Define a Polyhedral Arrangement.

    INPUT:

    - ``collection`` -- either a list or a tuple or a set of polyhedra, 
    or a dictionary of polyhedra.
    The keys of the dictionary are the dimensions of the polyhedra

    - **kwds -- keywords to pass on self._grid
 
    Attributes:

    - ``ambient_dim`` -- the ambient dim of each polyhedron
        
    - ``collection`` -- a dictionary of the polyhedra (not including sub-polyhedra).
      whose keys are the dimensions of the polyhedra
    
    - ``grid`` -- a :class:`Grid` that can be used for discrete binary search
      on polyhedra

    - ``grid_contents`` -- a dictionary whose keys are buckets and whose values are
      sets of polyhedra that intersect with the keys    
        
    - ``number_of_polyhedra`` -- the numbers of polyhedra (not including sub-polyhedra)
      when set up the arrangement
        
    - ``polyhedron_buckets`` -- a dictionary whose keys are elements and whose values are
      lists of buckets that intersect with the elements. This attribute is like an inverse of
      ``grid_contents``

    EXAMPLES::

        sage: p = Polyhedron(vertices=[[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]])
        sage: q = Polyhedron(vertices=[[10,30,50,60,70]])
        sage: collection = (p, q)
        sage: PA = PolyhedralArrangement(collection)
        sage: PA
        A polyhedral arrangement of 2 polyhedra in ambient dimension 5.
        It is arranged on a 5-dimensional grid with side length 70.
        The grid has 2^5 buckets. Each bucket has width 35.
        The lower left coordinate of the grid is [0, 0, 0, 0, 0]
    """
    def __init__(self, collection, **kwds):
        if isinstance(collection, (list, tuple, set)):
            coll_dict = {}
            for p in collection:
                d = p.dim()
                if d >= 0:
                    if d in coll_dict:
                        coll_dict[d].append(p)
                    else:
                        coll_dict[d] = [p]
        elif isinstance(collection, dict):
            coll_dict = copy(collection)
        else:
            raise ValueError
        if not coll_dict:
            return
        self._ambient_dim = coll_dict.values()[0][0].ambient_dim()
        self._number_of_polyhedra = sum([len(v) for v in coll_dict.values()])
        self._collection = coll_dict
        self._grid_setting = current_grid_setting
        if current_grid_setting == "one_bucket":
            self._grid = Grid(self._number_of_polyhedra, self._ambient_dim, grid_setting="one_bucket")
        else:
            grid_info = self._grid_info()
            if kwds is not None:
                for key, value in kwds.iteritems():
                    if key == "lower_left_coordinates":
                        grid_info[0] = value
                    elif key == "side_length":
                        grid_info[1] = value
                    elif key == "bucket_width":
                        grid_info[2] = value
                    elif key == "grid_setting":
                        self._grid_setting = value
            self._grid = Grid(self._number_of_polyhedra, self._ambient_dim,
                              lower_left_coordinates=grid_info[0],
                              side_length=grid_info[1],
                              bucket_width=grid_info[2], 
                              grid_setting=self._grid_setting)
        if self._grid_setting == "default":
            # if grid_setting is "one_bucket", then no dictionaries based
            # on grid is constructed
            self._grid_contents = self._grid._contents
            self._polyhedron_buckets = defaultdict(set)

    def __contains__(self, p):
        r"""
        Check ``p`` is an element of ``self``. 
        """
        if p in self.elements():
            return True
        else:
            return False

    def _grid_info(self):
        r"""
        Return a list of the default grid info on.

        INPUT::

        - ``poly`` -- a list, a tuple, or a set of polyhera that are
        in the grid we want to construct

        OUTPUT::

        - a list the first element of the tuple is the lower left integer
          coordinates of the bounding box of ``self``,
          the second element is the side length of the bounding box,
          the third element is the bucket width of the grid, default is "ideal" (see
          docstrings in :class:``Grid`` for more information)

        EXAMPLES::

            sage: p = Polyhedron(vertices=[[-1, 2], [2, -1]])
            sage: q = Polyhedron(vertices=[[2, 2], [6, 7]])
            sage: pa = PolyhedralArrangement([p, q])
            sage: pa._grid_info()
            [[-1, -1], 8, 'ideal']
        """
        elements = self.elements()
        dim = self.ambient_dim()
        ll_coord = [0 for i in range(dim)]
        ur_coord = [0 for i in range(dim)]
        for p in elements:
            min_max = _min_max_int_coordinates(p)
            ll_coord = [min_max[i][0] if min_max[i][0] <= ll_coord[i] else ll_coord[i] for i in range(dim)]
            ur_coord = [min_max[i][1] if min_max[i][1] >= ur_coord[i] else ur_coord [i] for i in range(dim)]
        return list((ll_coord, max([ur - ll for ur, ll in zip(ur_coord, ll_coord)]), "ideal"))

    def __iter__(self):
        for polyhedron in self.elements():
            yield polyhedron

    def __repr__(self):
        desc = 'A polyhedral arrangement of '
        desc += repr(self._number_of_polyhedra)
        desc += ' polyhedra in ambient dimension '
        desc += repr(self._ambient_dim)
        desc += '.\n'
        if self.grid_setting() == "one_bucket":
            desc += repr(self._grid)
        else:
            desc += 'It is arranged on a' + repr(self._grid)[1:]
        return desc

    def _add_constraints(self, lp, p, x):
        r""""
        Add the equations and inequalities of ``p`` as new equations and inequalities
        to ``lp`` with variables ``x``.
        """
        def linear_constraint(cstr, cstr_type):
                """
                Construct linear constraint (can be equation or inequality depends
                on constraint type 'cstr_type') by the given constraint 'cstr'
                """
                s = sum(x[index] * coef for index, coef in enumerate(cstr[1:]))
                if cstr_type == "==":
                    return s == -cstr[0]
                else:
                    return s >= -cstr[0]
        for eqn in p.equations_list():
            f = linear_constraint(eqn, "==")
            lp.add_constraint(f)
        for ieq in p.inequalities_list():
            g = linear_constraint(ieq, ">=")
            lp.add_constraint(g)
        return lp

    def _set_up_linear_programming(self, p):
        r"""
        Set up a linear programming for a polyhedron ``p``.
        """
        lp = MixedIntegerLinearProgram(solver="GLPK")
        x = lp.new_variable()
        lp.set_objective(None)
        self._add_constraints(lp, p, x)
        return lp, x

    def _set_up_linear_programming_for_intersection(self, p, q):
        r"""
        Construct a LP to solve if polyhedron ``p`` and polyhedron ``q`` intersect

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]], eqns = [[1,-1,-1,-1,-1]])
            sage: q = Polyhedron(eqns = [[10,30,50,60,70]])
            sage: PA = PolyhedralArrangement((p, q))
            sage: lp = PA._set_up_linear_programming_for_intersection(p, q)
            sage: lp.show()
            Maximization:
            <BLANKLINE>
            <BLANKLINE>
            Constraints:
              1.0 <= x_0 + x_1 + x_2 + x_3 <= 1.0
              - x_0 <= 0.0
              x_0 + x_1 + x_2 <= 1.0
              - x_2 <= 0.0
              - x_1 <= 0.0
              -1.0 <= 3.0 x_0 + 5.0 x_1 + 6.0 x_2 + 7.0 x_3 <= -1.0
            Variables:
              x_0 is a continuous variable (min=-oo, max=+oo)
              x_1 is a continuous variable (min=-oo, max=+oo)
              x_2 is a continuous variable (min=-oo, max=+oo)
              x_3 is a continuous variable (min=-oo, max=+oo)
            sage: # compare with the equations and inequalities from the polyhedra
            sage: p.equations(); q.equations(); p.inequalities(); q.inequalities()
            (An equation (1, 1, 1, 1) x - 1 == 0,)
            (An equation (3, 5, 6, 7) x + 1 == 0,)
            (An inequality (1, 0, 0, 0) x + 0 >= 0,
             An inequality (-1, -1, -1, 0) x + 1 >= 0,
             An inequality (0, 0, 1, 0) x + 0 >= 0,
             An inequality (0, 1, 0, 0) x + 0 >= 0)
            ()
        """
        if p.dim() == -1 or q.dim() == -1:
            raise ValueError("Cannot intersect with empty polyhedron")
        lp, x = self._set_up_linear_programming(p)
        self._add_constraints(lp, q, x)
        return lp

    def _small_perturbation(self, l):
        r"""
        Add small perturbations to elements of the list ``l``.
        """
        n = len(l)
        non_zero_elements = [s for s in l if s != 0]
        if len(non_zero_elements) == 0:
            minimum = 1
        else:
            minimum = min(non_zero_elements)
        pert = [randint(1, 10^3)/10^8 * minimum for i in range(n)]
        return [l[i] + pert[i] for i in range(n)]

    def _solve_linear_programming(self, lp):
        try:
            lp.solve()
            return True
        except MIPSolverException:
            return False
        else:
            print "Unexpected error"
            raise

    def _update_variables_bounds(self, lower_left, upper_right, lp, x):
        r"""
        Update the bounds of variables ``x`` of a linear programming
        problem ``lp`` by the new lower left ``lower_left`` and
        the upper right ``upper_right`` buckets for the bounding box
        of ``x``.
        """
        var_min, offset = self.grid().variables_min_offset(lower_left, upper_right)
        for index, minimum in enumerate(var_min):
            lp.set_min(x[index], minimum)
            lp.set_max(x[index], minimum + offset[index])
        lp.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
        lp.solver_parameter("primal_v_dual", "GLP_DUAL")
        return lp, x

    def ambient_dim(self):
        r"""
        The ambient dim of each polyhedron.
        """
        return self._ambient_dim

    def bounding_box_in_grid(self, p):
        r"""
        Return a set of buckets that form a bounding box of a
        polyhedron ``p`` in the grid.
        """
        g = self.grid()
        box = bounding_box(p, False)
        vertices_as_buckets = [g.point_in_bucket(v) for v in box.vertices_list()]
        buckets_vertices_matrix = matrix(vertices_as_buckets)
        min_max_int_coords = [(min(c), max(c)) for c in buckets_vertices_matrix.columns()]
        return set(product(*[[i[0], i[1]] for i in min_max_int_coords]))

    @cached_method
    def elements(self):
        r"""
        A list of maximal faces of the ``self``.
        """
        collection = set()
        coll_dict = self.collection_dict()
        for polyhedra_list in coll_dict.values():
            for polyhedron in polyhedra_list:
                collection.add(polyhedron)
        return list(collection)
 
    def collection_dict(self):
        r"""
        Return a dictionary of the polyhedra (not including sub-polyhedra).
        whose keys are the dimensions of the polyhedra
        """
        return self._collection

    @cached_method
    def faces_of_polyhedra(self):
        r"""
        Construct a dictionary whose keys are the elements of ``self``
        and the valeus are a set of faces of the elements.
        """
        faces_of_poly = defaultdict(set)
        for c in self:
            n = c.dim()
            for i in range(n+1):
                for f in c.faces(i):
                    faces_of_poly[c].add(f.as_polyhedron())
        return faces_of_poly

    def grid(self):
        r"""
        A :class:``Grid`` that is used for grid method
        on the elements of ``self``.
        """
        return self._grid

    def grid_contents(self):
        r"""
        Return a dictionary whose keys are buckets and whose values are
        sets of polyhedra that intersect with the keys.
        """
        if self.grid_setting() == "one_bucket":
            return self.grid()
        self.update_dictionaries()
        return self._grid_contents

    def grid_setting(self):
        return self._grid_setting

    def in_buckets(self, p):
        r"""
        Return a set of buckets that the element ``p`` intersect with.

        EXAMPLES::

            sage: logging.disable(logging.INFO)
            sage: p = Polyhedron(vertices=[[8/21, 1/4], [4/5, 1/5], [5/6, 4/5]])
            sage: q = Polyhedron(vertices=[[1/2, 1/2], [5/6, 1/2], [1/2, 1]])
            sage: x = Polyhedron(vertices=[[4/5, 4/5]])
            sage: pa = PolyhedralArrangement([p, q, x], lower_left_coordinates=None, side_length=1, bucket_width=None)
            sage: pa.in_buckets(p)
            {(1, 0), (1, 1), (2, 0), (2, 1), (2, 2)}
            sage: h = gmic()
            sage: poly = construct_2d_strip_by_FastPiecewise(h)
            sage: pa = PolyhedralArrangement(poly, lower_left_coordinates=None, side_length=1, bucket_width="ideal")
            sage: collection = pa.elements()
            sage: for p in collection:
            ....:         p.vertices()
            ....:         pa.in_buckets(p)
            (A vertex at (0, 0), A vertex at (0, 1))
            {(0, 0), (0, 1), (0, 2)}
            (A vertex at (0, 0),
             A vertex at (0, 1),
             A vertex at (4/5, 0),
             A vertex at (4/5, 1))
            {(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)}
            (A vertex at (1, 0),
             A vertex at (1, 1),
             A vertex at (4/5, 0),
             A vertex at (4/5, 1))
            {(2, 0), (2, 1), (2, 2)}
            (A vertex at (4/5, 0), A vertex at (4/5, 1))
            {(2, 0), (2, 1), (2, 2)}
            (A vertex at (1, 0), A vertex at (1, 1))
            {(2, 0), (2, 1), (2, 2)}
        """
        if self.grid_setting() == "one_bucket":
            return self.grid()
        return self.polyhedron_buckets()[p]

    def intersecting_pairs(self):
        r"""
        A set of pairs of intersecting elements of ``self``.
        """
        elements = self.elements()
        l = len(elements)
        pairs = set()
        for i in range(l):
            p = elements[i]
            for j in range(i+1, l):
                q = elements[j]
                if self.grid_setting() == "one_bucket":
                    if dim(p.intersection(q)) >= 0:
                        pairs.add((p, q))
                else:
                    if self.intersects(p, q):
                        pairs.add((p, q))
        return pairs

    @cached_method
    def intersects(self, p, q):
        r"""
        If ``grid_setting`` is ``deafult``, use grid method to check if
        two polyhedra ``p`` and ``q`` intersect. 
        
        If ``grid_setting`` is ``one_bucket``, use the dimension of the
        intsection of ``p`` and ``q to check if they intersect.
        
        Return True if they intersect, else return False.

        EXAMPLES::

            sage: a = Polyhedron(vertices=[[19/32, 14/16], [3/4, 11/16]])
            sage: b = Polyhedron(vertices=[[11/16, 6/16], [11/16, 11/16]])
            sage: c = Polyhedron(vertices=[[9/16, 1/8], [10/16, 1/8], [13/16, 14/16], [14/16, 14/16]])
            sage: pa = PolyhedralArrangement([a, b, c])
            sage: a.intersection(b)
            The empty polyhedron in QQ^2
            sage: pa.intersects(a, b)
            False
            sage: b.intersection(c)
            A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices
            sage: pa.intersects(b, c)
            True
            sage: a.intersection(c)
            A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex
            sage: pa.intersects(a, c)
            True

            sage: p = Polyhedron(vertices=[[1, 1/2, 0], [0, 1/2, 0], [1, 1/2, 1], [0, 1/2, 1]])
            sage: q = Polyhedron(vertices=[[1, 0, 0], [0, 1, 1]])
            sage: x = Polyhedron(vertices=[[1/2, 1/2, 1], [1/2, 1, 1]])
            sage: pa = PolyhedralArrangement([p, q, x])
            sage: p.intersection(q)
            A 0-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex
            sage: pa.intersects(p, q)
            True
            sage: p.intersection(x)
            A 0-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex
            sage: pa.intersects(p, x)
            True
            sage: q.intersection(x)
            The empty polyhedron in QQ^3
            sage: pa.intersects(q, x)
            False
        """
        if self.grid_setting() == "one_bucket":
            return dim(p.intersection(q)) >= 0
        p_buckets = self.in_buckets(p)
        q_buckets = self.in_buckets(q)
        # compare two sets of buckets that p and q belong to 
        if p_buckets.intersection(q_buckets) == set():
            return False
        if self.intersects_by_linear_programming(p, q):
            return True
        else:
            return False

    def intersects_by_linear_programming(self, p, q):
        r"""
        Use LP to check if two polyhedra intersect.
        Return True if polyhedron ``p`` and polyhedron ``q`` intersect.
        Return False it they do not intersect.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]], eqns = [[1,-1,-1,-1,-1]])
            sage: q = Polyhedron(eqns = [[10,30,50,60,70]])
            sage: PA = PolyhedralArrangement((p, q))
            sage: PA.intersects_by_linear_programming(p, q)
            False
            sage: p.intersection(q)
            The empty polyhedron in QQ^4

            sage: cube = polytopes.hypercube(3)
            sage: oct = polytopes.cross_polytope(3)
            sage: PA.intersects_by_linear_programming(cube, oct)
            True
            sage: cube.intersection(oct*2)
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 12 vertices

            sage: P = Polyhedron([(0,0),(1,1)], base_ring=ZZ)
            sage: PA.intersects_by_linear_programming(P, P)
            True
            sage: P.intersection(P)
            A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
        """
        lp = self._set_up_linear_programming_for_intersection(p, q)
        return self._solve_linear_programming(lp)

    def is_all_faces_simplices(self, full_dim=None):
        r"""
        Check if all ``faces`` are simplices

        EXAMPLES::

            sage: simplices = [polytopes.simplex(i) for i in range(1, 8)]
            sage: test_is_simplex = PolyhedralArrangement(simplices)
            sage: test_is_simplex.is_all_faces_simplices()
            True
            sage: square = Polyhedron(vertices=[[0, 0], [0, 1], [1, 0], [1, 1]])
            sage: pa = PolyhedralArrangement([square])
            sage: pa.is_all_faces_simplices()
            False
        """
        for f in self:
            if not self.is_simplex(f):
                return False
        return True

    def is_dissection(self):
        r"""
        Check if ``self`` is a polyhedral dissection of iself.
        
        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[[0, 0], [0, 1], [1, 0], [1, 1]])
            sage: p2 = Polyhedron(vertices=[[1/2, 1/2], [1/2, 3/2], [3/2, 1/2], [3/2, 3/2]])
            sage: pa = PolyhedralArrangement([p1, p2])
            sage: pa.is_dissection()
            False
            sage: p2 = Polyhedron(vertices=[[1, 0], [1, 1], [2, 0], [2, 1]])
            sage: pa = PolyhedralArrangement([p1, p2])
            sage: pa.is_dissection()
            True
            sage: pa = PolyhedralArrangement([p1])
            sage: pa.is_dissection()
            True
            sage: p3 = Polyhedron(vertices=[[0, 0]])
            sage: pa = PolyhedralArrangement([p1, p2, p3])
            sage: pa.is_dissection()
            True
            sage: p4 = Polyhedron(vertices=[[0, 0], [0, 1]])
            sage: pa = PolyhedralArrangement([p1, p2, p3, p4])
            sage: pa.is_dissection()
            True
        """
        elements = self.elements()
        n = len(elements)
        for i in range(n):
            for j in range(i+1, n):
                relint_i = self.relative_interior(elements[i])
                relint_j = self.relative_interior(elements[j])
                if not relint_i.is_disjoint_from(relint_j):
                    return False
        return True

    def is_face_to_face(self, show_errors=False):
        r"""
        Check if the intersection of any pair of elements of ``self``
        is the face of both of the elements.

        EXAMPLES::

            sage: # Polyhedral Arrangement, no need to be face-to-face
            sage: up = 2/10; space = 1/10
            sage: pa1 = Polyhedron(vertices=[[1/2, 0 + up], [1/2, 1 + up], [1, 0 + up], [1, 1 + up]])
            sage: pa2 = Polyhedron(vertices=[[3/2, 1/2 + up], [1, 0 + up], [1, 1 + up]])
            sage: pa3 = Polyhedron(vertices=[[1/2, 1/2 + up], [3/2, 1/2 + up]])
            sage: PA = PolyhedralArrangement([pa1, pa2, pa3])
            sage: PA.is_face_to_face()
            False
            sage: # Polyhedral Dissection, no need to be face-to-face
            sage: pd1 = Polyhedron(vertices=[[1/2, 0 + up], [1/2, 1 + up], [1, 0 + up], [1, 1 + up]])
            sage: pd2 = Polyhedron(vertices=[[3/2, 1/2 + up], [1, 1/2 + up], [1, 1 + up]])
            sage: pd3 = Polyhedron(vertices=[[3/2, 1/2 + up], [1, 1/2 + up], [1, 0 + up]])
            sage: PD = PolyhedralArrangement([pd1, pd2, pd3])
            sage: PD.is_face_to_face()
            False
            sage: # Polyhedral Complex
            sage: pc1 = Polyhedron(vertices=[[3/2, 1/2 + up], [1, 1/2 + up], [1, 1 + up]])
            sage: pc2 = Polyhedron(vertices=[[3/2, 1/2 + up], [1, 1/2 + up], [1, 0 + up]])
            sage: pc3 = Polyhedron(vertices=[[1/2, 0 + up], [1/2, 1 + up], [1, 1/2 + up], [1, 1 + up]])
            sage: pc4 = Polyhedron(vertices=[[1/2, 0 + up], [1, 1/2 + up], [1, 0 + up]])
            sage: PC = PolyhedralArrangement([pc1, pc2, pc3, pc4])
            sage: PC.is_face_to_face()
            True
            sage: # Polyhedral Triangulation
            sage: pt1 = Polyhedron(vertices=[[3/2, 1/2 + up], [1, 1/2 + up], [1, 1 + up]])
            sage: pt2 = Polyhedron(vertices=[[3/2, 1/2 + up], [1, 1/2 + up], [1, 0 + up]])
            sage: pt3 = Polyhedron(vertices=[[1/2, 1 + up], [1, 1/2 + up], [1, 1 + up]])
            sage: pt4 = Polyhedron(vertices=[[1/2, 0 + up], [1/2, 1 + up], [1, 1/2 + up]])
            sage: pt5 = Polyhedron(vertices=[[1/2, 0 + up], [1, 1/2 + up], [1, 0 + up]])
            sage: PT = PolyhedralArrangement([pt1, pt2, pt3, pt4, pt5])
            sage: PT.is_face_to_face()
            True
        """
        faces_of_poly = self.faces_of_polyhedra()
        elements = self.elements()
        n = self.number_of_polyhedra()
        for i in range(n):
            for j in range(i, n):
                p = elements[i]
                q = elements[j]
                int_of_p_q = p.intersection(q)
                if dim(int_of_p_q) >= 0:
                    if (not int_of_p_q in faces_of_poly[p]) or (not int_of_p_q in faces_of_poly[q]):
                        if show_errors:
                            print "Polyhedron with vertices", p.vertices_list(), \
                                  "and polyhedron with vertices", q.vertices_list(), \
                                  "have intersection with vertices", int_of_p_q.vertices_list(), \
                                  "but the intersection is not a faces of both of them \n"
                        return False
        return True

    def is_fine_triangulation(self, trig, pt_config):
        r"""
        Check if a triangulation is fine.

        INPUT::

        - ``trig`` -- a list of polyhedrons for the triangulation
        - ``pt_config`` -- a list of tuples to represent the point configuration

        EXAMPLES::

            sage: pa = PolyhedralArrangement([Polyhedron(vertices=[[0, 0], [0, 1], [1, 0], [1, 1]])])
            sage: pt_config = set([(0, 0), (0, 1), (1, 0), (1, 1)])
            sage: pa.is_fine_triangulation([Polyhedron(vertices=[[0, 0], [0, 1], [1, 0]])], pt_config)
            False
            sage: pa.is_fine_triangulation([Polyhedron(vertices=[[0, 0], [0, 1], [1, 0], [1, 1]])], pt_config)
            True
        """
        trig_vert = set([tuple(vertex) for t in trig for vertex in t.vertices_list()])
        if set(pt_config).issubset(trig_vert):
            return True
        else:
            return False

    def is_polyhedron_in_polyhedron(self, p, q):
        r"""
        (Required by Shuidie)
        
        Test if polyhedron ``p`` is inside polyhedron ``q``, i.e., every
        vertex of ``p`` is in ``q.
        
        Notice that ``p`` does not need to be an element of ``self``,
        and ``q`` is an element of ``self``.

        sage: p1 = Polyhedron(vertices=[[0, 4], [0, 3], [1, 4], [1, 3]])
        sage: p2 = Polyhedron(vertices=[[0, 4], [0, 2], [2, 4], [2, 2]])
        sage: p3 = Polyhedron(vertices=[[4, 4], [4, 2], [2, 4], [2, 2]])
        sage: p4 = Polyhedron(vertices=[[0, 0], [0, 4], [4, 0], [4, 4]])
        sage: pa = PolyhedralArrangement([p1, p2, p3, p4])
        sage: pa.is_polyhedron_in_polyhedron(p1, p2)
        True
        sage: pa.is_polyhedron_in_polyhedron(p1, p4)
        True
        sage: pa.is_polyhedron_in_polyhedron(p1, p3)
        False
        sage: all(pa.is_polyhedron_in_polyhedron(p, p4) for p in [p1, p2, p3, p4])
        True
        """
        if self.grid_setting() == "one_bucket":
            vertices = p.vertices_list()
            for v in vertices:
                if not q.contains(v):
                    return False
            return True
        else:
            p_buckets = self.bounding_box_in_grid(p)
            q_buckets = self.in_buckets(q)
            if not p_buckets.issubset(q_buckets):
                return False
            else:
                vertices = p.vertices_list()
                for v in vertices:
                    if not self.polyhedron_contains_point(q, v):
                        return False
            return True

    def is_simplex(self, face, check_full_dim=False):
        r"""
        Check if a face is a simplex.

        INPUT::

        - ``face`` -- a :class:`Polyhedron` instance

        - ``full_dim`` -- (default: False) check if the face
          is full dimensional, i.e. have the same dimension as
          its ambient dimension

        EXAMPLES::

            sage: triangle_in_3d = Polyhedron(vertices=[[0, 0, 0], [1, 0, 0], [1, 1, 0]])
            sage: tetrahedron_with_ray = Polyhedron(vertices=[[0, 0, 0], [1, 0, 0], [1, 1, 0], [1/2, 1/2, 1]], rays=[[0, 0, 1]])
            sage: point_in_3d = Polyhedron(vertices=[[0, 0, 0]])
            sage: square_in_3d = Polyhedron(vertices=[[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
            sage: simplex = polytopes.simplex(2)
            sage: test_is_simplex = PolyhedralArrangement([triangle_in_3d, tetrahedron_with_ray, point_in_3d, square_in_3d, simplex])
            sage: test_is_simplex.is_simplex(triangle_in_3d)
            True
            sage: test_is_simplex.is_simplex(tetrahedron_with_ray)
            False
            sage: test_is_simplex.is_simplex(point_in_3d)
            True
            sage: test_is_simplex.is_simplex(square_in_3d)
            False
            sage: test_is_simplex.is_simplex(simplex)
            True
        """
        # full_dim is the ambient dimension
        # 0) check if it is full dimensional
        affine_dim = face.affine_hull().dim()
        if check_full_dim:
            if face.ambient_dim() != affine_dim:
                return False
        # 1) compactness
        if not face.is_compact():
            return False
        # 3) has n+1 vertices, where n is the affine dimension
        if len(face.vertices_list()) != (affine_dim + 1):
            return False
        return True

    def lowest_dim_contains(self, p):
        r"""
        (Required by Shuidie)
        
        Return an element of ``self`` that contains ``p`` and
        has the lowest dimension among other elements that
        contains ``p``.

        Notice that ``p`` does not need to be an element of ``self``.

        INPUT::

        - ``p`` -- a :class:`Polyhedron`

        OUTPUT::

        - a :class:`Polyhedron` that is an element of ``self``
        containing ``p`` and has the lowest dimension among
        other elements that contains ``p``.

        sage: p1 = Polyhedron(vertices=[[0, 4], [0, 3], [1, 4], [1, 3]])
        sage: p2 = Polyhedron(vertices=[[0, 4], [0, 2], [2, 4], [2, 2]])
        sage: p3 = Polyhedron(vertices=[[4, 4], [4, 2], [2, 4], [2, 2]])
        sage: p4 = Polyhedron(vertices=[[0, 0], [0, 2], [2, 0], [2, 2]])
        sage: p5 = Polyhedron(vertices=[[0, 0], [0, 4], [4, 0], [4, 4]])
        sage: p6 = Polyhedron(vertices=[[0, 0], [0, 1]])
        sage: s = [p1, p2, p3, p4, p5]
        sage: pa = PolyhedralArrangement(s)
        sage: pa.lowest_dim_contains(p6).vertices_list()
        [[0, 0], [0, 2], [2, 0], [2, 2]]
        sage: pa.lowest_dim_contains(p1).vertices_list()
        [[0, 3], [0, 4], [1, 3], [1, 4]]
        """
        p_dim = p.dim()
        coll_dict = self.collection_dict()
        keys = coll_dict.keys()
        
        # if the dimension of p is less than the dimensions of all polyhedra in ``self``
        if p_dim < keys[0]:
            key_index = 0
        else:
            key_index = keys.index(p_dim)
        # if there exists higher dimension polyhedra than p
        for i in keys[key_index:]:
            for poly in coll_dict[i]:
                if self.is_polyhedron_in_polyhedron(p, poly):
                    return poly
        print "can't find a polyhedron that contains a polyhedron with vertices", p.vertices_list()
        return None

    @cached_method
    def relative_interior(self, p):
        r"""
        Construct a ``NNC_polyhedron`` to represent the relative interior of ``p``.
        
        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[[0, 0], [0, 1], [1, 0], [1, 1]])
            sage: p2 = Polyhedron(vertices=[[1/2, 1/2], [1/2, 3/2], [3/2, 1/2], [3/2, 3/2]])
            sage: p3 = Polyhedron(vertices=[[0, 0]])
            sage: p4 = Polyhedron(vertices=[[0, 0], [0, 1]])
            sage: p5 = Polyhedron(vertices=[[0, 0], [1, 0], [1/2, 1/2]])
            sage: pa = PolyhedralArrangement([p1, p2, p3, p4, p5])
            sage: pa.relative_interior(p1).constraints()
            Constraint_System {x0>0, x1>0, -x0+1>0, -x1+1>0}
            sage: p1.Hrepresentation()
            (An inequality (1, 0) x + 0 >= 0,
             An inequality (0, 1) x + 0 >= 0,
             An inequality (-1, 0) x + 1 >= 0,
             An inequality (0, -1) x + 1 >= 0)
            sage: pa.relative_interior(p2).constraints()
            Constraint_System {2*x0-1>0, 2*x1-1>0, -2*x0+3>0, -2*x1+3>0}
            sage: p2.Hrepresentation()
            (An inequality (2, 0) x - 1 >= 0,
             An inequality (0, 2) x - 1 >= 0,
             An inequality (-2, 0) x + 3 >= 0,
             An inequality (0, -2) x + 3 >= 0)
            sage: pa.relative_interior(p3).constraints()
            Constraint_System {x1==0, x0==0}
            sage: p3.Hrepresentation()
            (An equation (0, 1) x + 0 == 0, An equation (1, 0) x + 0 == 0)
            sage: pa.relative_interior(p4).constraints()
            Constraint_System {x0==0, x1>0, -x1+1>0}
            sage: p4.Hrepresentation()
            (An equation (1, 0) x + 0 == 0,
             An inequality (0, 1) x + 0 >= 0,
             An inequality (0, -1) x + 1 >= 0)
            sage: pa.relative_interior(p5).constraints()
            Constraint_System {x1>0, x0-x1>0, -x0-x1+1>0}
            sage: p5.Hrepresentation()
            (An inequality (0, 1) x + 0 >= 0,
             An inequality (1, -1) x + 0 >= 0,
             An inequality (-1, -1) x + 1 >= 0)
        """
        def construct_strict_inequality(cstr):
            s = sum(variables[index] * coef for index, coef in enumerate(cstr[1:]))
            return s > -cstr[0]
        def construct_equation(cstr):
            s = sum(variables[index] * coef for index, coef in enumerate(cstr[1:]))
            return s == -cstr[0]
        variables = [Variable(i) for i in range(p.ambient_dim())]
        equations = p.equations_list()
        inequalities = p.inequalities_list()
        if len(equations) > 0:
            cs = Constraint_System(construct_equation(equations[0]))
            if len(inequalities) > 0:
                cs.insert(construct_strict_inequality(inequalities[0]))
        else:
            cs = Constraint_System(construct_strict_inequality(inequalities[0]))
        for eqn in equations[1:]:
            cs.insert(construct_equation(eqn))
        for ieq in inequalities[1:]:
            cs.insert(construct_strict_inequality(ieq))
        return NNC_Polyhedron(cs)

    def number_of_polyhedra(self):
        r"""
        Return the numbers of polyhedra of arrangement.
        """
        return self._number_of_polyhedra

    def plot(self, polys=None, **kwds):
        r"""
        Plot some or all elements of ``self``.

        INPUT::

        - ``polys`` -- (default: None) if None, then assign ``polys``
          as the elements of ``self``, else take ``polys`` as a list,
          tuple, or set of polyhedra given from the user.

        - **kwds -- keywords to pass on :func:`plot`
        """
        if polys == None:
            polys = self.elements()
        return sum(p.plot(**kwds) for p in polys)

    def point_lookup(self, point):
        r"""
        Given a point, find where the bucket that the point belong to,
        and return the contents on that bucket as a set.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[[0, 0], [0, 1], [1, 0], [1, 1]])
            sage: q = Polyhedron(vertices=[[0, 0], [0, 2], [2, 0], [2, 2]])
            sage: pa = PolyhedralArrangement([p, q])
            sage: for poly in pa.point_lookup([1/5, 2/5]):
            ....:     poly.vertices_list()
            [[0, 0], [0, 1], [1, 0], [1, 1]]
            [[0, 0], [0, 2], [2, 0], [2, 2]]
            """
        if self.grid_setting() == "one_bucket":
            return self.grid()
        bucket = self.grid().point_in_bucket(point)
        return self.grid_contents()[bucket]

    def point_in_polyhedra(self, point):
        r"""
        (Required by Shuidie)

        Find which polyhedra the point locates on.
        Return a list of those polyhedra.

        EXAMPLES::

            sage: logging.disable(logging.INFO)
            sage: h1 = gmic()
            sage: poly = construct_2d_strip_by_FastPiecewise(h1)
            sage: pa = PolyhedralArrangement(poly)
            sage: point = (4/5, 1)
            sage: intersect_poly = pa.point_in_polyhedra(point)
            sage: for p in intersect_poly:
            ....:         p.vertices()
            (A vertex at (1, 0),
             A vertex at (1, 1),
             A vertex at (4/5, 0),
             A vertex at (4/5, 1))
            (A vertex at (0, 0),
             A vertex at (0, 1),
             A vertex at (4/5, 0),
             A vertex at (4/5, 1))
            (A vertex at (4/5, 0), A vertex at (4/5, 1))
            sage: pa_one_bucket = PolyhedralArrangement(poly, grid_setting="one_bucket")
            sage: intersect_poly_one_bucket = pa_one_bucket.point_in_polyhedra(point)
            sage: set(intersect_poly_one_bucket) == set(intersect_poly)
            True
            sage: h2 = rlm_dpl1_extreme_3a()
            sage: poly = construct_2d_strip_by_FastPiecewise(h2)
            sage: pa = PolyhedralArrangement(poly)
            sage: point = (5/8, 1/2)
            sage: intersect_poly = pa.point_in_polyhedra(point)
            sage: for p in intersect_poly:
            ....:         p.vertices()
            (A vertex at (1, 0),
             A vertex at (1, 1),
             A vertex at (5/8, 0),
             A vertex at (5/8, 1))
            (A vertex at (1/4, 0),
             A vertex at (1/4, 1),
             A vertex at (5/8, 0),
             A vertex at (5/8, 1))
            (A vertex at (5/8, 0), A vertex at (5/8, 1))
        """
        if self.grid_setting() == "one_bucket":
            return [p for p in self if p.contains(point)]
        else:
            possible_intersect = self.point_lookup(point)
            return [p for p in possible_intersect if p.contains(point)]

    def polyhedron_buckets(self):
        r"""
        Return a dictionary whose keys are elements and whose values are
        lists of buckets that intersect with the elements.
        """ 
        if self.grid_setting() == "one_bucket":
            return self.grid()
        self.update_dictionaries()
        return self._polyhedron_buckets

    def polyhedron_contains_point(self, polyhedron, point):
        r"""
        Return True if ``polyhedron`` contains ``point``.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[[0, 0], [0, 1], [1, 0], [1, 1]])
            sage: q = Polyhedron(vertices=[[0, 0], [0, 2], [2, 0], [2, 2]])
            sage: pa = PolyhedralArrangement([p, q])
            sage: all(pa.polyhedron_contains_point(poly, [1/5, 2/5]) for poly in [p, q])
            True
            sage: pa_one_bucket = PolyhedralArrangement([p, q], grid_setting="one_bucket")
            sage: all(pa_one_bucket.polyhedron_contains_point(poly, [1/5, 2/5]) for poly in [p, q])
            True
            sage: pa.polyhedron_contains_point(p, [3/2, 0]); pa.polyhedron_contains_point(q, [3/2, 0])
            False
            True
            sage: pa_one_bucket.polyhedron_contains_point(p, [3/2, 0]); pa_one_bucket.polyhedron_contains_point(q, [3/2, 0])
            False
            True
        """
        if self.grid_setting() == "one_bucket":
            return polyhedron.contains(point)
        else:
            pt_bucket = self.grid().point_in_bucket(point)
            buckets = self.in_buckets(polyhedron)
            if not (pt_bucket in buckets):
                return False
            else:
                return polyhedron.contains(point)

    def update_dictionaries_for_a_polyhedron(self, p, lower_left=None, upper_right=None, lp=None, x=None):
        r"""
        If grid_setting is "default", then use grid method to find what buckets
        the polyhedron ``p`` located in and update two dictionaries: ``grid_contents``
        and ``polyhedron_buckets``.

        Recall:

        - ``grid_contents`` -- a dictionary whose keys are buckets and whose values are
          sets of polyhedra that intersect with the keys  
        
        - ``polyhedron_buckets`` -- a dictionary whose keys are elements and whose values are
          lists of buckets that intersect with the elements. This attribute is like an inverse of
          ``grid_contents``

        INPUT:

        - ``p`` -- a :class:`Polyhedron`

        - ``selection`` -- (default: None) a tuple of two buckets
          (buckets are represented by its integer coordinates,
          see :class:`Grid` for examples)
          if ``selection`` is None, then choose the entire grid as
          selection

        EXAMPLES::

            sage: p = Polyhedron(vertices=[[8/21, 1/4], [4/5, 1/5], [5/6, 4/5]])
            sage: q = Polyhedron(vertices=[[1/2, 1/2], [4/5, 1/2], [1/2, 5/6]])
            sage: x = Polyhedron(vertices=[[4/5, 4/5]])
            sage: pa = PolyhedralArrangement([p, q, x], lower_left_coordinates=None, side_length=1, bucket_width=None)
            sage: pa.update_dictionaries_for_a_polyhedron(p)
            sage: pa._grid_contents
            defaultdict(<type 'set'>, {(2, 0): set([A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices]), (1, 0): set([A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices]), (1, 1): set([A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices]), (2, 1): set([A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices]), (2, 2): set([A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices])})
            sage: pa._polyhedron_buckets
            defaultdict(<type 'set'>, {A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices: set([(2, 0), (1, 0), (1, 1), (2, 1), (2, 2)])})

            sage: m1 = Polyhedron(vertices=[[-1, -1], [0, 0]])
            sage: m2 = Polyhedron(vertices=[[1, -1], [0, 0]])
            sage: m3 = Polyhedron(vertices=[[1, -1], [2, 0]])
            sage: m4 = Polyhedron(vertices=[[3, -1], [2, 0]])
            sage: k1 = Polyhedron(vertices=[[-1/2, -1/2], [-1/2, 2 + 1/2]])
            sage: k2 = Polyhedron(vertices=[[-1/2, 1 + 1/2], [1, 2]])
            sage: k3 = Polyhedron(vertices=[[-1/2, 1 + 1/2], [1 + 1/2, -1/2]])
            sage: o = Polyhedron(vertices=[[2, 2], [2, 3], [3, 2], [3, 3]])
            sage: s = [m1, m2, m3, m4, k1, k2, k3, o]
            sage: pa = PolyhedralArrangement(s, lower_left_coordinates=(-1, -1), side_length=4, bucket_width=1)
            sage: len(pa.polyhedron_buckets())
            8
            sage: len(pa.grid_contents())
            16
        """
        if self.grid_setting() == "one_bucket":
            return self.grid()

        grid = self.grid()
        if lower_left is None:
            dim = grid.dim()
            lower_left = tuple(0 for i in range(dim))
            upper_right = tuple(grid.size() - 1 for i in range(dim))
            lp, x = self._set_up_linear_programming(p)
        self._update_variables_bounds(lower_left, upper_right, lp, x)
        if self._solve_linear_programming(lp):
            if lower_left == upper_right:
                grid._contents[lower_left].add(p)
                self._polyhedron_buckets[p].add(lower_left)
            else:
                a, b = grid.cut_selection(lower_left, upper_right)
                self.update_dictionaries_for_a_polyhedron(p, lower_left=a[0], upper_right=a[1], lp=lp, x=x)
                self.update_dictionaries_for_a_polyhedron(p, lower_left=b[0], upper_right=b[1], lp=lp, x=x)

    @cached_method
    def update_dictionaries(self, lower_left=None, upper_right=None):
        r"""
        If grid_setting is "default", then, for each polyhedron in the collection, 
        use grid method to find what buckets the polyhedron located in and update
        two dictionaries: ``grid_contents``and ``polyhedron_buckets``.

        See :meth:`update_dictionaries_for_a_polyhedron` for more
        docstrings.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[[8/21, 1/4], [4/5, 1/5], [5/6, 4/5]])
            sage: q = Polyhedron(vertices=[[1/2, 1/2], [5/6, 1/2], [1/2, 1]])
            sage: x = Polyhedron(vertices=[[4/5, 4/5]])
            sage: pa = PolyhedralArrangement([p, q, x], lower_left_coordinates=None, side_length=1, bucket_width=None)
            sage: pa.update_dictionaries()
            sage: pa.grid_contents()
            defaultdict(<type 'set'>, {(1, 2): set([A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices]), (2, 1): set([A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices, A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices]), (2, 0): set([A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices]), (2, 2): set([A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices, A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices, A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex]), (1, 0): set([A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices]), (1, 1): set([A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices, A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices])})
            sage: pa.polyhedron_buckets()
            defaultdict(<type 'set'>, {A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices: set([(2, 0), (1, 0), (1, 1), (2, 1), (2, 2)]), A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices: set([(1, 2), (1, 1), (2, 1), (2, 2)]), A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex: set([(2, 2)])})
            sage: pa = PolyhedralArrangement([p, q, x], lower_left_coordinates=None, side_length=1, bucket_width=1)
            sage: pa.update_dictionaries()
            sage: pa.grid_contents()
            defaultdict(<type 'set'>, {(0, 0): set([A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices, A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices, A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex])})
        """
        if self.grid_setting() == "one_bucket":
            return self.grid()
        for p in self.elements():
            self.update_dictionaries_for_a_polyhedron(p, lower_left, upper_right)

    def random_vertex(self, digits):
        return tuple([randint(0, 10^digits) / 10^digits for i in range(self.ambient_dim())])

    def regular_triangulation_of_point_configuration(self, pt_config, dim):
        r"""
        Return a regular triangulation corresponding to the point
        configuration ``pt_config``

        INPUT::

        - ``pt_config`` -- a list of tuples to represent the points configuration
        - ``dim`` -- an integer indicating the dimension of the polyhedron we want to
          triangulate

        OUTPUT::

        a list of polyhedra that forms a triangulation
        """
        # dim is the dimension of the polyhedron
        # 1) construct a one dimension higher polyhedron by adding
        # the weights to all the points in the point configuration
        l = len(pt_config)
        # generate perturbations for the weights function
        old_weights = [self.weight(p) for p in pt_config]
        weights = self._small_perturbation(old_weights)
        new_vertices_list = [pt_config[i] + tuple(weights[i]) for i in range(l)]
        
        # 2) assign a ray in the direction of weights to the higher dimemsion polyhedron
        weight_direction = [[0] * self.ambient_dim() + [1]]
        higher_dim_poly = Polyhedron(vertices=new_vertices_list, rays=weight_direction)
        
        # 3) make a projection
        higher_dim_F = higher_dim_poly.faces(dim)

        # higher_maximal_faces is a set of faces with an extra dimension
        # 4) so we need to remove the extra dimension
        vert_lists = [[v[:-1] for v in face.as_polyhedron().vertices_list()] for face in higher_dim_F]
        
        # since the dimension of the face is not the same as the
        # dimension of the face as a polyhedron
        # 5) make a new polyhedron to get the true dimension
        max_dim = 0
        trigs = set()
        for vl in vert_lists:
            poly = Polyhedron(vertices=vl)
            true_dim = poly.dim()
            if true_dim > max_dim:
                trigs = set([poly])
                max_dim = true_dim
            elif true_dim == max_dim:
                trigs.add(poly)
        return trigs

    def regular_triangulation(self, p, digits=3):
        r"""
        Make a fine regular triangulation for the point configuration of ``p``,
        then return the triangulation as a list of :class:`Polyhedron`.
        """
        pt_config = p.vertices()
        trig = self.regular_triangulation_of_point_configuration(pt_config, p.dim())
        pa_trig = PolyhedralArrangement(trig, grid_setting='one_bucket')
        if not pa_trig.is_face_to_face(show_errors=True):
            raise ValueError("not face to face")
        if pa_trig.is_all_faces_simplices():
            return pa_trig.collection_of_polyhedra()
        else:
            raise ValueError("not all the face is simplex")

    def volume(self):
        r"""
        Return the volume as the sum of the polhedral triangulations of the
        polyhedral arrangement.

        sage: square_in_3D = Polyhedron(vertices=[[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]])
        sage: pa = PolyhedralArrangement([square_in_3D])
        sage: pa.volume()
        0
        sage: c3 = polytopes.hypercube(3)
        sage: c3.vertices_list()
        [[-1, -1, -1],
         [-1, -1, 1],
         [-1, 1, -1],
         [-1, 1, 1],
         [1, -1, -1],
         [1, -1, 1],
         [1, 1, -1],
         [1, 1, 1]]
        sage: pa = PolyhedralArrangement([c3])
        sage: pa.volume() == c3.volume()
        True
        """
        return inclusion_exclusion_on_full_dim_volumes(self.elements())

    def weight(self, p):
        r"""
        The weight of ``p`` is defined as the square of the Euclidean distance
        between ``p`` and the origin.

        EXAMPLES::

            sage: v = [1/2, 1/2, 1/2]
            sage: p = Polyhedron(vertices=[v])
            sage: pa = PolyhedralArrangement([p])
            sage: pa.weight(v)
            3/4
        """
        return sum([i^2 for i in p])

class Grid(object):
    r"""
    Define a grid that are divided into buckets. The grid
    will be used in the grid method to do calculations on elements
    of :class:`PolyhedralArrangement`.

    Definition of terms:
    1. bucket - a single bucket in the grid with integer coordinates.
    2. selection - a list of buckets which can be determined by 
       the lower left bucket and upper right bucket.

    :class:`Grid` is supposed to now nothing about the elements in
    :class:`PolyhedralArrangement`, but only the number of elements
    in :class:`PolyhedralArrangement`.

    INPUT:

    - ``number_of_polyhedra`` -- an integer, the number of polyhedra
      in the collection in a :class:`PolyhedralArrangement`

    - ``dim`` -- the dimension of the grid

    - ``side_length`` -- (default: 1) the length of the hypercube as ``self``

    - ``lower_left_coordinates`` -- (default: the origin) the real coordinates
      of the lower left bucket

    - ``bucket_width`` -- the width of the bucket

    Attributes:

    - ``_bucket_width`` -- the width of each bucket

    - ``_contents`` -- a dictionary, the keys are buckets
      and the values are the contents on the buckets. In
      this codes (so far), the contents are sets of polyhedra
      that intersect with the keys

    - ``_dim`` -- an integer that indicates the dimension
      of the grid

    - ``_lower_left_coordinates`` -- (default: None) a tuple or a list that
      indicates the real coordinates of the lower left point of the lower left bucket

    - ``_side_length`` -- (default: 1) an integer that 
      indicates the length of the hypercube as ``self``

    - ``_size`` -- an integer that indicates how many pieces
      are divided from 0 to 1

    EXAMPLES::

        sage: g = Grid(1, 2, lower_left_coordinates=(-1, -2), side_length=1, bucket_width=1/3)
        sage: g
        A 2-dimensional grid with side length 1.
        The grid has 3^2 buckets. Each bucket has width 1/3.
        The lower left coordinate of the grid is (-1, -2) 
    """ 
    def __init__(self, number_of_polyhedra, dim, grid_setting="default", lower_left_coordinates=None, side_length=1, bucket_width=None):
        if grid_setting == "default":
            if lower_left_coordinates is None:
                lower_left_coordinates = tuple(0 for i in range(dim))
            if side_length == 0:
                side_length = 1
            if bucket_width is None:
                self._bucket_width = side_length / number_of_polyhedra
            elif bucket_width == "ideal":
                self._bucket_width = ideal_width(number_of_polyhedra, dim, side_length=side_length)
            else:
                self._bucket_width = bucket_width
            # size is the number of buckets in each dimension
            # side_length = size * bucket_width
            self._size = Integer(side_length / self._bucket_width)
            self._dim = dim
            self._side_length = side_length
            self._lower_left_coordinates = lower_left_coordinates
            self._contents = defaultdict(set)
            self._grid_setting = grid_setting
        else:
            self._grid_setting = grid_setting

    def __repr__(self):
        if self._grid_setting == "default":
            desc = 'A '
            desc += repr(self._dim)
            desc += '-dimensional grid with side length '
            desc += repr(self._side_length)
            desc += '.\nThe grid has '
            desc += repr(self._size)
            desc += '^'
            desc += repr(self._dim)
            desc += ' buckets. Each bucket has width '
            desc += repr(self._bucket_width)
            desc += '.\n'
            desc += 'The lower left coordinate of the grid is '
            desc += repr(self._lower_left_coordinates)
            return desc
        else:
            desc = "The current grid setting is one bucket. \n"
            desc += "No functions on the grid is functional is this setting. \n"
            desc += "To avoid errors, don't call any functions on grid."
            return desc

    def bucket_width(self):
        r"""
        Return the width of the grid.

        EXAMPLES::

            sage: g = Grid(2, 3)
            sage: g.bucket_width()
            1/2
        """
        return self._bucket_width

    def change_coordinate(self, bucket, index, offset):
        r"""
        Change one coordinate of the given bucket.

        EXAMPLES::

            sage: g = Grid(3, 3)
            sage: g.change_coordinate((0, 1, 2), 2, -1)
            (0, 1, 1)
        """
        return bucket[:index] + (bucket[index] + offset,) + bucket[1 + index:]

    def contents(self):
        r"""
        A dictionary, the keys are buckets and the values are the contents
        on the buckets. In this codes (so far), the contents are sets of
        polyhedra that intersect with the keys.
        """
        return self._contents

    def cut_selection(self, lower_left, upper_right):
        r"""
        Given the lower left bucket and upper right bucket of the selection,
        find the 'longest' edges, and cut it into two new selections.

        If all edges are equal, cut the selection from its first coordinate.

        Return two selections with their lower left and upper right buckets.

        EXAMPLES::

            sage: g = Grid(3, 3)
            sage: g.cut_selection((0, 0, 0), (2, 2, 2))
            (((0, 0, 0), (0, 2, 2)), ((1, 0, 0), (2, 2, 2)))
            sage: g.cut_selection((1, 1, 1), (1, 1, 1))
            ((1, 1, 1), (1, 1, 1))
        """
        if lower_left == upper_right:
            return lower_left, lower_left
        index, length = self.longest_edge_of_selection(lower_left, upper_right)
        lower_left_offset = length // 2
        # cut the selection into two new selection
        # assume the lower left - is smaller if length is odd
        small = (lower_left, self.change_coordinate(upper_right, index, -(length - lower_left_offset)))
        large = (self.change_coordinate(lower_left, index, lower_left_offset), upper_right)
        return small, large

    def dim(self):
        r"""
        Return the dimension of the grid.
        """
        return self._dim

    def edges_of_selection(self, lower_left, upper_right):
        r"""
        Return the lengths of the edges of the selection.

        EXAMPLES::

            sage: g = Grid(3, 3)
            sage: g.edges_of_selection((0, 1, 1), (1, 2, 1))
            [2, 2, 1]
        """
        return [c1 - c0 + 1 for c0, c1 in zip(lower_left, upper_right)]

    def longest_edge_of_selection(self, lower_left, upper_right):
        r"""
        Find the longest edge of a selection.

        EXAMPLES::

            sage: g = Grid(3, 3)
            sage: g.longest_edge_of_selection((0, 1, 1), (1, 2, 1))
            (0, 2)
        """
        edges = self.edges_of_selection(lower_left, upper_right)
        maximal_edge = max(edges)
        index = edges.index(maximal_edge)
        return index, maximal_edge

    def point_in_bucket(self, point):
        r"""
        Return a bucket where the point are located in.

        EXAMPLES:

            sage: g = Grid(3, 3)
            sage: g.point_in_bucket((1/4, 1/2, 1))
            (0, 1, 2)
            sage: g = Grid(3, 2, side_length=3, bucket_width=3/6, lower_left_coordinates=(1, 3))
            sage: g.point_in_bucket((1, 3))
            (0, 0)
            sage: g.point_in_bucket((7/2, 3))
            (5, 0)
            sage: g.point_in_bucket((7/2, 7/2))
            (5, 1)
            sage: g.point_in_bucket((7/2, 11/2))
            (5, 5)
            sage: g.point_in_bucket((7/2, 11/2 - 1/10))
            (5, 4)
            sage: g.point_in_bucket((7/2, 11/2 + 1/10))
            (5, 5)
            sage: g.point_in_bucket((3, 6))
            (4, 5)
            sage: g.point_in_bucket((4, 6))
            (5, 5)
        """
        dim = self._dim
        if len(point) != dim:
            raise ValueError("the dimension of the point does not match with the dimension of the grid")
        ll_coord = self._lower_left_coordinates
        sl = self._side_length
        bw = self._bucket_width
        coord_max = [ll_coord[i] + sl for i in range(dim)]
        return tuple(((point[i] - ll_coord[i]) // bw).floor() if point[i] != coord_max[i]
                    else self._size - 1
                    for i in range(len(point)))

    def real_coordinate(self, bucket):
        r"""
        Get the actual coordinate of the bucket.

        EXAMPLES::

            sage: g = Grid(3, 3)
            sage: g.real_coordinate((0, 1, 2))
            (0, 1/3, 2/3)
            sage: g = Grid(3, 2, side_length=3, bucket_width=3/6, lower_left_coordinates=(1, 3))
            sage: g.real_coordinate((0, 0))
            (1, 3)
            sage: g.real_coordinate((5, 0))
            (7/2, 3)
            sage: g.real_coordinate((5, 5))
            (7/2, 11/2)
        """
        ll_coord = self._lower_left_coordinates
        return tuple(bucket[i] * self.bucket_width() + ll_coord[i] for i in range(len(bucket)))

    def variables_min_offset(self, lower_left, upper_right):
        r"""
        Consider using the selection in a linear programming,
        given the lower left and the upper right buckets of
        the selection, the method will return the miminum of
        the variables and the offsets (mimimum + offset = maximum)
        for a selection.

        EXAMPLES::

            sage: g = Grid(3, 3)
            sage: # build a flat box 
            sage: # (i.e. the height is half of the width or of the length)
            sage: g.variables_min_offset((1, 1, 2), (2, 2, 2))
            ((1/3, 1/3, 2/3), (2/3, 2/3, 1/3))
            sage: g.variables_min_offset((2, 2, 2), (2, 2, 2))
            ((2/3, 2/3, 2/3), (1/3, 1/3, 1/3))
            sage: g = Grid(1, 2, lower_left_coordinates=(-1, -2), side_length=1, bucket_width=1/2)
            sage: g.variables_min_offset((0, 0), (1, 1))
            ((-1, -2), (1, 1))
            sage: g.variables_min_offset((0, 0), (0, 0))
            ((-1, -2), (1/2, 1/2))
        """
        if lower_left == upper_right:
            offsets = tuple(self.bucket_width() for _ in range(self.dim()))
        else:
            edges = self.edges_of_selection(lower_left, upper_right)
            offsets = tuple([e * self.bucket_width() for e in edges])
        return self.real_coordinate(lower_left), offsets

    def size(self):
        r"""
        Return the size of the grid.

        EXAMPLES::

            sage: g = Grid(2, 3)
            sage: g.size()
            2
        """
        return self._size

###############################
# Unique Minimal Lifting
###############################

def _min_max_int_coordinates(p, use_integral_points=True):
    r"""
    Find the minimal and maiximal integers in each coordinates of Polyhedron ``p``
    
    EXAMPLES::

        sage: p = Polyhedron(vertices=[[-1 - 1/10], [0], [1 + 1/10]])
        sage: _min_max_int_coordinates(p)
        [(-2, 2)]
        sage: q = Polyhedron(vertices=[[1/10, 1, 2], [2 + 2/10, -1, 0], [1, 0 , 3 - 1/10]])
        sage: _min_max_int_coordinates(q)
        [(0, 3), (-1, 1), (0, 3)]
    """
    vmatrix = p.vertices_matrix()
    rows = vmatrix.rows()
    if use_integral_points:
        # use ceil and floor to get the boundary
        return [(floor(min(r)), ceil(max(r))) for r in rows]
    else:
        return [(min(r), max(r)) for r in rows]

def bounding_box(p, use_integral_points=True):
    r"""
    Return the bounding box that contains ``p``

    EXAMPLES::

        sage: p = Polyhedron(vertices=[[-1 - 1/10], [0], [1 + 1/10]])
        sage: bounding_box(p).vertices()
        (A vertex at (-2), A vertex at (2))
        sage: q = Polyhedron(vertices=[[1/10, 1, 2], [2 + 2/10, -1, 0], [1, 0 , 3 - 1/10]])
        sage: bounding_box(q).vertices()
        (A vertex at (0, -1, 0),
         A vertex at (0, -1, 3),
         A vertex at (0, 1, 0),
         A vertex at (0, 1, 3),
         A vertex at (3, -1, 0),
         A vertex at (3, -1, 3),
         A vertex at (3, 1, 0),
         A vertex at (3, 1, 3))
    """
    return Polyhedron(vertices=product(*[[i[0], i[1]] for i in _min_max_int_coordinates(p, use_integral_points)]))

def find_volme_of_B_with_f(B, f):
    r"""
    Find the volume of ``B`` with respect to ``f``.

    EXAMPLES::

        sage: triangle = Polyhedron(vertices=[[0, 0], [0, 2], [2, 0]])
        sage: find_volme_of_B_with_f(triangle, pick_f(triangle))
        1
    """
    return find_volume_of_lifting_region(B, lifting_region(B, f))

def find_volume_of_lifting_region(B, lr):
    return inclusion_exclusion_on_full_dim_volumes(translate_lifting_region(B, lr).elements())

def is_lattice_free(B):
    r"""
    INPUT::

    - ``B`` -- a :class:`Polyhedron`

    EXAMPLES::

        sage: thin_diamond = Polyhedron(vertices=[[-1/10, 3/2], [11/10, 3/2], [1/2, 5/2], [1/2, 1/2]])
        sage: square = Polyhedron(vertices=[[0, 0], [0, 1], [1, 0], [1, 1]])
        sage: square_minus_corner = Polyhedron(vertices=[[0, 0], [0, 1], [1, 0], [9/10, 1], [1, 9/10]])
        sage: triangle = Polyhedron(vertices=[[0, 0], [2, 0], [0, 2]])
        sage: contain_point = Polyhedron(vertices=[[-1/10, -1/10], [-1/10, 1/10], [1/10, -1/10], [1/10, 1/10]])
        sage: is_lattice_free(thin_diamond)
        True
        sage: is_lattice_free(square)
        True
        sage: is_lattice_free(square_minus_corner)
        True
        sage: is_lattice_free(triangle)
        True
        sage: is_lattice_free(contain_point)
        False
    """
    for pt in possible_lattices(B):
        if B.interior_contains(pt):
            return False
    return True

def is_maximal_lattice_free(B):
    r"""
    A polyhedron is maximal lattice tree if it is lattice free and every facet
    has at least one lattice point in its relative interior.

    EXAMPLES::

        sage: not_maximal_1 = Polyhedron(vertices=[[-1/5, 1/2], [2/5, 2], [13/10, 1/2], [2/5, -2/5]])
        sage: is_maximal_lattice_free(not_maximal_1)
        False
        sage: not_maximal_2 = Polyhedron(vertices=[[-1/5, 1/2], [2/5, 2], [13/10, 1/2]])
        sage: is_maximal_lattice_free(not_maximal_2)
        False
        sage: not_maximal_3 = Polyhedron(vertices=[[-1/5, 1/2], [2/5, 2], [13/10, 1/2], [1, 0]])
        sage: is_maximal_lattice_free(not_maximal_3)
        False
        sage: not_maximal_4 = Polyhedron(vertices=[[0, 0], [0, 1], [1, 0]])
        sage: is_maximal_lattice_free(not_maximal_4)
        False
        sage: point = Polyhedron(vertices=[[0, 0]])
        sage: is_maximal_lattice_free(point)
        True
        sage: triangle_1 = Polyhedron(vertices=[[-3/13, 21/13], [1 - 4/10, 3], [3/2, 3/4]])
        sage: is_maximal_lattice_free(triangle_1)
        True
        sage: triangle_2 = Polyhedron(vertices=[[0, 0], [0, 2], [2, 0]])
        sage: is_maximal_lattice_free(triangle_2)
        True
        sage: diamond = Polyhedron(vertices=[[-1/2, 1/2], [1/2, 3/2], [1/2, -1/2], [3/2, 1/2]])
        sage: is_maximal_lattice_free(diamond)
        True
    """
    if not is_lattice_free(B):
        return False
    else:
        for f in B.faces(dim(B) - 1):
            facet = f.as_polyhedron()
            integral_points = facet.integral_points()
            if not integral_points:
                # no lattice point contains in the facet
                return False
            for pt in integral_points:
                if facet.relative_interior_contains(pt):
                    break
            else:
                return False
        return True

def inclusion_exclusion_on_full_dim_volumes(polys):
    r"""
    Compute the volume of the union of polyhedra.

    Use inclusion and exclusion property to calculates the sum of the
    full dimension volumes for each polyhedron in ``polys``.

    sage: p1 = Polyhedron(vertices=[[0, 4], [0, 3], [1, 4], [1, 3]])
    sage: p2 = Polyhedron(vertices=[[0, 4], [0, 2], [2, 4], [2, 2]])
    sage: p3 = Polyhedron(vertices=[[4, 4], [4, 2], [2, 4], [2, 2]])
    sage: p4 = Polyhedron(vertices=[[0, 0], [0, 2], [2, 0], [2, 2]])
    sage: p5 = Polyhedron(vertices=[[4, 0], [4, 2], [2, 0], [2, 2]])
    sage: p6 = Polyhedron(vertices=[[3, 0], [3, 1], [4, 0], [4, 1]])
    sage: p7 = Polyhedron(vertices=[[0, 0], [0, 4], [4, 0], [4, 4]])
    sage: p8 = Polyhedron(vertices=[[0, 0], [0, 4]])
    sage: p9 = Polyhedron(vertices=[[0, 0], [1, 2], [3, 0]])
    sage: p10 = Polyhedron(vertices=[[2, 3]])
    sage: p11 = Polyhedron(vertices=[[0, 1], [2, 3], [3, 4]])
    sage: polys = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11]
    sage: inclusion_exclusion_on_full_dim_volumes(polys)
    16
    """
    return sum(inclusion_exclusion_recursion(value, 1, index, polys) for index, value in enumerate(polys))

def inclusion_exclusion_recursion(intersect, sign, first_index, polys):
    if intersect.is_full_dimensional():
        return intersect.volume() * sign + \
                sum(inclusion_exclusion_recursion(intersect.intersection(polys[index]), -sign, index, polys) \
                for index in range(len(polys)) if index > first_index)
    else:
        return 0

@cached_method
def intersection(p, q):
    # used for inclusion exclusion method
    # there will be repeated intersections in the calculasion
    # so we make this method as a cached method
    return p.intersection(q)

def lifting_region(B, f):
    r"""
    INPUT::

    - ``B`` -- a :class:`Polyhedron`

    - ``f`` -- a point in unit cube
    """
    # FIXME:: let's focus on polytopes rather than polyhedra at this point
    if not is_lattice_free(B):
        raise ValueError("B is not lattice free")
    if not f in B: ### WAS: B.interior_contains(f):    ### Must allow f on boundary. Needs some changes below!! --mkoeppe
        raise ValueError("f is not contained in B")
    vl = [tuple(v) for v in B.vertices_list()]
    dim = B.dim()
    f_rays = dict()
    for v in vl:
        f_rays[v] = make_ray(f, v)
    facets = B.faces(dim - 1)
    facets_as_poly = [facet.as_polyhedron() for facet in facets]
    lattices = possible_lattices(B)
    region = set()
    for lattice_pt in lattices:
        # Notice that B is lattice free, if B contains a lattice point,
        # then this point is on one of the facet
        if B.contains(lattice_pt):
            for facet in facets_as_poly:
                if facet.contains(lattice_pt):
                    v_on_facet = [tuple(v) for v in facet.vertices()]
                    facet_rays = [f_rays[v] for v in v_on_facet]
                    f_cone = Polyhedron(rays=facet_rays, vertices=[f])
                    lattice_pt_cone = Polyhedron(rays=[[-i for i in ray] for ray in facet_rays], vertices=[lattice_pt])
                    R_i = f_cone.intersection(lattice_pt_cone)
                    f_cone.plot()
                    lattice_pt_cone.plot()
                    R_i.plot()
                    region.add(R_i)
                    break
    return PolyhedralArrangement(list(region))

def make_ray(start, end):
    r"""
    Given two points, return the ray that points from the
    starting point to the end point

    EXAMPLES::

        sage: make_ray((0, 0, 0), (1, 2, 3))
        (1, 2, 3)
    """
    ls = len(start)
    le = len(end)
    if ls != le:
        raise ValueError("two points have different lengths")
    return tuple(end[i] - start[i] for i in range(ls))

def pick_f(B, digit=2):
    vl = B.vertices_list()
    while(True):
        random_v = vl[randint(0, len(vl)-1)]
        de = randint(1, 10^digit)
        nu_list = [randint(1, de) for i in range(len(vl[0]))]
        f = [(pair[1] - pair[0]) * nu_list[index] / de + random_v[index] for index, pair in enumerate(_min_max_int_coordinates(B))]
        if B.interior_contains(f):
            break
    return f

def possible_lattices(B):
    # return possible integral points in ``B``
    return bounding_box(B).integral_points()

def translate_lifting_region(B, lr, full_dimensional=True):
    r"""
    Given a lifting region, make translations on it.
    
    Return the translations that intersect with the unit hypercube
    as a :class:``PolyhedralArrangement`` instance.

    INPUT::

    - ``B`` -- a polyhedron

    - ``lr`` -- a :class:`PolyhedralArrangement` that is the lifting region of ``B``

    - ``full_dimensional`` -- (deafult: True) a boolean value. If true,
      then add only full dimensional translated lifting regions to the output

    OUTPUT::

    - a :class:`PolyhedralArrangement` consisting of the translated lifting
    regions that intersect with the unit cube
    """
    unit_cube = construct_hypercube(B.ambient_dim())
    translations_list = set([])
    for p in lr:
        bound_box_of_p = bounding_box(p, False)
        # randomly choose a vertex of the bounding box of p as the vector to translate p to the unit cube
        # use ceil, since we only make translation on integral vectors
        vector_to_unit_cube = [i.ceil() for i in bound_box_of_p.vertices_list()[0]]
        coord_offset = [ceil((i[1] - i[0] / 2)) for i in _min_max_int_coordinates(bound_box_of_p)]
        for t in product(*[range(-i, i+1) for i in coord_offset]):
            temp_int = (p + vector(QQ, t) - vector(QQ, vector_to_unit_cube)).intersection(unit_cube)
            if full_dimensional:
                if temp_int.is_full_dimensional():
                    translations_list.add(temp_int)
            else:
                translations_list.add(temp_int)
    return PolyhedralArrangement(translations_list, grid_setting='one_bucket')

def construct_hypercube(dim):
    r"""
    EXAMPLES::

        sage: construct_hypercube(2).vertices_list()
        [[0, 0], [0, 1], [1, 0], [1, 1]]
    """
    unit_grid_iter = product([0, 1], repeat=dim)
    unit_grid = [list(v) for v in unit_grid_iter]
    return Polyhedron(unit_grid)

def facets(dim):
    r"""
    Return the facets (as a list of :class:``Polyhedron``)
    of a hypercube in dimension ``dim``.
    
    EXAMPLES::

        sage: F = facets(3)
        sage: for f in F:
        ....:     f.vertices_list()
        [[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1]]
        [[0, 0, 0], [0, 0, 1], [1, 0, 0], [1, 0, 1]]
        [[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]]
        [[0, 0, 1], [0, 1, 1], [1, 0, 1], [1, 1, 1]]
        [[0, 1, 0], [0, 1, 1], [1, 1, 0], [1, 1, 1]]
        [[1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]]
    """
    F = construct_hypercube(dim).faces(dim-1)
    return [f.as_polyhedron() for f in F]

def lifting_graphics(B, f=None, full_dimensional=True, show_graphics=True):
    r"""
    If ``show_graphics`` is True, then return the lifting regions and
    translated lifting regions of ``B`` with respect to ``f``
    (or randly picked ``f`` if ``f`` is not assigned).

    If ``show_graphics`` is False, then return the sum of the volumes of the
    translated lifting regions.

    ``full_dimensional`` -- (deafult: True) a boolean value. If true,
      then add only full dimensional translated lifting regions to the output
    """
    if f is None:
        f = pick_f(B)
    lr = lifting_region(B, f)
    trl = translate_lifting_region(B, lr, full_dimensional=full_dimensional)
    volume = inclusion_exclusion_on_full_dim_volumes(trl.elements())
    if show_graphics:
        if B.ambient_dim() < 3:
            return graphics_array([lr.plot(color='red',alpha=0.4) + B.plot(alpha=0.4, color='yellow')
                                   + point(f, color='red')
                                   + text('f', [i + 1/9 for i in f], color='black', fontsize=15)
                                   + sum(point(pt, color='gray') for pt in bounding_box(B).integral_points()),
                                   trl.plot(alpha=0.3, color='red')
                                   + text(r'$\mathrm{vol}_{\mathbb{T}^{%s}} \bar{R}(f)/\mathbb{Z}^{%s} = %s$'
                                          % (B.ambient_dim(), B.ambient_dim(), latex(volume)), [1/2, 1.1],
                                          fontsize='x-large')
                                   + text(r'$\mathrm{f} = %s$'
                                          % (latex(f)), [1/2, 1.25],
                                          fontsize='x-large')
                                   ])
        elif B.ambient_dim() == 3:
            P = lr.plot(color='red',alpha=0.4)
            P += B.plot(alpha=0.8, color='blue')
            P += sum(point(pt, color='gray') for pt in bounding_box(B).integral_points())
            P.show()
            return trl.plot(alpha=0.2, color='red')
        else:
            raise ValueError("can't draw graphics for ambient dimension higher than 3!")
    else:
        return volume