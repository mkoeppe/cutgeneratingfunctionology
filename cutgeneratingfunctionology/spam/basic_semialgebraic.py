r"""
Mutable basic semialgebraic sets
"""

from __future__ import division, print_function, absolute_import

from sage.structure.element import Element
from sage.modules.free_module_element import vector
from sage.rings.all import QQ, ZZ
from sage.rings.real_double import RDF
from sage.misc.abstract_method import abstract_method

try:
    from ppl import Variable, Constraint, Linear_Expression, Constraint_System, NNC_Polyhedron, Poly_Con_Relation, Poly_Gen_Relation, Generator, MIP_Problem, point as ppl_point
except ImportError:
    # old Sage
    from sage.libs.ppl import Variable, Constraint, Linear_Expression, Constraint_System, NNC_Polyhedron, Poly_Con_Relation, Poly_Gen_Relation, Generator, MIP_Problem, point as ppl_point

from sage.numerical.mip import MixedIntegerLinearProgram, MIPVariable, MIPSolverException

class BasicSemialgebraicSet_base(Element):

    """
    Abstract base class of mutable basic semialgebraic sets.
    """

    def __init__(self, ambient_dim):
        """
        Initialize a basic semialgebraic set as the universe in ``ambient_dim``.
        """
        self._ambient_dim = ambient_dim

    def ambient_dim(self):
        return self._ambient_dim

    @abstract_method
    def __contains__(self, point):
        """
        Whether the set contains the ``point`` (vector).
        """

    @abstract_method
    def closure_polyhedron(self):
        """
        If the topological closure is a polyhedron, return it in
        the form of a Sage polyhedron.  Otherwise raise an error.
        """

    def add_space_dimensions_and_embed(self, space_dim_to_add):
        self._ambient_dim += space_dim_to_add

    @abstract_method
    def add_linear_constraint(lhs, op):
        """
        ``lhs`` should be a vector of length ambient_dim.
        Add the constraint ``lhs`` * x ``op`` 0,
        where ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.
        """

    @abstract_method
    def is_linear_constraint_valid(lhs, op):
        """
        Whether the constraint ``lhs`` * x ``op`` 0
        is satisfied for all points of ``self``.
        """

## (1) In the first step, we implement the following class.  Everything is linear.
## Rewrite all direct uses of PPL in ParametricRealFieldElement, ParametricRealField
## using method calls to this class.

class BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(BasicSemialgebraicSet_base):

    """
    A (possibly half-open) polyhedral basic semialgebraic set,
    represented by a PPL ``NNC_Polyhedron``
    """

    def __init__(self, ambient_dim):
        super(BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron, self).__init__(ambient_dim)
        self._polyhedron = NNC_Polyhedron(ambient_dim, 'universe')

    def __copy__(self):
        # to implement
        raise NotImplementedError()

    # override the abstract methods


## (2) Then
class BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram(BasicSemialgebraicSet_base):

    """
    A (possibly half-open) polyhedral basic semialgebraic set,
    represented by a PPL ``NNC_Polyhedron``
    """

    def __init__(self, ambient_dim, solver=None):
        """
        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram(3)

        """
        super(BasicSemialgebraicSet_polyhedral_MixedIntegerLinearProgram, self).__init__(ambient_dim)
        self._mip = MixedIntegerLinearProgram(solver=solver)

    def __copy__(self):
        raise NotImplementedError()


## (3) Then introduce the following class to simplify the code in parametric.sage
class BasicSemialgebraicSet_eq_lt_le_sets(BasicSemialgebraicSet_base):

    """
    A basic semialgebraic set, represented in a straightforward way
    as 3 sets of constraints.
    """

    def __init__(self, ambient_dim):
        """
        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_eq_lt_le_sets(3)

        """
        super(BasicSemialgebraicSet_eq_lt_le_sets, self).__init__(ambient_dim)
        self._eq = set([])
        self._lt = set([])
        self._le = set([])

    def __copy__(self):
        raise NotImplementedError()

## (4) Later... introduce a class that takes care of the monomial lifting etc.

class BasicSemialgebraicSet_veronese(BasicSemialgebraicSet_base):

    """
    A basic semialgebraic set that delegates to another semialgebraic set
    via a Veronese embedding (RLT).
    """

    pass
