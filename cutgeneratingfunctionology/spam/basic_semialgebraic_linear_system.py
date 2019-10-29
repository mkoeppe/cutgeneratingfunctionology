r"""
Basic polyhedral semialgebraic sets represented as linear systems
"""

from __future__ import division, print_function, absolute_import

# Implement by rewriting code from formulations.sage on branch symbolic_FM. (FourierSystem, ...)

class BasicSemialgebraicSet_polyhedral_linear_system(BasicSemialgebraicSet_base):

    """
    A closed polyhedral basic semialgebraic set.

    In contrast to ``BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron``, it does not
    eagerly compute the double description, so it is suitable for large linear systems.
    Also it is suitable for arbitrary real fields as the ``base_ring``, such as ``ParametricRealField``.
    """

    def __init__(self, M=None, base_ring=None, ambient_dim=None):
        r"""
        Initialize a closed polyhedral basic semialgebraic set.

        (What exactly is M??)

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic_linear_system import *
            sage: ls = BasicSemialgebraicSet_polyhedral_linear_system(QQ, 2)

        """
        # Compute base_ring, ambient_dim from M if they are None but M is provided
        # ......
        super(BasicSemialgebraicSet_polyhedral_linear_system, self).__init__(base_ring, ambient_dim)

    def coordinate_projection(self, coordinates, bsa_class='linear_system')
        r"""
        Compute the projection to ``coordinates`` (a list or tuple of
        indices or variables of ``self.poly_ring``) as a new instance of
        ``BasicSemialgebraicSet_polyhedral_linear_system`` or the given
        ``bsa_class``.

        The projection is a set in the space of those coordinates.
        """
        raise NotImplementedError()
