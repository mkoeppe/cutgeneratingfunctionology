r"""
Basic polyhedral semialgebraic sets represented as linear systems
"""

from __future__ import division, print_function, absolute_import

from cutgeneratingfunctionology.spam.basic_semialgebraic import BasicSemialgebraicSet_base

# Implement by rewriting code from formulations.sage on branch symbolic_FM. (FourierSystem, ...)

class BasicSemialgebraicSet_polyhedral_linear_system(BasicSemialgebraicSet_base):

    """
    A closed polyhedral basic semialgebraic set.

    In contrast to ``BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron``, it does not
    eagerly compute the double description, so it is suitable for large linear systems.
    Also it is suitable for arbitrary real fields as the ``base_ring``, such as ``ParametricRealField``.
    """

    def __init__(self, base_ring=None, ambient_dim=None, poly_ring=None, eq=[], lt=[], le=[]):
        r"""
        Initialize a closed polyhedral basic semialgebraic set.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic_linear_system import *
            sage: ls = BasicSemialgebraicSet_polyhedral_linear_system(QQ, 2)

        """
        # Compute base_ring, ambient_dim from M if they are None but M is provided
        # ......
        if poly_ring is None:
            polys = list(eq) + list(lt) + list(le)
            if polys:
                poly_ring = polys[0].parent()
        if (ambient_dim is None) and (poly_ring is not None):
            ambient_dim = poly_ring.ngens()
            if base_ring is None and poly_ring is not None:
                base_ring = poly_ring.base_ring()
        if ambient_dim is None or base_ring is None:
            raise ValueError("if eq, lt, and le are all empty, must provide either poly_ring or both of base_ring and ambient_dim")
        if poly_ring is None:
            poly_ring = PolynomialRing(base_ring, "x", ambient_dim)
        super(BasicSemialgebraicSet_polyhedral_linear_system, self).__init__(base_ring, ambient_dim)
        self._poly_ring = poly_ring
        self._eq = set(eq)
        self._lt = set(lt)
        self._le = set(le)

    def one_step_elimination(self, coordinate, degree = 1, bsa_class='linear_system'):
        r"""
        Compute the projection by eliminating ``coordinates^degree``  as a new instance of
        ``BasicSemialgebraicSet_polyhedral_linear_system``.
        """
        new_eq=[]
        new_lt=[]
        new_le=[]
        sub=None
        for e in self._eq:
            if self._base_ring(e.coefficient({coordinate: degree})) != 0:
                sub = coordinate^degree - e/e.coefficient({coordinate: degree})
                break
        if sub is None:
            new_eq=self._eq
            lt_lower=[]
            lt_upper=[]
            le_lower=[]
            le_upper=[]
            
            for lt in self._lt:
                if self._base_ring(lt.coefficient({coordinate: degree}))>0:
                    lt_upper.append(lt)
                elif self._base_ring(lt.coefficient({coordinate: degree}))<0:
                    lt_lower.append(lt)
                else:
                    new_lt.append(lt)
            for le in self._le:
                if self._base_ring(le.coefficient({coordinate: degree}))>0:
                    le_upper.append(le)
                elif self._base_ring(le.coefficient({coordinate: degree}))<0:
                    le_lower.append(le)
                else:
                    new_le.append(le)

            # compute less than or equal to inequality
            for l in le_lower:
                for u in le_upper:
                    new_le.append(l*u.coefficient({coordinate: degree})-(u*l.coefficient({coordinate: degree})))

            # compute strictly less than inequality
            for l in le_lower:
                for u in lt_upper:
                    new_lt.append(l*u.coefficient({coordinate: degree})-(u*l.coefficient({coordinate: degree})))
            for l in lt_lower:
                for u in le_upper:
                    new_lt.append(l*u.coefficient({coordinate: degree})-(u*l.coefficient({coordinate: degree})))
            for l in lt_lower:
                for u in lt_upper:
                    new_lt.append(l*u.coefficient({coordinate: degree})-(u*l.coefficient({coordinate: degree})))
        else:
            for e in self._eq:
                new_eq.append(e+e.coefficient({coordinate: degree})*(sub-coordinate^degree))
            for lt in self._lt:
                new_lt.append(lt+lt.coefficient({coordinate: degree})*(sub-coordinate^degree))
            for le in self._le:
                new_le.append(le+le.coefficient({coordinate: degree})*(sub-coordinate^degree))
        
        return BasicSemialgebraicSet_polyhedral_linear_system(base_ring=self.base_ring, ambient_dim=self.ambient_dim()-1, poly_ring=self.poly_ring, eq=new_eq, lt=new_lt, le=new_le)
    
                                                          

    def coordinate_projection(self, coordinates, bsa_class='linear_system'):
        r"""
        Compute the projection to ``coordinates`` (a list or tuple of
        indices or variables of ``self.poly_ring``) as a new instance of
        ``BasicSemialgebraicSet_polyhedral_linear_system`` or the given
        ``bsa_class``.

        The projection is a set in the space of those coordinates.
        """
        raise NotImplementedError()

    def eq_poly(self):
        r"""
        Return a list, set, or iterator of the polynomials `f` in equations `f(x) = 0`
        in the description of ``self``.
        
        Together, ``eq_poly``, ``lt_poly``, and ``le_poly`` describe ``self``.
        """
        return self._eq

    def lt_poly(self):
        r"""
        Return a list, set, or iterator of the polynomials `f` in strict inequalities `f(x) < 0`
            in the description of ``self``.
            
        Together, ``eq_poly``, ``lt_poly``, and ``le_poly`` describe ``self``.
        """
        return self._lt

    def le_poly(self):
        r"""
        Return a list, set, or iterator of the polynomials `f` in inequalities `f(x) \leq 0`
        in the description of ``self``.
        
        Together, ``eq_poly``, ``lt_poly``, and ``le_poly`` describe ``self``.
        """
        return self._le
