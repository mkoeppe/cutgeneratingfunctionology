r"""
Basic polyhedral semialgebraic sets represented as linear systems
"""

from __future__ import division, print_function, absolute_import

from cutgeneratingfunctionology.spam.basic_semialgebraic import BasicSemialgebraicSet_base, BasicSemialgebraicSet_polyhedral
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.geometry.polyhedron.all import *
import sage.structure.element
import operator

cm = sage.structure.element.get_coercion_model()

# Implement by rewriting code from formulations.sage on branch symbolic_FM. (FourierSystem, ...)

class BasicSemialgebraicSet_polyhedral_linear_system(BasicSemialgebraicSet_polyhedral):

    """
    A closed polyhedral basic semialgebraic set.

    In contrast to ``BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron``, it does not
    eagerly compute the double description, so it is suitable for large linear systems.
    Also it is suitable for arbitrary real fields as the ``base_ring``, such as ``ParametricRealField``.
    """

    def __init__(self, base_ring=None, ambient_dim=None, poly_ring=None, eq=[], lt=[], le=[], history_set=None):
        r"""
        Initialize a closed polyhedral basic semialgebraic set.

        EXAMPLES::

            sage: import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *
            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic_linear_system import *
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: # One dimensional relu formulation.
            sage: K.<L,U,W,b>=ParametricRealField([QQ(-2),QQ(2),QQ(1),QQ(1/2)])
            sage: Q.<x0,x1,x,y,z>=K[]
            sage: le = [x0+x1-x, -x0-x1+x, -y, W*x0-b*z+b, W*x1-y+b*z, -W*x1+y-b*z, x0+U*z-U, -x0-L*z+L, x1-U*z, -x1+L*z]
            sage: bsa = BasicSemialgebraicSet_polyhedral_linear_system(poly_ring=Q, le=le)
            sage: bsa_eliminated = bsa.coordinate_projection([x0,x1])    # eliminate variable x0,x1.
            sage: bsa_eliminated.le_poly()
            {((L - U)~)*z,
            ((-L*W - b)~)*z + (L*W + b)~,
            ((-L + U)~)*z + (L - U)~,
            -y,
            -y + ((L*W + b)~)*z,
            y + ((-U*W - b)~)*z,
            -x + L~,
            -x + y + ((-L*W - b)~)*z + (L*W)~,
            x + (-U)~,
            x + ((-U*W - b)~)*z + b~,
            x - y + (W*b)~,
            x - y + ((U*W + b)~)*z + (-U*W)~}
            
        A non-parametric example::
            sage: D = polytopes.dodecahedron()
            sage: D.ambient_dim()
            3
            sage: bsa = BasicSemialgebraicSet_polyhedral_linear_system.from_polyhedron_to_linear_system(D)
            sage: set(D.vertices()) == set(bsa.to_polyhedron().vertices())
            True
            sage: # project out the second variable
            sage: proj_mat = matrix([[1,0,0],[0,0,1]])
            sage: D_proj = Polyhedron(vertices = (proj_mat*D.vertices_matrix()).columns())
            sage: bsa_proj = bsa.coordinate_projection([bsa.poly_ring().gens()[1]])
            sage: set(D_proj.vertices()) == set(bsa_proj.to_polyhedron().vertices())
            True
            
        Another non-parametric example::
            sage: b3 = polytopes.Birkhoff_polytope(3)
            sage: b3.ambient_dim()
            9
            sage: bsa = BasicSemialgebraicSet_polyhedral_linear_system.from_polyhedron_to_linear_system(b3)
            sage: set(b3.vertices()) == set (bsa.to_polyhedron().vertices())
            True
            sage: # project out x0, x2, x4, x6, x8 variable
            sage: proj_mat=matrix([[0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0]])
            sage: b3_proj = Polyhedron(vertices = (proj_mat*b3.vertices_matrix()).columns())
            sage: bsa_proj = bsa.coordinate_projection([bsa.poly_ring().gens()[i] for i in [0,2,4,6,8]])
            sage: set(b3_proj.vertices()) == set(bsa_proj.to_polyhedron().vertices())
            True
        """
        if poly_ring is None:
            raise ValueError("must specify the poly_ring")
        polys = list(eq) + list(lt) + list(le)
        #check if the polynomials are actually linear.
        for e in polys:
            if e.degree()>1:
                raise ValueError("only suitable for linear system.")
        if len(polys)>0:
            proper_poly_ring = cm.common_parent(*polys)
            if poly_ring != proper_poly_ring:
                raise ValueError("not a proper poly_ring.")
        if ambient_dim is None:
            ambient_dim = poly_ring.ngens()
        if base_ring is None:
            base_ring = poly_ring.base_ring()
        if history_set is None:
            history_set={}
            i=1
            for p in polys:
                history_set[p]={i}
                i+=1
        super(BasicSemialgebraicSet_polyhedral_linear_system, self).__init__(base_ring, ambient_dim)
        self._poly_ring = poly_ring
        self.history_set = history_set
        self._eq = set(eq)
        self._lt = set(lt)
        self._le = set(le)

    def remove_redundant_constant_polynomial(self):
        r"""
        Remove redundent constant polynomial or use constant polynomial to check naive infeasibility.
        If the linear system is infeasible, return a system with one inequality 1<=0.
        """
        temp_eq = self._eq.copy()
        temp_lt = self._lt.copy()
        temp_le = self._le.copy()
        for e in self._eq:
            if e.degree()<1:
                if self._base_ring(e) != self._base_ring(0):
                    #replace with an invalid inequality.
                    self._le={poly_ring(1)}
                    self._eq={}
                    self._lt={}
                    self.history_set={poly_ring(1):1}
                    return
                else:
                    temp_eq.remove(e)
                    if e in self.history_set:
                        del self.history_set[e]
        for lt in self._lt:
            if lt.degree()<1:
                if self._base_ring(lt) >= self._base_ring(0):
                    self._le={poly_ring(1)}
                    self._eq={}
                    self._lt={}
                    self.history_set={poly_ring(1):1}
                    return
                else:
                    temp_lt.remove(lt)
                    if lt in self.history_set:
                        del self.history_set[lt]
        for le in self._le:
            if le.degree()<1:
                if self._base_ring(le) > self._base_ring(0):
                    self._le={poly_ring(1)}
                    self._eq={}
                    self._lt={}
                    self.history_set={poly_ring(1):1}
                    return
                else:
                    temp_le.remove(le)
                    if le in self.history_set:
                        del self.history_set[le]
        self._le=temp_le
        self._eq=temp_eq
        self._lt=temp_lt

    def remove_redundancy_non_minimal_history_set(self):
        r"""
        Remove redundant (in)equality if the history set of the (in)equality is not minimal.
        """
        history_set = self.history_set.copy()
        for key1 in history_set.keys():
            for key2 in history_set.keys():
                if key1 == key2:
                    continue
                if history_set[key2].issubset(history_set[key1]) and history_set[key1] != history_set[key2]:
                    # remove key1
                    self.history_set.discard(key1)
                    self._le.discard(key1)
                    self._lt.discard(key1)
                    self._eq.discard(key1)
                    break

    def one_step_elimination(self, coordinate_index, bsa_class='linear_system'):
        r"""
        Compute the projection by eliminating ``coordinates``  as a new instance of
        ``BasicSemialgebraicSet_polyhedral_linear_system``.
        """
        # create a new poly_ring with one less generator (coordinate).
        if coordinate_index >= self._poly_ring.ngens():
            raise ValueError("doesn't exist the elimination variable")
        coordinate = self._poly_ring.gens()[coordinate_index]
        variables_names = [str(self._poly_ring.gens()[i]) for i in range(self._poly_ring.ngens()) if i != coordinate_index]
        new_poly_ring = PolynomialRing(self._base_ring, variables_names, len(variables_names))
        # create the ring hommorphism
        polynomial_map = [new_poly_ring.gens()[i] for i in range(new_poly_ring.ngens())]
        polynomial_map.insert(coordinate_index,new_poly_ring(0))
        
        new_eq=[]
        new_lt=[]
        new_le=[]
        new_history_set={}
        # try to find a substitution of coordinate in equalities.
        sub=None
        for e in self._eq:
            if e.monomial_coefficient(coordinate) != self.base_ring()(0):
                sub = coordinate - e/e.monomial_coefficient(coordinate)
                sub_history=self.history_set[e]
                break
        if sub is None:
            new_eq=self._eq
            for e in self._eq:
                new_history_set[e]=self.history_set[e]
            
            lt_lower=[]
            lt_upper=[]
            le_lower=[]
            le_upper=[]
            
            for lt in self._lt:
                if self._base_ring(lt.monomial_coefficient(coordinate))>self.base_ring()(0):
                    lt_upper.append(lt)
                elif self._base_ring(lt.monomial_coefficient(coordinate))<self.base_ring()(0):
                    lt_lower.append(lt)
                else:
                    new_lt.append(lt)
                    new_history_set[lt]=self.history_set[lt]
            for le in self._le:
                if self._base_ring(le.monomial_coefficient(coordinate))>self.base_ring()(0):
                    le_upper.append(le)
                elif self._base_ring(le.monomial_coefficient(coordinate))<self.base_ring()(0):
                    le_lower.append(le)
                else:
                    new_le.append(le)
                    new_history_set[le]=self.history_set[le]

            # compute less than or equal to inequality
            for l in le_lower:
                for u in le_upper:
                    polynomial=l*u.monomial_coefficient(coordinate)-(u*l.monomial_coefficient(coordinate))
                    new_le.append(polynomial)
                    new_history_set[polynomial]=self.history_set[l].union(self.history_set[u])

            # compute strictly less than inequality
            for l in le_lower:
                for u in lt_upper:
                    polynomial=l*u.monomial_coefficient(coordinate)-(u*l.monomial_coefficient(coordinate))
                    new_lt.append(polynomial)
                    new_history_set[polynomial]=self.history_set[l].union(self.history_set[u])
            for l in lt_lower:
                for u in le_upper:
                    polynomial=l*u.monomial_coefficient(coordinate)-(u*l.monomial_coefficient(coordinate))
                    new_lt.append(polynomial)
                    new_history_set[polynomial]=self.history_set[l].union(self.history_set[u])
            for l in lt_lower:
                for u in lt_upper:
                    polynomial=l*u.monomial_coefficient(coordinate)-(u*l.monomial_coefficient(coordinate))
                    new_lt.append(polynomial)
                    new_history_set[polynomial]=self.history_set[l].union(self.history_set[u])
        else:
            for e in self._eq:
                polynomial=e+e.monomial_coefficient(coordinate)*(sub-coordinate)
                new_eq.append(polynomial)
                if e.monomial_coefficient(coordinate) != self.base_ring()(0):
                    new_history_set[polynomial]=self.history_set[e].union(sub_history)
                else:
                    new_history_set[polynomial]=self.history_set[e]
            for lt in self._lt:
                polynomial=lt+lt.monomial_coefficient(coordinate)*(sub-coordinate)
                new_lt.append(polynomial)
                if lt.monomial_coefficient(coordinate) != self.base_ring()(0):
                    new_history_set[polynomial]=self.history_set[lt].union(sub_history)
                else:
                    new_history_set[polynomial]=self.history_set[lt]
            for le in self._le:
                polynomial=le+le.monomial_coefficient(coordinate)*(sub-coordinate)
                new_le.append(polynomial)
                if le.monomial_coefficient(coordinate) != self.base_ring()(0):
                    new_history_set[polynomial]=self.history_set[le].union(sub_history)
                else:
                    new_history_set[polynomial]=self.history_set[le]

        bsa = BasicSemialgebraicSet_polyhedral_linear_system(base_ring=self._base_ring, ambient_dim=self.ambient_dim(), poly_ring=self._poly_ring, eq=new_eq, lt=new_lt, le=new_le, history_set=new_history_set)
        new_bsa = bsa.section(polynomial_map,bsa_class=bsa_class,poly_ring=new_poly_ring,history_set=new_history_set)
        down_stair_history_set={}
        for key in new_bsa.history_set.keys():
            down_stair_history_set[key(polynomial_map)]=new_bsa.history_set[key]
        new_bsa.history_set=down_stair_history_set
        # remove constant polynomial
        new_bsa.remove_redundant_constant_polynomial()
        # remove polynomial with non_minimal_history_set
        new_bsa.remove_redundancy_non_minimal_history_set()
        return new_bsa

    def coordinate_projection(self, coordinates, bsa_class='linear_system'):
        r"""
        Compute the projection after projecting out the ``coordinates`` (a list or tuple of
        indices or variables of ``self.poly_ring``) as a new instance of
        ``BasicSemialgebraicSet_polyhedral_linear_system`` or the given ``bsa_class``.
        """
        res=self
        for c in coordinates:
            if not c in self._poly_ring.gens():
                raise ValueError("Coordinate not found in the polynomial ring")
        elimination_positions=[self._poly_ring.gens().index(c) for c in coordinates]
        for i in range(len(elimination_positions)):
            coordinate_index=elimination_positions[i]
            res = res.one_step_elimination(coordinate_index)
            # update indices
            for j in range(i+1,len(elimination_positions)):
                if elimination_positions[j]>coordinate_index:
                    elimination_positions[j]-=1
        return res
    
    def poly_ring(self):
        r"""
        Return the polynomial ring of self.
        """
        return self._poly_ring

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

    def add_polynomial_constraint(self, lhs, op):
        """
        ``lhs`` should be a polynomial.
        Add the constraint ``lhs``(x) ``op`` 0,
        where ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.
        """
        if lhs.parent() != self._poly_ring:
            try:
                lhs = self._poly_ring(lhs)
            except TypeError:
                raise TypeError("can not convert lhs into self.poly_ring")
        if lhs.degree() > 1:
            raise ValueError("{} is not a valid linear polynomial.".format(lhs))
        if op == operator.lt:
            self._lt.add(lhs)
        elif op == operator.gt:
            self._lt.add(-lhs)
        elif op == operator.eq:
            self._eq.add(lhs)
        elif op == operator.le:
            self._le.add(lhs)
        elif op == operator.ge:
            self._le.add(-lhs)
        else:
            raise ValueError("{} is not a supported operator".format(op))

    def add_linear_constraint(self, lhs_vector, cst, op):
        """
        Add the constraint ``lhs`` * x + cst ``op`` 0,
        where ``lhs`` is a vector of length ``ambient_dim`` and
        ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``.
        """
        if len(lhs_vector) != self._ambient_dim:
            raise ValueError("length of lhs_vector and ambient_dim do not match.")
        lhs=sum(lhs_vector[i]*self.poly_ring().gens()[i] for i in range(len(lhs_vector)))+cst
        self.add_polynomial_constraint(lhs, op)

    def to_polyhedron(self, **kwds):
        # not suitable for Paramatric field
        if len(self._lt)>0:
            raise ValueError("Contain strict inequalities.")
        ieqs=[]
        eqns=[]
        for le in self._le:
            temp_le = [le.constant_coefficient()] + [le.monomial_coefficient(v) for v in self.poly_ring().gens()]
            ieqs.append([-x for x in temp_le])
        for eq in self._eq:
            eqns.append([eq.constant_coefficient()] + [eq.monomial_coefficient(v) for v in self.poly_ring().gens()])
        return Polyhedron(ieqs=ieqs, eqns=eqns, **kwds)

    @classmethod
    def from_polyhedron_to_linear_system(cls, p, poly_ring=None):
        """
        Convert a instance of Polyhedron to a instance of BasicSemialgebraicSet_polyhedral_linear_system
        """
        if p.__class__ == cls:
            return p
        base_ring = p.base_ring().fraction_field()
        ambient_dim = p.ambient_dim()
        if poly_ring is None:
            poly_ring = PolynomialRing(base_ring, "x", ambient_dim)
        else:
            if ambient_dim != poly_ring.ambient_dim():
                raise ValueError("Ambient dimensions of p and poly_ring do not match.")
        self = cls(base_ring=base_ring, ambient_dim=ambient_dim, poly_ring=poly_ring)
        i=1
        for lhs in p.inequalities_list():
            self.add_linear_constraint(lhs[1:],lhs[0],operator.ge)
            self.history_set[-sum(lhs[1:][i]*self.poly_ring().gens()[i] for i in range(len(lhs)-1))-lhs[0]]={i}
            i+=1
        for lhs in p.equations_list():
            self.add_linear_constraint(lhs[1:],lhs[0],operator.eq)
            self.history_set[sum(lhs[1:][i]*self.poly_ring().gens()[i] for i in range(len(lhs)-1))+lhs[0]]={i}
            i+=1
        return self

