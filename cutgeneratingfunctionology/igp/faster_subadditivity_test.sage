from itertools import chain

import queue as queue
import itertools
import numpy as np
import gc

class SubadditivityTestTreeNode(object):

    r"""
    Class for the node in the spatial branch and bound tree for subadditivity testing.

    Each node ``N`` has the following attributes:

     - ......

     - ``N.affine_estimators = [[slope_I, intercept_I], [slope_J, intercept_J], [slope_K, intercept_K]]``
       describes affine linear underestimators for `\pi` on the restrictions to intervals `I` and `J`
       and an overestimator for `pi` on the restriction to the interval `K`.


    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = kzh_7_slope_1()
        sage: T = SubadditivityTestTree(h)
        sage: T.is_subadditive()
        True
        sage: # check affine_estimators are exact functions if the node is indivisible.
        sage: for node in T.leaf_set:
        ....:     if not node.is_divisible():
        ....:         break
        sage: slope_I, intercept_I = node.affine_estimators[0]
        sage: slope_I*node.I_bkpts()[0]+intercept_I == node.I_values()[0]
        True
        sage: slope_I*node.I_bkpts()[1]+intercept_I == node.I_values()[1]
        True
        sage: slope_J, intercept_J = node.affine_estimators[1]
        sage: slope_J*node.J_bkpts()[0]+intercept_J == node.J_values()[0]
        True
        sage: slope_J*node.J_bkpts()[1]+intercept_J == node.J_values()[1]
        True
        sage: slope_K, intercept_K = node.affine_estimators[2]
        sage: slope_K*node.K_bkpts()[0]+intercept_K == node.K_values()[0]
        True
        sage: slope_K*node.K_bkpts()[1]+intercept_K == node.K_values()[1]
        True
    """

    def __init__(self, fn, level, intervals, use_symmetry = False):

        self.intervals=tuple(tuple(I) for I in intervals)
        self.level=level
        self.function=fn
        self.use_symmetry=use_symmetry
        if use_symmetry:
            self.vertices=lower_triangle_vertices(verts(*intervals))
        else:
            self.vertices=verts(*intervals)
        self.projections=projections(self.vertices)

        # function values on the two endpoints
        self.I_end_points_values=[fn(self.projections[0][0]),fn(self.projections[0][1])]
        self.J_end_points_values=[fn(self.projections[1][0]),fn(self.projections[1][1])]
        self.K_end_points_values=[fn(fractional(self.projections[2][0])),fn(fractional(self.projections[2][1]))]

        # start and end index for breakpoints contained in the interval other than two endpoints.
        # for example, use fn.end_points()[i:j].
        self.I_bkpts_index=find_bkpts_index(fn.end_points(),self.projections[0])
        self.J_bkpts_index=find_bkpts_index(fn.end_points(),self.projections[1])
        self.K_bkpts_index=find_bkpts_index(fn.extended_end_points,self.projections[2])

        self.left_child=None
        self.right_child=None
        self.parent=None
        self.affine_estimators=None

    # SubadditivityTestTreeNode sort lexicographically by intervals.
    # We implement this so that PriorityQueue can be used with it in Py3.
    def __lt__(self, other):
        return self.intervals < other.intervals
    def __eq__(self, other):
        return self.intervals == other.intervals
    def __le__(self, other):
        return self.intervals <= other.intervals
    def __ge__(self, other):
        return self.intervals >= other.intervals
    def __gt__(self, other):
        return self.intervals > other.intervals
    def __ne__(self, other):
        return self.intervals != other.intervals
    def __hash__(self):
        return hash(self.intervals)

    def I_bkpts(self):
        if hasattr(self,'_I_bkpts'):
            return self._I_bkpts
        new_I=self.projections[0]
        self._I_bkpts = tuple(find_all_bkpts_in_the_interval(self.function.end_points(),new_I))
        return self._I_bkpts

    def J_bkpts(self):
        if hasattr(self,'_J_bkpts'):
            return self._J_bkpts
        new_J=self.projections[1]
        self._J_bkpts = tuple(find_all_bkpts_in_the_interval(self.function.end_points(),new_J))
        return self._J_bkpts

    def K_bkpts(self):
        if hasattr(self,'_K_bkpts'):
            return self._K_bkpts
        new_K=self.projections[2]
        self._K_bkpts = tuple(find_all_bkpts_in_the_interval(self.function.end_points(),new_K))
        return self._K_bkpts

    def I_values(self):
        if hasattr(self,'_I_values'):
            return self._I_values
        self._I_values = tuple(self.function(bkpt) for bkpt in self.I_bkpts())
        return self._I_values

    def J_values(self):
        if hasattr(self,'_J_values'):
            return self._J_values
        self._J_values = tuple(self.function(bkpt) for bkpt in self.J_bkpts())
        return self._J_values

    def K_values(self):
        if hasattr(self,'_K_values'):
            return self._K_values
        self._K_values = tuple(self.function(bkpt if bkpt<=1 else bkpt-1) for bkpt in self.K_bkpts())
        return self._K_values

    def I_slope_set(self):
        if hasattr(self,'_I_slope_set'):
            return self._I_slope_set
        self._I_slope_set = set(self.function.slope_values[self.I_bkpts_index[0]-1:self.I_bkpts_index[1]])
        return self._I_slope_set

    def J_slope_set(self):
        if hasattr(self,'_J_slope_set'):
            return self._J_slope_set
        self._J_slope_set = set(self.function.slope_values[self.J_bkpts_index[0]-1:self.J_bkpts_index[1]])
        return self._J_slope_set

    def K_slope_set(self):
        if hasattr(self,'_K_slope_set'):
            return self._K_slope_set
        self._K_slope_set = set(self.function.extended_slope_values[self.K_bkpts_index[0]-1:self.K_bkpts_index[1]])
        return self._K_slope_set

    def delta_pi_constant_lower_bound(self):
        r"""
        Compute the constant lower bound of delta pi in the current region defined by self.intervals. The bound is (min of fn on proj(I)) + (min of fn on proj(J)) - (max of fn on proj(K)).
        The second output is the slope and intercept values of three estimators.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)
            sage: h = kzh_7_slope_1()
            sage: T = SubadditivityTestTree(h)
            sage: T.is_subadditive()
            True
            sage: T.height
            16
            sage: N = T.root.left_child.right_child.left_child.right_child.left_child.right_child
            sage: N.intervals
            ((13/33, 16/33), (14/33, 1), (0, 1))
            sage: N.projections
            [[13/33, 16/33], [14/33, 20/33], [9/11, 1]]
            sage: lb, estimators = N.delta_pi_constant_lower_bound() # min([13/33, 16/33]) + min([14/33, 20/33]) - max([9/11, 1])
            sage: lb
            -9/22
            sage: estimators
            [[0, 13/66], [0, 13/66], [0, 53/66]]
        """
        alpha_I=min(self.I_end_points_values + self.function.values_at_end_points()[self.I_bkpts_index[0]:self.I_bkpts_index[1]])
        alpha_J=min(self.J_end_points_values + self.function.values_at_end_points()[self.J_bkpts_index[0]:self.J_bkpts_index[1]])
        beta_K=max(self.K_end_points_values + self.function.extended_values_at_end_points[self.K_bkpts_index[0]:self.K_bkpts_index[1]])
        return alpha_I+alpha_J-beta_K, [[0,alpha_I],[0,alpha_J],[0,beta_K]]

    def delta_pi_trivial_affine_lower_bound(self):
        r"""
        One heuristic affine lower bound of delta pi. Each affine estimator's slope is fixed and depends on the two endpoints of the corresponding interval.
        Time complexity is the same as constant lower bound.
        No guarantee that the bound is better than constant lower bound.
        The second output is the slope and intercept values of three estimators.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)
            sage: h = kzh_7_slope_1()
            sage: T = SubadditivityTestTree(h)
            sage: T.is_subadditive()
            True
            sage: N = T.root.left_child.right_child.left_child.right_child.left_child.right_child
            sage: lb, estimators = N.delta_pi_trivial_affine_lower_bound()
            sage: lb
            -13/18
            sage: estimators
            [[7/6, -1/3], [-19/12, 11/12], [-41/12, 35/9]]
        """
        slope_I=(self.I_end_points_values[0]-self.I_end_points_values[1])/(self.projections[0][0]-self.projections[0][1])
        slope_J=(self.J_end_points_values[0]-self.J_end_points_values[1])/(self.projections[1][0]-self.projections[1][1])
        slope_K=(self.K_end_points_values[0]-self.K_end_points_values[1])/(self.projections[2][0]-self.projections[2][1])
        return self.delta_pi_fast_affine_lower_bound(slope_I,slope_J,slope_K)

    def delta_pi_fast_lower_bound(self):
        r"""
        Choose the better bound between constant and trivial affine lower bounds.
        The second output is the slope and intercept values of three estimators.
        """
        constant_bound, constant_estimators = self.delta_pi_constant_lower_bound()
        trivial_affine_bound, affine_estimators = self.delta_pi_trivial_affine_lower_bound()
        if constant_bound>=trivial_affine_bound:
            return constant_bound, constant_estimators
        else:
            return trivial_affine_bound, affine_estimators

    def delta_pi_affine_lower_bound(self,solver='Coin'):
        r"""
        Compute the best lower bound of delta pi if using affine estimators, by solving an LP.
        The bound is guaranteed to be no worse than any other bounds in this code.
        The second output is the slope and intercept values of three estimators.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)
            sage: h = kzh_7_slope_1()
            sage: T = SubadditivityTestTree(h)
            sage: T.is_subadditive()
            True
            sage: N = T.root.left_child.right_child.left_child.right_child.left_child.right_child
            sage: lb, estimators = N.delta_pi_affine_lower_bound('PPL')
            sage: lb
            -1/3
            sage: estimators
            [[-5/2, 4/3], [-5/2, 4/3], [-5/2, 3]]
            sage: N.delta_pi_affine_lower_bound('PPL')[0] > N.delta_pi_fast_lower_bound()[0]
            True
        """

        p = MixedIntegerLinearProgram(maximization=True, solver=solver)
        v = p.new_variable()
        m1, b1, m2, b2, m3, b3, deltamin= v['m1'], v['b1'], v['m2'], v['b2'], v['m3'], v['b3'], v['deltamin']
        p.set_objective(deltamin)

        p.add_constraint(m1*self.projections[0][0]+b1<=self.I_end_points_values[0])
        p.add_constraint(m1*self.projections[0][1]+b1<=self.I_end_points_values[1])
        p.add_constraint(m2*self.projections[1][0]+b2<=self.J_end_points_values[0])
        p.add_constraint(m2*self.projections[1][1]+b2<=self.J_end_points_values[1])
        p.add_constraint(m3*self.projections[2][0]+b3>=self.K_end_points_values[0])
        p.add_constraint(m3*self.projections[2][1]+b3>=self.K_end_points_values[1])

        for i in range(self.I_bkpts_index[0],self.I_bkpts_index[1]):
            p.add_constraint(m1*self.function.end_points()[i]+b1<=self.function.values_at_end_points()[i])
        for j in range(self.J_bkpts_index[0],self.J_bkpts_index[1]):
            p.add_constraint(m2*self.function.end_points()[j]+b2<=self.function.values_at_end_points()[j])
        for k in range(self.K_bkpts_index[0],self.K_bkpts_index[1]):
            p.add_constraint(m3*self.function.extended_end_points[k]+b3>=self.function.extended_values_at_end_points[k])
        for v in self.vertices:
            x, y, z=v[0], v[1], v[0]+v[1]
            p.add_constraint(deltamin<=m1*x+b1+m2*y+b2-m3*z-b3)

        p.solve()
        # Due to floating point arithematic, the optimal solution is not exact. 
        # We round the solutions for slope values and recompute the intercepts, making sure the estimators are valid.
        m_I=QQ(p.get_values(m1))
        m_J=QQ(p.get_values(m2))
        m_K=QQ(p.get_values(m3))
        b_I=min([self.I_end_points_values[0] - m_I*self.projections[0][0], self.I_end_points_values[1] - m_I*self.projections[0][1]] + [self.function.values_at_end_points()[i]-m_I*self.function.end_points()[i] for i in range(self.I_bkpts_index[0],self.I_bkpts_index[1])])
        b_J=min([self.J_end_points_values[0] - m_J*self.projections[1][0], self.J_end_points_values[1] - m_J*self.projections[1][1]] + [self.function.values_at_end_points()[j]-m_J*self.function.end_points()[j] for j in range(self.J_bkpts_index[0],self.J_bkpts_index[1])])
        b_K=max([self.K_end_points_values[0] - m_K*self.projections[2][0], self.K_end_points_values[1] - m_K*self.projections[2][1]] + [self.function.extended_values_at_end_points[k]-m_K*self.function.extended_end_points[k] for k in range(self.K_bkpts_index[0],self.K_bkpts_index[1])])

        # it is possible that the bound is worse than fast_bound due to rounding, so we choose the better one.
        affine_bound = min(m_I*v[0]+b_I+m_J*v[1]+b_J-m_K*(v[0]+v[1])-b_K for v in self.vertices)
        fast_bound, fast_estimators = self.delta_pi_fast_lower_bound()
        if affine_bound >= fast_bound:
            return affine_bound, [[m_I,b_I],[m_J,b_J],[m_K,b_K]]
        else:
            return fast_bound, fast_estimators

    def delta_pi_fast_affine_lower_bound(self,slope_I,slope_J,slope_K):
        r"""
        If the slopes of the affine estimators are fixed, then compute the best intersepts first. Then compute the affine lower bound of delta pi.
        Note that if slope_I=slope_J=slope_K=0, then the bound is the constant bound.
        The second output is the slope and intercept values of three estimators.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)
            sage: h = kzh_7_slope_1()
            sage: T = SubadditivityTestTree(h)
            sage: T.is_subadditive()
            True
            sage: N = T.root.left_child.right_child.left_child.right_child.left_child.right_child
            sage: N.delta_pi_fast_affine_lower_bound(0,0,0) == N.delta_pi_constant_lower_bound()
            True
        """
        intercept_I=find_best_intercept([self.projections[0][0],self.projections[0][1]]+self.function.end_points()[self.I_bkpts_index[0]:self.I_bkpts_index[1]],self.I_end_points_values+self.function.values_at_end_points()[self.I_bkpts_index[0]:self.I_bkpts_index[1]],slope_I,lower_bound=True)
        intercept_J=find_best_intercept([self.projections[1][0],self.projections[1][1]]+self.function.end_points()[self.J_bkpts_index[0]:self.J_bkpts_index[1]],self.J_end_points_values+self.function.values_at_end_points()[self.J_bkpts_index[0]:self.J_bkpts_index[1]],slope_J,lower_bound=True)
        intercept_K=find_best_intercept([self.projections[2][0],self.projections[2][1]]+self.function.extended_end_points[self.K_bkpts_index[0]:self.K_bkpts_index[1]],self.K_end_points_values+self.function.extended_values_at_end_points[self.K_bkpts_index[0]:self.K_bkpts_index[1]],slope_K,lower_bound=False)
        lower_bound=min((slope_I*vertex[0]+intercept_I)+(slope_J*vertex[1]+intercept_J)-(slope_K*(vertex[0]+vertex[1])+intercept_K) for vertex in self.vertices)
        return lower_bound, [[slope_I,intercept_I],[slope_J,intercept_J],[slope_K,intercept_K]]

    def delta_pi_lower_bound(self,max_number_of_bkpts=0,solver='Coin'):
        r"""
        Strategic lower bound of delta pi. If the total number of bkpts in I,J,K is small, use affine bound by solving the LP. Use constant bound if max_number_of_bkpts=0. Otherwise use delta_pi_fast_lower_bound.
        """
        if hasattr(self,'_delta_pi_lb'):
            return self._delta_pi_lb
        if not self.is_divisible():
            return self.delta_pi_upper_bound()
        if max_number_of_bkpts==0:
            lower_bound, estimators=self.delta_pi_constant_lower_bound()
        elif self.I_bkpts_index[1]-self.I_bkpts_index[0] + self.J_bkpts_index[1]-self.J_bkpts_index[0] + self.K_bkpts_index[1]-self.K_bkpts_index[0] + 6 <=max_number_of_bkpts:
            lower_bound, estimators=self.delta_pi_affine_lower_bound(solver=solver)
        else:
            lower_bound, estimators=self.delta_pi_fast_lower_bound()
        self._delta_pi_lb=lower_bound
        # cache the three estimators.
        self.affine_estimators=estimators
        return lower_bound

    def delta_pi_upper_bound(self):
        r"""
        Compute the upper bound of delta pi based on the values of delta pi on the vertices of the region.
        """
        if hasattr(self,'_delta_pi_ub'):
            return self._delta_pi_ub
        self._delta_pi_ub=min(self.function(v[0])+self.function(v[1])-self.function(fractional(v[0]+v[1])) for v in self.vertices)
        return self._delta_pi_ub

    def branching_direction(self):
        new_I,new_J,new_K=self.projections
        lenI=new_I[1]-new_I[0]
        lenJ=new_J[1]-new_J[0]
        lenK=new_K[1]-new_K[0]
        candidates=[]
        if self.I_bkpts_index[0] != self.I_bkpts_index[1]:
            candidates.append(('I',lenI))
        if self.J_bkpts_index[0] != self.J_bkpts_index[1]:
            candidates.append(('J',lenJ))
        if self.K_bkpts_index[0] != self.K_bkpts_index[1]:
            candidates.append(('K',lenK))
        if not candidates:
            return None
        else:
            ind=np.argmax([candidate[1] for candidate in candidates])
            return candidates[ind][0]

    def new_intervals(self):
        dir=self.branching_direction()
        if dir=='I':
            new_bkpt=self.function.end_points()[(self.I_bkpts_index[0] + self.I_bkpts_index[1])//2]
            return (((self.intervals[0][0],new_bkpt),self.intervals[1],self.intervals[2]),
                    ((new_bkpt,self.intervals[0][1]),self.intervals[1],self.intervals[2]))
        elif dir=='J':
            new_bkpt=self.function.end_points()[(self.J_bkpts_index[0] + self.J_bkpts_index[1])//2]
            return ((self.intervals[0],(self.intervals[1][0],new_bkpt),self.intervals[2]),
                    (self.intervals[0],(new_bkpt,self.intervals[1][1]),self.intervals[2]))
        elif dir=='K':
            new_bkpt=self.function.extended_end_points[(self.K_bkpts_index[0] + self.K_bkpts_index[1])//2]
            return ((self.intervals[0],self.intervals[1],(self.intervals[2][0],new_bkpt)),
                    (self.intervals[0],self.intervals[1],(new_bkpt,self.intervals[2][1])))
        else:
            raise ValueError("Indivisible Region.")

    def is_divisible(self):
        if self.I_bkpts_index[0] == self.I_bkpts_index[1] and self.J_bkpts_index[0] == self.J_bkpts_index[1] and self.K_bkpts_index[0] == self.K_bkpts_index[1]:
            # if the node is indivisible, then affine_estimators equals the exact function.
            if self.affine_estimators is None:
                slope_I=(self.I_end_points_values[0]-self.I_end_points_values[1])/(self.projections[0][0]-self.projections[0][1])
                slope_J=(self.J_end_points_values[0]-self.J_end_points_values[1])/(self.projections[1][0]-self.projections[1][1])
                slope_K=(self.K_end_points_values[0]-self.K_end_points_values[1])/(self.projections[2][0]-self.projections[2][1])
                intercept_I=self.I_end_points_values[0]-slope_I*self.projections[0][0]
                intercept_J=self.J_end_points_values[0]-slope_J*self.projections[1][0]
                intercept_K=self.K_end_points_values[0]-slope_K*self.projections[2][0]
                self.affine_estimators=[[slope_I,intercept_I],[slope_J,intercept_J],[slope_K,intercept_K]]
            return False
        else:
            return True

    def generate_children(self,upper_bound=0,stop_only_if_strict=True,**kwds):
        if not self.is_fathomed(upper_bound,stop_only_if_strict,**kwds):
            I1,I2=self.new_intervals()
            self.left_child=SubadditivityTestTreeNode(self.function,self.level+1,I1,self.use_symmetry)
            self.right_child=SubadditivityTestTreeNode(self.function,self.level+1,I2,self.use_symmetry)
            self.left_child.parent=self
            self.right_child.parent=self

    def is_fathomed(self,value=0,stop_only_if_strict=True,**kwds):
        if self.is_divisible():
            if stop_only_if_strict:
                if self.delta_pi_lower_bound(**kwds)>value:
                    return True
                else:
                    return False
            else:
                # stop branching early.
                if self.delta_pi_lower_bound(**kwds)>=value:
                    return True
                else:
                    return False
        else:
            return True

    def verify_estimators(self):
        r"""
        If affine estimators have been computed, then verirify that it is indeed lower/upper estimator.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)
            sage: h = kzh_7_slope_1()
            sage: T = SubadditivityTestTree(h)
            sage: T.is_subadditive()
            True
            sage: all(node.verify_estimators() for node in T.complete_node_set)
            True
        """
        if self.affine_estimators is None:
            raise ValueError("Affine estimators have not been computed.")

        for i in range(len(self.I_bkpts())):
            bkpt = self.I_bkpts()[i]
            value = self.I_values()[i]
            slope_I, intercept_I = self.affine_estimators[0]
            if slope_I*bkpt+intercept_I>value:
                return False

        for j in range(len(self.J_bkpts())):
            bkpt = self.J_bkpts()[j]
            value = self.J_values()[j]
            slope_J, intercept_J = self.affine_estimators[1]
            if slope_J*bkpt+intercept_J>value:
                return False

        for k in range(len(self.K_bkpts())):
            bkpt = self.K_bkpts()[k]
            value = self.K_values()[k]
            slope_K, intercept_K = self.affine_estimators[2]
            if slope_K*bkpt+intercept_K<value:
                return False

        return True

    def verify_vertices(self,**kwds):
        r"""
        If the lower bound of delta pi has been computed, verify delta pi value is no smaller than lower bound.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)
            sage: h = kzh_7_slope_1()
            sage: T = SubadditivityTestTree(h)
            sage: T.is_subadditive()
            True
            sage: all(node.verify_vertices() for node in T.complete_node_set)
            True
        """
        for v in self.vertices:
            delta_pi = self.function(v[0]) + self.function(v[1]) - self.function(fractional(v[0]+v[1]))
            if delta_pi<self.delta_pi_lower_bound(**kwds):
                return False
        return True

    def plot(self, colorful=False, **bound_kwds):
        v=[ver for ver in self.vertices]
        region=Polyhedron(vertices=v)
        p = region.projection().render_outline_2d()
        if colorful:
            if not self.is_divisible():
                p+=Face(self.intervals).plot(rgbcolor="yellow", fill_color="yellow")
            elif self.delta_pi_lower_bound(**bound_kwds)>0:
                p+=Face(self.intervals).plot(rgbcolor="red", fill_color="red")
        return p

class SubadditivityTestTree:

    r"""
    Class for spatial branch and bound for subadditivity testing.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = not_minimal_1()
        sage: T = SubadditivityTestTree(h, objective_limit=-1000)
        sage: T.minimum(max_number_of_bkpts=0, search_method='BB')
        -1/5
        sage: T = SubadditivityTestTree(h)
        sage: T.is_subadditive()
        False

        sage: h = kzh_7_slope_1()
        sage: T1 = SubadditivityTestTree(h, use_symmetry = False)
        sage: T1.is_subadditive()
        True
        sage: T1.height, len(T1.complete_node_set)
        (16, 1447)

        sage: T2 = SubadditivityTestTree(h, use_symmetry = True)
        sage: T2.is_subadditive()
        True
        sage: T2.height, len(T2.complete_node_set)
        (13, 749)

    """

    def __init__(self,fn,intervals=((0,1), (0,1), (0,2)),global_upper_bound=0,objective_limit=0, use_symmetry = False):
        fn.extended_end_points=fn.end_points() + [bkpt+1 for bkpt in fn.end_points() if bkpt != 0]
        fn.extended_values_at_end_points=fn.values_at_end_points() + fn.values_at_end_points()[1:]
        fn.slope_values=[(fn.values_at_end_points()[i+1]-fn.values_at_end_points()[i])/(fn.end_points()[i+1]-fn.end_points()[i]) for i in range(len(fn.end_points())-1)]
        fn.extended_slope_values=fn.slope_values+fn.slope_values
        self.function=fn
        self.intervals=intervals
        self.global_upper_bound=global_upper_bound
        self.objective_limit=objective_limit
        self.root=SubadditivityTestTreeNode(fn,0,intervals,use_symmetry)
        self.height=0
        self.complete_node_set=set([self.root])
        self.leaf_set=set([self.root])
        # the order of unfathomed nodes matters, like DFS or BFS
        self.unfathomed_node_list=queue.PriorityQueue()
        self.nonsubadditive_vertices=set()
        self.additive_vertices=set()
        self.maximal_additive_faces=set()

    def number_of_nodes(self):
        return len(self.complete_node_set)

    def number_of_leaves(self):
        return len(self.leaf_set)

    def node_branching(self,node,search_method='BB',find_min=True,stop_only_if_strict=True,**kwds):
        r"""
        Branch on a given node, dependent on current self.global_upper_bound.
        """
        if not node.left_child:
            if find_min:
                node.generate_children(self.global_upper_bound,stop_only_if_strict,**kwds)
            else:
                node.generate_children(self.objective_limit,stop_only_if_strict,**kwds)
            if node.left_child:
                if search_method=='BFS':
                    self.unfathomed_node_list.put((len(self.complete_node_set),node.left_child))
                    self.unfathomed_node_list.put((len(self.complete_node_set)+1,node.right_child))
                elif search_method=='DFS':
                    self.unfathomed_node_list.put((-len(self.complete_node_set),node.left_child))
                    self.unfathomed_node_list.put((-len(self.complete_node_set)-1,node.right_child))
                elif search_method=='BB':
                    if node.left_child.is_divisible():
                        self.unfathomed_node_list.put((node.left_child.delta_pi_lower_bound(**kwds),node.left_child))
                    else:
                        self.unfathomed_node_list.put((node.left_child.delta_pi_upper_bound(),node.left_child))
                    if node.right_child.is_divisible():
                        self.unfathomed_node_list.put((node.right_child.delta_pi_lower_bound(**kwds),node.right_child))
                    else:
                        self.unfathomed_node_list.put((node.right_child.delta_pi_upper_bound(),node.right_child))
                else:
                    raise ValueError("Can't recognize search_method.")
                self.height=max(self.height,node.left_child.level)
                self.complete_node_set.update({node.left_child,node.right_child})
                self.leaf_set.discard(node)
                self.leaf_set.update({node.left_child,node.right_child})

    def next_level(self, node_set,search_method='BB',find_min=True,stop_only_if_strict=True,**kwds):
        r"""
        Generate nodes in the next level.
        """
        next_level=set()
        for node in node_set:
            self.node_branching(node,search_method=search_method,find_min=find_min,stop_only_if_strict=stop_only_if_strict,**kwds)
            if node.left_child:
                next_level.update({node.left_child,node.right_child})
        return next_level

    def is_subadditive(self,stop_if_fail=False,cache_additive_vertices=False,search_method='BB',**kwds):
        r"""
        Check whether the mininum of delta pi is no smaller than objective_limit.
        If the parameter cache_additive_vertices is True, then we cache vertices where delta pi equals objective_limit.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)
            sage: h = not_minimal_1()
            sage: T = SubadditivityTestTree(h)
            sage: T.minimum()
            -1/5
            sage: T = SubadditivityTestTree(h, objective_limit=-1/5)
            sage: T.is_subadditive()
            True
            sage: T = SubadditivityTestTree(h, objective_limit=0)
            sage: T.is_subadditive()
            False

            sage: T = SubadditivityTestTree(h, objective_limit=-1/5)
            sage: T.is_subadditive(cache_additive_vertices=True)
            True
            sage: T.additive_vertices
            {(1/5, 1/5)}
            sage: h(1/5)+h(1/5)-h(2/5) == T.objective_limit
            True
        """
        if search_method=='BFS' or search_method=='DFS':
            self.unfathomed_node_list.put((0,self.root))
        elif search_method=='BB':
            self.unfathomed_node_list.put((self.root.delta_pi_lower_bound(**kwds),self.root))
        else:
            raise ValueError("Can't recognize search_method.")
        while not self.unfathomed_node_list.empty():
            current_node=self.unfathomed_node_list.get()[1]
            if len(current_node.vertices)==0:
                continue
            for v in current_node.vertices:
                delta=self.function(v[0])+self.function(v[1])-self.function(fractional(v[0]+v[1]))
                if delta<self.objective_limit:
                    self._is_subadditive=False
                    # can stop early
                    if stop_if_fail:
                        return False
                    self.nonsubadditive_vertices.add(v)
                if delta==self.objective_limit and cache_additive_vertices:
                    self.additive_vertices.add(v)
            self.node_branching(current_node,search_method=search_method,find_min=False,stop_only_if_strict=cache_additive_vertices,**kwds)
        if not hasattr(self,'_is_subadditive'):
            self._is_subadditive=True
        return self._is_subadditive

    def minimum(self,search_method='BB',**kwds):
        if search_method=='BFS' or search_method=='DFS':
            self.unfathomed_node_list.put((0,self.root))
        elif search_method=='BB':
            self.unfathomed_node_list.put((self.root.delta_pi_lower_bound(**kwds),self.root))
        else:
            raise ValueError("Can't recognize search_method.")
        while not self.unfathomed_node_list.empty():
            current_node=self.unfathomed_node_list.get()[1]
            if len(current_node.vertices)==0:
                continue
            upper_bound=current_node.delta_pi_upper_bound()
            if upper_bound<self.global_upper_bound:
                self.global_upper_bound=upper_bound
            self.node_branching(current_node,search_method=search_method,find_min=True,stop_only_if_strict=False,**kwds)
        self.min=self.global_upper_bound
        return self.min

    def generate_maximal_additive_faces(self,search_method='BB',**kwds):
        epsilon=QQ(1)/100000000
        if search_method=='BFS' or search_method=='DFS':
            self.unfathomed_node_list.put((0,self.root))
        elif search_method=='BB':
            self.unfathomed_node_list.put((self.root.delta_pi_lower_bound(**kwds),self.root))
        else:
            raise ValueError("Can't recognize search_method.")
        while not self.unfathomed_node_list.empty():
            current_node=self.unfathomed_node_list.get()[1]
            if len(current_node.vertices)==0:
                continue
            # prune if no additive set exists based on slope value set.
            if current_node.I_slope_set().isdisjoint(current_node.J_slope_set()) and current_node.I_slope_set().isdisjoint(current_node.K_slope_set()) and current_node.J_slope_set().isdisjoint(current_node.K_slope_set()):
                continue
            if current_node.is_divisible():
                self.node_branching(current_node,search_method=search_method,find_min=False,stop_only_if_strict=True,**kwds)
            else:
                temp = []
                keep = False
                for vertex in current_node.vertices:
                    if delta_pi(self.function, vertex[0],vertex[1]) == 0:
                        temp.append(vertex)
                        keep = True
                if len(temp) == 2:
                    if temp[0][0] == temp[1][0]:
                        if temp[1][0] == current_node.projections[0][0]:
                            if delta_pi(self.function, temp[0][0]-epsilon,(temp[0][1]+temp[1][1])/2)==0:
                                keep = False
                        else:
                            keep = False
                    elif temp[0][1] == temp[1][1]:
                        if temp[1][1] == current_node.projections[1][0]:
                            if delta_pi(self.function, (temp[0][0]+temp[1][0])/2,temp[0][1]-epsilon)==0:
                                keep = False
                        else:
                            keep = False
                    elif temp[0][0] + temp[0][1] == temp[1][0] + temp[1][1]:
                        if temp[1][0] + temp[1][1] == current_node.projections[2][0]:
                            if delta_pi(self.function, (temp[0][0]+temp[1][0])/2-epsilon,(temp[0][1]+temp[1][1])/2-epsilon)==0:
                                keep = False
                        else:
                            keep = False
                    else:
                        keep = False
                elif len(temp) == 1:
                    x, y = temp[0]
                    if delta_pi(self.function,x+epsilon,y)==0 or delta_pi(self.function,x-epsilon,y)==0 or delta_pi(self.function,x,y+epsilon)==0 or delta_pi(self.function,x,y-epsilon)==0 or delta_pi(self.function,x+epsilon,y-epsilon)==0 or delta_pi(self.function,x-epsilon,y+epsilon)==0:
                        keep = False
                if keep:
                    trip = projections(temp)
                    self.maximal_additive_faces.add(Face(trip, vertices=temp, is_known_to_be_minimal=True))
        return self.maximal_additive_faces

    def generate_covered_components_big_cell(self,search_method='BB',**kwds):
        r"""
        Generate covered components stratigically using big cell.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)
            sage: h = kzh_7_slope_1()
            sage: T = SubadditivityTestTree(h,use_symmetry = True)
            sage: covered_components_1 = T.generate_covered_components_big_cell()
            sage: covered_components_2 = generate_covered_components(h)
            sage: covered_components_3 = generate_covered_components_strategically(h)
            sage: set(tuple(c) for c in covered_components_1) == set(tuple(c) for c in covered_components_2)
            False
            sage: set(tuple(c) for c in covered_components_1) == set(tuple(c) for c in covered_components_3)
            True
        """
        if hasattr(self,'covered_components'):
            return self.covered_components
        step = 1
        epsilon=QQ(1)/100000000
        additive_edges = []
        # covered_components will always be mutually exclusive throughout the algorithm.
        covered_components = []
        if search_method=='BFS' or search_method=='DFS':
            self.unfathomed_node_list.put((0,self.root))
        elif search_method=='BB':
            self.unfathomed_node_list.put((self.root.delta_pi_lower_bound(**kwds),self.root))
        else:
            raise ValueError("Can't recognize search_method.")
        while not self.unfathomed_node_list.empty():
            current_node=self.unfathomed_node_list.get()[1]
            if len(current_node.vertices)==0:
                continue
            # prune if no additive set exists based on slope value set.
            if current_node.I_slope_set().isdisjoint(current_node.J_slope_set()) and current_node.I_slope_set().isdisjoint(current_node.K_slope_set()) and current_node.J_slope_set().isdisjoint(current_node.K_slope_set()):
                continue
            if current_node.is_divisible():
                self.node_branching(current_node,search_method=search_method,find_min=False,stop_only_if_strict=True,**kwds)
            else:
                temp = []
                for vertex in current_node.vertices:
                    if delta_pi(self.function, vertex[0],vertex[1]) == 0:
                        temp.append(vertex)
                if len(temp) == 2:
                    # possible additive edges
                    if temp[0][0] == temp[1][0] and temp[1][0] == current_node.projections[0][0] and delta_pi(self.function, temp[0][0]-epsilon,(temp[0][1]+temp[1][1])/2) != 0:
                        # vertical additive edge.
                        additive_edges.append(Face(projections(temp), vertices=temp, is_known_to_be_minimal=True))
                    elif temp[0][1] == temp[1][1] and temp[1][1] == current_node.projections[1][0] and delta_pi(self.function, (temp[0][0]+temp[1][0])/2,temp[0][1]-epsilon) != 0:
                        # horizontal additive edge.
                        additive_edges.append(Face(projections(temp), vertices=temp, is_known_to_be_minimal=True))
                    elif temp[0][0] + temp[0][1] == temp[1][0] + temp[1][1] and temp[1][0] + temp[1][1] == current_node.projections[1][0] and delta_pi(self.function, (temp[0][0]+temp[1][0])/2-epsilon,(temp[0][1]+temp[1][1])/2-epsilon) != 0:
                        # diagonal additive edge.
                        additive_edges.append(Face(projections(temp), vertices=temp, is_known_to_be_minimal=True))
                elif len(temp) >=3:
                    # 2D additive faces
                    I,J,K = current_node.projections
                    K_mod_1 = interval_mod_1(K)
                    component = union_of_coho_intervals_minus_union_of_coho_intervals([[open_interval(* I)], [open_interval(* J)], [open_interval(* K_mod_1)]],[])
                    new_component, remaining_components = merge_components_with_given_component_strategically(component, covered_components)
                    if new_component:
                        covered_components = remaining_components + [new_component]
                        if logging.getLogger().isEnabledFor(logging.DEBUG):
                            if new_component == component:
                                logging.debug("Step %s: Consider the 2d additive %s.\n%s is directly covered." % (step, Face(current_node.projections, vertices=current_node.vertices, is_known_to_be_minimal=True), component))
                            else:
                                logging.debug("Step %s: By merging components that overlap with projections of the 2d additive %s, we obtain a larger covered component %s" % (step, Face(current_node.projections, vertices=current_node.vertices, is_known_to_be_minimal=True), new_component))
                            step += 1
        for edge in additive_edges:
            fdm = edge.functional_directed_move()
            sym_fdm = [fdm]
            if fdm.sign() == 1:
                backward_fdm = FunctionalDirectedMove(fdm.range_intervals(), (1, -fdm[1]))
                sym_fdm.append(backward_fdm)
            for covered_component in covered_components:
                component = []
                for fdm in sym_fdm:
                    overlapped_ints = list(intersection_of_coho_intervals([covered_component, fdm.intervals()]))
                    moved_intervals = [[fdm.apply_to_coho_interval(overlapped_int)] for overlapped_int in overlapped_ints ]
                    component = union_of_coho_intervals_minus_union_of_coho_intervals(moved_intervals + [overlapped_ints] + [component], [])
                if component:
                    new_component, remaining_components = merge_components_with_given_component_strategically(component, covered_components)
                    if new_component:
                        if logging.getLogger().isEnabledFor(logging.DEBUG):
                            if len(remaining_components) == len(covered_components)-1:
                                # only extend one component
                                newly_covered = union_of_coho_intervals_minus_union_of_coho_intervals([new_component], covered_components)
                                logging.debug("Step %s: %s is indirectly covered." % (step, newly_covered))
                            else:
                                # extend and merge
                                logging.debug("Step %s: By merging components that are connected by the 1d additive %s, we obtain a larger covered component %s" % (step, edge, new_component))
                            step += 1
                        covered_components = remaining_components + [new_component]
        self.covered_components = covered_components
        self.proof_step = step - 1
        return self.covered_components

    def plot_current_regions(self,colorful=False,**kwds):
        p=Graphics()
        legend_kwds = { 'legend_label1': "indivisible face" , 'legend_label2': "strict subadditive divisible face"}
        legend=[0,0]
        for node in self.leaf_set:
            v=[ver for ver in node.vertices]
            region=Polyhedron(vertices=v)
            p+=region.projection().render_outline_2d()
            if colorful:
                if not node.is_divisible():
                    if legend[0]==0:
                        p+=Face(node.intervals).plot(rgbcolor = "yellow",fill_color = "yellow",legend_label=legend_kwds['legend_label1'])
                        legend[0]=1
                    else:
                        p+=Face(node.intervals).plot(rgbcolor = "yellow",fill_color = "yellow")
                elif node.delta_pi_lower_bound(**kwds)>self.global_upper_bound:
                    if legend[1]==0:
                        p+=Face(node.intervals).plot(rgbcolor = "red",fill_color = "red",legend_label=legend_kwds['legend_label2'])
                        legend[1]=1
                    else:
                        p+=Face(node.intervals).plot(rgbcolor = "red",fill_color = "red")
        return p

def test_maximal_additive_faces(fn,**kwds):
    r"""
    Check the difference of maximal additive faces generated by naive and faster algorithms.
    """
    naive=set(generate_maximal_additive_faces(fn))
    faster=SubadditivityTestTree(fn).generate_maximal_additive_faces(**kwds)
    return naive.difference(faster),faster.difference(naive)

def find_all_bkpts_in_the_interval(bkpts,interval):
    r"""
    Return a list of breakpoints contained in interval and bkpts+Z. bkpts=fn.end_points(), and interval can be any closed interval contained in [0,2].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: fn=gj_2_slope()
        sage: bkpts=fn.end_points()
        sage: bkpts
        [0, 4/15, 1/3, 3/5, 1]
        sage: find_all_bkpts_in_the_interval(bkpts,[4/15,3/5])
        [4/15, 1/3, 3/5]
        sage: find_all_bkpts_in_the_interval(bkpts,[1/15,4/5])
        [1/15, 4/15, 1/3, 3/5, 4/5]
        sage: find_all_bkpts_in_the_interval(bkpts,[1/3,2])
        [1/3, 3/5, 1, 19/15, 4/3, 8/5, 2]
        sage: find_all_bkpts_in_the_interval(bkpts,[9/8,19/10])
        [9/8, 19/15, 4/3, 8/5, 19/10]
    """
    if interval[1]<=1:
        i,j=find_bkpts_index(bkpts,interval)
        return [interval[0]]+bkpts[i:j]+[interval[1]]
    if interval[0]>=1:
        new_interval=[interval[0]-1,interval[1]-1]
        i,j=find_bkpts_index(bkpts,new_interval)
        return [interval[0]]+[bkpt+1 for bkpt in bkpts[i:j]]+[interval[1]]
    interval1=[interval[0],1]
    interval2=[0,interval[1]-1]
    i1,j1=find_bkpts_index(bkpts,interval1)
    i2,j2=find_bkpts_index(bkpts,interval2)
    return [interval[0]]+bkpts[i1:j1+1]+[bkpt+1 for bkpt in bkpts[i2:j2]]+[interval[1]]

def find_bkpts_index(bkpts,interval):
    i=bisect_left(bkpts, interval[0])
    j=bisect_left(bkpts, interval[1])
    if interval[0]==bkpts[i]:
        i=i+1
    return i,j

def plot_2d_regions(fn,colorful=False,search_method='BB',find_min=True,stop_only_if_strict=True, show_plots=True, **kwds):
    # FIXME: This should really be integrated to the class via 'show_plots'.
    r"""
    Visualize the subadditivity test based on spatial branch and bound.

    Cells of each diagrams correspond to leaves of the branch and bound tree.

    If ``colorful=True``, the plot is done in fried eggs and tomato style:
    the yellow regions are pruned by indivisibility (i.e.,
    they are faces of `\Delta\mathcal{P}`), and the red regions are pruned by bounds.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = kzh_7_slope_1()
        sage: plot_2d_regions(h, colorful=True)      # not tested
    """
    T=SubadditivityTestTree(fn)
    p=Graphics()
    # current_level=T.complete_node_set results in error: Set changed size during iteration
    current_level=set([T.root])
    level_index = 0
    while current_level:
        p+=plot_2d_regions_in_one_level(current_level,colorful=colorful)   ## FIXME: Should this be passing kwds?
        show_plot(p, show_plots=show_plots, tag='bb-{}'.format(level_index))
        next_level=T.next_level(current_level,search_method=search_method,find_min=find_min,stop_only_if_strict=stop_only_if_strict,**kwds)
        level_index += 1
        current_level=next_level

def plot_2d_regions_in_one_level(node_set,colorful=False,**bound_kwds):
    if node_set:
        return sum(node.plot(colorful=colorful, **bound_kwds)
                   for node in node_set)
    else:
        return Graphics()

def values_of_delta_pi_over_grid(fn,q):
    r"""
    Return a matrix representing the values of delta pi over the grid (1/q)Z*(1/q)Z.
    """
    res=np.zeros((q+1,q+1))
    for i in range(q+1):
        for j in range(q+1):
            res[q-j][i]=delta_pi(fn,i/q,j/q)
    return res

def strategic_delta_pi_min(fn,approximation='constant',norm='one',branching_point_selection='median'):
    r"""
    Return the minimum of delta_pi.
    """
    f=find_f(fn)
    a=delta_pi_min(fn,[0,f],[0,f],[0,f],approximation=approximation,norm=norm,branching_point_selection=branching_point_selection)
    b=delta_pi_min(fn,[0,1],[0,1],[f,1+f],approximation=approximation,norm=norm,branching_point_selection=branching_point_selection)
    c=delta_pi_min(fn,[f,1],[f,1],[1+f,2],approximation=approximation,norm=norm,branching_point_selection=branching_point_selection)
    return min(a,b,c)
 
def find_slope_intercept_trivial(X,Y):
    m=(Y[1]-Y[0])/(X[1]-X[0])
    b=Y[1]-m*X[1] 
    return m,b

def find_best_intercept(X,Y,slope,lower_bound=True):
    r"""
    Find the intercept of the affine lower/upper estimator with fixed slope.
    """
    min_b=infinity
    max_b=-infinity
    for i in range(len(X)):
        b=Y[i]-X[i]*slope
        if b<min_b:
            min_b=b
        if b>max_b:
            max_b=b
    if lower_bound:
        return min_b
    else:
        return max_b

                           
def find_best_slope_intercept(X,Y,lower_bound=True,solver='Coin',norm='one'):
    r"""
    Find the slope and intercept of the affine lower/upper estimator.
    """
    p = MixedIntegerLinearProgram(maximization=False, solver=solver)
    v = p.new_variable()
    if norm=='one':
        m, b= v['m'], v['b']
        if lower_bound:
            p.set_objective(-m*sum(X)-b*len(X))
            for i in range(len(X)):
                p.add_constraint(m*X[i]+b<=Y[i])
        else:
            p.set_objective(m*sum(X)+b*len(X))
            for i in range(len(X)):
                p.add_constraint(m*X[i]+b>=Y[i])
    elif norm=='inf':
        m, b, z=v['m'], v['b'], v['z']
        p.set_objective(z)
        if lower_bound:
            for i in range(len(X)):
                p.add_constraint(z>=Y[i]-m*X[i]-b>=0)
        else:
            for i in range(len(X)):
                p.add_constraint(z>=m*X[i]+b-Y[i]>=0)
    else:
        raise ValueError("Can't recognize norm.")
    p.solve()
    return p.get_values(m),p.get_values(b)

def find_affine_bounds_new(fn,X,Y,lower_bound=True,norm='one'):
    if len(X)==2:
        slope,intercept=find_slope_intercept_trivial(X,Y)
    else:
        slope,intercept=find_best_slope_intercept(X,Y,lower_bound=lower_bound,norm=norm)
    return slope, intercept


def find_affine_bounds(fn,I,lower_bound=True,extended=False,norm='one'):
    r"""
    Find the affine lower/upper estimator of the function fn over the closed interval I.
    """
    if not extended:
        bkpts=fn.end_points()
    else:
        bkpts=fn.end_points()[:-1]+[1+bkpt for bkpt in fn.end_points()]
    I_index=find_possible_branching_bkpts_index(bkpts,I)
    X,Y=[I[0],I[1]],[fn(fractional(I[0])),fn(fractional(I[1]))]
    if I_index:
        for i in range(I_index[0], I_index[1]+1):
            X.append(bkpts[i])
            Y.append(fn(fractional(bkpts[i])))
        slope,intercept=find_best_slope_intercept(X,Y,lower_bound=lower_bound,norm=norm)
    else:
        slope,intercept=find_slope_intercept_trivial(X,Y)    
    return slope, intercept

def find_possible_branching_bkpts_index(bkpts,I):
    r"""
    Return [i,j] such that I[0]<bkpts[k]<I[1] for any i<=k<=j. Return None if empty.
    """
    i=bisect_left(bkpts, I[0])
    j=bisect_left(bkpts, I[1])-1
    if I[0]==bkpts[i]:
        i=i+1
    if i>j:
        return None
    else:
        return [i,j]

def fn_constant_bounds(fn,I,lower_bound=True,extended=False):
    r"""
    Find the constant lower/upper bound of the function fn over the closed interval I.
    """
    if not extended:
        bkpts=fn.end_points()
    else:
        bkpts=fn.end_points()[:-1]+[1+bkpt for bkpt in fn.end_points()]
    I_index=find_possible_branching_bkpts_index(bkpts,I)
    if lower_bound:
        if not I_index:
            return min(fn(fractional(I[0])), fn(fractional(I[1])))
        else:
            return min(fn(fractional(I[0])), fn(fractional(I[1])), min(fn(fractional(bkpts[i])) for i in range(I_index[0], I_index[1]+1)))
    else:
        if not I_index:
            return max(fn(fractional(I[0])), fn(fractional(I[1])))
        else:
            return max(fn(fractional(I[0])), fn(fractional(I[1])), max(fn(fractional(bkpts[i])) for i in range(I_index[0], I_index[1]+1)))

def find_branching_bkpt(fn,I_index,extended=False,branching_point_selection='median'):
    r"""
    Select the branching point.
    """
    if extended:
        bkpts=fn.end_points()[:-1]+[1+bkpt for bkpt in fn.end_points()]
    else:
        bkpts=fn.end_points()
    if branching_point_selection=='median':
        return bkpts[floor((I_index[0]+I_index[1])/2)]
    elif branching_point_selection=='sparse':
        distance=0
        for k in range(I_index[0],I_index[1]+1):
        # choose large interval. Could be more flexible?
            new_distance=bkpts[k+1]-bkpts[k-1]
            if new_distance>distance:
                distance=new_distance
                best_bkpt=bkpts[k]
        return best_bkpt
    else:
        raise ValueError("Can't recognize branching_point_selection.")

def find_branching_direction(fn,I,J,K,I_index,J_index,K_index):
    lenI=I[1]-I[0]
    lenJ=J[1]-J[0]
    lenK=K[1]-K[0]
    candidates=[]
    if I_index:
        candidates.append(['I',lenI])
    if J_index:
        candidates.append(['J',lenJ])
    if K_index:
        candidates.append(['K',lenK])
    if not candidates:
        return None
    else:
        ind=np.argmax([candidate[1] for candidate in candidates])
        return candidates[ind][0]

def delta_pi_min(fn,I,J,K,approximation='constant',norm='one',branching_point_selection='median'):
    r"""
    Return the minimum of the function delta fn in the region whose projections are I,J,K.
    """
    bkpts1=fn.end_points()
    bkpts2=fn.end_points()[:-1]+[1+bkpt for bkpt in fn.end_points()]
    vertices=verts(I,J,K)
    if len(vertices)==0:
        return 100
    elif len(vertices)==1:
        return delta_pi(fn,vertices[0][0],vertices[0][1])
    new_I,new_J,new_K=projections(vertices)
    I_index=find_possible_branching_bkpts_index(bkpts1,new_I)
    J_index=find_possible_branching_bkpts_index(bkpts1,new_J)
    K_index=find_possible_branching_bkpts_index(bkpts2,new_K)
    if approximation=='constant':
        alpha_I=fn_constant_bounds(fn,new_I,lower_bound=True,extended=False)
        alpha_J=fn_constant_bounds(fn,new_J,lower_bound=True,extended=False)
        beta_K=fn_constant_bounds(fn,new_K,lower_bound=False,extended=True)
        if alpha_I+alpha_J-beta_K>0:
            return alpha_I+alpha_J-beta_K
    elif approximation=='affine':
        m_I, b_I=find_affine_bounds(fn,new_I,lower_bound=True,extended=False,norm=norm)
        m_J, b_J=find_affine_bounds(fn,new_J,lower_bound=True,extended=False,norm=norm)
        m_K, b_K=find_affine_bounds(fn,new_K,lower_bound=False,extended=True,norm=norm)
        lower_bound=min((m_I*vertex[0]+b_I)+(m_J*vertex[1]+b_J)-(m_K*fractional(vertex[0]+vertex[1])+b_K) for vertex in vertices)
        if lower_bound>0:
            return lower_bound
    else:
        raise ValueError("Can't recognize approximation.")
    upper_bound=min(delta_pi(fn,v[0],v[1]) for v in vertices)
    if upper_bound<0:
        return upper_bound
    branch_direction=find_branching_direction(fn,new_I,new_J,new_K,I_index,J_index,K_index)
    if not branch_direction:
        return upper_bound
    elif branch_direction=='I':
        bkpt=find_branching_bkpt(fn,I_index,extended=False,branching_point_selection=branching_point_selection)
        return min(delta_pi_min(fn,[I[0],bkpt],J,K,approximation=approximation,norm=norm,branching_point_selection=branching_point_selection),delta_pi_min(fn,[bkpt,I[1]],J,K,approximation=approximation,norm=norm,branching_point_selection=branching_point_selection))
    elif branch_direction=='J':
        bkpt=find_branching_bkpt(fn,J_index,extended=False,branching_point_selection=branching_point_selection)
        return min(delta_pi_min(fn,I,[J[0],bkpt],K,approximation=approximation,norm=norm,branching_point_selection=branching_point_selection),delta_pi_min(fn,I,[bkpt,J[1]],K,approximation=approximation,norm=norm,branching_point_selection=branching_point_selection))
    elif branch_direction=='K':
        bkpt=find_branching_bkpt(fn,K_index,extended=True,branching_point_selection=branching_point_selection)
        return min(delta_pi_min(fn,I,J,[K[0],bkpt],approximation=approximation,norm=norm,branching_point_selection=branching_point_selection),delta_pi_min(fn,I,J,[bkpt,K[1]],approximation=approximation,norm=norm,branching_point_selection=branching_point_selection))
    
def lower_triangle_vertices(vertices):
    r"""
    Given the polytope P defined by vertices, return the vertices of the new polytope Q,
    where Q = P \cap {x>=y}, by computing vertices on x=y.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: I=[0,1]
        sage: J=[0,1]
        sage: K=[0,2]
        sage: lower_triangle_vertices(verts(I,J,K))
        [(1, 0), (0, 0), (1, 1)]
    """
    upper_diagonal_vertices=[]
    lower_diagonal_vertices=[]
    on_diagonal_vertices=[]
    for v in vertices:
        if v[0]>v[1]:
            lower_diagonal_vertices.append(v)
        elif v[0]<v[1]:
            upper_diagonal_vertices.append(v)
        else:
            on_diagonal_vertices.append(v[0])
    for v_upper in upper_diagonal_vertices:
        for v_lower in lower_diagonal_vertices:
            new_point_on_diagonal = find_intersection_with_diagonal(v_upper, v_lower)
            on_diagonal_vertices.append(new_point_on_diagonal[0])
    if len(on_diagonal_vertices)==0:
        return lower_diagonal_vertices
    l=min(on_diagonal_vertices)
    u=max(on_diagonal_vertices)
    if l<u:
        return lower_diagonal_vertices+[(l,l),(u,u)]
    else:
        return lower_diagonal_vertices+[(u,u)]

def find_intersection_with_diagonal(v1,v2):
    r"""
    Find a point (x,x) lying on the line segment between v1,v2.
    """
    x1,y1=v1
    x2,y2=v2
    x=(x2*y1-y2*x1)/(y1-y2-x1+x2)
    return (x,x)

def merge_components_with_given_component_strategically(given_component, other_components):
    remaining_components = []
    for component in other_components:
        component_intersection = list(intersection_of_coho_intervals([given_component, component]))
        if component_intersection == given_component:
            # stop early if given_component is already a subset of component.
            return None, other_components
        elif component_intersection:
            # merge if have nontrivial intersection.
            given_component = union_of_coho_intervals_minus_union_of_coho_intervals([given_component, component], [])
        else:
            # if component is not merged into the given_component, it remains in remaining_components.
            remaining_components.append(component)
    return given_component, remaining_components
