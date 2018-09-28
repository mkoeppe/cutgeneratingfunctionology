from itertools import chain

import itertools
import numpy as np

class SubadditivityTestTreeNode :

    def __init__(self, fn, level, intervals):

        self.intervals=intervals
        self.level=level
        self.function=fn
        self.vertices=verts(*intervals)
        self.projections=projections(self.vertices)
        self.left_child=None
        self.right_child=None
        self.parent=None

    def I_bkpts(self):
        if hasattr(self,'_I_bkpts'):
            return self._I_bkpts
        new_I=self.projections[0]
        self._I_bkpts=find_all_bkpts_in_the_interval(self.function.end_points(),new_I)
        return self._I_bkpts

    def J_bkpts(self):
        if hasattr(self,'_J_bkpts'):
            return self._J_bkpts
        new_J=self.projections[1]
        self._J_bkpts=find_all_bkpts_in_the_interval(self.function.end_points(),new_J)
        return self._J_bkpts

    def K_bkpts(self):
        if hasattr(self,'_K_bkpts'):
            return self._K_bkpts
        new_K=self.projections[2]
        self._K_bkpts=find_all_bkpts_in_the_interval(self.function.end_points(),new_K)
        return self._K_bkpts

# affine bounds don't always dominate constant bound

    def delta_pi_constant_lower_bound(self):
        """
        EXAMPLE::

            sage: logging.disable(logging.INFO)
            sage: N=SubadditivityTestTreeNode(gj_forward_3_slope(),0,intervals=[[28/45,4/5],[28/45,4/5],[1+8/45,1+28/45]])
            sage: N.delta_pi_constant_lower_bound()
            5/36
            sage: N.delta_pi_affine_lower_bound(norm='one')
            -0.08333333333333326
            sage: N.delta_pi_affine_lower_bound(norm='inf')
            0.41666666666666674
        """
        alpha_I=min(self.function(bkpt) for bkpt in self.I_bkpts())
        alpha_J=min(self.function(bkpt) for bkpt in self.J_bkpts())
        beta_K=max(self.function(fractional(bkpt)) for bkpt in self.K_bkpts())
        return alpha_I+alpha_J-beta_K

    def delta_pi_affine_lower_bound(self,norm='one'):
        I_values=[self.function(bkpt) for bkpt in self.I_bkpts()]
        J_values=[self.function(bkpt) for bkpt in self.J_bkpts()]
        K_values=[self.function(fractional(bkpt)) for bkpt in self.K_bkpts()]
        self.I_slope, self.I_intercept=find_affine_bounds_new(self.function,self.I_bkpts(),I_values,lower_bound=True,norm=norm)
        self.J_slope, self.J_intercept=find_affine_bounds_new(self.function,self.J_bkpts(),J_values,lower_bound=True,norm=norm)
        self.K_slope, self.K_intercept=find_affine_bounds_new(self.function,self.K_bkpts(),K_values,lower_bound=False,norm=norm)
        delta_lower_bound=min((self.I_slope*vertex[0]+self.I_intercept)+(self.J_slope*vertex[1]+self.J_intercept)-(self.K_slope*(vertex[0]+vertex[1])+self.K_intercept) for vertex in self.vertices)
        return delta_lower_bound

    def delta_pi_fast_affine_lower_bound(self,slope_I,slope_J,slope_K):
        I_values=[self.function(bkpt) for bkpt in self.I_bkpts()]
        J_values=[self.function(bkpt) for bkpt in self.J_bkpts()]
        K_values=[self.function(fractional(bkpt)) for bkpt in self.K_bkpts()]
        intercept_I=find_best_intercept(self.I_bkpts(),I_values,slope_I,lower_bound=True)
        intercept_J=find_best_intercept(self.J_bkpts(),J_values,slope_J,lower_bound=True)
        intercept_K=find_best_intercept(self.K_bkpts(),K_values,slope_K,lower_bound=False)
        delta_lower_bound=min((slope_I*vertex[0]+intercept_I)+(slope_J*vertex[1]+intercept_J)-(slope_K*(vertex[0]+vertex[1])+intercept_K) for vertex in self.vertices)
        return delta_lower_bound

    def delta_pi_lower_bound(self,max_number_of_bkpts=10,norm='inf'):
        """
        Stratigic lower bound of delta pi. If the number of bkpts is small, use affine bound. Use constant bound otherwise.
        """
        I_values=[self.function(bkpt) for bkpt in self.I_bkpts()]
        J_values=[self.function(bkpt) for bkpt in self.J_bkpts()]
        K_values=[self.function(fractional(bkpt)) for bkpt in self.K_bkpts()]
        if len(self.I_bkpts())<=max_number_of_bkpts:
            slope_I, intercept_I=find_affine_bounds_new(self.function,self.I_bkpts(),I_values,lower_bound=True,norm=norm)
        else:
            slope_I=0
            intercept_I=min(I_values)
        if len(self.J_bkpts())<=max_number_of_bkpts:
            slope_J, intercept_J=find_affine_bounds_new(self.function,self.J_bkpts(),J_values,lower_bound=True,norm=norm)
        else:
            slope_J=0
            intercept_J=min(J_values)
        if len(self.K_bkpts())<=max_number_of_bkpts:
            slope_K, intercept_K=find_affine_bounds_new(self.function,self.K_bkpts(),K_values,lower_bound=False,norm=norm)
        else:
            slope_K=0
            intercept_K=max(K_values)
        delta_lower_bound=min((slope_I*vertex[0]+intercept_I)+(slope_J*vertex[1]+intercept_J)-(slope_K*(vertex[0]+vertex[1])+intercept_K) for vertex in self.vertices)
        return delta_lower_bound

    def delta_pi_upper_bound(self):
        return min(delta_pi(self.function,v[0],v[1]) for v in self.vertices)

    def branching_direction(self):
        new_I,new_J,new_K=self.projections
        lenI=new_I[1]-new_I[0]
        lenJ=new_J[1]-new_J[0]
        lenK=new_K[1]-new_K[0]
        candidates=[]
        if len(self.I_bkpts())>2:
            candidates.append(['I',lenI])
        if len(self.J_bkpts())>2:
            candidates.append(['J',lenJ])
        if len(self.K_bkpts())>2:
            candidates.append(['K',lenK])
        if not candidates:
            return None
        else:
            ind=np.argmax([candidate[1] for candidate in candidates])
            return candidates[ind][0]

    def new_intervals(self):
        dir=self.branching_direction()
        if dir=='I':
            new_bkpt=self.I_bkpts()[floor(len(self.I_bkpts())/2)]
            return [[self.intervals[0][0],new_bkpt],self.intervals[1],self.intervals[2]],[[new_bkpt,self.intervals[0][1]],self.intervals[1],self.intervals[2]]
        elif dir=='J':
            new_bkpt=self.J_bkpts()[floor(len(self.J_bkpts())/2)]
            return [self.intervals[0],[self.intervals[1][0],new_bkpt],self.intervals[2]],[self.intervals[0],[new_bkpt,self.intervals[1][1]],self.intervals[2]]
        elif dir=='K':
            new_bkpt=self.K_bkpts()[floor(len(self.K_bkpts())/2)]
            return [self.intervals[0],self.intervals[1],[self.intervals[2][0],new_bkpt]],[self.intervals[0],self.intervals[1],[new_bkpt,self.intervals[2][1]]]
        else:
            raise ValueError, "Indivisible Region."

    def is_divisible(self):
        if len(self.I_bkpts())==2 and len(self.J_bkpts())==2 and len(self.K_bkpts())==2:
            return False
        else:
            return True

    def generate_children(self,upper_bound=0,stop_only_if_strict=True):
        if not self.is_fathomed(upper_bound,stop_only_if_strict):
            I1,I2=self.new_intervals()
            self.left_child=SubadditivityTestTreeNode(self.function,self.level+1,I1)
            self.right_child=SubadditivityTestTreeNode(self.function,self.level+1,I2)
            self.left_child.parent=self
            self.right_child.parent=self

    def is_fathomed(self,upper_bound=0,stop_only_if_strict=True):
        if self.is_divisible():
            if stop_only_if_strict:
                if self.delta_pi_lower_bound()>upper_bound:
                    return True
                else:
                    return False
            else:
                # stop branching early.
                if self.delta_pi_lower_bound()>=upper_bound:
                    return True
                else:
                    return False
        else:
            return True


class SubadditivityTestTree :

    def __init__(self,fn,intervals=[[0,1],[0,1],[0,2]],global_upper_bound=0,objective_limit=0):
        self.function=fn
        self.intervals=intervals
        self.global_upper_bound=global_upper_bound
        self.objective_limit=objective_limit
        self.root=SubadditivityTestTreeNode(fn,0,intervals)
        self.height=0
        self.complete_node_set=set([self.root])
        self.leaf_set=set([self.root])
        # the order of unfathomed nodes matters, like DFS or BFS
        self.unfathomed_node_list=[self.root]
        self.nonsubadditive_vertices=set()
        self.additive_vertices=set()


    def number_of_nodes(self):
        return len(self.complete_node_set)

    def number_of_leaves(self):
        return len(self.leaf_set)

    def node_branching(self,node,find_min=True,stop_only_if_strict=True):
        """
        Branch on a given node, dependent on current self.global_upper_bound.
        """
        if not node.left_child:
            if find_min:
                node.generate_children(self.global_upper_bound,stop_only_if_strict)
            else:
                node.generate_children(self.objective_limit,stop_only_if_strict)
            if node.left_child:
                self.height=max(self.height,node.left_child.level)
                self.complete_node_set.update({node.left_child,node.right_child})
                self.leaf_set.discard(node)
                self.leaf_set.update({node.left_child,node.right_child})
                self.unfathomed_node_list=self.unfathomed_node_list+[node.left_child,node.right_child]

    def next_level(self, node_set,find_min=True,stop_only_if_strict=True):
        """
        Generate nodes in the next level.
        """
        next_level=set()
        for node in node_set:
            self.node_branching(node,find_min,stop_only_if_strict)
            if node.left_child:
                next_level.update({node.left_child,node.right_child})
        return next_level

    def is_subadditive(self,stop_if_fail=False,cache_additive_vertices=True,search_method='BFS'):
        self.unfathomed_node_list=[self.root]
        while self.unfathomed_node_list:
            if search_method=='BFS':
                current_node=self.unfathomed_node_list.pop(0)
            elif search_method=='DFS':
                current_node=self.unfathomed_node_list.pop()
            else:
                raise ValueError, "Can't recognize search_method."
            for v in current_node.vertices:
                delta=delta_pi(self.function,v[0],v[1])
                if delta<self.objective_limit:
                    self.nonsubadditive_vertices.add(v)
                    # can stop early
                    if stop_if_fail:
                        return False
                if delta==self.objective_limit and cache_additive_vertices:
                    self.additive_vertices.add(v)
            self.node_branching(current_node,find_min=False,stop_only_if_strict=cache_additive_vertices)
        return True

    def minimum(self,search_method='BFS'):
        self.unfathomed_node_list=[self.root]
        while self.unfathomed_node_list:
            if search_method=='BFS':
                current_node=self.unfathomed_node_list.pop(0)
            elif search_method=='DFS':
                current_node=self.unfathomed_node_list.pop()
            else:
                raise ValueError, "Can't recognize search_method."
            upper_bound=current_node.delta_pi_upper_bound()
            if upper_bound<self.global_upper_bound:
                self.global_upper_bound=upper_bound
            self.node_branching(current_node,find_min=True,stop_only_if_strict=False)
        return self.global_upper_bound

    def plot_current_regions(self,colorful=False):
        p=Graphics()
        kwds = { 'legend_label1': "indivisible face" , 'legend_label2': "strict subadditive divisible face"}
        legend=[0,0]
        for node in self.leaf_set:
            region=Polyhedron(vertices=node.vertices)
            p+=region.projection().render_outline_2d()
            if colorful:
                if not node.is_divisible():
                    if legend[0]==0:
                        p+=Face(node.intervals).plot(rgbcolor = "yellow",fill_color = "yellow",legend_label=kwds['legend_label1'])
                        legend[0]=1
                    else:
                        p+=Face(node.intervals).plot(rgbcolor = "yellow",fill_color = "yellow")
                elif node.delta_pi_lower_bound()>self.global_upper_bound:
                    if legend[1]==0:
                        p+=Face(node.intervals).plot(rgbcolor = "red",fill_color = "red",legend_label=kwds['legend_label2'])
                        legend[1]=1
                    else:
                        p+=Face(node.intervals).plot(rgbcolor = "red",fill_color = "red")
        return p

def find_all_bkpts_in_the_interval(bkpts,interval):
    """
    Return a list of breakpoints contained in interval and bkpts+Z. bkpts=fn.end_points(), and interval can be any closed interval contained in [0,2].

    EXAMPLES::

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
        i,j=find_bkpts_index_from_zero_to_one(bkpts,interval)
        return [interval[0]]+bkpts[i:j+1]+[interval[1]]
    if interval[0]>=1:
        new_interval=[interval[0]-1,interval[1]-1]
        i,j=find_bkpts_index_from_zero_to_one(bkpts,new_interval)
        return [interval[0]]+[bkpt+1 for bkpt in bkpts[i:j+1]]+[interval[1]]
    interval1=[interval[0],1]
    interval2=[0,interval[1]-1]
    i1,j1=find_bkpts_index_from_zero_to_one(bkpts,interval1)
    i2,j2=find_bkpts_index_from_zero_to_one(bkpts,interval2)
    return [interval[0]]+bkpts[i1:j1+2]+[bkpt+1 for bkpt in bkpts[i2:j2+1]]+[interval[1]]

def find_bkpts_index_from_zero_to_one(bkpts,interval):
    i=bisect_left(bkpts, interval[0])
    j=bisect_left(bkpts, interval[1])-1
    if interval[0]==bkpts[i]:
        i=i+1
    return i,j

def plot_2d_regions(fn,colorful=False,find_min=True,stop_only_if_strict=True):
    T=SubadditivityTestTree(fn)
    p=Graphics()
    # current_level=T.complete_node_set results in error: Set changed size during iteration
    current_level=set([T.root])
    while current_level:
        p+=plot_2d_regions_in_one_level(current_level,colorful=colorful)
        p.show()
        next_level=T.next_level(current_level,find_min,stop_only_if_strict)
        current_level=next_level

def plot_2d_regions_in_one_level(node_set,colorful=False):
    p=Graphics()
    for node in node_set:
        region=Polyhedron(vertices=node.vertices)
        p+=region.projection().render_outline_2d()
        if colorful:
            if not node.is_divisible():
                p+=Face(node.intervals).plot(rgbcolor = "yellow", fill_color = "yellow")
            elif node.delta_pi_lower_bound()>0:
                p+=Face(node.intervals).plot(rgbcolor = "red", fill_color = "red")
    return p

def values_of_delta_pi_over_grid(fn,q):
    """
    Return a matrix representing the values of delta pi over the grid (1/q)Z*(1/q)Z.
    """
    res=np.zeros((q+1,q+1))
    for i in range(q+1):
        for j in range(q+1):
            res[q-j][i]=delta_pi(fn,i/q,j/q)
    return res

def strategic_delta_pi_min(fn,approximation='constant',norm='one',branching_point_selection='median'):
    """
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
    """
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

                           
def find_best_slope_intercept(X,Y,lower_bound=True,solver='GLPK',norm='one'):
    """
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
        raise ValueError, "Can't recognize norm."
    p.solve()
    return p.get_values(m),p.get_values(b)

def find_affine_bounds_new(fn,X,Y,lower_bound=True,norm='one'):
    if len(X)==2:
        slope,intercept=find_slope_intercept_trivial(X,Y)
    else:
        slope,intercept=find_best_slope_intercept(X,Y,lower_bound=lower_bound,norm=norm)
    return slope, intercept


def find_affine_bounds(fn,I,lower_bound=True,extended=False,norm='one'):
    """
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
    """
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
    """
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
    """
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
        raise ValueError, "Can't recognize branching_point_selection."

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
    """
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
        delta_lower_bound=min((m_I*vertex[0]+b_I)+(m_J*vertex[1]+b_J)-(m_K*fractional(vertex[0]+vertex[1])+b_K) for vertex in vertices)
        if delta_lower_bound>0:
            return delta_lower_bound
    else:
        raise ValueError, "Can't recognize approximation."
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
    
