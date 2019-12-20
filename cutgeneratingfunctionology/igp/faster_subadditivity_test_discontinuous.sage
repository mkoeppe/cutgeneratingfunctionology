from itertools import chain

import queue as queue
import itertools
import numpy as np

class SubadditivityTestTreeNodeGeneral(object):

    def __init__(self, fn, level, I, J, K):

        """
        Interval I = [(x1,epsilon1),(x2,epsilon2)]
        epsilon = 0 means at x, epsilon = 1 means right limit, epsilon = -1 means left limit.
        """

        self.intervals=tuple(tuple(I[0][0],I[1][0]),tuple(J[0][0],J[1][0]),tuple(K[0][0],K[1][0]))
        self.I=I
        self.J=J
        self.K=K
        self.level=level
        self.function=fn
        self.vertices=verts(*intervals)
        self.projections=projections(self.vertices)
        self.left_child=None
        self.right_child=None
        self.parent=None
        self.affine_estimators=None

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
        self._I_bkpts = tuple(find_all_bkpts_in_the_interval_general(self.function.end_points(),new_I))
        return self._I_bkpts

    def J_bkpts(self):
        if hasattr(self,'_J_bkpts'):
            return self._J_bkpts
        new_J=self.projections[1]
        self._J_bkpts = tuple(find_all_bkpts_in_the_interval_general(self.function.end_points(),new_J))
        return self._J_bkpts

    def K_bkpts(self):
        if hasattr(self,'_K_bkpts'):
            return self._K_bkpts
        new_K=self.projections[2]
        self._K_bkpts = tuple(find_all_bkpts_in_the_interval_general(self.function.end_points(),new_K))
        return self._K_bkpts

    def I_values(self):
        if hasattr(self,'_I_values'):
            return self._I_values
        self._I_values = tuple(self.function.limit(bkpt[0],bkpt[1]) for bkpt in self.I_bkpts())
        return self._I_values

    def J_values(self):
        if hasattr(self,'_J_values'):
            return self._J_values
        self._J_values = tuple(self.function.limit(bkpt[0],bkpt[1]) for bkpt in self.J_bkpts())
        return self._J_values

    def K_values(self):
        if hasattr(self,'_K_values'):
            return self._K_values
        fractional_K_bkpts=[tuple(fractional(b[0]),b[1]) for b in self.K_bkpts()]
        self._K_values = tuple(self.function.limit(bkpt[0],bkpt[1]) for bkpt in fractional_K_bkpts)
        return self._K_values
