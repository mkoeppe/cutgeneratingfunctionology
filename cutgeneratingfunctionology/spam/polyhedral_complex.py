r"""
PolyhedralComplex
"""

from __future__ import division, print_function, absolute_import

from copy import copy
try:
    from sage.topology.cell_complex import GenericCellComplex
except ImportError:
    from sage.homology.cell_complex import GenericCellComplex
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.modules.free_module_element import vector

class PolyhedralComplex(GenericCellComplex):
    r"""
    Define a PolyhedralComplex.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: pc = PolyhedralComplex([Polyhedron(base_ring=QQ, vertices=[(1/3, 1/3), (QQ(0), QQ(0)), (1/7, 2/7)]), Polyhedron(base_ring=QQ, vertices=[(1/7, 2/7), (QQ(0), QQ(0)), (QQ(0), 1/4)])])
        sage: list(sorted( sorted(pol.vertices_list()) for pol in pc.cells_list() ))
        [[[0, 0]],
         [[0, 0], [0, 1/4]],
         [[0, 0], [0, 1/4], [1/7, 2/7]],
         [[0, 0], [1/7, 2/7]],
         [[0, 0], [1/7, 2/7], [1/3, 1/3]],
         [[0, 0], [1/3, 1/3]],
         [[0, 1/4]],
         [[0, 1/4], [1/7, 2/7]],
         [[1/7, 2/7]],
         [[1/7, 2/7], [1/3, 1/3]],
         [[1/3, 1/3]]]
        sage: pc.is_convex()
        True
        sage: p = pc.union_as_polyhedron()
        sage: p.Hrepresentation()
        (An inequality (1, -4) x + 1 >= 0,
         An inequality (-1, 1) x + 0 >= 0,
         An inequality (1, 0) x + 0 >= 0)

    TESTS::

        sage: pc.face_poset()     # known bug

    """
    def __init__(self, maximal_cells, maximality_check=True, **kwds):
        r"""
        INPUT:

        - ``maximal_cells`` -- a list, a tuple, or a dictionary (indexed by dimension) of cells of the Complex. Each cell is of class :class:`Polyhedron` of the same ambient dimension. To set up a :class:PolyhedralComplex, it is sufficient to provide the maximal faces. Use keyword argument partial=``True`` to set up a partial polyhedral complex, which is a subset of the faces (viewed as relatively open) of a polyhedral complex that is not necessarily closed under taking faces.
        - ``maximality_check`` -- boolean; default ``False``; if is ``True``, check that each given maximal cells are indeed maximal. In this case, when producing the internal representation of the polyheral complex, omit those that are not (move them from ``self._maximal_cells`` to ``self._non_maximal_cells_given``).

        - ``partial`` -- boolean; default ``False``;

        - ``check_face_to_face`` -- boolean; default ``False``;
        """
        if isinstance(maximal_cells, (list, tuple)):
            cells_dict = {}
            for cell in maximal_cells:
                d = cell.dimension()
                if d in cells_dict:
                    cells_dict[d].append(cell)
                else:
                    cells_dict[d] = [ cell ]
        elif isinstance(maximal_cells, dict):
            cells_dict = copy(maximal_cells)
        else:
            raise ValueError
        if not cells_dict:
            return
        self._dim = list(cells_dict.keys())[-1]
        self._ambient_dim = list(cells_dict.values())[0][0].ambient_dim()

        #TODO: partial
        self._partial = kwds.get('partial', False)
        #TODO: check_face_to_face
        check_face_to_face = kwds.pop('check_face_to_face', False)

        if not maximality_check:
            self._maximal_cells = cells_dict
            self._non_maximal_cells_given = {}
        else:
            self._maximal_cells = {}
            self._non_maximal_cells_given = {}
            for k in range(self._dim, -1, -1):
                k_faces = []
                k_faces_not_maximal = []
                for face in cells_dict.get(k, []):
                    maximal_face = True
                    for (l, l_faces) in self._maximal_cells.items():
                        for c in l_faces:
                            for f in c.faces(k):
                                if face == f.as_polyhedron():
                                    maximal_face = False
                                    break
                            if not maximal_face:
                                break
                        if not maximal_face:
                            break
                    if maximal_face:
                        k_faces.append(face)
                    else:
                        k_faces_not_maximal.append(face)
                if k_faces:
                    self._maximal_cells[k] = k_faces
                if k_faces_not_maximal:
                    self._non_maximal_cells_given[k] = k_faces_not_maximal

    def cells(self, subcomplex=None):
        if hasattr(self, '_cells'):
            return self._cells
        if subcomplex is not None:
            raise NotImplementedError
        cells = copy(self._maximal_cells)
        for k in range(self._dim, 0, -1):
            k_cells =  cells.get(k,[])
            k_minus_1_cells = set(cells.get(k-1,[]))
            for c in k_cells:
                k_minus_1_cells.update([f.as_polyhedron() for f in c.faces(k-1)])
            if k_minus_1_cells:
                cells[k-1] = list(k_minus_1_cells)
        self._cells = cells
        return cells

    def cells_list(self):
        cells_dict = self.cells()
        cells = []
        for kcells in cells_dict.values():
            cells += kcells
        return cells

    def maximal_cells(self, k=None):
        if k is None:
            return copy(self._maximal_cells)
        return self._maximal_cells.get(k, [])

    def maximal_cells_list(self):
        cells_dict = self.maximal_cells()
        cells = []
        for kcells in cells_dict.values():
            cells += kcells
        return cells

    def dimension(self):
        r"""
        The dimension of this cell complex: the maximum dimension of its cells.
        """
        return self._dim

    def ambient_dimension(self):
        return self._ambient_dim

    def is_pure(self):
        return len(self._maximal_cells) == 1

    def is_full_dimensional(self):
        return self._dim == self._ambient_dim

    def stratify(self, k):
        k_faces = self.maximal_cells(k)
        if k_faces:
            return PolyhedralComplex(k_faces, maximality_check=False)

    def boundary_cells(self):
        r"""
        Returns a list of co-dimension one cells on the boundary.
        """
        if not self.is_pure():
            raise NotImplementedError
        k = self._dim-1
        k_faces = [f.as_polyhedron() for face in self._maximal_cells[k+1] for f in face.faces(k)]
        return [face for face in k_faces if k_faces.count(face)==1]

    def is_convex(self):
        """
        sage: from cutgeneratingfunctionology.spam.polyhedral_complex import PolyhedralComplex
        sage: from sage.geometry.polyhedron.constructor import Polyhedron
        sage: pc = PolyhedralComplex([Polyhedron(base_ring=QQ, vertices=[[1,0],[0,1]],
        ....: rays=[[1,0],[0,1]])])
        sage: pc.is_convex()
        True
        sage: pc = PolyhedralComplex([Polyhedron(base_ring=QQ, vertices=[[1,0,0],[0,1,0]],
        ....: rays=[[1,0,0],[0,1,0]])])
        sage: pc.is_convex()
        True
        sage: pc = PolyhedralComplex([Polyhedron(base_ring=QQ, vertices=[[-1,0],[1,0]],
        ....: lines=[[0,1]])])
        sage: pc.is_convex()
        True
        sage: pc = PolyhedralComplex([Polyhedron(base_ring=QQ, vertices=[[1,0],[0,1]],
        ....: rays=[[1,0],[0,1]]),Polyhedron(base_ring=QQ, vertices=[[1,0],[0,-1]],
        ....: rays=[[1,0],[0,-1]])])
        sage: pc.is_convex()
        False
        sage: pc = PolyhedralComplex([Polyhedron(base_ring=QQ, vertices=[[0,0]],
        ....: rays=[[1,0],[-1,1]]),Polyhedron(base_ring=QQ, vertices=[[0,0]],
        ....: rays=[[1,0],[-1,-1]])])
        sage: pc.is_convex()
        False
        sage: pc = PolyhedralComplex([Polyhedron(base_ring=QQ, vertices=[[0,0,0]],
        ....: rays=[[1,0,0],[-1,1,0]]),Polyhedron(base_ring=QQ, vertices=[[0,0,0]],
        ....: rays=[[1,0,0],[-1,-1,0]])])
        sage: pc.is_convex()
        False
        sage: pc = PolyhedralComplex([Polyhedron(base_ring=QQ, vertices=[[0,0,0]],rays=[[1,0,0],[0,1,0],[0,0,-1]]),
        ....: Polyhedron(base_ring=QQ, vertices=[[0,0,0]],rays=[[1,0,0],[0,-1,0],[0,0,-1]]),
        ....: Polyhedron(base_ring=QQ, vertices=[[0,0,0]],rays=[[1,0,0],[0,-1,0],[0,0,1]]),
        ....: Polyhedron(base_ring=QQ, vertices=[[0,0,0]],rays=[[-1,0,0],[0,-1,0],[0,0,-1]]),
        ....: Polyhedron(base_ring=QQ, vertices=[[0,0,0]],rays=[[-1,0,0],[0,-1,0],[0,0,1]]),
        ....: Polyhedron(base_ring=QQ, vertices=[[0,0,0]],rays=[[-1,0,0],[0,1,0],[0,0,-1]]),
        ....: Polyhedron(base_ring=QQ, vertices=[[0,0,0]],rays=[[-1,0,0],[0,1,0],[0,0,1]])])
        sage: pc.is_convex()
        False
        """
        if hasattr(self, '_is_convex'): # FIXME: bad! _is_convex can not be changed later.
            return self._is_convex
        if not self.is_pure():
            self._is_convex = False
            return False
        if not self.is_full_dimensional():
            from sage.modules.free_module import span
            face = self._maximal_cells[self._dim][0]
            affine_space = span(face.equations_list(), face.base_ring())
            for face in self._maximal_cells[self._dim][1::]:
                if span(face.equations_list(), face.base_ring()) != affine_space:
                    self._is_convex = False
                    return False
            # # If they lie in different subspaces, can't be convex. When they all lie in the same subspace, then you orient the boundary halfspaces toward a strict convex combination of the vertices. Then you check whether all vertices are contained. After you made sure that the affine hulls of the cells are the same, it does not matter that is not full dimensional.
        boundaries = self.boundary_cells()
        vertices = set([])
        rays = set([])
        lines = set([])
        # lines are useless, because they are in the affine space of each boundary cell.
        for cell in boundaries: # is that enough, or need vertices of all cells? I think that is enough.
            for v in cell.vertices_list():
                vv = vector(v)
                vv.set_immutable()
                vertices.add(vv)
        for cell in self._maximal_cells[self._dim]:
            for r in cell.rays_list():
                rr = vector(r)
                rr.set_immutable()
                rays.add(rr)
            for l in cell.lines_list():
                ll = vector(l)
                ll.set_immutable()
                lines.add(ll)
        center = sum(vertices) / len(vertices)
        for cell in boundaries:
            for equation in cell.equations_list(): # if not full-dim, cell has more than one equaiton.
                coeff = vector(equation[1::])
                const = equation[0]
                if const + coeff * center == 0:
                    sign = 0
                elif const + coeff * center > 0:
                    sign = 1
                    for v in vertices:
                        if const + coeff * v < 0:
                            self._is_convex = False
                            return False
                elif const + coeff * center < 0:
                    sign = -1
                    for v in vertices:
                        if const + coeff * v > 0:
                            self._is_convex = False
                            return False
                for r in rays:
                    if sign == 0:
                        sign = coeff * r
                    else:
                        if sign * (coeff * r) < 0:
                            self._is_convex = False
                            return False
        self._is_convex = True
        self._polyhedron = Polyhedron(vertices=vertices,rays=rays,lines=lines)
        return True

    def union_as_polyhedron(self):
        if not self.is_convex():
            raise ValueError("The polyhedral complex is not convex.")
        return self._polyhedron

    def has_maximal_cell(self, c):
        d = c.dimension()
        if not d in self._maximal_cells:
            return False
        if c in self._maximal_cells[d]:
            return True
        return False

    def has_cell(self, c):
        d = c.dimension()
        cells = self.cells()
        if not d in cells:
            return False
        if c in cells[d]:
            return True
        return False

    def has_non_maximal_cell_given(self, c):
        d = c.dimension()
        if not d in self._non_maximal_cells_given:
            return False
        if c in self._non_maximal_cells_given[d]:
            return True
        return False

    def plot(self, **kwds):
        return sum(cell.plot()
                   for stratum in self.maximal_cells().values()
                   for cell in stratum)
