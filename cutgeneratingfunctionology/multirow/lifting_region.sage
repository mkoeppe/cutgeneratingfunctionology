import itertools
from sage.geometry.integral_points import rectangular_box_points
from six.moves import range
    
def random_point_in_polytope(B):
    vertices= B.vertices_list()
    n = len(vertices)
    random_ints = [0] + sorted([randint(0,100) for i in range(n - 1)]) + [100]
    coeffs = [(random_ints[i+1]-random_ints[i])/100 for i in range(n)]
    pt = sum(vector(vertices[i])*coeffs[i] for i in range(n))
    return tuple(pt)

def lifting_regions(B, f, shifted=False):
    # assume f is in B.
    # assume B is lattice-free
    d = B.dim()
    regions = []
    for facet in B.faces(d - 1):
        p = facet.as_polyhedron()
        if f in p:
            continue
        ray_vectors = [vector(f)-vector(v) for v in p.vertices()]
        cone_f = Polyhedron(rays=[tuple(-ray) for ray in ray_vectors], vertices=[tuple(f)])
        for v in p.integral_points():
            cone_v = Polyhedron(rays=[tuple(ray) for ray in ray_vectors], vertices=[tuple(v)])
            region_v = cone_f.intersection(cone_v)
            if not region_v.is_full_dimensional():
                continue
            if shifted:
                regions.append(region_v - v) # shift to origin
            else:
                regions.append(region_v)
    return regions

def lifting_regions_tile_box(regions, box=None):
    # box is a tuple describing the lowerleft and upperright corners, or a polyhedron type, or default = unit cube.
    # assume box is integral
    d = regions[0].dim() #assume that regions are full-dimensional
    if hasattr(box, 'bounding_box'):
        box_min, box_max = box.bounding_box()
    else:
        if box is None:
            box_min, box_max = [0] * d, [1] * d
        else:
            box_min, box_max = box[0], box[1]
        box_vertices = [[box_min[i] if v[i]==0 else box_max[i] for i in range(d)] for v in itertools.product([0, 1], repeat=d)]
        box = Polyhedron(vertices = box_vertices)
    regions_tile = []
    for r in regions:
        r_box = r.bounding_box(integral=True)
        shift_min = vector(box_min) - vector(r_box[1]) + vector([1]*d)
        shift_max = vector(box_max) - vector(r_box[0]) - vector([1]*d)
        shift_vectors = rectangular_box_points(list(shift_min), list(shift_max), None)
        for v in shift_vectors:
            r_v = r + v
            r_v_in_box = r_v.intersection(box)
            if r_v_in_box.is_full_dimensional():
                regions_tile.append(r_v_in_box)
    return regions_tile

def volume_of_union_of_polytopes(polyhedra):
    cache_vol = dict()
    polys = []
    signs = []
    vol = 0
    for p in polyhedra:
        polys.append(p)
        signs.append(1)
        p_shift = tuple(sorted((tuple(vector(v) - p.center()) for v in p.vertices()))) # caching vertices, avoid setting up the shifted polytope.
        if p_shift in cache_vol:
            p_vol = cache_vol[p_shift]
        else:
            p_vol = p.volume()
            cache_vol[p_shift] = p_vol
        vol += p_vol #vol += p.volume()
        n = len(polys)
        for i in range(n-1):
            q = p.intersection(polys[i])
            if q.is_full_dimensional():
                sign = -signs[i]
                q_shift = tuple(sorted((tuple(vector(v) - q.center()) for v in q.vertices())))
                if q_shift in cache_vol:
                    q_vol = cache_vol[q_shift]
                else:
                    q_vol = q.volume()
                    cache_vol[q_shift] = q_vol
                vol += sign * q_vol #q.volume()
                polys.append(q)
                signs.append(sign)
    return vol
        
def volume_of_lifting_region(B, f=None, show_plots=False, return_plots=False):
    r"""
    Examples::

        sage: from cutgeneratingfunctionology.multirow import *
        sage: M11 = Polyhedron(vertices=[[-1, 0, 0], [0, 0, 2], [0, 2, 0], [1, 0, 0], [1, 2, 2], [2, 0, 2]])
        sage: volume_of_lifting_region(M11)
        11/12
        sage: B = Polyhedron(vertices=[[-3/13, 21/13], [1 - 4/10, 3], [3/2, 3/4]])
        sage: f = B.center()
        sage: volume_of_lifting_region(B, f)
        1763/2600
        sage: B = Polyhedron(vertices=[[-3/5, 6/5], [-3, 0], [3, 0]])
        sage: f = (-39/25,27/50)
        sage: volume_of_lifting_region(B, f) # long time # 12 s
        1
    """
    if f is None:
        f = random_point_in_polytope(B)
    regions = lifting_regions(B, f)
    regions_tile = lifting_regions_tile_box(regions)
    vol = volume_of_union_of_polytopes(regions_tile)
    result = vol
    if show_plots or return_plots:
        d = B.dim()
        box = B.bounding_box(integral=True)
        lattices = [tuple(pt) for pt in rectangular_box_points(list(box[0]), list(box[1]), None)]
        if d < 3:
            gb = B.plot(alpha=0.4, color='yellow', title="$f=%s$" % (latex(f))) + point(f, color='red') + text("$f$", f, vertical_alignment='bottom',horizontal_alignment='left', color='black')   
            gb += sum(p.plot(color='red',alpha=0.4) for p in regions)
            gb += points(lattices, color='gray')
            gtile = sum(p.plot(color='red',alpha=0.3, title = r"$f=%s,\quad \mathrm{vol}_{\mathbb{T}^{%s}} \bar{R}(f)/\mathbb{Z}^{%s} = %s$" % ((latex(f)), latex(d), latex(d), latex(vol))) for p in regions_tile)
            #gtile += points(box, color='white')
            if return_plots:
                result = vol, gb, gtile
            if show_plots:
                g = graphics_array([gb, gtile])
                g.show()
        elif d == 3:
            gb = B.plot(point=False, line='blue', polygon=False)
            #gb += sum(p.plot(alpha=0.4, color='red') for p in regions)
            # workaround. polygon plot forgets about keyword argument opacity. 
            for p in regions:
                for facet in p.faces(2):
                    vertices = [list(v) for v in facet.vertices()]
                    face = polygon3d(vertices, color='red')
                    face._extra_kwds['opacity'] = 0.4
                    gb += face
            gb += points(lattices, color='gray')
            gtile = sum(p.plot(color='red',alpha=0.3) for p in regions_tile)
            if return_plots:
                result = vol, gb, gtile
            if show_plots:
                gb.show(viewer='threejs')
                gtile.show(viewer='threejs')
    return result

def is_lattice_free_polyhedron(B):
    for pt in B.integral_points():
        if B.relative_interior_contains(pt):
            return False
    return True

def is_maximal_lattice_free_polyhedron(B):
    return is_lattice_free_polyhedron(B) and all(is_lattice_free_polyhedron(p.as_polyhedron()) == False for p in B.faces(B.dim()-1))

def is_volume_affine_in_f(B, num_tests=10):
    # assume B is full-dimensional
    d = B.dim()
    vertices = B.vertices_list()
    A = Matrix([v+[1] for v in vertices])
    Y = vector([volume_of_lifting_region(B, v, False) for v in vertices])
    X = A.solve_right(Y)
    print(X)
    for i in range(num_tests):
        f = random_point_in_polytope(B)
        vol = volume_of_lifting_region(B, f, False)
        print(i, f, vol)
        if not vector(f+(1,)) * X == vol:
            return False
    return True

def plot_polyhedron_and_f(B,f):
    d = B.dim()
    box = B.bounding_box(integral=True)
    lattices = [tuple(pt) for pt in rectangular_box_points(list(box[0]), list(box[1]), None)]
    if d < 3:
        gb = B.plot(alpha=0.5, color='yellow') + point(f, color='red', size=30,zorder=1)
        gb += points(lattices, color='blue', size=30, zorder=1)
        #gb.show(figsize=4)
    elif d == 3:
        gb = point3d(f, color='red', size=15)
        gb += B.plot(point=False, line='orange', polygon=False)
        #gb += sum(p.plot(alpha=0.4, color='red') for p in regions)
        # workaround. polygon plot forgets about keyword argument opacity. 
        for facet in B.faces(2):
            vertices = [list(v) for v in facet.vertices()]
            face = polygon3d(vertices, color='yellow')
            face._extra_kwds['opacity'] = 0.5
            gb += face
        gb += points(lattices, color='blue', size=15)
        #gb.show(viewer='threejs')
    return gb
