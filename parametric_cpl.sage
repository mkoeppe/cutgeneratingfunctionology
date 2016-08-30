# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

def cpl_n_group_function(n, cpleq=False, merge=True):
    return CPLFunctionsFactory(n, cpleq, merge)

class CPLFunctionsFactory:
    """
    A Factory of CPL functions.

    EXAMPLES::

        sage: logging.disable(logging.INFO)
        sage: cpl3 = cpl_n_group_function(3)
        sage: h = cpl3(f=1/6, z=(1/12, 1/12), o=(1/5, 0))
        sage: extremality_test(h)
        True

        Omit the argument o, then it takes the mapping self._theta.
        With the default value of self._theta, it creates a gmic function.
        sage: g = cpl3(1/6, (1/12, 1/12))
        sage: g == gmic(1/6)

        Change self._theta:
        sage: cpl3._theta = lambda f, z: (1/5, 0)
        sage: hh = cpl3(1/6, (1/12, 1/12))
        sage: h == hh
        True
    """
    def __init__(self, n, cpleq=False, merge=True):
        self._n = n
        self._theta = lambda f, z: tuple([z[i]/(1-f) for i in range(self._n-1)])
        self._cpleq = cpleq # in cpl3eq, z1=z2
        self._merge = merge # need to set merge=False in cpl_regions_from_arrangement_of_bkpts
    def __call__(self, f, z, o=None, field=None):
        if o is None:
            o = self._theta(f, z)
        if not (bool(0 < f < 1) & bool(all(0 < zi for zi in z)) & bool (sum(z) <= (1-f)/2)):
            raise ValueError, "Bad parameters. Unable to construct the function."
        if not (bool(0 <= oi for oi in o) & bool(sum(o) <= 1/2)):
            raise ValueError, "Bad thetas parameters. function value outside [0,1]."
        if sum(z) < (1-f)/2:
            m = self._n
        else:
            m = self._n - 1 
        bkpts = [0] + [f + sum(z[0:i]) for i in range(m)] + \
                [1 - sum(z[0:i]) for i in range(m - 1, -1, -1)]
        phi_values = [0] + [sum(o[0:i]) for i in range(m)] + \
                     [1 - sum(o[0:i]) for i in range(m - 1, -1, -1)]
        values = [(bkpts[i] - phi_values[i])/f for i in range(len(bkpts))]
        return piecewise_function_from_breakpoints_and_values(bkpts, values, field=field, merge=self._merge)

def cpl_regions_from_arrangement_of_bkpts(n=3, cpleq=True, max_iter=0, flip_ineq_step=1/3000, check_completion=False, wall_crossing_method='heuristic', goto_lower_dim=True):
    """
    Got regions[0:143]: 3-dim; regions[143:420]: 2-dim; regions[420:581]: 1-dim; regions[581:607]: 0-dim.

    sage: logging.disable(logging.INFO)                      # not tested
    sage: regions = cpl_regions_from_arrangement_of_bkpts(3) # not tested
    sage: len(regions)                                       #not tested
    607
    """
    cpln = cpl_n_group_function(n, cpleq, merge=False)
    if cpleq:
        var_name = ['f', 'z']
        var_value = [1/5, 1/9/n]
    else:
        var_name = ['f']+['z%s' % i for i in range(1, n)]
        var_value = [1/5]+[1/(5 * 2^i) for i in range(1, n)]
    def frt(K, h):
        return find_region_type_igp(K, h, region_level='minimal', is_minimal=None)
    arr_complex=SemialgebraicComplex(cpln, var_name, max_iter=max_iter, find_region_type=frt)
    if flip_ineq_step != 0:
        arr_complex.bfs_completion(var_value, \
                                   flip_ineq_step=flip_ineq_step, \
                                   check_completion=check_completion, \
                                   wall_crossing_method=wall_crossing_method, \
                                   goto_lower_dim=goto_lower_dim)
    else:
        arr_complex.shoot_random_points(1000)
    regions = [c for c in arr_complex.components if c.region_type != 'not_constructible']
    regions.sort(key=lambda r: len(r.leq))
    return regions

def plot_cpl_components(components, show_testpoints=False):
    """
    sage: regions = cpl_regions_from_arrangement_of_bkpts(3)
    sage: g = plot_cpl_components(regions)
    sage: g.show(xmin=0, xmax=1, ymin=0, ymax=1/4) #not tested
    """
    g = Graphics()
    for c in components:
        if not c.leq:
            g += c.plot(show_testpoints=show_testpoints)
    for c in components:
        if len(c.leq)==1:
            g += c.plot(show_testpoints=show_testpoints)
    for c in components:
        if len(c.leq)==2:
            ptcolor = find_region_color(c.region_type)
            g += point(c.var_value, color = ptcolor, zorder=10)
    return g
