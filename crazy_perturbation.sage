# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

# Reminder: need coerce all input to common RNF.

class CrazyPiece:
    def __init__(self, interval, generators=None, cosets=[]):
        # assume interval is open, represented as (a, b)
        self.interval = interval
        # assume that generators QQ linearly independent
        self.generators = generators
        if not is_QQ_linearly_independent(*generators):
            logging.warn("Generators are not linearly independent over Q.")
        if len(generators) < 2:
            logging.warn("The group is not dense.")
        self.hermite_form_generators = find_hermite_form_generators(generators)
        # cosets is a list of (coset_repr, shift),
        # assume coset_repr represent distinct cosets;
        # assume shifts are non-zero
        self.cosets = cosets
        for i in range(len(cosets)):
            (r, s) = cosets[i]
            if s == 0:
                logging.warn("Have shift = 0.")
            for j in range(i):
                rr = cosets[j][0]
                if is_in_group_over_ZZ(r-rr, generators):
                    logging.warn("Not unique coset representative.")

    def __call__(self, x):
        ###assert self.interval[0] < x < self.interval[1]
        for (r,s) in self.cosets:
            if is_in_group_over_ZZ(x-r, self.generators):
                return s
        return x.parent().zero()

class CrazyPerturbation:
    # assume that all inputs are elements from a same RNF.
    def __init__(self, pwl, crazy_pieces):
        self.pwl = pwl
        # assume that crazy pieces' intervals are disjoint.
        self.crazy_pieces = crazy_pieces

    def find_crazy_piece(self, x, xeps=0):
        x = fractional(x)
        if x == 1 and xeps == 1:
            x = x.parent().zero()
        elif x == 0 and xeps == -1:
            x = x.parent().one()
        for cp in self.crazy_pieces:
            if (cp.interval[0] < x < cp.interval[1]) or \
               (xeps == 1) and (x == cp.interval[0]) or \
               (xeps == -1) and (x == cp.interval[1]):
                return cp
        return None

    @cached_method
    def __call__(self, x):
        crazy_piece = self.find_crazy_piece(x, 0)
        if crazy_piece is None:
            return self.pwl(x)
        else:
            return self.pwl(x) + crazy_piece(x)

    def plot(self):
        pwl = self.pwl
        g = pwl.plot(color='magenta')
        for crazy_piece in self.crazy_pieces:
            a, b = crazy_piece.interval[0], crazy_piece.interval[1]
            pwla = pwl.limit(a, 1)
            pwlb = pwl.limit(b, -1)
            max_shift = max([s for (r, s) in crazy_piece.cosets])
            min_shift = min([s for (r, s) in crazy_piece.cosets])
            g += polygon([(a, pwla+min_shift), (a, pwla+max_shift), (b, pwlb+max_shift), (b, pwlb+min_shift)], color='magenta')
        return g

    def range(self):
        pwl = self.pwl
        lb = min(flatten(pwl.limits_at_end_points()))
        ub = max(flatten(pwl.limits_at_end_points()))
        for crazy_piece in self.crazy_pieces:
            pwla = pwl.limit(crazy_piece.interval[0],1)
            pwlb = pwl.limit(crazy_piece.interval[1],-1)
            max_shift = max([s for (r, s) in crazy_piece.cosets])
            min_shift = min([s for (r, s) in crazy_piece.cosets])
            if min_shift + min(pwla, pwlb) < lb:
                lb = min_shift + min(pwla, pwlb)
            if max_shift + max(pwla, pwlb) > ub:
                ub = max_shift + max(pwla, pwlb)
        return (lb, ub)

    def end_points(self):
        bkpts = copy(self.pwl.end_points())
        for crazy_piece in self.crazy_pieces:
            bkpts += [crazy_piece.interval[0], crazy_piece.interval[1]]
        return uniq(bkpts)

def is_in_group_over_ZZ(x, generators):
    # assume that all inputs are elements from a same RNF.
    # generators are linearly independent over Q
    lgens = [g.list() for g in generators]
    lx = x.list()
    #if rank(matrix(QQ,lgens+[lx])) != len(generators):
    #    return False
    try:
        s = matrix(QQ,lgens).solve_left(matrix(QQ,lx))
        return all((si in ZZ) for si in s[0])
    except ValueError:
        return False

def find_hermite_form_generators(generators):
    lgens = [g.list() for g in generators]
    return (matrix(QQ,lgens).hermite_form())

def find_epsilon_for_crazy_perturbation(fn, cp, show_plots=False):
    """
    EXAMPLE::

        sage: logging.disable(logging.INFO)
        sage: fn = discontinuous_bhk_irrational_dense_move_not_affine()
        sage: bkpt = fn.end_points()
        sage: generators = [bkpt[4]-bkpt[2], bkpt[6]-bkpt[4]]
        sage: pwl = piecewise_function_from_breakpoints_and_slopes([0,1],[0])
        sage: crazy_piece_1 = CrazyPiece((bkpt[8], bkpt[9]), generators, [(bkpt[8], 1), (bkpt[9], -1)])
        sage: crazy_piece_2 = CrazyPiece((bkpt[10], bkpt[11]), generators, [(bkpt[10], 1), (bkpt[11], -1)])
        sage: cp = CrazyPerturbation(pwl, [crazy_piece_1, crazy_piece_2])
        sage: find_epsilon_for_crazy_perturbation(fn, cp)
        RNF0.0001237724214802864?

        sage: h = bhk_discontinuous_irrational()
        sage: find_epsilon_for_crazy_perturbation(h, cp)
        0
    """
    # assume fn is a subadditive pwl function, cp is a non_zero CrazyPerturbation
    bkpt = uniq(copy(fn.end_points())+cp.end_points())
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]
    type_1_vertices = [(x, y, x+y) for x in bkpt for y in bkpt if x <= y]
    type_2_vertices = [(x, z-x, z) for x in bkpt for z in bkpt2 if x < z < 1+x]
    vertices = set(type_1_vertices + type_2_vertices)
    g = plot_2d_complex(fn)
    epsilon = 3
    for (x, y, z) in vertices:
        for (xeps, yeps, zeps) in [(0,0,0)]+list(nonzero_eps):
            deltafn = delta_pi_general(fn, x, y, (xeps, yeps, zeps))
            if deltafn > 0:
                possible =  True
                if deltafn < epsilon:
                    epsilon = deltafn
            elif (xeps, yeps, zeps) == (0, 0, 0): # 0-d face
                possible = (delta_pi(cp, x, y) == 0)
            elif xeps!= 0 and yeps!=0 and zeps!=0: # 2-d face, 6 cases
                possible = ((delta_pi_general(cp.pwl, x, y, (xeps, yeps, zeps))==0)\
                            and (cp.find_crazy_piece(x, xeps) is None) \
                            and (cp.find_crazy_piece(y, yeps) is None) \
                            and (cp.find_crazy_piece(z, zeps) is None) )
            elif xeps == 0: # vertical 1-d face, 2 cases
                possible = (cp(x) + cp.pwl.limit(y, yeps) == cp.pwl.limit(fractional(z), zeps)) \
                           and check_move_on_crazy_pieces((1, x), cp.find_crazy_piece(y, yeps), cp.find_crazy_piece(z, zeps))
            elif yeps == 0: # horizontal 1-d face, 2 cases
                possible = (cp(y) + cp.pwl.limit(x, xeps) == cp.pwl.limit(fractional(z), zeps)) \
                           and check_move_on_crazy_pieces((1, y), cp.find_crazy_piece(x, xeps), cp.find_crazy_piece(z, zeps))
            elif zeps == 0: # diagonal 1-d face, 2 cases
                possible = (cp.pwl.limit(x, xeps) + cp.pwl.limit(y, yeps) == cp(fractional(z))) \
                           and check_move_on_crazy_pieces((-1, z), cp.find_crazy_piece(x, xeps), cp.find_crazy_piece(y, yeps))
            else:
                raise ValueError, "Something Wrong in enumerating limits"
            if not possible:
                #return (x, y, z, xeps, yeps, zeps)
                if not show_plots:
                    return 0
                else:
                    epsilon = 0
                    g += plot_limit_cone_of_vertex(x, y, epstriple_to_cone((xeps, yeps, zeps)), r=0.01)
                    g += plot_limit_cone_of_vertex(y, x, epstriple_to_cone((yeps, xeps, zeps)), r=0.01)
    range_cp = cp.range()
    epsilon = epsilon / max([abs(range_cp[0]), abs(range_cp[1])]) / 3
    if epsilon == 0 and show_plots:
        g.show(figsize=40)
    return epsilon

def check_move_on_crazy_pieces((move_sign, move_dist), cp1, cp2):
    if (cp1 is None) and (cp2 is not None) or (cp1 is not None) and (cp2 is None):
        return False
    elif (cp1 is not None) and (cp2 is not None):
        #compare if the groups on cp1 and cp2 are the same
        # TODO: set up these subgroups in a high-level way in Sage, and compare.
        if not cp1.hermite_form_generators == cp2.hermite_form_generators:
            logging.warn("Different groups. Open question.")
            return False
        if move_sign == 1:
            return (all((s == cp2(r + move_dist)) for (r, s) in cp1.cosets) \
                    and all((s == cp1(r - move_dist)) for (r, s) in cp2.cosets))
        else: # move_sign == -1:
            return (all((s == - cp2(move_dist - r)) for (r, s) in cp1.cosets) \
                    and all((s == - cp1(move_dist - r)) for (r, s) in cp2.cosets))
    else: # (cp1 is None) and (cp2 is None)
        return True
