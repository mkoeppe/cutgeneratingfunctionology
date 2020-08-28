from six.moves import range
# Reminder: need coerce all input to common RNF.

class CrazyPiece:
    def __init__(self, interval, generators=None, cosets=[]):
        # assume interval is open, represented as (a, b)
        self.interval = interval
        # assume that generators QQ linearly independent
        self.generators = generators
        if not is_QQ_linearly_independent(*generators):
            logging.warning("Generators are not linearly independent over Q.")
        if len(generators) < 2:
            logging.warning("The group is not dense.")
        self.hermite_form_generators = find_hermite_form_generators(generators)
        # cosets is a list of (coset_repr, shift),
        # assume coset_repr represent distinct cosets;
        # assume shifts are non-zero
        self.cosets = cosets
        for i in range(len(cosets)):
            (r, s) = cosets[i]
            if s == 0:
                logging.warning("Have shift = 0.")
            for j in range(i):
                rr = cosets[j][0]
                if is_in_ZZ_span(r-rr, generators):
                    logging.warning("Not unique coset representative.")

    def __call__(self, x):
        x = fractional(x)
        ###assert self.interval[0] < x < self.interval[1]
        for (r,s) in self.cosets:
            if is_in_ZZ_span(x-r, self.generators):
                return s
        return x.parent().zero()

    def __neg__(self):
        new_cosets = [(r, -s) for (r, s) in self.cosets]
        return CrazyPiece(self.interval, self.generators, new_cosets)

    def __mul__(self, other):
        r"""
        Multiply self by a scalar.
        """
        # assume scalar multiplication
        new_cosets = [(r, s*other) for (r, s) in self.cosets]
        return CrazyPiece(self.interval, self.generators, new_cosets)

    __rmul__ = __mul__

    def __div__(self, other):
        return self * (1 / other)

    def range(self):
        shifts =  [s for (r, s) in self.cosets]
        return unique_list(shifts+[0])

from sage.structure.element import Element, ModuleElement

class PiecewiseCrazyFunction(ModuleElement):

    # assume that all inputs are elements from a same RNF.
    def __init__(self, pwl, crazy_pieces, parent=None):
        self.pwl = pwl
        # assume that crazy pieces' intervals are disjoint.
        self.crazy_pieces = crazy_pieces
        if parent is None:
            from sage.structure.element import get_coercion_model
            cm = get_coercion_model()
            crazy_generators = [generator.parent()
                                for piece in crazy_pieces
                                for generator in piece.generators]
            domain_parent = cm.common_parent(pwl.parent().domain(), *crazy_generators)
            from sage.categories.pushout import pushout
            codomain_parent = pushout(domain_parent, pwl.parent().codomain())
            parent = PiecewiseCrazyFunctionsSpace(domain_parent, codomain_parent)
        ModuleElement.__init__(self, parent)

    def _repr_(self):
        return "PiecewiseCrazyFunction({}, crazy_pieces={})".format(self.pwl, self.crazy_pieces)

    def _richcmp_(self, other, op):
        if not isinstance(other, PiecewiseCrazyFunction):
            return NotImplemented
        def is_equal():
            # FIXME: Need better comparison -- there are many ways to express the same function
            return self.pwl == other.pwl and self.crazy_pieces == other.crazy_pieces
        if op == op_EQ:
            return is_equal()
        elif op == op_NE:
            return not is_equal()
        else:
            # FIXME: Implement partial order of pointwise comparison
            return NotImplemented

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

    def _add_(self,other):
        # assume that intervals of crazy pieces from self and from other are disjoint.
        # FIXME: This needs to be checked
        return self.__class__(self.pwl + other.pwl, self.crazy_pieces + other.crazy_pieces, parent=self.parent())

    def _neg_(self):
        return self.__class__(-self.pwl, [-cp  for cp in self.crazy_pieces], parent=self.parent())

    def _lmul_(self, scalar):
        r"""
        Multiply self by a scalar.
        """
        return self.__class__(self.pwl * scalar, [cp * scalar for cp in self.crazy_pieces], parent=self.parent())

    _rmul_ = _lmul_

    # Adding the following method helps the coercion system
    # discover the action of multiplication with a scalar that is not the base ring.
    def _acted_upon_(self, actor, self_on_left):
        r"""
        Multiply self by a matrix from left or right.
        """
        return self.__class__(self.pwl * actor, [cp * actor for cp in self.crazy_pieces],
                              parent=self.parent())

    def __truediv__(self, other):
        return self * (QQ(1) / other)

    __div__ = __truediv__

    def _sub_(self, other):
        return self + (-other)

    @cached_method
    def __call__(self, x):
        x = fractional(x)
        crazy_piece = self.find_crazy_piece(x, 0)
        if crazy_piece is None:
            return self.pwl(x)
        else:
            return self.pwl(x) + crazy_piece(x)

    def plot(self, rgbcolor='magenta'):
        pwl = self.pwl
        g = pwl.plot(rgbcolor=rgbcolor)
        for crazy_piece in self.crazy_pieces:
            a, b = crazy_piece.interval[0], crazy_piece.interval[1]
            pwla = pwl.limit(a, 1)
            pwlb = pwl.limit(b, -1)
            shifts = [s for (r, s) in crazy_piece.cosets]
            max_shift = max(shifts)
            min_shift = min(shifts)
            g += polygon([(a, pwla+min_shift), (a, pwla+max_shift), (b, pwlb+max_shift), (b, pwlb+min_shift)], color=rgbcolor, alpha=0.1)
            for s in shifts:
                g += line([(a, pwla+s), (b, pwlb+s)], color=rgbcolor)
        return g

    def end_points(self):
        bkpts = copy(self.pwl.end_points())
        for crazy_piece in self.crazy_pieces:
            bkpts += [crazy_piece.interval[0], crazy_piece.interval[1]]
        return sorted(set(bkpts))

    def limit(self, x0, epsilon):
        x0 = fractional(x0)
        if epsilon != 0:
            if not (self.find_crazy_piece(x0, epsilon) is None):
                raise ValueError("The crazy function has no limit at {}{}".format(x0, print_sign(epsilon)))
            else:
                return self.pwl.limit(x0, epsilon)
        else:
            return self(x0)

    def limit_range(self, x0, epsilon):
        x0 = fractional(x0)
        crazy_piece = self.find_crazy_piece(x0, epsilon)
        limit_pwl = self.pwl.limit(x0, epsilon)
        if crazy_piece is None:
            return [limit_pwl]
        else:
            return [limit_pwl + s for s in crazy_piece.range()]

from .fast_piecewise import PiecewiseFunctionsSpace, PiecewiseLinearFunctionsSpace

class PiecewiseCrazyFunctionsSpace(UniqueRepresentation, PiecewiseFunctionsSpace):

    """
    Space of piecewise crazy functions.

    TESTS::

        sage: from cutgeneratingfunctionology.igp import *
        sage: pushout(PiecewiseCrazyFunctionsSpace(QQ, QQ), PiecewiseLinearFunctionsSpace(AA, AA))  # known bug -- need to implement construction functor!

        sage: TestSuite(PiecewiseCrazyFunctionsSpace(AA, AA)).run(skip=['_test_zero'])
        sage: TestSuite(PiecewiseCrazyFunctionsSpace(AA, AA)).run(verbose=True, catch=False)  # not tested - for interactive use

    """
    ## FIXME: We skip _test_zero - we need to improve equality testing to get it to work.

    Element = PiecewiseCrazyFunction

    def _repr_(self):
        return "Vector space of piecewise crazy partial functions from {} to {} over {}".format(
            self.domain(), self.codomain(), self.base_ring())

    def _coerce_map_from_(self, S):
        """
        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: PiecewiseCrazyFunctionsSpace(AA, AA).has_coerce_map_from(PiecewiseLinearFunctionsSpace(QQ, QQ))
            True
        """
        if isinstance(S, PiecewiseLinearFunctionsSpace) or isinstance(S, PiecewiseCrazyFunctionsSpace):
            return self.codomain().has_coerce_map_from(S.codomain())
        return False

    def _element_constructor_(self, x, *args, **kwds):
        element_class = self.element_class
        if isinstance(x, PiecewiseCrazyFunction):
            if args or kwds:
                raise TypeError("no extra args allowed")
            return element_class(x.pwl, x.crazy_pieces, parent=self)
        if isinstance(x, FastPiecewise):
            if not args and not kwds:
                args = [[]]
            return element_class(x, parent=self, *args, **kwds)
        if x in self.base_ring():
            pwl = PiecewiseLinearFunctionsSpace(self.domain(), self.codomain())(x)
            return element_class(pwl, [], parent=self)
        return element_class(x, parent=self, *args, **kwds)

    def _an_element_(self):
        pwl = PiecewiseLinearFunctionsSpace(self.domain(), self.codomain())._an_element_()
        return self(pwl, [])

def is_in_ZZ_span(x, generators):
    # assume that all inputs are elements from a same RNF.
    # generators are linearly independent over Q
    numbers = [x]+generators
    if not is_all_the_same_number_field_fastpath(numbers):
        numbers = nice_field_values(numbers, RealNumberField)
        if not is_NumberFieldElement(numbers[0]):
            if is_all_QQ(numbers):
                raise ValueError("generators are not Q-linear independent")
            raise ValueError("Q-linear independence test only implemented for algebraic numbers")
    x = numbers[0]
    generators = numbers[1::]
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
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = kzh_minimal_has_only_crazy_perturbation_1()
        sage: bkpts = h.end_points()
        sage: t1 = bkpts[10]-bkpts[6]
        sage: t2 = bkpts[13]-bkpts[6]
        sage: f = bkpts[37]
        sage: ucl = bkpts[17]
        sage: ucr = bkpts[18]
        sage: generators = [t1, t2]
        sage: pwl = piecewise_function_from_breakpoints_and_slopes([0,1],[0])
        sage: crazy_piece_1 = CrazyPiece((ucl, ucr), generators, [(ucl, 1), (ucr, -1)])
        sage: crazy_piece_2 = CrazyPiece((f-ucr, f-ucl), generators, [(f-ucr, 1), (f-ucl, -1)])
        sage: cp = PiecewiseCrazyFunction(pwl, [crazy_piece_1, crazy_piece_2])
        sage: find_epsilon_for_crazy_perturbation(h, cp)
        0.0003958663221935161?

        sage: h = minimal_no_covered_interval()
        sage: cp = equiv7_example_2_crazy_perturbation()
        sage: find_epsilon_for_crazy_perturbation(h, cp)
        1/6
    """
    # assume fn is a subadditive pwl function, cp (crazy perturbation) is a non_zero PiecewiseCrazyFunction with cp(0)=cp(f)=0.
    bkpt = sorted(set(copy(fn.end_points())+cp.end_points()))
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]
    type_1_vertices = [(x, y, x+y) for x in bkpt for y in bkpt if x <= y]
    type_2_vertices = [(x, z-x, z) for x in bkpt for z in bkpt2 if x < z < 1+x]
    vertices = set(type_1_vertices + type_2_vertices)
    g = plot_2d_complex(fn)
    m = 3
    M = 0
    for (x, y, z) in vertices:
        for (xeps, yeps, zeps) in [(0,0,0)]+list(nonzero_eps):
            deltafn = delta_pi_general(fn, x, y, (xeps, yeps, zeps))
            if deltafn > 0:
                possible =  True
                if deltafn < m:
                    m = deltafn
                xcp_range = cp.limit_range(x, xeps)
                ycp_range = cp.limit_range(y, yeps)
                zcp_range = cp.limit_range(z, zeps)
                deltacp_min = min(xcp_range) + min(ycp_range) - max(zcp_range)
                deltacp_max = max(xcp_range) + max(ycp_range) - min(zcp_range)
                deltacp_abs = max([abs(deltacp_min), abs(deltacp_max)])
                if deltacp_abs > M:
                    M = deltacp_abs
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
                raise ValueError("Something Wrong in enumerating limits")
            if not possible:
                #return (x, y, z, xeps, yeps, zeps)
                if not show_plots:
                    return 0
                else:
                    m = 0
                    g += plot_limit_cone_of_vertex(x, y, epstriple_to_cone((xeps, yeps, zeps)), r=0.01)
                    g += plot_limit_cone_of_vertex(y, x, epstriple_to_cone((yeps, xeps, zeps)), r=0.01)
    if M == 0:
        return 1
    epsilon = m / M # not a maximal value.
    if epsilon == 0 and show_plots:
        g.show(figsize=40)
    return epsilon

def check_move_on_crazy_pieces(xxx_todo_changeme, cp1, cp2):
    (move_sign, move_dist) = xxx_todo_changeme
    if (cp1 is None) and (cp2 is not None) or (cp1 is not None) and (cp2 is None):
        return False
    elif (cp1 is not None) and (cp2 is not None):
        #compare if the groups on cp1 and cp2 are the same
        # TODO: set up these subgroups in a high-level way in Sage, and compare.
        if not cp1.hermite_form_generators == cp2.hermite_form_generators:
            logging.warning("Different groups. Open question.")
            return False
        if move_sign == 1:
            return (all((s == cp2(r + move_dist)) for (r, s) in cp1.cosets) \
                    and all((s == cp1(r - move_dist)) for (r, s) in cp2.cosets))
        else: # move_sign == -1:
            return (all((s == - cp2(move_dist - r)) for (r, s) in cp1.cosets) \
                    and all((s == - cp1(move_dist - r)) for (r, s) in cp2.cosets))
    else: # (cp1 is None) and (cp2 is None)
        return True


def random_test_number(fn, extra_cosets=[]):
    if randint(0, 5)==0:
        # Pick f
        try:
            return find_f(fn.pwl)
        except AttributeError:
            return find_f(fn)
    breakpoints = fn.end_points()
    if randint(0, 5) == 0:
        # Pick a breakpoint
        return breakpoints[randint(0, len(breakpoints)-1)]
    # Pick a point from the interior of some interval
    crazy_pieces = []
    intervals = []
    try:
        crazy_pieces = fn.crazy_pieces
        intervals = fn.pwl.intervals()
    except AttributeError:
        intervals = fn.intervals()
    if crazy_pieces and randint(0, 1) == 0:
        # Pick from crazy piece
        crazy_piece = crazy_pieces[randint(0, len(crazy_pieces)-1)]
        if randint(0, 0) == 0: # Always...
            # Pick from support of microperiodic
            cosets = [ coset for (coset, value) in crazy_piece.cosets ] + extra_cosets
            coset = cosets[randint(0, len(cosets)-1)]
        else:
            coset = ZZ(randint(0, 12345678)) / ZZ(randint(1, 1234567))
        generators = crazy_piece.generators
        x = coset + sum(randint(0, 12345678) * gen for gen in generators)
        assert generators[0] < 1
        x = x - floor(x / generators[0] * generators[0])
        i = floor((1 - x) / generators[0])
        x = x + randint(0, i) * generators[0]
        return x
    interval = intervals[randint(0, len(intervals)-1)]
    denom = 12345678
    x = interval[0] + ZZ(randint(0, denom)) / denom * (interval[1] - interval[0])
    return x

def random_6_tuple(fn, extra_cosets=[]):
    eps_tuples = list(dic_eps_to_cone.keys()) # 13 possibilities.
    while True:
        (xeps, yeps, zeps)=eps_tuples[randint(0, 12)]
        if randint(0, 1) == 0:
            x = random_test_number(fn, extra_cosets=extra_cosets)
            y = random_test_number(fn, extra_cosets=extra_cosets)
            z = x + y
        else:
            x = random_test_number(fn, extra_cosets=extra_cosets)
            z = randint(0, 1) + random_test_number(fn, extra_cosets=extra_cosets)
            y = fractional(z - x)
        try:
            fn_x = fn.limit(x, xeps)
            fn_y = fn.limit(y, yeps)
            fn_z = fn.limit(z, zeps)
            return x, y, z, xeps, yeps, zeps
        except ValueError:
            # if fn is a crazy function such that not all of fn_x, fn_y, fn_z exist, then restart the random generator.
            pass

def minimality_test_randomized(fn, orig_function=None, testpoint_function=None, extra_cosets=[], max_iterations=None, limits=True, lost_additivity_is_error=False):
    r"""
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = kzh_minimal_has_only_crazy_perturbation_1()
        sage: bkpts = h.end_points()
        sage: t1 = bkpts[10]-bkpts[6]
        sage: t2 = bkpts[13]-bkpts[6]
        sage: f = bkpts[37]
        sage: ucl = bkpts[17]
        sage: ucr = bkpts[18]
        sage: generators = [t1, t2]
        sage: pwl = piecewise_function_from_breakpoints_and_slopes([0,1],[0])
        sage: crazy_piece_1 = CrazyPiece((ucl, ucr), generators, [(ucl, 1), (ucr, -1)])
        sage: crazy_piece_2 = CrazyPiece((f-ucr, f-ucl), generators, [(f-ucr, 1), (f-ucl, -1)])
        sage: cp = PiecewiseCrazyFunction(pwl, [crazy_piece_1, crazy_piece_2])
        sage: eps = find_epsilon_for_crazy_perturbation(h, cp)
        sage: hcp = PiecewiseCrazyFunction(h+eps*pwl, [eps*crazy_piece_1, eps*crazy_piece_2]) # In fact, hcp = h + eps * cp, but it fails in doctest.
        sage: minimality_test_randomized(hcp, h, max_iterations=10)
        True
    """
    smallest_delta = 10
    num_it = 0
    if testpoint_function is None:
        testpoint_function = fn
    f = find_f(testpoint_function)
    while max_iterations is None or num_it < max_iterations:
        num_it = num_it + 1
        x, y, z, xeps, yeps, zeps = random_6_tuple(testpoint_function, extra_cosets=extra_cosets)
        if limits:
            delta = delta_pi_general(fn, x, y, (xeps, yeps, zeps))
            if orig_function is not None:
                delta_orig = delta_pi_general(orig_function, x, y, (xeps, yeps, zeps))
        else:
            xeps = yeps = zeps = 0
            delta = delta_pi(fn, x, y)
            if orig_function is not None:
                delta_orig = delta_pi(orig_function, x, y)
        if delta < 0:
            logging.warning("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(xeps), y, print_sign(yeps), z, print_sign(zeps)))
            return False
        if 0 < delta and (x + y == f or x + y == 1 + f) and zeps == 0:
            # Symmetry violated
            logging.warning("pi(%s%s) + pi(%s%s) - pi(%s%s) != 0" % (x, print_sign(xeps), y, print_sign(yeps), z, print_sign(zeps)))
            return False
        if 0 < delta and orig_function is not None:
            if delta_orig == 0:
                logging.warning("Lost additivity: pi(%s%s) + pi(%s%s) - pi(%s%s) > 0" % (x, print_sign(xeps), y, print_sign(yeps), z, print_sign(zeps)))
                if lost_additivity_is_error:
                    return False
        if 0 == delta and orig_function is not None:
            if delta_orig != 0:
                logging.info("New additivity: pi(%s%s) + pi(%s%s) - pi(%s%s) = 0" % (x, print_sign(xeps), y, print_sign(yeps), z, print_sign(zeps)))
        if 0 < delta < smallest_delta:
            smallest_delta = delta
            logging.info("After {} tries, smallest Delta pi now: {}".format(num_it, delta))
    return True
