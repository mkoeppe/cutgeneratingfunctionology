## Old code for random walks and BFS.
## This is superseded by directed_move_composition_completion.

def find_possible_directed_moves(x,directed_moves):
    """Find the directed moves applicable at x.
    """
    possible_directed_moves = []
    for directed_move in directed_moves:
        if is_directed_move_possible(x, directed_move):
            possible_directed_moves.append(directed_move)
    return possible_directed_moves

def find_impossible_directed_moves(x,directed_moves):
    """
    Find the directed moves NOT applicable at x.
    """
    impossible_directed_moves = []
    for directed_move in directed_moves:
        if not is_directed_move_possible(x, directed_move):
            impossible_directed_moves.append(directed_move)
    return impossible_directed_moves

def apply_directed_move(x, directed_move):
    return directed_move(x)


def random_walk(seed, moves, fn, num_it):
    """
    Find the orbit of a seed through a random walk.  It may not find
    all the possible elements in the orbit (if unlucky).  It may not
    be very efficient if the size of orbit is small because it will
    continue to run until the predetermined number of iterations is
    reached.

    deterministic_walk is preferable in most cases.
    """
    import random
    x = seed
    xlist = {x:[1,None,None]}
    walk_sign = 1
    # points_plot = point((x,0))
    for i in range(num_it):
        possible_directed_moves = find_possible_directed_moves(x,moves)
        move = possible_directed_moves[random.randint(0,len(possible_directed_moves)-1)]
        move_sign = move.sign()
        directed_move = move
        
        next_x = move(x)

        next_walk_sign = move_sign * walk_sign
        if next_x not in xlist:
            pass
            # points_plot += point((next_x,0))
        else:
            if xlist[next_x][2] != None:
                x = xlist[next_x][1]
                directed_move = xlist[next_x][2]
            if xlist[next_x][0] != next_walk_sign:
                next_walk_sign = 0
                print next_x, ","
                previous = xlist[next_x][1]
                while previous != next_x and previous != None:
                    print previous, ","
                    previous = xlist[previous][1]
                print next_x, "."
        # Prep for next
        xlist[next_x] = [next_walk_sign, x, directed_move] 
        walk_sign = next_walk_sign
        x = next_x
    # show(points_plot + plot(-0.1+x-x, [A,A0]) + plot(-0.1+x-x, [f-A0,f-A])\
    #     + plot(-0.2+x-x, [A,A+a0], color = "green") + plot(-0.3+x-x, [A,A+a1], color = "green") + \
    #     plot(-0.4+x-x, [A,A+a2], color = "green"), ymin = -1, ymax = 0)
    return xlist

class MaximumNumberOfIterationsReached(Exception):
    pass

class SignContradiction(Exception):
    pass

def deterministic_walk(seed, moves, fn=None, max_num_it = 1000, error_if_sign_contradiction=False, error_if_max_num_it_exceeded=True):
    """
    Compute the orbit of a given seed. (Done by a breadth-first search.)
    To avoid infinite computations in the case of a dense orbit,
    there is a maximum number of iterations (by default, it is 1000).

    The moves know the domains on which they are allowed.

    Returns a dictionary:
    - keys are elements of the orbit
    - values are lists of the form [walk_sign, predecessor, directed_move_from_predecessor].

    If `error_if_max_num_it_exceeded` is `False`, then 
    a secondary result is the to_do list.
    """
    ## FIXME: If `fn` is provided, store the dictionary in `fn.walk_dict`
    ## and the to_do list in `fn.walk_to_do`, to allow us to resume
    ## an interrupted BFS.  OR: Turn it into an iterator/generator (yield).

    logging.info("Breadth-first search to discover the reachable orbit...")
    to_do = [seed]
    # xlist is actually a dictionary.
    
    xlist = {seed:[1,None,None]}
    walk_sign = 1
    # points_plot = point((seed,0))
    contradiction_reached = False
    num_it = 0
    while to_do and num_it < max_num_it:
        if (num_it > 0 and num_it % 100 == 0):
            logging.info("(Iteration %d, to do list has %d items)" % (num_it, len(to_do)))
        x = to_do.pop(0)
        possible_directed_moves = find_possible_directed_moves(x, moves)
        for directed_move in possible_directed_moves:
            walk_sign = xlist[x][0]
            move_sign = directed_move.sign()
            move_for_next_x = directed_move
            next_x = directed_move(x)
            
            next_walk_sign = move_sign * walk_sign
            
            if next_x not in xlist:
                # points_plot += point((next_x,0))
                xlist[next_x] = [next_walk_sign, x, move_for_next_x]
                #if next_x not in to_do:
                to_do.append(next_x)
                
            else:
                # If next_x is already in xlist, we do not need to modify anything, except probably the sign.
                if xlist[next_x][0] != next_walk_sign:
                    xlist[next_x][0] = 0
                    if contradiction_reached == False:
                        logging.info('A contradiction of signs was reached. All the elements in the reachable orbit take the value 0.')
                        if error_if_sign_contradiction:
                            raise SignContradiction
                        else:
                            contradiction_reached = True
                            for pt in xlist.keys():
                                xlist[pt][0] = 0
        num_it = num_it + 1

    if error_if_max_num_it_exceeded:
        if num_it == max_num_it:
            raise MaximumNumberOfIterationsReached, "Reached %d iterations, to do list has still %d items" % (num_it, len(to_do))
        logging.info("Breadth-first search to discover the reachable orbit... done")
        return xlist
    else:
        if num_it == max_num_it:
            logging.info("Breadth-first search to discover the reachable orbit reached %d iterations, to do list has still %d items" % (num_it, len(to_do)))
        return xlist, to_do

def find_generic_seed(fn, max_num_it = 1000):
    intervals = generate_uncovered_intervals(fn)
    if not intervals:
        raise ValueError, "Need an uncovered interval"
    moves = generate_functional_directed_moves(fn)
    seed = intervals[0][1]
    while True:
        seed = 2/3 * intervals[0][0] + 1/3 * seed
        logging.info("Trying seed %s" % seed)
        try:
            stab_int, walk_list = find_stability_interval_with_deterministic_walk_list \
                                  (seed, intervals, moves, fn, max_num_it = max_num_it, \
                                   error_if_sign_contradiction = True)
            if not stab_int[0] < 0 < stab_int[1]:
                logging.info("Stability interval does not contain 0 in its interior, continuing search.")
                continue
            logging.info("Seed %s has a proper stability interval %s, reachable orbit has %s elements" % (seed, stab_int, len(walk_list)))
            return (seed, stab_int, walk_list)
        except SignContradiction:
            continue

def find_stability_interval_with_deterministic_walk_list(seed, intervals, moves, fn, max_num_it = 1000, error_if_sign_contradiction=False):
    """
    Returns the stability interval (an open, half-open, or closed interval)
    and the deterministic_walk_list.
    """
    walk_dict = deterministic_walk(seed, moves, fn, max_num_it, \
                                   error_if_sign_contradiction=error_if_sign_contradiction)
    return (compute_stability_interval(walk_dict, intervals, moves, fn), walk_dict)

def generate_shifted_stability_intervals(stab_int, walk_list):
    """The result is sorted."""
    orbit = sorted(walk_list.keys())
    intervals = []
    for i in orbit:
        sign = walk_list[i][0]
        if sign == 1:
            intervals.append(closed_or_open_or_halfopen_interval(stab_int.a + i, stab_int.b + i, stab_int.left_closed, stab_int.right_closed))
        elif sign == -1:
            intervals.append(closed_or_open_or_halfopen_interval(i - stab_int.b, i - stab_int.a, stab_int.right_closed, stab_int.left_closed))
        elif sign == 0:
            assert(stab_int.a == -stab_int.b)
            assert(stab_int.left_closed == stab_int.right_closed)
            intervals.append(closed_or_open_or_halfopen_interval(stab_int.a + i, stab_int.b + i, stab_int.left_closed, stab_int.right_closed))
    return intervals

def find_decomposition_into_stability_intervals(fn):
    ## experimental.
    intervals = generate_uncovered_intervals(fn)
    moves = generate_functional_directed_moves(fn)
    return iterative_stability_refinement(intervals, moves)

def find_decomposition_into_stability_intervals(fn, show_plots=False, max_num_it=1000):
    fn._stability_orbits = []
    uncovered_intervals = generate_uncovered_intervals(fn)
    intervals = uncovered_intervals
    moves = generate_functional_directed_moves(fn)
    orbits = []
    while intervals:
        #print "Intervals: ", intervals
        seed = intervals[0][0] + (intervals[0][1] - intervals[0][0]) / 3
        print "Seed: ", seed
        walk_dict, to_do = deterministic_walk(seed, moves, fn, max_num_it = max_num_it, \
                                              error_if_sign_contradiction=False, 
                                              error_if_max_num_it_exceeded=False)
        # When to_do is nonempty, the BFS did not finish.
        # We compute the stability interval anyway -- it gives an interval of points that ALL have such a long orbit!
        # But we don't try to make the stability orbits disjoint in this case (this is not well-defined);
        # instead we count on our scan functions to deal with the "multiset" case properly.
        make_disjoint = not to_do
        stab_int = compute_stability_interval(walk_dict, uncovered_intervals, moves, fn, make_disjoint = make_disjoint)
        orbits.append((stab_int, walk_dict))
        if show_plots:
            plot_orbit_comparison(fn, orbits).show(dpi=500)
        #print "Stability interval: ", stab_int
        shifted_stability_intervals = generate_shifted_stability_intervals(stab_int, walk_dict)
        print "Stability orbit: ", shifted_stability_intervals[0], ", ... (length ", len(walk_dict), ")"
        fn._stability_orbits.append((shifted_stability_intervals, walk_dict, to_do))
        remaining = union_of_coho_intervals_minus_union_of_coho_intervals([intervals], [shifted_stability_intervals])
        intervals = remaining
        
    logging.info("Total: %s stability orbits, lengths: %s" \
                 % (len(fn._stability_orbits), \
                    [ ("%s+" if to_do else "%s") % len(shifted_stability_intervals) \
                      for (shifted_stability_intervals, walk_dict, to_do) in fn._stability_orbits ]))



def compute_stability_interval(deterministic_walk_list, intervals, moves, fn, make_disjoint = True):
    ## FIXME: Refactor using above.
    a = -10
    b = 10
    left_closed = True
    right_closed = True
    for pt in deterministic_walk_list.keys():
        if deterministic_walk_list[pt][0] == 0:
            return closed_or_open_or_halfopen_interval(0, 0, True, True)
        for interval in intervals:
            if element_of_int(pt, interval):
                if deterministic_walk_list[pt][0] == 1:
                    if (interval[0] - pt) > a:
                        a = interval[0] - pt
                        left_closed = True
                    if (interval[1] - pt) < b:
                        b = interval[1] - pt
                        right_closed = True
                elif deterministic_walk_list[pt][0] == -1:
                    if (pt - interval[0]) < b:
                        b = pt - interval[0]
                        right_closed = True
                    if (interval[1] - pt) < -a:
                        a = pt - interval[1]
                        left_closed = True      
        impossible_directed_moves = find_impossible_directed_moves(pt, moves)
        ### We now take the set difference of
        ### the __directed__ moves with the possible directed moves.
        ### This was not done in the old code: --Matthias
        # impossible_directed_moves = []
        # for move in moves:
        #     if move not in possible_directed_moves:
        #         impossible_directed_moves.append(move)
        
        for move in impossible_directed_moves:
            if move.sign() == 1:
                impossible_next_x = fractional(pt + move[1])
                for interval in intervals:
                    temp = interval[0] - impossible_next_x
                    if 0 < temp <= b:
                        b = temp
                        right_closed = False
                    temp2 = interval[1] - impossible_next_x
                    if a <= temp2 < 0:
                        a = temp2
                        left_closed = False
            elif move.sign() == -1:
                impossible_next_x = fractional(move[1] - pt)
                for interval in intervals:
                    temp = interval[0] - impossible_next_x
                    if 0 < temp <= -a:
                        a = -temp
                        left_closed = False
                    temp2 = impossible_next_x - interval[1]
                    if 0 < temp2 <= b:
                        b = temp2
                        right_closed = False
    if make_disjoint:
        ### Now we make the stability orbit disjoint by looking at adjacent intervals
        ### and shrinking if necessary.
        orbit = sorted(deterministic_walk_list.keys())
        min_midpt_dist = 1
        for i in range(len(orbit)-1):
            sign_1 = deterministic_walk_list[orbit[i]][0]
            sign_2 = deterministic_walk_list[orbit[i+1]][0]
            half_distance = (orbit[i+1]-orbit[i])/2
            ## Two intervals of different signs might overlap.
            if sign_1 == 1 and sign_2 == -1:
                if b > half_distance:
                    b = half_distance
                    logging.info("Making stability intervals disjoint: Reducing right to %s to separate %s and %s" % (b, orbit[i], orbit[i+1]))
                    right_closed = False
            elif sign_1 == -1 and sign_2 == 1:
                if -a > half_distance:
                    a = -half_distance
                    logging.info("Making stability intervals disjoint: Reducing left to %s to separate %s and %s" % (a, orbit[i], orbit[i+1]))
                    left_closed = False
            else:
                assert(-a + b <= 2 * half_distance)
    return closed_or_open_or_halfopen_interval(a, b, left_closed, right_closed)

def one_step_stability_interval(x, intervals, moves):
    """Returns the stability interval, i.e., an open, half-open, or closed
    interval I), such that for all moves, `move`(`x`) lies in
    `intervals` if and only if `move`(`x` + t) lies in `intervals` for
    t in I."""
    a = -10
    b = 10
    left_closed = True
    right_closed = True
    for move in moves:
        next_x = directed_move(x)
        for interval in intervals:
            if (interval[0] <= next_x <= interval[1]):
                # Move was possible, so:
                if move[0] == 1:
                    if (interval[0] - next_x) > a:
                        a = interval[0] - next_x
                        left_closed = True
                    if (interval[1] - next_x) < b:
                        b = interval[1] - next_x
                        right_closed = True
                elif move[0] == -1:
                    if (next_x - interval[0]) < b:
                        b = next_x - interval[0]
                        right_closed = True
                    if (interval[1] - next_x) < -a:
                        a = next_x - interval[1]
                        left_closed = True      
                else:
                    raise ValueError, "Move not valid: %s" % list(move)
                break # next_x can only lie in one interval 
        else:
            # Move was not possible.
            if move[0] == 1:
                for interval in intervals:
                    temp = interval[0] - next_x
                    if 0 < temp <= b:
                        b = temp
                        right_closed = False
                    temp2 = interval[1] - next_x
                    if a <= temp2 < 0:
                        a = temp2
                        left_closed = False
            elif move[0] == -1:
                for interval in intervals:
                    temp = interval[0] - next_x
                    if 0 < temp <= -a:
                        a = -temp
                        left_closed = False
                    temp2 = next_x - interval[1]
                    if 0 < temp2 <= b:
                        b = temp2
                        right_closed = False
            else:
                raise ValueError, "Move not valid: %s" % list(move)
    return closed_or_open_or_halfopen_interval(a, b, left_closed, right_closed)

def one_step_stability_refinement(interval, intervals, moves):
    if len(interval) == 2:
        interval = closed_or_open_or_halfopen_interval(interval[0], interval[1], True, True)
    seed = (interval[0] + interval[1]) / 2
    (stab_a, stab_b, stab_left_closed, stab_right_closed) = one_step_stability_interval(seed, intervals, moves)
    stab_a += seed
    stab_b += seed
    left = None
    right = None
    #print "  Shifted stability interval: ", stab_a, stab_b, stab_left_closed, stab_right_closed
    #print interval.a < stab_a, interval.left_closed
    if interval.a < stab_a \
       or (interval.a == stab_a and interval.left_closed and not stab_left_closed):
        left = closed_or_open_or_halfopen_interval(interval.a, stab_a, interval.left_closed, not stab_left_closed)
        #print "  Left: ", left
    else:
        stab_a = interval.a
        stab_left_closed = interval.left_closed
    if stab_b < interval.b \
       or (stab_b == interval.b and interval.right_closed and not stab_right_closed):
        right = closed_or_open_or_halfopen_interval(stab_b, interval.b, not stab_right_closed, interval.right_closed)
        #print "  Right: ", right
    else:
        stab_b = interval.b
        stab_right_closed = interval.right_closed
    #print "  Modified stability interval: ", stab_a, stab_b, stab_left_closed, stab_right_closed
    refinement = []
    if left != None:
        refinement.append(left)
    refinement.append(closed_or_open_or_halfopen_interval(stab_a, stab_b, stab_left_closed, stab_right_closed))
    if right != None:
        refinement.append(right)
    return refinement

from heapq import *

def iterative_stability_refinement(intervals, moves):
    interval_pq = [ (-interval_length(interval), interval) \
                    for interval in intervals ]
    heapify(interval_pq)
    finished_list = []
    while (interval_pq):
        priority, int = heappop(interval_pq)
        length = -priority
        logging.info("%s of length %s " % (int, float(length)))
        refinement = one_step_stability_refinement(int, intervals, moves)
        logging.debug("  Refinement: %s" % refinement)
        if len(refinement) == 1:
            finished_list.append(refinement[0])
        else:
            for new_int in refinement:
                heappush(interval_pq, (-interval_length(new_int), new_int))
    return sorted(finished_list + \
                  [ interval for (priority, interval) in interval_pq ])

def plot_orbit_comparison(fn, orbits):
    """
    `orbits` can be a list of:
    tuples (stab_int, walk_dict)
    or numbers (seeds)
    """
    g = plot_intervals(generate_uncovered_intervals(fn))
    ymin = 0.5
    ymax = 1
    for orbit, color in itertools.izip(orbits, rainbow(len(orbits), 'rgbtuple')):
        ymax -= 0.2
        ymin -= 0.2
        stab_color = Color(color).lighter() 
        try:
            stab_int, walk_dict = orbit
        except TypeError:
            seed = orbit
            stab_int, walk_dict = find_stability_interval_with_deterministic_walk_list(seed, generate_uncovered_intervals(fn), generate_functional_directed_moves(fn), fn)
        g += plot_walk_with_shifted_stability_intervals(stab_int, walk_dict, stab_color=stab_color, color=color, ymin=ymin, ymax=ymax)
    return g

def plot_shifted_stability_intervals(shifted_stability_intervals, color="cyan", ymin=0, ymax=1, **options):
    g = Graphics()
    for i in shifted_stability_intervals:
        g += polygon2d([(i.a, ymin), (i.b, ymin), (i.b, ymax), (i.a, ymax)], \
                       color=color, zorder = -5, **options)
    return g

def plot_walk_with_shifted_stability_intervals(stab_int, walk_dict, color="black", stab_color="cyan", **options):
    shifted_stability_intervals = generate_shifted_stability_intervals(stab_int, walk_dict)
    return plot_shifted_stability_intervals(shifted_stability_intervals, color=stab_color, **options) \
        + plot_walk(walk_dict, color=color, **options)

