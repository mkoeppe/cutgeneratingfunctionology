def generate_type_1_vertices_general(fn, comparison, continuity=True):
    """A generator...
    "...'general' refers to the fact that it outputs 6-tuples (x,y,z,xeps,yeps,zeps).
    """
    bkpt = fn.end_points()
    if not continuity:
        limits = fn.limits_at_end_points()
        limits[0][-1] = limits[-1][-1]
        limits[-1][1] = limits[0][1]
    for i in range(len(bkpt)):
        for j in range(i,len(bkpt)):
            x = bkpt[i]
            y = bkpt[j]
            z = fractional(x + y)
            if comparison(fn.values_at_end_points()[i] + fn.values_at_end_points()[j], fn(z)):
                yield (x, y, x+y, 0, 0, 0)
            if not continuity:
                limits_x = limits[i]
                limits_y = limits[j]
                if z == 0:
                    # take care of limit at 0-
                    limits_z = limits[0]
                elif z == 1:
                    # take care of limit at 1+
                    limits_z = limits[-1]
                else:
                    limits_z = fn.limits(z)
                zeps_set = [0, 1, -1]
                if limits_x[0] == limits_x[1] == limits_x[-1] and limits_y[0] == limits_y[1] == limits_y[-1]:
                    # continuous at x and y
                    for zeps in [1, -1]:
                        # note that (0 0 0) has alreadly been checked.
                        if comparison(limits_x[0] + limits_y[0] - limits_z[eps], 0):
                            yield (x, y, x+y, zeps, zeps, zeps)
                else:
                    for eps in [1, -1]:
                        # (+ + +), (- - -). note that (0 0 0) has alreadly been checked.
                        if comparison(limits_x[eps] + limits_y[eps] - limits_z[eps], 0):
                            yield (x, y, x+y, eps, eps, eps)
                    for eps in [1, -1]:
                        for zeps in zeps_set:
                            # (+ - 0), (+ - +), (+ - -), (- + 0), (- + +), (- + -)
                            if comparison(limits_x[eps] + limits_y[-eps] - limits_z[zeps], 0):
                                yield (x, y, x+y, eps, -eps, zeps)
                        if comparison(limits_x[0] + limits_y[eps] - limits_z[eps], 0):
                            # (0 + +), (0 - -)
                            yield (x, y, x+y, 0, eps, eps)
                        if comparison(limits_x[eps] + limits_y[0] - limits_z[eps], 0):
                            # (+ 0 +), (- 0 -)
                            yield (x, y, x+y, eps, 0, eps)

def generate_type_2_vertices_general(fn, comparison, continuity=True):
    bkpt = fn.end_points()
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]
    if not continuity:
        limits = fn.limits_at_end_points()
        limits[0][-1] = limits[-1][-1]
        limits[-1][1] = limits[0][1]
    for i in range(len(bkpt)):
        for k2 in range(i + 1, i + len(bkpt) - 1):
            # only need to check for 0 < y < 1. and note that bkpt2[i + len(bkpt) - 1] == bkpt[i] + 1.
            x = bkpt[i]
            z = bkpt2[k2]
            y = z - x
            if k2 < len(bkpt):
                k = k2
            else:
                k = k2 - len(bkpt) + 1
            if comparison(fn.values_at_end_points()[i] + fn(y), fn.values_at_end_points()[k]):
                yield (x, y, z, 0, 0, 0)
            if not continuity:
                limits_x = limits[i]
                limits_z = limits[k]
                limits_y = fn.limits(y)
                # no trouble at 0- and 1+ since 0 < y < 1.
                if not (limits_y[0] == limits_y[1] == limits_y[-1]):
                    # then y is a in bkpt. this is done in type1check_general.
                    continue
                for xeps in [0, 1, -1]:
                    for zeps in [0, 1, -1]:
                        if comparison(limits_x[xeps] + limits_y[0], limits_z[zeps]):
                            yield (x, y, z, xeps, 0, zeps)

@cached_function
def generate_nonsubadditive_vertices_general(fn, continuity=True):
    """
    We are returning a set of 6-tuples (x, y, z, xeps, yeps, zeps),
    so that duplicates are removed, and so the result can be cached for later use.
    """
    return { (x, y, z, xeps, yeps, zeps)
             for (x, y, z, xeps, yeps, zeps) in itertools.chain(generate_type_1_vertices_general(fn, operator.lt, continuity),
                                                                generate_type_2_vertices_general(fn, operator.lt, continuity))
           }

def subadditivity_check_general(fn, continuity=None):
    """
    Check if fn is subadditive. Works for discontinuous functions as well.
    """
    bkpt = fn.end_points()
    if continuity == None:
        continuity = fn.is_continuous_defined()
    result = True
    for (x, y, z, xeps, yeps, zeps) in generate_nonsubadditive_vertices_general(fn, continuity):
        logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(xeps), y, print_sign(yeps), z, print_sign(zeps)))
        result = False
    if result:
        logging.info("pi is subadditive.")
    else:
        logging.info("Thus pi is not subadditive.")
    return result

def symmetric_check_general(fn, f, continuity=None):
    """
    Check if fn is symmetric. Works for discontinuous functions as well.
    """
    result = True
    if fn(f) != 1:
        logging.info('pi(f) is not equal to 1.')
        result = False
    bkpt = fn.end_points()
    if continuity == None:
        continuity = fn.is_continuous_defined()
    if not continuity:
        limits = fn.limits_at_end_points()
        limits[0][-1] = limits[-1][-1]
        limits[-1][1] = limits[0][1]
    for i in range(len(bkpt)):
        x = bkpt[i]
        if x == f:
            continue
        if x < f:
            y = f - x
        else:
            y = 1 + f - x
        if fn.values_at_end_points()[i] + fn(y) != 1:
            logging.info('pi(%s) + pi(%s) is not equal to 1' % (x, y))
            result = False
        if not continuity:
            limits_x = limits[i]
            limits_y = fn.limits(y)
            # no trouble at 0- and 1+ since 0 < y < 1.
            if limits_x[-1] + limits_y[1] != 1:
                logging.info('pi(%s-) + pi(%s+) is not equal to 1' % (x, y))
                result = False
            if limits_x[1] + limits_y[-1] != 1:
                logging.info('pi(%s+) + pi(%s-) is not equal to 1' % (x, y))
                result = False
    if result:
        logging.info("pi is symmetric.")
    else:
        logging.info("Thus pi is not symmetric.")
    return result

def minimality_test_general(fn, f=None):
    """
    Check if fn is minimal with respect to f. Works for discontinuous functions as well.
    """
    for x in fn.values_at_end_points():
        if (x < 0) or (x > 1):
            logging.info('pi is not minimal because it does not stay in the range of [0, 1].')
            return False
    if f==None:
        f = find_f(fn, no_error_if_not_minimal_anyway=True)
        if f==None:
            return False
    if fn(0) != 0:
        logging.info('pi is NOT minimal because pi(0) is not equal to 0.')
        return False
    logging.info('pi(0) = 0')
    bkpt = fn.end_points()
    continuity = fn.is_continuous_defined()
    if not continuity:
        limits = fn.limits_at_end_points()
        limits[0][-1] = limits[-1][-1]
        limits[-1][1] = limits[0][1]
        for x in limits:
            if not ((0 <= x[-1] <=1) and (0 <= x[1] <=1)):
                logging.info('pi is not minimal because it does not stay in the range of [0, 1].')
                return False
    if subadditivity_check_general(fn, continuity) and symmetric_check_general(fn, f, continuity):
        logging.info('Thus pi is minimal.')
        return True
    logging.info('Thus pi is NOT minimal.')
    return False
