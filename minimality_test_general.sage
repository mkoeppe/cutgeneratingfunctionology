def type1check_general(fn, continuity=True, limits=None):
    endpts = fn.end_points()
    for i in range(len(endpts)):
        for j in range(i,len(endpts)):
            x = endpts[i]
            y = endpts[j]
            z = fractional(x + y)
            if fn._values_at_end_points[i] + fn._values_at_end_points[j] < fn(z):
                logging.info("pi(%s) + pi(%s) - pi(%s) < 0" % (x, y, x + y))
                return False
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
                epsz_set = [0, 1, -1]
                if limits_x[0] == limits_x[1] == limits_x[2] and limits_y[0] == limits_y[1] == limits_y[2]:
                    # continuous at x and y
                    for epsz in epsz_set:
                        if limits_x[1] + limits_y[1] - limits_z[eps + 1] < 0:
                            logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(epsz), y, print_sign(epsz), x + y, print_sign(epsz)))
                            return False
                else:
                    for eps in [1, -1]:
                        # (+ + +), (- - -). note that (0 0 0) has alreadly been checked.
                        if limits_x[eps + 1] + limits_y[eps + 1] - limits_z[eps + 1] < 0:
                            logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(eps), y, print_sign(eps), x + y, print_sign(eps)))
                            return False
                    for eps in [1, -1]:
                        for epsz in epsz_set:
                            # (+ - 0), (+ - +), (+ - -), (- + 0), (- + +), (- + -)
                            if limits_x[eps + 1] + limits_y[-eps + 1] - limits_z[epsz + 1] < 0:
                                logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(eps), y, print_sign(-eps), x + y, print_sign(epsz)))
                                return False
                        if limits_x[1] + limits_y[eps + 1] - limits_z[eps + 1] < 0:
                            # (0 + +), (0 - -)
                            logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(0), y, print_sign(eps), x + y, print_sign(eps)))
                            return False
                        if limits_x[eps + 1] + limits_y[1] - limits_z[eps + 1] < 0:
                            # (+ 0 +), (- 0 -)
                            logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(eps), y, print_sign(0), x + y, print_sign(eps)))
                            return False
    return True

def type2check_general(fn, continuity=True, limits=None):
    endpts = fn.end_points()
    endpts2 = endpts + [(endpts[i] + 1) for i in range(1, len(endpts))]
    for i in range(len(endpts)):
        for k2 in range(i + 1, i + len(endpts) - 1):
            x = endpts[i]
            z = endpts2[k2]
            y = z - x
            if k2 < len(endpts):
                k = k2
            else:
                k = k2 - len(endpts) + 1
            if fn._values_at_end_points[i] + fn(y) < fn._values_at_end_points[k]:
                logging.info("pi(%s) + pi(%s) - pi(%s) < 0" % (x, y, z))
                return False
            if not continuity:
                limits_x = limits[i]
                limits_z = limits[k]
                limits_y = fn.limits(y)
                # no trouble at 0- and 1+ since 0 < y < 1.
                if not (limits_y[0] == limits_y[1] == limits_y[2]):
                    # then y is a in endpts. this is done in type1check_general.
                    continue
                for epsx in [0, 1, -1]:
                    for epsz in [0, 1, -1]:
                        if limits_x[epsx + 1] + limits_y[1] < limits_z[epsz + 1]:
                            logging.info("pi(%s%s) + pi(%s) - pi(%s%s) < 0" % (x, print_sign(epsx), y, z, print_sign(epsz)))
                            return False
    return True
    
def subadditivity_check_general(fn, continuity=True, limits=None):
    """
    Check if fn is subadditive. Works for discontinuous functions as well.
    """
    endpts = fn.end_points()
    if continuity == None:
        continuity = fn.is_continuous_defined()
    if not continuity and limits == None:
        limits = [fn.limits(x) for x in endpts]
        limits[0][0] = limits[-1][0]
        limits[-1][2] = limits[0][2]
    if type1check_general(fn, continuity, limits) and type2check_general(fn, continuity, limits):
        logging.info("pi is subadditive.")
        return True
    else:
        logging.info("Thus pi is NOT subadditive.")
        return False

def symmetric_check_general(fn, f, continuity=True, limits=None):
    """
    Check if fn is symmetric. Works for discontinuous functions as well.
    """
    if fn(f) != 1:
        logging.info('pi(f) is not equal to 1.')
        logging.info('Thus pi is NOT symmetric.')
        return False
    endpts = fn.end_points()
    if continuity == None:
        continuity = fn.is_continuous_defined()
    if not continuity and limits == None:
        limits = [fn.limits(x) for x in endpts]
        limits[0][0] = limits[-1][0]
        limits[-1][2] = limits[0][2]
    for i in range(len(endpts)):
        x = endpts[i]
        if x == f:
            continue
        if x < f:
            y = f - x
        else:
            y = 1 + f - x
        if fn._values_at_end_points[i] + fn(y) != 1:
            logging.info('pi(%s) + pi(%s) is not equal to 1' % (x, y))
            logging.info('Thus pi is NOT symmetric.')
            return False
        if not continuity:
            limits_x = limits[i]
            limits_y = fn.limits(y)
            # no trouble at 0- and 1+ since 0 < y < 1.
            if limits_x[0] + limits_y[2] != 1:
                logging.info('pi(%s-) + pi(%s+) is not equal to 1' % (x, y))
                logging.info('Thus pi is NOT symmetric.')
                return False
            if limits_x[2] + limits_y[0] != 1:
                logging.info('pi(%s+) + pi(%s-) is not equal to 1' % (x, y))
                logging.info('Thus pi is NOT symmetric.')
                return False
    logging.info('pi is symmetric.')
    return True

def minimality_test_general(fn, f=None):
    """
    Check if fn is minimal with respect to f. Works for discontinuous functions as well.
    """
    for x in fn._values_at_end_points:
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
    endpts = fn.end_points()
    if fn.is_continuous_defined():
        continuity = True
        limits = None
    else:
        continuity = False
        limits = [fn.limits(x) for x in endpts]
        limits[0][0] = limits[-1][0]
        limits[-1][2] = limits[0][2]
        for x in limits:
            if not ((0 <= x[0] <=1) and (0 <= x[2] <=1)):
                logging.info('pi is not minimal because it does not stay in the range of [0, 1].')
                return False
    if subadditivity_check_general(fn, continuity, limits) and symmetric_check_general(fn, f, continuity, limits):
        logging.info('Thus pi is minimal.')
        return True
    logging.info('Thus pi is NOT minimal.')
    return False
