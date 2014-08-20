def print_sign(epsilon):
    if epsilon > 0:
        return "+"
    elif epsilon < 0:
        return "-"
    else:
        return ""

def range_zero_one_general(fn):
    """
    Check if fn has range in [0,1]. Works for discontinuous functions as well.
    """
    for x in fn.end_points():
        for epsilon in [0,1,-1]:
            if not (0 <= limit_periodic_one(fn, x, epsilon) <= 1):
                return False
    return True

def limit_periodic_one(fn, x, epsilon):
    """
    Similar as fn.limit(x, epsilon).
    Take care of left_limit at 0 and right_limit at 1.
    """
    if x == 0 and epsilon < 0:
        return fn.limit(1, epsilon)
    if x == 1 and epsilon > 0:
        return fn.limit(0, epsilon)
    return fn.limit(x, epsilon)
    
def delta_pi_limit(fn,x,y, epsx=0, epsy=0, epsz=0):
    """
    Compute the limit slack in subaddivity.
    """
    value_x = limit_periodic_one(fn, fractional(x), epsx)
    value_y = limit_periodic_one(fn, fractional(y), epsy)
    value_z = limit_periodic_one(fn, fractional(x+y), epsz)
    return value_x + value_y - value_z
  
def type1check_general(fn):
    endpts = fn.end_points()
    for x in endpts:
        for y in endpts:
            for eps in [0, 1, -1]:
            # (0 0 0), (+ + +), (- - -)
                if delta_pi_limit(fn, x, y, eps, eps, eps) < 0:
                    logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(eps), y, print_sign(eps), x + y, print_sign(eps)))
                    return False
            for eps in [1, -1]:
                for epsz in [0, 1, -1]:
                # (+ - 0), (+ - +), (+ - -), (- + 0), (- + +), (- + -)
                    if delta_pi_limit(fn, x, y, eps, -eps, epsz) < 0:
                        logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(eps), y, print_sign(-eps), x + y, print_sign(epsz)))
                        return False
                if delta_pi_limit(fn, x, y, 0, eps, eps) < 0:
                # (0 + +), (0 - -)
                        logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(0), y, print_sign(eps), x + y, print_sign(eps)))
                        return False
                if delta_pi_limit(fn, x, y, eps, 0, eps) < 0:
                # (+ 0 +), (- 0 -)
                        logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(eps), y, print_sign(0), x + y, print_sign(eps)))
                        return False
    return True

def type2check_general(fn):
    endpts = fn.end_points()
    endpts2 = endpts + [endpts[i] for i in range(1, len(endpts))]
    for x in endpts:
        for z in endpts2:
            y = z - x
            if (0 < y < 1) and not (y in endpts):
                for epsx in [0, 1, -1]:
                    for epsz in [0, 1, -1]:
                        if delta_pi_limit(fn, x, y, epsx, 0, epsz) < 0:
                            logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(epsx), y, print_sign(0), x + y, print_sign(epsz)))
                            return False
    return True
    
def subadditivity_check_general(fn):
    """
    Check if fn is subadditive. Works for discontinuous functions as well.
    """
    if type1check_general(fn) and type2check_general(fn):
        logging.info("pi is subadditive.")
        return True
    else:
        logging.info("Thus pi is NOT subadditive.")
        return False

from bisect import bisect_left
def symmetric_check_general(fn, f):
    """
    Check if fn is symmetric. Works for discontinuous functions as well.
    """
    if fn(f) != 1:
        logging.info('pi(f) is not equal to 1.')
        logging.info('Thus pi is NOT symmetric.')
        return False
    if fn.limit(0, 1) + fn.limit(f, -1) != 1:
        logging.info('pi(0+) + pi(f-) is not equal to 1.')
        logging.info('Thus pi is NOT symmetric.')
        return False
    if fn(1) != 0:
        logging.info('pi(1) is not equal to 0.')
        logging.info('Thus pi is NOT symmetric.')
        return False
    if fn.limit(f, 1) + fn.limit(1, -1) != 1:
        logging.info('pi(f+) + pi(1-) is not equal to 1.')
        logging.info('Thus pi is NOT symmetric.')
        return False    
    endpts = fn.end_points()
    k = bisect_left(endpts, f)
    for i in range(1, len(endpts) - 1):
        x = endpts[i]
        if i <= k:
            y = f - x
        else:
            y = 1 + f - x
        if i != k:
            for epsilon in [0, 1, -1]:
                if limit_periodic_one(fn, x, epsilon) + limit_periodic_one(fn, y, -epsilon) != 1:
                    logging.info('For x = %s, y = %s, pi(x%s) + pi(y%s) is not equal to 1' % (x, y, print_sign(epsilon), print_sign(epsilon)))
                    logging.info('Thus pi is NOT symmetric.')
                    return False
    logging.info('pi is symmetric.')
    return True

def minimality_test_general(fn, f=None):
    """
    Check if fn is minimal with respect to f. Works for discontinuous functions as well.
    """
    if not range_zero_one_general(fn):
        raise ValueError, "The given function does not stay in the range of [0, 1]."
    if f==None:
        f = find_f(fn)
    if fn(0) != 0:
        logging.info('pi is NOT minimal because pi(0) is not equal to 0.')
        return False
    logging.info('pi(0) = 0')
    if fn(1) != 0:
        logging.info('pi is NOT minimal because pi(1) is not equal to 0.')
        return False
    if  subadditivity_check_general(fn) and symmetric_check_general(fn, f):
        logging.info('Thus pi is minimal.')
        return True
    logging.info('Thus pi is NOT minimal.')
    return False
