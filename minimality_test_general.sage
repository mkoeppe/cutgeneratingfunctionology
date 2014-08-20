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
  
def type1check_general(fn, endpts=None, continuity=None):
    if endpts == None:
        endpts = fn.end_points()
    if continuity == None:
        if fn.is_continuous_defined():
            continuity = True
        else:
            continuity = [(limit_periodic_one(fn, x, -1) == limit_periodic_one(fn, x, 0) == limit_periodic_one(fn, x, 1)) for x in endpts]
    for i in range(len(endpts)):
        for j in range(i,len(endpts)):
            x = endpts[i]
            y = endpts[j]
            if continuity == True or not ((x + y) in endpts):
                epsz_set = [0]
            else:
                # Want to check if fn is continuous at x+y. But finding k s.t. endpts[k] = x + y is O(log n). Worse than
                # checking directly (limit_periodic_one(fn, z, -1) == limit_periodic_one(fn, z, 0) == limit_periodic_one(fn, z, 1)).
                # It's not worth distinguishing the continuity at x+y.
                epsz_set = [0, 1, -1]
            if continuity == True or continuity[i] == True and continuity[j] == True:
                # continuous at x and y
                for epsz in epsz_set:
                    if delta_pi_limit(fn, x, y, 0, 0, epsz) < 0:
                        logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(epsz), y, print_sign(epsz), x + y, print_sign(epsz)))
                        return False
            else:
                # if continuity[i] == True and continuity[j] == False:
                # if continuity[i] == False and continuity[j] == True:
                # if continuity[i] == False and continuity[j] == False:
                # Too many cases... Just do them all in else:
                for eps in [0, 1, -1]:
                    # (0 0 0), (+ + +), (- - -)
                    if delta_pi_limit(fn, x, y, eps, eps, eps) < 0:
                        logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(eps), y, print_sign(eps), x + y, print_sign(eps)))
                        return False
                for eps in [1, -1]:
                    for epsz in epsz_set:
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

def type2check_general(fn, endpts=None, continuity=None):
    if endpts == None:
        endpts = fn.end_points()
    endpts2 = endpts + [(endpts[i] + 1) for i in range(1, len(endpts))]
    if continuity == None:
        if fn.is_continuous_defined():
            continuity = True
        else:
            continuity = [(limit_periodic_one(fn, x, -1) == limit_periodic_one(fn, x, 0) == limit_periodic_one(fn, x, 1)) for x in endpts]
    if not continuity == True:
        continuity2 = continuity + [continuity[i] for i in range(1, len(continuity))]
    for i in range(len(endpts)):
        for k in range(i + 1, i + len(endpts) -1):
            x = endpts[i]
            z = endpts2[k]
            y = z - x
            if not (y in endpts):
                # the case y in endpts is done in type1check_general.
                if continuity == True or continuity[i] == True and continuity2[k] == True:
                    # continuous at x and z (and y since y is not in endpts)
                    if delta_pi_limit(fn, x, y, 0, 0, 0) < 0:
                        logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(epsx), y, print_sign(0), x + y, print_sign(epsz)))
                        return False
                else:
                    for epsx in [0, 1, -1]:
                        for epsz in [0, 1, -1]:
                            if delta_pi_limit(fn, x, y, epsx, 0, epsz) < 0:
                                logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) < 0" % (x, print_sign(epsx), y, print_sign(0), x + y, print_sign(epsz)))
                                return False
    return True
    
def subadditivity_check_general(fn, endpts=None, continuity=None):
    """
    Check if fn is subadditive. Works for discontinuous functions as well.
    """
    if endpts == None:
        endpts = fn.end_points()
    if continuity == None:
        if fn.is_continuous_defined():
            continuity = True
        else:
            continuity = [(limit_periodic_one(fn, x, -1) == limit_periodic_one(fn, x, 0) == limit_periodic_one(fn, x, 1)) for x in endpts]
    if type1check_general(fn, endpts, continuity) and type2check_general(fn, endpts, continuity):
        logging.info("pi is subadditive.")
        return True
    else:
        logging.info("Thus pi is NOT subadditive.")
        return False

def symmetric_check_general(fn, f, endpts=None, continuity=None):
    """
    Check if fn is symmetric. Works for discontinuous functions as well.
    """
    if fn(f) != 1:
        logging.info('pi(f) is not equal to 1.')
        logging.info('Thus pi is NOT symmetric.')
        return False
    if fn(1) != 0:
        logging.info('pi(1) is not equal to 0.')
        logging.info('Thus pi is NOT symmetric.')
        return False
    if endpts == None:
        endpts = fn.end_points()
    if continuity == None:
        if fn.is_continuous_defined():
            continuity = True
        else:
            continuity = [(limit_periodic_one(fn, x, -1) == limit_periodic_one(fn, x, 0) == limit_periodic_one(fn, x, 1)) for x in endpts]
    for i in range(len(endpts)):
        x = endpts[i]
        if x <= f:
            y = f - x
        else:
            y = 1 + f - x
        if continuity == True or continuity[i] == True:
            eps_set = [0]
        else:
            eps_set = [0, 1, -1]
        for eps in eps_set:
            if limit_periodic_one(fn, x, eps) + limit_periodic_one(fn, y, -eps) != 1:
                logging.info('pi(%s%s) + pi(%s%s) is not equal to 1' % (x, print_sign(eps), y, print_sign(eps)))
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
    endpts = fn.end_points()
    # if the function fn is continuous of domain, set continuity = True.
    # Otherwise continuity is a list recording whether fn is continuous at each endpts.
    if fn.is_continuous_defined():
        continuity = True
    else:
        continuity = [(limit_periodic_one(fn, x, -1) == limit_periodic_one(fn, x, 0) == limit_periodic_one(fn, x, 1)) for x in endpts]
    if  subadditivity_check_general(fn, endpts, continuity) and symmetric_check_general(fn, f, endpts, continuity):
        logging.info('Thus pi is minimal.')
        return True
    logging.info('Thus pi is NOT minimal.')
    return False
