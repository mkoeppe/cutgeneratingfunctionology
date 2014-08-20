def print_sign(epsilon):
    if epsilon > 0:
        return "+"
    elif epsilon < 0:
        return "-"
    else:
        return ""

def nonnegative_check_general(fn):
    """
    Check if fn is nonnegative. Works for discontinuous functions as well.
    """
    if fn(0) < 0:
        logging.info('pi(0) is negative.')
        return False
    if fn.limit(0, 1) < 0:
        logging.info('pi(0+) is negative.')
        return False
    if fn(1) < 0:
        logging.info('pi(1) is negative.')
        return False
    if fn.limit(1, -1) < 0:
        logging.info('pi(1-) is negative.')
        return False
    endpts = fn.end_points()
    for i in range(1, len(endpts) - 1):
        x = endpts[i]
        for epsilon in [0,1,-1]:
            if fn.limit(x, epsilon) < 0:
                logging.info('pi(x%s) is negative.' % print_sign(epsilon))
                return False
    logging.info('pi is nonnegative.')
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
    
def delta_pi_limit(fn,x,y,epsilon):
    """
    Compute the right limit slack (if epsilon > 0) or left limit slack (if epsilon < 0) 
    or slack (if epsilon == 0) in subaddivity.
    """
    value_x = limit_periodic_one(fn, fractional(x), epsilon)
    value_y = limit_periodic_one(fn, fractional(y), epsilon)
    value_z = limit_periodic_one(fn, fractional(x+y), epsilon)
    return value_x + value_y - value_z
  
def type1check_general(fn):
    endpts = fn.end_points()
    for x in endpts:
        for y in endpts:
            for epsilon in [0, 1, -1]:
                if delta_pi_limit(fn, x, y, epsilon) < 0:
                    logging.info("For x = %s, y = %s, x+y = %s, Delta pi(x%s, y%s) = %s < 0" % (x, y, x + y, \
                                    print_sign(epsilon), print_sign(epsilon), delta_pi_limit(fn, x, y, epsilon)))
                    return False        
    return True

def type2check_general(fn):
    endpts = fn.end_points()
    endpts2 = []
    for i in range(len(endpts)-1):
        endpts2.append(endpts[i])
    for i in range(len(endpts)):
        endpts2.append(endpts[i]+1)
    for x in endpts:
        for z in endpts2:
            if z - x > 0:
                y = z - x
                for epsilon in [0, 1, -1]:
                    if delta_pi_limit(fn, x, y, epsilon) < 0:
                        logging.info("For x = %s, y = %s, x+y = %s, Delta pi(x%s, y%s) = %s < 0" % (x, y, x + y, \
                                        print_sign(epsilon), print_sign(epsilon), delta_pi_limit(fn, x, y, epsilon)))
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
    if f==None:
        f = find_f(fn)
    if fn(0) != 0:
        logging.info('pi is NOT minimal because pi(0) is not equal to 0.')
        return False
    logging.info('pi(0) = 0')
    if nonnegative_check_general(fn) and subadditivity_check_general(fn) and symmetric_check_general(fn, f):
        logging.info('Thus pi is minimal.')
        return True
    logging.info('Thus pi is NOT minimal.')
    return False
