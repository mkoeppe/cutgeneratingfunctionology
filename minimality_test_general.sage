def generate_type_1_vertices_general(fn, comparison, continuity=True, reduced=True):
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
                if reduced and limits_x[0] == limits_x[1] == limits_x[-1] and limits_y[0] == limits_y[1] == limits_y[-1]:
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
                        for zeps in [0, 1, -1]:
                            # (+ - 0), (+ - +), (+ - -), (- + 0), (- + +), (- + -)
                            if comparison(limits_x[eps] + limits_y[-eps] - limits_z[zeps], 0):
                                yield (x, y, x+y, eps, -eps, zeps)
                        if comparison(limits_x[0] + limits_y[eps] - limits_z[eps], 0):
                            # (0 + +), (0 - -)
                            yield (x, y, x+y, 0, eps, eps)
                        if comparison(limits_x[eps] + limits_y[0] - limits_z[eps], 0):
                            # (+ 0 +), (- 0 -)
                            yield (x, y, x+y, eps, 0, eps)

def generate_type_2_vertices_general(fn, comparison, continuity=True, reduced=True):
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
                if reduced:
                    for xeps in [0, 1, -1]:
                        for zeps in [0, 1, -1]:
                            if comparison(limits_x[xeps] + limits_y[0], limits_z[zeps]):
                                yield (x, y, z, xeps, 0, zeps)
                else:
                    for eps in [1, -1]:
                        # (+ + +), (- - -). note that (0 0 0) has alreadly been checked.
                        if comparison(limits_x[eps] + limits_y[eps] - limits_z[eps], 0):
                            yield (x, y, x+y, eps, eps, eps)
                    for eps in [1, -1]:
                        for zeps in [0, 1, -1]:
                            # (+ - 0), (+ - +), (+ - -), (- + 0), (- + +), (- + -)
                            if comparison(limits_x[eps] + limits_y[-eps] - limits_z[zeps], 0):
                                yield (x, y, x+y, eps, -eps, zeps)
                        if comparison(limits_x[0] + limits_y[eps] - limits_z[eps], 0):
                            # (0 + +), (0 - -)
                            yield (x, y, x+y, 0, eps, eps)
                        if comparison(limits_x[eps] + limits_y[0] - limits_z[eps], 0):
                            # (+ 0 +), (- 0 -)
                            yield (x, y, x+y, eps, 0, eps)
@cached_function
def generate_nonsubadditive_vertices_general(fn, continuity=True, reduced=True):
    """
    We are returning a set of 6-tuples (x, y, z, xeps, yeps, zeps),
    so that duplicates are removed, and so the result can be cached for later use.
    """
    return { (x, y, z, xeps, yeps, zeps)
             for (x, y, z, xeps, yeps, zeps) in itertools.chain( \
                generate_type_1_vertices_general(fn, operator.lt, continuity=continuity, reduced=reduced),\
                generate_type_2_vertices_general(fn, operator.lt, continuity=continuity, reduced=reduced)) }

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

def plot_2d_diagram_general(fn, show_function=True, continuity=None):
    """
    To show only a part of it, use
    `show(plot_2d_diagram(h), xmin=0.25, xmax=0.35, ymin=0.25, ymax=0.35)`

    EXAMPLES::
    sage: h = FastPiecewise([[closed_interval(0,1/4), FastLinearFunction(4, 0)], \
    ...          [open_interval(1/4, 1), FastLinearFunction(4/3, -1/3)], \
    ...          [singleton_interval(1), FastLinearFunction(0,0)]], merge=False)
    sage: plot_2d_diagram_general(h)

    sage: h = FastPiecewise([[closed_interval(0,1/4), FastLinearFunction(4, 0)], \
    ...           [open_interval(1/4,1/2), FastLinearFunction(3, -3/4)], \
    ...           [closed_interval(1/2, 3/4), FastLinearFunction(-2, 7/4)], \
    ...           [open_interval(3/4,1), FastLinearFunction(3, -2)], \
    ...           [singleton_interval(1), FastLinearFunction(0,0)]], merge=False)
    """
    if continuity == None:
        continuity = fn.is_continuous_defined()
    p = plot_2d_complex(fn)

    ## FIXME: Need to develope code for Face of discontinuous functions
    #y = var('y')
    #faces = generate_maximal_additive_faces(function)
    #kwds = { 'legend_label': "Additive face" }
    #for face in faces:
    #    p += face.plot(**kwds)
    #    if 'legend_label' in kwds:
    #        del kwds['legend_label']

    ### For non-subadditive functions, show the points where delta_pi
    ### is negative.


    def plot_nonsubadditive_fan(x, y, fan):
        r = 0.03
        if fan == (-1,-1,-1): #13
            p = disk((x,y), r, (pi, 3*pi/2), color='red', zorder=-5)
            p += line([(x,y), (x, y-r)], color='white', zorder=-4, thickness=3) # -3
            p += line([(x,y), (x-r, y)], color='white', zorder=-4, thickness=3) # -9
        if fan == (-1, 1,-1): #12
            p = disk((x,y), r, (3*pi/4, pi), color='red', zorder=-5)
            p += line([(x,y), (x-r, y)], color='white', zorder=-4, thickness=3) # -9
            p += line([(x,y), (x-r/sqrt(2), y+r/sqrt(2))], color='white', zorder=-4, thickness=3) #-10
        if fan == (-1, 1, 1): #11
            p = disk((x,y), r, (pi/2, 3*pi/4), color='red', zorder=-5)
            p += line([(x,y), (x, y+r)], color='white', zorder=-4, thickness=3) # -2
            p += line([(x,y), (x-r/sqrt(2), y+r/sqrt(2))], color='white', zorder=-4, thickness=3) #-10
        if fan == (-1, 1, 0): #10
            p = line([(x,y), (x-r/sqrt(2), y+r/sqrt(2))], color='red', zorder=-3, thickness=3)
            p += point([(x,y)], color = 'white', size=20, zorder=-2) # -1
        if fan == (-1, 0,-1): #9
            p = line([(x,y), (x-r, y)], color='red', zorder=-3, thickness=3)
            p += point([(x,y)], color = 'white', size=20, zorder=-2) # -1
        if fan == ( 1,-1,-1): #8
            p = disk((x,y), r, (3*pi/2, 7*pi/4), color='red', zorder=-5)
            p += line([(x,y), (x, y-r)], color='white', zorder=-4, thickness=3) # -3
            p += line([(x,y), (x+r/sqrt(2), y-r/sqrt(2))], color='white', zorder=-4, thickness=3) # -6
        if fan == ( 1,-1, 1): #7
            p = disk((x,y), r, (7*pi/4, 2*pi), color='red', zorder=-5)
            p += line([(x,y), (x+r, y)], color='white', zorder=-4, thickness=3) # -4
            p += line([(x,y), (x+r/sqrt(2), y-r/sqrt(2))], color='white', zorder=-4, thickness=3) # -6
        if fan == ( 1,-1, 0): #6
            p = line([(x,y), (x+r/sqrt(2), y-r/sqrt(2))], color='red', zorder=-3, thickness=3)
            p += point([(x,y)], color = 'white', size=20, zorder=-2) # -1
        if fan == ( 1, 1, 1): #5
            p = disk((x,y), r, (0, pi/2), color='red', zorder=-5)
            p += line([(x,y), (x+r, y)], color='white', zorder=-4, thickness=3) # -4
            p += line([(x,y), (x, y+r)], color='white', zorder=-4, thickness=3) # -2
        if fan == ( 1, 0, 1): #4
            p = line([(x,y), (x+r, y)], color='red', zorder=-3, thickness=3)
            p += point([(x,y)], color = 'white', size=20, zorder=-2) # -1
        if fan == ( 0,-1,-1): #3
            p = line([(x,y), (x, y-r)], color='red', zorder=-3, thickness=3)
            p += point([(x,y)], color = 'white', size=20, zorder=-2) # -1
        if fan == ( 0, 1, 1): #2
            p = line([(x,y), (x, y+r)], color='red', zorder=-3, thickness=3)
            p += point([(x,y)], color = 'white', size=20, zorder=-2) # -1
        if fan == ( 0, 0, 0): #1
            p = point([(x,y)], color = 'red', size=20, zorder=-1)
        return p

    nonsubadditive_vertices = generate_nonsubadditive_vertices_general(fn, continuity=continuity, reduced=False)
    if continuity:
        nonsubadditive_vertices = {(x,y) for (x, y, z, xeps, yeps, zeps) in nonsubadditive_vertices}
        p += point(list(nonsubadditive_vertices), \
                            color = "red", size = 50, legend_label="Subadditivity violated", zorder=-1)
        p += point([ (y,x) for (x,y) in nonsubadditive_vertices ], color = "red", size = 50, zorder=-1)
    elif nonsubadditive_vertices != {}:
        for (x, y, z, xeps, yeps, zeps) in nonsubadditive_vertices:
            p += plot_nonsubadditive_fan(x, y, (xeps, yeps, zeps))
            if x != y:
                p += plot_nonsubadditive_fan(y, x, (yeps, xeps, zeps))
        # add legend_label
        p += point([(0,0)], color = "red", size = 50, legend_label="Subadditivity violated", zorder=-10)
        p += point([(0,0)], color = "white", size = 50, zorder=-9)
    if show_function:
        x = var('x')
        #FIXME parametric_plot doesn't work for discontinuous functions.
        p += parametric_plot((lambda x: x, lambda x: 0.3 * float(fn(x)) + 1), \
                                                (x, 0, 1), color='blue', legend_label="Function pi")
        p += parametric_plot((lambda x: - 0.3 * float(fn(x)), lambda x: x), \
                                                (x, 0, 1), color='blue')
    return p

