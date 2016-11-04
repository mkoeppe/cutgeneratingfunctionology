nonzero_eps = { (-1,-1,-1), (-1, 1,-1), (-1, 1, 1), (-1, 1, 0), (-1, 0,-1), ( 1,-1,-1), \
                ( 1,-1, 1), ( 1,-1, 0), ( 1, 1, 1), ( 1, 0, 1), ( 0,-1,-1), ( 0, 1, 1) }
continuous_xy_eps = { (-1,-1,-1), (1, 1, 1) }
type2_reduced_eps = { (0,-1,-1), (0, 1, 1), (1,-1,-1), (1, 1, 1), (-1,-1,-1), (-1, 1, 1), \
                      (1,-1, 0), (-1, 1, 0) }
dic_eps_to_cone = { (-1,-1,-1): [(-1, 0), (0, -1)], \
                    (-1, 1,-1): [(-1, 1), (-1, 0)], \
                    (-1, 1, 1): [(0, 1), (-1, 1)], \
                    (-1, 1, 0): [(-1, 1)], \
                    (-1, 0,-1): [(-1, 0)], \
                    ( 1,-1,-1): [(0, -1), (1, -1)], \
                    ( 1,-1, 1): [(1, -1), (1, 0)], \
                    ( 1,-1, 0): [(1, -1)], \
                    ( 1, 1, 1): [(1, 0), (0, 1)], \
                    ( 1, 0, 1): [(1, 0)], \
                    ( 0,-1,-1): [(0, -1)], \
                    ( 0, 1, 1): [(0, 1)], \
                    ( 0, 0, 0): [] \
                  }

def generate_type_1_vertices_general_dff(fn, comparison):
    bkpt=fn.end_points()
    limits=[fn.limits(x) for x in bkpt]
    for i in range(len(bkpt)):
        for j in range(i,len(bkpt)):
            x=bkpt[i]
            y=bkpt[j]
            z=x+y
            if z<=1:
                if comparison(limits[i][0] + limits[j][0], fn(z)):
                    yield (x, y, x+y, 0, 0, 0)
                limits_x = limits[i]
                limits_y = limits[j]
                limits_z = fn.limits(z)
                for (xeps, yeps, zeps) in nonzero_eps:
                    if x==0 and xeps==-1:
                        continue
                    if y==0 and yeps==-1:
                        continue
                    if z==0 and zeps==-1:
                        continue
                    if x==1 and xeps==1:
                        continue
                    if y==1 and yeps==1:
                        continue
                    if z==1 and zeps==1:
                        continue
                    if comparison(limits_x[xeps] + limits_y[yeps] - limits_z[zeps], 0):
                        yield (x, y, x+y, xeps, yeps, zeps)

def generate_type_2_vertices_general_dff(fn, comparison):
    bkpt=fn.end_points()
    limits=[fn.limits(x) for x in bkpt]
    for i in range(len(bkpt)):
        for j in range(i,len(bkpt)):
            x=bkpt[i]
            z=bkpt[j]
            y=z-x
            if comparison(limits[i][0] + fn(y), limits[j][0]):
                yield (x, y, x+y, 0, 0, 0)
            limits_x = limits[i]
            limits_y = fn.limits(y)
            limits_z = limits[j]
            for (xeps, yeps, zeps) in nonzero_eps:
                if x==0 and xeps==-1:
                    continue
                if y==0 and yeps==-1:
                    continue
                if z==0 and zeps==-1:
                    continue
                if x==1 and xeps==1:
                    continue
                if y==1 and yeps==1:
                    continue
                if z==1 and zeps==1:
                    continue
                if comparison(limits_x[xeps] + limits_y[yeps] - limits_z[zeps], 0):
                    yield (x, y, x+y, xeps, yeps, zeps)

def generate_nonsymmetric_vertices_general_dff(fn):
    bkpt = fn.end_points()
    limits = fn.limits_at_end_points()
    for i in range(len(bkpt)):
        x = bkpt[i]
        y=1-x
        if limits[i][0] + fn(y) != 1:
            yield (x, y, 0, 0)
        limits_x = limits[i]
        limits_y = fn.limits(y)
        if x!=0 and y!=1:
            if limits_x[-1] + limits_y[1] != 1:
                yield (x, y, -1, 1)
        if x!=1 and y!=0:
            if limits_x[1] + limits_y[-1] != 1:
                yield (x, y, 1, -1)

def generate_nonsuperadditive_vertices_general(fn):
    return set(itertools.chain( \
                generate_type_1_vertices_general_dff(fn, operator.gt),\
                generate_type_2_vertices_general_dff(fn, operator.gt)))

def generate_additive_vertices_general_dff(fn):
    return set(itertools.chain( \
                generate_type_1_vertices_general_dff(fn, operator.eq),\
                generate_type_2_vertices_general_dff(fn, operator.eq)))

def symmetry_test_general_dff(fn):
    result = True
    for (x, y, xeps, yeps) in generate_nonsymmetric_vertices_general_dff(fn):
        logging.info('f(%s%s) + f(%s%s) is not equal to 1' % (x, print_sign(xeps), y, print_sign(yeps)))
        result = False
    if result:
        logging.info('f is symmetric.')
    else:
        logging.info('Thus f is not symmetric.')
    return result

def superadditivity_test_general(fn):
    """
    Check if `fn` is superadditive.
    """
    result = True
    for (x, y, z, xeps, yeps, zeps) in generate_nonsuperadditive_vertices_general(fn):
        logging.info("pi(%s%s) + pi(%s%s) - pi(%s%s) > 0" % (x, print_sign(xeps), y, print_sign(yeps), z, print_sign(zeps)))
        result = False
    if result:
        logging.info("f is superadditive.")
    else:
        logging.info("Thus f is not superadditive.")
    return result

def maximality_test_general_dff(fn):
    bkpt = fn.end_points()
    limits=[fn.limits(x) for x in bkpt]
    for i in range(len(limits)):
        for j in range(3):
            if (limits[i][j] < 0) or (limits[i][j] > 1):
                logging.info('f is not maximal because it does not stay in the range of [0, 1].')
                return False
    if fn(0) != 0:
        logging.info('f is NOT maximal because f(0) is not equal to 0.')
        return False
    logging.info('f(0) = 0')
    if superadditivity_test_general(fn) and symmetry_test_general_dff(fn):
        logging.info('Thus f is maximal.')
        is_maximal = True
    else:
        logging.info('Thus f is NOT maximal.')
        is_maximal = False
    return is_maximal


