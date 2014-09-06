########## Code for Continuous Case ###########

def generate_nonsymmetric_vertices_continuous(fn, f):
    bkpt = fn.end_points()
    for i in range(len(bkpt)):
        x = bkpt[i]
        if x == f:
            continue
        if x < f:
            y = f - x
        else:
            y = 1 + f - x
        if fn.values_at_end_points()[i] + fn(y) != 1:
            yield (x, y, 0, 0)

def generate_type_1_vertices_continuous(fn, comparison):
    """A generator...
    "...'general' refers to the fact that it outputs 6-tuples (x,xeps,y,yeps,z,zeps).
    FIXME: it currently does not take care of any discontinuities at all.
    """
    bkpt = fn.end_points()
    return ( (x, y, x+y, 0, 0, 0) for x in bkpt for y in bkpt if x <= y and comparison(delta_pi(fn,x,y), 0) ) # generator comprehension

def generate_type_2_vertices_continuous(fn, comparison):
    bkpt = fn.end_points()
    bkpt2 = bkpt[:-1] + [ x+1 for x in bkpt ]
    return ( (x, z-x, z, 0, 0, 0) for x in bkpt for z in bkpt2 if x < z < 1+x and comparison(delta_pi(fn, x, z-x), 0) ) # generator comprehension

