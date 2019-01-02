
def plot_beams_of_one_face(face):
    g = Graphics()
    I, J, K = face.minimal_triple
    if len(I) == 2:
        g += polygon([(I[0], 1), (I[1], 1), (I[1], 0), (I[0], 0)], color='yellow', fill=True, alpha=0.35, zorder=-5)
    if len(J) == 2:
        g += polygon([(0, J[0]), (0, J[1]), (1, J[1]), (1, J[0])], color='yellow', fill=True, alpha=0.35, zorder=-5)
    if len(K) == 2:
        if coho_interval_contained_in_coho_interval(K, [0,1]):
            g += polygon([(K[0], 0), (K[1], 0), (0, K[1]), (0, K[0])], color='yellow', fill=True, alpha=0.35, zorder=-5)
        elif coho_interval_contained_in_coho_interval(K, [1,2]):
            g += polygon([(1, K[0]-1), (1, K[1]-1), (K[1] - 1, 1), (K[0] - 1, 1)], color='yellow', fill=True, alpha=0.35, zorder=-5)
        else:
            raise ValueError("Bad face: %s" % face)
    return g
