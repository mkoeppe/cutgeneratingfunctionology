def generate_additive_faces_discrete(fn):
    """
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: points = [i/13 for i in range(14)]
        sage: values = [0, 6/7, 1/6, 5/7, 1/3, 4/7, 1/2, 3/7, 2/3, 2/7, 5/6, 1/7, 1, 0]
        sage: fn = discrete_function_from_points_and_values(points, values)
        sage: len(generate_additive_faces_discrete(fn))
        61
    """
    logging.info("Computing maximal additive faces...")
    faces = []
    additive_vertices = set()
    for (x, y, z, xeps, yeps, zeps) in generate_additive_vertices(fn):
        if x != 1 and y != 1:
            additive_vertices.add((x,y))
            additive_vertices.add((y,x))
    for (x, y) in additive_vertices:
        face = Face(([x], [y], [x+y]))
        faces.append(face)
    logging.info("Computing maximal additive faces... done")
    return faces
