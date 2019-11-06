def generate_additive_faces_discrete(fn):
    logging.info("Computing maximal additive faces...")
    faces = []
    additive_vertices = {(x,y) for (x, y, z, xeps, yeps, zeps) in generate_additive_vertices(fn) if x != 1 and y != 1}
    for (x, y) in additive_vertices:
        face = Face(([x], [y], [x+y]))
        faces.append(face)
        if x != y:
            face = Face(([y], [x], [x+y]))
            faces.append(face)
    logging.info("Computing maximal additive faces... done")
    return faces
