### Uses global variable h

def format_interval(interval):
    l = latex(realset_from_interval(interval))
    return l.replace(r'-', r'\compactop-').replace(r'+', r'\compactop+')  # reduce spacing

def boldmath(s):
    return r'\text{\boldmath$' + s + '$}'

def format_interval_extra_fancy(nominal_I, minimal_I, number=None):
    if len(nominal_I) < 2:
        s = format_interval(nominal_I)
        if nominal_I[0] == number:
            s = boldmath(s)
        return s
    minimal_I = list(interval_to_endpoints(minimal_I))
    left = ""
    if nominal_I[0] < minimal_I[0]:
        left += r'\langle '
    elif h.limit(fractional(nominal_I[0]), 0) == h.limit(fractional(nominal_I[0]), 1):
        left += r'['
    else:
        left += r'('
    left += latex(nominal_I[0]).replace(r'-', r'\compactop-').replace(r'+', r'\compactop+')  # reduce spacing
    if nominal_I[0] == number:
        left = boldmath(left)
    right = latex(nominal_I[1]).replace(r'-', r'\compactop-').replace(r'+', r'\compactop+')  # reduce spacing
    if minimal_I[1] < nominal_I[1]:
        right += r'\rangle '
    elif h.limit(fractional(nominal_I[1]), 0) == h.limit(fractional(nominal_I[1]), -1):
        right += r']'
    else:
        right += r')'
    if nominal_I[1] == number:
        right = boldmath(right)
    s = left + r', ' + right
    ## if not h.special_intervals.is_disjoint_from(realset_from_interval(interval_mod_1(minimal_I))):
    ##     s = specialinterval(s)
    return s

def format_face_triple_extra_fancy(triple, vector=None):
    face = Face(triple)
    xyz = [vector[0], vector[1], vector[0] + vector[1]]
    return [ format_interval_extra_fancy(nominal_I, minimal_I, number)
             for nominal_I, minimal_I, number in zip(triple, face.minimal_triple, xyz) ]
