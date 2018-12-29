from six.moves import range
#######################
# interface mathematica
#######################

def write_mathematica_constraints(eqs, ineqs, strict=True):
    r"""
    Write polynomial constraints in the mathematica format.
    Notice that the string ends with ' && '; in practice, often take condstr[:-4]

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: P.<x,y,z>=QQ[]
        sage: eqs = [z]
        sage: ineqs = [-x, x-1, -y, y-1]
        sage: write_mathematica_constraints(eqs, ineqs, strict=True)
        'z == 0 && -x < 0 && x - 1 < 0 && y - 1 < 0 && -y < 0 && '
        sage: write_mathematica_constraints(eqs, ineqs, strict=False)
        'z == 0 && -x <= 0 && x - 1 <= 0 && y - 1 <= 0 && -y <= 0 && '
    """
    condstr = ''
    for l in set(eqs):
        condstr += str(l) + ' == 0 && '
    for l in set(ineqs):
        if strict:
            condstr += str(l) + ' < 0 && '
        else:
            condstr += str(l) + ' <= 0 && '
    return condstr

def write_mathematica_variables(var_name):
    r"""
    Write the variables in the Mathematica format.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: var_name = ['x','y','z']
        sage: write_mathematica_variables(var_name)
        '{x, y, z}'
    """
    varstr = var_name[0]
    for v in var_name[1::]:
        varstr = varstr + ', ' + v
    return '{' + varstr + '}'

def find_instance_mathematica(condstr, var_name):
    r"""
    Call the Mathematica's ``FindInstance`` to get a point that satisfies the given conditions.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: condstr = 'z == 0 && -x < 0 && x - 1 < 0 && y - 1 < 0 && -y < 0'
        sage: var_name = ['x','y','z']
        sage: find_instance_mathematica(condstr, var_name)     # optional - mathematica
        (1/2, 1/2, 0)
    """
    varstr =  write_mathematica_variables(var_name)
    newvarstr = varstr.replace("_", "")
    newcondstr = condstr.replace("_", "")
    pt_math = mathematica.FindInstance(newcondstr, newvarstr)
    if len(pt_math) == 0:
        return None
    n = len(var_name)
    pt = []
    for i in range(n):
        try:
            pt_i = QQ(pt_math[1][i+1][2])
        except TypeError:
            pt_i = pt_math[1][i+1][2]
        pt.append(pt_i)
    return tuple(pt)
