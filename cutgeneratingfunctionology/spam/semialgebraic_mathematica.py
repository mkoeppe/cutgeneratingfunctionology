r"""
Basic semialgebraic sets using the Mathematica interface
"""

from __future__ import division, print_function, absolute_import

from cutgeneratingfunctionology.spam.basic_semialgebraic import  BasicSemialgebraicSet_eq_lt_le_sets
from sage.interfaces.mathematica import mathematica
import operator

class BasicSemialgebraicSet_mathematica(BasicSemialgebraicSet_eq_lt_le_sets):
    """
    EXAMPLES::

        sage: from cutgeneratingfunctionology.spam.semialgebraic_mathematica import BasicSemialgebraicSet_mathematica
        sage: Q.<x,y,z> = QQ[]
        sage: bsa = BasicSemialgebraicSet_mathematica(QQ, poly_ring=Q)
        sage: bsa.add_polynomial_constraint(z-1, operator.eq)
        sage: bsa.add_polynomial_constraint(x^2 + y^2 - 4, operator.lt)
        sage: bsa.add_polynomial_constraint(y-1, operator.ge)
        sage: bsa
        BasicSemialgebraicSet_mathematica(eq=[z - 1], lt=[x^2 + y^2 - 4], le=[-y + 1])
        sage: bsa.variables_string()
        '{x, y, z}'
        sage: bsa.constraints_string()
        'z - 1 == 0 && x^2 + y^2 - 4 < 0 && -y + 1 <= 0'
        sage: bsa.find_point()  # optional - mathematica
        (0, 3/2, 1)
        sage: bsa.add_polynomial_constraint(y-x^2/3, operator.gt)
        sage: bsa.add_polynomial_constraint(x-3, operator.le)
        sage: bsa.add_polynomial_constraint(x+1, operator.ge)
        sage: bsa
        BasicSemialgebraicSet_mathematica(eq=[z - 1], lt=[1/3*x^2 - y, x^2 + y^2 - 4], le=[-y + 1, -x - 1, x - 3])
        sage: bsa.remove_redundant_constraints()  # optional - mathematica
        sage: bsa                                 # optional - mathematica
        BasicSemialgebraicSet_mathematica(eq=[z - 1], lt=[x^2 + y^2 - 4], le=[-y + 1, -x - 1])

    Test that the underscores have been replaced properly::

        sage: P.<x_0, x_1, x_2> = QQ[]
        sage: bsa = BasicSemialgebraicSet_mathematica(eq=[2*x_0*x_1 + 27*x_1 - 1/2], lt=[-27*x_0 - 2*x_1 - 1, 27/113*x_0^2 + x_1*x_2 + 1/2], le=[x_1^3 + x_0])
        sage: bsa.variables_string()
        '{x@0, x@1, x@2}'
        sage: bsa.constraints_string()
        '2*x@0*x@1 + 27*x@1 - 1/2 == 0 && 27/113*x@0^2 + x@1*x@2 + 1/2 < 0 && -27*x@0 - 2*x@1 - 1 < 0 && x@1^3 + x@0 <= 0'
        sage: bsa.find_point()  # optional - mathematica
         (-5/256, 64/3451, -512)
    """

    def _repr_(self):
        return 'BasicSemialgebraicSet_mathematica(eq={}, lt={}, le={})'.format(sorted(self._eq), sorted(self._lt), sorted(self._le))

    # def is_polynomial_constraint_valid(self, lhs, op):
    #     try:
    #         return super(BasicSemialgebraicSet_mathematica, self).is_polynomial_constraint_valid(lhs, op)
    #     except NotImplementedError:
    #         #TODO: Return true if FindInstance bsa cap lhs !op 0 is empty.
    #         #FindInstance is slow, and is_polynomial_constraint_valid and is_factor_known are used a lot. Is it worth doing the check?
    #         pass

    # def add_polynomial_constraint(self, lhs, op):
    #     pass # TODO: check is_polynomial_constraint_valid first. Worth doing?

    def find_point(self):
        """
        Use Mathematica's ``FindInstance``
        """
        pt_math = mathematica.FindInstance(self.constraints_string(), self.variables_string())  # For the first doctest, pt_math is of the form {{x -> 0, y -> 3/2, z -> 1}}.
        if len(pt_math) == 0:
            return None
        pt = []
        for i in range(self.ambient_dim()):
            try:
                pt_i = self.base_ring()(pt_math[1][i+1][2])
            except TypeError:
                pt_i = pt_math[1][i+1][2] #<class 'sage.interfaces.mathematica.MathematicaElement'>
            pt.append(pt_i)
        return tuple(pt)

    def variables_string(self):
        names = self.poly_ring().gens()
        varstr = str(names[0]).replace("_", "@")  # underscore is not allowed in variables names in mathematica. Replace x_1 by x@1, which will be treated as x[1] in mathematica.
        for v in names[1::]:
            varstr = varstr + ', ' + str(v).replace("_", "@")
        return '{' + varstr + '}'

    def constraints_string(self):
        constr = ''
        for l in self.eq_poly():
            constr += str(l) + ' == 0 && '
        for l in self.lt_poly():
            constr += str(l) + ' < 0 && '
        for l in self.le_poly():
            constr += str(l) + ' <= 0 && '
        return constr[:-4].replace("_", "@")

    def remove_redundant_constraints(self):
        r"""
        Use Mathematica's ``FindInstance`` to discard redundant inequalities or equations in the bsa description. Update self._eq, self._lt and self._le. The result is order dependent.
        """
        eqs = sorted(self.eq_poly())
        lts = sorted(self.lt_poly())
        les = sorted(self.le_poly())
        eq = []; lt = []; le = []
        for i in range(len(les)):
            bsa = BasicSemialgebraicSet_mathematica(eq=eqs, lt=lts, le=le+les[i+1::])
            bsa.add_polynomial_constraint(les[i], operator.gt)
            pt = bsa.find_point()
            if pt is not None:
                le.append(les[i])
        for i in range(len(lts)):
            bsa = BasicSemialgebraicSet_mathematica(eq=eqs, lt=lt+lts[i+1::], le=le)
            bsa.add_polynomial_constraint(lts[i], operator.ge)
            pt = bsa.find_point()
            if pt is not None:
                lt.append(lts[i])
        for i in range(len(eqs)):
            bsa = BasicSemialgebraicSet_mathematica(eq=eq+eqs[i+1::], lt=lt, le=le)
            bsa.add_polynomial_constraint(eqs[i], operator.lt)
            pt = bsa.find_point()
            if pt is not None:
                eq.append(eqs[i])
            else:
                bsa = BasicSemialgebraicSet_mathematica(eq=eq+eqs[i+1::], lt=lt, le=le)
                bsa.add_polynomial_constraint(eqs[i], operator.gt)
                pt = bsa.find_point()
                if pt is not None:
                    eq.append(eqs[i])
        self._eq = set(eq)
        self._lt = set(lt)
        self._le = set(le)
