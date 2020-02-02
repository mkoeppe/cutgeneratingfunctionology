r"""
Basic semialgebraic sets using the Mathematica interface
"""

from __future__ import division, print_function, absolute_import

from cutgeneratingfunctionology.spam.basic_semialgebraic import  BasicSemialgebraicSet_eq_lt_le_sets
from sage.interfaces.mathematica import mathematica
import operator
from sage.modules.free_module_element import vector
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.all import QQ, AA, RIF

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

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.semialgebraic_mathematica import BasicSemialgebraicSet_mathematica
            sage: P.<x> = QQ[]
            sage: BasicSemialgebraicSet_mathematica(QQ, poly_ring=P).find_point()   # optional - mathematica
            (0)
            sage: BasicSemialgebraicSet_mathematica(eq=[x-1]).find_point()   # optional - mathematica
            (1)
            sage: BasicSemialgebraicSet_mathematica(eq=[x^3-4]).find_point()   # optional - mathematica
            (1.587401051968200?)
            sage: Q.<y,z> = QQ[]
            sage: BasicSemialgebraicSet_mathematica(eq=[y], le=[z]).find_point()
            (0, -75)

       Was a bug example, whose second component -12/5 - I  is of <class 'sage.interfaces.mathematica.MathematicaElement'>. Now fixed by passing 'Reals' to FindInstance::

            sage: BasicSemialgebraicSet_mathematica(eq=[y]).find_point()   # optional - mathematica
            (0, -12/5)
        """
        try:
            # treat the cases where constraints_string is empty. Can also set self.constraints_string() to 'True', but returned point would have complex number component.
            if super(BasicSemialgebraicSet_mathematica, self).is_empty():
                return None
            if super(BasicSemialgebraicSet_mathematica, self).is_universe():
                return self.ambient_space().zero()
        except NotImplementedError:
            pass
        pt_math = mathematica.FindInstance(self.constraints_string(), self.variables_string(), 'Reals')  # For the first doctest, pt_math is of the form {{x -> 0, y -> 3/2, z -> 1}}.
        if len(pt_math) == 0:
            return None
        pt = vector(from_mathematica(pt_math[1][i+1][2]) for i in range(self.ambient_dim()))
        return self.ambient_space(field=pt.parent().base_ring())(pt)

    def is_empty(self):
        """
        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.semialgebraic_mathematica import BasicSemialgebraicSet_mathematica
            sage: P.<x,y> = QQ[]
            sage: BasicSemialgebraicSet_mathematica(QQ, poly_ring=P).is_empty()  # optional - mathematica
            False
            sage: bsa = BasicSemialgebraicSet_mathematica(eq=[x], le=[x-y])
            sage: bsa.is_empty()   # optional - mathematica
            False
            sage: bsa.add_polynomial_constraint(x-1, operator.eq)
            sage: bsa.is_empty()   # optional - mathematica
            True
        """
        try:
            return super(BasicSemialgebraicSet_mathematica, self).is_empty()
        except NotImplementedError:
            return (self.find_point() is None)

    def is_universe(self):
        """
        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.semialgebraic_mathematica import BasicSemialgebraicSet_mathematica
            sage: P.<x,y> = QQ[]
            sage: BasicSemialgebraicSet_mathematica(QQ, 2).is_universe()
            True
            sage: BasicSemialgebraicSet_mathematica(lt=[x-x]).is_universe()   # optional - mathematica
            False
            sage: BasicSemialgebraicSet_mathematica(eq=[y-y], lt=[x-x-1], le=[x-x]).is_universe()
            True
            sage: BasicSemialgebraicSet_mathematica(le=[-x^2-y^2]).is_universe()  # optional - mathematica
            True
        """
        try:
            return super(BasicSemialgebraicSet_mathematica, self).is_universe()
        except NotImplementedError:
            pt_math = mathematica.FindInstance('!(' + self.constraints_string() + ')', self.variables_string(), 'Reals')
            return len(pt_math) == 0

    def variables_string(self):
        names = self.poly_ring().gens()
        varstr = str(names[0]).replace("_", "@")  # underscore is not allowed in variables names in mathematica. Replace x_1 by x@1, which will be treated as x[1] in mathematica.
        for v in names[1::]:
            varstr = varstr + ', ' + str(v).replace("_", "@")
        if len(names) == 1:
            return varstr
        else:
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

def from_mathematica(a):
    try:
        return QQ(a.sage())
    except Exception:
        pass
    try:
        return AA(a.sage())
    except Exception:
            coefficients = mathematica.CoefficientList(mathematica.MinimalPolynomial(a, 'x'), 'x').sage()
            x = polygen(QQ)
            minpoly = x.parent()(coefficients)
            interval = mathematica.IsolatingInterval(a).sage()
            rif_interval = RIF(interval)
            return AA.polynomial_root(minpoly, rif_interval)
