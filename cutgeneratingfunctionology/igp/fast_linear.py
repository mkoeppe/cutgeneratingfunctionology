"""
Linear functions of 1 variable
"""
from __future__ import division, print_function, absolute_import

try:
    from sage.misc.repr import repr_lincomb
except ImportError:  # Sage < 9.2
    from sage.misc.misc import repr_lincomb

## FIXME: Its __name__ is "Fast..." but nobody so far has timed
## its performance against the other options. --Matthias
class FastLinearFunction :

    def __init__(self, slope, intercept):
        self._slope = slope
        self._intercept = intercept

    def __call__(self, x):
        if type(x) == float:
            # FIXME: There must be a better way.
            return float(self._slope) * x + float(self._intercept)
        else:
            return self._slope * x + self._intercept


    def __float__(self):
        return self

    def __add__(self, other):
        return FastLinearFunction(self._slope + other._slope,
                                  self._intercept + other._intercept)

    def __mul__(self, other):
        # scalar multiplication
        return FastLinearFunction(self._slope * other,
                                  self._intercept * other)


    def __neg__(self):
        return FastLinearFunction(-self._slope,
                                  -self._intercept)

    __rmul__ = __mul__

    def __eq__(self, other):
        if not isinstance(other, FastLinearFunction):
            return False
        return self._slope == other._slope and self._intercept == other._intercept

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        from cutgeneratingfunctionology.spam.parametric_real_field_element import is_parametric_element
        # Following the Sage convention of returning a pretty-printed
        # expression in __repr__ (rather than __str__).
        if not(is_parametric_element(self._slope) or is_parametric_element(self._intercept)):
            # repr_lincomb tests for 0, don't do this for parametric elements.
            try:
                return '<FastLinearFunction ' + repr_lincomb([('x', self._slope), (1, self._intercept)], strip_one = True) + '>'
            except TypeError:
                pass
        return '<FastLinearFunction (%s)*x + (%s)>' % (self._slope, self._intercept)

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.
        """
        return sib.name('FastLinearFunction')(sib(self._slope), sib(self._intercept))

    ## FIXME: To be continued.

fast_linear_function = FastLinearFunction

def linear_function_through_points(p, q):
    slope = (q[1] - p[1]) / (q[0] - p[0])
    intercept = p[1] - slope * p[0]
    return FastLinearFunction(slope, intercept) 

