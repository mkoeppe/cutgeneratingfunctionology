r"""
Piecewise-defined Functions

Sage implements a very simple class of piecewise-defined functions.
Functions may be any type of symbolic expression. Infinite
intervals are not supported. The endpoints of each interval must
line up.

TODO:

- Implement max/min location and values,

- Need: parent object - ring of piecewise functions

- This class should derive from an element-type class, and should
  define ``_add_``, ``_mul_``, etc. That will automatically take care
  of left multiplication and proper coercion. The coercion mentioned
  below for scalar mult on right is bad, since it only allows ints and
  rationals. The right way is to use an element class and only define
  ``_mul_``, and have a parent, so anything gets coerced properly.

AUTHORS:

- David Joyner (2006-04): initial version

- David Joyner (2006-09): added __eq__, extend_by_zero_to, unextend,
  convolution, trapezoid, trapezoid_integral_approximation,
  riemann_sum, riemann_sum_integral_approximation, tangent_line fixed
  bugs in __mul__, __add__

- David Joyner (2007-03): adding Hann filter for FS, added general FS
  filter methods for computing and plotting, added options to plotting
  of FS (eg, specifying rgb values are now allowed). Fixed bug in
  documentation reported by Pablo De Napoli.

- David Joyner (2007-09): bug fixes due to behaviour of
  SymbolicArithmetic

- David Joyner (2008-04): fixed docstring bugs reported by J Morrow; added
  support for Laplace transform of functions with infinite support.

- David Joyner (2008-07): fixed a left multiplication bug reported by
  C. Boncelet (by defining __rmul__ = __mul__).

- Paul Butler (2009-01): added indefinite integration and default_variable

TESTS::

    sage: R.<x> = QQ[]
    sage: f = Piecewise([[(0,1),1*x^0]])
    doctest:...: DeprecationWarning: use lower-case piecewise instead
    See http://trac.sagemath.org/14801 for details.
    sage: 2*f
    doctest:...: DeprecationWarning: use lower-case piecewise instead
    See http://trac.sagemath.org/14801 for details.
    Piecewise defined function with 1 parts, [[(0, 1), 2]]
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 David Joyner <wdjoyner@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from six.moves import zip

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.sage_eval import sage_eval
from sage.rings.all import QQ, RR, Integer, Rational, infinity
from sage.calculus.functional import derivative
from sage.symbolic.expression import is_Expression
from sage.symbolic.assumptions import assume, forget

from sage.calculus.calculus import SR, maxima
from sage.calculus.all import var

from sage.structure.sage_object import SageObject

def Piecewise(list_of_pairs, var=None):
    """
    Deprecated spelling of :func:`sage.functions.piecewise`.

    Return a piecewise function from a list of (interval, function)
    pairs.
    
    ``list_of_pairs`` is a list of pairs (I, fcn), where
    fcn is a Sage function (such as a polynomial over RR, or functions
    using the lambda notation), and I is an interval such as I = (1,3).
    Two consecutive intervals must share a common endpoint.
    
    If the optional ``var`` is specified, then any symbolic expressions
    in the list will be converted to symbolic functions using 
    ``fcn.function(var)``.  (This says which variable is considered to
    be "piecewise".)

    We assume that these definitions are consistent (ie, no checking is
    done).
    
    EXAMPLES::
    
        sage: f1(x) = -1
        sage: f2(x) = 2
        sage: f = Piecewise([[(0,pi/2),f1],[(pi/2,pi),f2]])
        doctest:...: DeprecationWarning: use lower-case piecewise instead
        See http://trac.sagemath.org/14801 for details.
        sage: f(1)
        -1
        sage: f(3)
        2
        sage: f = Piecewise([[(0,1),x], [(1,2),x^2]], x); f
        Piecewise defined function with 2 parts, [[(0, 1), x |--> x], [(1, 2), x |--> x^2]]
        sage: f(0.9)
        0.900000000000000
        sage: f(1.1)
        1.21000000000000
    """
    from sage.misc.superseded import deprecation
    deprecation(14801, 'use lower-case piecewise instead')
    return PiecewisePolynomial(list_of_pairs, var=var)

class PiecewisePolynomial (SageObject):
    """
    Returns a piecewise function from a list of (interval, function)
    pairs.
    
    EXAMPLES::
    
        sage: f1(x) = -1
        sage: f2(x) = 2
        sage: f = Piecewise([[(0,pi/2),f1],[(pi/2,pi),f2]])
        doctest:...: DeprecationWarning: use lower-case piecewise instead
        See http://trac.sagemath.org/14801 for details.
        sage: f(1)
        -1
        sage: f(3)
        2
    """
    def __init__(self, list_of_pairs, var=None):
        r"""
        ``list_of_pairs`` is a list of pairs (I, fcn), where
        fcn is a Sage function (such as a polynomial over RR, or functions
        using the lambda notation), and I is an interval such as I = (1,3).
        Two consecutive intervals must share a common endpoint.
        
        If the optional ``var`` is specified, then any symbolic
        expressions in the list will be converted to symbolic
        functions using ``fcn.function(var)``.  (This says which
        variable is considered to be "piecewise".)

        We assume that these definitions are consistent (ie, no checking is
        done).

        EXAMPLES::

            sage: f1(x) = 1
            sage: f2(x) = 1 - x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.list()
            [[(0, 1), x |--> 1], [(1, 2), x |--> -x + 1]]
            sage: f.length()
            2
        """
        self._length = len(list_of_pairs)
        self._intervals = [x[0] for x in list_of_pairs]
        functions = [x[1] for x in list_of_pairs]
        if var is not None:
            for i in range(len(functions)):
                if is_Expression(functions[i]):
                    functions[i] = functions[i].function(var)
        self._functions = functions
        # We regenerate self._list in case self._functions was modified
        # above.  This also protects us in case somebody mutates a list
        # after they use it as an argument to piecewise().
        self._list = [[self._intervals[i], self._functions[i]] for i in range(self._length)]
 
    def list(self):
        """
        Returns the pieces of this function as a list of functions.
        
        EXAMPLES::

            sage: f1(x) = 1
            sage: f2(x) = 1 - x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.list()
            [[(0, 1), x |--> 1], [(1, 2), x |--> -x + 1]]
        """
        return self._list
 
    def length(self):
        """
        Returns the number of pieces of this function.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1 - x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.length()
            2
        """
        return self._length
 
    def __repr__(self):
        """
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1 - x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]]); f
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            Piecewise defined function with 2 parts, [[(0, 1), x |--> 1], [(1, 2), x |--> -x + 1]]
        """
        return 'Piecewise defined function with %s parts, %s'%(
            self.length(),self.list())
 
    def _latex_(self):
        r"""
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1 - x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: latex(f)
            \begin{cases}
            x \ {\mapsto}\ 1 &\text{on $(0, 1)$}\cr
            x \ {\mapsto}\ -x + 1 &\text{on $(1, 2)$}\cr
            \end{cases}
        
        ::
        
            sage: f(x) = sin(x*pi/2)
            sage: g(x) = 1-(x-1)^2
            sage: h(x) = -x
            sage: P = Piecewise([[(0,1), f], [(1,3),g], [(3,5), h]])
            sage: latex(P)
            \begin{cases}
            x \ {\mapsto}\ \sin\left(\frac{1}{2} \, \pi x\right) &\text{on $(0, 1)$}\cr
            x \ {\mapsto}\ -{\left(x - 1\right)}^{2} + 1 &\text{on $(1, 3)$}\cr
            x \ {\mapsto}\ -x &\text{on $(3, 5)$}\cr
            \end{cases}
        """
        from sage.misc.latex import latex
        tex = ['\\begin{cases}\n']
        for (left, right), f in self.list():
            tex.append('%s &\\text{on $(%s, %s)$}\\cr\n' % (latex(f), left, right))
        tex.append(r'\end{cases}')
        return ''.join(tex)

    def intervals(self):
        """
        A piecewise non-polynomial example.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.intervals()
            [(0, 1), (1, 2), (2, 3), (3, 10)]
        """
        return self._intervals
 
    def domain(self):
        """
        Returns the domain of the function.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.domain()
            (0, 10)
        """
        endpoints = sum(self.intervals(), ())
        return (min(endpoints), max(endpoints))

    def functions(self):
        """
        Returns the list of functions (the "pieces").
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.functions()
            [x |--> 1, x |--> -x + 1, x |--> e^x, x |--> sin(2*x)]
        """
        return self._functions
        
    def extend_by_zero_to(self,xmin=-1000,xmax=1000):
        """
        This function simply returns the piecewise defined function which
        is extended by 0 so it is defined on all of (xmin,xmax). This is
        needed to add two piecewise functions in a reasonable way.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1 - x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.extend_by_zero_to(-1, 3)
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            Piecewise defined function with 4 parts, [[(-1, 0), 0], [(0, 1), x |--> 1], [(1, 2), x |--> -x + 1], [(2, 3), 0]]
        """
        zero = QQ['x'](0)
        list_of_pairs = self.list()
        a, b = self.domain()
        if xmin < a:
            list_of_pairs = [[(xmin, a), zero]] + list_of_pairs
        if xmax > b:
            list_of_pairs = list_of_pairs + [[(b, xmax), zero]]
        return Piecewise(list_of_pairs)

    def unextend(self):
        """
        This removes any parts in the front or back of the function which
        is zero (the inverse to extend_by_zero_to).
        
        EXAMPLES::
        
            sage: R.<x> = QQ[]
            sage: f = Piecewise([[(-3,-1),1+2+x],[(-1,1),1-x^2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: e = f.extend_by_zero_to(-10,10); e
            Piecewise defined function with 4 parts, [[(-10, -3), 0], [(-3, -1), x + 3], [(-1, 1), -x^2 + 1], [(1, 10), 0]]
            sage: d = e.unextend(); d
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            Piecewise defined function with 2 parts, [[(-3, -1), x + 3], [(-1, 1), -x^2 + 1]]
            sage: d==f
            True
        """
        list_of_pairs = self.list()
        funcs = self.functions()
        if funcs[0] == 0:
            list_of_pairs = list_of_pairs[1:]
        if funcs[-1] == 0:
            list_of_pairs = list_of_pairs[:-1]
        return Piecewise(list_of_pairs)

    def base_ring(self):
        """
        Returns the base ring of the function pieces.   This
        is useful when this class is extended.

        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = x^2-5
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: base_ring(f)
            Symbolic Ring

        ::

            sage: R.<x> = QQ[]
            sage: f1 = x^0
            sage: f2 = 10*x - x^2
            sage: f3 = 3*x^4 - 156*x^3 + 3036*x^2 - 26208*x
            sage: f = Piecewise([[(0,3),f1],[(3,10),f2],[(10,20),f3]])
            sage: f.base_ring()
            Rational Field
        """
        return (self.functions()[0]).base_ring()

    def end_points(self):
        """
        Returns a list of all interval endpoints for this function.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = x^2-5
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.end_points()
            [0, 1, 2, 3]
        """
        intervals = self.intervals()
        return [ intervals[0][0] ] + [b for a,b in intervals]

    def __call__(self,x0):
        """
        Evaluates self at x0. Returns the average value of the jump if x0
        is an interior endpoint of one of the intervals of self and the
        usual value otherwise.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f(0.5)
            1
            sage: f(5/2)
            e^(5/2)
            sage: f(5/2).n()
            12.1824939607035
            sage: f(1)
            1/2
        """
        #x0 = QQ(x0) ## does not allow for evaluation at pi
        n = self.length()
        endpts = self.end_points()
        for i in range(1,n):
            if x0 == endpts[i]:
                return (self.functions()[i-1](x0) + self.functions()[i](x0))/2
        if x0 == endpts[0]:
            return self.functions()[0](x0)
        if x0 == endpts[n]:
            return self.functions()[n-1](x0)
        for i in range(n):
            if endpts[i] < x0 < endpts[i+1]:
                return self.functions()[i](x0)
        raise ValueError("Value not defined outside of domain.")

    def which_function(self,x0):
        """
        Returns the function piece used to evaluate self at x0.
        
        EXAMPLES::
        
            sage: f1(z) = z
            sage: f2(x) = 1-x
            sage: f3(y) = exp(y)
            sage: f4(t) = sin(2*t)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f.which_function(3/2)
            x |--> -x + 1
        """
        for (a,b), f in self.list():
            if a <= x0 <= b:
                return f
        raise ValueError("Function not defined outside of domain.")

    def default_variable(self):
        r"""
        Return the default variable. The default variable is defined as the
        first variable in the first piece that has a variable. If no pieces have
        a variable (each piece is a constant value), `x` is returned.

        The result is cached.

        AUTHOR: Paul Butler

        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 5*x
            sage: p = Piecewise([[(0,1),f1],[(1,4),f2]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: p.default_variable()
            x

            sage: f1 = 3*var('y')
            sage: p = Piecewise([[(0,1),4],[(1,4),f1]])
            sage: p.default_variable()
            y

        """
        try:
            return self.__default_variable
        except AttributeError:
            pass
        for _, fun in self._list:
            try:
                fun = SR(fun)
                if fun.variables():
                    v = fun.variables()[0]
                    self.__default_variable = v
                    return v
            except TypeError:
                # pass if fun is lambda function
                pass
        # default to x
        v = var('x')
        self.__default_value = v
        return v

    def plot(self, *args, **kwds):
        """
        Returns the plot of self.
        
        Keyword arguments are passed onto the plot command for each piece
        of the function. E.g., the plot_points keyword affects each
        segment of the plot.
        
        EXAMPLES::
        
            sage: f1(x) = 1
            sage: f2(x) = 1-x
            sage: f3(x) = exp(x)
            sage: f4(x) = sin(2*x)
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: P = f.plot(rgbcolor=(0.7,0.1,0), plot_points=40)
            sage: P
            Graphics object consisting of 4 graphics primitives

        Remember: to view this, type show(P) or P.save("path/myplot.png")
        and then open it in a graphics viewer such as GIMP.

        TESTS:

        We should not add each piece to the legend individually, since
        this creates duplicates (:trac:`12651`). This tests that only
        one of the graphics objects in the plot has a non-``None``
        ``legend_label``::

            sage: f1 = sin(x)
            sage: f2 = cos(x)
            sage: f = piecewise([[(-1,0), f1],[(0,1), f2]])
            sage: p = f.plot(legend_label='$f(x)$')
            sage: lines = [
            ....:   line
            ....:   for line in p._objects
            ....:   if line.options()['legend_label'] is not None ]
            sage: len(lines)
            1
        """
        from sage.plot.all import plot, Graphics

        g = Graphics()

        for i, ((a,b), f) in enumerate(self.list()):
            # If it's the first piece, pass all arguments. Otherwise,
            # filter out 'legend_label' so that we don't add each
            # piece to the legend separately (trac #12651).
            if i != 0 and 'legend_label' in kwds:
                del kwds['legend_label']

            g += plot(f, a, b, *args, **kwds)

        return g

    def _make_compatible(self, other):
        """
        Returns self and other extended to be defined on the same domain as
        well as a refinement of their intervals. This is used for adding
        and multiplying piecewise functions.
        
        EXAMPLES::
        
            sage: R.<x> = QQ[]
            sage: f1 = Piecewise([[(0, 2), x]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: f2 = Piecewise([[(1, 3), x^2]])
            sage: f1._make_compatible(f2)
            (Piecewise defined function with 2 parts, [[(0, 2), x], [(2, 3), 0]],
            Piecewise defined function with 2 parts, [[(0, 1), 0], [(1, 3), x^2]],
            [(0, 1), (1, 2), (2, 3)])
        """
        a1, b1 = self.domain()
        a2, b2 = other.domain()
        a = min(a1, a2)
        b = max(b1, b2)
        F = self.extend_by_zero_to(a,b)
        G = other.extend_by_zero_to(a,b)
        endpts = list(set(F.end_points()).union(set(G.end_points())))
        endpts.sort()
        return F, G, list(zip(endpts, endpts[1:]))

    def __add__(self,other):
        """
        Returns the piecewise defined function which is the sum of self and
        other. Does not require both domains be the same.
        
        EXAMPLES::
        
            sage: x = PolynomialRing(QQ,'x').gen()
            sage: f1 = x^0
            sage: f2 = 1-x
            sage: f3 = 2*x
            sage: f4 = 10-x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: g1 = x-2
            sage: g2 = x-5
            sage: g = Piecewise([[(0,5),g1],[(5,10),g2]])
            sage: h = f+g
            sage: h
            Piecewise defined function with 5 parts, [[(0, 1), x - 1], [(1, 2), -1], [(2, 3), 3*x - 2], [(3, 5), 8], [(5, 10), 5]]
        
        Note that in this case the functions must be defined using
        polynomial expressions *not* using the lambda notation.
        """
        F, G, intervals = self._make_compatible(other)
        fcn = []
        for a,b in intervals:
            fcn.append([(a,b), F.which_function(b)+G.which_function(b)])        
        return Piecewise(fcn)
        
    def __mul__(self,other):
        r"""
        Returns the piecewise defined function which is the product of one
        piecewise function (self) with another one (other).
        
        EXAMPLES::
        
            sage: x = PolynomialRing(QQ,'x').gen()
            sage: f1 = x^0
            sage: f2 = 1-x
            sage: f3 = 2*x
            sage: f4 = 10-x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: g1 = x-2
            sage: g2 = x-5
            sage: g = Piecewise([[(0,5),g1],[(5,10),g2]])
            sage: h = f*g
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: h
            Piecewise defined function with 5 parts, [[(0, 1), x - 2], [(1, 2), -x^2 + 3*x - 2], [(2, 3), 2*x^2 - 4*x], [(3, 5), -x^2 + 12*x - 20], [(5, 10), -x^2 + 15*x - 50]]
            sage: g*(11/2)
            Piecewise defined function with 2 parts, [[(0, 5), 11/2*x - 11], [(5, 10), 11/2*x - 55/2]]
            
        Note that in this method the functions must be defined using
        polynomial expressions *not* using the lambda notation.
        """
        ## needed for scalar multiplication
        if isinstance(other,Rational) or isinstance(other,Integer):
            return Piecewise([[(a,b), other*f] for (a,b),f in self.list()])
        else:
            F, G, intervals = self._make_compatible(other)
            fcn = []
            for a,b in intervals:
                fcn.append([(a,b),F.which_function(b)*G.which_function(b)])     
            return Piecewise(fcn)

    __rmul__ = __mul__

    def __eq__(self,other):
        """
        Implements Boolean == operator.

        EXAMPLES::
        
            sage: f1 = x^0
            sage: f2 = 1-x
            sage: f3 = 2*x
            sage: f4 = 10-x
            sage: f = Piecewise([[(0,1),f1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            doctest:...: DeprecationWarning: use lower-case piecewise instead
            See http://trac.sagemath.org/14801 for details.
            sage: g = Piecewise([[(0,1),1],[(1,2),f2],[(2,3),f3],[(3,10),f4]])
            sage: f==g
            True
        """
        return self.list()==other.list()
