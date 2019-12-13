"""
Parametric families of functions.
"""
from __future__ import print_function, absolute_import

from sage.misc.abstract_method import abstract_method
from sage.structure.unique_representation import UniqueRepresentation
from inspect import isclass
from sage.misc.sageinspect import sage_getargspec, sage_getvariablename
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.sage_object import SageObject
import logging
from .class_call import Classcall
from copy import copy
from collections import OrderedDict
from cutgeneratingfunctionology.spam.basic_semialgebraic import BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron

    ## @staticmethod
    ## def __classcall__(cls, *args, **options):
    ##     """
    ##     Construct a new object of this class and store the (normalized) init parameters.

    ##     This is like :meth:`~sage.structure.unique_representation.CachedRepresentation.__classcall__`,
    ##     but there is no caching.
    ##     """
    ##     if options.pop('conditioncheck', True):
    ##         cls.check_conditions(*args, **options)
    ##     if options.pop('compute_args_only', False):
    ##         return (cls, args, options)
    ##     return super(ParametricFamilyElement, cls).__classcall__(cls, *args, **options)

    # FIXME: centralize default argument handling:  Currently in derived classes, default arguments are provided both in __classcall__ and again in claimed_parameter_attribute.

    ## @classmethod
    ## def check_conditions(cls, *args, **kwargs):
    ##     try:
    ##         c = cls.claimed_parameter_attribute(*args, **kwargs)
    ##     except NotImplementedError:
    ##         return
    ##     if c == 'not_constructible':
    ##         raise ValueError("Bad parameters. Unable to construct the function.")
    ##     elif c == 'constructible':
    ##         logging.info("Conditions for extremality are NOT satisfied.")
    ##     else:   # FIXME: What about minimal, not extreme?
    ##         logging.info("Conditions for extremality are satisfied.")

class ParametricFamily_base(UniqueRepresentation, Parent):

    """
    Base class of parametric families of functions.

    EXAMPLES:

        sage: import logging; logging.disable(logging.INFO)
        sage: from cutgeneratingfunctionology.igp import *

    Parametric families contain their functions::

        sage: ll_strong_fractional() in ll_strong_fractional
        True

    """

    def __init__(self, default_values, names):
        from .fast_piecewise import FastPiecewise, PiecewiseLinearFunctionsSpace
        from sage.rings.rational_field import QQ
        Parent.__init__(self) #, facade=PiecewiseLinearFunctionsSpace(QQ, QQ))    ### how does it help??
        self._names = names
        self._default_values = OrderedDict(default_values)

    def parameter_space(self):
        """
        Return the parameter space of ``self`` as a basic semialgebraic set.

        EXAMPLE::

            sage: from cutgeneratingfunctionology.igp import *
            sage: ParametricFamily(chen_4_slope).parameter_space()
            BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(Constraint_System {}, names=[f, s_pos, s_neg, lam1, lam2])

        """
        return BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(names=self.names())

    @staticmethod
    def _format_function_call(fn_name, *v, **k):
        """
        Return a Python function call as a string.

        Keywords are sorted.
        """
        args = [ "%s" % a for a in v ] + [ "%s=%r" % (arg, val) for arg, val in sorted(k.items()) ]
        return "{}({})".format(fn_name, ", ".join(args))

    def _repr_(self):
        cls, args, kwds = self._reduction
        if self.__class__.__init__ == ParametricFamily_base.__init__:
            return self._format_function_call(cls.__name__, default_values=args[0], names=args[1])
        else:
            return self._format_function_call(cls.__name__, *args, **kwds)

    def names(self):
        """
        Return a list of the names of the parameters.
        """
        return self._names

    def default_point(self):
        """
        Return a list of the default values of the parameters corresponding to ``names``.
        """
        default_values = self.default_values()
        return [ default_values[name] for name in self.names() ]

    def default_values(self):
        """
        Return an ``OrderedDict`` of all arguments and their defaults
        that will be passed to ``constructor``.
        """
        return self._default_values

    def args_from_point(self, point, **options):
        """
        Return a dictionary of arguments suitable for ``_element_constructor_``,
        corresponding to ``point`` in the ``parameter_space``.
        """
        args = dict(zip(self.names(), point))
        args.update(options)
        return args

    def __call__(self, *args, **kwds):
        # Override to get rid of default x=0
        if not args:
            return self._element_constructor_(**kwds)
        #return super(ParametricFamily, self).__call__(*args, **kwds)
        ##  Above does not work because FastPiecewise is not an Element.
        return self._element_constructor_(*args, **kwds)

    def _args_with_defaults(self, *args, **options):
        names = self.names()
        assert len(args) <= len(names)
        values = copy(self.default_values())
        values.update(zip(names, args))
        values.update(options)
        return values

    @abstract_method
    def _construct_function(self, **values):
        """
        Construct the function.
        """

    def _element_constructor_(self, *args, **options):
        """
        The preferred calling method of this parametric family.
        """
        if len(args) == 1 and not options:
            try:
                if self in args[0]._parametric_families:
                    return args[0]
            except AttributeError:
                pass
        values = self._args_with_defaults(*args, **options)
        element = self._construct_function(**values)
        try:
            if not hasattr(element, "_parametric_families"):
                element._parametric_families = dict()
            element._parametric_families[self] = values
        except AttributeError:
            # Some legacy "constructors" actually just return an object of some built-in type.
            # Can't store an attribute in those.
            pass
        return element

    def _claimed_parameter_attribute(self, *args, **options):
        """
        Describe what the literature claims about the function.

        Return 'not_constructible', 'constructible', 'minimal', or 'extreme'.
        """
        try:
            h = self._element_constructor_(*args, **options)
            return h._claimed_parameter_attribute
        except Exception:
            # Function is non-contructible at this point.
            return 'not_constructible'
        return None

    def parameter_attribute(self, according_to='literature', *args, **options):
        """
        Determine whether the function obtained by ``_element_constructor_``, applied
        to the given parameters, is 'not_constructible', 'constructible', or 'extreme'.

        This information is by default given according to the published literature
        regarding this function (``according_to='literature'``), or according to an
        algorithmic test (``according_to='algorithm'``).

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO)
            sage: F_chen_4_slope = ParametricFamily(chen_4_slope)
            sage: F_chen_4_slope.parameter_attribute()
            'extreme'
            sage: F_chen_4_slope.parameter_attribute(according_to='literature', f=1/2, s_pos=5, s_neg=-5, lam1=1/5, lam2=1/5)
            'constructible'
            sage: F_chen_4_slope.parameter_attribute(according_to='literature', f=7/10, s_pos=2, s_neg=-4, lam1=1/100, lam2=49/100)
            'extreme'
            sage: F_chen_4_slope.parameter_attribute(according_to='best', f=7/10, s_pos=2, s_neg=-4, lam1=1/100, lam2=49/100)
            'constructible'

        """
        if according_to != 'algorithm':
            if 'condition_according_to_literature' in self.default_values():
                options['condition_according_to_literature'] = according_to == 'literature'
            return self._claimed_parameter_attribute(*args, **options)
        try:
            h = self._element_constructor_(*args, **options)
        except Exception:
            # Function is non-contructible at this point.
            return 'not_constructible'
        from . import extremality_test       # FIXME: This duplicates find_region_type
        if extremality_test(h):
            return 'extreme'
        else:
            return 'constructible'

    #def subfamily()    ...


class ParametricFamily(ParametricFamily_base):
    r"""A parametric family of functions.

    It is defined in the traditional way by defining a ``constructor`` and then
    introspecting it for argument names and default arguments.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *

    Some parametric families are already predefined::

        sage: isinstance(ll_strong_fractional, ParametricFamily)
        True
        sage: ll_strong_fractional
        ParametricFamily_ll_strong_fractional(default_values=(('f', 2/3), ('field', None), ('conditioncheck', True)), names=('f',))

    Making a parametric family is an idempotent operation::

        sage: F_ll_strong_fractional = ParametricFamily(ll_strong_fractional)
        sage: F_ll_strong_fractional is ll_strong_fractional
        True

    More examples of parametric families::

        sage: ParametricFamily(gj_2_slope)
        ParametricFamily_gj_2_slope(default_values=(('f', 3/5), ('lambda_1', 1/6), ('field', None), ('conditioncheck', True)), names=('f', 'lambda_1'))
        sage: ParametricFamily(bcdsp_arbitrary_slope, names=['f'])
        ParametricFamily_bcdsp_arbitrary_slope(ParametricFamily_bcdsp_arbitrary_slope(default_values=(('f', 1/2), ('k', 4), ('field', None), ('conditioncheck', True)), names=('f', 'k')), default_values=(('f', 1/2), ('k', 4), ('field', None), ('conditioncheck', True)), names=('f',))

        sage: from cutgeneratingfunctionology.igp import *
        sage: F_gj_2_slope = ParametricFamily(gj_2_slope)
        sage: h = F_gj_2_slope(); h
        <FastPiecewise ...>
        sage: h(3/5)
        1
        sage: h2 = F_gj_2_slope(3/4, 1/7)
        sage: h2(3/4)
        1

        sage: F_dg_2_step_mir = ParametricFamily(dg_2_step_mir)
        sage: h3 = F_dg_2_step_mir()
        sage: h3(4/5)
        1
        sage: h4 = F_dg_2_step_mir(5/6, 3/10)
        sage: h4(5/6)
        1

    Making a new parametric family by overriding defaults::

        sage: F_dg_2_step_mir_new_default = ParametricFamily(F_dg_2_step_mir, default_values={'f': 5/6, 'alpha': 3/10})
        sage: F_dg_2_step_mir_new_default() == h4
        True

    TESTS:

    Parametric families have unique representation behavior::

        sage: F_dg_2_step_mir_f_default = ParametricFamily(dg_2_step_mir, names=['alpha']); F_dg_2_step_mir_f_default
        ParametricFamily_dg_2_step_mir(ParametricFamily_dg_2_step_mir(default_values=(('f', 4/5), ('alpha', 3/10), ('field', None), ('conditioncheck', True)), names=('f', 'alpha')), default_values=(('f', 4/5), ('alpha', 3/10), ('field', None), ('conditioncheck', True)), names=('alpha',))
        sage: F_dg_2_step_mir_f_default is ParametricFamily(dg_2_step_mir, names=['alpha'])
        True
        sage: F_dg_2_step_mir_f_default is ParametricFamily(dg_2_step_mir, names=['alpha'], default_values={'f': 4/5})
        True
        sage: F_dg_2_step_mir_f_default is ParametricFamily(dg_2_step_mir, names=['alpha'], default_values={'f': 2/3})
        False

    Parametric families with globally defined constructors can be pickled::

        sage: p_F_ll_strong_fractional = dumps(F_ll_strong_fractional)
        sage: explain_pickle(p_F_ll_strong_fractional)                  # not tested
        pg_unreduce = unpickle_global('sage.structure.unique_representation', 'unreduce')
        pg_ParametricFamily = unpickle_global('cutgeneratingfunctionology.igp.parametric_family', 'ParametricFamily')
        pg_ll_strong_fractional = unpickle_global('cutgeneratingfunctionology.igp', 'll_strong_fractional')
        pg_make_rational = unpickle_global('sage.rings.rational', 'make_rational')
        pg_unreduce(pg_ParametricFamily, (pg_ll_strong_fractional, (('f', pg_make_rational('2/3'))), ('f',)), {})
        sage: loads(p_F_ll_strong_fractional) is F_ll_strong_fractional
        True

    Elements of parametric families can be pickled and remember their construction
    but also their computed attributes, but there is no unique representation behavior::

        sage: from cutgeneratingfunctionology.igp import *
        sage: import logging; logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = ParametricFamily(ll_strong_fractional)()
        sage: h._parametric_families
        {ParametricFamily_ll_strong_fractional(default_values=(('f', 2/3), ('field', None), ('conditioncheck', True)), names=('f',)): OrderedDict([('f', 2/3), ('field', None), ('conditioncheck', True)])}
        sage: h._difficult_computational_result = 42
        sage: p_h = dumps(h)
        sage: explain_pickle(p_h)              # not tested
        sage: h_copy = loads(p_h)
        sage: h_copy == h
        True
        sage: h_copy is h
        False
        sage: h_copy._parametric_families
        {ParametricFamily_ll_strong_fractional(default_values=(('f', 2/3), ('field', None), ('conditioncheck', True)), names=('f',)): OrderedDict([('f', 2/3), ('field', None), ('conditioncheck', True)])}
        sage: h_copy._difficult_computational_result
        42

    """

    @staticmethod
    def __classcall__(cls, constructor=None, default_values=None, names=None,
                      ignore_special_args=('cls', 'self'),
                      option_only_args=('conditioncheck', 'field', 'merge', 'condition_according_to_literature')):
        # For classes, sage_getargspec uses the argspec of __call__, which is not useful for us.
        if isinstance(constructor, ParametricFamily_base):
            if default_values is None and names is None:
                return constructor
            arg_default_dict = copy(constructor.default_values())
            args = constructor.names()
            cls = constructor._reduction[0]    # FIXME
        else:
            if constructor is None:
                c = cls._construct_function
            else:
                c = constructor
            args, varargs, keywords, defaults = sage_getargspec(c)
            #import pdb; pdb.set_trace()
            args = [ name for name in args if name not in ignore_special_args ]     # cls and self do not appear in the default list, remove them first
            if defaults:
                arg_default_dict = OrderedDict(zip(args, defaults))
            else:
                arg_default_dict = OrderedDict.fromkeys(args)
        if default_values is not None:
            default_values = dict(default_values)
            arg_default_dict.update(default_values)
        if names is None:
            names = [ name for name in args
                     if name not in option_only_args ]
        arg_default_dict = tuple((name, value) for name, value in arg_default_dict.items())
        return super(ParametricFamily, cls).__classcall__(
            cls, constructor, arg_default_dict, tuple(names))

    def __init__(self, constructor, default_values, names):
        from cutgeneratingfunctionology.igp import FastPiecewise
        super(ParametricFamily, self).__init__(default_values, names)
        if constructor is not None:
            self._construct_function = constructor

    def _repr_(self):
        cls, args, kwds = self._reduction
        def name(x):
            try:
                return x.__name__
            except AttributeError:
                return x
        if self.__class__.__init__ == ParametricFamily.__init__:
            constructor, default_values, names = args
            if constructor is None:
                return self._format_function_call(name(cls), default_values=default_values, names=names)
            else:
                return self._format_function_call(name(cls), name(constructor), default_values=default_values, names=names)
        else:
            return self._format_function_call(name(cls), *args, **kwds)

    def constructor(self):
        return self._construct_function
