"""
Parametric families of functions.
"""
from __future__ import print_function, absolute_import

from sage.misc.abstract_method import abstract_method
from sage.structure.unique_representation import UniqueRepresentation
from inspect import isclass
from sage.misc.sageinspect import sage_getargspec, sage_getvariablename
from sage.structure.parent import Parent
import logging
from .class_call import Classcall
from collections import OrderedDict
from cutgeneratingfunctionology.spam.basic_semialgebraic import BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron

class ParametricFamilyElement(Classcall):

    r"""
    Class of an element of a parametric family of functions.

    TESTS:

    Elements of parametric families can be pickled and remember their construction
    but also their computed attributes, but there is no unique representation behavior::

        sage: from cutgeneratingfunctionology.igp import *
        sage: import logging; logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = ll_strong_fractional()
        sage: h._init_args
        (<class 'cutgeneratingfunctionology.igp.ll_strong_fractional'>,
        (2/3,),
        {'field': None})
        sage: h._difficult_computational_result = 42
        sage: p_h = dumps(h)
        sage: explain_pickle(p_h)                                        # not tested
        pg_ll_strong_fractional = unpickle_global('cutgeneratingfunctionology.igp', 'll_strong_fractional')
        si1 = unpickle_newobj(pg_ll_strong_fractional, ())
        ...
        sage: h_copy = loads(p_h)
        sage: h_copy == h
        True
        sage: h_copy is h
        False
        sage: h_copy._init_args
        (<class 'cutgeneratingfunctionology.igp.ll_strong_fractional'>,
         (2/3,),
         {'field': None})
        sage: h_copy._difficult_computational_result
        42
    """

    @staticmethod
    def __classcall__(cls, *args, **options):
        """
        Construct a new object of this class and store the (normalized) init parameters.

        This is like :meth:`~sage.structure.unique_representation.CachedRepresentation.__classcall__`,
        but there is no caching.
        """
        if options.pop('conditioncheck', True):
            cls.check_conditions(*args, **options)
        if options.pop('compute_args_only', False):
            return (cls, args, options)
        return super(ParametricFamilyElement, cls).__classcall__(cls, *args, **options)

    @classmethod
    def claimed_parameter_attribute(cls, *args, **kwargs):
        """
        Describe what the literature claims about the function.

        Return 'not_constructible', 'constructible', 'minimal', or 'extreme'.
        """
        raise NotImplementedError

    # FIXME: centralize default argument handling:  Currently in derived classes, default arguments are provided both in __classcall__ and again in claimed_parameter_attribute.

    @classmethod
    def check_conditions(cls, *args, **kwargs):
        try:
            c = cls.claimed_parameter_attribute(*args, **kwargs)
        except NotImplementedError:
            return
        if c == 'not_constructible':
            raise ValueError("Bad parameters. Unable to construct the function.")
        elif c == 'constructible':
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")


class ParametricFamily(UniqueRepresentation, Parent):
    r"""A parametric family of functions.

    It is defined in the traditional way by defining a ``constructor`` and then
    introspecting it for argument names and default arguments.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: F_ll_strong_fractional = ParametricFamily(ll_strong_fractional); F_ll_strong_fractional
        ParametricFamily(ll_strong_fractional, names=['f'], default_values={'field': None, 'conditioncheck': True, 'f': 2/3})
        sage: ParametricFamily(gj_2_slope)
        ParametricFamily(gj_2_slope, names=['f', 'lambda_1'], default_values={'field': None, 'conditioncheck': True, 'lambda_1': 1/6, 'f': 3/5})
        sage: ParametricFamily(bcdsp_arbitrary_slope, names=['f'])
        ParametricFamily(bcdsp_arbitrary_slope, names=['f'], default_values={'field': None, 'k': 4, 'conditioncheck': True, 'f': 1/2})

    TESTS:

    Parametric families have unique representation behavior::

        sage: F_dg_2_step_mir_f_default = ParametricFamily(dg_2_step_mir, names=['alpha']); F_dg_2_step_mir_f_default
        ParametricFamily(dg_2_step_mir, names=['alpha'], default_values={'alpha': 3/10, 'conditioncheck': True, 'f': 4/5, 'field': None})
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
        pg_unreduce(pg_ParametricFamily, (pg_ll_strong_fractional, (('f', pg_make_rational('2/3')), ('field', None), ('conditioncheck', True)), ('f',), 'conditioncheck'), {})
        sage: loads(p_F_ll_strong_fractional) is F_ll_strong_fractional
        True

    """

    @staticmethod
    def __classcall__(cls, constructor, default_values=None, names=None,
                      ignore_special_args=('cls', 'self'),
                      ignore_args=('conditioncheck', 'field', 'merge', 'condition_according_to_literature')):
        # For classes, sage_getargspec uses the argspec of __call__, which is not useful for us.
        c = constructor
        if isclass(c):
            if hasattr(constructor, "__classcall__"):
                c = constructor.__classcall__
            elif hasattr(constructor, "__init__"):
                c = constructor.__init__
        args, varargs, keywords, defaults = sage_getargspec(c)
        #import pdb; pdb.set_trace()
        args = [ name for name in args if name not in ignore_special_args ]     # cls and self do not appear in the default list, remove them first
        if defaults:
            arg_default_dict = OrderedDict(zip(args, defaults))
        else:
            arg_default_dict = OrderedDict((arg, None) for arg in args)
        if default_values is not None:
            default_values = dict(default_values)
            arg_default_dict.update(default_values)
        if names is None:
            names = [ name for name in args
                     if name not in ignore_args ]
        return super(ParametricFamily, cls).__classcall__(
            cls, constructor, tuple(arg_default_dict.items()), tuple(names))

    def __init__(self, constructor, default_values, names):
        from cutgeneratingfunctionology.igp import FastPiecewise
        Parent.__init__(self, facade=(FastPiecewise,))
        self._constructor = constructor
        self._names = names
        self._default_values = OrderedDict(default_values)

    def _repr_(self):
        return "ParametricFamily({}, names={}, default_values={})".format(
            self.constructor().__name__, list(self.names()), dict(self.default_values()))

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

    def constructor(self):
        return self._constructor

    def parameter_space(self):
        """
        Return the parameter space of ``self`` as a basic semialgebraic set.

        EXAMPLE::

            sage: from cutgeneratingfunctionology.igp import *
            sage: ParametricFamily(chen_4_slope).parameter_space()
            BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(Constraint_System {}, names=[f, s_pos, s_neg, lam1, lam2])

        """
        return BasicSemialgebraicSet_polyhedral_ppl_NNC_Polyhedron(names=self.names())

    #def subfamily()    ...
