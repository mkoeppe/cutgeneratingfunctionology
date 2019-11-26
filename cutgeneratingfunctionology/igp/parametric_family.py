"""
Parametric families of functions.
"""

from sage.misc import six
from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.misc.abstract_method import abstract_method
import logging

class ParametricFamily(six.with_metaclass(ClasscallMetaclass)):
    r"""
    A parametric family of functions.

    TESTS:

    Globally defined parametric families can be pickled::

        sage: from cutgeneratingfunctionology.igp import *
        sage: p_ll_strong_fractional = dumps(ll_strong_fractional)
        sage: explain_pickle(p_ll_strong_fractional)                     # not tested
        unpickle_global('cutgeneratingfunctionology.igp', 'll_strong_fractional')
        sage: loads(p_ll_strong_fractional) is ll_strong_fractional
        True

    Members of parametric families can be pickled and remember their construction
    but also their computed attributes, but there is no unique representation behavior::

        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
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
        instance = typecall(cls, *args, **options)
        assert isinstance( instance, cls )
        # if instance.__class__.__reduce__ == ParametricFamily.__reduce__:
        #     instance._reduction = (cls, args, options)
        instance._init_args = (cls, args, options)
        return instance

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

class ParametricFamily_introspect(ParametricFamily):
    r"""
    A parametric family of functions, defined in the traditional way by
    defining a function and then introspecting it for argument names and
    default arguments.
    """

    @staticmethod
    def __classcall__(cls, *args, **options):
        raise NotImplementedError
