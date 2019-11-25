"""
Parametric families of functions.
"""

from sage.misc import six
from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.misc.abstract_method import abstract_method
import logging

class ParametricFamily(six.with_metaclass(ClasscallMetaclass)):

    @staticmethod
    def __classcall__(cls, *args, **options):
        """
        Construct a new object of this class and store the (normalized) init parameters.

        This is like :meth:`~sage.structure.unique_representation.CachedRepresentation.__classcall__`,
        but there is no caching.
        """
        if options.pop('conditioncheck', True):
            cls.check_conditions(*args, **options)
        instance = typecall(cls, *args, **options)
        assert isinstance( instance, cls )
        # if instance.__class__.__reduce__ == ParametricFamily.__reduce__:
        #     instance._reduction = (cls, args, options)
        instance._init_args = (cls, args, options)
        return instance

    @classmethod
    def claimed_parameter_attributes(cls, *args, **kwargs):
        """
        Describe what the literature claims about the function.

        Return 'not_constructible', 'constructible', 'minimal', or 'extreme'.
        """
        raise NotImplementedError

    # FIXME: centralize default argument handling in something called by

    @classmethod
    def check_conditions(cls, *args, **kwargs):
        try:
            c = cls.claimed_parameter_attributes(*args, **kwargs)
        except NotImplementedError:
            return
        if c == 'not_constructible':
            raise ValueError("Bad parameters. Unable to construct the function.")
        elif c == 'constructible':
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    # What should reduce do?? We still want to save the function's computed attributes.

    ## @classmethod
    ## def from_constructor(cls, constructor):   ### or a separate class.

class ParametricFamily_introspect(ParametricFamily):

    @staticmethod
    def __classcall__(cls, *args, **options):


