from __future__ import print_function, absolute_import

try:
    from sage.misc import six
except ImportError:
    import six

from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall

class Classcall(six.with_metaclass(ClasscallMetaclass)):

    @staticmethod
    def __classcall__(cls, *args, **options):
        instance = typecall(cls, *args, **options)
        assert isinstance( instance, cls )
        instance._init_args = (cls, args, options)
        return instance
