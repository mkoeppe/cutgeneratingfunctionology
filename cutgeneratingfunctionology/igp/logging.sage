import logging

# this is a version of Python's StreamHandler which prints log
# messages to the stream *currently* pointed to by sys.stderr (not the
# one when StreamHandler is set up).  This is useful in a Sage notebook, 
# where every cell has its own set of streams.

class DynamicStdErrStreamHandler(logging.StreamHandler):
    """
    A handler class which writes logging records, appropriately formatted,
    to a stream. Note that this class does not close the stream, as
    sys.stdout or sys.stderr may be used.
    """

    def __init__(self):
        logging.StreamHandler.__init__(self, sys.stderr)
        self.parent_class = logging.StreamHandler           # save in object because name logging.StreamHandler is not available at exit

    def flush(self):
        try:
            self.stream = sys.stderr
        except (NameError, AttributeError):                                   # happens at exit in terminal
            pass
        self.parent_class.flush(self)

    def emit(self, record):
        try:
            self.stream = sys.stderr
        except (NameError, AttributeError):                                   # happens at exit in terminal
            pass
        self.parent_class.emit(self, record)

#logging.basicConfig(format='%(levelname)s: %(asctime)s %(message)s', level=logging.INFO)

# Set up stream handler for sage notebook.
# We always want to add the sage_hdlr to the root handler to enable logging.
# Other handlers could be added to the root logger such as file handlers
# before loading of the igp module.
# We don't want assume that igps sage stream logging should be dependent
# on no other handlers existing.
# Disable sage stream by calling logging.getLogger().removeHandler(sage_hdlr).
fmt = logging.Formatter('%(levelname)s: %(asctime)s %(message)s', None)
sage_hdlr = DynamicStdErrStreamHandler()
sage_hdlr.setFormatter(fmt)
logging.getLogger().addHandler(sage_hdlr) # Add the sage steam handler to the root logger. This is necessary.

# Note to future person looking at this code:
# Default logging levels are set by module loggers.
# .sage files use __name__ to be the python module they are called by not the full module path as would be in pure python.
# We name all module loggers by cutgeneratingfunctionology.igp.+file name in anticpation for a pure python file change.

def logger(func):
    def inner(*args, **kwargs): #1
        print("Arguments to %s were: %s, %s" % (func, args, kwargs))
        result = func(*args, **kwargs) #2
        print("Result is: %s" % (result))
        return result
        
    return inner

