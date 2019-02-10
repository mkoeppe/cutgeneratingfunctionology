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

if not logging.getLogger().handlers:
    fmt = logging.Formatter('%(levelname)s: %(asctime)s %(message)s', None)
    hdlr = DynamicStdErrStreamHandler()
    hdlr.setFormatter(fmt)
    logging.getLogger().addHandler(hdlr)
    logging.getLogger().setLevel(logging.INFO)


def logger(func):
    def inner(*args, **kwargs): #1
        print("Arguments to %s were: %s, %s" % (func, args, kwargs))
        result = func(*args, **kwargs) #2
        print("Result is: %s" % (result))
        return result
        
    return inner

