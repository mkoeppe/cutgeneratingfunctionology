## Module

from __future__ import print_function, absolute_import

from sage.all import *
del SetPartitionsAk
del SetPartitionsBk
del SetPartitionsIk
del SetPartitionsPRk
del SetPartitionsPk
del SetPartitionsRk
del SetPartitionsSk
del SetPartitionsTk


igp_dir = os.path.dirname(__file__)
if igp_dir:
    igp_dir += "/"

import os

## We define a (simplified version of) Sage's load here
## so that exec uses the __future__ statement of this module,
## rather than that of the module where load is defined.
def igp_load(fpath):
    the_globals = globals()
    ext = os.path.splitext(fpath)[1].lower()
    if ext == '.py':
        with open(fpath) as f:
            code = compile(f.read(), fpath, 'exec')
            exec(code, the_globals)
    elif ext == '.sage':
        from sage.repl.attach import load_attach_mode
        from sage.repl.preparse import preparse_file_named, preparse_file
        # Preparse to a file to enable tracebacks with
        # code snippets. Use preparse_file_named to make
        # the file name appear in the traceback as well.
        # See Trac 11812.
        with open(preparse_file_named(fpath)) as f:
            code = compile(f.read(), preparse_file_named(fpath), 'exec')
            exec(code, the_globals)
    else:
        raise ValueError('unknown file extension %r for load or attach (supported extensions: .py, .pyx, .sage, .spyx, .f, .f90, .m)' % ext)

igp_load(igp_dir + "logging.sage")

from .intervals import *

igp_load(igp_dir + "real_number_field.sage")
igp_load(igp_dir + "functions.sage")
igp_load(igp_dir + "continuous_case.sage")
igp_load(igp_dir + "discontinuous_case.sage")
igp_load(igp_dir + "discrete_case.sage")
igp_load(igp_dir + "compendium_procedures.sage")
igp_load(igp_dir + "extreme_functions_in_literature.sage")
igp_load(igp_dir + "extreme_functions_sporadic.sage")
igp_load(igp_dir + "survey_examples.sage")
igp_load(igp_dir + "simple_extremality_test.sage")
igp_load(igp_dir + "parametric.sage")
igp_load(igp_dir + "parametric_cpl.sage")
igp_load(igp_dir + "kslope_ppl_mip.py")
igp_load(igp_dir + "vertex_enumeration.py")
igp_load(igp_dir + "quasi_periodic.sage")
igp_load(igp_dir + "extreme_functions_mlr_cpl3.sage")
igp_load(igp_dir + "kslope_pattern.sage")
igp_load(igp_dir + "2q_mip.sage")
igp_load(igp_dir + "kslope_mip.sage")
#igp_load("old_walk_and_stability_interval_code.sage")
igp_load(igp_dir + "animation_2d_diagram.sage")
igp_load(igp_dir + "crazy_perturbation.sage")
igp_load(igp_dir + "crazy_perturbation_examples.sage")
igp_load(igp_dir + "bug_examples.sage")
igp_load(igp_dir + "lifting_project.sage")
igp_load(igp_dir + "plot_options.sage")
igp_load(igp_dir + "faster_subadditivity_test.sage")
igp_load(igp_dir + "faster_subadditivity_test_discontinuous.sage")

from . import extreme_functions, procedures

try:
    igp_load(igp_dir + "config.sage")
except IOError:
    pass

## logging.info("Welcome to the infinite-group-relaxation-code. DON'T PANIC. See demo.sage for instructions.")
