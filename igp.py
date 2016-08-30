## Module

from sage.all import *

load("logging.sage")
load("intervals.sage")
load("real_number_field.sage")
load("fast_linear.sage")
load("functions.sage")
load("continuous_case.sage")
load("discontinuous_case.sage")
load("discrete_case.sage")
load("compendium_procedures.sage")
load("extreme_functions_in_literature.sage")
load("extreme_functions_sporadic.sage")
load("survey_examples.sage")
load("simple_extremality_test.sage")
load("parametric.sage")
load("parametric_cpl.sage")
load("quasi_periodic.sage")
load("extreme_functions_mlr_cpl3.sage")
#load("old_walk_and_stability_interval_code.sage")

load("bug_examples.sage")

import extreme_functions, procedures

try:
    load("config.sage")
except IOError:
    pass

logging.info("Welcome to the infinite-group-relaxation-code. DON'T PANIC. See demo.sage for instructions.")
