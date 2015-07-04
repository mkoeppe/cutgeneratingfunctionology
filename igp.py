## Module

from sage.all import *

load("logging.sage")
load("real_number_field.sage")
load("functions.sage")
load("continuous_case.sage")
load("discontinuous_case.sage")
load("discrete_case.sage")
load("compendium_procedures.sage")
load("extreme_functions_in_literature.sage")
load("survey_examples.sage")
load("simple_extremality_test.sage")
load("parametric.sage")
load("kslope_ppl_mip.py")
load("quasi_periodic.sage")
load("extreme_functions_mlr_cpl3.sage")
load("kslope_pattern.sage")
load("2q_mip.sage")
load("kslope_mip.sage")
#load("old_walk_and_stability_interval_code.sage")

load("bug_examples.sage")

try:
    load("config.sage")
except IOError:
    pass

logging.info("Welcome to the infinite-group-relaxation-code. DON'T PANIC. See demo.sage for instructions.")
