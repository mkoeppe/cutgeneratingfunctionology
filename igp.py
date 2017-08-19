## Module

from sage.all import *

igp_dir = os.path.dirname(__file__)
if igp_dir:
    igp_dir += "/"

load(igp_dir + "logging.sage")
load(igp_dir + "intervals.sage")
load(igp_dir + "real_number_field.sage")
load(igp_dir + "fast_linear.sage")
load(igp_dir + "functions.sage")
load(igp_dir + "continuous_case.sage")
load(igp_dir + "discontinuous_case.sage")
load(igp_dir + "discrete_case.sage")
load(igp_dir + "compendium_procedures.sage")
load(igp_dir + "extreme_functions_in_literature.sage")
load(igp_dir + "extreme_functions_sporadic.sage")
load(igp_dir + "survey_examples.sage")
load(igp_dir + "simple_extremality_test.sage")
load(igp_dir + "kslope_ppl_mip.py")
load(igp_dir + "vertex_enumeration.py")
load(igp_dir + "quasi_periodic.sage")
load(igp_dir + "extreme_functions_mlr_cpl3.sage")
load(igp_dir + "kslope_pattern.sage")
load(igp_dir + "2q_mip.sage")
load(igp_dir + "kslope_mip.sage")
#load("old_walk_and_stability_interval_code.sage")
load(igp_dir + "animation_2d_diagram.sage")

load(igp_dir + "bug_examples.sage")

import extreme_functions, procedures

try:
    load(igp_dir + "config.sage")
except IOError:
    pass

logging.info("Welcome to the infinite-group-relaxation-code. DON'T PANIC. See demo.sage for instructions.")
