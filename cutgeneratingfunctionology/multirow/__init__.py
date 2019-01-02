## Module

from __future__ import absolute_import
from sage.all import *
del SetPartitionsAk
del SetPartitionsBk
del SetPartitionsIk
del SetPartitionsPRk
del SetPartitionsPk
del SetPartitionsRk
del SetPartitionsSk
del SetPartitionsTk

from cutgeneratingfunctionology.igp import *

multirow_dir = os.path.dirname(__file__)
if multirow_dir:
    multirow_dir += "/"

# multirow
load(multirow_dir + "piecewise_functions.sage")
load(multirow_dir + "lifting_region.sage")
