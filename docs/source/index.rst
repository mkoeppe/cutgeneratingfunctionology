====================================================================================================
cutgeneratingfunctionology: Code for cut-generating functions in the Gomory-Johnson model and beyond
====================================================================================================

.. plot::

    from cutgeneratingfunctionology.igp import *
    h = equiv7_example_xyz_2()
    g = plot_completion_diagram(h)
    g += point((2.5, 1), color='white', figsize=(8,4))
    g.set_legend_options(title="Moves closure of equiv7_example_xyz_2()")
    sphinx_plot(g)

To use this module, you need to import it:: 

    import cutgeneratingfunctionology.igp as igp
    from cutgeneratingfunctionology.igp import *

Code for the 1-row Gomory--Johnson Infinite Group Relaxation
============================================================

.. toctree::
   :maxdepth: 1

   igp
   extreme_functions
   procedures

Code for Classical and General Dual-Feasible Functions
======================================================

.. toctree::
   :maxdepth: 1

   dff

Code for Multi-Row Models
=========================

.. toctree::
   :maxdepth: 1

   multirow

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
