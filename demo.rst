=============================================
 Demonstration of cutgeneratingfunctionology
=============================================

We interact with Sage using commands that basically follow standard
Python syntax.  

See http://www.sagemath.org/doc/tutorial/ for information on how to use Sage.

Copy these commands into your Sage terminal session or Jupyter notebook.


Basic operations
================

First load the code::

    sage: import cutgeneratingfunctionology.igp as igp; from cutgeneratingfunctionology.igp import *

First we load a function and store it in variable ``h``.

We start with the easiest function, the GMIC::

    sage: h = gmic()


Plot the function::

    sage: plot_with_colored_slopes(h)


Test its extremality; this will create informative output and
plots::

    sage: extremality_test(h, show_plots=True)


The documentation string of each function reveals its optional
arguments, usage examples, and bibliographic information.  ::

    sage: gmic?


The docstring tells us that we can set the `f' value using an
optional argument::

    sage: h = gmic(1/5)


Of course, we know it will still be extreme; but let's test it to
see all the informative graphs::

    sage: extremality_test(h, show_plots=True)


Let's learn about a different function.  
It's from Dey--Richard--Li--Milller -- hence the prefix of the function::

    sage: drlm_backward_3_slope?


Let's change the parameters a little bit, so that they do NOT
satisfy the known sufficient conditions from the literature about
this class of functions::

    sage: h = drlm_backward_3_slope(f=1/12, bkpt=4/12)


Let's run the extremality test::

    sage: extremality_test(h, show_plots=True)

Indeed, it's not extreme.  We see a perturbation in magenta and the
two perturbed functions in blue and red, whose average is the
original function (black).

The extremality test stops when it has found one perturbation.  To
see more perturbations, we use the following: ::

    sage: extremality_test(h, show_plots=True, show_all_perturbations=True)


Here's the Gomory fractional cut.   ::

    sage: gomory_fractional?

::

    sage: h = gomory_fractional()

::

    sage: plot_with_colored_slopes(h)


It is not even minimal: ::

    sage: minimality_test(h, True)


Let's consider an interesting discontinuous function.
It was defined by Letchford and Lodi::

    sage: ll_strong_fractional?

The docstring suggests a few things to try with this function.
In particular it tells us that it fails to be minimal if 0 < f < 1/2.
Let's verify that::

    sage: h = ll_strong_fractional(1/3)

::

    sage: extremality_test(h, True)

There's many more functions to explore.  Try some of the following
to get more information.

Also consult the survey, which shows the graphs
of these functions for default parameters next to their names::

    sage: gj_2_slope?

::

    sage: gj_2_slope_repeat?

::

    sage: dg_2_step_mir?

::

    sage: kf_n_step_mir?

::

    sage: gj_forward_3_slope?

::

    sage: drlm_backward_3_slope?

::

    sage: dg_2_step_mir_limit?

::

    sage: drlm_2_slope_limit?

::

    sage: drlm_3_slope_limit?

::

    sage: bccz_counterexample?

::

    sage: bhk_irrational?

::

    sage: chen_4_slope?

::

    sage: rlm_dpl1_extreme_3a?

::

    sage: mlr_cpl3_g_3_slope?

Many more ``mlr_cpl3_...`` functions::

    sage: not_extreme_1?

::

    sage: drlm_not_extreme_2?

::

    sage: hildebrand_5_slope_28_1?

::

    sage: hildebrand_2_sided_discont_1_slope_1?

::

    sage: hildebrand_2_sided_discont_2_slope_1?

::

    sage: hildebrand_discont_3_slope_1?

::

    sage: gomory_fractional?


See the extreme function with the world-record number of different slopes::

    sage: extreme_function_with_world_record_number_of_slopes?

::

    sage: h = extreme_function_with_world_record_number_of_slopes()

::

    sage: plot_with_colored_slopes(h, thickness=4).show(figsize=30)


Use tab completion to see more example functions in the module
``extreme_functions``::

    sage: h = igp.extreme_functions.

Also see the live manual at
http://mkoeppe.github.io/cutgeneratingfunctionology/doc/html/extreme_functions.html


Procedures
==========

There are also "procedures", operations that can be applied to
extreme functions.  

First, the multiplicative homomorphism::

    sage: multiplicative_homomorphism?

::

    sage: h = multiplicative_homomorphism(gmic(f=4/5), 3)


Note, this function has several points where it takes value 1,
and hence several candidates for "f".  If we don't provide the f value
that we mean, it will warn and pick the first one from the left.
So let's provide the f value::

    sage: extremality_test(h, True, f=4/15)


A special case of the above::

    sage: automorphism?

::

    sage: h = automorphism(gmic(f=4/5))

::

    sage: extremality_test(h, True)


We can restrict to a finite group problem::

    sage: restrict_to_finite_group?

::

    sage: hr = restrict_to_finite_group(gmic(f=4/5))

::

    sage: extremality_test(hr, True)


For the finite group problems, automorphisms are more interesting!
 ::

    sage: ha = automorphism(hr, 2)

::

    sage: extremality_test(ha, True)


We can interpolate to get a function for the infinite group
problem::

    sage: hi = interpolate_to_infinite_group(ha)

::

    sage: extremality_test(hi, True)


The docstring has more interesting examples::

    sage: interpolate_to_infinite_group?


There's also this: ::

    sage: projected_sequential_merge?

The module ``procedures`` has a catalog of all procedures.  Use tab
completion to explore them::


    sage: igp.procedures.

Also see the live manual at
http://mkoeppe.github.io/cutgeneratingfunctionology/doc/html/procedures.html

    
Customizing the graphics
========================


Sometimes, for complicated functions, the graphics do not show
enough detail.   ::

    sage: h = bhk_irrational(delta=[1/200, 3/200])

::

    sage: extremality_test(h, True)


:-(

There are two ways to see more.

Approach 1: Use specific plotting functions to zoom in to some
areas of interest.

The 2d diagram, showing non-subadditive vertices and additive faces::

    sage: plot_2d_diagram(h).show(xmin=0.15, xmax=0.35, ymin=0.15, ymax=0.35)


The diagram of covered intervals::

    sage: plot_covered_intervals(h).show(xmin=0.15, xmax=0.35, ymin=0.22, ymax=0.55)


The completion diagram (to be explained in a forthcoming paper)::

    sage: plot_completion_diagram(h).show(xmin=0.28, xmax=0.52, ymin=0.25, ymax=0.35)


The perturbation diagram.  1 is the index of the perturbation shown::

    sage: plot_perturbation_diagram(h, 1).show(xmin=0.28, xmax=0.35, ymin=0.33, ymax=0.4)


Approach 2: Increase the plotting figure size (the default is 10)::

    sage: igp.show_plots_figsize = 40

::

    sage: h = bhk_irrational(delta=[1/200, 3/200])

::

    sage: extremality_test(h, True)


Defining new functions
======================

Of course, we can define functions from scratch::

    sage: h = piecewise_function_from_breakpoints_and_values([0, 1/5, 2/5, 4/5, 1], [0, 1/5, 3/5, 1, 0]);

::

    sage: extremality_test(h, True)


Here's another way::

    sage: slopes = [10/3,0,10/3,0,10/3,-10/3,0,-10/3,0,-10/3]

::

    sage: interval_lengths = [1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10]

::

    sage: h = piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes)

::

    sage: extremality_test(h, True, show_all_perturbations=True)
