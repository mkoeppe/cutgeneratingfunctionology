## Demonstration of how to use infinite-group-relaxation-code.
##
## We interact with Sage using commands that basically follow standard
## Python syntax.  
##
## See http://www.sagemath.org/doc/tutorial/ for information on how to use Sage.
##
## Copy these commands into your Sage terminal session or notebook session.

## First load the code.
import igp; from igp import *

## First we load a function and store it in variable h.
##
## The functions are defined in the file
## extreme_functions_in_literature.sage
##
## Additional functions are defined in the file
## survey_examples.sage
##
## We start with the easiest function, the GMIC.
h = gmic()

## Plot the function.
plot_with_colored_slopes(h)

## Test its extremality; this will create informative output and
## plots.
extremality_test(h, show_plots=True)

## The documentation string of each function reveals its optional
## arguments, usage examples, and bibliographic information. 
gmic?

## The docstring tells us that we can set the `f' value using an
## optional argument.
h = gmic(1/5)

## Of course, we know it will still be extreme; but let's test it to
## see all the informative graphs.
extremality_test(h, show_plots=True)

## Let's learn about a different function.  
## It's from Dey--Richard--Li--Milller -- hence the prefix of the function.
drlm_backward_3_slope?

## Let's change the parameters a little bit, so that they do NOT
## satisfy the known sufficient conditions from the literature about
## this class of functions.
h = drlm_backward_3_slope(f=1/12, bkpt=4/12)

## Let's run the extremality test.
extremality_test(h, show_plots=True)
## Indeed, it's not extreme.  We see a perturbation in magenta and the
## two perturbed functions in blue and red, whose average is the
## original function (black).

## Here's the Gomory fractional cut.  
gomory_fractional?
h = gomory_fractional()
plot_with_colored_slopes(h)

## It is not even minimal:
minimality_test(h, True)

## Let's consider an interesting discontinuous function.
## It was defined by Letchford and Lodi.
ll_strong_fractional?

## The docstring suggests a few things to try with this function.
## In particular it tells us that it fails to be minimal if 0 < f < 1/2.
## Let's verify that.
h = ll_strong_fractional(1/3)
extremality_test(h, True)

## There's many more functions to explore.  Try some of the following
## to get more information.
## 
## Also consult the survey, which shows the graphs
## of these functions for default parameters next to their names.
gj_2_slope?
gj_2_slope_repeat?
dg_2_step_mir?
kf_n_step_mir?
gj_forward_3_slope?
drlm_backward_3_slope?
dg_2_step_mir_limit?
drlm_2_slope_limit?
drlm_3_slope_limit?
bccz_counterexample?
bhk_irrational?
chen_4_slope?
rlm_dpl1_extreme_3a?
not_extreme_1?
drlm_not_extreme_2?
hildebrand_5_slope_28_1?
hildebrand_2_sided_discont_1_slope_1?
hildebrand_2_sided_discont_2_slope_1?
hildebrand_discont_3_slope_1?
gomory_fractional?

## There are also "procedures" operations that can be applied to
## extreme functions.  They are defined in compendium_procedures.sage.

## First, the multiplicative homomorphism.
multiplicative_homomorphism?

h = multiplicative_homomorphism(gmic(f=4/5), 3)

## Note, this function has several points where it takes value 1,
## and hence several candidates for "f".  If we don't provide the f value
## that we mean, it will warn and pick the first one from the left.
## So let's provide the f value.
extremality_test(h, True, f=4/15)

## A special case of the above.
automorphism?

h = automorphism(gmic(f=4/5))
extremality_test(h, True)

## We can restrict to a finite group problem.
restrict_to_finite_group?

hr = restrict_to_finite_group(gmic(f=4/5))
extremality_test(hr, True)

## For the finite group problems, automorphisms are more interesting!

ha = automorphism(hr, 2)
extremality_test(ha, True)

## We can interpolate to get a function for the infinite group
## problem.
hi = interpolate_to_infinite_group(ha)
extremality_test(hi, True)

## The docstring has more interesting examples.
interpolate_to_infinite_group?

## There's also this:
projected_sequential_merge?

## Sometimes, for complicated functions, the graphics do not show
## enough detail.  
h = bhk_irrational(delta=[1/200, 3/200])
extremality_test(h, True)

## :-(

## There are two ways to see more.
##
## Approach 1: Use specific plotting functions to zoom in to some
## areas of interest.

## The 2d diagram, showing non-subadditive vertices and additive faces.
show(plot_2d_diagram(h), xmin=0.15, xmax=0.35, ymin=0.15, ymax=0.35)

## The diagram of covered intervals.
show(plot_covered_intervals(h), xmin=0.15, xmax=0.35, ymin=0.22, ymax=0.55)

## The completion diagram (to be explained in a forthcoming paper).
show(plot_completion_diagram(h), xmin=0.28, xmax=0.52, ymin=0.25, ymax=0.35)

## Approach 2: Increase the plotting figure size (the default is 10).
igp.show_plots_figsize = 40
h = bhk_irrational(delta=[1/200, 3/200])
extremality_test(h, True)

## Of course, we can define functions from scratch.

h = piecewise_function_from_breakpoints_and_values([0, 1/5, 2/5, 4/5, 1], [0, 1/5, 3/5, 1, 0]);
extremality_test(h, True)

## Here's another way.

slopes = [10/3,0,10/3,0,10/3,-10/3,0,-10/3,0,-10/3]
interval_lengths = [1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10]
h = piecewise_function_from_interval_lengths_and_slopes(interval_lengths, slopes)
extremality_test(h, True)

## See extreme_functions_in_literature.sage and survey_examples.sage for many more examples.

