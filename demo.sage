## Demonstrate what functions.sage can do.
##
## Load first one of the example functions to define `h'.

## Plot the function.
plot(h)

# The minimality test.
minimality_test(h, f)

# The 2d complex with all bells and whistles.
plot_2d_diagram(h)

# The covered (or connected-to-covered) intervals
generate_covered_intervals(h)

# A nice colorful plot of those
plot_covered_intervals(h)


# Or everything in one go:
extremality_test(h, show_plots=True)


try:
    lattice_plot(A, A0, t1, t2, 10)
except NameError:
    pass


