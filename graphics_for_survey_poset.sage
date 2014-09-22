f = 6/7
# Base
v = dict()
v[0] = gomory_fractional(f)
v[1] = gmic(f)
assert extremality_test(v[1])
v[2] = multiplicative_homomorphism(gmic(f=fractional(2*f)),2)
#FastPiecewise([[(QQ(0), 5/14), FastLinearFunction(14/5, QQ(0))], [left_open_interval(5/14, 1/2), FastLinearFunction(-QQ(7), 7/2)], [left_open_interval(1/2, 6/7), FastLinearFunction(14/5, -7/5)], [left_open_interval(6/7, QQ(1)), FastLinearFunction(-QQ(7), QQ(7))]])
# obtained by sage_input(lift(drlm_not_extreme_1()))
assert extremality_test(v[2])
v[3] = ll_strong_fractional(f)
assert extremality_test(v[3])
v[12] = (v[1] + v[2]) / 2
v[13] = (v[1] + v[3]) / 2
v[23] = (v[2] + v[3]) / 2
v[123] = (v[1] + v[2] + v[3]) / 3

## Just the envelope of v[2] and v[1]
v['b'] = piecewise_function_from_breakpoints_and_values([0, 5/14, 6/7, 1], [0, 1, 1, 0])

## # from Wolfram mathworld
## w(x) = 3 * sqrt(1-(x/7)^2)
## l(x) = 1/2*(x+3) - 3/7*sqrt(10) * sqrt(4-(x+1)^2) + 6/7*sqrt(10)
## h(x) = 1/2*(3*(abs(x-1/2) + abs(x+1/2) + 6) - 11 * (abs(x-3/4) + abs(x+3/4)))
## r(x) = (6/7)*sqrt(10) + (3 - x)/2 - (3/7)*sqrt(10)*sqrt(4 - (x - 1)^2)
## b(x) = w(x) + (l(x) - w(x))*heaviside(x + 3) + (h(x) - l(x))*heaviside(x + 1) + (r(x) - h(x))*heaviside(x - 1) + (w(x) - r(x))*heaviside(x - 3)
## b2(x) = (1/2)*(3*sqrt(1 - (x/7)^2) + sqrt(1 - (abs(abs(x) - 2) - 1)^2) + abs(x/2) - ((3*sqrt(33) - 7)/112)*x^2 - 3)*((x + 4)/abs(x + 4) - (x - 4)/abs(x - 4)) - 3*sqrt(1 - (x/7)^2)
# {x, -7, 7},

v['b'] = california_ip()

def graph_ticks_keywords(function, y_ticks_for_breakpoints=False):
    #xticks = [0, 1] 
    xticks = v[2].end_points()
    xtick_formatter = [ "" for x in xticks ]
    #xtick_formatter = 'latex'  # would not show rationals as fractions
    ytick_formatter = None
    if y_ticks_for_breakpoints:
        yticks = xticks
        ytick_formatter = xtick_formatter
    else:
        #yticks = [0, 1]
        yticks = [1]
        ytick_formatter = [ "" for y in yticks ]
    return {'ticks': [xticks, yticks],
            'gridlines': True,
            'tick_formatter': [xtick_formatter, ytick_formatter],
            'axes': True}

igp.ticks_keywords = graph_ticks_keywords

for name, func in v.items():
    #g = plot_something(func)
    if name == 'b':
        uncovered_color
        g = plot(func, xmin=0, ymin=0, xmax=1, ymax=1.2, rgbcolor=(192/255, 54/255, 44/255), **graph_ticks_keywords(func))
        g += line([(0, func(0)), (5/14, func(5/14))], color='black')
        g += line([(6/7, func(6/7)), (1, func(1))], color='black')
        #### the following does not work; gives unwanted discontinuity markers
        #g += plot(func, xmin=0, ymin=0, xmax=5/14, ymax=1.2, color='black')
        #g += plot(func, xmin=6/7, ymin=0, xmax=1, ymax=1.2, color='black')
    else:
        g = plot_covered_intervals(func, xmin=0, ymin=0, xmax=1, ymax=1.2, color=uncovered_color)
    g.save(destdir + "graph-%s.pdf" % name, transparent=True, aspect_ratio=0.8, figsize=compendium_figsize)


