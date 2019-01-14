from cutgeneratingfunctionology.igp import *
h = gj_forward_3_slope()
g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), ticks=[h.end_points(),[]], tick_formatter=[["$0$","$a'$","$a$","$b$","$b'$","$f$","$1$"], []], thickness=2)
slope_formatter = ["$s^+$", "$s^-$", "$\\frac{1}{f}$", "$s^-$", "$s^+$", "$s^-$"]
for i in range(len(h.end_points())-1):
    x = (h.end_points()[i] + h.end_points()[i+1])/2 -1/50
    y = h(x) +1/10
    g += text(slope_formatter[i], (x, y), axis_coords=False, vertical_alignment='bottom',horizontal_alignment='left',color='black')
sphinx_plot(g)