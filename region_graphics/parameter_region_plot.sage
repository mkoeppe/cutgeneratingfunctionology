def linesplot(K, fill=False, color="blue", linewidth=1, legend_label=None, xmin=-0.1, xmax=1.1, ymin=-0.1, ymax=1.1):
    x,y = var('x,y')
    leq, lin = simplify_eq_lt_poly_via_ppl(K.get_eq_factor(), K.get_lt_factor())
    if leq:
        print "WARNING: equation list is not empty!"
    if fill:
        g = region_plot([ lhs(x, y) < 0 for lhs in lin ], (x, xmin, xmax), (y, ymin, ymax), incol=color, plot_points=1000)
        return g
    g = Graphics()
    for l in lin:
        g += implicit_plot(l(x, y) == 0, (x, xmin, xmax), (y, ymin, ymax), color=color, linewidth=linewidth)
    g += line([(0,0),(0,1)], color = color, legend_label=legend_label, zorder=-10)
    return g
    
def regionplot_drlm_backward_3_slope(f_val=1/12-1/1000, b_val=2/12):
    K.<f, b> = SymbolicRealNumberField([f_val, b_val])
    h = drlm_backward_3_slope(f, b, field=K, conditioncheck=False)
    reg_cf = linesplot(K, color="orange", linewidth=5, legend_label="construction")

    minimality_test(h)
    reg_cfm = linesplot(K, color="green", linewidth=3, legend_label="min_test")

    extremality_test(h)
    reg_cfe = linesplot(K, color="blue", linewidth=1, legend_label="ext_test")
    reg = linesplot(K, fill=True, color="blue")

    K.<f, b> = SymbolicRealNumberField([f_val, b_val])
    h = drlm_backward_3_slope(f, b, field=K, conditioncheck=True)
    reg_ct  = linesplot(K, color="red", linewidth=4, legend_label="conditioncheck")

    p = point([K._values],color = "white", size = 20, zorder=5)

    g = reg+reg_cf+reg_ct+reg_cfm+reg_cfe+p
    return g
    
#sage: logging.disable(logging.INFO)
#sage: g1 = regionplot_drlm_backward_3_slope(f_val=1/12-1/30, b_val=2/12)
#sage: g1.save("drlm_backward_3_slope_1.pdf", title="drlm_backward_3_slope(1/12-1/30, 2/12)")
#sage: g2= regionplot_drlm_backward_3_slope(f_val=1/12+1/30, b_val=2/12)
#sage: g2.save("drlm_backward_3_slope_2.pdf", title="drlm_backward_3_slope(1/12+1/30, 2/12)")
