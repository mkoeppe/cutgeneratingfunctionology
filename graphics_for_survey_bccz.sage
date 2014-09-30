# f = 1/2
# q = 4
#eta = 0
#f = 4/5
### Standard parameters of bccz_counterexample:
f = 2/3 
q = 4
eta = 1/1000
for n in [0, 1, 2, 3, Infinity]:
    if n == Infinity:
        h = bccz_counterexample(f=f, q=q, eta=eta)
    else:
        e = generate_example_e_for_psi_n(f=f, n=n, q=q, eta=eta)
        h = psi_n_in_bccz_counterexample_construction(f=f, e=e)
    xticks = [f, 1]
    xtick_formatter = ["$f$", "$1$"]
    if n >= 1:
        #xticks += [f/2-e[0]/2, f/2+e[0]/2]
        if n != Infinity: 
            xticks += h.end_points()
            xtick_formatter += [ "" for x in h.end_points() ]
            #xtick_formatter += ["$\\frac{1}{2} f-\\frac{1}{2} \\epsilon_1$", "$\\frac{1}{2} f+\\frac{1}{2} \\epsilon_1$"]
            #xtick_formatter += ["", ""]
    yticks = [1]
    ytick_formatter = ["$1$"]
    g = plot(h, ticks=[xticks, yticks],
             tick_formatter=[xtick_formatter, ytick_formatter],
             gridlines=True, color='black')
    if n != Infinity:
        for i, fn in h.list():
            if fn._slope < 0:
                try: 
                    k = 1 + e.index(interval_length(i))
                    g += text("$\\epsilon_{%s}$" % k, ((i[0] + i[1])/2, 0), axis_coords=True, vertical_alignment='top', horizontal_alignment='left', color='black')
                except ValueError:
                    pass
    if n == Infinity:
        fname = "bccz_counterexample_with_ticks.pdf"
    else:
        fname = "bccz_counterexample_with_ticks-%s.pdf" % n
    save(g, destdir + fname, figsize=3.5, aspect_ratio=0.9)
    
