
def plot_no_legend(f, *args, **kwds):
    # really should rather use plot_kwds_hook everywhere in functions.sage
    plot_kwds_hook_no_legend(kwds)
    return sage.plot.plot.plot(f, *args, **kwds)

def plot_kwds_hook_no_legend(kwds):
    if 'legend_label' in kwds:
        del kwds['legend_label']
    if 'legend_loc' in kwds:
        del kwds['legend_loc']
    if 'legend_title' in kwds:
        del kwds['legend_title']

def dont_plot_rescaled_perturbation(*args, **kwds):
    return Graphics()
