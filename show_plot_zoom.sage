### Only works in terminal, not in worksheets.

# Set to True to have plots replace older plots with same tag and object.
# But FIXME: Preview on Mac OS X does not update reliably.
show_plots_update = False

global_plots_by_tag = dict()

def show_plot_zoomable(graphics, show_plots, tag, object=None, **show_kwds):
    """
    Display or save `graphics`.

    `show_plots` should be `False` (do nothing), 
    `True` (use `show` to display on screen),
    a string (file name format such as "FILENAME-%s.pdf", 
    where %s is replaced by `tag`.
    """
    if isinstance(show_plots, str):
        graphics.save(show_plots % tag, figsize=show_plots_figsize, **show_kwds)
    elif show_plots:
        fname_graphics_kwds = [None, None, None]
        if object is not None:
            if not hasattr(object, '_plots_by_tag'):
                object._plots_by_tag = dict()
            plots_by_tag = object._plots_by_tag
            if not show_plots_update or not tag in plots_by_tag:
                plots_by_tag[tag] = fname_graphics_kwds
            else:
                fname_graphics_kwds = plots_by_tag[tag]
        global_plots_by_tag[tag] = fname_graphics_kwds
        if fname_graphics_kwds[0] is None:
            fname_graphics_kwds[0] = tmp_filename('igp_%s_' % tag, ext='.png')
        fname = fname_graphics_kwds[0]
        fname_graphics_kwds[1] = graphics
        fname_graphics_kwds[2] = show_kwds
        graphics.show(figsize=show_plots_figsize, filename=fname, **show_kwds)

def zoom_plot(tag, object=None, **new_show_kwds):
    """
    After `minimality_test' or `extremality_test' has been run with show_plots=True,
    use zoom_plot to zoom in to a graph.

    EXAMPLE::
    sage: logging.disable(logging.INFO)
    sage: h = bhk_irrational_extreme_limit_to_rational_nonextreme()
    sage: extremality_test(h, True)
    False
    sage: zoom_plot('2d_diagram', xmin=0.1, xmax=0.4, ymin=0.1, ymax=0.4)
    sage: zoom_plot('completion', xmin=0.3, xmax=0.5, ymin=0.3, ymax=0.5)
    sage: zoom_plot('perturbation-1', xmin=0.25, xmax=0.4, ymin=0.25, ymax=0.5)
    """
    if object is not None:
        fname_graphics_kwds = object._plots_by_tag[tag]
    else:
        fname_graphics_kwds = global_plots_by_tag[tag]
    (fname, graphics, kwds) = fname_graphics_kwds
    kwds = copy(kwds)
    kwds.update(new_show_kwds)
    graphics.show(figsize=show_plots_figsize, filename=fname, **kwds)
