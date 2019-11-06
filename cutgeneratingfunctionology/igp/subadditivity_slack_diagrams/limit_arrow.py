"""
A version of arrow2d.
"""

from sage.plot.arrow import Arrow

class LimitArrow(Arrow):

    def _render_on_subplot(self, subplot):
        r"""
        Render this arrow in a subplot.

        This version of the method uses a narrower arrow head,
        which is not customizable by parameters in the Sage class.
        """

        from sage.plot.misc import get_matplotlib_linestyle

        options = self.options()
        head = options.pop('head')
        if head == 0: style = '<|-'
        elif head == 1: style = '-|>'
        elif head == 2: style = '<|-|>'
        else: raise KeyError('head parameter must be one of 0 (start), 1 (end) or 2 (both).')
        style='fancy'
        width = float(options['width'])
        arrowshorten_end = float(options.get('arrowshorten', 0)) / 2.0
        arrowsize = float(options.get('arrowsize', 5))
        #head_width = arrowsize * 0.5
        head_width = arrowsize * 0.7
        tail_width = arrowsize * 0.7
        head_length = arrowsize * 2.0
        color = to_mpl_color(options['rgbcolor'])
        from matplotlib.patches import FancyArrowPatch
        p = FancyArrowPatch((self.xtail, self.ytail), (self.xhead, self.yhead),
                            lw=width,
                            arrowstyle='%s,head_width=%s,head_length=%s,tail_width=%s' % (style, head_width, head_length, tail_width),
                            shrinkA=arrowshorten_end, shrinkB=arrowshorten_end,
                            fc=color, ec=color,
                            linestyle=get_matplotlib_linestyle(options['linestyle'], return_type='long'))
        p.set_zorder(options['zorder'])
        p.set_label(options['legend_label'])

        subplot.add_patch(p)
        return p

from sage.misc.decorators import options, rename_keyword
from sage.plot.colors import to_mpl_color
@rename_keyword(color='rgbcolor')
@options(width=1, rgbcolor=(0,0,1), zorder=2, head=1, linestyle='solid', legend_label=None)
def limit_arrow(tailpoint=None, headpoint=None, path=None, **options):
    from sage.plot.all import Graphics
    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options))

    if headpoint is not None and tailpoint is not None:
        xtail, ytail = tailpoint
        xhead, yhead = headpoint
        g.add_primitive(LimitArrow(xtail, ytail, xhead, yhead, options=options))
    elif path is not None:
        g.add_primitive(CurveArrow(path, options=options))
    elif tailpoint is None and headpoint is None:
        return g
    else:
        raise TypeError('Arrow requires either both headpoint and tailpoint or a path parameter.')
    if options['legend_label']:
        g.legend(True)
        g._legend_colors = [options['legend_color']]
    return g
