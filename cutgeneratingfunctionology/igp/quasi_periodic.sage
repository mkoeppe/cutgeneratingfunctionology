from six.moves import range
class PiecewiseQuasiPeriodic(FastPiecewise):
    r"""
    Returns a piecewise quasi-periodic function from a list of (interval, function) pairs.
    Decomposes a piecewise quasi-periodic function into a periodic function and a linear function.

    It is assumed 'interval' from the first list of (interval, function) pairs is in the form (0,x1) for some x1 > 0.
    """

    def __init__(self, list_of_pairs, quasiperiodic_extension=True):
        r"""
        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO) # to disable output in automatic tests.
            sage: q = PiecewiseQuasiPeriodic([[(0,1/2), FastLinearFunction(3,0)],[(1/2,3/2), FastLinearFunction(-1,2)]])
            sage: q.period
            3/2
            sage: q.linear_term
            <FastLinearFunction 1/3*x>
            sage: q.periodic_term
            <FastPiecewise with 2 parts,
             (0, 1/2)   <FastLinearFunction 8/3*x>       values: [0, 4/3]
             (1/2, 3/2) <FastLinearFunction -4/3*x + 2>  values: [4/3, 0]>
        """

        FastPiecewise.__init__(self, list_of_pairs, var=None, periodic_extension=quasiperiodic_extension, merge=True)

        bkpt = self._end_points
        values = self._values_at_end_points

        self.period = bkpt[-1]-bkpt[0]

        # Linear function from quasi-periodic function
        self.linear_term = linear_function_through_points([bkpt[0],values[0]],[bkpt[-1],values[-1]])

        # Periodic function from quasi-periodic function
        values_of_periodic_function = []
        for i in range(len(bkpt)):
            values_of_periodic_function.append(values[i]-self.linear_term(bkpt[i]))

        self.periodic_term = piecewise_function_from_breakpoints_and_values(bkpt, values_of_periodic_function, field=None)
        self._quasiperiodic_extension = quasiperiodic_extension

    def __call__(self, x0):
        r"""
        Evaluats self at x0.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: logging.disable(logging.INFO) # to disable output in automatic tests.
            sage: q = PiecewiseQuasiPeriodic([[(0,1/2), FastLinearFunction(3,0)],[(1/2,3/2), FastLinearFunction(-1,2)]])
            sage: q(3/2)
            1/2
            sage: q(-1)
            1
        """

        bkpt = self._end_points
        values = self._values_at_end_points

        # can use "x0 - floor(interval / self.period) * self.period" to replace line 67-76
        if bkpt[0] <= x0 <= bkpt[-1]:
            return FastPiecewise.__call__(self,x0)
        else:

            if x0 < bkpt[0]:
                interval = bkpt[0] - x0
                shifted_x0 = x0 + ceil(interval/self.period)*self.period
            else:
                interval = x0 - bkpt[-1]
                shifted_x0 = x0 - floor(x0/self.period)*self.period

            value_of_periodic_term_at_shifted_x0 = self.periodic_term.__call__(shifted_x0)
            value_of_linear_term_at_shifted_x0 = self.linear_term(x0)

            value_at_shifted_x0 = value_of_periodic_term_at_shifted_x0 + value_of_linear_term_at_shifted_x0

        return value_at_shifted_x0

    def plot(self, *args, **kwds):
        r"""
        Returns the plot of self.

        Keyword arguments are passed onto the plot command for each piece
        of the function. E.g., the plot_points keyword affects each
        segment of the plot.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: q = PiecewiseQuasiPeriodic([[(0,1/2), FastLinearFunction(3,0)],[(1/2,3/2), FastLinearFunction(-1,2)]])
            sage: p = plot(q)

        The implementation of the plot method in Sage 5.11 piecewise.py
        is incompatible with the use of the xmin and xmax arguments.  Test that
        this has been fixed::

            sage: p = q.plot(xmin=0, xmax=3)
            sage: p = plot(q, xmin=0, xmax=3)
            sage: p = plot(q, 0, 3)
            sage: p = plot(q, 0, 3, color='red')

        Also the following plot syntax should be accepted::

            sage: p = plot(q, [0, 3])
        """

        from sage.plot.all import plot, Graphics

        g = Graphics()
        if not 'rgbcolor' in kwds:
            color = 'blue'
        else:
            color = kwds['rgbcolor']
        ### Code duplication with xmin/xmax code in plot.py.
        n = len(args)
        xmin = None
        xmax = None
        if n == 0:
            # if there are no extra args, try to get xmin,xmax from
            # keyword arguments
            xmin = kwds.pop('xmin', None)
            xmax = kwds.pop('xmax', None)
        elif n == 1:
            # if there is one extra arg, then it had better be a tuple
            xmin, xmax = args[0]
            args = []
        elif n == 2:
            # if there are two extra args, they should be xmin and xmax
            xmin = args[0]
            xmax = args[1]
            args = []
        point_kwds = dict()
        if 'alpha' in kwds:
            point_kwds['alpha'] = kwds['alpha']
        if 'legend_label' in kwds and self.is_discrete():
            point_kwds['legend_label'] = kwds['legend_label']
        # pop discontinuity markers.
        # Since quasiperiodic functions are continuous, discontinuity_markers are not used.
        discontinuity_markers = kwds.pop('discontinuity_markers', True)
        # if self._quasiperiodic_extension is True, plot the quasiperiodic extension of the function.

        if (xmin is None) or (not self._quasiperiodic_extension):
            xmin = 0
        if (xmax is None) or (not self._quasiperiodic_extension):
            xmax = xmin+self.period

        diff = self._values_at_end_points[-1]-self._values_at_end_points[0]
        number_of_repetition = floor(xmin/self.period)

        temp_list = []
        for (i, f) in self.list():

            a = i[0] + number_of_repetition*self.period
            b = i[1] + number_of_repetition*self.period

            temp_list.append([(a,b), linear_function_through_points([a, f(i[0]) + number_of_repetition*diff], [b, f(i[1]) + number_of_repetition*diff])])

        repetition = self._end_points[0] + number_of_repetition*self.period
        factor = 0
        while repetition <= xmax:

            for (i, f) in temp_list:

                a = i[0] + factor*self.period
                b = i[1] + factor*self.period

                if (xmin is not None) and (a < xmin):
                    a = xmin
                if (xmax is not None) and (b > xmax):
                    b = xmax
                if a < b:
                    g += plot(linear_function_through_points([a, f(i[0]) + factor*diff], [b, f(i[1]) + factor*diff]), *args, xmin=a, xmax=b, zorder=-1, **kwds)
                    delete_one_time_plot_kwds(kwds)

            factor += 1
            repetition += self.period
        return g

