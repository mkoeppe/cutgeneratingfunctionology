# Make sure current directory is in path.
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

class PiecewiseQuasiPeriodic(FastPiecewise):
    """
    Returns a piecewise quasi-periodic function from a list of (interval, function) pairs.
    Decomposes a piecewise quasi-periodic function into a periodic function and a linear function.

    It is assumed 'interval' from the first list of (interval, function) pairs is in the form (0,x1) for some x1 > 0.
    """

    def __init__(self, list_of_pairs):
        """
        EXAMPLE::

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

        FastPiecewise.__init__(self, list_of_pairs, var=None, periodic_extension=False, merge=True)

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

    def __call__(self, x0):
        """
        Evaluats self at x0.

        EXAMPLE::

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

