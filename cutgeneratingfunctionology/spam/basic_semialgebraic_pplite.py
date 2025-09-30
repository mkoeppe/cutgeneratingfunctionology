from cutgeneratingfunctionology.spam.basic_semialgebraic import BasicSemialgebraicSet_polyhedral

from pplite import Variable as pplite_Var, Constraint as pplite_Con, Linear_Expression as pplite_Lin_Expr, Affine_Expression as pplite_Aff_expr, NNC_Polyhedron as pplite_NNC_Polyhedron, PPliteGenerator, Polyhedron_Constraint_Rel, Polyhedron_Generator_Rel

poly_is_included_pplite = Polyhedron_Constraint_Rel.is_included()
point_is_included_pplite = Polyhedron_Generator_Rel.subsumes()

class BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron(BasicSemialgebraicSet_polyhedral):

    r"""
    A (possibly half-open) polyhedral basic semialgebraic set,
    represented by a PPLite ``NNC_Polyhedron``

    """

    def __init__(self, ambient_dim=None, polyhedron=None, base_ring=None, poly_ring=None, **options):
        r"""
        Initialize a basic semialgebraic set as the universe in
        ``ambient_dim``, or, if ``polyhedron`` (an ``NNC_Polyhedron``,
        which after that belongs to this object) is provided, as
        that.

        TEST::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P = BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron(2)
            sage: P.add_linear_constraint([0,1],0,operator.ge)
            sage: P.add_linear_constraint([1,0],0,operator.ge)
            sage: P.add_linear_constraint([2,3],-6,operator.lt)
            sage: P
            BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron([x1>=0, x0>=0, -2*x0-3*x1+6>0], names=[x0, x1])
            sage: sorted(P.eq_poly())
            []
            sage: sorted(P.lt_poly())
            [2*x0 + 3*x1 - 6]
            sage: sorted(P.le_poly())
            [-x1, -x0]
        """
        if ambient_dim is None and polyhedron is not None:
            ambient_dim = polyhedron.space_dimension()
        if base_ring is None and poly_ring is None:
            base_ring = QQ
        poly_ring, base_ring, ambient_dim, names = self._poly_ring_from_options(
            ambient_dim=ambient_dim, base_ring=base_ring, poly_ring=poly_ring, **options)
        if base_ring is not QQ:
            raise ValueError("only base_ring=QQ is supported")
        super(BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron, self).__init__(poly_ring=poly_ring)
        if polyhedron is None:
            self._polyhedron = pplite_NNC_Polyhedron(dim_type=int(ambient_dim), spec_elem='universe', topology="nnc")  # To work with pplite, ambient_dim is required to be type int
        else:
            self._polyhedron = polyhedron

    @staticmethod
    def _pplite_constraint(lhs, cst, op):
        r"""
        Make a PPL ``Constraint`` ``lhs`` * x + cst ``op`` 0,
        where ``lhs`` is be a vector of length ambient_dim.
        """
        lcd = lcm(lcm(x.denominator() for x in lhs), cst.denominator())
        lin_expr = sum([int((lcd*lhs[i]))*pplite_Var(i) for i in range(len(lhs))])
        aff_expr = pplite_Aff_expr(lin_expr, int(lcd*cst))
        #linexpr = pplite_Lin_Expr(lhs * lcd, cst * lcd)
        if op == operator.lt:
            return (aff_expr < 0)
        elif op == operator.gt:
            return (aff_expr > 0)
        elif op == operator.eq:
            return (aff_expr == 0)
        elif op == operator.le:
            return (aff_expr <= 0)
        elif op == operator.ge:
            return (aff_expr >= 0)
        else:
            raise ValueError("{} is not a supported operator".format(op))

    def __copy__(self):
        r"""
        Make a copy of ``self``.

        TESTS:

        Test that it is actually making a copy of the (mutable!) NNC_Polyhedron::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P = BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron(2)
            sage: P._polyhedron is copy(P)._polyhedron
            False
        """
        return self.__class__(polyhedron=pplite_NNC_Polyhedron(nnc_poly=self._polyhedron), poly_ring=self.poly_ring())

    def _repr_(self):
        constraints = self._polyhedron.constraints()
        names = list(self.poly_ring().gens())
        return 'BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron({}, names={})'.format(
            constraints, names)

    def closure(self, bsa_class='formal_closure'):
        r"""
        Return the basic semialgebraic set that is the topological closure
        of ``self``.
        """
        # Because our description consists of minimized constraints, the closure is
        # just the formal closure.
        return self.formal_closure(bsa_class=bsa_class)

    def relint(self, bsa_class='formal_relint'):
        r"""
        Return the basic semialgebraic set that is the topological relative interior
        of ``self``.
        """
        # Because our description consists of minimized constraints, the relint is
        # just the formal relint.
        return self.formal_relint(bsa_class=bsa_class)

    def eq_poly(self):
        r"""
        Return a list of the polynomials `f` in equations `f(x) = 0`
        in the description of ``self``.

        Together, ``eq_poly``, ``lt_poly``, and ``le_poly`` describe ``self``.
        """
        # add tests
        for c in self._polyhedron.constraints():
            if c.is_equality():
                coeff = [c.coefficient(pplite_Var(i)) for i in range(c.space_dimension())]
                # observe: coeffients in a constraint of NNC_Polyhedron could have gcd != 1.
                gcd_c = gcd(gcd(coeff), c.inhomogeneous_term())
                t = sum(QQ(x)/gcd_c*y for x, y in zip(coeff, self.poly_ring().gens())) + QQ(c.inhomogeneous_term())/gcd_c  # not type stable, make it type stable
                yield self.poly_ring()(t)

    def lt_poly(self):
        r"""
        Return a list of the polynomials `f` in strict inequalities `f(x) < 0`
        in the description of ``self``.

        Together, ``eq_poly``, ``lt_poly``, and ``le_poly`` describe ``self``.
        """
        for c in self._polyhedron.constraints():
            if c.is_strict_inequality():
                coeff = [c.coefficient(pplite_Var(i)) for i in range(c.space_dimension())]
                gcd_c = gcd(gcd(coeff), c.inhomogeneous_term())
                # constraint is written with '>', while lt_poly records '<' relation
                t = sum(-QQ(x)/gcd_c*y for x, y in zip(coeff, self.poly_ring().gens())) - QQ(c.inhomogeneous_term())/gcd_c
                yield self.poly_ring()(t)

    def le_poly(self):
        r"""
        Return a list of the polynomials `f` in inequalities `f(x) \leq 0`
        in the description of ``self``.

        Together, ``eq_poly``, ``lt_poly``, and ``le_poly`` describe ``self``.
        """
        for c in self._polyhedron.constraints():
            if c.is_nonstrict_inequality():
                coeff = [c.coefficient(pplite_Var(i)) for i in range(c.space_dimension())]
                gcd_c = gcd(gcd(coeff), c.inhomogeneous_term())
                # constraint is written with '>=', while lt_poly records '<=' relation
                t = sum(-QQ(x)/gcd_c*y for x, y in zip(coeff, self.poly_ring().gens())) - QQ(c.inhomogeneous_term())/gcd_c
                yield self.poly_ring()(t)

    # override the default implementation
    def __contains__(self, point):
        r"""
        Whether the set contains the ``point`` (vector).
        """
        rational_list = [ QQ(x) for x in point ]
        num_list = [x.numerator() for x in rational_list]
        den_list = [x.denominator() for x in rational_list]
        common_den = lcm(den_list)
        coef = [common_den // den_list[i] * num_list[i] for i in range(len(rational_list))]
        pt = ppl_point(Linear_Expression(coef, 0), common_den)
        return self._polyhedron.relation_with(pt).implies(point_is_included_pplite)

    # override the abstract methods
    def find_point(self):
        r"""
        Find a point in ``self``.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P = BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron(2)
            sage: P.add_linear_constraint([0,1],0,operator.ge)
            sage: P.add_linear_constraint([1,0],0,operator.ge)
            sage: P.add_linear_constraint([2,3],-6,operator.lt)
            sage: P.find_point()
            (1, 2/3)
        """
        # pplite has a different representation points, closure points, of NNC polys compared to ppl
        # so the find_point method yields different results

        def to_point(g, ambient_dim):
            den = g.divisor()
            # g.set_space_dimension(ambient_dim) # PPlite generators have space dim of largest dimension of variables in point expression.
            # To sum points as vectors in sagemath vectors need to have the same dimension.
            # To fix, update the points space dim to be a defined ambient dimension.
            # TODO: update after this gets fixed in pplite.
            return vector(QQ, (QQ(x)/den for x in [g.coefficient(v) for v in range(ambient_dim)])) # based on email this should in theory works

        def to_vector(g, ambient_dim):
            den = g.divisor()
            # g.set_space_dimension(ambient_dim)
            return vector(QQ, (QQ(x)/den for x in [g.coefficient(v) for v in range(ambient_dim)]))
        points = [to_point(g, self._polyhedron.space_dimension()) for g in self._polyhedron.generators()
                    if g.is_point() or g.is_closure_point()]
        rays = [to_vector(g, self._polyhedron.space_dimension()) for g in self._polyhedron.generators()
                    if g.is_ray()]
        if points:
            p = sum(points) / len(points)
            if rays:
                p += sum(rays) / len(rays)
            return p
        raise NotImplementedError("find_test_point implementation cannot handle this case")

    def add_space_dimensions_and_embed(self, space_dim_to_add):
        r"""
        Mutate ``self`` by injecting it into a higher dimensional space.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P = BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron(2)
            sage: P.add_linear_constraint([0,1],0,operator.ge)
            sage: P.add_linear_constraint([1,0],0,operator.ge)
            sage: P.add_linear_constraint([2,3],-6,operator.lt)
            sage: P.add_space_dimensions_and_embed(2)
            sage: P
            BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron([x1>=0, x0>=0, -2*x0-3*x1+6>0], names=[x0, x1, x2, x3])
            sage: P.ambient_dim()
            4
        """
        super(BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron, self).add_space_dimensions_and_embed(space_dim_to_add)
        self._polyhedron.add_space_dimensions(int(space_dim_to_add), False)

    @staticmethod
    def _pplite_constraint(lhs, cst, op):
        r"""
        Make a PPLite ``Constraint`` ``lhs`` * x + cst ``op`` 0,
        where ``lhs`` is be a vector of length ambient_dim.

        TESTS::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: from pplite import Constraint as pplite_Con
            sage: test_constraint = BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron._pplite_constraint([0,1], 0, operator.ge)
            sage: test_constraint
            x1>=0
            sage: isinstance(test_constraint, pplite_Con)
            True

        """
        lcd = lcm(lcm(x.denominator() for x in lhs), cst.denominator())
        aff_expr = sum([int((lcd*lhs[i]))*pplite_Var(i) for i in range(len(lhs))]) + int(lcd*cst)
        if op == operator.lt:
            return (aff_expr < 0)
        elif op == operator.gt:
            return (aff_expr > 0)
        elif op == operator.eq:
            return (aff_expr == 0)
        elif op == operator.le:
            return (aff_expr <= 0)
        elif op == operator.ge:
            return (aff_expr >= 0)
        else:
            raise ValueError("{} is not a supported operator".format(op))

    def linear_function_upper_bound(self, form):
        r"""
        Find an upper bound for ``form`` (a vector) on ``self``.
        This upper bound is the supremum.

        If ``self`` is empty, it returns -oo
        """

        def to_point(g):
            den = g.divisor()
            return vector(QQ, (QQ(x)/den for x in [g.coefficient(v) for v in range(g.space_dimision())]))

        def to_vector(g):
            return vector(QQ, (QQ(x) for x in [g.coefficient(v) for v in range(g.space_dimision())]))
        if self._polyhedron.is_empty():
            return -Infinity
        form = vector(form)
        for g in self._polyhedron.generators():
            if g.is_line():
                if to_vector(g) * form != 0:
                    return +Infinity
            if g.is_ray():
                if to_vector(g) * form > 0:
                    return +Infinity
        points = [to_point(g) for g in self._polyhedron.generators()
                   if g.is_point() or g.is_closure_point()]
        return max(p * form for p in points)

    def is_linear_constraint_valid(self, lhs, cst, op):
        r"""
        Whether the constraint ``lhs`` * x + cst ``op`` 0
        is satisfied for all points of ``self``,
        where ``lhs`` is be a vector of length ambient_dim.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P = BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron(2)
            sage: P.add_linear_constraint([0,1],0,operator.ge)
            sage: P.add_linear_constraint([1,0],0,operator.ge)
            sage: P.add_linear_constraint([2,3],-6,operator.lt)
            sage: P.is_linear_constraint_valid([1,1],-3,operator.lt)
            True
            sage: P.is_linear_constraint_valid([0,1],0,operator.gt)
            False
        """
        lhs = vector(lhs)
        constraint = self._pplite_constraint(lhs, cst, op)
        return self._polyhedron.relation_with(constraint).implies(poly_is_included_pplite)

    def add_linear_constraint(self, lhs, cst, op):
        r"""
        Add the constraint ``lhs`` * x + cst ``op`` 0,
        where ``lhs`` is a vector of length ambient_dim, and
        ``op`` is one of ``operator.lt``, ``operator.gt``, ``operator.eq``,
        ``operator.le``, ``operator.ge``

        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: P = BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron(2)
            sage: P.add_linear_constraint([2,3],-6,operator.gt)
            sage: sorted(P.lt_poly())
            [-2*x0 - 3*x1 + 6]
        """
        lhs = vector(lhs)
        constraint = self._pplite_constraint(lhs, cst, op)
        self._polyhedron.add_constraint(constraint)

    def is_empty(self):
        """
        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron(1)
            sage: S.add_linear_constraint([1], -1, operator.ge)
            sage: S.is_empty()
            False
            sage: S.add_linear_constraint([1], +1, operator.le)
            sage: S.is_empty()
            True
        """
        return self._polyhedron.is_empty()

    def is_universe(self):
        """
        EXAMPLES::

            sage: from cutgeneratingfunctionology.spam.basic_semialgebraic import *
            sage: S = BasicSemialgebraicSet_polyhedral_pplite_NNC_Polyhedron(1)
            sage: S.add_linear_constraint([0], 0, operator.eq)
            sage: S.is_universe()
            True
            sage: S.add_linear_constraint([1], 1, operator.le)
            sage: S.is_universe()
            False
        """
        self._polyhedron.minimize()
        return self._polyhedron.is_universe()