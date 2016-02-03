# Make sure current directory is in path.  
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

def cpl3_group_function(r0=1/6, z1=1/12, o1=1/5, o2=0):
    if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    #if not (bool(0 <= o1) & bool(0 <= o2) & bool(o1+o2 <= 1/2)):
    #    logging.info("Conditions for a CPL-3 function are NOT satisfied.")

    if z1 < (1-r0)/4:
        bkpts = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
        phi_values = [0, 0, o1, o1+o2, 1-(o1+o2), 1-o1, 1]
    else:
        bkpts = [0, r0, r0+z1, 1-z1, 1]
        phi_values = [0, 0, o1, 1-o1, 1]
    values = [(bkpts[i] - phi_values[i])/r0 for i in range(len(bkpts))]
    return piecewise_function_from_breakpoints_and_values(bkpts, values)

def cpl3_bkpts(r0=1/6, z1=1/12):
    if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    # construct the function on [0, 2]
    if z1 < (1-r0)/4:
        bkpts = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
        values = [0, 0, 0, 0, 0, 0, 0]
    else:
        bkpts = [0, r0, r0+z1, 1-z1, 1]
        values = [0, 0, 0, 0, 0]
    # Cannot return piecewise_function_from_breakpoints_and_values(bkpts, values)
    # because we don't want to merge slopes
    # for arbitrary values of o1, o2 such as 0, 0
    list_of_pairs = [[(bkpts[i], bkpts[i+1]), linear_function_through_points([bkpts[i], values[i]], [bkpts[i+1], values[i+1]])] for i in range(len(bkpts)-1)]
    return FastPiecewise(list_of_pairs, periodic_extension=False, merge=False)

def cpl3_lifting_function(r0=1/6, z1=1/12, o1=1/5, o2=0):
    if not (bool(0 < r0 < 1) & bool(0 < z1 <= (1-r0)/4)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    # construct the function on [0, 2]
    if z1 < (1-r0)/4:
        phi_bkpts = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
        phi_values = [0, 0, o1, o1+o2, 1-(o1+o2), 1-o1, 1]
    else:
        phi_bkpts = [0, r0, r0+z1, 1-z1, 1]
        phi_values = [0, 0, o1, 1-o1, 1]
    bkpts = phi_bkpts + [x+1 for x in phi_bkpts[1::]]
    values = phi_values + [y+1 for y in phi_values[1::]]
    # Cannot return piecewise_function_from_breakpoints_and_values(bkpts, values)
    # because we don't want to merge slopes
    # for arbitrary values of o1, o2 such as 0, 0
    list_of_pairs = [[(bkpts[i], bkpts[i+1]), linear_function_through_points([bkpts[i], values[i]], [bkpts[i+1], values[i+1]])] for i in range(len(bkpts)-1)]
    return FastPiecewise(list_of_pairs, periodic_extension=False, merge=False)

class Cpl3Complex(SageObject):

    def __init__(self, var_name, theta=None, max_iter=8):
        #self.num_components = 0
        self.components = []
        self.d = len(var_name)
        self.var_name = var_name
        self.monomial_list = []
        self.v_dict = {}
        K = SymbolicRealNumberField([0]*self.d, var_name)
        for i in range(self.d):
            v = K.gens()[i].sym().numerator()
            self.monomial_list.append(v)
            self.v_dict[v] = i
        self.graph = Graphics()
        self.num_plotted_components = 0
        self.points_to_test = set()
        self.theta = theta
        # cpl region graphic looks like linear case,
        # but higher degree terms appear in monomials during the computation
        # so let's try mccormiks
        self.max_iter = max_iter

    def generate_random_var_value(self, var_bounds=None):
        var_value = []
        for i in range(self.d):
            if not var_bounds:
                x = QQ(uniform(-0.1, 1.1))
            else:
                if hasattr(var_bounds[i][0], '__call__'):
                    l =  var_bounds[i][0](*var_value)
                else:
                    l = var_bounds[i][0]
                if hasattr(var_bounds[i][1], '__call__'):
                    u =  var_bounds[i][1](*var_value)
                else:
                    u = var_bounds[i][1]
                x = QQ(uniform(l, u))
            var_value.append(x)
        return var_value

    def is_point_covered(self, var_value):
        monomial_value = [m(var_value) for m in self.monomial_list]
        # coefficients in ppl point must be integers.
        lcm_monomial_value = lcm([x.denominator() for x in monomial_value])
        #print [x * lcm_monomial_value for x in monomial_value]
        pt = Generator.point(Linear_Expression([x * lcm_monomial_value for x in monomial_value], 0), lcm_monomial_value)
        for c in self.components:
            # Check if the random_point is contained in the box.
            if c.region_type == 'not_constructible' and c.leq == [] and c.lin == []:
                continue
            if is_point_in_box(monomial_value, c.bounds):
                # Check if all eqns/ineqs are satisfied.
                if c.polyhedron.relation_with(pt).implies(point_is_included):
                    return True
        return False

    def find_uncovered_random_point(self, var_bounds=None, max_failings=1000):
        num_failings = 0
        while not max_failings or num_failings < max_failings:
            if self.points_to_test:
                var_value = list(self.points_to_test.pop())
            else:
                var_value = self.generate_random_var_value(var_bounds=var_bounds)
            # This point is not already covered.
            if self.is_point_covered(var_value):
                num_failings += 1
            else:
                return var_value
        logging.warn("The graph has %s components. Cannot find one more uncovered point by shooting %s random points" % (len(self.components), max_failings))
        return False

    def add_new_component(self, var_value, flip_ineq_step=0):
        # Remark: the sign of flip_ineq_step indicates how to search for neighbour testpoints:
        # if flip_ineq_step = 0, don't search for neighbour testpoints. Used in shoot_random_points().
        # if flip_ineq_step < 0, we assume that the walls of the cell are linear eqn/ineq over original parameters. (So, gradient is constant; easy to find a new testpoint on the wall and another testpoint (-flip_ineq_step away) across the wall.) Used in bfs.
        # if flip_ineq_step > 0, we don't assume the walls are linear. Apply generate_one_point_by_flipping_inequality() with flip_ineq_step to find new testpoints across the wall only. Used in bfs.

        unlifted_space_dim =  len(self.monomial_list)
        K = SymbolicRealNumberField(var_value, self.var_name)
        K.monomial_list = self.monomial_list # change simultaneously while lifting
        K.v_dict = self.v_dict # change simultaneously while lifting
        K.polyhedron.add_space_dimensions_and_embed(len(K.monomial_list))
        if self.theta is None:
            try:
                h = cpl3_bkpts(*K.gens())
                t = subadditivity_test(h) # always True
                region_type = 'not_minimal'
            except:
                # Function is non-contructible at this random point.
                region_type = 'not_constructible'
        else:
            o1 = self.theta[0](*K.gens())
            o2 = self.theta[1](*K.gens())
            try:
                h = cpl3_group_function(K.gens()[0], K.gens()[1], o1, o2)
                # t = minimality_test(h) # no need, since cpl3 satisfies minimality
                t = extremality_test(h)
                if t:
                    region_type = 'is_extreme'
                else:
                    region_type = 'not_extreme'
            except:
                region_type = 'not_constructible'
        new_component = SemialgebraicComplexComponent(self, K, var_value, region_type)
        # Temporary code to check if everything is linear:
        if flip_ineq_step < 0:
            for l in (new_component.lin+new_component.leq):
                if l.degree() > 1:
                    raise NotImplementedError, "Alas, non-linear term appeared."
        #if see new monomial, lift polyhedrons of the previously computed components.
        dim_to_add = len(self.monomial_list) - unlifted_space_dim
        if dim_to_add > 0:
            for c in self.components:
                c.polyhedron.add_space_dimensions_and_embed(dim_to_add)
        self.components.append(new_component)
        if (region_type != 'not_constructible') and (flip_ineq_step != 0):
            #(self.points_to_test).update(new_component.generate_neighbour_points(flip_ineq_step))
            neighbour_points = new_component.generate_neighbour_points(flip_ineq_step)
            (self.points_to_test).update(neighbour_points)
            ## debug
            #(self.plot()+point(neighbour_points)).show(xmin=0, xmax=1, ymin=0, ymax=1/4)
            #print new_component.var_value, new_component.region_type
            #print new_component.lin, new_component.leq
            

    def shoot_random_points(self, num, var_bounds=None, max_failings=1000):
        for i in range(num):
            var_value = self.find_uncovered_random_point(var_bounds=var_bounds, max_failings=max_failings)
            if var_value is False:
                return
            else:
                self.add_new_component(var_value, flip_ineq_step=0)

    def plot(self, alpha=0.5, plot_points=300, slice_value=None, restart=False):
        if restart:
            self.graph = Graphics()
            self.num_plotted_components = 0
        for c in self.components[self.num_plotted_components::]:
            self.graph += c.plot(alpha=alpha, plot_points=plot_points, slice_value=slice_value)
        self.num_plotted_components = len(self.components)
        return self.graph

    def bfs_completion(self, var_value=None, var_bounds=None, max_failings=1000, flip_ineq_step=-1/100):
        # See remark about flip_ineq_step in def add_new_component().
        if var_value:
            (self.points_to_test).add(tuple(var_value))
        while self.points_to_test:
            var_value = list(self.points_to_test.pop())
            if not self.is_point_covered(var_value) and point_satisfies_var_bounds(var_value, var_bounds):
                self.add_new_component(var_value, flip_ineq_step=flip_ineq_step)

def regions_r0_z1_from_arrangement_of_bkpts(show_plots=False):
    """
    sage: logging.disable(logging.INFO)  # not tested
    sage: regions = regions_r0_z1_from_arrangement_of_bkpts() # not tested
    sage: len(regions) #not tested
    30

    Figure clp_30_regions is obtained by setting show_plots=True
    """
    complex=Cpl3Complex(['r0','z1'], None)
    complex.bfs_completion(var_value=[6/10,4/100])
    regions=[]
    for c in complex.components:
        x, y = c.var_value
        if x > 0 and x < 1 and y > 0 and y<=1/4-x/4:
            regions.append(c)
    if show_plots:
        complex.plot().show(xmin=0, xmax=1, ymin=0, ymax=1/4)
    regions.sort(key=lambda r: len(r.leq))
    return regions

def treat_constraint_of_PTheta3(rnf_c):
    """
    The emaxple follows subcase1 of case 1, section 3.2, CPL3 paper Page 175,
    We assume that 0 < r0 < z1 and 2*r0 + 6*z1 <= 1.
    sage: logging.disable(logging.WARNING)
    sage: K.<r0,z1,o1,o2>=SymbolicRealNumberField([1/20, 1/12, 1/5, 0])
    sage: phi = cpl3_lifting_function(r0, z1, o1, o2)

    Let rnf_c be the second constraint of PTheta3= in section 3.1:
    sage: rnf_c = 2*o1 - phi(2*r0+2*z1)

    We look for the coefficients of o1 and o2 in equation(4), and the rhs of eqn(4):
    sage: c_o1, c_o2, c_rhs = treat_constraint_of_PTheta3(rnf_c)
    sage: c_o1
    ((-r0 + 4*z1 - 1)/(r0 + 4*z1 - 1))~
    sage: c_o2
    ((-3*r0 - 4*z1 + 1)/(r0 + 4*z1 - 1))~
    sage: c_rhs
    ((-r0)/(r0 + 4*z1 - 1))~
    """
    sym_c = rnf_c.sym()
    c = sym_c.numerator()
    r0_sym, z1_sym, o1_sym, o2_sym = c.args()
    c_rhs_sym_numerator = - (c.reduce([o1_sym, o2_sym]))
    c_rhs_sym = c_rhs_sym_numerator / sym_c.denominator()
    c_o1_sym = (c + c_rhs_sym_numerator).mod(o2_sym) / o1_sym / sym_c.denominator()
    c_o2_sym = (c + c_rhs_sym_numerator).mod(o1_sym) / o2_sym / sym_c.denominator()
    magic_field = rnf_c.parent()
    v = magic_field._values
    c_rhs = SymbolicRNFElement(c_rhs_sym(v), c_rhs_sym, magic_field)
    c_o1 = SymbolicRNFElement(c_o1_sym(v), c_o1_sym, magic_field)
    c_o2 = SymbolicRNFElement(c_o2_sym(v), c_o2_sym, magic_field)
    return c_o1, c_o2, c_rhs

# Approach that uses simplified description of PTheta3 stated on page 175 sec 3.1 of the paper.
# This approach allows plotting the colorful diagrams cpl_i_j_thetas.

def constraints_PTheta3(r0, z1, o1, o2):
    """
    sage: logging.disable(logging.WARNING)
    sage: K.<r0,z1,o1,o2>=SymbolicRealNumberField([1/20, 1/12, 1/5, 0])
    sage: constraints_PTheta3(r0, z1, o1, o2)
    [(((-r0 + z1)/z1)~, (-1)~, 0~),
     (((-r0 + 4*z1 - 1)/(r0 + 4*z1 - 1))~,
      ((-3*r0 - 4*z1 + 1)/(r0 + 4*z1 - 1))~,
      ((-r0)/(r0 + 4*z1 - 1))~),
     (((-r0 - 2*z1 + 1)/(r0 + 4*z1 - 1))~,
      (2*z1/(r0 + 4*z1 - 1))~,
      (z1/(r0 + 4*z1 - 1))~),
     (((-r0 + 2*z1 - 1)/(r0 + 4*z1 - 1))~,
      ((-2*r0 - 2*z1)/(r0 + 4*z1 - 1))~,
      ((-r0 - z1)/(r0 + 4*z1 - 1))~),
     (((-r0 + 1)/(r0 + 4*z1 - 1))~,
      ((-r0 + 1)/(r0 + 4*z1 - 1))~,
      (2*z1/(r0 + 4*z1 - 1))~),
     (((-r0 - 1)/(r0 + 4*z1 - 1))~,
      ((-r0 - 1)/(r0 + 4*z1 - 1))~,
      ((-r0 - 2*z1)/(r0 + 4*z1 - 1))~),
     ((-1)~, 1~, 0~),
     ((-1)~, 0~, 0~),
     (0~, (-1)~, 0~)]
    """
    phi = cpl3_lifting_function(r0, z1, o1, o2)
    # here we copy the constraints from page 175 sec 3.1
    # to generate the constraints from the superadditivity of phi,
    # see generate_regions_and_theta_ext()
    rnf_constraints = [ -o2 - phi(r0-z1+1) + 1,\
                        2*o1 - phi(2*r0 + 2*z1), \
                        -2*o1 - o2 + phi(r0 + 3*z1), \
                        2*o1 + o2 - phi(2*r0 + 3*z1), \
                        -2*o1 - 2*o2 + phi(r0 + 4*z1), \
                        2*o1 + 2*o2 - phi(2*r0 + 4*z1), \
                        -o1 + o2, \
                        -o1, \
                        -o2 ]
    coeff_o_rhs = [treat_constraint_of_PTheta3(rnf_c) for rnf_c in rnf_constraints]
    return coeff_o_rhs

def regions_r0_z1_with_thetas_and_feasibility(regions=None):
    """
    sage: regions = regions_r0_z1_with_thetas_and_feasibility() # not tested
    # long time
    sage: len(regions) #not tested
    30
    """
    if regions is None:
        regions = regions_r0_z1_from_arrangement_of_bkpts()
    for r in regions:
        r.theta_ij={}
        r.feas_ij={}
        for i in range(9):
            for j in range(i+1, 9):
                r0_val, z1_val = r.var_value
                K.<r0,z1,o1,o2>=SymbolicRealNumberField([r0_val, z1_val, 0, 0])
                coeff_o_rhs = constraints_PTheta3(r0,z1,o1,o2)
                a11, a12, b1 = coeff_o_rhs[i]
                a21, a22, b2 = coeff_o_rhs[j]
                d = a11 * a22 - a12 * a21
                if d == 0:
                    r.theta_ij[i,j] = (None, None)
                    r.feas_ij[i,j] = None
                    continue
                theta1 = (b1 * a22 - a12 * b2) / d
                theta2 = (a11 * b2 - b1 * a21) / d
                feasibility = True
                for (c_o1, c_o2, c_rhs) in coeff_o_rhs:
                    if c_o1 * theta1 + c_o2 * theta2 > c_rhs:
                        feasibility = False
                        break
                r.theta_ij[i,j] = (theta1, theta2)
                r.feas_ij[i,j] = feasibility
    return regions

def plot_cpl_i_j_thetas_diagram(regions, i, j):
    """
    sage: regions = regions_r0_z1_with_thetas_and_feasibility() # not tested
    sage: i = 0; j = 1; # not tested
    sage: g = plot_cpl_i_j_thetas_diagram(regions, i, j) # not tested
    sage: g.save("cpl_%s_%s_thetas.pdf" %(i,j)) # not tested
    """
    assert (0 <= i < j < 9)
    # uniq doesn't work for the case (i,j) = (1,3), (1,5), (2,4), (2,6), (3,5), (4,6)
    #for (i,j) in [(1,3), (1,5), (2,4), (2,6), (3,5), (4,6)]:
    thetaij_dup = uniq([(r.theta_ij[i,j][0].sym(), r.theta_ij[i,j][1].sym()) for r in regions if r.feas_ij[i,j]])
    thetaij = []
    for d in thetaij_dup:
        to_add = True
        for t in thetaij:
            if t == d:
                to_add = False
                break
        if to_add:
            thetaij.append(d)
    n = len(thetaij)
    for r in regions:
        feas = r.feas_ij[i,j]
        if not feas:
            r.region_type = 'gray'
        else:
            theta1, theta2 = r.theta_ij[i,j][0].sym(), r.theta_ij[i,j][1].sym()
            for k in range(n):
                if (theta1, theta2) == thetaij[k]:
                    break
            r.region_type=rainbow(n)[k]
    g = Graphics()
    for r in regions:
        g += r.plot()
    for k in range(n):
        t =  text(thetaij[k], (0.5, 1/4 + k/25), color = rainbow(n)[k])
        g += t
    return g

def retrieve_theta_ext_from_regions(regions):
    """
    sage: regions = regions_r0_z1_with_thetas_and_feasibility() # not tested
    sage: theta_ext = retrieve_theta_ext_from_regions(regions); # not tested
    sage: len(theta_ext) #not tested
    18
    sage: theta_ext[0] # not tested
    (z1/(-r0 - 2*z1 + 1), 0)
    """
    thetaij = uniq([(r.theta_ij[i,j][0].sym(), r.theta_ij[i,j][1].sym()) \
                    for i in range(9) for j in range(i+1,9) for r in regions if r.feas_ij[i,j]])
    theta_ext = []
    for r in regions:
        r.thetas = set([])
    K.<r0,z1>=QQ[]
    for (t1, t2) in thetaij:
        d1 = t1(r0, z1, 0, 0)
        d2 = t2(r0, z1, 0, 0)
        d = (d1, d2)
        to_add = True
        for t in theta_ext:
            if t == d:
                to_add = False
                break
        if to_add:
            theta_ext.append(d)
            for r in regions:
                for i in range(9):
                    for j in range(i+1, 9):
                        if r.feas_ij[i,j] and (r.theta_ij[i,j][0].sym(), r.theta_ij[i,j][1].sym()) == d:
                            (r.thetas).add(d)
    return theta_ext

# Approach that does not use simplified description of PTheta3 stated on page 175 sec 3.1 of the paper.
# Instead, it derives everything starting from the superadditivity of phi.
def generate_coefficients_of_constraints_PTheta3(r0_val, z1_val):
    """
    sage: generate_coefficients_of_constraints_PTheta3(3/5, 1/25) # not tested
    [(2~, 2~, 1~),
     (2~, 0~, 1~),
     ((-1)~, (-1)~, 0~),
     ((-1)~, 0~, 0~),
     ((-1)~, 1~, 0~),
     (0~, 0~, 0~),
     (((-r0 + 1)/(r0 + 4*z1 - 1))~,
      ((-r0 + 1)/(r0 + 4*z1 - 1))~,
      (2*z1/(r0 + 4*z1 - 1))~),
     (0~, 0~, 1~),
     (0~, (-1)~, 0~),
     (((-r0 - 2*z1 + 1)/(r0 + 4*z1 - 1))~,
      (2*z1/(r0 + 4*z1 - 1))~,
      (z1/(r0 + 4*z1 - 1))~),
     (1~, 1~, 1~),
     (2~, 1~, 1~),
     (1~, 0~, 1~)]
    """
    K.<r0,z1,o1,o2>=SymbolicRealNumberField([r0_val, z1_val, 0, 0])
    phi = cpl3_lifting_function(r0, z1, o1, o2)
    rnf_constraints = set([])
    bkpts = phi.end_points()
    for x in bkpts:
        for y in bkpts:
            if x <= y < 1:
                rnf_c = phi(x) + phi(y) - phi(x+y)
                rnf_constraints.add(rnf_c)
    for x in bkpts:
        for z in bkpts:
            if x < 1 and x < z < 1+x:
                rnf_c = phi(x) + phi(z-x) - phi(z)
                rnf_constraints.add(rnf_c)
    coeff_o_rhs = [treat_constraint_of_PTheta3(rnf_c) for rnf_c in rnf_constraints]
    return coeff_o_rhs

def generate_regions_and_theta_ext():
    """
    sage: regions, theta_ext = generate_regions_and_theta_ext() # not tested
    sage: len(regions) # not tested
    30  -> 96
    sage: len(theta_ext) # not tested
    18
    """
    regions = regions_r0_z1_from_arrangement_of_bkpts()
    theta_ext = []
    K.<r0,z1>=QQ[]
    for r in regions:
        r.thetas = set([])
        r0_val, z1_val = r.var_value
        coeff_o_rhs = generate_coefficients_of_constraints_PTheta3(r0_val, z1_val)
        n = len(coeff_o_rhs)
        for i in range(n):
            for j in range(i+1, n):
                a11, a12, b1 = coeff_o_rhs[i]
                a21, a22, b2 = coeff_o_rhs[j]
                d = a11 * a22 - a12 * a21
                if d == 0:
                    continue
                theta1 = (b1 * a22 - a12 * b2) / d
                theta2 = (a11 * b2 - b1 * a21) / d
                feasibility = True
                for (c_o1, c_o2, c_rhs) in coeff_o_rhs:
                    if c_o1 * theta1 + c_o2 * theta2 > c_rhs:
                        feasibility = False
                        break
                if not feasibility:
                    continue
                d = (theta1.sym()(r0, z1, 0, 0), theta2.sym()(r0, z1, 0, 0))
                to_add = (not r.leq) ### too many theta_ext for lower dim case.
                for t in theta_ext:
                    if t == d:
                        to_add = False
                        break
                if to_add:
                    theta_ext.append(d)
                    (r.thetas).add(d)
                else:
                    (r.thetas).add(t)
    return regions, theta_ext

# Plotting diagrams
def plot_cpl_thetas_ext_diagram(regions, t, k):
    """
    # either approach in paper:
    sage: regions = regions_r0_z1_with_thetas_and_feasibility() # not tested
    sage: theta_ext = retrieve_theta_ext_from_regions(regions); # not tested
    # or better direct approach:
    sage: regions, theta_ext = generate_regions_and_theta_ext() # not tested

    sage: k = 0; t = theta_ext[k]; # not tested
    sage: g = plot_cpl_thetas_ext_diagram(regions, t, k) # not tested
    sage: g.save("cpl_thetas_ext_%s.pdf" %k) # not tested
    """
    g = Graphics()
    if len(str(t)) > 100:
        g += text("extreme point %s:\ntheta = (%s,\n %s)" %(k,t[0], t[1]), (0.5, 1/4), color='black')
    else:
        g += text("extreme point %s:  theta = %s" %(k,t), (0.5, 1/4), color='black')
    for r in regions:
        if t in r.thetas:
            r.region_type = "red"  #is feasible vertex theta
        else:
            r.region_type = "lightgrey"  #is not feasible vertex theta
        g += r.plot()
    return g

def complex_of_cpl_extreme_case_k(regions, t):
    """
    # either approach in paper:
    sage: regions = regions_r0_z1_with_thetas_and_feasibility() # not tested
    sage: theta_ext = retrieve_theta_ext_from_regions(regions); # not tested
    # or better direct approach:
    sage: regions, theta_ext = generate_regions_and_theta_ext() # not tested

    sage: k = 0; t = theta_ext[k]; # not tested
    sage: complex = complex_of_cpl_extreme_case_k(regions, t) # not tested
    """
    complex = Cpl3Complex(['r0','z1'], t)
    for r in regions:
        if t in r.thetas:
            (complex.points_to_test).add(tuple(r.var_value))
        else:
            r.region_type = "lightgrey"  #is not feasible vertex theta
            complex.components.append(r)
    complex.bfs_completion()
    #var_bounds = [(0,1), (0, (lambda x: 1/4-x/4))]
    #complex.shoot_random_points(1000, var_bounds=var_bounds, max_failings=10000)
    return complex

def plot_cpl_extreme_case_k_diagram(complex, t, k):
    """
    # either approach in paper:
    sage: regions = regions_r0_z1_with_thetas_and_feasibility() # not tested
    sage: theta_ext = retrieve_theta_ext_from_regions(regions); # not tested
    # or better direct approach:
    sage: regions, theta_ext = generate_regions_and_theta_ext() # not tested

    sage: k = 0; t = theta_ext[k]; # not tested
    sage: complex = complex_of_cpl_extreme_case_k(regions, t) # not tested
    sage: g = plot_cpl_extreme_case_k_diagram(complex, t, k) # not tested
    sage: g.save("cpl_extreme_case_%s.pdf" %k, xmin=0, xmax=1, ymin=0, ymax=1/4)
    """
    g = Graphics()
    if len(str(t)) > 100:
        g += text("extreme point %s:\ntheta = (%s,\n %s)" %(k,t[0], t[1]), (0.5, 1/4), color='black')
    else:
        g += text("extreme point %s:  theta = %s" %(k,t), (0.5, 1/4), color='black')
    g += complex.plot()
    return g
