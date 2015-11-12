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
    # should generate the constraints from the superadditivity of phi
    # here we copy the constraints from page 175 sec 3.1
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

def copy_magic_field_elements(magic_field, elements):
    return [SymbolicRNFElement(x.val(), x.sym(), magic_field) for x in elements]
    
def find_bases(coeff_o_rhs, K):
    """
    sage: logging.disable(logging.WARNING)
    sage: K.<r0,z1,o1,o2>=SymbolicRealNumberField([1/20, 1/12, 1/5, 0])
    sage: coeff_o_rhs = constraints_PTheta3(r0, z1, o1, o2)
    sage: theta_and_magic_field = find_bases(coeff_o_rhs, K)
    sage: theta_and_magic_field
    [(((-r0 - z1)/(-r0 - 1))~,
      ((-z1)/(-r0 - 1))~,
      SymbolicRNF[(r0)~, (z1)~, (o1)~, (o2)~]),
     (((r0 + 2*z1)/(-2*r0 + 2))~,
      ((-r0 + 2*z1)/(-2*r0 + 2))~,
      SymbolicRNF[(r0)~, (z1)~, (o1)~, (o2)~]),
     (((-z1)/(r0 - 1))~,
      ((-z1)/(r0 - 1))~,
      SymbolicRNF[(r0)~, (z1)~, (o1)~, (o2)~]),
     (((-r0 - 2*z1)/(-2*r0 - 2))~,
      ((-r0 - 2*z1)/(-2*r0 - 2))~,
      SymbolicRNF[(r0)~, (z1)~, (o1)~, (o2)~])]
    """
    n = len(coeff_o_rhs)
    theta_and_magic_field = []
    seen_thetas = []
    for i in range(n):
        for j in range(i, n):
            magic_field = copy(K)
            a11_K, a12_K, b1_K = coeff_o_rhs[i]
            a21_K, a22_K, b2_K = coeff_o_rhs[j]
            [a11, a12, b1, a21, a22, b2] = copy_magic_field_elements(magic_field, \
                                           [a11_K, a12_K, b1_K, a21_K, a22_K, b2_K])
            d = a11 * a22 - a12 * a21
            if d == 0:
                continue
            theta1 = (b1 * a22 - a12 * b2) / d
            theta2 = (a11 * b2 - b1 * a21) / d
            thetas = (theta1.sym(), theta2.sym())
            # I tried to use seen_thetas as a set and check if thetas in seen_thetas.
            # But it didn't work.
            # Because '(z1/r0_=).sym() == (2*z1/(2*r0)).sym()' but not 'is'?
            # See also [sage-trac] #16993: Broken fraction field of rational polynomial ring
            # and [sage-trac] #15297: Elements from a Field of Fractions that compare equal
            #                         should have equal hashes
            seen = False
            for t in seen_thetas:
                if t == thetas:
                    seen = True
                    break
            if seen:
                continue
            seen_thetas.append(thetas)
            feasibility = True
            for (c_o1_K, c_o2_K, c_rhs_K) in coeff_o_rhs:
                [c_o1, c_o2, c_rhs] = copy_magic_field_elements(magic_field, \
                                      [c_o1_K, c_o2_K, c_rhs_K])
                if c_o1 * theta1 + c_o2 * theta2 > c_rhs:
                    feasibility = False
                    break
            if feasibility:
                #print i, j, theta1, theta2
                theta_and_magic_field.append((theta1, theta2, magic_field))
    return theta_and_magic_field
                
def setup_vertex_group_function(theta1, theta2, r0_val, z1_val):
    """
    sage: logging.disable(logging.WARNING)
    sage: K.<r0,z1,o1,o2>=SymbolicRealNumberField([1/20, 1/12, 1/5, 0])
    sage: coeff_o_rhs = constraints_PTheta3(r0, z1, o1, o2)
    sage: theta_and_magic_field = find_bases(coeff_o_rhs, K)
    sage: theta1, theta2, magic_field = theta_and_magic_field[0]
    sage: r0_val = 1/20 # = magic_field.gens()[0].val()
    sage: z1_val = 1/12 # = magic_field.gens()[1].val()
    sage: K2, h = setup_vertex_group_function(theta1, theta2, r0_val, z1_val)
    sage: read_simplified_leq_lin(K2)
    ([], [r0 + 4*z1 - 1, -z1, -r0])
    """
    K2.<r0, z1, o1, o2> = SymbolicRealNumberField([r0_val, z1_val, 0, 0])
    vertex_o1, vertex_o2 = copy_magic_field_elements(K2, [theta1, theta2])
    h  = cpl3_group_function(r0, z1, vertex_o1, vertex_o2)
    return K2, h

def testCPL3(K):
    """
    sage: logging.disable(logging.WARNING)
    sage: K.<r0,z1,o1,o2>=SymbolicRealNumberField([1/20, 1/12, 1/5, 0])
    sage: fields, thetas, chambres, extremality = testCPL3(K)

    The following four PTheta3= vertices (theta1, theta2) correspond to
    extreme points (c), (d), (a), (b) respectively (cf. page 177 table 2)
    sage: thetas
    [(((-r0 - z1)/(-r0 - 1))~, ((-z1)/(-r0 - 1))~),
     (((r0 + 2*z1)/(-2*r0 + 2))~, ((-r0 + 2*z1)/(-2*r0 + 2))~),
     (((-z1)/(r0 - 1))~, ((-z1)/(r0 - 1))~),
     (((-r0 - 2*z1)/(-2*r0 - 2))~, ((-r0 - 2*z1)/(-2*r0 - 2))~)]

    The extremality of the funciton in the chambre around the testpoint
    sage: extremality
    [True, False, True, True]

    Compare the returned chambres around testpoint with 
    the known conditions for extremality [page 188, table 5]
    Ext pt (c), expected 3*r0 + 4*z1 <=1, got (-r0*z1<1,) r0<z1, 0<r0, 3*r0 + 4*z1 < 1
    Ext pt (d), expected r0=2*z1, r0+8*z1<=1, got a chambre where the function is not extreme
    Ext pt (a), expected 0<r0<1, 0<=z1<1, got r0 + 4*z1 < 1, 0 < z1, 0 < r0 < 1/2
    Ext pt (b), expected 3*r0 + 8*z1 <=1, got 3*r0 + 8*z1 < 1, r0<2*z1, 0<r0, (-2*r0*z1< 1)
    sage: chambres
    [([], [-r0*z1 - 1, r0 - z1, -r0, 3*r0 + 4*z1 - 1]),
     ([-5*r0 + 3*z1],
      [43*r0 - 3,
       10*r0^2*z1 - 16*r0*z1^2 + 6*z1^3 + r0^2 - 6*r0*z1 + z1^2,
       -40000*r0^2*z1 + 40000*r0*z1^2 - r0^2 + 40001*r0*z1 - 39998*z1^2,
       -2*r0,
       30*r0*z1 - 18*z1^2 - 28*r0,
       79990*r0*z1 + 6*z1^2 - 133324*r0,
       3*r0^2 - 6*r0*z1 - z1^2]),
     ([], [r0 + 4*z1 - 1, -z1, -r0, 2*r0 - 1]),
     ([], [3*r0 + 8*z1 - 1, r0 - 2*z1, -r0, -2*r0*z1 - 1])]

    Let's try another testpoint.
    In the following example, the third PTheta3= vertex (theta1, theta2) is extreme point (c).
    The function is extreme in the chambre defined by (-r0*z1<1,) 0<z1, z1<r0, 3*r0 + 4*z1 < 1
    Together with the chambre ext pt (c) in the example above,
    the conditions for extreme function for ext point (c) in the paper are completely recovered
    sage: K.<r0,z1,o1,o2>=SymbolicRealNumberField([1/6, 1/12, 0, 0])
    sage: fields, thetas, chambres, extremality = testCPL3(K)
    sage: thetas[2]
    (((-r0 - z1)/(-r0 - 1))~, ((-z1)/(-r0 - 1))~)
    sage: extremality[2]
    True
    sage: chambres[2]
    ([], [-r0*z1 - 1, -z1, -r0 + z1, 3*r0 + 4*z1 - 1])
    """
    fields = []
    thetas = []
    chambres = []
    extremality = []
    coeff_o_rhs = constraints_PTheta3(K.gens()[0], K.gens()[1], K.gens()[2], K.gens()[3])
    r0_val = K.gens()[0].val()
    z1_val = K.gens()[1].val()
    theta_and_magic_field = find_bases(coeff_o_rhs, K)
    for theta1, theta2, magic_field in  theta_and_magic_field:
        K2, h = setup_vertex_group_function(theta1, theta2, r0_val, z1_val)
        is_extreme =  extremality_test(h)
        fields.append(K2)
        thetas.append((theta1, theta2))
        leq, lin = read_simplified_leq_lin(K2)
        chambres.append((leq, lin))
        extremality.append(is_extreme)
    return fields, thetas, chambres, extremality

def find_all_vertex_thetas(sampling=20):
    """
    Finds too many vertex thetas, to check!
    sage: logging.disable(logging.WARNING)
    sage: vertex_thetas = find_all_vertex_thetas()
    sage: len(vertex_thetas)
    27
    sage: vertex_thetas
    [((-r0 - z1)/(-r0 - 1), (-z1)/(-r0 - 1)),
     ((r0 + 2*z1)/(-2*r0 + 2), (-r0 + 2*z1)/(-2*r0 + 2)),
     ((-z1)/(r0 - 1), (-z1)/(r0 - 1)),
     ((-r0 - 2*z1)/(-2*r0 - 2), (-r0 - 2*z1)/(-2*r0 - 2)),
     ((4*z1 - 1)/(6*r0 + 24*z1 - 6), (2*r0 + 4*z1 - 1)/(6*r0 + 24*z1 - 6)),
     ((-r0 - z1)/(3*r0 + 12*z1 - 3), (2*r0 + 5*z1 - 1)/(3*r0 + 12*z1 - 3)),
     (1/6, 1/6),
     ((-z1)/(3*r0 - 1), (r0 - z1)/(3*r0 - 1)),
     (z1/(-2*r0 + 1), z1/(-2*r0 + 1)),
     ((-r0)/(-r0 + 4*z1 - 1), 0),
     (2*z1/(-r0 + 1), 0),
     ((-z1)/(3*r0 + 12*z1 - 3), (r0 + 5*z1 - 1)/(3*r0 + 12*z1 - 3)),
     (z1/(-r0 + 4*z1), (-r0 + z1)/(-r0 + 4*z1)),
     ((-z1)/(-r0 - 8*z1 + 1), (-r0 - 5*z1 + 1)/(-r0 - 8*z1 + 1)),
     (1/4, 1/4),
     ((-r0^2 - 6*r0*z1 - 4*z1^2 + r0 + z1)/(-r0^2 - 8*r0*z1 - 4*z1 + 1),
      (-2*r0*z1 - 4*z1^2 + z1)/(-r0^2 - 8*r0*z1 - 4*z1 + 1)),
     ((-2*r0*z1 - 5*z1^2 + z1)/(2*r0^2 + 7*r0*z1 - 3*r0 - 5*z1 + 1),
      (-r0*z1 - 5*z1^2 + z1)/(2*r0^2 + 7*r0*z1 - 3*r0 - 5*z1 + 1)),
     ((2*r0 + 5*z1 - 1)/(4*r0 + 12*z1 - 2), z1/(4*r0 + 12*z1 - 2)),
     (1/3, 0),
     ((r0 + 4*z1 - 1)/(7*r0 + 20*z1 - 5), (2*r0 + 4*z1 - 1)/(7*r0 + 20*z1 - 5)),
     (1/5, 1/5),
     (z1/(-2*r0 + 1), 0),
     ((-z1)/(r0 + 2*z1 - 1), 0),
     ((-2*z1)/(4*r0 - 2), (r0 - 2*z1)/(4*r0 - 2)),
     (2*z1/(-4*r0 + 2), (-2*r0 - 2*z1 + 1)/(-4*r0 + 2)),
     (1/3, 1/6),
     (1/2, 0)]
    """
    # cannot use set because of
    #sage: P.<x,y,z> = QQ[]
    #sage: (-x)/(-y) in set([x/y])
    #False
    vertex_thetas = []
    for r0_num in range(1, sampling):
        for z1_num in range(1, ceil((sampling-r0_num)/4)):
            r0_val = r0_num / sampling
            z1_val = z1_num / sampling
            K.<r0,z1,o1,o2>=SymbolicRealNumberField([r0_val, z1_val, 0, 0])
            coeff_o_rhs = constraints_PTheta3(r0,z1,o1,o2)
            for theta1, theta2, magic_field in find_bases(coeff_o_rhs, K):
                thetas = (theta1.sym(), theta2.sym())
                seen = false
                for t in vertex_thetas:
                    if t == thetas:
                        seen = True
                        break
                if not seen:
                    vertex_thetas.append(thetas)
    return vertex_thetas

def random_r0z1_chambres_on_diagram_ij(i=1, j=5, random_points=100, plot_points=500, starting_graph=None):
    """
    0 <= i < j <=8
    sage: diagram15 = random_r0z1_chambres_on_diagram_ij() # not tested
    sage: diagram15.save("cpl3_diagram_1_5.pdf", title="cpl3 diagram, basis=(1, 5), yellow=not vertex thetas, red = non extreme, blue = extreme") # not tested
    """
    logging.disable(logging.WARNING)
    if starting_graph is None:
        g = Graphics()
        g += line([(0,1/4),(1,0)], color='black')
    else:
        g = starting_graph
    for num_test_points in range(random_points):
        r0_val = QQ(uniform(0, 1))
        z1_val = QQ(uniform(0, 1))
        while not z1_val <= (1-r0_val)/4:
            r0_val = QQ(uniform(0, 1))
            z1_val = QQ(uniform(0, 1))
        g += draw_r0z1_chambre_on_diagram_ij(r0_val, z1_val, i, j, plot_points)
    return g

def find_basis_thetas_on_diagram_ij(r0_val=1/6, z1_val=1/12, i=1, j=5):
    K.<r0,z1,o1,o2>=SymbolicRealNumberField([r0_val, z1_val, 0, 0])
    try:
        phi = cpl3_lifting_function(r0, z1, o1, o2)
    except:
        return K, (None, None), None
    rnf_constraints = [ -o2 - phi(r0-z1+1) + 1,\
                        2*o1 - phi(2*r0 + 2*z1), \
                        -2*o1 - o2 + phi(r0 + 3*z1), \
                        2*o1 + o2 - phi(2*r0 + 3*z1), \
                        -2*o1 - 2*o2 + phi(r0 + 4*z1), \
                        2*o1 + 2*o2 - phi(2*r0 + 4*z1), \
                        -o1 + o2, \
                        -o1, \
                        -o2 ]
    coeff_o_rhs = constraints_PTheta3(r0,z1,o1,o2)
    a11, a12, b1 = coeff_o_rhs[i]
    a21, a22, b2 = coeff_o_rhs[j]
    d = a11 * a22 - a12 * a21
    if d == 0:
        return K, (None, None), None
    theta1 = (b1 * a22 - a12 * b2) / d
    theta2 = (a11 * b2 - b1 * a21) / d
    feasibility = True
    for (c_o1, c_o2, c_rhs) in coeff_o_rhs:
        if c_o1 * theta1 + c_o2 * theta2 > c_rhs:
            feasibility = False
            break
    return K, (theta1, theta2), feasibility

def draw_r0z1_chambre_on_diagram_ij(r0_val=1/6, z1_val=1/12, i=1, j=5, plot_points=500):
    K, (theta1, theta2), feasibility = find_basis_thetas_on_diagram_ij(r0_val, z1_val, i, j)
    if feasibility is None:
        # not basis
        return point([(r0_val, z1_val)], color = 'black', size = 2, zorder=10)
    g = point([(r0_val, z1_val)], color = 'white', size = 2, zorder=10)
    if not feasibility:
        g += plot_r0z1_chambre(K, 'orange', plot_points=plot_points)
        return g
    # ?? Cannot use  K2, h = setup_vertex_group_function(theta1, theta2, r0_val, z1_val)
    # because K2 forgest the chambre for getting theta1,theta2 as functions of r0 and z1.
    h  = cpl3_group_function(K.gens()[0], K.gens()[1], theta1, theta2)
    is_extreme =  extremality_test(h)
    if is_extreme:
        g += plot_r0z1_chambre(K, 'blue', plot_points=plot_points)
    else:
        g += plot_r0z1_chambre(K, 'red', plot_points=plot_points)
    return g

def plot_r0z1_chambre(K, color, alpha=0.5, xmin=-0.1, xmax=1.1, ymin=-0.1, ymax=0.3, plot_points=1000):
    x, y = var('x, y')
    leq, lin = read_simplified_leq_lin(K)
    g =region_plot([ lhs(x, y, 0, 0) == 0 for lhs in leq ] + [ lhs(x, y, 0 , 0) < 0 for lhs in lin ], \
        (x, xmin, xmax), (y, ymin, ymax), incol=color, alpha=alpha, plot_points=plot_points, bordercol=color)
    return g

def save_random_r0z1_chambres_diagrams(random_points=100, plot_points=500):
    for i in range(9):
        for j in range(i+1, 9):
            g = random_r0z1_chambres_on_diagram_ij(i,j,random_points, plot_points)
            g.save("region_graphics/cpl3_diagram_%s_%s.pdf" %(i,j), title="cpl3 diagram, basis=(%s, %s)" % (i,j))
    return

def cpl3_r0z1_on_diagram_ij(r0_val=1/6, z1_val=1/12, i=1, j=5):
    """
    h = cpl3_r0z1_on_diagram_ij(r0_val=1/6, z1_val=1/12, i=1, j=5)
    extremality_test(h)
    """
    K, (theta1, theta2), feasibility = find_basis_thetas_on_diagram_ij(r0_val, z1_val, i, j)
    if feasibility is None:
        print "not a basis!"
        return
    o1_val = theta1.val()
    o2_val = theta2.val()
    return cpl3_group_function(r0_val, z1_val, o1_val, o2_val)
    
def random_r0z1_chambres_for_fixed_ext_pt(ext_pt='c', random_points=100, plot_points=500, starting_graph=None):
    logging.disable(logging.WARNING)
    if starting_graph is None:
        g = Graphics()
        g += line([(0,1/4),(1,0)], color='black')
    else:
        g = starting_graph
    for num_test_points in range(random_points):
        r0_val = QQ(uniform(0, 1))
        z1_val = QQ(uniform(0, 1))
        while not z1_val <= (1-r0_val)/4:
            r0_val = QQ(uniform(0, 1))
            z1_val = QQ(uniform(0, 1))
        K.<r0,z1,o1,o2>=SymbolicRealNumberField([r0_val, z1_val, 0, 0])
        thetas_of_ext_pt = {'a': (z1/(1-r0), z1/(1-r0)), \
                            'b': ((r0+2*z1)/(2+2*r0), (r0+2*z1)/(2+2*r0)), \
                            'c': ((r0+z1)/(1+r0), z1/(1+r0)), \
                            'd': ((r0+2*z1)/(2-2*r0), (2*z1-r0)/(2-2*r0)), \
                            'e': (z1/(1-2*r0), z1/(1-2*r0)), \
                            'f': ((-z1-r0+6*z1*r0+r0*r0+4*z1*z1)/(4*z1+8*z1*r0-1-r0*r0), \
                                  (z1*(2*r0+4*z1-1)/(4*z1+8*z1*r0-1+r0*r0))), \
                            'g': (z1/(1-3*r0), (z1-r0)/(1-3*r0)), \
                            'h': (1/4, 1/4),\
                            'i': ((1-2*r0-5*z1)/(2-4*r0-12*z1), (-z1/(2-4*r0-12*z1))), \
                            'j': (z1*(2*r0+5*z1-1)/(-1+3*r0+5*z1-2*r0*r0-7*z1*r0), \
                                  z1*(r0+5*z1-1)/(-1+3*r0+5*z1-2*r0*r0-7*z1*r0)), \
                            'k': (1/3, 0), \
                            'l': (z1/(1-2*r0), (2*z1-r0)/(2-4*r0)), \
                            'm': (r0/(1+r0-4*z1), 0), \
                            'n': (2*z1/(1-r0), 0), \
                            'o': (z1/(1-2*r0), (1-2*r0-2*z1)/(2-4*r0)), \
                            'p': (z1/(1-2*r0), 0), \
                            'q': (z1/(1-r0-2*z1), 0), \
                            'r': (1/2, 0), \
                           }
        (theta1, theta2) = thetas_of_ext_pt[ext_pt]
        h  = cpl3_group_function(r0, z1, theta1, theta2)
        is_extreme =  extremality_test(h)
        if is_extreme:
            g += plot_r0z1_chambre(K, 'blue', plot_points=plot_points)
        else:
            g += plot_r0z1_chambre(K, 'red', plot_points=plot_points)
    return g

def save_random_r0z1_for_ext_pts(random_points=100, plot_points=500):
    for ext_pt in ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r']:
        g = random_r0z1_chambres_for_fixed_ext_pt(ext_pt, random_points, plot_points)
        g.save("region_graphics/cpl3_ext_pt_%s.pdf" % ext_pt, title="cpl3 diagram for extreme point (%s)" % ext_pt)
    return




class Cpl3Complex(SageObject):

    def __init__(self, var_name, theta=None):
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
            self.v_dict[v] = Variable(i)
        self.graph = Graphics()
        self.num_plotted_components = 0
        self.points_to_test = set()
        self.theta = theta

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
            if is_point_in_box(var_value, c.bounds):
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
        unlifted_space_dim =  len(self.monomial_list)
        K = SymbolicRealNumberField(var_value, self.var_name)
        K.monomial_list = self.monomial_list
        K.v_dict = self.v_dict
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
        leq, lin = read_simplified_leq_lin(K)
        #if see new monomial, lift polyhedrons of the previously computed components.
        dim_to_add = len(self.monomial_list) - unlifted_space_dim
        if dim_to_add > 0:
            for c in self.components:
                c.polyhedron.add_space_dimensions_and_embed(dim_to_add)
        new_component = SemialgebraicComplexComponent(K, leq, lin, var_value, region_type)
        self.components.append(new_component)
        if flip_ineq_step > 0:
            (self.points_to_test).update(new_component.generate_neighbour_points(flip_ineq_step))

    def shoot_random_points(self, num, var_bounds=None, max_failings=1000):
        for i in range(num):
            var_value = self.find_uncovered_random_point(var_bounds=var_bounds, max_failings=max_failings)
            if not var_value is False:
                self.add_new_component(var_value, flip_ineq_step=0)
            else:
                return

    def plot(self, alpha=0.5, plot_points=300, slice_value=None, restart=False):
        if restart:
            self.graph = Graphics()
            self.num_plotted_components = 0
        for c in self.components[self.num_plotted_components::]:
            self.graph += c.plot(alpha=alpha, plot_points=plot_points, slice_value=slice_value)
        self.num_plotted_components = len(self.components)
        return self.graph

    def bfs_completion(self, var_value=None, var_bounds=None, max_failings=1000, flip_ineq_step=0.01):
        if var_value:
            (self.points_to_test).add(tuple(var_value))
        while self.points_to_test:
            var_value = list(self.points_to_test.pop())
            if not self.is_point_covered(var_value) and point_satisfies_var_bounds(var_value, var_bounds):
                self.add_new_component(var_value, flip_ineq_step=flip_ineq_step)
