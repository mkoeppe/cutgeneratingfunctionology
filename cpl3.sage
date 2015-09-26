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
    return piecewise_function_from_breakpoints_and_values(bkpts, values)

def treat_constraint_of_PTheta3(rnf_c):
    """
    The emaxple follows subcase1 of case 1, section 3.2, CPL3 paper Page 175,
    We assume that 0 < r0 < z1 and 2*r0 + 6*z1 <= 1.
    sage: K.<r0,z1,o1,o2>=SymbolicRealNumberField([1/20, 1/12, 1/5, 0])
    sage: phi = cpl3_lifting_function(r0, z1, o1, o2)
    Let rnf_c be the second constraint of PTheta3= in section 3.1:
    sage: rnf_c = 2*o1 - phi(2*r0+2*z1)
    We look for the coefficients of o1 and o2 in equation(4), and the rhs of eqn(4):
    sage: c_o1, c_o2, c_rhs = treat_constraint_of_PTheta3(rnf_c)
    sage: c_o1
    (r0 - 4*z1 + 1)~
    sage: c_o2
    (3*r0 + 4*z1 - 1)~
    sage: c_rhs
    (r0)~
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
    sage: K.<r0,z1,o1,o2>=SymbolicRealNumberField([1/20, 1/12, 1/5, 0])
    sage: coeff_o_rhs = constraints_PTheta3(r0, z1, o1, o2)
    sage: theta_and_magic_field = find_bases(coeff_o_rhs, K)
    sage: theta_and_magic_field
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
    sage: K.<r0,z1,o1,o2>=SymbolicRealNumberField([1/20, 1/12, 1/5, 0])
    sage: fields, thetas, chambres, extremality = testCPL3(K)
    # The following four PTheta3= vertices (theta1, theta2) correspond to
    # extreme points (c), (d), (a), (b) respectively (cf. page 177 table 2)
    sage: thetas
    [(((-r0 - z1)/(-r0 - 1))~, ((-z1)/(-r0 - 1))~),
     (((r0 + 2*z1)/(-2*r0 + 2))~, ((-r0 + 2*z1)/(-2*r0 + 2))~),
     (((-z1)/(r0 - 1))~, ((-z1)/(r0 - 1))~),
     (((-r0 - 2*z1)/(-2*r0 - 2))~, ((-r0 - 2*z1)/(-2*r0 - 2))~)]
    # the extremality of the funciton in the chambre around the testpoint
    sage: extremality
    [True, False, True, True]
    # page 188, table 5 shows the conditions for extremality for igp. 
    # compare with chambres around testpoint.  
    # ext pt (c), expected 3*r0 + 4*z1 <=1, got (-r0*z1<1,) r0<z1, 0<r0, 3*r0 + 4*z1 < 1
    # ext pt (d), expected r0=2*z1, r0+8*z1<=1, got a chambre where the function is not extreme
    # ext pt (a), expected 0<r0<1, 0<=z1<1, got r0 + 4*z1 < 1, 0 < z1, 0 < r0 < 1/2
    # ext pt (b), expected 3*r0 + 8*z1 <=1, got 3*r0 + 8*z1 < 1, r0<2*z1, 0<r0, (-2*r0*z1< 1)
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
    In the following example, the first PTheta3= vertex (theta1, theta2) is extreme point (c).
    The function is extreme in the chambre defined by (-r0*z1<1,) 0<z1, z1<r0, 3*r0 + 4*z1 < 1
    Together with the chambre ext pt (c) in the example above,
    the conditions for extreme function for ext point (c) in the paper are completely recovered
    sage: K.<r0,z1,o1,o2>=SymbolicRealNumberField([1/6, 1/12, 0, 0])
    sage: fields, thetas, chambres, extremality = testCPL3(K)
    sage: thetas
    [(((-r0 - z1)/(-r0 - 1))~, ((-z1)/(-r0 - 1))~),
     (((2*r0*z1 + 2*z1^2)/(-r0^2 - 2*r0*z1 + r0 + 2*z1))~,
      (2*z1^2/(-r0^2 - 2*r0*z1 + r0 + 2*z1))~),
     (((-z1)/(r0 - 1))~, ((-z1)/(r0 - 1))~),
     (((-r0 - 2*z1)/(-2*r0 - 2))~, ((-r0 - 2*z1)/(-2*r0 - 2))~)]
    sage: extremality
    [True, False, True, False]
    sage: chambres
    [([], [-r0*z1 - 1, -z1, -r0 + z1, 3*r0 + 4*z1 - 1]),
     ([-r0 + 2*z1],
      [5*r0 - 1,
       10000*r0^3 + 20000*r0^2*z1 - 10000*r0^2 - 19998*r0*z1 + 3*z1^2,
       -r0,
       40000*r0^2 + 80000*r0*z1 - 79993*r0,
       -2*r0^2 - r0*z1 + 2*z1^2]),
     ([], [r0 + 4*z1 - 1, -z1, -r0, 2*r0 - 1]),
     ([-12*z1 + 1, -6*r0 + 1], [])]

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

        

#### case d #####
# sage: K.<r0,z1,o1,o2>=SymbolicRealNumberField([1/6, 1/12, 1/5, 0])
# sage: h = cpl3_group_function(r0,z1,o1,o2)
# sage: leq, lin = read_simplified_leq_lin(K)
# sage: leq, lin
# ([],
#  [2*z1*o1 - r0*o2 - 2*z1*o2 - z1 + o2, -z1, -r0, r0 + 4*z1 - 1, -o1 + o2, -o1])
# sage: extremality_test(h)
# True
# sage: leq2, lin2 = read_simplified_leq_lin(K)
# WARNING: 2015-09-18 02:19:24,693 equation list [-r0*o1 - r0*o2 - r0 + o1, -2*z1*o1 + r0*o2 + 2*z1*o2 - r0 + o1, -o2, -r0 + 2*z1] is not empty!
# sage: leq2
# [-r0*o1 - r0*o2 - r0 + o1,
#  -2*z1*o1 + r0*o2 + 2*z1*o2 - r0 + o1,
#  -o2,
#  -r0 + 2*z1]
# sage: lin2
# [5*r0 - 1,
#  r0 - o1,
#  -3*r0 + 2*o1,
#  2*r0*z1 - 8*z1^2 - r0,
#  2*r0*z1 - 3*r0,
#  r0^2 - 2*r0 - 1]
#### after simplification: o1 = r/(1-r); o2 = 0; r = 2*z; 0 < r < 1/5 ####
#### identique to the known condition for case d in the paper ###





