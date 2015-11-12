# from cpl3.sage import Cpl3Complex, cpl3_lifting_function, constraints_PTheta3

complex=Cpl3Complex(['r0','z1'], None)
complex.bfs_completion(var_value=[6/10,4/100])

regions=[]
for c in complex.components:
    x, y = c.var_value
    if x >= 0 and x <=1 and y>=0 and y<=1/4-x/4:
        regions.append(c)

sage: len(regions)
30
sage: cpl3c = copy(complex)
sage: cpl3c.components = copy(regions)
sage: cpl3c.plot(restart=True)

for r in regions:
    r.theta_ij={}
    r.feas_ij={}
    for i in range(9):
        for j in range(i+1, 9):
            r0_val, z1_val = r.var_value
            K.<r0,z1,o1,o2>=SymbolicRealNumberField([r0_val, z1_val, 0, 0])
            phi = cpl3_lifting_function(r0, z1, o1, o2)
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


sage: [r.var_value for r in regions]
[[173929566/281939773, 36007817/382248036],
 [80767979/131663654, 72372091/941212003],
 [115267641/231256003, 13535997/158825152],
 [93353921/152718865, 15818191/247470619],
 [155978156/314143601, 61635301/606584617],
 [130874027/262708804, 35849890/428805121],
 [33958816/77010563, 29583910/262451343],
 [11442929/1320139567, 80074554/325164023],
 [56160945/151392271, 26071232/204125627],
 [76061882/230283697, 30317099/181118528],
 [62152306/138203933, 22008572/237467581],
 [164768233/353449855, 6278125/92865963],
 [62286115/157600803, 191111951/1870340064],
 [10661906/1390447175, 893892844/3659611281],
 [71617359/424605569, 65515849/394059058],
 [602731061/1827400053, 34467283/306548244],
 [130244120/290288403, 11553922/230598079],
 [59855629/212167067, 80789373/562680065],
 [31270622/186503255, 18106887/110233786],
 [23297588/2640906353, 84815005/428297632],
 [164830655/500154363, 23960570/216221143],
 [46483979/284577824, 28789343/208146307],
 [11535925/1418556061, 74410489/379057325],
 [45699005/282154742, 36360980/269610727],
 [36467573/260815246, 43147179/302316497],
 [121710915/382066426, 69202949/779179895],
 [65018857/551842147, 4095187/32316337],
 [39645285/196298446, 39269765/393229541],
 [38712983/249818664, 53855146/472974751],
 [42352143/247725019, 25506811/311572622]]



sage: thetaij
[(z1^2/(2*r0^2 + r0*z1 - r0 + z1), (-r0*z1 + z1^2)/(2*r0^2 + r0*z1 - r0 + z1)),
 ((-z1)/(3*r0 - 1), (r0 - z1)/(3*r0 - 1)),
 (z1/(-2*r0 + 1), 0),
 ((-r0)/(-r0 + 4*z1 - 1), 0),
 (1/3, 0),
 (1/2, 0)]


i=0; j=1;

for i in range(9):
    for j in range(i+1, 9):
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
        g.save("cpl_%s_%s_thetas.pdf" %(i,j))
        
g.show()

# uniq doesn't work for the case (i,j) = (1,3), (1,5), (2,4), (2,6), (3,5), (4,6) ;

thetaij = uniq([(r.theta_ij[i,j][0].sym(), r.theta_ij[i,j][1].sym()) for i in range(9) for j in range(i+1,9) for r in regions if r.feas_ij[i,j]])
theta_ext = []
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
sage: len(theta_ext)
18
sage: print theta_ext
[(z1/(-r0 - 2*z1 + 1), 0), (2*z1/(-2*r0 + 2), 2*z1/(-2*r0 + 2)), (z1/(-2*r0 + 1), 0), (2*z1/(-4*r0 + 2), (-2*r0 - 2*z1 + 1)/(-4*r0 + 2)), ((-2*z1)/(r0 - 1), 0), ((-2*r0*z1 - 5*z1^2 + z1)/(2*r0^2 + 7*r0*z1 - 3*r0 - 5*z1 + 1), (-r0*z1 - 5*z1^2 + z1)/(2*r0^2 + 7*r0*z1 - 3*r0 - 5*z1 + 1)), ((-2*z1)/(4*r0 - 2), (r0 - 2*z1)/(4*r0 - 2)), ((-z1)/(3*r0 - 1), (r0 - z1)/(3*r0 - 1)), (z1/(-2*r0 + 1), z1/(-2*r0 + 1)), ((-r0)/(-r0 + 4*z1 - 1), 0), (1/4, 1/4), (1/3, 0), ((2*r0 + 5*z1 - 1)/(4*r0 + 12*z1 - 2), z1/(4*r0 + 12*z1 - 2)), (1/2, 0), ((r0 + z1)/(r0 + 1), z1/(r0 + 1)), ((r0 + 2*z1)/(-2*r0 + 2), (-r0 + 2*z1)/(-2*r0 + 2)), ((-r0 - 2*z1)/(-2*r0 - 2), (-r0 - 2*z1)/(-2*r0 - 2)), ((-r0^2 - 6*r0*z1 - 4*z1^2 + r0 + z1)/(-r0^2 - 8*r0*z1 - 4*z1 + 1), (-2*r0*z1 - 4*z1^2 + z1)/(-r0^2 - 8*r0*z1 - 4*z1 + 1))]


sage: cpl01=Cpl3Complex(['r0','z1'],0,1)
sage: cpl01.bfs_completion(var_value=[6/10,5/100])
KeyboardInterrupt: 
sage: len(cpl01.components)
35
sage: cpl01.plot()
Launched png viewer for Graphics object consisting of 105 graphics primitives
sage: cpl01.monomial_list
[r0, z1, r0*z1, r0^2, r0^3, r0^2*z1, z1^2, r0*z1^2]
sage: g = cpl01.plot()
sage: g.show(xmin=0, xmax=1, ymin=0, ymax=1/4)
Launched png viewer for Graphics object consisting of 105 graphics primitives
sage: c01 = cpl01.components[0]
sage: cpl01.thetas
{(1/2, 0): 0}

sage: cpl02=Cpl3Complex(['r0','z1'],0,2)
sage: cpl02.bfs_completion(var_value=[6/10,5/100],var_bounds=[(-0.1, 1.1), (-0.1, 0.4)])
sage: cpl02.thetas
{}
sage: len(cpl02.components)
26

sage: cpl03=Cpl3Complex(['r0','z1'],0,3)
sage: cpl03.bfs_completion(var_value=[6/10,5/100])
sage: cpl03.thetas
{(1/2, 0): 0}
sage: len(cpl03.components)
25
sage: g = cpl03.plot()
sage: g.show(xmin=0, xmax=1, ymin=0, ymax=1/4)
Launched png viewer for Graphics object consisting of 75 graphics primitives


sage: cpl04=Cpl3Complex(['r0','z1'],0,4)
sage: cpl04.bfs_completion(var_value=[6/10,1/100])
WARNING: 2015-11-11 16:14:05,296 equation list [-r0 - 8*z1 + 1] is not empty!
KeyboardInterrupt:
sage: cpl04.thetas
{(2*z1/(-r0 + 1), 0): 0}
sage: len(cpl04.components)
13


sage: thetaij[0] == theta_ext[0]
True
sage: thetaij[0][0].parent()
Fraction Field of Multivariate Polynomial Ring in r0, z1, o1, o2 over Rational Field
sage: theta_ext[0][0].parent()
Fraction Field of Multivariate Polynomial Ring in r0, z1 over Rational Field

for k in range(len(theta_ext)):
    t = theta_ext[k]
    g = Graphics()
    g += text("extreme point %s:  theta = %s" %(k,t), (0.5, 1/4), color='black')
    #g += text("extreme point %s:\ntheta = (%s,\n %s)" %(k,t[0], t[1]), (0.5, 1/4), color='black') # for k =5 and 17
    for r in regions:
        to_add = False
        for i in range(9):
            for j in range(i+1, 9):
                if r.feas_ij[i,j] and (r.theta_ij[i,j][0].sym(), r.theta_ij[i,j][1].sym()) == t:
                    to_add = True
                    break
            if to_add:
                break
        if to_add:
            r.region_type = "red"  #is feasible vertex theta
        else:
            r.region_type = "lightgrey"  #is not feasible vertex theta
        g += r.plot()
    g.save("cpl_thetas_ext_%s.pdf" %k)


cpl3s = []
for k in range(len(theta_ext)):
    print "extreme case %s" %k
    t = theta_ext[k]
    g = Graphics()
    if k == 5 or k == 17:
        g += text("extreme point %s:\ntheta = (%s,\n %s)" %(k,t[0], t[1]), (0.5, 1/4), color='black') # for k =5 and 17
    else:
        g += text("extreme point %s:  theta = %s" %(k,t), (0.5, 1/4), color='black')
    cpl = Cpl3Complex(['r0','z1'], t)
    for r in regions:
        possible_region = False
        for i in range(9):
            for j in range(i+1, 9):
                if r.feas_ij[i,j] and (r.theta_ij[i,j][0].sym(), r.theta_ij[i,j][1].sym()) == t:
                    possible_region = True
                    break
            if possible_region:
                break
        if possible_region:
            (cpl.points_to_test).add(tuple(r.var_value))
        else:
            r.region_type = "lightgrey"  #is not feasible vertex theta
            cpl.components.append(r)

    #cpl.bfs_completion()  # this is much slower than randomly shooting points

    var_bounds = [(0,1), (0, (lambda x: 1/4-x/4))]
    cpl.shoot_random_points(1000, var_bounds=var_bounds, max_failings=10000)
    cpl3s.append(cpl)  # in case random shooting is not complete, can continue from here.
    print "extreme case %s finish" %k
    g += cpl.plot()
    g.save("cpl_extreme_case_%s.pdf" %k, xmin=0, xmax=1, ymin=0, ymax=1/4)

