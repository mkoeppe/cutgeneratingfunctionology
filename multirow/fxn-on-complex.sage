if '' not in sys.path:
    sys.path = [''] + sys.path
from igp import *
import random
import itertools

def FastPiecewise_to_2d_diagonal_to_function_on_complex(f):
    """
    Take a FastPiecewise function, convert it to diagonal 2d with pi_tilda(x,y) = pi(x+y), and then to function_on_complex
    Convert the 1d interval where pi is defined to a 2d [0,1]*[0,1] area where pi_tilda is defined
    
    Example::
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = gmic()
        sage: F = FastPiecewise_to_2d_diagonal_to_function_on_complex(h)
        sage: F.minimality_test()
        pi is minimal

        sage: h = dg_2_step_mir_limit()
        sage: F = FastPiecewise_to_2d_diagonal_to_function_on_complex(h)
        sage: F.minimality_test()
        pi is minimal

        sage: h = drlm_3_slope_limit()
        sage: F = FastPiecewise_to_2d_diagonal_to_function_on_complex(h)
        sage: F.minimality_test()
        pi is minimal

        sage: h = not_minimal_1()
        sage: F = FastPiecewise_to_2d_diagonal_to_function_on_complex(h)
        sage: F.minimality_test()
        Minimality test fails because subadditivity test fails.

    """      
    R.<x,y> = PolynomialRing(QQ)

    #construct the basic trapezoids(or triangles) with regard to the intervals of f
    #for the lower triangle
    polyhedrons = [Polyhedron(vertices = [[intv[0],0],[0,intv[0]],[intv[1],0],[0,intv[1]]]) for intv in f.intervals()]
    #for the upper triangle
    polyhedrons.extend([Polyhedron(vertices = [[intv[0],1],[1,intv[0]],[intv[1],1],[1,intv[1]]]) for intv in f.intervals()])

    #the diagonal function pi_tilda(x,y) = pi(x+y)
    functions = [fxn._slope * (x+y) + fxn._intercept for fxn in f.functions()]
    functions2 = [fxn._slope * (x+y-1) + fxn._intercept for fxn in f.functions()]
    functions.extend(functions2)
    l = zip(polyhedrons, functions)

    #the lower triangle
    for v in f.end_points()[1:len(f.end_points())-1]:
        if f.limit(v,-1) != f.limit(v,1):   
            l.append((Polyhedron(vertices = [[v,0],[0,v]]), f(v)))

    for v in f.end_points()[1:len(f.end_points())-1]:
        if f.limit(v,-1) != f.limit(v,1):   
            l.append((Polyhedron(vertices = [[v,1],[1,v]]), f(v)))
           
    #printPolyFxnPair(l)
    return function_on_complex(l)


def FastPiecewise_to_2d_strip_to_function_on_complex(f):
    """
    Take a FastPiecewise function, convert it to straight 2d, and then to function_on_complex
    Points are extended to lines, intervals are converted to rectangles
    
    Example::
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = gmic()
        sage: F = FastPiecewise_to_2d_strip_to_function_on_complex(h)
        sage: F.minimality_test()
        pi is minimal

        sage: h = dg_2_step_mir_limit()
        sage: F = FastPiecewise_to_2d_strip_to_function_on_complex(h)
        sage: F.minimality_test()
        pi is minimal

        sage: h = drlm_3_slope_limit()
        sage: F = FastPiecewise_to_2d_strip_to_function_on_complex(h)
        sage: F.minimality_test()
        pi is minimal

        sage: h = hildebrand_2_sided_discont_1_slope_1 ()
        sage: F = FastPiecewise_to_2d_strip_to_function_on_complex(h)
        sage: F.minimality_test()
        pi is minimal

        sage: h = not_minimal_1()
        sage: F = FastPiecewise_to_2d_strip_to_function_on_complex(h)
        sage: F.minimality_test()
        Minimality test fails because subadditivity test fails.
    """  
    R.<x,y> = PolynomialRing(QQ)
    polyhedrons = [Polyhedron(vertices = [[intv[0],0],[intv[1],0],[intv[0],1],[intv[1],1]]) for intv in f.intervals()]
    functions = [fxn._slope * x + fxn._intercept + 0 * y for fxn in f.functions()]
    l = zip(polyhedrons, functions)

    for v in f.end_points()[1:len(f.end_points())-1]:
        if f.limit(v,-1) != f.limit(v,1):   
            l.append((Polyhedron(vertices = [[v,0],[v,1]]), f(v)))

    #printPolyFxnPair(l)
    return function_on_complex(l)

def FastPiecewise_to_function_on_complex(f,cont=True):
    """
    Take a FastPiecewise function and convert it to function_on_complex
    
    Example::
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = gmic()
        sage: F = FastPiecewise_to_function_on_complex(h)
        sage: F(1/2)
        5/8
        sage: F.minimality_test(4/5)
        pi is minimal

        sage: h = dg_2_step_mir_limit()
        sage: F = FastPiecewise_to_function_on_complex(h,cont=False)
        sage: F.minimality_test()
        pi is minimal

        sage: h = drlm_3_slope_limit()
        sage: F = FastPiecewise_to_function_on_complex(h)
        sage: F.minimality_test()
        pi is minimal

        sage: h = hildebrand_2_sided_discont_1_slope_1()
        sage: F = FastPiecewise_to_function_on_complex(h)
        sage: F.minimality_test()
        pi is minimal

        sage: h = not_minimal_1()
        sage: F = FastPiecewise_to_function_on_complex(h)
        sage: F.minimality_test()
        Minimality test fails because subadditivity test fails.

        sage: h = ll_strong_fractional()
        sage: F = FastPiecewise_to_function_on_complex(h)
        sage: printPolyFxnPair(F.polyhedronFxnPair)
        [[1]] 0
        [[0], [2/3]] 3/2*x
        [[2/3]] 1
        [[1], [2/3]] 3/2*x - 3/4

    """                                                                                                                       
    
    R.<x> = PolynomialRing(QQ)
    polyhedrons = [Polyhedron(vertices = [[intv[0]],[intv[1]]]) for intv in f.intervals()]
    functions = [fxn._slope * x + fxn._intercept for fxn in f.functions()]
    #print functions
    l = zip(polyhedrons, functions)

    #Add lowest dimension polyhedron(points) where f is not continuous
    for v in f.end_points()[1:len(f.end_points())-1]:
        if f.limit(v,-1) != f.limit(v,1):   
            l.append((Polyhedron(vertices = [[v]]), f(v)))

    #printPolyFxnPair(l)
    return function_on_complex(l)

class function_on_complex :

    def __init__(self, l,cont=True):
        """
        Take a list of pairs (polyhedron, function) and construct them on complex.
        Assume that the polyhedrons don't overlap and cover the whole complex.
        
        Raise error when:
        -The input list doesn't follow the format (polyhedron, function)
        -The ambient dimension of input polyhedron is inconsistent.
        -At some point, the value of pi on the polyhedrons which has the lowest dimesion is not the some.
        
        List of attributes:
        -self.polyhedronFxnPair:       (input) A list of pairs (polyhedron, function)
        -self.polyhedronFxnDict:       A dictonary version of matching polyhedron to functions
        -self.ambientDim:              The ambientDim() of input polyhedrons
        -self.vertices:                A list of vertices of all polyhedrons
        -self._nonnegative:             True if the function pi is nonnegative
        -self.continuous:              True if the function pi is continuous
        
        Example::
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: h = ll_strong_fractional()
            sage: F = FastPiecewise_to_function_on_complex(h)

            sage: R.<x1,x2> = PolynomialRing(QQ)
            sage: f1 = 2*x1+x2+1
            sage: f2 = 2*x1+x2+1
            sage: p1 = Polyhedron(vertices = [[0,0],[0,1],[1,0]])
            sage: p2 = Polyhedron(vertices = [[0,1],[1,0],[1,1]])
            sage: l = [(p1, f1),(p2,f2)]
            sage: F = function_on_complex(l)
            Traceback (most recent call last):
            ...
            ValueError: Error: Value of pi does not match at ... 
        """
        self.polyhedronFxnPair = ((pair[0],pair[1]) for pair in l)
        self.polyhedronFxnPair = tuple(set(self.polyhedronFxnPair))
        self.polyhedronFxnDict = {}
        self.continuous = cont
        ##check input for consistent dimension
        try:
            self.ambientDim = self.polyhedronFxnPair[0][0].ambient_dim()
        except:
            raise ConstructError("Error: The input list must be a list of pairs (polyhedron, function).Example: list = [(polyhedron1, function1), (polyhedron2, function2)]")

        for pair in self.polyhedronFxnPair :
            if pair[0].ambient_dim() != self.ambientDim:
                raise ConstructError('Error: All polyhedrons should be in the same ambient dimension')
            self.polyhedronFxnDict[pair[0]] = pair[1]
        
        ##Note: vector saves time for other operations, but disabled use of set() so that an additional if-statement is required to takes away duplicates.
        ##create the list of all vertices - self.vertices
        self.vertices = []
        for pair in self.polyhedronFxnPair  :
            for v in pair[0].vertices():
                if v.vector() not in self.vertices:
                    self.vertices.extend([v.vector()])     
        #print "vertices:",self.vertices
        
        self._nonnegative = True

        for vertex in self.vertices :
            #Check if every vertex is nonnegative
            if self.get_nonnegative() == true and self.testNeg(vertex) == true :
                self._nonnegative == false
            #Check if the function is continuous and raise construct error if, for any vertex, the value of pi at the lowest dimension polyhedron that contains this vertex don't agree
            if not self.is_continuous_at(vertex) :
                self.continuous = false

        
    def get_nonnegative(self):
        return self._nonnegative


    def testNeg(self,p) :
        """
        For every element(real number) in p, if any of them is negative, return true.
        """
        for number in p :
            if n < 0 :
                return true
        return false

              
    def evaluatePairs(self,polypairs,point) :
        """
        Given a list of pair of (polyhedron,function) and a point p
        Return a list of pair of (polyhedron,value at p)
        Help to build self.verticesLimitPair
        """
        limitValuePair=[]
        for pair in polypairs:
            #print 'pair',pair
            f = pair[1]
            val = self.evaluatePoint(f,point)
            tmp = (pair[0],val)
            limitValuePair.append(tmp)
        return limitValuePair


    def evaluatePoint(self,f,point):
        """
        Given function f and point p
        Return the value of f(p)
        """
        #If the point is a vector, turn it to a tuple
        try :
            p = tuple(point)
        except TypeError:
            p = point
        
        #If f is a constant, return that value
        try :
            val = f(p)
        except TypeError:
            val = f

        return val


    def is_continuous_at(self, point) :
        """
        Translate the given point into the space [0,1]^k
        Example::
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: h = ll_strong_fractional()
            sage: F = FastPiecewise_to_function_on_complex(h)
            sage: F.is_continuous_at([0])
            False
            sage: F.is_continuous_at([2/3])
            False  
            sage: F.is_continuous_at([1/2])
            True
        """
        point = tuple(mod(x) for x in point)
        return self._is_continuous_at(point) 


    def _is_continuous_at(self, point) :
        """
        Given a point, return false if pi is not continuous at this point
        Raise error if the value of pi at the lowest dimension polyhedron at this point don't agree
        """
        l = self.polyhedron_value_pairs_containing(point)
        lowestPolyValuePairs = l[0]
        higherPolyValuePairs = l[1]

        #print 'lowestPolyPairs:'
        #printPolyFxnPair(lowestPolyPairs)
        #print 'higherPolyPairs:'
        #printPolyFxnPair(higherPolyPairs)
        higherlist = list(higherPolyValuePairs)
        higherlist.append(lowestPolyValuePairs[0])
        if len(higherlist) > 0:
            if not self.value_agree_at(higherlist):
                return false
        return true 


    def value_agree_at(self,list):
        """
        Given a list of (polyhedra, value) pairs
        Return true if all values agree 
        """
        if len(list) > 1:
            val = list[0][1]
            for t in list[1:]:
                if val != t[1]:
                    return false
                
        #if self.tmpVertex != [] :
        #    self.verticesValuePair[tuple(self.tmpVertex)] = val
            #self.verticesValuePair.append((self.tmpVertex, val))
        return true 

    
    def polyhedron_value_pairs_containing(self,point):
        """
        Translate the given point into the space [0,1]^k
        Due to periodic, translate 1s to 0s

        Example::
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: h = ll_strong_fractional()
            sage: F = FastPiecewise_to_2d_strip_to_function_on_complex(h)
            sage: l = F.polyhedron_value_pairs_containing([0,0])
            sage: printPolyFxnPair(l[0])
            [[1, 0], [1, 1]] 0
            sage: printPolyFxnPair(l[1])
            [[0, 0], [0, 1], [2/3, 0], [2/3, 1]] 0
            [[1, 0], [1, 1], [2/3, 0], [2/3, 1]] 3/4

        """
        point = tuple(mod(x) for x in point)
        return self._polyhedron_value_pairs_containing(point) 


    @cached_method
    def _polyhedron_value_pairs_containing(self,point):
        """
        Given a point p, 
        Return a list of [lowestPolyValuePairs, higherPolyValuePairs], where lowestPolyValuePairs has lowest dimension polyhedron 
        and higherPolyValuePairs has higher dimension polyhedron.
        Raise error if the value of pi in lowestPolyValuePairs at this point don't agree
        """
        plist = addPeriodicPoints(point)
        #print point,plist
        #Use brute force to find all polyhedrons that contain this point
        l = []   # A list of such polyhedron function pairs
        for pair in self.polyhedronFxnPair :
            for p in plist:
                if p in pair[0] :
                    evaluatedPair = (pair[0],self.evaluatePoint(pair[1],p))
                    l.append(evaluatedPair)

        lset = list(set(l))
        #printPolyFxnPair(lset)

        if len(lset) == 0 :
            raise RangeError("Error: The given point {} is not in any polyhedron.".format(point))

        ##get a list of polyhedron pairs of current lowest dimensions (lowestPolyPairs)
        lowestPolyValuePairs = [lset[0]]
        higherPolyValuePairs = []
        if len(lset) > 1:
            for pair in lset[1:] :
                if pair[0].dim() < lowestPolyValuePairs[0][0].dim():
                    higherPolyValuePairs.extend(lowestPolyValuePairs)
                    lowestPolyValuePairs = [pair]
                elif pair[0].dim() == lowestPolyValuePairs[0][0].dim():
                    lowestPolyValuePairs.append(pair)
                else:
                    higherPolyValuePairs.append(pair)

        if not self.value_agree_at(lowestPolyValuePairs) :
            #print 'Problem:',point
            #for v in lowestPolyValuePairs:
            #    print v[0].vertices_list(),v[1]
            raise ValueError('Error: Value of pi does not match at {}'.format(point))


        return (tuple(lowestPolyValuePairs), tuple(higherPolyValuePairs))


    def __call__(self, point):
        """
        Return the value of the function corresponding to the given point(tuple).

        Example:
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: h = ll_strong_fractional()
            sage: F = FastPiecewise_to_function_on_complex(h)
            sage: F([1/2])
            3/4
            sage: F([2/3])
            1
        """
        point = tuple([mod(x) for x in point])

        list = self.polyhedron_value_pairs_containing(point)
        lowestPolyValuePairs = list[0]
        f = lowestPolyValuePairs[0][1]
        return self.evaluatePoint(f,point)


    def limits(self, point):
        """
        Translate the given point into the space [0,1]^k
        """
        point = tuple([mod(x) for x in point])
        return self._limits(point)         

    
    def _limits(self, point):
        """
        Given a point (list)
        Return a list of all (polyhedra, value) pairs where the polyhedra contains this point
        Example::
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: h = ll_strong_fractional()
            sage: F = FastPiecewise_to_2d_strip_to_function_on_complex(h)
            sage: l = F.limits([0,0])
            sage: printPolyFxnPair(l)
            [[1, 0], [1, 1]] 0
            [[0, 0], [0, 1], [2/3, 0], [2/3, 1]] 0
            [[1, 0], [1, 1], [2/3, 0], [2/3, 1]] 3/4

        """
        l = self.polyhedron_value_pairs_containing(point)
        lowestPolyValuePairs = list(l[0])
        #print "low:", lowestPolyValuePairs
        higherPolyValuePairs = list(l[1]) 
        #print "high:",higherPolyValuePairs

        lowestPolyValuePairs.extend(higherPolyValuePairs)
        #print "all:",lowestPolyValuePairs
        return lowestPolyValuePairs


    def limit(self, point,T):
        """
        Translate the given point into the space [0,1]^k
        """
        point = tuple(mod(x) for x in point)
        return self._limit(point,T)  


    @cached_method
    def _limit(self,point,T): 
        """
        Return the limit value of point in the direction of a polyhedron T
        T is assumed to be contained in one of the polyhedron that construct this function_on_complex
        Raise error if:
            -point is not in T
            -T is not contained in one of the polyhedron that construct this function_on_complex

        Example::
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: h = ll_strong_fractional()
            sage: F = FastPiecewise_to_function_on_complex(h)
            sage: p1 = Polyhedron(vertices = [[1/2],[2/3]])
            sage: F.limit([2/3],p1)
            1
            sage: p2 = Polyhedron(vertices = [[2/3]])
            sage: F.limit([2/3],p2)
            1
            sage: p3 = Polyhedron(vertices = [[2/3],[4/5]])
            sage: F.limit([2/3],p3)
            1/4

            sage: p = Polyhedron(vertices = [[0],[1/2]])
            sage: F.limit([2/3],p)
            Traceback (most recent call last):
            ...
            ConstructError: 'The limit is undefined because the point ... is not in the polyhedra [[0], [1/2]]'

            sage: p = Polyhedron(vertices = [[1/2],[4/5]])
            sage: F.limit([2/3],p)
            Traceback (most recent call last):
            ...
            NotImplementedError: The directional polyhedra [[1/2], [4/5]] is not inside any polyhedron that construct this function_on_complex instance
        """
        #Check if point is in T
        plist = addPeriodicPoints(point)
        inside = False
        for pt in plist:
            if inside: 
                break
            elif pt in T:
                inside = True
        if not inside:
            raise ConstructError('The limit is undefined because the point {} is not in the polyhedra {}'.format(point, T.vertices_list()))


        list = self.polyhedron_value_pairs_containing(point)
        lowestPolyValuePairs = list[0]
        #print 'low:'
        #printPolyFxnPair(lowestPolyValuePairs)
        higherPolyValuePairs = list[1] 
        #print 'high'
        #printPolyFxnPair(higherPolyValuePairs)

        if self.is_continuous_at(point):  
            return lowestPolyValuePairs[0][1]  #continuous imples that all values are the same

        if T.dim() == 0:          #T is a point
            return lowestPolyValuePairs[0][1]

        for pair in lowestPolyValuePairs:
            if poly_isInsidePoly(T,pair[0]):
                return pair[1]

        #Find the lowest dimension polyhedra from all polyhedron that contains the point
        lowestDimPair = None
        lowestDim = higherPolyValuePairs[0][0].dim()
        for pair in higherPolyValuePairs:
            #print pair[0].vertices_list(), pair[0].dim(), lowestDim
            if poly_isInsidePoly(T,pair[0]) and pair[0].dim() <= lowestDim:
                lowestDim = pair[0].dim()
                lowestDimPair = pair
        
        if lowestDimPair != None:
            return lowestDimPair[1]
        else:
            raise NotImplementedError, "The directional polyhedra {} is not inside any polyhedron that construct this function_on_complex instance".format(T.vertices_list())


    def polyhedron_containing_line(self,line):
        """
        Given a line (as Polyhedron instance), assume the line is in the PolyhedronFxnPair.
        Return the list of polyhedron in PolyhedronFxnPair that contains this line.
        """
        vert = line.vertices_list()
        #find list of all 2d polyhedrons that contains this line(two points)

        polyVal1List = self.polyhedron_value_pairs_containing(vert[0])
        if len(polyVal1List[1]) == 0:
            poly1List = [p[0] for p in polyVal1List[0]]
        else:
            poly1List = [p[0] for p in polyVal1List[1]]

        polyVal2List = self.polyhedron_value_pairs_containing(vert[1])
        if len(polyVal2List[1]) == 0:
            poly2List = [p[0] for p in polyVal2List[0]]
        else:
            poly2List = [p[0] for p in polyVal2List[1]]   
        
        commom_polyhedra = set(poly1List).intersection(set(poly2List))

        #get rid of periodic
        commom_polyhedra_real = []
        for polyhedra in commom_polyhedra:
            if vert[1] in polyhedra:
                commom_polyhedra_real.append(polyhedra)

        return commom_polyhedra_real


    def plot_2d_complex(self, discontWidth = 1/60):
        """
        Plot the 2d complex
        DiscontWidth is the width of the white strip which shows discontinuous in the graph.
        Example::
            sage: h = hildebrand_discont_3_slope_1()
            sage: F = FastPiecewise_to_2d_diagonal_to_function_on_complex(h)
            sage: F.plot_2d_complex()
            Graphics object consisting of ... graphics primitives

            sage: L = F.polyhedronFxnPair
            sage: p2 = Polyhedron(vertices = [[1/2,1/2]])
            sage: L2 = list(L)
            sage: L2.append([p2,1/2])
            sage: F2 = function_on_complex(L2)
            sage: F2.plot_2d_complex()          
            Graphics object consisting of ... graphics primitives

            sage: p3 = Polyhedron(vertices = [[1/32,1/32]])
            sage: L3 = list(L)
            sage: L3.append([p3,1/8])
            sage: F3 = function_on_complex(L3)
            sage: F3.plot_2d_complex()          
            Graphics object consisting of ... graphics primitives            

        """
        if self.ambientDim != 2:
            raise NotImplementedError, "Currently there's no plot method except 2 dimensional complex"

        #create a dict. key: gradient vector, value: index 
        self.gradientDic = {}
        index = 0
        for pair in self.polyhedronFxnPair:
            if isinstance(pair[1],Rational):
                g = tuple([0] * self.ambientDim)
            else:
                g = tuple((pair[1].gradient()))

            if g not in self.gradientDic:
                self.gradientDic[g] = index
                index += 1

        #colors list
        self.colors = rainbow(len(self.gradientDic))
        #print self.gradientDic

        image = Graphics()

        #Get all ticks on axis that the graph needs to show. 
        t1 = {v[0] for v in self.vertices}  #x-axis
        t2 = {v[1] for v in self.vertices}  #y-axis 
        self.xtick_formatter = [ "$%s$" % latex(x) for x in t1 ]
        self.ytick_formatter = [ "$%s$" % latex(y) for y in t2 ]

        self.labeled = []   #labeled gradients on the graph
        #the complex, plot part of the graph by self.polyhedronFxnPair[0:n], for some number n
        for pair in self.polyhedronFxnPair:
            #print pair[0].vertices_list()
            #2d polyhedrons
            if pair[0].dim() == 2:
                plots = self.add_2d_polyhedron_plots(pair,discontWidth,t1,t2)
            #1d polyhedrons(lines)
            elif pair[0].dim() == 1:
                plots = self.add_1d_lines_plots(pair)
            #0d polyhedrons(points)
            elif not self._is_continuous_at(pair[0].vertices_list()[0]):
                plots = self.add_0d_points_plots(pair,discontWidth)
            for element in plots:
                image += element

        return image

    def add_2d_polyhedron_plots(self,pair,discontWidth,t1,t2):
        l = []
        g = tuple(pair[1].gradient())
        color = self.colors[self.gradientDic[g]]
        draw = self.calculateWhiteStrip(pair,discontWidth)
        shrinkedPolyhedron = draw[0]

        #add polyhedrons
        if g not in set(self.labeled):
            l.append(plot(shrinkedPolyhedron,wireframe=False, fill=color, alpha=0.5,ticks=[list(t1),list(t2)],tick_formatter = [self.xtick_formatter, self.ytick_formatter],zorder = 1, legend_label = "Gradient %s" % str(g)))
            self.labeled.append(g)
        else:
            l.append(plot(shrinkedPolyhedron,wireframe=False, fill=color, alpha=0.5,zorder = 1))

        #add sectors
        for s in draw[1]:
            l.append(s)  
        return l    

    def add_1d_lines_plots(self,pair):
        l = []
        vert = pair[0].vertices_list()
        linecolor = []
        neighborPoly = self.polyhedron_containing_line(pair[0])

        for polyhedra in neighborPoly:
            if polyhedra.dim() == 1:
                continue
            f = self.polyhedronFxnDict[polyhedra]
            if self.value_matches_at_line(vert[0],vert[1],f,pair[1]):
                linecolor.append(self.colors[self.gradientDic[tuple(f.gradient())]])

        if len(linecolor) == 0:
            linecolor.append('black')
        if len(linecolor) < 2:
            l.append(line2d(pair[0].vertices_list(),rgbcolor=linecolor[0],thickness=2, zorder = 2))
        #print linecolor

        return l            

    def add_0d_points_plots(self,pair,discontWidth):
        l = []
        point = pair[0].vertices_list()[0]
        #find the list of polyhedron that contains the point
        polyValList = self.polyhedron_value_pairs_containing(point)[1]
        polyValList = [polyVal for polyVal in polyValList if point in polyVal[0]]
        point_value = self.evaluatePoint(pair[1],point)
        for polyVal in polyValList:
            #if the point is discontinuous for the 2d polyehedra
            if polyVal[0].dim() == 2 and self.evaluatePoint(polyVal[1],point) != point_value:
                #the case where the point is on the vertices has been done 
                if list(point) not in polyVal[0].vertices_list():
                    #if the point is on the edge of the polyhedra, add a semicircle; else a circle.
                    semicircle = self.calculateSemicircle(point,polyVal[0],discontWidth)
                    if not semicircle:
                        calculatedWidth = min(stripWidth_upper_bound_given_point(point,polyVal[0]),discontWidth)
                        l.append(circle(point,calculatedWidth,rgbcolor = 'white',fill = True, facecolor = 'white',zorder = 2))
                    else:
                        l.append(semicircle)

        l.append(point2d(pair[0].vertices_list(),rgbcolor='black',size = 10,zorder=3))
        return l        

    def calculateWhiteStrip(self,pair,width):
        """
        Given a pair of polyhedra, function
        Return a list of items to draw.
        The first one is a shrinked polyhedra which is equivalent to:
         adding white strips to the edges where the original polyhedron is not continuous
        Others are sectors at the discontinuous vertices of the original polyhedra
        """
        if pair[0].dim() < 2:
            raise NotImplementedError
        edgeIneqDict = {}   #key: edge : (vertex1, vertex2), value: inequality
        vertices = pair[0].vertices_list()
        #print "polyhedra", vertices
        inequalities = pair[0].inequalities_list()

        #construct the edgeIneqDict by brute force
        for ineq in pair[0].inequality_generator():
            edge = []
            for v in vertices:
                sum = 0
                for i in range(self.ambientDim):
                    sum += ineq.A()[i] * v[i]
                if sum == ineq.b()*(-1):
                    edge.append(tuple(v))
            edgeIneqDict[tuple(edge)] = ineq
        #print "edgeIneqDict",edgeIneqDict

        #construct the shrinked polyhedra and vertEdgeDict
        vertEdgeDict = {}   #key: vertex, value: (edge1, edge2)
        for key, value in edgeIneqDict.iteritems():
            #check if the polyhedra is generally continuous at this edge
            if not self.is_continuous_at_line(key[0],key[1],pair[1]):
                #make sure that the width is not too large
                calculatedWidth = min(stripWidth_upper_bound_given_edge(key[0],key[1],pair[0]),width)
                #normalize v
                v = value.vector() / RR(value.vector().norm(2))
                
                if v[2] == 0 and v[1] > 0: #edge is vertical, goes right
                    stripVector = [v[0]-calculatedWidth*v[1],v[1],0]
                elif v[2] == 0 and v[1] < 0: #edge is vertical, goes left
                    stripVector = [v[0]+calculatedWidth*v[1],v[1],0]
                elif v[2] < 0:  #edge goes down
                    c = shiftDistance(key[0],key[1],calculatedWidth)
                    stripVector = [v[0]+c*v[2],v[1],v[2]]
                else:  #edge goes up
                    c = shiftDistance(key[0],key[1],calculatedWidth)
                    stripVector = [v[0]-c*v[2],v[1],v[2]]

                #print "stripVector",stripVector
                inequalities.append(stripVector)

            #add key[0] to vertEdgeDict
            if not key[0] in vertEdgeDict:
                vertEdgeDict[key[0]] = [key[1]]
            else:
                edges = vertEdgeDict[key[0]]
                if len(edges) == 1 and edges[0] != key[1]:
                    edges.append(key[1])
                    vertEdgeDict[key[0]] = tuple(edges)
            #add key[1] to vertEdgeDict
            if not key[1] in vertEdgeDict:
                vertEdgeDict[key[1]] = [key[0]]
            else:
                edges = vertEdgeDict[key[1]]
                if len(edges) == 1 and edges[0] != key[0]:
                    edges.append(key[0])
                    vertEdgeDict[key[1]] = tuple(edges)

        #print vertEdgeDict

        #add sectors at discontinuous points
        sectorL = self.addSectorsInPolyhedron(vertEdgeDict,pair,width)
        return (Polyhedron(ieqs = inequalities),sectorL)


    def calculateSemicircle(self,point,polyhedra,width):
        """
        If the point is on any edge of the polyhedra, return the semicircle
        Else return false
        """
        edgeList = [e for e in polyhedra.bounded_edges()]
        for e in edgeList:
            edge1 = tuple(e[0])
            edge2 = tuple(e[1])
            if is_on_the_line(edge1,edge2,point):
                dist1 = sqrt((edge1[1]-point[1])^2+(edge1[0]-point[0])^2)
                dist2 = sqrt((edge2[1]-point[1])^2+(edge2[0]-point[0])^2)
                calculatedWidth = min(stripWidth_upper_bound_given_edge(edge1,edge2,polyhedra),width,dist1,dist2)
                return self.calculateSectorInPolyhedra(point,edge1,edge2,calculatedWidth,polyhedra)

        return False


    def is_continuous_at_line(self,v1,v2,f):
        """
        Given a line [v1,v2] and a function f,
        Return true if the line is generally continuous(besides some discontinous points on the line)
        """
        if self(v1) == self.evaluatePoint(f,v1) and self(v2) == self.evaluatePoint(f,v2):
            return True
        test1 = random.randint(1,99) / 100
        p1 = point_on_line_from_percentage(test1,v1,v2)
        test2 = random.randint(1,99) / 100
        p2 = point_on_line_from_percentage(test2,v1,v2)
        test3 = random.randint(1,99) / 100
        p3 = point_on_line_from_percentage(test3,v1,v2)

        if self(p1) == self.evaluatePoint(f,p1) or self(p2) == self.evaluatePoint(f,p2) or self(p3) == self.evaluatePoint(f,p3):   
            return True
        return False


    def value_matches_at_line(self,v1,v2,f1,f2):
        """
        Given a line [v1,v2] and two functions f,
        Return true if the line has the same value under f1 and f2
        """
        if self.evaluatePoint(f1,v1) != self.evaluatePoint(f2,v1) or self.evaluatePoint(f1,v2) != self.evaluatePoint(f2,v2):
            return False
        test1 = random.randint(1,99) / 100
        p1 = point_on_line_from_percentage(test1,v1,v2)
        test2 = random.randint(1,99) / 100
        p2 = point_on_line_from_percentage(test2,v1,v2)   

        if self.evaluatePoint(f1,p1) != self.evaluatePoint(f2,p1) or self.evaluatePoint(f1,p2) != self.evaluatePoint(f2,p2):   
            return False

        return True


    def addSectorsInPolyhedron(self,vertEdgeDict,polyFxnPair,width):
        """
        Return a list of sectors (as disk instances) with given radius, which fits exactly 
        to the corners of the polyhedra, for all discontinuous vertices of the polyhedra.
        """
        sectorL = []
        for key, value in vertEdgeDict.iteritems():
            #print key,value
            if self(key) != self.evaluatePoint(polyFxnPair[1],key):
                edge1 = value[0]
                edge2 = value[1]

                width1 = stripWidth_upper_bound_given_edge(key,edge1,polyFxnPair[0])
                width2 = stripWidth_upper_bound_given_edge(key,edge2,polyFxnPair[0])
                calculatedWidth = min([width1,width2,width])

                s = self.calculateSectorInPolyhedra(key,edge1,edge2,calculatedWidth,polyFxnPair[0])
                sectorL.append(s)
        return sectorL


    def calculateSectorInPolyhedra(self,orig,edge1,edge2,radius,polyhedra):
        """
        Return a sector (as a disk instance) with given radius, where its origin is a given vertex of a polyhedra,
        and its two edges are the same edges of the polyhedra from the vertex.
        """
        vector1 = ((edge1[0]-orig[0]),(edge1[1]-orig[1]))
        vector2 = ((edge2[0]-orig[0]),(edge2[1]-orig[1]))
                
        angle1 = arctan2(vector1[1],vector1[0])
        angle2 = arctan2(vector2[1],vector2[0])
        if angle1 < 0:
            positive_angle1 = 2*pi + angle1
        else:
            positive_angle1 = angle1
        if angle2 < 0:
            positive_angle2 = 2*pi + angle2
        else:
            positive_angle2 = angle2
                
        mid_angle = (positive_angle1+positive_angle2)/2
        testPoint = findPointFrom(orig,mid_angle,radius)

        if testPoint in polyhedra:
            if mid_angle <= positive_angle1:
                s = disk(orig,radius, (angle2,angle1),color='white',zorder=2)
            else:
                s = disk(orig,radius, (angle1,angle2),color='white',zorder=2)
        else:
            if mid_angle <= positive_angle1:
                s = disk(orig,radius, (angle1,angle2),color='white',zorder=2)
            else:
                s = disk(orig,radius, (angle2,angle1),color='white',zorder=2)    
        return s               


    def minimality_test(self,f=None) :
        """
        Check if pi is minimal on the complex

        List of attributes:
        -self.delta_p_vertices:             A list of vertices of delta p
        -self.delta_p_list:                 A list of polyhedrons (FIJK) that constructs delta p
        -self.pi_on_delta_p:                A function_on_complex instance

        Example::
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: h = ll_strong_fractional()
            sage: F = FastPiecewise_to_function_on_complex(h)
            sage: F.minimality_test()
            pi is minimal
            sage: F = FastPiecewise_to_2d_strip_to_function_on_complex(h)
            sage: F.minimality_test()
            pi is minimal
            sage: F = FastPiecewise_to_2d_diagonal_to_function_on_complex(h)
            sage: F.minimality_test()
            pi is minimal

        """
        if self.get_nonnegative() == false :
            print "Minimality test fails because the function pi is negative"
            return

        ##Find f
        if f == None :
            f = self.find_f()
        else :
            if self(f) != 1 :
                print "Minimality test fails because symmetry test fails: pi(f) is not equal to 1."
                return

        if f == None :
            print "Minimality test fails because there's no point that the function pi evaluates to 1."
            return
        self.f = f
    
        ##check pi(0) = 0
        zeroVector = [0] * self.ambientDim
        if self((zeroVector)) != 0 :
            #print zeroVector, self((zeroVector))
            print "Minimality test fails because pi(0) is not equal to 0."
            return

        ##construct a delta p and perform subadditivity and symmetry test
        self.construct_delta_p()
        if self.subadditivity_test() == false:
            print "Minimality test fails because subadditivity test fails."
            return          
        if self.symmetry_test() == false:
            print "Minimality test fails because symmertry test fails."
            return
        print "pi is minimal"


    def symmetry_test(self) :
        for a in self.delta_p_vertices :
            if self.modUVF(a) and self.delta_pi(a) != 0 :                
                return false
        return true     

        
    def modUVF(self,a):
        """
        Return true if u+v-f is an integer
        """
        u = a[:(len(a)/2)]
        v = a[(len(a)/2):]
        for i in range (len(u)) :
            if (u[i] + v[i] - self.f[i]) not in ZZ :
                return false
        return true


    def subadditivity_test(self) :
        if self.continuous:
            for a in self.delta_p_vertices :
                if self.delta_pi(a) < 0 :
                    return false
            return true
        else:
            comb2 = combinationTwo(len(self.delta_p_list))
            for comb in comb2:
                E = self.delta_p_list[comb[0]]
                F = self.delta_p_list[comb[1]]
                #F contains E
                if poly_isStrictlyInsidePoly(E.polyhedron,F.polyhedron):
                    for a in E.polyhedron.vertices_list():
                        if self.delta_pi_with_limit(a,F) < 0 :
                            print 'u:',a[:(len(a)/2)], 'v:', a[(len(a)/2):]
                            print 'p1:', F.p1.vertices_list(), 'p2:',F.p2.vertices_list(),'p3:',F.p3.vertices_list()
                            print E.polyhedron.vertices_list(), "is in", F.polyhedron.vertices_list()
                            return false
            return true


    def find_f(self) :
        """
        Find a vertex that its value is 1. Return 0 if no such vertex exists
        """
        for v in self.vertices:
            if self(v) == 1:
                return v
        return None
        

    def delta_pi(self,a) :
        u = a[:(len(a)/2)]
        v = a[(len(a)/2):]
        sum = [ u[i] + v[i] for i in range(len(u))]
        return self(u) + self(v) - self(sum)


    def delta_pi_with_limit(self,a,F) :
        u = a[:(len(a)/2)]
        v = a[(len(a)/2):]
        sum = [ u[i] + v[i] for i in range(len(u))]
        #if (self.limit(u,F.I_pi)+ self.limit(v,F.J_pi) - self.limit(sum,F.K_pi)) < 0:
            #print u, v,sum
            #print F.I.vertices_list(), F.J.vertices_list(),F.K.vertices_list()
            #print self.limit(u,F.I_pi), self.limit(v,F.J_pi), self.limit(sum,F.K_pi)
        return self.limit(u,F.p1) + self.limit(v,F.p2) - self.limit(sum,F.p3)


    def construct_delta_pi_on_delta_p(self):
        """
        Return a function_on_complex instance of function pi on complex delta_p
        Attention: Seems to work only for 1d cases. 
        Bug(value don't match) occurs for 2d cases, such as 2d diagonal and 2d strip.

        Example::
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: h = ll_strong_fractional()
            sage: F = FastPiecewise_to_function_on_complex(h)
            sage: F2 = F.construct_delta_pi_on_delta_p()
            sage: F2.get_nonnegative()
            True

            sage: h = hildebrand_discont_3_slope_1()
            sage: F3 = FastPiecewise_to_function_on_complex(h)
            sage: F4 = F.construct_delta_pi_on_delta_p()
            sage: F5 = FastPiecewise_to_2d_strip_to_function_on_complex(h)
            sage: F6 = F5.construct_delta_pi_on_delta_p()

        """
        l = len(self.polyhedronFxnPair)
        comb2 = list(product(range(l),repeat = 2))
        l_length_list = [i for i in range(l*2)]
        F_I_J_K_comb = list(product(comb2,l_length_list))
        
        delta_p_function_pair = []
        d = self.ambientDim
        #create variable name for ring
        R = PolynomialRing(QQ,d,var_array = ['u','v'])
        #get variables
        var = R.gens()
        translate = False

        for t in F_I_J_K_comb:
            #print t
            if t[1] >= l:
                translate = True  #k is in [0,1]^k
            else:
                translate = False  #k is in [1,2]^k

            F_I_J_K = self.createFIJK(t[0][0],t[0][1],t[1])
            vertList = F_I_J_K.polyhedron.vertices_list()
            if vertList:
                #print 'FIJK:',vertList
                pi_I = self.findFxn(F_I_J_K.p1)
                #print "I",F_I_J_K.I.vertices_list(), "p1", F_I_J_K.p1.vertices_list()
                #print "before:",pi_I
                if not isinstance(pi_I,Rational): 
                    pi_I = pi_I([var[2*i] for i in range(d)])
                #print "after:",pi_I


                pi_J = self.findFxn(F_I_J_K.p2)
                #print "J",F_I_J_K.J.vertices_list(), "p2", F_I_J_K.p2.vertices_list()
                #print "before:",pi_J
                if not isinstance(pi_J,Rational): 
                    pi_J = pi_J([var[2*i+1] for i in range(d)])
                #print "after:",pi_J


                if translate:
                    pi_K = self.find_translated_fxn(F_I_J_K.p3)
                else:
                    pi_K = self.findFxn(F_I_J_K.p3)
                #print "K",F_I_J_K.K.vertices_list(), "p3", F_I_J_K.p3.vertices_list()
                #print "before:",pi_K
                if not isinstance(pi_K,Rational): 
                    pi_K = pi_K([var[2*i]+var[2*i+1] for i in range(d)])
                #print "after:",pi_K

                pi = pi_I + pi_J - pi_K
                #print "pi",pi,'\n'
                pair = (F_I_J_K.polyhedron,pi)
                delta_p_function_pair.append(pair)

        s = list(set(delta_p_function_pair))

        #printPolyFxnPair(delta_p_function_pair)
        return function_on_complex(s)

    def find_translated_fxn(self, p):
        """
        Translate the given polyhedra into the space [0,1]^k
        """
        l = p.vertices_list()
        
        translated_vertices = []
        for v in l:
            vertex = [mod(p) for p in v]
            translated_vertices.append(vertex)
        return self.findFxn(Polyhedron(vertices = translated_vertices))  


    @cached_method
    def findFxn(self,p):
        """
        Given polyhedra p (which is contained in some polyhedronFxnPair)
        Return the function on the lowest dimension polyhedra that contains it.
        """
        l = list(set(pair[0] for pair in self.polyhedronFxnPair))
        periodicPolyList = self.addPeriodicPolyhedron(p)
        polyContainDict = {} #key: input p, value: the lowest dim polyhedra that contains it.

        for poly in periodicPolyList:
            for original_polyhedron in l:
                #print "test:", original_polyhedron.vertices_list()
                if poly_isInsidePoly(poly,original_polyhedron):
                    if not poly in polyContainDict.keys() or polyContainDict[poly].dim() > original_polyhedron.dim():
                        polyContainDict[poly] = original_polyhedron

        if not polyContainDict:
            raise ValueError('Could not find a polyhedra from the polyhedronFxnPair which contains the input polyhedra {}'.format(p.vertices_list()))

        #find the answer with the lowest dimension
        for key, value in polyContainDict.iteritems(): #get first in the dict
            lowestDim = value.dim()
            lowestDimPoly = key
            break
      
        for key, value in polyContainDict.iteritems():
            if value.dim() < lowestDim:
                lowestDimPoly = key
            #print key.vertices_list(), value.vertices_list()

        #find the function and modify it if we are using a periodic polyhedron different from p
        original_fxn = self.polyhedronFxnDict[polyContainDict[lowestDimPoly]]

        #find the different index & difference
        for i in range(self.ambientDim):
            diff = lowestDimPoly.vertices_list()[0][i] - p.vertices_list()[0][i]
            if diff != 0:
                break
        different_index = i


        d = self.ambientDim
        R = PolynomialRing(QQ,d,var_array = ['x'])
        #get variables
        var = R.gens()

        modified_fxn = original_fxn
        if not isinstance(original_fxn,Rational):
            modified_fxn = original_fxn([var[i]+diff if i == different_index else var[i] for i in range(d)] )

        #print p.vertices_list(), "is in", polyContainDict[lowestDimPoly].vertices_list()
        #print 'original:',original_fxn
        #print 'modified:',modified_fxn
        return modified_fxn

    
    def addPeriodicPolyhedron(self,p):
        """
        Given a polyhedra p, if p is on the boundary, return a list of its periodic polyhedra(including p)
        Else return [p]

        Example::
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: p = Polyhedron(vertices = [[0,0,0]])
            sage: F1 = function_on_complex([(p,0)])
            sage: p1 = Polyhedron(vertices = [[0,0,0],[0,0,1/2],[1/2,0,0],[1/2,0,1/2]])
            sage: printPoly(F1.addPeriodicPolyhedron(p1))
            [[0, 0, 0], [0, 0, 1/2], [1/2, 0, 0], [1/2, 0, 1/2]]
            [[0, 1, 0], [0, 1, 1/2], [1/2, 1, 0], [1/2, 1, 1/2]]

            sage: p1 = Polyhedron(vertices = [[0,0,0],[1,0,0]])
            sage: printPoly(F1.addPeriodicPolyhedron(p1))
            [[0, 0, 0], [1, 0, 0]]
            [[0, 1, 0], [1, 1, 0]]
            [[0, 0, 1], [1, 0, 1]]
            [[0, 1, 1], [1, 1, 1]]
        """        
        if p.ambient_dim() != self.ambientDim:
            raise ConstructError('The ambient dimension of the input polyhedra should be'.format(self.ambientDim))
        if p.dim() > self.ambientDim:
            #print p.vertices_list()
            raise ConstructError('The dimension of the input polyhedra should be less than or equal to the ambient dimension of the polyhedrons: {}'.format(self.ambientDim))

        if not self.is_on_boundary(p):
            return [p]

        if p.dim() == 0:
            l = addPeriodicPoints(p.vertices_list()[0])
            return [Polyhedron(vertices = [vertList]) for vertList in l]

        sameVector = [1]*self.ambientDim  #1 if all vertices of p is the same on this index
        for i in range(len(p.vertices_list())):
            if sameVector == [0]*self.ambientDim:
                break
            for j in range(i+1,len(p.vertices_list())):
                for k in range(self.ambientDim):
                    if p.vertices_list()[i][k] != p.vertices_list()[j][k]:
                        sameVector[k] = 0
        #print sameVector

        fixedVector = []  #get the part of the original vector that is fixed
        for i in range(self.ambientDim):
            if sameVector[i] == 1:
                fixedVector.append(p.vertices_list()[0][i])
        #print fixedVector

        l = addPeriodicPoints(fixedVector)
        #print l

        polyhedronList = []
        for periodicVector in l:
            vertList = []
            for original_vertex in p.vertices_list():
                changedVertex = original_vertex
                index = 0  #keep record of the working index
                for i in range(len(sameVector)):
                    if sameVector[i] == 1:
                        changedVertex[i] = periodicVector[index]
                        index += 1
                vertList.append(changedVertex)
            #print vertList
            polyhedronList.append(Polyhedron(vertices = vertList))

        return polyhedronList        


    def is_on_boundary(self,p):
        """
        Given a polyhedra p
        Return True if p is on the boundary of [0,1]^k
        """
        bList = facets(self.ambientDim)

        for poly in bList:
            if poly_isInsidePoly(p,poly):
                return True
        return False


    def construct_delta_p(self):

        l = len(self.polyhedronFxnPair)
        comb3 = self.combinationThree(l)

        self.delta_p_vertices = []
        self.delta_p_list = []
        for t in comb3 :
            F_I_J_K = self.createFIJK(t[0],t[1],t[2])
            vertList = F_I_J_K.polyhedron.vertices_list()
            if vertList:
                self.delta_p_vertices.extend(vertList)
                self.delta_p_list.append(F_I_J_K)

        #print "delta_p:"
        #for f in self.delta_p:
        #    print f.polyhedron.vertices_list()
        #Note! delta_p has duplicates!
        #print 'delta_p',self.delta_p

    def createFIJK(self, a, b, c) :
        """
        Given integers a,b,c from combinationThree
        Find the corresponding polyhedronFxnPair according to the index a,b,c
        """
        A = self.polyhedronFxnPair[a]
        B = self.polyhedronFxnPair[b]
        l = len(self.polyhedronFxnPair)
        #k in [0,1]^k
        if c < l:
            C = self.polyhedronFxnPair[c]
            return FIJK(A,B,C)

        #k in [1,2]^k
        C = self.polyhedronFxnPair[c-l]
        vertList = C[0].vertices_list()
        
        translated_vertices = []
        for v in vertList:
            vertex = [p+1 for p in v]
            translated_vertices.append(vertex)
        translated_pair = (Polyhedron(vertices = translated_vertices),C[1])
        return FIJK(A,B,translated_pair)


    def combinationThree(self,n):
        """
        Given a list of n different integers from 0 to n-1
        Return a list of 3-tuples, which are all possible combinations of these integers.
        Note: For use of FIJK, since I and J is interchangable, if there's [1 2 3], then [2 1 3] is not included
        """
        comb3 = []
        for i in range (n):
            a = []
            a.append(i)
            for j in range (i, n):
                if len(a) >= 2:
                    a = a[0:1]
                a.append(j)
                for k in range (n) :
                    if len(a) >= 3 :
                        a = a[0:2]
                    a.append(k)
                    comb3.append(a)
        #print 'comb:',comb
        return comb3

def count_0_and_1(l):
    """
    Given list l, count the number of 0s and 1s in the list
    """
    count = 0
    for i in l:
        if i == 1 or i == 0:
            count += 1

    return count

def stripWidth_upper_bound_given_edge(v1,v2,polyhedra):
    """
    Given a line fro
    m v1 to v2, which is an edge of the polyhedra
    calculate the shortest distance between the line and all vertices of polyhedra
    Return 1/3 of such distance.
    """
    vlist = polyhedra.vertices_list()
    disList = []

    for v in vlist:
        dis = distance_from_point_to_line(v,v1,v2)
        if dis != 0:
            disList.append(dis)
    return RR(min(disList)/3)

def stripWidth_upper_bound_given_point(p,polyhedra):
    """
    Given a point p which is not on any edge of the polyhedra
    calculate the shortest distance between the point and all vertices of polyhedra
    Return 1/3 of such distance.
    """
    edgeList = [e for e in polyhedra.bounded_edges()]
    disList = [distance_from_point_to_line(p,tuple(e[0]),tuple(e[1])) for e in edgeList]

    return RR(min(disList)/3)

def distance_from_point_to_line(p,v1,v2):
    """
    Given line (v1,v2)
    Return the distance from point p to this line
    """
    dis_numerator = abs((v2[0]-v1[0])*(v1[1]-p[1])-(v1[0]-p[0])*(v2[1]-v1[1]))
    dis_denomiator = sqrt((v2[0]-v1[0])^2+(v2[1]-v1[1])^2)
    if dis_denomiator == 0 or dis_numerator == 0:
        return 0;

    return RR(dis_numerator/dis_denomiator)

def point_on_line_from_percentage(p,v1,v2):
    """
    Given line (v1,v2), 0<p<1
    Return a point that is p from v1.
    """
    if v1[0] < v2[0]:
        x = v1[0] + p * abs(v1[0] - v2[0])
    else:
        x = v1[0] - p * abs(v1[0] - v2[0])

    if v1[1] < v2[1]:
        y = v1[1] + p * abs(v1[1] - v2[1])
    else:
        y = v1[1] - p * abs(v1[1] - v2[1])
    return [x,y]

def is_on_the_line(p1,p2,p):
    """
    Given 3 points,
    Return true if p is on the line between p1 and p2
    """
    #check if p is in the middle
    if p[1] > max(p1[1],p2[1]) or p[1] < min(p1[1],p2[1]):
        return false
    if p[0] > max(p1[0],p2[0]) or p[0] < min(p1[0],p2[0]):
        return false        

    #if the line is vertical
    if p1[0] == p2[0]:
        return p == p1[0]

    slope1 = (p1[1]-p[1])/(p1[0]-p[0])
    slope2 = (p[1]-p2[1])/(p[0]-p2[0])
    return slope1==slope2

def shiftDistance(edge1,edge2,width):
    """
    Given an edge which vertices are edge1 and edge2
    Find a line parallel to this edge such that the shortest distance between them is width
    Return their distance on y-axis
    """
    v = (edge1[0]-edge2[0],edge1[1]-edge2[1])
    
    phi1 = arctan2(v[1],v[0])
    phi2 = arctan2(-1,0)

    alpha = abs(phi1-phi2)
    #make alpha an acute angle
    if alpha > pi:
        alpha = alpha - pi
    elif alpha == pi:
        raise NotImplementedError
    elif alpha > (pi/2) :
        alpha = pi - alpha

    return RR(width/sin(alpha))

def findPointFrom(orig,mid_angle,width):
    """
    Assume 2d.
    Return a point that is width(value) from the orig(point) in direction mid_angle
    """
    return (orig[0]+cos(mid_angle)*width,orig[1]+sin(mid_angle)*width)

def mod(x):
    """
    Given a number x, return its mod(1)
    """
    if x == 1:
        return 0
    else:
        return fractional(x)


class FIJK:

    def __init__(self,I_pair,J_pair,K_pair):

        self.I = I_pair[0]
        self.I_pi = I_pair[1] 

        self.J = J_pair[0]
        self.J_pi = J_pair[1] 
        self.K = K_pair[0]
        self.K_pi = K_pair[1] 

        self.buildPolyhedron()
        self.computeProjections()


    def buildPolyhedron(self):
        """
        Given I, J, K, calculate the Polyehedron FIJK
        """
        zero = [0] * self.I.ambient_dim()
        equationList = []
        inequationList = []

        for t in self.I.inequalities_list() :
            inequationList.append([t[0]] + t[1:] + zero)

        for t in self.I.equations_list() :
            equationList.append([t[0]] + t[1:] + zero)
        
        for t in self.J.inequalities_list() :
            inequationList.append([t[0]] + zero + t[1:])

        for t in self.J.equations_list() :
            equationList.append([t[0]] + zero + t[1:])    

        for t in self.K.inequalities_list() :
            inequationList.append([t[0]] + t[1:] + t[1:])
 
        for t in self.K.equations_list() :
            equationList.append([t[0]] + t[1:] + t[1:])  
        
        self.polyhedron = Polyhedron(ieqs = inequationList, eqns = equationList)

    def computeProjections(self):
        """
        Compute the projection of FIJK(self.polyhedron) to I, J, and K
        """
        self.p1 = Polyhedron(vertices = [v[:len(v)/2] for v in self.polyhedron.vertices_list()])
        self.p2 = Polyhedron(vertices = [v[len(v)/2:] for v in self.polyhedron.vertices_list()])

        #mod 1
        self.p3 = Polyhedron(vertices = [vector(v[:len(v)/2])+vector(v[len(v)/2:]) for v in self.polyhedron.vertices_list()])





class function_on_standard_2d_triangulation(function_on_complex):

    def __init__(self,matrix): 
        """
        A special case.
        The matrix has to be q x q.
        
        EXAMPLES::
            sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
            sage: m = matrix([[0,1/4,1/2],[1/4,1/2,3/4],[1/2,3/4,1]])
            sage: F = function_on_standard_2d_triangulation(m)
            sage: F.minimality_test()
            pi is minimal

            m = matrix([[0,2,2,2,2],[2,2,2,3,1],[2,2,4,2,2],[2,2,2,1,2],[2,2,2,2,3]])
            sage: F = function_on_standard_2d_triangulation(m)
            sage: F.minimality_test()
            pi is minimal            
        """
        #sage: m = matrix([[0,2,2,2,2],[2,2,2,3,1],[2,2,4,2,2],[2,2,2,1,2],[2,2,2,2,3]])
      
        if matrix.nrows() != matrix.ncols():
            raise ConstructError('Error: The matrix has to have the same number of rows and columns')  
        self.matrix = matrix
        self.q = matrix.nrows()
        self.num = self.q ^ 2 * 2  #of triangles

        ##Ring
        R.<x,y> = PolynomialRing(QQ,2)
        self.x = x
        self.y = y

        self.vertices_trangle_calc()
        self.polyhedron_fxn_list = [self.polyFxnPair(self.Tri3d(i)) for i in range (self.num)]
        function_on_complex.__init__(self,self.polyhedron_fxn_list)   

    def vertices_trangle_calc(self) :
        """
        Given the input matrix and the number q, calculate vertices, triangle pairs and the pair vertNumPair: (vertex, function value on this vertex)
        
        """
        q = self.q
        matrix = self.matrix
        ##find all the vertices
        list = [Rational(i/float(q)) for i in range (q + 1)]
        self.trig_vertices = zip([(0)] * (q + 1), list)
        for i in range (q) :
            v_temp = [list[i + 1] for j in range (q + 1)]
            self.trig_vertices.extend(zip(v_temp, list))
        #print self.trig_vertices

        ##compute all triangulation pair list t
        self.triangles = [(i, i+1, i+1+q) for i in range (q)]
        for i in range (q-1) :
            t_temp = [(j+(i+1)*(q+1), j+1+(i+1)*(q+1), j+q+1+(i+1)*(q+1))  for j in range (q)]
            self.triangles.extend(t_temp)
        for i in range (q) :
            t_temp = [(j+1+i*(q+1), j+1+q+i*(q+1), j+2+q+i*(q+1))  for j in range (q)]
            self.triangles.extend(t_temp)
        ##print "triangles:",self.triangles

        ##add numbers to the vertices
        list = []
        for i in range (q) :
            for j in range (q) :
                list += [matrix[j][i]]
            list += [matrix[0][i]]
        for j in range (q) :
            list += [matrix[j][0]]
        list += [0]
        ##print list
        self.vertNumPair = zip(self.trig_vertices, list)
        #print "vertNumPair:",self.vertNumPair
    

    def Tri3d(self,i) :
        """
        Given the number i corresponding to the ith triangle, return a list of the 3d vertices of this triangle.
        This fxn is used in __init__().        
        """
        list = [(self.vertNumPair[self.triangles[i][0]][0][0],self.vertNumPair[self.triangles[i][0]][0][1],self.vertNumPair[self.triangles[i][0]][1]),(self.vertNumPair[self.triangles[i][1]][0][0],self.vertNumPair[self.triangles[i][1]][0][1],self.vertNumPair[self.triangles[i][1]][1]),(self.vertNumPair[self.triangles[i][2]][0][0],self.vertNumPair[self.triangles[i][2]][0][1],self.vertNumPair[self.triangles[i][2]][1])]
        return list


    def Tri2d(self,i) :
        """
        Given the number i corresponding to the ith triangle, return a list of the 2d vertices of this triangle.
        """
        list = [(self.vertNumPair[self.triangles[i][0]][0][0],self.vertNumPair[self.triangles[i][0]][0][1]),(self.vertNumPair[self.triangles[i][1]][0][0],self.vertNumPair[self.triangles[i][1]][0][1]),(self.vertNumPair[self.triangles[i][2]][0][0],self.vertNumPair[self.triangles[i][2]][0][1])]
        return list    


    def polyFxnPair(self,l) :
        """
        Given a list of 3 points, return the pair (polyhedron, function)
        Polyhedron is the plane ax+by+cz=d 
        Function is a polynomial z = (d - ax - by)/c
        """
        a = l[0]
        b = l[1]
        c = l[2]
        para_a = (b[1] - a[1])*(c[2] - a[2]) - (c[1] - a[1])*(b[2] - a[2])
        para_b = (b[2] - a[2])*(c[0] - a[0]) - (c[2] - a[2])*(b[0] - a[0])
        para_c = (b[0] - a[0])*(c[1] - a[1]) - (c[0] - a[0])*(b[1] - a[1])
        para_d = para_a * a[0] + para_b * a[1] + para_c * a[2]

        return (Polyhedron(vertices=[[a[0],a[1]],[b[0],b[1]],[c[0],c[1]]],base_ring=QQ),(para_d - para_a * self.x - para_b * self.y) / para_c)


    def plot_complex(self) :
        """
        Plot the 2d complex with numbers on vertices and different colors for different gradients
        """
        ##list of pair(gradient vector, integer(which triangle))
        self.gradient_list = []
        for i in range (self.num) :
            f = self.polyhedron_fxn_list[i][1]
            temp = ((f.derivative(self.x), f.derivative(self.y)),i)
            self.gradient_list.append(temp)
        self.gradient_list.sort()

        ##calculate the number of different gradients
        self.grad_num = 1
        for i in range (self.num - 1) :
            if self.gradient_list[i][0] !=  self.gradient_list[i+1][0] :
                self.grad_num = self.grad_num + 1
        
        ##colors
        self.colors = rainbow(self.grad_num)

        image = Graphics()
        count = 0
        
        for i in range(self.num) :
            l = self.Tri2d(self.gradient_list[i][1])    
            if i == 0 : ##first color
                image += polygon2d(l, edgecolor="black", rgbcolor=self.colors[count])
            else :  ##others
                if self.gradient_list[i][0] !=  self.gradient_list[i-1][0] :
                    count += 1
                    image += polygon2d(l, edgecolor="black", rgbcolor=self.colors[count])
                else:  ##the same gradient
                     image += polygon2d(l, edgecolor="black", rgbcolor=self.colors[count])

        #print self.trig_vertices
        ##add boxes at the position of numbers
        for i in range ((self.q+1)^2) :
            #print self.trig_vertices[i]
            image += self.box(self.trig_vertices[i])

        ##add numbers to the vertices
        for i in range (self.q) :
            for j in range (self.q) :
                image += text(self.matrix[i][j], (j/float(self.q),i/float(self.q)), fontsize=16, color='black')
        for j in range (self.q) :
            image += text(self.matrix[0][j], (j/float(self.q),1), fontsize=16, color='black')
        for i in range (self.q) :
            image += text(self.matrix[i][0], (1,i/float(self.q)), fontsize=16, color='black')
        image += text('0', (1,1), fontsize=16, color='black')
        return image    

    def box(self,a) :
        """
        Given the coordinate in the center, return a white square
        """
        l = [(a[0]-1/5*(1/self.q), a[1]-1/5*(1/self.q)), (a[0]-1/5*(1/self.q), a[1]+1/5*(1/self.q)), (a[0]+1/5*(1/self.q), a[1]+1/5*(1/self.q)),(a[0]+1/5*(1/self.q), a[1]-1/5*(1/self.q))] 
        return polygon2d(l,color = 'white', fill = true, edgecolor = 'black', axes = False, transparent=False)



def addPeriodicPoints(point):
    """
    Given a point(list)
    Return a list of all points in [0,1]^k which is periodically equivalent to this point

    Example::
        sage: addPeriodicPoints([0,1])
        [[0, 0], [1, 0], [0, 1], [1, 1]]
        sage: addPeriodicPoints([0,1/2,1])
        [[0, 1/2, 0], [1, 1/2, 0], [0, 1/2, 1], [1, 1/2, 1]]
    """
    l = []
    n = len(point) 
    index_of_periodic_numbers = [i for i in range(n) if point[i] == 0 or point[i] == 1]
    k = len(index_of_periodic_numbers)
    #print index_of_periodic_numbers
    orderedList = [i+1 for i in range (k)]
    subsetList = Subsets(orderedList).list()

    for subset in subsetList:
        sub = [int(i+1 in subset) for i in range(k)]
        p = list(point)
        for i in range(len(sub)):
            p[index_of_periodic_numbers[i]] = sub[i]
        l.append(p)

    return l

def concatenate_ring_variables(R):
    """
    Example::
        sage: R = PolynomialRing(QQ,'x',4)
        sage: concatenate_ring_variables(R)
        'x0,x1,x2,x3'
    """
    s = R.variable_names()[0]

    for n in R.variable_names()[1:]:
        s = s + ','
        s = s + n
    return s

def gradientToRational(l):
    """
    Given a gradient vector from polynomial ring
    return a vector such that each element is rational
    """
    newl = [Rational(n) for n in l]
    return newl

    
def printPolyFxnPair(l):
    #l2 = tuple(l)
    #l2 = list(set(l2))
    for pair in l:
        print pair[0].vertices_list(), pair[1]
    return

def printPoly(l):
    #count = 0
    for poly in l:
        #print count
        print poly.vertices_list()
        #count = count+1
    return

def pointEqual(l1,l2):
    try:
        if len(l1) != len(l2):
            return false
        else:
            for i in range (len(l1)):
                if l1[i] != l2[i]:
                    return false
        return true
    except TypeError:
        return l1 == l2

def poly_isStrictlyInsidePoly(E,F):
    """
    Give polyhedrons E and F, return true if E is contained in F
    E is not equal to F
    check if every vertice of E is in F
    """
    #the dimension of F has to be higher than that of E
    if E.dim() > F.dim():
        #print E.polyhedron.vertices_list(), "is not in", F.polyhedron.vertices_list()
        return false

    #if every point of E is in F
    for v in E.vertices_list():
        if v not in F:
            #print v, "of", E.polyhedron.vertices_list(), "is not in", F.polyhedron.vertices_list()
            return false

    #if there exist a point in F not in E
    for v in F.vertices_list():
        if v not in E:
            #print E.polyhedron.vertices_list(), "is in", F.polyhedron.vertices_list()
            return true

    #print E.polyhedron.vertices_list(), "is not in", F.polyhedron.vertices_list()
    return false

def poly_isInsidePoly(E,F):
    """
    Give polyhedrons E and F, return true if E is contained in F
    E can be equal to F
    check if every vertice of E is in F
    """
    #the dimension of F has to be higher than that of E
    if E.dim() > F.dim():
        #print E.polyhedron.vertices_list(), "is not in", F.polyhedron.vertices_list()
        return false

    #if every point of E is in F
    for v in E.vertices_list():
        if v not in F:
            #print v, "of", E.polyhedron.vertices_list(), "is not in", F.polyhedron.vertices_list()
            return false

    return true

def combinationTwo(n):
    """
    Given a list of n different integers from 0 to n-1
    Return a list of 2-tuples, which are all possible combinations of these integers.
    Note: The result list includes [1,1], [2,1], and [1,2]
    """
    comb2 = []
    for i in range (n):
        a = []
        a.append(i)
        for j in range (n):
            a.append(j)
            comb2.append(a)
            a = a[0:1]
    #print 'comb2:',comb2
    return comb2


class ConstructError(Exception):
    def __init__(self,output):
        self.output = output
    def __str__(self):
        return repr(self.output)

class RangeError(Exception):
    def __init__(self,output):
        self.output = output
    def __str__(self):
        return repr(self.output)
