
class function_on_2dcomplex_matrix :

    def __init__(self, q, matrix):
        """
        The matrix has to be q x q.
        EXAMPLES::
        
            sage: m = matrix([[1,2],[3,4]])
            sage: F = function_on_2dcomplex_matrix(2,m)
        """
        self.q = q
        self.matrix = matrix

        ##Ring
        R.<x,y> = PolynomialRing(QQ,2)
        self.x = x
        self.y = y

        ##calculate verticies and triangles. 
        ##num is the number of triangles
        self.vertices_trangle_calc(q,matrix)
        self.num = q ^ 2 * 2

        ##list of pairs (plane, function)
        self.face_fxn_list = [self.faceFxnPair(self.Tri3d(i)) for i in range (self.num)]
        ##print self.face_fxn_list

        ##list of pair(gradient vector, integer(which triangle))
        self.gradient_list = []
        for i in range (self.num) :
            f = self.face_fxn_list[i][1]
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



    def __call__(self,x0,y0) :
        """
        Evaluates self at (x0,y0).

        EXAMPLES:
            sage: m = matrix([[1,2],[3,4]])
            sage: F = function_on_2dcomplex_matrix(2,m)
            sage: F(1/2,1/2)
            4
        
        """
        n = self.find(x0,y0)
        f = self.face_fxn_list[n][1]
        a = f(x = x0, y = y0)
        return a

    def find(self,x0,y0) :
        """
        Given a pair of coordinates(x0,y0), return i, where the point is in the ith triangle. 
        This fxn is used in __call__().
        """
        x = x0
        y = y0
        a = 0
        b = 0 
        while x > (1/self.q) :
            a += 1
            x -= 1/self.q

        while y > (1/self.q) :
            b += 1
            y -= 1/self.q

        if x + y < (1/self.q) :
            num = self.q * a + b
        else :
            num = self.q * a + b + self.q * self.q
        return num



    def Tri3d(self,i) :
        """
        Given the number i corresponding to the ith triangle, return a list of the 3d vertices of this triangle.
        This fxn is used in __init__().        
        """
        list  = [(self.vertNumPair[self.triangles[i][0]][0][0],self.vertNumPair[self.triangles[i][0]][0][1],self.vertNumPair[self.triangles[i][0]][1]),(self.vertNumPair[self.triangles[i][1]][0][0],self.vertNumPair[self.triangles[i][1]][0][1],self.vertNumPair[self.triangles[i][1]][1]),(self.vertNumPair[self.triangles[i][2]][0][0],self.vertNumPair[self.triangles[i][2]][0][1],self.vertNumPair[self.triangles[i][2]][1])]
        return list

    def Tri2d(self,i) :
        """
        Given the number i corresponding to the ith triangle, return a list of the 2d vertices of this triangle.
        """
        list  = [(self.vertNumPair[self.triangles[i][0]][0][0],self.vertNumPair[self.triangles[i][0]][0][1]),(self.vertNumPair[self.triangles[i][1]][0][0],self.vertNumPair[self.triangles[i][1]][0][1]),(self.vertNumPair[self.triangles[i][2]][0][0],self.vertNumPair[self.triangles[i][2]][0][1])]
        return list    


    def faceFxnPair(self,l) :
        """
        Given a list of 3 points, return the pair (face, function)
        Face is the plane ax+by+cz=d as a class of polyhedra
        Function is a polynomial z = (d - ax - by)/c
        """
         a = l[0]
        b = l[1]
        c = l[2]
        para_a = (b[1] - a[1])*(c[2] - a[2]) - (c[1] - a[1])*(b[2] - a[2])
        para_b = (b[2] - a[2])*(c[0] - a[0]) - (c[2] - a[2])*(b[0] - a[0])
        para_c = (b[0] - a[0])*(c[1] - a[1]) - (c[0] - a[0])*(b[1] - a[1])
        para_d = para_a * a[0] + para_b * a[1] + para_c * a[2]
        ##print "%f %f %f %f" %(para_a,para_b,para_c,para_d)
        ##print (para_d - para_a * self.x - para_b * self.y) / para_c
        return (Polyhedron(eqns=[(-para_d, para_a, para_b, para_c)],base_ring=QQ),(para_d - para_a * self.x - para_b * self.y) / para_c)
        


    def vertices_trangle_calc(self,q,matrix) :
        """
        Given the input matrix and the number q, calculate vertices, triangle pairs and the pair vertNumPair: (vertex, function value on this vertex)
        
        """
        ##find all the vertices
        list = [i/float(q) for i in range(q + 1)]
        self.vertices = zip([0] * (q + 1), list)
        for i in range (q) :
            v_temp = [list[i + 1] for j in range (q + 1)]
            self.vertices.extend(zip(v_temp, list))
        ##print self.vertices

        ##compute all triangulation pair list t
        self.triangles = [(i, i+1, i+1+q) for i in range (q)]
        for i in range (q-1) :
            t_temp = [(j+(i+1)*(q+1), j+1+(i+1)*(q+1), j+q+1+(i+1)*(q+1))  for j in range (q)]
            self.triangles.extend(t_temp)
        for i in range (q) :
            t_temp = [(j+1+i*(q+1), j+1+q+i*(q+1), j+2+q+i*(q+1))  for j in range (q)]
            self.triangles.extend(t_temp)
        ##print "triangles:"
        ##print self.triangles

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
        self.vertNumPair = zip(self.vertices, list)
        ##print "vertex, function(vertex):"
        ##print self.vertNumPair

    def box(self,a) :
        """
        Given the coordinate in the center, return a white square
        """
        l = [(a[0]-1/5*(1/self.q), a[1]-1/5*(1/self.q)), (a[0]-1/5*(1/self.q), a[1]+1/5*(1/self.q)), (a[0]+1/5*(1/self.q), a[1]+1/5*(1/self.q)),(a[0]+1/5*(1/self.q), a[1]-1/5*(1/self.q))] 
        return polygon2d(l,color = 'white', fill = true, edgecolor = 'black', axes = False, transparent=False)


    def plot_complex(self) :
        """
        Plot the 2d complex with numbers on vertices and different colors for different gradients
        """
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

        ##add boxes at the position of numbers
        for i in range ((self.q+1)^2) :
            image += self.box(self.vertices[i])

        ##add numbers to the vertices
        for i in range (self.q) :
            for j in range (self.q) :
                image += text(self.matrix[i][j], (j/self.q,i/self.q), fontsize=16, color='black')
        for j in range (self.q) :
            image += text(self.matrix[0][j], (j/self.q,1), fontsize=16, color='black')
        for i in range (self.q) :
            image += text(self.matrix[i][0], (1,i/self.q), fontsize=16, color='black')
        image += text('0', (1,1), fontsize=16, color='black')
        return image



    def plot_function(self)    :
        """
        Plot the function on the complex (3d).
        Haven't tried the code since there was a bug in the 3d plot software.
        """
        image = Graphics()
    
        count = 0
        ##triangles
        for i in range (self.num) :
            l = self.Tri3d(self.gradient_list[i][1])
                        
            if i == 0 : ##first color
                image += polygon3d(l, color=self.colors[count])
            else :  ##others
                if self.gradient_list[i][0] !=  self.gradient_list[i-1][0] :
                    count += 1
                    image += polygon3d(l, color=self.colors[count])
                else:  ##the same gradient
                     image += polygon3d(l, color=self.colors[count])

        vertNumPair = self.vertNumPair
        ##triangle edges
        ##horizontal    
        ##for i in range (self.q+1) :
        ##    for j in range (self.q) :
        ##        image += line([(vertNumPair[i+3*j][0][0],vertNumPair[i+3*j][0][1],vertNumPair[i+3*j][1]), (vertNumPair[3*j+i+1][0][0],vertNumPair[i+1+3*j][0][1],vertNumPair[i+1+3*j][1])], color = 'black')

        ##vertical
        ##for i in range (self.q+1) :
        ##    for j in range (self.q) :
        ##        image += line([(vertNumPair[j+3*i][0][0],vertNumPair[j+3*i][0][1],vertNumPair[j+3*i][1]), (vertNumPair[3*i+j+1][0][0],vertNumPair[j+1+3*i][0][1],vertNumPair[j+1+3*i][1])], color = 'black')

        ##diagonal
        ##for i in range (self.q-1) : ##lower triangle
        ##    for j in range (i+1) :
        ##        image += line([(vertNumPair[i+1+self.q*j][0][0],vertNumPair[i+1+self.q*j][0][1],vertNumPair[i+1+self.q*j][1]), (vertNumPair[i+1+self.q*(j+1)][0][0],vertNumPair[i+1+self.q*(j+1)][0][1],vertNumPair[i+1+self.q*(j+1)][1])], color = 'black')
        ##for i in range (self.q) :  ##middile diagonal
        ##    image += line([(vertNumPair[self.q*(i+1)][0][0],vertNumPair[self.q*(i+1)][0][1],vertNumPair[self.q*(i+1)][1]), (vertNumPair[self.q*(i+2)][0][0],vertNumPair[self.q*(i+2)][0][1],vertNumPair[self.q*(i+2)][1])], color = 'black')
        ##c = self.q^2 - 1
        ##for i in range (self.q-1) : ##upper triangle
        ##    for j in range (i+1) :
        ##        image += line([(vertNumPair[c-(self.q+1)*i+self.q*j][0][0],vertNumPair[c-(self.q+1)*i+self.q*j][0][1],vertNumPair[c-(self.q+1)*i+self.q*j][1]), (vertNumPair[c-(self.q+1)*i+self.q*(j+1)][0][0],vertNumPair[c-(self.q+1)*i+self.q*(j+1)][0][1],vertNumPair[c-(self.q+1)*i+self.q*(j+1)][1])], color = 'black')    
            



        ##bottom lines
        for i in range (self.q-1):    ##vertical    
            image += line([((i+1)/self.q,0,0),((i+1)/self.q,1,0)], color = 'black')
        for i in range (self.q-1):    ##horizontal    
            image += line([(0,(i+1)/self.q,0),(1,(i+1)/self.q,0)], color = 'black')
        ##diagonal
        for i in range (self.q-1):        
            image += line([((i+1)/self.q,0,0),(0,(i+1)/self.q,0)], color = 'black')
        image += line([(0,1,0),(1,0,0)], color = 'black')
        for i in range (self.q-1):        
            image += line([((i+1)/self.q,1,0),(1,(i+1)/self.q,0)], color = 'black')
        return image.show(viewer = "tachyon")

    













