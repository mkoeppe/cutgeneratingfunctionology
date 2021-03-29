from six.moves import range
def vertex_enumeration(polytope, exp_dim=-1, vetime=False):
    r"""
    Returns the vertices of the polytope.

    - Do preprocessing if exp_dim >= igp.exp_dim_prep, i.e., call the function redund provided by lrslib to remove redundant inequalities.
    - If normaliz is installed, use it to enumerate vertices after preprocessing, else use lrs vertex enumeration if exp_dim >= exp_dim_lrs.
    - Use ppl vertex enumeration for small dimensions.
    - Print the vertex enumeration running time if vetime is ``True``.
    
    EXAMPLE::

        sage: from cutgeneratingfunctionology.igp import *
        sage: x = Variable(0)
        sage: y = Variable(1)
        sage: cs = Constraint_System()
        sage: cs.insert(x >= 0)
        sage: cs.insert(y >= 0)
        sage: cs.insert(x <= 1)
        sage: cs.insert(y <= 1)
        sage: polytope = C_Polyhedron(cs)
        sage: vertex_enumeration(polytope)
        Generator_System {point(1/1, 0/1), point(1/1, 1/1), point(0/1, 1/1), point(0/1, 0/1)}
    """
    if vetime:
        st = os.times();
    try:
        import PyNormaliz
        if exp_dim < exp_dim_prep:
            extreme_points = polytope.minimized_generators()
        else:
            # do preprocessing by lrslib redund
            cs = polytope.constraints()
            cs_prep = remove_redundancy_from_cs(cs)
            ieqs = []; eqns=[]
            for c in cs_prep:
                coef = [c.inhomogeneous_term()] + list(c.coefficients())
                if c.is_equality():
                    eqns.append(coef)
                else:
                    ieqs.append(coef)       
            sage_polyhedron = Polyhedron(ieqs=ieqs, eqns=eqns, backend='normaliz')
            extreme_points = sage_polyhedron.vertices()
    except ImportError:
        if exp_dim >= exp_dim_prep:
            # do preprocessing
            cs = polytope.constraints()
            if exp_dim >= exp_dim_lrs:
                # preprocessing and vertex enumertation using redund + lrs
                cs_prep_lrs_str = remove_redundancy_from_cs(cs, return_lrs=True)
                extreme_points = lrs_lrsinput_pploutput(cs_prep_lrs_str)
            else:
                # preprocessing and vertex enumertation using redund + ppl
                cs_prep = remove_redundancy_from_cs(cs)
                polytope = C_Polyhedron(cs_prep)
                extreme_points = polytope.minimized_generators()
        else:
            # no preprocessing
            if exp_dim >= exp_dim_lrs:
                # vertex enumertation using lrs
                cs = polytope.constraints()
                cs_lrs_str = convert_pplcs_to_lrs(cs)
                extreme_points = lrs_lrsinput_pploutput(cs_lrs_str) 
            else:
                # vertex enumertation using ppl
                extreme_points = polytope.minimized_generators()
    if vetime:
        et = os.times(); 
        logging.info("user=%s, sys=%s, child user=%s, child sys=%s" %(et[0]-st[0], et[1]-st[1], et[2]-st[2], et[3]-st[3]))
        t = sum([et[i]-st[i] for i in range(4)]);
        logging.info("Vertex enumeration time = %s" % t)
    return extreme_points

#output_dir = "./"
## override this in a file named "config.sage" if necessary

def write_panda_format_cs(cs, fname=None, newcode=True):
    if fname:
        dir = output_dir+"profiler/panda/"
        mkdir_p(dir)
        filename = open(dir+fname, "w")
    else:
        filename = sys.stdout
    if newcode:
        panda_string = convert_pplcs_to_panda(cs)
        print(panda_string, file=filename)
    else:
        print('Inequalities:', file=filename)
        for c in cs:
            for x in c.coefficients():
                print(-x, end=' ', file=filename)
            print(-c.inhomogeneous_term(), file=filename)
            if c.is_equality():
                for x in c.coefficients():
                    print(x, end=' ', file=filename)
                print(c.inhomogeneous_term(), file=filename)
    if fname:
        filename.close()
    return

def convert_pplcs_to_panda(cs):
    s = "Names\n"
    q = cs.space_dimension() - 1;
    for i in range(1, q+2):
        s += "x%s " % i
    s += "\n"
    s += 'Equations\n'
    for c in cs:
        if c.is_equality():
            coefs = c.coefficients()
            for i in range(q+1):
                x = coefs[i]
                if x > 0:
                    s += '+%sx%s ' % (x, i+1)
                elif x < 0:
                    s += '%sx%s ' % (x, i+1)
            s += '= '
            s += '%s\n' % -c.inhomogeneous_term()
    s += 'Inequalities\n'
    for c in cs:
        if not c.is_equality():
            coefs = c.coefficients()
            for i in range(q+1):
                x = coefs[i]
                if x > 0:
                    s += '+%sx%s ' % (x, i+1)
                elif x < 0:
                    s += '%sx%s ' % (x, i+1)
            s += '>= '
            s += '%s\n' % -c.inhomogeneous_term()
    return s

def write_porta_ieq(q, f, destdir=None):
    vertices_color = initial_vertices_color(q, f);
    cs = initial_cs(q, f, vertices_color)
    if destdir is None:
        fname = None
    else:
        dir = output_dir+"profiler/porta/"
        mkdir_p(dir)
        fname = dir + "porta_q%sf%s.ieq" % (q, f)
    write_porta_format_cs(cs, q=q, f=f, fname=fname)
    return

def write_porta_format_cs(cs, q=None, f=None, fname=None):
    # q, f are used to find a valid point -- gmic function
    # this valid point is required by 'traf' as the 0 is not in the feasible region.
    if fname:
        filename = open(fname, "w")
    else:
        filename = sys.stdout
    porta_string = convert_pplcs_to_porta(cs, q=q, f=f)
    print(porta_string, file=filename)
    if fname:
        filename.close()
    return

def convert_pplcs_to_porta(cs, q=None, f=None):
    # q, f are used to find a valid point -- gmic function
    # this valid point is required by 'traf' as the 0 is not in the feasible region.
    s = ""
    s += "DIM = %s\n" % cs.space_dimension()
    s += "\n"
    if not q is None:
        s += 'VALID\n'
        s += '0 '
        for i in range(1, f):
            s += '%s/%s ' % (i, f)
        s += '1 '
        for i in range(q - f - 1, 0, -1):
            s += '%s/%s ' % (i, q - f)
        s += '0\n'
        s += '\n'
    s += 'INEQUALITIES_SECTION\n'
    for c in cs:
        coefs = c.coefficients()
        for i in range(q+1):
            x = coefs[i]
            if x > 0:
                s += '+%sx%s ' % (x, i+1)
            elif x < 0:
                s += '%sx%s ' % (x, i+1)
        if c.is_equality():
            s += '== '
        else:
            s += '>= '
        s += '%s\n' % -c.inhomogeneous_term()
    s += '\n'
    s += 'END\n'
    return s

def write_lrs_ine(q, f, destdir=None):
    vertices_color = initial_vertices_color(q, f);
    cs = initial_cs(q, f, vertices_color)
    if destdir is None:
        fname = None
    else:
        dir = output_dir+"profiler/lrs/"
        mkdir_p(dir)
        fname = dir + "lrs_q%sf%s.ine" % (q, f)
    write_lrs_format_cs(cs, fname=fname)
    return

def write_lrs_format_cs(cs, fname=None):
    if fname:
        filename = open(fname, "w")
    else:
        filename = sys.stdout
    lrs_string = convert_pplcs_to_lrs(cs, fname=fname)
    print(lrs_string, file=filename)
    if fname:
        filename.close()
    return

def convert_pplcs_to_lrs(cs, fname=None):
    if fname:
        s = fname + '\n'
    else:
        s = ""
    s += "H-representation" + '\n'
    m = len(cs)
    n = cs.space_dimension() + 1
    k = 0
    linearities = []
    for i in range(m):
        c = cs[i]
        if c.is_equality():
            k += 1
            linearities.append(i)
    if k > 0:
        s += "linearity %s " % k
    for i in linearities:
        s += str(i + 1)  + ' '
    s += '\n'
    s += "begin\n"
    s += "%s %s rational\n" %(m, n)
    for c in cs:
        s += str(c.inhomogeneous_term()) + ' '
        for x in c.coefficients():
            s += str(x) + ' '
        s += '\n'
    s += 'end\n'
    return s

def read_lrs_to_cs_or_gs(fname):
    myfile = open(fname, "r")
    lrs_string = myfile.read()
    myfile.close()
    return convert_lrs_to_ppl(lrs_string)

def convert_lrs_to_ppl(lrs_string):
    r"""
    Convert lrs format H-representation to ppl cs;
    or lrs format V-representation (of a polytope) to ppl gs.

    COPY from src/geometry/polyhedron/backend_cdd.py and edit the function
    ``Polyhedron_cdd._init_from_cdd_output(self, cdd_output_string)``


    EXAMPLE::

        sage: from cutgeneratingfunctionology.igp import *
        sage: x = Variable(0)
        sage: y = Variable(1)
        sage: cs = Constraint_System()
        sage: cs.insert(x >= 0)
        sage: cs.insert(y >= 0)
        sage: cs.insert(x <= 1)
        sage: cs.insert(y <= 1)
        sage: in_str = convert_pplcs_to_lrs(cs)
        sage: out_str = lrs_redund(in_str)     #      optional - lrslib
        sage: type(out_str)                    # py3, optional - lrslib
        <class 'str'>
        sage: convert_lrs_to_ppl(out_str)      #      optional - lrslib
        Constraint_System {x0>=0, x1>=0, -x0+1>=0, -x1+1>=0}
    """
    cddout=lrs_string.splitlines()

    # nested function
    def expect_in_cddout(expected_string):
        l = cddout.pop(0).strip()
        if l != expected_string:
            raise ValueError('Error while parsing cdd output: expected "'
                               +expected_string+'" but got "'+l+'".\n')
    # nested function
    def cdd_linearities():
        l = cddout[0].split()
        if l[0] != "linearity":
            return []
        cddout.pop(0)
        assert len(l) == int(l[1])+2, "Not enough linearities given"
        return [int(i)-1 for i in l[2:]]  # make indices pythonic

    # nested function
    def cdd_convert(string, field=QQ):
        r"""
        Converts the cdd output string to a QQ numerical value.
        """
        return [field(x) for x in string.split()]

    # nested function
    def find_in_cddout(expected_string):
        r"""
        Find the expected string in a list of strings, and
        truncates ``cddout`` to start at that point. Returns
        ``False`` if search fails.
        """
        for pos in range(0,len(cddout)):
            l = cddout[pos].strip();
            if l==expected_string:
                # must not assign to cddout in nested function
                for i in range(0,pos+1):
                    cddout.pop(0)
                return True
        return False

    def lrs_row_to_linear_expression(l):
        rational_list = cdd_convert(l)
        num_list = [x.numerator() for x in rational_list]
        den_list = [x.denominator() for x in rational_list]
        common_den = lcm(den_list)
        ihom = int(common_den / den_list[0]) * num_list[0]
        coef = [int(common_den / den_list[i]) * num_list[i] for i in range(1, len(rational_list))]
        return Linear_Expression(coef, ihom)

    def lrs_row_to_point(l):
        rational_list = cdd_convert(l)
        num_list = [x.numerator() for x in rational_list]
        den_list = [x.denominator() for x in rational_list]
        common_den = lcm(den_list)
        coef = [int(common_den / den_list[i]) * num_list[i] for i in range(len(rational_list))]
        return ppl_point(Linear_Expression(coef, 0), common_den)

    if find_in_cddout('V-representation'):
        # Suppose it's the V-representation of a polytope.
        # Return the ppl generator_system that contains the vertices.
        # Sometimes the number of vertices is unknown from the input file, so read till 'end'
        # ex: in the output of lrs vertex enumertaion file.ext.
        #raise NotImplementedError, "V-representation Not implemented."
        gs = Generator_System()
        equations = cdd_linearities()
        expect_in_cddout('begin')
        l = cddout.pop(0).split()
        n = int(l[1]) - 1
        l = cddout.pop(0).strip()
        while l != 'end':
            l_type = l[0]
            l = l[1:]
            #if (i in equations) or l_type == '0':
            #    raise NotImplementedError, "V-representation of line or ray is NOT implemented."
            vertex = lrs_row_to_point(l)
            gs.insert(vertex)
            l = cddout.pop(0).strip()
        return gs

    if find_in_cddout('H-representation'):
        cs = Constraint_System()
        equations = cdd_linearities()
        expect_in_cddout('begin')
        l = cddout.pop(0).split()
        m = int(l[0])
        n = int(l[1]) - 1
        fn = [ Variable(i) for i in range(n) ]
        for i in range(m):
            l = cddout.pop(0).strip()
            lin_exp = lrs_row_to_linear_expression(l)
            if i in equations:
                cs.insert(lin_exp == 0)
            else:
                cs.insert(lin_exp >= 0)
        expect_in_cddout('end')
        return cs

    raise ValueError('Input is neither a V-representation nor H-representation')

from sage.misc.temporary_file import tmp_filename
from subprocess import Popen, PIPE

def lrs_redund(in_str, verbose=False):
    r"""
    To remove redundant inequalities from an H-representation or
    input points that are not vertices from a V-representation
    use the command 'redund' from lrslib.

    Input: lrs format in_str; Output: lrs format out_str;

    Copy and edit from ``def _volume_lrs(self, verbose=False)``,
    http://www.sagenb.org/src/geometry/polyhedron/base.py

    EXAMPLE::

        sage: from cutgeneratingfunctionology.igp import *
        sage: x = Variable(0)
        sage: y = Variable(1)
        sage: cs = Constraint_System()
        sage: cs.insert(x >= 0)
        sage: cs.insert(y >= 0)
        sage: cs.insert(x <= 1)
        sage: cs.insert(y <= 1)
        sage: in_str = convert_pplcs_to_lrs(cs)
        sage: out_str = lrs_redund(in_str)      #      optional - lrslib
        sage: type(out_str)                     # py3, optional - lrslib
        <class 'str'>
    """
    #if is_package_installed('lrslib') != True:
    #    print 'You must install the optional lrs package ' \
    #          'for this function to work'
    #    raise NotImplementedError

    in_filename = tmp_filename()
    in_file = open(in_filename,'w')
    in_file.write(in_str)
    in_file.close()
    if verbose: print(in_str)
    redund_procs = Popen(['redund',in_filename],stdin = PIPE, stdout=PIPE, stderr=PIPE)
    out_bytes, err = redund_procs.communicate()
    out_str = out_bytes.decode("utf-8") # lrslib changed, output is of type bytes
    if verbose:
        print(out_str)
    return out_str

def remove_redundancy_from_cs(cs, verbose=False, return_lrs=False):
    r"""
    Remove redundant inequalities from cs using the command 'redund' from lrslib.
    Return the new cs without redundancy
    in ppl format if return_lrs==False (default), or in lrs format if return_lrs==True 
    """
    in_str = convert_pplcs_to_lrs(cs, fname=None)
    out_str = lrs_redund(in_str, verbose=verbose)
    if return_lrs:
        return out_str
    else:
        return convert_lrs_to_ppl(out_str)

def lrs_lrs(in_str, verbose=False):
    r"""
    Use the command 'lrs' from lrslib.
    Input: lrs format in_str; Output: lrs format out_str;
    """
    #if is_package_installed('lrslib') != True:
    #    print 'You must install the optional lrs package ' \
    #          'for this function to work'
    #    raise NotImplementedError
    in_filename = tmp_filename()
    in_file = open(in_filename,'w')
    in_file.write(in_str)
    in_file.close()
    if verbose: print(in_str)

    redund_procs = Popen(['lrs',in_filename],stdin = PIPE, stdout=PIPE, stderr=PIPE)
    out_bytes, err = redund_procs.communicate()
    out_str = out_bytes.decode("utf-8") # lrslib changed, output is of type bytes
    if verbose:
        print(out_str)
    return out_str
    if verbose:
        print(out_str)
    return out_str

def lrs_lrsinput_pploutput(in_str):
    r"""
    Use the command lrs from lrslib.

    Input: 
        lrs format in_str; 

    Output: 
        ppl format extreme_points;
    
    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: cube_in_str = "cube\n*cube of side 2 centred at origin\nH-representation\nbegin\n6  4 rational" + "\n1 1 0 0\n1 0 1 0\n1 0 0 1\n1 -1 0 0\n1 0 -1 0\n1 0 0 -1\nend"
        sage: lrs_lrsinput_pploutput(cube_in_str)  # optional - lrslib
        Generator_System {point(1/1, 1/1, 1/1), point(-1/1, 1/1, 1/1), point(1/1, -1/1, 1/1), point(-1/1, -1/1, 1/1), point(1/1, 1/1, -1/1), point(-1/1, 1/1, -1/1), point(1/1, -1/1, -1/1), point(-1/1, -1/1, -1/1)}
        sage: lrs_q5f3_str = "lrs_q5f3\nH-representation\nlinearity 5 1 2 3 13 21\nbegin\n21 7 rational" + "\n0 1 0 0 0 0 0\n0 0 0 0 0 0 1\n-1 0 0 0 1 0 0\n0 0 1 0 0 0 0\n1 0 -1 0 0 0 0\n0 0 0 1 0 0 0\n1 0 0 -1 0 0 0\n0 0 0 0 1 0 0" + "\n1 0 0 0 -1 0 0\n0 0 0 0 0 1 0\n1 0 0 0 0 -1 0\n0 0 2 -1 0 0 0\n0 0 1 1 -1 0 0\n0 0 1 0 1 -1 0\n0 -1 1 0 0 1 0" + "\n0 0 0 2 0 -1 0\n0 -1 0 1 1 0 0\n0 0 -1 1 0 1 0\n0 0 -1 0 2 0 0\n0 0 0 -1 1 1 0\n0 0 0 0 1 -2 0\nend"
        sage: lrs_lrsinput_pploutput(lrs_q5f3_str)  # optional - lrslib
        Generator_System {point(0/6, 2/6, 4/6, 6/6, 3/6, 0/6), point(0/4, 3/4, 1/4, 4/4, 2/4, 0/4)}    
    r"""
    v_lrs_str = lrs_lrs(in_str)
    extreme_points = convert_lrs_to_ppl(v_lrs_str)
    return extreme_points
   
def lcdd_rational(in_str, verbose=False):
    r"""
    Use the command ``lcdd_gmp`` from cddlib.

    Input: cdd format in_str; Output: cdd format out_str;
    """
    in_filename = tmp_filename()
    in_file = open(in_filename,'w')
    in_file.write(in_str)
    in_file.close()
    if verbose: 
        print(in_str)
    redund_procs = Popen(['lcdd_gmp',in_filename],stdin = PIPE, stdout=PIPE, stderr=PIPE)
    out_str, err = redund_procs.communicate()
    if verbose:
        print(out_str)
    return out_str


def write_normaliz_in(q, f, destdir=None):
    vertices_color = initial_vertices_color(q, f);
    cs = initial_cs(q, f, vertices_color)
    if destdir is None:
        fname = None
    else:
        dir = output_dir+"profiler/normaliz/"
        mkdir_p(dir)
        fname = dir + "normaliz_q%sf%s.in" % (q, f)
    write_normaliz_format_cs(cs, fname=fname)
    return

def write_normaliz_remove_redund_in(q, f):
    infname = output_dir+"profiler/lrs/lrs_remove_redund_q%sf%s.ine" % (q, f)
    cs = read_lrs_to_cs_or_gs(infname)
    fname = output_dir+"profiler/normaliz/normaliz_remove_redund_q%sf%s.in" % (q, f)
    write_normaliz_format_cs(cs, fname=fname)
    return

def write_normaliz_format_cs(cs, fname=None):
    if fname:
        filename = open(fname, "w")
    else:
        filename = sys.stdout
    normaliz_string = convert_pplcs_to_normaliz(cs)
    print(normaliz_string, file=filename)
    if fname:
        filename.close()
    return

def convert_pplcs_to_normaliz(cs):
    m = len(cs)
    n = cs.space_dimension() + 1
    s = 'amb_space %s\n' % n
    s += 'hom_constraints %s\n' % m
    for c in cs:
        for x in c.coefficients():
            s += str(x) + ' '
        if c.is_equality():
            s += '= '
        else:
            s += '>= '
        s += '%s\n' % -c.inhomogeneous_term()
    s += 'grading\n'
    s += 'unit_vector %s\n' % n
    s += 'ExtremeRays\n'
    return s
