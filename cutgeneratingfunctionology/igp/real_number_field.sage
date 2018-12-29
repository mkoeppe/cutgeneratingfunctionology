
import sage.rings.number_field.number_field
from sage.rings.number_field.number_field import NumberField_absolute, NumberField_quadratic

import sage.rings.number_field.number_field_element
from sage.rings.number_field.number_field_element import NumberFieldElement_absolute

from sage.structure.richcmp import op_EQ, op_NE, op_LE, op_GE, op_LT, richcmp_not_equal

class RealNumberFieldElement(NumberFieldElement_absolute):

    def embedded(self):
        e = getattr(self, '_embedded', None)
        if e is None:
            parent = self.parent()
            embedding = parent.coerce_embedding()
            self._embedded = e = embedding(self)
        return e

    def _richcmp_(left, right, op):    # In Sage >= 7.1, need to specialize this function.
        #print "_richcmp_", left, right, op
        if NumberFieldElement._richcmp_(left, right, op_EQ):
            return op == op_EQ or op == op_LE or op == op_GE
        elif op == op_NE:
            return True
        elif op == op_EQ:
            return False

        if left.embedded() == right.embedded():
            raise NotImplementedError("Precision of real interval field not sufficient to continue")
        return richcmp_not_equal(left.embedded(), right.embedded(), op)
    
    def __abs__(self):
        if self.sign() >= 0:
            return self
        else:
            return -self

    def sign(self):
        parent = self.parent()
        return cmp(self, parent._zero_element)

    def __float__(self):
        embedded = self.embedded()
        result = (float(embedded.lower()) + float(embedded.upper())) / 2
        return result

    def __repr__(self):
        global show_RNFElement_by_embedding
        if show_RNFElement_by_embedding:
            embedded = self.embedded()
            symbolic = getattr(self, '_symbolic', None)
            if symbolic is None:
                return '%s' % embedded
            else:
                return '<%s=%s>' % (symbolic, embedded)
        else:
            return NumberFieldElement_absolute._repr_(self)

    def _latex_(self):
        symbolic = getattr(self, '_symbolic', None)
        if symbolic is None:
            parent = self.parent()
            embedding = parent._exact_embedding
            symbolic = embedding(self)
        return "%s" % (latex(symbolic))

    def _maxima_(self, session=None):
        symbolic = getattr(self, '_symbolic', None)
        if symbolic is None:
            raise NotImplementedError("Cannot make %s a Maxima number" % self)
        else:
            return symbolic._maxima_()

    def _add_(self, other):
        result = NumberFieldElement_absolute._add_(self, other)
        self_e = self.embedded()
        other_e = other.embedded()
        result._embedded = self_e + other_e
        return result

    def _sub_(self, other):
        result = NumberFieldElement_absolute._sub_(self, other)
        self_e = self.embedded()
        other_e = other.embedded()
        result._embedded = self_e - other_e
        return result

    def __neg__(self):
        result = NumberFieldElement_absolute._neg_(self)
        self_e = self.embedded()
        result._embedded = -self_e
        return result

    def _mul_(self, other):
        result = NumberFieldElement_absolute._mul_(self, other)
        self_e = self.embedded()
        other_e = other.embedded()
        result._embedded = self_e * other_e
        return result
        
    def _div_(self, other):
        result = NumberFieldElement_absolute._div_(self, other)
        self_e = self.embedded()
        other_e = other.embedded()
        result._embedded = self_e / other_e
        return result

    def __hash__(self):
        if not hasattr(self, '_hash'):
            self._hash = NumberFieldElement_absolute.__hash__(self)
        return self._hash

class RealNumberField_absolute(NumberField_absolute):
    r"""
    A ``RealNumberField`` knows its embedding into a ``RealIntervalField``
    and use that for ``<``, ``>`` comparisons.
    ``==`` comparison is exact, using the underlying numberfield.
    
    A ``RealNumberField`` also knows an embedding into an exact field (``SR``)
    for the purpose of latexing.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: field, field_values, morphism = number_field_elements_from_algebraics((sqrt(2), sqrt(3)))
        sage: emb_field = RealNumberField(field.polynomial(), 'a', 
        ....:                 embedding=morphism(field.gen(0)),
        ....:                 exact_embedding=SR(morphism(field.gen(0))))
        sage: hom = field.hom([emb_field.gen(0)])
        sage: Integer(7)/5 < hom(field_values[0])
        True
        sage: hom(field_values[0]) < Integer(3)/2
        True
        sage: hom(field_values[0]) < hom(field_values[1])
        True
        sage: Integer(3)/2 < hom(field_values[1])
        True
        sage: hom(field_values[1]) < 2
        True
        sage: logging.disable(logging.INFO)
        sage: field = nice_field_values([sqrt(2),sqrt(3)])[0].parent()
        sage: field
        Real Number Field in `a` as the root of the defining polynomial y^4 - 4*y^2 + 1 near 0.5176380902050415?
        sage: field1 = NumberField(field.polynomial(),"a")
        sage: field1
        Number Field in a with defining polynomial y^4 - 4*y^2 + 1
        sage: field1 == field or field == field1
        False
        sage: field2 = VectorSpace(field, 2).base_field()
        sage: field2
        Real Number Field in `a` as the root of the defining polynomial y^4 - 4*y^2 + 1 near 0.5176380902050415?
        sage: field == field2 and field2 == field
        True
    """

    def __init__(self, polynomial, name=None, latex_name=None, check=True, embedding=None,
                 assume_disc_small=False, maximize_at_primes=None, exact_embedding=None):
        r"""
        Create a real number field.
        """
        NumberField_absolute.__init__(self, polynomial, name=name, check=check,
                                      embedding=embedding, latex_name=latex_name,
                                      assume_disc_small=assume_disc_small, maximize_at_primes=maximize_at_primes)
        self._standard_embedding = True
        self._element_class = RealNumberFieldElement
        self._zero_element = self(0)
        self._one_element =  self(1)
        self._exact_embedding = self.hom([exact_embedding])

    def _repr_(self):
        return "Real Number Field in `%s` as the root of the defining polynomial %s near %s"%(
                   self.variable_name(), self.polynomial(), self.gen(0).embedded())
    def __hash__(self):
        return NumberField_absolute.__hash__(self)

from sage.rings.number_field.number_field_element_quadratic import NumberFieldElement_quadratic

# global variable to controle the format of repr(RealNumberFieldElement_quadratic)
show_RNFElement_by_embedding = True

class RealNumberFieldElement_quadratic(NumberFieldElement_quadratic):
    
    def embedded(self):
        parent = self.parent()
        embedding = parent.coerce_embedding()
        e = embedding(self)
        return e
    
    def __repr__(self):
        global show_RNFElement_by_embedding
        if show_RNFElement_by_embedding:
            embedded = self.embedded()
            return '%s' % embedded
        else:
            return NumberFieldElement_quadratic._repr_(self)


    def __hash__(self):
        if not hasattr(self, '_hash'):
            self._hash = NumberFieldElement_quadratic.__hash__(self)
        return self._hash

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: from cutgeneratingfunctionology.igp import *
            sage: x, = nice_field_values([2^(1/2)])
            sage: sage_input(x)
            R.<y> = QQ[]
            RealNumberField(y^2 - 2, embedding=RR(1.4142135623730949), name='a')([0, 1])
            sage: sage_input((x, x))
            R.<y> = QQ[]
            K = RealNumberField(y^2 - 2, embedding=RR(1.4142135623730949), name='a')
            (K([0, 1]), K([0, 1]))

        """
        def maybe_ZZ(x):
            try:
                return ZZ(x)
            except (TypeError, ValueError):
                return x

        try:
            return sib(self.parent())(maybe_ZZ(QQ(self)))
        except (TypeError, ValueError):
            return sib(self.parent())([maybe_ZZ(x) for x in self.list()])

class RealNumberField_quadratic(NumberField_quadratic):
    def __init__(self, polynomial, name=None, latex_name=None, check=True, embedding=None,
                 assume_disc_small=False, maximize_at_primes=None, exact_embedding=None):
        r"""
        Create a real number field.
        """
        #### This is copy/paste from number_field.py, Sage 5.11,
        #### modified to change the element type.
        NumberField_quadratic.__init__(self, polynomial, name=name, check=check,
                                      embedding=embedding, latex_name=latex_name,
                                      assume_disc_small=assume_disc_small, maximize_at_primes=maximize_at_primes)
        #### Note: Modify "NumberField_absolute.__init__()" to "NumberField_quadratic.__init__()" as the former causes TypeError in Sage 5.12. --Yuan
        self._standard_embedding = True
        self._element_class = RealNumberFieldElement_quadratic
        c, b, a = [Rational(t) for t in self.defining_polynomial().list()]
        # set the generator
        Dpoly = b*b - 4*a*c
        D = (Dpoly.numer() * Dpoly.denom()).squarefree_part(bound=10000)
        self._D = D
        parts = -b/(2*a), (Dpoly/D).sqrt()/(2*a)
        self._NumberField_generic__gen = self._element_class(self, parts)

        # we must set the flag _standard_embedding *before* any element creation
        # Note that in the following code, no element is built.
        emb = self.coerce_embedding()
        if emb is not None:
            rootD = RealNumberFieldElement_quadratic(self, (QQ(0),QQ(1)))
            if D > 0:
                from sage.rings.real_double import RDF
                self._standard_embedding = RDF.has_coerce_map_from(self) and RDF(rootD) > 0
            else:
                raise DataError("RealNumberField_quadratic needs a positive discriminant")

        # we reset _NumberField_generic__gen has the flag standard_embedding
        # might be modified
        self._NumberField_generic__gen = self._element_class(self, parts)

        # NumberField_absolute.__init__(...) set _zero_element and
        # _one_element to NumberFieldElement_absolute values, which is
        # wrong (and dangerous; such elements can actually be used to
        # crash Sage: see #5316).  Overwrite them with correct values.
        #### Note: In my Sage 5.12, self(0) and self(1) return NumberFieldElement_quadratic. Should be RNFE_quadratic --Yuan
        #### self._zero_element = self(0) 
        #### self._one_element =  self(1)
        self._zero_element = RealNumberFieldElement_quadratic(self, QQ(0)) 
        self._one_element =  RealNumberFieldElement_quadratic(self, QQ(1))

    def _repr_(self):
        return "Real Number Field in `%s` as the root of the defining polynomial %s near %s"%(
                   self.variable_name(), self.polynomial(), self.gen(0).embedded())

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.
        """
        p = sib.name('RealNumberField')(self.polynomial(), embedding=RR(self.gen(0)), name='a')
        sib.id_cache(self, p, 'K')
        return p

# The factory.

def RealNumberField(polynomial, name=None, latex_name=None, names=None, check=True, embedding=None,
                    assume_disc_small=False, maximize_at_primes=None, exact_embedding=None):
    r"""
    A ``NumberField`` embedded into the real numbers, for which comparisons work according tothe embedding.

    This may not be necessary any more in Sage 7.1 after :trac:`17830`; but see :trac:`20184`.

    Some special tricks for speed and for pretty printing.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: x=polygen(QQ)
        sage: K.<cbrt2> = RealNumberField(x^3 - 2, embedding=RIF(AA.polynomial_root(x^3-2, RIF(0,3))), exact_embedding=AA.polynomial_root(x^3-2, RIF(0,3)))
        sage: 6064/4813 < cbrt2 < 90325/71691
        True

        sage: R.<y> = QQ[]; polynomial=y^2 - 2; embedding=RR(1.4142135623730949); exact_embedding=SR(AA(sqrt(2)))
        sage: RealNumberField(polynomial=polynomial, name='a', embedding=embedding, exact_embedding=exact_embedding)
        Real Number Field in `a` as the root of the defining polynomial y^2 - 2 near 1.414213562373095?


    """

    if names is not None:
        name = names
    if polynomial.degree() == 2:
        if embedding is not None:
            # In Sage 7.5, need to replace given embedding by exact root
            # before calling the constructors.
            try:
                from sage.rings.number_field.number_field_morphisms import root_from_approx
                embedding = root_from_approx(polynomial, embedding)
            except ImportError:
                pass

        K = RealNumberField_quadratic(polynomial, name, latex_name, check, embedding,
                                      assume_disc_small=assume_disc_small, 
                                      maximize_at_primes=maximize_at_primes, 
                                      exact_embedding=exact_embedding)
    else:
        K = RealNumberField_absolute(polynomial, name, latex_name, check, embedding,
                                     assume_disc_small=assume_disc_small, 
                                     maximize_at_primes=maximize_at_primes, 
                                     exact_embedding=exact_embedding)
    return K

def is_real_number_field_element(x):
    try:
        x.embedded # FIXME: this is a hack
        return True
    except AttributeError:
        return False

def is_all_the_same_real_number_field(values):
    real_number_field_seen = None
    for x in values:
        if is_real_number_field_element(x):
            if real_number_field_seen:
                if real_number_field_seen != x.parent():
                    return False, values
            else:
                real_number_field_seen = x.parent()
        elif not can_coerce_to_QQ(x):
            return False, values
    if real_number_field_seen:
        # Coerce any rationals to the number field.
        return True, [ real_number_field_seen(x) for x in values ]
    else:
        return False, values

