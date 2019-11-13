# Monkey patch the Sage number field repr   (does not work)

from sage.rings.number_field.number_field_element import NumberFieldElement, NumberFieldElement_absolute
from sage.rings.number_field.number_field_element_quadratic import NumberFieldElement_quadratic

# global variable to control the format of repr
show_RNFElement_by_embedding = True

class _RNF_Patch:

    def _repr_(self):
        if show_RNFElement_by_embedding:
            return repr(AA(self))
        else:
            return NumberFieldElement._repr_(self)

for cls in NumberFieldElement, NumberFieldElement_quadratic, NumberFieldElement_absolute:
    cls._repr_ = _RNF_Patch._repr_    ## Well, doesn't work, that would have been too nice.

