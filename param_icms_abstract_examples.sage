## read_default_args crashes if there are no default args...
## construct_field_and_test_point expects to pass 'field' and 'conditioncheck' arguments -- hence the **args 
def test_posdef(a=5, b=7, **args):
    return Matrix([[a, b], [b, 1/4]]).is_positive_definite()

def find_region_type(field, result):
    # Return these label because that's what the plotting code expects.
    # TODO: Make mapping customizable
    if result:
        return 'is_extreme'
    else:
        return 'not_extreme'
    
complex = SemialgebraicComplex(test_posdef, ['a', 'b'], find_region_type=find_region_type)

complex.shoot_random_points(10)
plot(complex)

