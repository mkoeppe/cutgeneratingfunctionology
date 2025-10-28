class FactorUndetermined(Exception): # FactorUndetermined is raised when an expression can not be evaluated with a test point.
    pass


class ParametricRealFieldFrozenError(ValueError):
    pass


class ParametricRealFieldInconsistencyError(ValueError):
    pass


class ParametricRealFieldRefinementError(ValueError):
    pass
