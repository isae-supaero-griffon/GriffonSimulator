class Material:
    """

    """
    def __init__(self, identifier, density=None, yield_stress=None, *args):
        self.name = identifier
        self.density = density
        self.yield_stress = yield_stress

