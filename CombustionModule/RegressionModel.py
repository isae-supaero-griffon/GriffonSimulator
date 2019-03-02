# RegressionModel.py allows for different hybrid fuel regression laws to be implemented.
# It defines the required model outputs through an abstract class definition,
# followed by several implementations.
# Author: Maxime Sicat
# Date: 26/02/2019

# -------------------------- IMPORT MODULES -------------------------------

from abc import ABC, abstractmethod

# -------------------------- CLASS DEFINITIONS -----------------------------


class RegressionModel(ABC):

    """
    The functionality of the RegressionModel object is to encapsulate a model
    of solid fuel regression and its specific regimes.

    There are no master attributes for increased generality.
    """

    def __init__(self, **kwargs):
        """
        class initializer
        """

        super().__init__()

    @abstractmethod
    def computeRegressionRate(self, geometry, ox_flow):
        """
        Return the instantenous solid fuel regression rate
        """
        pass

    def __str__(self):
        """ return a string of the objects representation """

        # Return the objects type
        class_name = self.__class__.__name__
        # Loop through the objects properties and values
        return "{0}, ".format(class_name) + \
               ", ".join(("{prop}, {val}".format(prop=prop, val=value) for prop, value in vars(self).items()))


class SingleRegimeMarxmanModel(RegressionModel):

    """
    Classical regression model in which there exists a single
    regression regime, governed my Marxman's classical equation
    dependant on oxidizer mass flux through the combustion port:
    r = a.(L^-0.2)(Go^n)

    # Attributes:
        1. a: Constant factor dependant on fuel choice (experimental)
        2. n: Oxidizer mass flux exponent dependant on fuel choice (experimental, around 0.8)
        3. m: Length exponent dependant on fuel choice (experimental, around -0.2)
    """

    def __init__(self, a, n, m, **kwargs):
        """
        class initializer
        """

        super().__init__(**kwargs)

        self.a = a
        self.n = n
        self.m = m

    def computeRegressionRate(self, geometry, ox_flow):
        """
        Return the instantenous solid fuel regression rate
        according to Marxman's equation
        """

        oxidizer_mass_flux = ox_flow / geometry.totalCrossSectionArea()

        return self.a * (geometry.get_length() ** self.m) * (oxidizer_mass_flux ** self.n)


class MarxmanAndConstantFloodingRegimeModel(RegressionModel):

    """
    Model with two different regression regimes:
    Marxman's model below a threshold oxidizer mass flux and
    above it a flooding regime where regression rate is constant.

    # Attributes:
        1. marxman: SingleRegimeMarxmanModel object modeling the regression
        below the threshold
        2. threshold: Threshold oxidizer mass flux
    """

    def __init__(self, a, n, m, maxGo, **kwargs):
        """
        class initializer
        """

        super().__init__(**kwargs)

        self.marxman = SingleRegimeMarxmanModel(a, n, m)
        self.threshold = maxGo

    def getMaxRegressionRate(self, geometry):
        """
        Return the maximum regression rate possible
        within this model
        """
        return self.marxman.a * (geometry.get_length() ** self.marxman.m) * (self.threshold ** self.marxman.n)

    def computeRegressionRate(self, geometry, ox_flow):
        """
        Return the instantenous solid fuel regression rate
        depending on the current regime
        """

        oxidizer_mass_flux = ox_flow / geometry.totalCrossSectionArea()

        if oxidizer_mass_flux < self.threshold:
            return self.marxman.computeRegressionRate(geometry, ox_flow)
        else:
            return self.getMaxRegressionRate(geometry)


class TwoRegimesMarxmanAndFloodedModel(RegressionModel):

    """
    Model with two different regression regimes:
    Marxman's model with n around 0.8 below a threshold oxidizer mass flux and
    above it a flooding regime where regression rate has both an oxidizer flux
    and an implicit pressure dependence, as we assume constant pressure.

    # Attributes:
        1. marxman: SingleRegimeMarxmanModel object modeling the regression
        below the threshold
        2. threshold: Threshold oxidizer mass flux
        3. flooded: SingleRegimeMarxmanModel object modeling the flooded regime
        4. continuity: Boolean indicating whether the regression rate continuously
        goes from one regime to the other
    """

    def __init__(self, a1, n1, m1, n2, m2, maxGo, **kwargs):
        """
        class initializer
        Warning : it is insufficient by itself, as the
        makeContinuous() has to be run before using this model.
        """

        super().__init__(**kwargs)

        self.marxman = SingleRegimeMarxmanModel(a1, n1, m1)
        self.threshold = maxGo
        self.flooded = SingleRegimeMarxmanModel(0, n2, m2)
        self.continuous = False

    def makeContinuous(self, geometry):
        """
        Sets the factor in the flooded equation to assure
        continuity between the two regimes
        :param geometry: the equations are length-dependent
        so continuity has to take fuel geometry into account.
        """

        a2 = self.marxman.a * (geometry.get_length() ** (self.marxman.m - self.flooded.m)) * \
             (self.threshold ** self.marxman.n - self.flooded.n)

        self.flooded.a = a2
        self.continuous = True

    def computeRegressionRate(self, geometry, ox_flow):
        """
        Return the instantenous solid fuel regression rate
        depending on the current regime
        """

        if not self.continuous:

            print("WARNING model is not continuous, use at your own risk")

        oxidizer_mass_flux = ox_flow / geometry.totalCrossSectionArea()

        if oxidizer_mass_flux < self.threshold:

            return self.marxman.ComputeRegressionRate(geometry, ox_flow)

        else:

            return self.flooded.ComputeRegressionRate(geometry, ox_flow)
