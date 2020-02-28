# RegressionModel.py allows for different hybrid fuel regression laws to be implemented.
# It defines the required model outputs through an abstract class definition,
# followed by several implementations.
# Author: Maxime Sicat / Jose Felix Zapata Usandivaras
# Date: 26/02/2019

# -------------------------- IMPORT MODULES -------------------------------

from abc import ABC, abstractmethod
import numpy as np

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
    def compute_regression_rate(self, geometry, ox_flow):
        """
        Return the instantenous solid fuel regression rate
        """
        pass

    @abstractmethod
    def compute_regression_rate_haltman(self, x, m_flux):
        """
        Return the instantaneous solid fuel regression rate
        :param x: position of the grain
        :param m_flux: total mass flow
        :return: regression rate
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
    r = a.(length^-0.2)(Go^n)

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

    def compute_regression_rate(self, geometry, ox_flow):
        """
        Return the instantenous solid fuel regression rate
        according to Marxman's equation
        """

        oxidizer_mass_flux = ox_flow / geometry.total_cross_section_area()

        return self.a * (geometry.get_length() ** self.m) * (oxidizer_mass_flux ** self.n)

    def compute_regression_rate_haltman(self, x, m_flux):
        """
        Return the instantaneous solid fuel regression rate
        according to Marxman's equation
        """
        if x != 0:
            return self.a * x ** self.m * m_flux ** self.n
        else:
            return 0


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

    def get_max_regression_rate(self, geometry):
        """
        Return the maximum regression rate possible
        within this model
        """
        return self.marxman.a * (geometry.get_length() ** self.marxman.m) * (self.threshold ** self.marxman.n)

    def get_max_regression_rate_local(self, x):
        """
        Return the maximum regression rate possible within this model
        :param x: x parameter along the grain
        :return: maximum regression rate
        """
        if x != 0:
            return self.marxman.a * (x ** self.marxman.m) * (self.threshold ** self.marxman.n)
        else:
            # TODO: Solve singularity
            return 0

    def compute_regression_rate(self, geometry, ox_flow):
        """
        Return the instantenous solid fuel regression rate
        depending on the current regime
        """

        oxidizer_mass_flux = ox_flow / geometry.total_cross_section_area()

        if oxidizer_mass_flux < self.threshold:
            return self.marxman.compute_regression_rate(geometry, ox_flow)
        else:
            return self.get_max_regression_rate(geometry)

    def compute_regression_rate_haltman(self, x, m_flux):
        """
        Return the instantaneous solid fuel regression rate
        depending on the current regime
        """

        if m_flux < self.threshold:
            return self.marxman.compute_regression_rate_haltman(x, m_flux)
        else:
            return self.get_max_regression_rate_local(x)


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

    def __init__(self, a, n, m, n2, maxGo, **kwargs):
        """
        class initializer
        Warning : it is insufficient by itself, as the
        make_continuous() has to be run before using this model.
        """

        super().__init__(**kwargs)

        self.marxman = SingleRegimeMarxmanModel(a, n, m)
        self.threshold = maxGo
        self.flooded = SingleRegimeMarxmanModel(0, n2, m)
        self.continuous = False

    def make_continuous(self):
        """
        Sets the factor in the flooded equation to assure
        continuity between the two regimes
        """

        a2 = self.marxman.a * (self.threshold ** (self.marxman.n - self.flooded.n))

        self.flooded.a = a2
        self.continuous = True

    def compute_regression_rate(self, geometry, ox_flow):
        """
        Return the instantenous solid fuel regression rate
        depending on the current regime
        :return regression rate
        """

        if not self.continuous:
            # Make the model continuous if it's not already continuous
            self.make_continuous()

        oxidizer_mass_flux = ox_flow / geometry.total_cross_section_area()

        if oxidizer_mass_flux < self.threshold:
            regression_rate = self.marxman.compute_regression_rate(geometry, ox_flow)
        else:
            regression_rate = self.flooded.compute_regression_rate(geometry, ox_flow)

        # Return the regression rate
        return regression_rate

    def compute_regression_rate_haltman(self, x, m_flux):
        """
        Return the instantaneous solid fuel regression rate depending on the current
        regime.
        :return regression rate
        """

        if not self.continuous:
            # Make the model continuous if it's not already continuous
            self.make_continuous()

        if m_flux < self.threshold:
            regression_rate = self.marxman.compute_regression_rate_haltman(x, m_flux)
        else:
            regression_rate = self.flooded.compute_regression_rate_haltman(x, m_flux)

        # Return the regression rate
        return regression_rate
