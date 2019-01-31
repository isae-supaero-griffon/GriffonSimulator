# Nozzle class for defining nozzle configuration
# Author: Maxime Sicat
# 28/12/2018

# ------------------------- IMPORT MODULES ----------------------

import math
import CombustionModule.Isentropic as iso

# -------------------------- CLASS DEFINITIONS ------------------


class Nozzle:
    """
    Nozzle class implements nozzle characteristics and evolution.

    Attributes:
        throat_area: current nozzle throat section area
        expansion_ratio: current nozzle expansion ratio (ie exit area divided by throat area)
        exit area: nozzle exit section area (set)
        lambda_efficiency: nozzle efficiency coefficient
        erosion_rate: fixed throat area regression rate
    """

    def __init__(self, At, expansion, lambda_e, erosion):
        """
        class initializer
        """
        self.throat_area = At
        self.expansion_ratio = expansion
        self.exit_area = expansion * At
        self.lambda_efficiency = lambda_e
        self.erosion_rate = erosion

    def __str__(self):
        """ return a string of the objects representation """

        # Return the objects type
        class_name = self.__class__.__name__
        # Loop through the objects properties and values
        return "\n{0} \n\t".format(class_name) + "\t\n\t".join(("{prop}, {val}".format(prop=prop,
                                                                val=value) for prop, value in vars(self).items()))

    @staticmethod
    def set_design_static(gamma, p_chamber, p_exit, c_star, isp, thrust):
        """
        Static method, return the main nozzle parameters for de given design conditions.
        :param gamma: flow isentropic parameter
        :param p_chamber: chamber pressure in Pascals
        :param p_exit: desired exit pressure (ie atmospheric pressure) in Pascals
        :param c_star: expected combustion c*
        :param isp: expected rocket specific impulse in seconds
        :param thrust: expected rocket thrust
        """

        # Determination of parameters based on isentropic flow and c* definition

        throat_area = c_star * thrust / isp / 9.81 / p_chamber

        exit_mach = iso.mach_via_pressure(gamma, p_chamber, p_exit)

        expansion_ratio = iso.expansion_via_mach(gamma, exit_mach)

        exit_area = throat_area * expansion_ratio

        return throat_area, expansion_ratio, exit_area

    def set_design(self, gamma, p_chamber, p_exit, c_star, isp, thrust):
        """
        Set the nozzle parameters as those defined by the previous static method. Inputs are identical.
        :param gamma: flow isentropic parameter
        :param p_chamber: chamber pressure in Pascals
        :param p_exit: desired exit pressure (ie atmospheric pressure) in Pascals
        :param c_star: expected combustion c*
        :param isp: expected rocket specific impulse in seconds
        :param thrust: expected rocket thrust
        """

        self.throat_area, self.expansion_ratio, self.exit_area = Nozzle.set_design_static(gamma, p_chamber, p_exit, c_star, isp, thrust)

    def get_throat_area(self):
        """
        Return the nozzle throat section area.
        """
        return self.throat_area

    def get_nozzle_effeciency(self):
        """
        Return the nozzle efficiency.
        """
        return self.lambda_efficiency

    def get_expansion_ratio(self):
        """
        Return the nozzle expansion ratio.
        """
        return self.expansion_ratio

    def get_exit_area(self):
        """
        Return the nozzle exit section area.
        """
        return self.exit_area

    def erode(self, dt):
        """
        Apply erosion to the throat
        :param dt: time increment
        """

        radius = math.sqrt(self.throat_area / math.pi)
        radius += self.erosion_rate*dt
        self.throat_area = math.pi * (radius**2)
        self.expansion_ratio = self.exit_area / self.throat_area



