# Density.py contains the declaration of the DensityLaw class used in the calculation
# of the density variation with altitude used to asses the rockets performance in flight.
# Author: Jose Felix Zapata Usandivaras
# Date: 30/12/2018

# ------------------------ IMPORT MODULES -----------------------

import math                             # Import math library

# ----------------------- CLASS DEFINITIONS ---------------------


class DensityLaw:
    """ DensityLaw class is in charge of calculating the density change with altitude """

    def __init__(self, density_parameters):
        """
        class initializer
        :param density_parameters: dictionary containing the parameters required for density calculation
        """

        # Set the attributes [unpack the dictionary for speed gain]
        self.g0 = density_parameters['g0']
        self._r = density_parameters['R']
        self._m = density_parameters['M']
        self._rho0 = density_parameters['rho_0']
        self._t_0 = density_parameters['T_0']
        self.h_0 = density_parameters['h_0']

    def compute_density(self, altitude):
        """
        compute the density at a given altitude
        :param altitude: float-altitude at which we want the density to be calculated [m]
        :return: density [kg/m^3]
        """
        return self._rho0 * math.exp(-self.g0 * self._m * (altitude - self.h_0) / (self._r * self._t_0))
