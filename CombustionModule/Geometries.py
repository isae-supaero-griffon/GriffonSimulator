# Geometry abstract class and its various implementations for hybrid rocket burn calculations
# Author: Maxime Sicat
# 28/12/2018

# ------------------------- IMPORT MODULES ----------------------

from abc import ABC, abstractmethod                 # Import abc module for abstract classes definitions
import matplotlib.pyplot as plt                     # Import matplotlib
import math as m                                    # Import the math module
import numpy as np                                  # Import numpy
import cv2                                          # Import cv2
from CombustionModule.RegressionModel import *      # Import the regression models
from scipy.integrate import solve_ivp               # Import solve_ivp method from scipy
import os.path                                      # Import the OS library
from Libraries.Interpolator import Interpolator     # Import the interpolator

# ------------------------ FUNCTION DEFINITIONS -----------------


# -------------------------- CLASS DEFINITIONS ------------------

class Geometry(ABC):
    """
    Geometry abstract class sets the common methods and implementation
    for all geometries meant to be tested in the Griffon hybrid burn simulation
    code.

    Attributes:
        1. length: Length of the fuel slab
        2. regression_model: RegressionModel object for regression
        3. r_ext: external radius
        rate calculation
    """

    def __init__(self, length, regression_model, r_ext):
        """
        class initializer
        :param length: float with length of the grain [m]
        :param regression_model: RegressionModel instance
        :param r_ext: float with external radius [m]
        """
        # Check the inputs
        assert length > 0, "Grain length (length) must be greater than 0 \n"
        assert r_ext > 0, "Grain external radius (r_ext) must be greater than 0 \n"
        assert isinstance(regression_model, RegressionModel), "Please insert a valid RegressionModel instance. \n"

        # Call superclass initializer
        super(Geometry, self).__init__()

        # Set the attributes
        self.length = length
        self.regression_model = regression_model
        self.r_ext = r_ext

    def __str__(self):
        """ return a string of the objects representation """

        # Return the objects type
        class_name = self.__class__.__name__
        # Loop through the objects properties and values
        return "Geometry" + "\t\n\ttype, {0}\t\n\t".format(class_name) + \
               "\t\n\t".join(("{prop}, {val}".format(prop=prop, val=value) for prop, value in vars(self).items()))

    def get_length(self):
        """
        Return the fuel slab length.
        """
        return self.length

    def return_external_radius(self):
        """
        Return external radius of the geometry
        :return: external radius of the geometry
        """
        return self.r_ext

    @abstractmethod
    def compute_fuel_rate(self, rho, ox_flow):
        """
        Return the instantaneous fuel mass flow rate for the geometry
        :param rho: fuel density
        :param ox_flow: oxidizer mass flow rate
        :return float with mass flow rate of fuel
        """
        pass

    @abstractmethod
    def solve_mass_flux(self, m_ox, rho_f):
        """
        solve_mass_flux solves for the mass flux of oxidizer and fuel
        :param m_ox: mass flow of oxidizer [kg/s]
        :param rho_f: fuel density in [kg/m^3]
        :return: solution_vector, mass_flows tupple, of_ratio, mass_fluxes tupple
        """
        pass

    @abstractmethod
    def regress(self, ox_flow, dt, *args):
        """
        Apply regression to the geometry
        :param ox_flow: instantaneous oxidizer flow
        :param dt: time increment
        :param args: additional arguments (depending on the class)
        """
        pass

    @abstractmethod
    def min_bloc_thickness(self):
        """
        Return the smallest fuel space between ports or fuel edge.
        Useful for setting an burn termination safety margin.
        """
        pass

    @abstractmethod
    def get_fuel_mass(self, fuel_density):
        """
        Return the initial fuel mass based on fuel_density
        :param fuel_density: density of the fuel [kg/m^3]
        """
        pass

    @abstractmethod
    def total_cross_section_area(self):
        """
        return the total_cross_section area of the grain. For 1D case it shall be the area at the end of the grain
        :return: total cross section area at the end of the grain
        """
