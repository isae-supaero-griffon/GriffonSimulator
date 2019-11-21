# Reynolds.py helps determine the pressure loss along the pipeline by means of the colebrook equation.
# @author: Jose Felix Zapata Usandivaras
# Date: 21/11/2019
# ISAE-SUPAERO, Space Section Project "Griffon", Combustion team

# --------------------- IMPORT MODULES ----------------------

import numpy as np                          # Import numpy
from HydraulicModule.Fluid import Fluid     # Import the Fluid class
from scipy.optimize import fixed_point      # Import fixed point solution

# ------------------------- FUNCTION DEFINITIONS ------------

def calculate_flow_speed(mass_flow, area, rho):
    """
    calculate_flow_speed determines the speed of the flow along the pipe
    :param mass_flow: mass flow of fluid going through the pipe [kg/sec]
    :param area: area of the pipe [m^2]
    :param rho: density of the fluid [kg/m^3]
    :return: flow speed [m/sec]
    """
    return mass_flow / (area * rho)

def calculate_reynolds_number(rho, viscosity, hydraulic_diameter, speed):
    """
    calculate_reynolds_number determines the value of the Re number based on the
    diameter
    :param rho: density of the fluid [kg/m^3]
    :param viscosity: dynamic viscosity of the fluid [Pa.sec]
    :param hydraulic_diameter: hydraulic diameter of the pipe [m]
    :param speed: speed of the flow [m/s]
    :return: Reynolds number
    """
    return rho * speed * hydraulic_diameter / viscosity

def colebrook_equation(x, roughness, hydraulic_diameter, reynolds_number):
    """
    colebrook_equation is the equation implemented for fixed point algorithm
    in the form (fun(x0) = x0)
    :param roughness: pipe relative roughness
    :param hydraulic_diameter: hydraulic diameter of the pipe [m]
    :param reynolds_number: reynolds number of the flow
    :return: x0
    """
    a, b = 2.51 / reynolds_number, roughness / (3.7 * hydraulic_diameter)
    return - 2 * np.log(a * x + b)


def calculate_darcy_friction_factor(roughness, hydraulic_diameter, reynolds_number):
    """
    calculate_darcy_friction_factor determines the friction factor along the pipe based
    on the colebrook equation.
    :param roughness: pipe relative roughness
    :param hydraulic_diameter: hydraulic diameter of the pipe [m]
    :param reynolds_number: reynolds number of the flow
    :return: darcy friction factor
    """

    # Define an initial value, error and tolerance
    x0 = 0.000001
    sol = fixed_point(colebrook_equation, [x0], args=(roughness, hydraulic_diameter, reynolds_number),
                      xtol=1e-4, maxiter=300)

    # Return result
    return 1 / sol ** 2

# --------------------- CLASSES DEFINITIONS -----------------

class Reynolds:
    """
    Reynolds class is in charge of providing a means to calculate the darcy friction factor
    for a given fluid.

        # Attributes:
            1. roughness: relative roughness of the pipe
            2. hydraulic_diameter: hydraulic diameter of the pipe
            3. area: area of the pipe (to be considered circular pipe)

    """

    def __init__(self, roughness, hydraulic_diameter):
        """
        class initializer
        :param roughness: relative roughness of the pipe
        :param hydraulic_diameter: hydraulic diameter of the pipe
        """
        # Check the inputs
        assert 0 <= roughness, "Pipe roughness must be greater or equal than 0 \n"
        assert 0 < hydraulic_diameter, "Pipe diameter has to be greater than 0 \n"

        # Set the attributes
        self.hydraulic_diameter = hydraulic_diameter
        self.roughness = roughness
        self.area = np.pi * self.hydraulic_diameter ** 2 / 4

    def solve_for_friction_factor(self, fluid, mass_flow):
        """
        solve_for_friction factor determines the friction factor of the pipe for the flow condition
        :param fluid: Fluid object of interest
        :param mass_flow: mass flow of fluid flowing through the pipe [kg/sec]
        :return: darcy friction factor
        """
        # Check the inputs
        assert isinstance(fluid, Fluid), "fluid input must be a Fluid instance \n"
        flow_speed = calculate_flow_speed(mass_flow, self.area, fluid.get_density())
        r_number = calculate_reynolds_number(fluid.get_density(), fluid.get_viscosity(), self.hydraulic_diameter,
                                             flow_speed)
        # Return the darcy friction factor
        return calculate_darcy_friction_factor(self.roughness, self.hydraulic_diameter, r_number)
