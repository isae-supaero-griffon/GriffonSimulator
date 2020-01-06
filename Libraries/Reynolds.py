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
    :param x: value of colebrook variable
    :param roughness: pipe relative roughness
    :param hydraulic_diameter: hydraulic diameter of the pipe [m]
    :param reynolds_number: reynolds number of the flow
    :return: x0
    """
    a, b = 2.51 / reynolds_number, roughness / (3.7 * hydraulic_diameter)
    return - 2 * np.log(a * x + b)

def estimate_friction_factor(roughness, hydraulic_diameter):
    """
    estimate_friction_factor determines the value of the friction factor for high reynolds numbers
    :param roughness: pipe roughness
    :param hydraulic_diameter: hydraulic diameter of the pipe
    :return: estimate for high reynolds number of the friction factor
    """
    return  1 / (2 * np.log(3.7*hydraulic_diameter/roughness)) ** 2

def calculate_darcy_friction_factor(f0, roughness, hydraulic_diameter, reynolds_number):
    """
    calculate_darcy_friction_factor determines the friction factor along the pipe based
    on the colebrook equation.
    :param f0: initial guess for friction factor
    :param roughness: pipe relative roughness
    :param hydraulic_diameter: hydraulic diameter of the pipe [m]
    :param reynolds_number: reynolds number of the flow
    :return: darcy friction factor
    """

    # Define an initial value, error and tolerance
    x0 = f0
    # print("Re:{0} \n".format(reynolds_number))
    sol = fixed_point(colebrook_equation, [x0], args=(roughness, hydraulic_diameter, np.abs(reynolds_number)),
                      xtol=1e-4, maxiter=300)

    # Return result
    return 1 / sol[0] ** 2

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

    def __init__(self, roughness, hydraulic_diameter, friction_estimate = 0.0001):
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
        self.friction_estimate = friction_estimate

    def solve_for_friction_factor(self, density, viscosity, mass_flow):
        """
        solve_for_friction factor determines the friction factor of the pipe for the flow condition
        :param density: density of fluid in [kg/m^3]
        :param viscosity: viscosity of fluid
        :param mass_flow: mass flow of fluid flowing through the pipe [kg/sec]
        :return: darcy friction factor
        """
        # Check the inputs
        assert isinstance(fluid, Fluid), "fluid input must be a Fluid instance \n"
        flow_speed = calculate_flow_speed(mass_flow, self.area, density)
        r_number = calculate_reynolds_number(fluid.get_density(), viscosity, self.hydraulic_diameter,
                                             flow_speed)

        self.friction_estimate = calculate_darcy_friction_factor(self.friction_estimate, self.roughness,
                                                                 self.hydraulic_diameter, r_number)
        # Return the darcy friction factor
        return self.friction_estimate
