# Fluid.py file is used to contain the declaration of the classes associated to the
# different fluids present in the HydraulicModule. So far only gaseous fluids and
# Liquid are considered with a predefined standard behavior. This design decision has
# been taken on the basis that two different non-mixing fluids are interacting in the
# hydraulic network. To simplify the problem no multiphase fluid is considered.
# @author: Jose Felix Zapata Usandivaras
# Date: 13/11/2019
# ISAE - SUPAERO Space Section, Project "Griffon", Combustion Team

# -------------------------- IMPORT MODULES -------------------------------

from abc import ABC, abstractmethod             # Import abc components to define abstract classes

# ------------------------- FUNCTION DEFINITIONS ---------------------------

def calculate_liquid_density(pure_relative_density, purity):
    """
    calculate_liquid_density determines the liquid density based on purity and actual density
    of the pure component.
    :param pure_relative_density: relative density of the pure liquid (relative to water)
    :param purity: purity of the liquid (mixed with water).
    :return: density of the oxidizer [kg/m^3]
    """
    water_density = 1000  # Define the water density
    return water_density / (1 - purity * (1 - 1 / pure_relative_density))

def calculate_gas_density(pressure, temperature, gas_constant):
    """
    calculate_gas_density determines the density of the gas based on the ideal gas law
    :param pressure: pressure at which the gas is stored  [Pa]
    :param temperature: temperature at which the gas is stored [K]
    :param gas_constant: specific constant of the gas  [J/kgK]
    :return: gas density [kg/m^3]
    """
    return pressure / (temperature * gas_constant)


# -------------------------- CLASS DEFINITIONS -----------------------------


class Fluid(ABC):
    """
    Fluid class helps define the density properties of the different fluids that appear
    in the hydraulic system
    """
    def __init__(self, name, **kwargs):
        """
        class initializer
        :param name: string containing the name of the fluid
        """

        # Set the attributes
        self.name = name
        self.density = 0
        self.viscosity = 0

        # Set remaining properties
        for prop, value in kwargs.items():
            setattr(self, prop, value)

    def __str__(self):
        return "Fluid:: \n\t" + "\t\n".join(("{0}:{1}".format(name, value)
                                             for name, value in self.__dict__))

    def set_density(self, density):
        self.density = density

    def get_density(self):
        return self.density

    def set_viscosity(self, viscosity):
        self.viscosity = viscosity

    def get_viscosity(self):
        return self.viscosity

    @abstractmethod
    def compute_density(self, *args, **kwargs):
        """
        compute_density is an abstract method that can have multiple implementations
        depending on the inherited class, thus its multiple input definition.
        :param args: multiple arguments accepted
        :param kwargs: multiple key, value pair arguments accepted
        :return: nothing
        """
        pass


class Liquid(Fluid):
    """
    Liquid class inherits from the Fluid class, main characteristic
    is based on the procedure to calculate the density
    """
    def __init__(self, name, **kwargs):
        """ class initializer """
        # Call superclass initializer
        super(Liquid, self).__init__(name, **kwargs)
        # Compute the density at initialization
        self.compute_density()

    def compute_density(self, *args, **kwargs):
        """ compute the density of the liquid """
        self.density = calculate_liquid_density(self.pure_density, self.purity)


class Gas(Fluid):
    """
    Gas class inherits from the Fluid class, main characteristic is based on the procedure to
    calculate the density. Since the external pressure and temperatures at which the gas is
    operating are need to know its density, then in compute density this inputs are required.

    The ideal gas law is used to calculate the density of the gas
    """
    def __init__(self, name, **kwargs):
        """ class initializer """
        # Call superclass initializer
        super(Gas, self).__init__(name, **kwargs)

    def compute_density(self, *args, **kwargs):
        """ compute the density of the gas
        :param args: args[0] must be the pressure, args[1] the temperature
        :param kwargs: dictionary containing the keywords pressure, temperature
        :return nothing """
        if args:
            pressure, temperature = args[0], args[1]
        elif {'pressure', 'temperature'} <= set(kwargs):
            pressure, temperature = kwargs['pressure'], kwargs['temperature']
        else:
            raise ValueError("Illegal input error in compute density - Gas \n")

        # Compute the density
        self.density = calculate_gas_density(pressure, temperature, self.gas_constant)

    def compute_local_gas_temperature(self, pressure):
        """ compute_local_gas_temperature determines the instantaneous
        temperature of the gas given the pressure condition, and it's current density
         :return temperature in [K]"""
        return pressure / (self.gas_constant * self.density)


class FluidCatalogue:
    """ FluidCatalogue is a static class intended to serve as catalogue for
    the different types of fluid in place """
    catalogue = {'Liquid': Liquid, 'Gas': Gas}

    @staticmethod
    def return_initializer(fluid_type):
        return FluidCatalogue.catalogue[fluid_type]



