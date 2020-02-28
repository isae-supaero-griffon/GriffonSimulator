# Fuel.py contains the declaration of the Fuel class which handles the communication with CEA through
# the RocketCEA library
# @ author: Jose Felix Zapata Usandivaras
# Date: 8/05/2019
# ISAE - SUPAERO, Scube GRIFFON project, Combustion Team

# ---------------------------- IMPORT MODULES ----------------------------

import rocketcea.cea_obj as cea_obj                     # Import RocketCEA
import numpy as np                                      # Import numpy
from DataLayer.JsonInterpreter import JsonInterpreter   # Import the JsonInterpreter class

# -------------------------- FUNCTION DEFINITIONS ------------------------


def bar2psia(a):
    """
    bar2psia converts from bar to psia
    :param a: float with pressure in bar
    :return: pressure in psia
    """
    return (a + 1.013) * 14.50377


def rankine2kelvin(a):
    """
    rankine2Kelvin converts from Rankine to Kelvin
    :param a: float with temperature in Rankine
    :return: temperature in Kelvin
    """
    return a / 1.8


def pascal2bar(a):
    """
        pascal2bar converts from pascals to bars
        :param a: float with pressure in pascals
        :return: pressure in bars
        """
    return a / 10 ** 5


def feet2meters(a):
    """
        feet2meters converts from feet to meters
        :param a: float with length/speed in feet
        :return: length/speed in meters
        """
    return a / 3.28083


def lbmfeet2kgm(a):
    """
    convert from lbm/ft^3 to kg/m^3
    :param a: float with lbm/ft^3
    :return: density in kg/m^3
    """
    return a * 16.0185


# ------------------------- CLASS DEFINITIONS ------------------------

class Fuel:
    """
    The class Fuel serves as wrapper of the RocketCEA module
    and obtain the thermochemical data associated to the fuel.
    All units handled are in SI, even though RocketCEA works with
    imperial units. Thus, proper conversion must be made in every output
    of RocketCEA.

    Attrs:
        cea_obj: CEA_Obj instance with the definition of fuel and oxydizer already implemented

    Note: The Fuel class is implemented as a singleton as the JsonInterpreter class.
    """

    # ----------------------- NESTED CLASS ---------------------------

    class __Fuel:
        """ Define the nested class for Fuel to implement it as a Singleton """

        def __init__(self, json_interpreter):
            """
            class initializer
            :param json_interpreter: JsonInterpreter instance from which the fuel data will be read
            """
            # Assert nature of the inputs - JsonInterpreter is a singleton so it is instantiated once
            assert json_interpreter == JsonInterpreter.instance, "Please insert a valid JsonInterpreter. \n"

            # Extract the propellants table
            propellant_table = json_interpreter.return_propellant_table()
            self.define_fuel(fuel_table=propellant_table['fuel'])
            self.define_oxidizer(oxidizer_table=propellant_table['oxidizer'])

            # Define the ispObj
            self.ispObj = cea_obj.CEA_Obj(oxName=propellant_table['oxidizer']['name'],
                                          fuelName=propellant_table['fuel']['name'])

        @staticmethod
        def define_fuel(fuel_table):
            """
            define_fuel constructs the propellant card (str) that's associated to
            the fuel.
            :param fuel_table: dict containing the fuel table as specified in the data-layer
            :return: Nothing
            """

            # Initialize the fuel card
            card = []

            # Loop through the fuels list and generate the string associated to it
            for fuel_type in fuel_table['components']:
                # Define the basic fields
                card.append("fuel={name}".format(name=fuel_type['fuel']))
                card.append("wt={wt}".format(wt=fuel_type['wt']))
                card.append("t{units}={value}".format(units=fuel_type['t']['units'], value=fuel_type['t']['value']))

                # Check if more complex fields are present for the fuel
                if "h" in fuel_type:
                    card.append("h,{units}={value}".format(units=fuel_type['h']['units'],
                                                           value=fuel_type['h']['value']))

                if "composition" in fuel_type:
                    card.append(" ".join(("{name} {value}".format(name=key, value=val)
                                          for key, val in fuel_type['composition'].items())))

            # Join the card
            card = " ".join(card)

            # Add the fuel to the cea_obj
            cea_obj.add_new_fuel(fuel_table['name'], card)

        @staticmethod
        def define_oxidizer(oxidizer_table):
            """
            define_oxidizer constructs the propellant card associated to the oxidizer
            :param oxidizer_table: dict containing the oxidizer table as specified in the data-layer
            :return: nothing
            """

            # Initialize the oxidizer card
            card = []

            # Loop through the oxidizer list and generate the string associated to it
            for ox_type in oxidizer_table['components']:
                # Define the basic fields
                card.append("oxid={name}".format(name=ox_type['oxid']))
                card.append("wt={wt}".format(wt=ox_type['wt']))
                card.append("t{units}={value}".format(units=ox_type['t']['units'], value=ox_type['t']['value']))

                # Check if more complex fields are present for the oxidizer
                if "h" in ox_type:
                    card.append("h,{units}={value}".format(units=ox_type['h']['units'],
                                                           value=ox_type['h']['value']))

                if "composition" in ox_type:
                    card.append(" ".join(("{name} {value}".format(name=key, value=val)
                                          for key, val in ox_type['composition'].items())))

            # Join the card
            card = " ".join(card)

            # Add the fuel to the cea_obj
            cea_obj.add_new_oxidizer(oxidizer_table['name'], card)

        def return_chamber_c_star(self, Pc, MR):
            """
            get_c_star_values returns the c_star values
            :param Pc: Chamber Pressure (in bars)
            :param MR: Mixture Ratio
            :return: float c_star value in SI units
            """

            # Get from the ispObj the c_star request
            return self.ispObj.get_Cstar(Pc=bar2psia(Pc), MR=MR)

        def return_combustion_variables(self, Pc, MR, eps):
            """
            return_combustion_variables returns the variables requested by
            the combustion module.
            :param Pc: Chamber Pressure (in bars)
            :param MR: Mixture Ratio
            :param eps: Expansion Ratio of the burn.
            :return: Tcomb, M, Gamma, Cstar, sonvel
            """

            # Call the process from cea_obj
            _, c_star, t_comb, mw, gam = self.ispObj.get_IvacCstrTc_ChmMwGam(Pc=bar2psia(Pc),
                                                                             MR=MR,
                                                                             eps=eps)

            # Get the sonvel
            sonvel = self.ispObj.get_Chamber_SonicVel(Pc=bar2psia(Pc),
                                                      MR=MR,
                                                      eps=eps)

            # Get the chamber density
            rho = self.ispObj.get_Chamber_Density(Pc=bar2psia(Pc), MR=MR, eps=eps)

            # Return the output
            return rankine2kelvin(t_comb), gam, mw/1000, feet2meters(c_star), \
                   feet2meters(sonvel), lbmfeet2kgm(rho)

    # ------------------- Instance:

    instance = None

    # ----------------------------- WRAPPER CLASS -----------------------------

    def __new__(cls, *args, **kwargs):
        """ __new__ always as a class method
        depends if we are instantiating the a new instance of we are calling an
        already created object and returning it """

        # Check if instance is None, if so, create the new object
        if not Fuel.instance:
            Fuel.instance = Fuel.__Fuel(*args, **kwargs)
        return Fuel.instance

    def __getattr__(self, name):
        return getattr(self.instance, name)

    def __instancecheck__(self, obj):
        return isinstance(obj, self.__Fuel)
