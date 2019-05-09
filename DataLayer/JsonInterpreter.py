# JsonInterpreter.py contains the declaration and implementation of the JsonInterpreter class which is designed to
# enable the extraction of static data concerning the Combustion and Mass simulation modules of Griffon project.
# To keep consistency and simplify the handling of the code, JsonInterpreter has been conceived as
# a singleton, that is, a class that can only be instantiated once in the lifespan of the program, all the
# other calls return the same object reference. The approach used to implement the pattern is the one of a nested class.
# Author: Jose Felix Zapata Usandivaras
# Date: 16/12/2018

# ------------------- IMPORT MODULES -----------------------

import json                     # Import the json module


# ------------------ CLASS DEFINITIONS ---------------------


class JsonInterpreter:
    """
    JsonInterpreter translates the json file which contains all of the static-data
    of the program into consumable data. Not designed to be used inside a loop.

    Attributes:
        1. json_data
    """

    # ---------------------------- NESTED CLASS --------------------------------

    class __JsonInterpreter:
        """ Define the nested class of JsonInterpreter to implement the singleton """

        def __init__(self, file_name):
            """ nested class initializer
             :param file_name: file_name of the json containing the data """

            # Load the file into memory
            with open(file_name, 'r') as f:
                json_str = f.read()
                self.json_data = json.loads(json_str)

            # Check whether the data is in proper form
            if not isinstance(self.json_data, dict):
                raise TypeError("Data loaded of not proper type, check your data. \n")

        def return_cea_data(self, simulation_type):
            """
            return_cea_data extracts the CEA
            :param simulation_type: str containing the simulation type
            :return: dict with requested data
            """
            # Get the right key from the dictionary and return it
            cea_table = self.json_data['CEA_Data']
            return cea_table[simulation_type]

        def return_combustion_table(self):
            """
            return_combustion_table returns the combustion table with all the static parameters
            :return: combustion table as a dict
            """
            return self.json_data['combustion_table']

        def return_mass_simulator_table(self):
            """
            return_mass_simulator_table returns the mass simulator table with all the static parameters
            :return: mass simulator table as a dict
            """
            return self.json_data['mass_simulator_table']

        def return_trajectory_table(self):
            """
            return_trajectory_table returns the trajectory module table with all the static parameters
            :return: trajectory module table as a dict
            """
            return self.json_data['trajectory_table']

        def return_propellant_table(self):
            """
            return_propellant_table returns the propellants data to be used to initialize the RocketCEA wrapper
            :return: propellant table as a dict
            """
            return self.json_data['propellant_table']

    # ------------------- Instance:

    instance = None

    # ----------------------------- WRAPPER CLASS -----------------------------

    def __new__(cls, *args, **kwargs):
        """ __new__ always as a class method
        depends if we are instantiating the a new instance of we are calling an
        already created object and returning it """

        # Check if instance is None, if so, create the new object
        if not JsonInterpreter.instance:
            JsonInterpreter.instance = JsonInterpreter.__JsonInterpreter(*args, **kwargs)
        return JsonInterpreter.instance

    def __getattr__(self, name):
        return getattr(self.instance, name)

    def __instancecheck__(self, obj):
        return isinstance(obj, self.__JsonInterpreter)
