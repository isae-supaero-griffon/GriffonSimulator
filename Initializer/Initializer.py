# Initializer.py contains the Initializer class implementation in order to achieve the initialization of the different
# modules of the Griffon program, that is Combustion and MassSimulator module in principle.
# Author: Jose Felix Zapata Usandivaras
# Date: 12/12/2018


# ------------------- IMPORT MODULES ----------------------------

from math import pi, sqrt                                   # Import math functions
from CombustionModule.Geometries import *                   # Import the Geometries
from CombustionModule.Geometries1D import *                 # Import the 1D Geometries
from CombustionModule.Geometries0D import *                 # Import the 0D Geometries
from CombustionModule.Nozzle import *                       # Import the nozzles
from DataLayer.JsonInterpreter import JsonInterpreter       # Import the JsonInterpreter class
from TrajectoryModule.Drag import *                         # Import the Drag Module (DensityLaw in it)
from Libraries.Collections import Collections               # Import Collections class
from MassEstimationModule.system import *                   # Import the Mass Module classes
import CombustionModule.RegressionModel as Reg              # Import the RegressionModel module
from HydraulicModule.Hydraulic import HydraulicModule       # Import the Hydraulic module

# ------------------ FUNCTIONS/STATIC CLASSES -------------------


class InitializerCalculator:
    """
    InitializerCalculator is an auxiliary class which helps decoupling the different
    calculations present in the code and handle then separately to ease debugging.
    It acts as a library of methods/function. All the presented methods are static and it
    does not have and initializer.
    """

    @staticmethod
    def calculate_fuel_mass(geometry_obj, fuel_density):
        """ calculate_fuel_mass determines the mass of fuel contained in the combustion chamber.
        :param geometry_obj: Geometry instance used to calculate initial mass
        :param fuel_density: density of the fuel in the combustion chamber [kg/m^3] """

        # Check the inputs
        assert isinstance(geometry_obj, Geometry) or isinstance(geometry_obj, Geometry1D), \
            "Please insert a valid Geometry type. \n"

        # Return the output
        return geometry_obj.get_fuel_mass(fuel_density)

    @staticmethod
    def calculate_oxidiser_flow(oxidant_pressure, chamber_pressure, injection_loss, area_injection, ox_density):
        """ calculate_oxidiser_flow determines the flow of oxidizer based on chamber and oxidant_tank pressure.
        :param oxidant_pressure: pressure in [Pa] at the oxidant tank.
        :param chamber_pressure: pressure in [Pa] at the combustion chamber.
        :param injection_loss: coefficient for injection 0 < K* < 1, as a proportion of Delta P
        :param area_injection: area of injection [m^2]
        :param ox_density: density of the oxidiser [kg/m^3]
        :return oxidiser flow
        """

        # Return the oxidizer flow
        return area_injection * sqrt(2 * ox_density * (oxidant_pressure - chamber_pressure) * (1 - injection_loss))

    @staticmethod
    def calculate_oxidizer_tank_pressure(ox_flow, chamber_pressure, injection_loss, area_injection, ox_density):
        """
        calculate_oxidizer_tank_pressure determines the pressure at the oxidizer tank based on ox_flow, etc
        :param ox_flow: oxidizer flow based (float) [kg/sec]
        :param chamber_pressure: pressure in [Pa] at the combustion chamber.
        :param injection_loss: coefficient for injection 0 < K* < 1, as a proportion of Delta P
        :param area_injection: area of injection [m^2]
        :param ox_density: density of the oxidiser [kg/m^3]
        :return: pressure at the oxidiser tank
        """

        # Return the pressure at the oxidiser tank
        return chamber_pressure / (1 - injection_loss) + (ox_flow / area_injection) ** 2 / \
                                                         (2 * ox_density * (1 - injection_loss))

    @staticmethod
    def calculate_oxidiser_mass_based_on_tank(tank_radius, tank_height, oxidiser_density, tank_filling=1):
        """
        calculate_oxidiser_mass_based_on_tank determines the mass of oxidizer given the size of the tank and
        a tank_filling coefficient.
        :param tank_radius: radius of the tank
        :param tank_height: height of the tank
        :param oxidiser_density: density of the oxidizer
        :param tank_filling: % of the tank that is filled. It's default value is 1.
        :return: total oxidiser mass
        """

        # TODO: this issue should be contemplated in the part.py module with an appropiate class.
        return pi * tank_radius * tank_radius * tank_height * tank_filling * oxidiser_density

    @staticmethod
    def calculate_oxidiser_density(ox_density_relative, ox_purity):
        """
        calculate_oxidiser_density determines the oxidiser density based on purity and actual density
        of the pure component.
        :param ox_density_relative: relative density of the pure oxidizer (relative to water)
        :param ox_purity: purity of the oxidizer (mixed with water).
        :return: density of the oxidizer [kg/m^3]
        """
        water_density = 1000            # Define the water density
        return water_density / (1 - ox_purity * (1 - 1 / ox_density_relative))

    @staticmethod
    def calculate_oxidiser_mass_based_on_burn(ox_flow, burn_time, extra_filling):
        """
        calculate_oxidiser_mass_based_on_burn determines the total mass of oxidiser based on the size of the oxidiser
        tank.
        :param ox_flow: oxidiser flow defined as an input parameter of the simulation
        :param burn_time: simulation burn time as obtained from the combustion module
        :param extra_filling: percentage of extra filling for the tank with oxidiser
        :return: oxidiser mass in kg.
        """
        # Return the oxidiser mass
        return ox_flow * burn_time * (1 + extra_filling)

    @staticmethod
    def calculate_oxidiser_tank_height(ox_mass, ox_density, tank_radius):
        """
        calculate_oxidiser_tank_height determines the height of the oxidiser tank based on the total mass required
        to complete the burn of the rocket.
        :param ox_mass: required oxidiser mass in kg
        :param ox_density: density of the oxidiser [kg/m^3]
        :param tank_radius: proposed radius for the tank (external constraint)
        :return: tank height
        """
        # Return the calculated tank length
        return ox_mass / (ox_density * pi * tank_radius ** 2)


def pretty_print(dictionary):
    """ output the dictionary on a string format that is more convenient """
    if isinstance(dictionary, dict):
        return ", ".join(("{0}, {1}".format(key, pretty_print(value)) for key, value in dictionary.items()))
    else:
        return dictionary

# --------------------- CLASS DEFINITIONS -----------------------


class Initializer:
    """
    The Initializer class pretends to generate the inputs required by the combustion
    and the MassSimulator module.

    Attributes:
        1. combustion_parameters: store the calculated combustion_parameters dictionary
        2. mass_simulator_parameters: store the calculated mass_simulator_parameters
        3. trajectory_parameters: store the calculated trajectory_parameters
        4. simulation_parameters: store the simulation parameters in order to command the simulation run
    """

    def __init__(self, init_parameters, simulation_parameters, json_interpreter):
        """
        class initializer
        :param init_parameters: list of parameters defined to initialize the complete
        object as a hole.
        :param simulation_parameters: parameters associated to the run of the simulation
        :param json_interpreter: JsonInterpreter instance (singleton)
        """

        # Check the inputs
        assert isinstance(init_parameters, dict), "The parameters of the Simulation must be passed as a dict \"" \
                                                        "instance \n"
        assert json_interpreter == JsonInterpreter.instance, "Please insert a valid JsonInterpreter instance. \n"

        # Store the input as an attribute
        self.simulation_parameters = simulation_parameters

        self.combustion_parameters = self.generate_combustion_parameters(init_parameters['combustion'],
                                                                          json_interpreter)
        self.mass_simulator_parameters = self.generate_mass_simulator_parameters(json_interpreter)

        self.trajectory_parameters = self.generate_trajectory_parameters(json_interpreter)

    def __str__(self):
        """ return a string representation of the object based on the geometry, nozzle and simulation parameters
        implemented. """

        # Get the geometry and nozzle strings
        geometry_obj = self.combustion_parameters['geometry_object']
        nozzle_obj = self.combustion_parameters['nozzle_object']

        # Introduction string
        intro_str = "\nInitializing simulation w/parameters:"

        # Generate a string for the geometry and nozzle
        inputs_str = "\n{geom} \n {nozzle}".format(geom=geometry_obj, nozzle=nozzle_obj)

        # Simulation parameters string
        data_str = "\nSimulation Parameters: \n\t" + \
                   "\t\n\t".join(("{name}, {value}".format(name=variable, value=pretty_print(val)) for variable, val
                                  in self.simulation_parameters.items()))

        return "\n".join([intro_str, inputs_str, data_str])

    @staticmethod
    def generate_combustion_parameters(combustion_init_dict, json_interpreter):
        """
        generate_combustion_parameters performs all of the necessary calculations (if there are )required to set the
        inputs for the combustion module
        :param combustion_init_dict: dictionary containing the data required to prepare the inputs for the
        combustion module.
        :param json_interpreter: JsonInterpreter instance
        :return: combustion_parameters dictionary
        """

        # ------------- Generate the objects:

        # Geometry
        geometry_type = combustion_init_dict['geometric_params'].pop('type')
        regression_type = combustion_init_dict['geometric_params'].pop('regression_model')
        regression_obj = regression_type(**json_interpreter.return_combustion_table())
        combustion_init_dict['geometric_params']['regression_model'] = regression_obj
        geometry_obj = geometry_type(**combustion_init_dict['geometric_params'])
        if geometry_type == SinglePortImageGeometry:
            geometry_obj.generate_fourier(**combustion_init_dict['shape_params'])

        # Nozzle
        nozzle_obj = Nozzle(**combustion_init_dict['nozzle_params'])

        # Set the nozzle design
        if combustion_init_dict['set_nozzle_design']:
            nozzle_obj.set_design(**combustion_init_dict['design_params'])

        # Return output
        return {'json_interpreter': json_interpreter, 'geometry_object': geometry_obj, 'nozzle_object': nozzle_obj}

    def generate_mass_simulator_parameters(self, json_interpreter):
        """
        generate_mass_simulator_parameters simply outputs the system_dict
        :param json_interpreter: JsonInterpreter instance
        :return: mass estimation parameters dictionary
        """

        # Generate the output dictionary based on the class of mass simulator module to be implemented
        if self.simulation_parameters['mass_simulator']['system'] == SystemDynamic:
            output = {'system_dict': json_interpreter.return_mass_simulator_table(),
                      'combustion_table': json_interpreter.return_combustion_table()}
        elif self.simulation_parameters['mass_simulator']['system'] == SystemStatic:
            output = {'system_dict': self.simulation_parameters['mass_simulator']['system_dict'],
                      'combustion_table': json_interpreter.return_combustion_table()}
        else:
            raise TypeError("System type not recognized, please define a correct system type. \n")

        # Return the output
        return output

    def update_mass_simulator_parameters(self):
        """
        update_mass_simulator_parameters calss the appropiate method depending on the mass simulator class
        selected
        """
        # Generate the output dictionary based on the class of mass simulator module to be implemented
        if self.simulation_parameters['mass_simulator']['system'] == SystemDynamic:
            self.update_mass_simulator_parameters_dynamic()
        elif self.simulation_parameters['mass_simulator']['system'] == SystemStatic:
            self.update_mass_simulator_parameters_static()
        else:
            raise TypeError("System type not recognized, please define a correct system type. \n")

    def update_mass_simulator_parameters_static(self):
        """
        update_mass_simulator_parameters_static performs all of the necessary calculations (if there are) required to
        set the inputs for the mass simulation module when using a SystemStatic class for mass simulator
        """

        # Get the system dict
        system_dict = self.mass_simulator_parameters['system_dict']

        # Get combustion table
        combustion_table = self.mass_simulator_parameters['combustion_table']

        # ----------------------- Calculate propellant mass

        propellant_mass = InitializerCalculator.calculate_fuel_mass(self.combustion_parameters['geometry_object'],
                                                                    combustion_table['rho_fuel'])

        # ----------------------- Calculate oxidiser mass

        ox_flow = self.simulation_parameters['mass_simulator']['ox_flow']
        burn_time = self.simulation_parameters['mass_simulator']['burn_time']
        extra_filling = self.simulation_parameters['mass_simulator']['extra_filling']

        ox_mass = InitializerCalculator.calculate_oxidiser_mass_based_on_burn(ox_flow, burn_time, extra_filling)

        # ----------------------- Allocate the values into the dictionary

        system_dict['oxidizer_mass'] = ox_mass
        system_dict['fuel_mass'] = propellant_mass

    def update_mass_simulator_parameters_dynamic(self):
        """
        update_mass_simulator_parameters_dynamic performs all of the necessary calculations (if there are) required to
        set the inputs for the mass simulation module when using a SystemDynamic class for mass simulator
        """

        # Get the system dict
        system_dict = self.mass_simulator_parameters['system_dict']

        # Get combustion table
        combustion_table = self.mass_simulator_parameters['combustion_table']

        # ----------------------- Calculate propellant mass

        # chamber_length = self.combustion_parameters['geometry_object'].length
        # chamber_radius = self.combustion_parameters['geometry_object'].return_external_radius()
        propellant_mass = InitializerCalculator.calculate_fuel_mass(self.combustion_parameters['geometry_object'],
                                                                    combustion_table['rho_fuel'])

        # ----------------------- Calculate oxidiser mass

        ox_flow = self.simulation_parameters['mass_simulator']['ox_flow']
        # burn_time = self.simulation_parameters['mass_simulator']['burn_time']
        # extra_filling = self.simulation_parameters['mass_simulator']['extra_filling']
        #
        # ox_mass = InitializerCalculator.calculate_oxidiser_mass_based_on_burn(ox_flow, burn_time, extra_filling)

        # Calculate oxidizer mass based on the size of the tank
        ox_tank_radius = system_dict['subsystems']['oxidiser']['parts']['oxidant_tank']['radius']
        ox_tank_height = system_dict['subsystems']['oxidiser']['parts']['oxidant_tank']['height']

        ox_density = InitializerCalculator.calculate_oxidiser_density(combustion_table['rho_ox_pure'],
                                                                      combustion_table['ox_purity'])

        ox_mass = InitializerCalculator.calculate_oxidiser_mass_based_on_tank(tank_radius=ox_tank_radius,
                                                                              tank_height=ox_tank_height,
                                                                              oxidiser_density=ox_density)

        # -------------------------- Determine the max burn time based on oxidizer mass and allocated

        max_burn_time = ox_mass / ox_flow
        self.simulation_parameters['combustion']['max_burn_time'] = max_burn_time

        # ----------------------- Calculate oxidiser tank length

        # ox_density = InitializerCalculator.calculate_oxidiser_density(combustion_table['rho_ox_pure'],
        #                                                               combustion_table['ox_purity'])
        #
        # ox_tank_radius = system_dict['subsystems']['oxidiser']['parts']['oxidant_tank']['radius']
        #
        # # noinspection PyTypeChecker
        # ox_tank_height = InitializerCalculator.calculate_oxidiser_tank_height(ox_mass, ox_density, ox_tank_radius)

        # Calculate the oxidiser pressure in the tank and redefine the dictionary
        # chamber_pressure_mass = system_dict['subsystems']['combustion']['parts']['chamber']['mass']
        # chamber_pressure_combustion = combustion_table['P_chamber_bar'] * (10 ** 5)
        #
        # assert chamber_pressure_combustion == chamber_pressure_mass, "Chamber pressure in both systems are not the " \
        #                                                              "same. \n"

        # Injection loss
        injection_loss = self.simulation_parameters['mass_simulator']['injection_loss']

        # Area injection
        area_injection = self.simulation_parameters['mass_simulator']['area_injection']

        # Calculate the tank pressure
        # ox_pressure = InitializerCalculator.calculate_oxidizer_tank_pressure(ox_flow=ox_flow,
        #                                                                      chamber_pressure=chamber_pressure_mass,
        #                                                                      injection_loss=injection_loss,
        #                                                                      area_injection=area_injection,
        #                                                                      ox_density=ox_density)

        # Set the system_dict to its proper values
        # system_dict['subsystems']['combustion']['parts']['chamber']['height'] = chamber_length
        # system_dict['subsystems']['combustion']['parts']['chamber']['radius'] = chamber_radius
        system_dict['subsystems']['combustion']['parts']['chamber']['propellant_mass'] = propellant_mass
        system_dict['subsystems']['oxidiser']['parts']['oxidant_tank']['propellant_mass'] = ox_mass
        # system_dict['subsystems']['oxidiser']['parts']['oxidant_tank']['height'] = ox_tank_height
        # system_dict['subsystems']['oxidiser']['parts']['oxidant_tank']['pressure'] = ox_pressure

        # Update the system_dict
        self.mass_simulator_parameters['system_dict'] = system_dict

    @staticmethod
    def generate_trajectory_parameters(json_interpreter):
        """
        generate_trajectory_parameters generates the inputs for the TrajectoryModule
        :param json_interpreter: JsonInterpreter instance (singleton)
        :return: TrajectoryModule inputs dictionary
        """

        # TODO: review if possible to implement calculation of aerodynamic area of rocket based on system dict

        # Extract the trajectory data and generate the objects
        trajectory_data = json_interpreter.return_trajectory_table()
        density_obj = DensityLaw(trajectory_data['density'])

        drag_parameters = {'drag_coefficient': trajectory_data['drag']['Cd'],
                           'area_ref': trajectory_data['drag']['area_ref'],
                           'density': density_obj}
        drag_obj = DragCatalogue.return_drag_object(trajectory_data['drag']['type'], drag_parameters)

        # Return the dictionary
        return {'density_obj': density_obj, 'drag_obj': drag_obj}

    @staticmethod
    def state_simulation_parameters(simulation_parameters):
        """
        state_simulation_parameters helps determine the simulation parameters that are linked and calculate
        them in order to run the simulation. In the present case we are talking about oxidizer flow.
        :param simulation_parameters: simulation_parameters dictionary
        :return: modified simulation_parameters dictionary
        """

        # Return the dictionary
        return simulation_parameters


class InitializerCollection(Collections):
    """
    InitializerCollection class works as a container of Initializer
    required to perform general operations to multiple Initializer
    Attributes:
        0. elements: list containing the constellation
    """

    def __init__(self, initializer_obj):
        """
        class initializer
        :param initializer_obj: Initializer instance, iterable of Initializer, or empty list
        """
        # Call superclass method
        super(InitializerCollection, self).__init__(initializer_obj, Initializer)

    def add_element(self, element, *args):
        """ Override method from parent class"""
        super().add_element(element, Initializer)