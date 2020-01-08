# SimulationObject.py includes the implementation of the Simulation object.
# This class is in charge of instantiating the Combustion and the MassSimulation
# modules and command their operation.
# Author: Jose Felix Zapata Usandivaras
# Date: 13/12/2018

# ------------------------ IMPORT MODULES -------------------------

from Initializer.Initializer import *                   # Import the Initializer module
import warnings                                         # Import the warnings module
from TrajectoryModule.Trajectory import *               # Import the Trajectory class
from Results.Results import *                           # Import Results classes
import CombustionModule.Combustion as Comb              # Import the Combustion module

# ---------------------- FUNCTION DEFINITIONS ---------------------


def prepare_inputs_for_trajectory_module(thrust_vector, isp_vector, simulation_time, dt):
    """
    prepare_inputs_for_trajectory_module rearranges the inputs in order to run the TrajectoryModule
    :param thrust_vector: thrust_vector as obtained from the combustion module
    :param isp_vector: isp_vector as obtained from the combustion module
    :param simulation_time: simulation_time [secs] for which we want the trajectory module to be run
    :param dt: time-spacing [secs] (same as used in combustion to make time-series coincident
    :return: time, thrust, isp,
    """

    # Get the time-array and complete the vectors
    time = dt * np.arange(start=0, stop=int(simulation_time / dt) + 1)

    # Generate thrust and isp
    thrust = np.zeros(time.shape)
    isp = np.zeros(time.shape)
    isp.fill(np.nan)

    if len(thrust_vector) < len(thrust):
        thrust[0:len(thrust_vector)] = thrust_vector
        isp[0:len(thrust_vector)] = isp_vector
    else:
        thrust = thrust_vector[0:len(thrust)]
        isp = isp_vector[0:len(thrust)]

    # Return the outputs
    return time, thrust, isp

# ----------------------- CLASS DEFINITIONS -----------------------


class SimulationObject:
    """
    SimulationObject is a master object which commands the operation of all the separate
    modules [Combustion, MassSimulator, Altitude].

    Attributes:
        1. initialization_object: Initializer instance with the data required for the simulation, it can be a list
        2. combustion_object: Combustion Module instance
        3. mass_simulator_object: MassSimulator instance
        4. altitude_object: Altitude instance
        5. results_dictionary: dictionary which contains the results of the code.
    """

    def __init__(self, initializer_collection=None):
        """
        class initializer
        :param initializer_collection: Initializer instance object used to set up the simulation
        """

        # Check the inputs, if initializer_obj is None then notify user that initializer must be set
        if isinstance(initializer_collection, InitializerCollection):

            # Set the value of the associated attribute
            self.initialization_collection = initializer_collection
            self.initialization_object = initializer_collection.elements_list[0]

        elif isinstance(initializer_collection, Initializer):

            # Create the collection and redo steps
            self.initialization_collection = InitializerCollection(initializer_obj=initializer_collection)
            self.initialization_object = initializer_collection

        elif initializer_collection is not None:

            # Raise a type error
            raise TypeError("Initializer collection must be a Initializer or InitializerCollection type. \n")

        else:
            # If the value passed for the initializer is None then notify user with a warning
            warnings.warn("Initialized simulation with None value \n ", UserWarning)
            self.initialization_collection = initializer_collection
            self.initialization_object = None

        # Set the remaining attributes
        self.combustion_module = self._generate_combustion_module()
        self.mass_simulator_module = None
        self.trajectory_module = self._generate_trajectory_module()

        # Initialize an empty ResultsCollection object
        self.results_collection = ResultsCollection([])

    def _generate_combustion_module(self):
        """
        _generate_combustion_module calls the constructor of the Combustion class and initializes it
        :return: Combustion object
        """

        # Check the validity of the initializer
        if self.initialization_object:
            # Initialize the combustion object if initializer is present
            obj = self.initialization_object.simulation_parameters['CombustionModel'](**self.initialization_object.combustion_parameters)
            # else:
                # obj = Comb.CombustionObjectClassic(**self.initialization_object.combustion_parameters)
        else:
            # Set the object to None value, issue warning
            warnings.warn("Set the CombustionObject to None, not representative object. \n")
            obj = None

        # Return the object
        return obj

    def update_mass_simulator_module(self):
        """
        _generate_mass_simulator_module calls the constructor of the MassSimulator class and initializes it
        :return: MassSimulator object
        """

        # Collect the object part
        object_class = self.initialization_object.simulation_parameters['mass_simulator']['system']

        # Check the validity of the initializer
        if self.initialization_object:
            # Initialize the Mass Simulation module if initializer is present
            obj = object_class(self.initialization_object.mass_simulator_parameters['system_dict'])
        else:
            # Set the object to None value, issue warning
            warnings.warn("Set the System object (Mass Simulator) to None, not representative object. \n")
            obj = None

        # Update the module
        self.mass_simulator_module = obj

    def _generate_trajectory_module(self):
        """
        _generate_trajectory_module calls the constructor of the TrajectoryModule class, this one without initialization
        :return: TrajectoryModule class instance with no initialization
        """

        # Check the validity of the initializer
        if self.initialization_object:
            # Initialize the Mass Simulation module if initializer is present
            obj = TrajectoryObject(**self.initialization_object.trajectory_parameters)

        else:
            # Set the object to None value, issue warning
            warnings.warn("Set the Trajectory Object to None, not representative object. \n")
            obj = None

        # Return the module
        return obj

    def update_modules(self, initialization_object):
        """
        update_modules re-runs the modules generators to update them to a new initializer object
        :param initialization_object: new initializer object
        :return: nothing
        """

        # Update the initialization_object
        self.initialization_object = initialization_object

        # Regenerate the modules
        self.combustion_module = self._generate_combustion_module()
        self.mass_simulator_module = None
        self.trajectory_module = self._generate_trajectory_module()

    def store_current_simulation_outputs(self):
        """
        store the current simulation outputs into the ResultsCollection attribute
        :return: nothing
        """

        # Instantiate a new Results object and add it to the collection
        result_obj = Results(initializer=self.initialization_object,
                             combustion_module=self.combustion_module,
                             mass_module=self.mass_simulator_module,
                             trajectory_module=self.trajectory_module)
        self.results_collection.add_element(result_obj)

    def export_results_to_file(self, file_name_expression):
        """
        export_results_to_file
        :return: nothing
        """
        if not self.results_collection.is_empty():
            self.results_collection.export_results_to_file(file_name=file_name_expression)
        else:
            print("\n-------- NO RESULTS DETECTED ON SIMULATION OBJECT ----------\n")

    def run_simulation_in_batch(self):
        """
        run_simulation_in_batch runs the simulations contained in the InitializerCollection and stores the output
        :return: nothing
        """

        # Run the simulations within a try-catch block
        for initializer_obj in self.initialization_collection.elements_list:

            # Update the simulation object to the initialization obj
            self.update_modules(initializer_obj)

            try:
                self.run_simulation()

            except (TypeError, ValueError, ArithmeticError, ZeroDivisionError, KeyError) as e:
                print("\n-----------ERROR ENCOUNTERED WHILE RUNNING THE SIMULATION---------\n")
                print(e)

    def run_simulation(self):
        """ run the simulation to obtain the data """

        # Start the simulation by calling the CombustionModule run_simulation method
        print("\n Running simulation on \n {initialization}".format(initialization=self.initialization_object))

        # ------------------------ RUN COMBUSTION MODULE ----------------------

        self.combustion_module.run_balanced_nozzle_analysis(
            **self.initialization_object.simulation_parameters['combustion'])

        # ------------------------ RUN MASS ESTIMATION ------------------------

        # Update the burn time
        if self.initialization_object.simulation_parameters['mass_simulator']['burn_time'] == 'TBD':
            self.initialization_object.simulation_parameters['mass_simulator']['burn_time'] = \
                self.combustion_module.results['magnitudes']['burn_time']

        # Update the mass simulator parameters
        self.initialization_object.update_mass_simulator_parameters()

        # Update the System object
        self.update_mass_simulator_module()

        # Calculate the mass
        initial_wet_mass = self.mass_simulator_module.get_mass()

        # Update the initial conditions for trajectory
        self.initialization_object.simulation_parameters['trajectory']['initial_conditions']['m0'] = initial_wet_mass

        # -------------------- PREPARE INPUTS FOR TRAJECTORY ------------------

        simulation_time = self.initialization_object.simulation_parameters['trajectory']['simulation_time']
        dt = self.initialization_object.simulation_parameters['combustion']['dt']
        initial_conditions = self.initialization_object.simulation_parameters['trajectory']['initial_conditions']

        # Get the results from combustion and unpack them
        combustion_results = self.combustion_module.return_results()
        thrust_vector = combustion_results['run_values']['thrust']
        isp = combustion_results['run_values']['isp']

        # Get time, thrust and isp
        time, thrust, isp = prepare_inputs_for_trajectory_module(thrust_vector, isp, simulation_time, dt)

        # ------------------------ RUN TRAJECTORY MODULE -------------------------

        self.trajectory_module.run_simulation_on_trajectory(time, thrust, isp, initial_conditions)

        # -------------------------- STORE THE RESULTS ----------------------------

        self.store_current_simulation_outputs()

        # print
        print("\n --- Simulation successfully performed \n")
