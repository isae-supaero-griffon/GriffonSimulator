# Results.py contains the declaration of the Results class and the ResultsCollection class which will enable to handle
# in an easier way the results obtained from the simulations of the rocket engine.
# Author: Jose Felix Zapata Usandivaras
# Date: 3/1/2018

# ------------------------ IMPORT MODULES -----------------------

from Libraries.Collections import Collections                   # Import Collections class to define results collections
import csv                                                      # Import csv to allow the use of csv files
import datetime                                                 # Import the datetime to record the date
import numpy as np                                              # Import numpy to handle numpy arrays
import pandas as pd                                             # Import pandas to handle time-series
from Initializer.Initializer import Initializer                 # Import the Initializer class
from CombustionModule.Combustion import CombustionObject        # Import CombustionObject class
from MassEstimationModule.system import System                  # Import the System class
from TrajectoryModule.Trajectory import TrajectoryObject        # Import the Trajectory class

# ----------------------- FUNCTION DEFINITIONS -------------------


def merge_dictionaries(combustion_dictionary, trajectory_dictionary):
    """
    merge_dictionaries generates a list of dictionaries that is later used to generate the table,
    since not necessarily both dictionaries have the same length, then it is imperative that for
    merging they are adjusted in time-series length
    :param combustion_dictionary: dictionary containing the combustion results
    :param trajectory_dictionary: dictionary containing the trajectory results
    :return: list of dictionaries
    """

    # Compare time-series length
    combustion_len = len(combustion_dictionary['time'])
    trajectory_len = len(trajectory_dictionary['time'])

    if combustion_len < trajectory_len:
        aux_array = np.zeros((1, trajectory_len - combustion_len))
        combustion_dictionary = {key: np.append(initial_array, aux_array)
                                 for key, initial_array in combustion_dictionary.items()
                                 if key != 'time'}
    else:
        aux_array = np.zeros((1, combustion_len - trajectory_len))
        trajectory_dictionary = {key: np.append(initial_array, aux_array)
                                 for key, initial_array in trajectory_dictionary.items()
                                 if key != 'time'}

    # Merge dictionaries and sort keys
    output_dictionary = {**combustion_dictionary, **trajectory_dictionary}

    # Return the list of dictionaries
    return pd.DataFrame(output_dictionary).to_dict('records')

# ------------------------ CLASS DEFINITIONS ---------------------


class Results:
    """
    Results class is intended to facilitate the handling of the results outputted
    by the GriffonSimulator program.

    Its attributes are defined by the results given by each module of the GriffonSimulator:

    1. initializer_obj: InitializerObject instance which was used to set up the simulation
    2. combustion: Combustion module results (includes time evolving and averages)
    3. mass: Mass estimation module results (includes only the total mass)
    4. trajectory: Trajectory module results

    """

    def __init__(self, initializer, combustion_module, mass_module, trajectory_module):
        """
        class initializer
        :param initializer: Initializer object instance
        :param combustion_module: CombustionObject instance with results
        :param mass_module: System instance
        :param trajectory_module: TrajectoryObject instance
        """

        # Check in inputs
        assert isinstance(initializer, Initializer), "Please insert a valid Initializer object.\n"
        assert isinstance(combustion_module, CombustionObject), "Please insert a valid CombustionObject.\n"
        assert isinstance(mass_module, System), "Please insert a valid System object.\n"
        assert isinstance(trajectory_module, TrajectoryObject), "Please insert a valid TrajectoryObject.\n"

        # Set the attributes values
        self.initializer = initializer
        self.combustion = combustion_module
        self.mass_estimation = mass_module
        self.trajectory = trajectory_module

    # noinspection PyTypeChecker
    def export_results_to_file(self, file_name):
        """
        export the data to a .txt that is formatted
        :param file_name: name of the file where to output the files
        :return: nothing
        """

        # Open the file
        with open(file_name, 'w', newline='') as file:

            # Write the introduction of the file and the initializer description
            date_time = datetime.datetime.now().strftime("%I:%M%p on %B %d - %Y")
            intro_text = "\n".join(["Griffon Simulation Results - {date}".format(date=date_time),
                                    str(self.initializer),
                                    str(self.mass_estimation),
                                    str(self.combustion),
                                    str(self.trajectory), "\n"])

            # Write initializer preamble
            file.write(intro_text)

            # Return data
            combustion_results = self.combustion.return_results()['run_values']
            trajectory_results = self.trajectory.return_results()['run_values']

            # Write the columns
            columns_names = sorted(set(list(combustion_results.keys()) + list(trajectory_results.keys())))

            # Get the root of the file
            writer = csv.DictWriter(file, fieldnames=columns_names, delimiter=',')

            # write the header
            writer.writeheader()

            # Collect the table
            table = merge_dictionaries(combustion_dictionary=combustion_results,
                                       trajectory_dictionary=trajectory_results)
            # Loop through the table and output the result
            for row in table:
                writer.writerow(row)


class ResultsCollection(Collections):
    """
    ResultsCollection class works as a container of Results
    required to perform general operations to multiple Results
    Attributes:
        0. elements: list containing the constellation
    """

    def __init__(self, results_obj):
        """
        class initializer
        :param results_obj: Results instance, iterable of Results, or empty list
        """
        # Call superclass method
        super().__init__(results_obj, Results)

    def add_element(self, element, *args):
        """ Override method from parent class"""
        super().add_element(element, Results)

    def export_results_to_file(self, file_name):
        """
        export_results_to_file will call on every element of the collection the export_result_to_file method
        :param file_name: name of the file (generic expression)
        :return: nothing
        """

        # Loop through the list and execute the method
        count = 0
        for element in self.elements_list:
            count += 1
            # Execute the statement in a try-except block
            try:

                f_name = file_name.format(number=count)
                element.export_results_to_file(f_name)
                print("\nSuccessfully saved file: {name}".format(name=f_name))

            finally:
                pass

        # Print to screen that process has finished successfully
        print("-----------------  Process finished --------------------\n")
