# GriffonModule.py includes the declaration of the GriffonModule abstract class. The
# role of the GriffonModule abstract class is to nucleate the behavior of the modules
# regarding results storing and their subsequent post-processing.
# @author Jose Felix Zapata Usandivaras
# Date: 28/02/2020
# ISAE-SUPAERO Space Section, Project "Griffon", Combustion Team

# -------------------------- IMPORT MODULES -------------------------------

from abc import ABC, abstractmethod                             # Import abstract classes module
import numpy as np                                              # Import numpy
import datetime                                                 # Import the datetime to record the date
import csv                                                      # Import csv to allow the use of csv files
import pandas as pd                                             # Import pandas to handle time-series


# -------------------------- CLASS DEFINITIONS -----------------------------


class Clock:
    """
    Clock class is a concrete class implemented as a Singleton whose's aim is to
    keep record of simulation time, thus simplifying time management.
    """

    # ---------------------------- NESTED CLASS --------------------------------

    class __Clock:
        """
        _Clock is the nested class used to implement it as a singleton.
        """
        def __init__(self, start_time=0.0, limit=None):
            """ nested class initializer
            :param start_time: float with start_time of the clock
            :param limit: limit if necessary to contrast with on the flag method
            """
            # Set the time attribute and the time_log
            self.time = start_time
            self.log = np.empty(shape=(0, 0))
            self.limit = limit

        def update(self, dt):
            """
            update the time by a time-step and append the value to the log
            :param dt: float with time-step
            :return: nothing
            """
            self._update_log()
            self.time += dt

        def _update_log(self):
            self.log = np.append(self.log, self.time)

        def reset(self):
            """
            reset the values of the clock to their original values
            :return: nothing
            """
            self.time = 0
            self.log = np.empty(shape=(0, 0))
            self.limit = None

        def get_time(self):
            """
            get_time returns the current time at the object
            :return: float with the recorded time
            """
            return self.time

        def set_time(self, new_time):
            """
            set_time sets the time value of the clock to a given time
            :param new_time: float with new_time in seconds
            :return: nothing
            """
            self.time = new_time

        def get_log(self):
            """
            get_log returns the history of the clock
            :return: np array with the log
            """
            return self.log

        def set_log(self, new_log):
            """
            set_log will set the history of the clock (used only in particular
            cases where time-resampling might be necessary)
            :param new_log: np.array
            :return: nothing
            """
            self.log = new_log

        def set_limit(self, limit):
            self.limit = limit

        def is_limit_reached(self):
            return self.limit <= self.time if self.limit is not None else False

    # ------------------- Instance:

    instance = None

    # ----------------------------- WRAPPER CLASS -----------------------------

    def __new__(cls, *args, **kwargs):
        """ __new__ always as a class method
        depends if we are instantiating the a new instance of we are calling an
        already created object and returning it """

        # Check if instance is None, if so, create the new object
        if not Clock.instance:
            Clock.instance = Clock.__Clock(*args, **kwargs)
        return Clock.instance

    def __getattr__(self, name):
        return getattr(self.instance, name)

    def __instancecheck__(self, obj):
        return isinstance(obj, self.__Clock)


class GriffonModule(ABC):
    """
    GriffonModule is an abstract class which serves as a base class for all of the
    modules implemented in the Griffon Simulator code. These modules are:
        1. Combustion Module
        2. Hydraulic Module
        3. Trajectory Module

    Mass simulator module is not considered at this point.

        Attributes:
            1. clock: single instance of Clock for simulation management
            2. results: dictionary containing the keys 'run_values' and 'magnitudes'
    """
    def __init__(self):
        """
        Class initializer
        """
        super(GriffonModule, self).__init__()
        # Set the attributes
        self.clock = Clock()
        self.results = {'time': [], 'run_values': {}, 'magnitudes': {}}

    def return_results(self):
        """ return the results from the simulation """
        return self.results

    def export_results_file(self, file_name):
        """
        export_results_file will output the data to the file specified
        :return: nothing
        """
        date_time = datetime.datetime.now().strftime("%I:%M%p on %B %d - %Y")
        intro_text = "\n".join(["{module} Results - {date}".format(module=type(self).__name__, date=date_time),
                                str(self), "\n"])
        columns_names = ['time'] + list(self.results['run_values'].keys())
        output_dictionary = dict(self.results['run_values'])
        output_dictionary.update(time=self.results['time'])
        data = pd.DataFrame(output_dictionary).to_dict('records')

        # Create the file and output the data
        with open(file_name, 'w', newline='') as file:
            # Write initializer preamble
            file.write(intro_text)
            # Get the root of the file
            writer = csv.DictWriter(file, fieldnames=columns_names, delimiter=',')
            # Write the header
            writer.writeheader()
            # Dump data
            for row in data:
                writer.writerow(row)

    @abstractmethod
    def _initialize_results(self):
        """
        _initialize_results is in charge of populating the run_values dictionary with the
        proper keys to be stored
        :return: nothing
        """
        pass

    def _populate_run_values(self, keys):
        """
        _populate_run_values will create the lists for run_values with their proper values
        of initialization.
        :param keys: iterable of keys
        :return: nothing
        """
        for field in keys:
            self.results['run_values'][field] = np.empty(shape=(0, 0))

    def _allocate_result(self, **kwargs):
        """
        _allocate_results is in charge of allocating the key,value pairs into run-values
        :param kwargs: key,value pair (key must be contained in run_values)
        :return: nothing
        """
        for key, value in kwargs.items():
            self.results['run_values'][key] = np.append(self.results['run_values'][key], value)

    def _calculate_maximums(self):
        """ determine the maximum values of run_values and allocate to magnitudes"""
        self.results["magnitudes"] = {"max_{key}".format(key=key): np.nanmax(value) for key, value
                                      in self.results['run_values'].items()}

    def _calculate_averages(self):
        """ determine the averages of run_values and allocate to magnitudes """
        # Calculate the
        self.results["magnitudes"] = {key: np.nanmean(value) for key, value in self.results["run_values"].items()}

    def post_process_data(self, *args):
        """
        post_process_data is a protected member of the GriffonModule class, it is in
        charge of handling the population of magnitudes
        :param args: might differ on the child class implementation (time-step, etc)
        :return: nothing
        """
        # Store the time and replace the initial value
        self.results['time'] = self.clock.get_log()

    def remap_run_results(self, new_time_array):
        """
        remap_run_results has for function to remap the results obtained by an asynchronous simulation
        to a one that is synchronous with a fixed time-step
        :param new_time_array: resampled time array for which to perform the interpolation
        :return: nothing
        """
        # Loop through the variables on run_results and remap by doing interpolation of the results
        # to increase the density of the array
        time_interp = self.clock.get_log()

        for key, value in self.results['run_values'].items():
            self.results['run_values'][key] = np.interp(new_time_array, time_interp, value)
