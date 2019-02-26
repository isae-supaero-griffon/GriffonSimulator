# Combustion.py nucleates the combustion functionality in order to ease integration of results
# It performs the required integration of the Regression rate in order to obtain the results
# used for modelling of the rocket.
# Author: Maxime Sicat/Jose Felix Zapata Usandivaras
# Date: 29/12/2018

# -------------------------- IMPORT MODULES -------------------------------

from CombustionModule.Geometries import *                   # Import the geometry classes
from CombustionModule.Nozzle import *                       # Import the nozzle class
import CombustionModule.Isentropic as iso                   # Import the isoentropic module
import math                                                 # Import math
import matplotlib.pyplot as plt                             # Import matplotlib to plot results
from numpy import mean, trapz, asarray                      # Import mean from numpy for list arrays
from Libraries.Interpolator import *                        # Import interpolator
from DataLayer.JsonInterpreter import JsonInterpreter       # Import the JsonInterpreter



# -------------------------- CLASS DEFINITIONS -----------------------------

class CombustionObject:
    """
    The functionality of the combustion object is to perform the integration of the regression
    rate based on the parameters defined for the rocket engine.

    # Global attributes:
        simulation_type: simulation type used to interpolate from the CEA table

    # Attributes:
        1. data_dictionary: combustion table that contains all of the static data required for
        object to work properly, (mainly constants stored in the json file and uploaded by the
        JsonInterpreter.
        2. interpolator: Interpolator instance used to process CEA data table
        3. geometry: Geometry instance used to define the geometry of the port (hybrid rocket)
        4. nozzle: Nozzle instance to define the parameters of the nozzle and convert regression rate
        into thrust and Isp.
        5. results
    """

    # Global attribute
    simulation_type = "THEORETICAL ROCKET PERFORMANCE ASSUMING EQUILIBRIUM"

    def __init__(self, json_interpreter, geometry_object, nozzle_object):
        """
        class initializer
        :param json_interpreter: JsonInterpreter instance used to collect the data
        :param geometry_object: Geometry instance
        :param nozzle_object: Nozzle instance
        """

        # Assert nature of the inputs - JsonInterpreter is a singleton so it is instantiated once
        assert json_interpreter == JsonInterpreter.instance, "Please insert a valid JsonInterpreter. \n"
        assert isinstance(geometry_object, Geometry), " Please insert a valid Geometry instance. \n"
        assert isinstance(nozzle_object, Nozzle), " Please insert a valid Nozzle instance. \n"

        # Allocate the attributes
        self.data_dictionary = json_interpreter.return_combustion_table()
        self.interpolator = Interpolator(simulation_type=self.simulation_type,
                                         json_interpreter=json_interpreter)
        self.geometry = geometry_object
        self.nozzle = nozzle_object
        # TODO: convert the results lists to numpy array to be consistent with format all across the program.
        self.results = {
                        "run_values": {"time": [],
                                       "thrust": [],
                                       "isp": [],
                                       "pressure": [],
                                       "temperature": [],
                                       "radius": [],
                                       "of": []},
                        "magnitudes": {}
                       }

    def __str__(self):
        """ return a string that represents the Object """

        # Generate a string for the data dictionary
        data_str = "Combustion Table \n\t" + "\t\n\t".join(("{name}, {value}".format(name=variable, value=val) for
                                                            variable, val in self.data_dictionary.items()))

        # Generate a string for the geometry and nozzle
        inputs_str = "\n{geom} \n {nozzle}".format(geom=self.geometry, nozzle=self.nozzle)

        # Generate a string for the results
        if self.results['magnitudes']:
            results_str = "\nAverages values: \n\t" + "\t\n\t".\
                join(("{name}, {value:5.2f}".format(name=name, value=value) for
                      name, value in self.results['magnitudes'].items()))
        else:
            results_str = "\nNo results have been yet calculated. \n"

        # Return the string
        return "\n".join([data_str, inputs_str, results_str]) + "\n"

    def lookup_from_cea(self, desired_of_ratio):
        """
        lookup_from_cea go and collects the data from the CEA file contained in the interpolator object
        :param desired_of_ratio: float containing the desired O/F ratio used for the simulation
        :return: tuple with the variables of interest in their right locations ()
        """

        # Define the variables to interpolate
        vars_ = ["t", "m", "gammas", "cstar"]
        # Call out the interpolator with the required variables and get the output
        output = self.interpolator.interpolate_data(o_f_desired_value=desired_of_ratio,
                                                    variables=vars_,
                                                    mole_fractions="")

        # Return the values
        return output['variables']['t']['CHAMBER'], output['variables']['gammas']['CHAMBER'],\
               output['variables']['m']['CHAMBER']/1000, output['variables']['cstar']['THROAT']

    def unpack_data_dictionary(self):
        """ unpack the data dictionary into different variables to reduce code size
        order: (g0, R, pression_atmo, a, n, m, rho_fuel, initial_chamber_pressure) """

        return self.data_dictionary['g0'],          \
               self.data_dictionary['R'],           \
               self.data_dictionary['Pa'],          \
               self.data_dictionary['a'],           \
               self.data_dictionary['n'],           \
               self.data_dictionary['m'],           \
               self.data_dictionary['rho_fuel'],    \
               self.data_dictionary['P_chamber_bar'] * (10 ** 5)

    def post_process_data(self, dt):
        """ Compute results averages and total impulse"""
        self.results["magnitudes"] = {key: mean(value) for key, value in self.results["run_values"].items()
                                      if key != "time"}
        self.results["magnitudes"]["impulse"] = trapz(self.results["run_values"]["thrust"],
                                                      self.results["run_values"]["time"],
                                                      dt)

        self.results["magnitudes"]["burn_time"] = self.results["run_values"]["time"][-1]

    def return_results(self):
        """ return the results from the simulation """
        return self.results

    def run_simulation_constant_fuel_sliver(self, ox_flow, safety_thickness, dt, max_burn_time=None):
        """
        Simulate the hybrid rocket burn.
        Every value is expressed in ISU unless specified otherwise. Simulation is only valid for normal ground
        atmosphere conditions (1 bar, 293 K).
        :param ox_flow: value of input oxidizer flow
        :param safety_thickness: minimum allowed thickness for fuel slivers, which conditions burn termination
        :param dt: time increment
        :param max_burn_time: maximum burn time (if inputted)
        :return: TBD
        """

        # Check the inputs
        assert ox_flow > 0, "Oxidizer flow not greater than 0, check your inputs to CombustionModule. \n"
        assert safety_thickness > 0, "Safety thickness not greater than 0, check your inputs to CombustionModule. \n"
        assert dt > 0, "Simulation time-increment not greater than 0, check your inputs to CombustionModule. \n"

        # Call for workspace definition methods
        g0, R, pression_atmo, a, n, m, rho_fuel, initial_chamber_pressure = self.unpack_data_dictionary()

        # -------------------- Initialize simulation_variables if required -------------

        # regression_rate, fuel_flow, total_mass_flow, OF, T_chambre, gamma, masse_molaire_gaz, c_star = \
        #     self.initialize_variables()

        pression_chambre = initial_chamber_pressure

        # ---------------------------- MAIN SIMULATION LOOP -----------------------------

        # Set a counter to keep-track of the loop
        k = 0

        # Set the flag for the maximum burn-time
        flag_burn = True

        while self.geometry.min_bloc_thickness() > safety_thickness and flag_burn:

            # Regression rate and mass flow calculations
            fuel_flow = self.geometry.compute_fuel_rate(rho=rho_fuel, ox_flow=ox_flow)
            total_mass_flow = ox_flow + fuel_flow
            of_ratio = ox_flow / fuel_flow

            # Data extraction using the interpolator and OF ratio from already defined variables
            # Check lookup_from_cea method which variables have been chosen

            t_chamber, gamma, gas_molar_mass, cea_c_star = self.lookup_from_cea(desired_of_ratio=of_ratio)
            r = R / gas_molar_mass  # Specific gaz constant

            # Calculate chamber conditions and motor performance
            # Chamber pressure is not recalculated as it is assumed constant

            mach_exit = iso.exit_mach_via_expansion(gamma=gamma, expansion=self.nozzle.get_expansion_ratio(),
                                                    supersonic=True)
            exit_pressure = pression_chambre / iso.pressure_ratio(gamma, mach_exit)
            t_exit = t_chamber / iso.temperature_ratio(gamma, mach_exit)
            v_exit = math.sqrt(gamma * r * t_exit) * mach_exit

            thrust = self.nozzle.get_nozzle_effeciency() * (total_mass_flow*v_exit + (exit_pressure-pression_atmo) *
                                                            self.nozzle.get_exit_area())

            isp = thrust / total_mass_flow / g0

            # Results updates
            self.results['run_values']['time'].append(k * dt)
            self.results["run_values"]["thrust"].append(thrust)
            self.results["run_values"]["isp"].append(isp)
            self.results["run_values"]["pressure"].append(pression_chambre)
            self.results["run_values"]["temperature"].append(t_chamber)
            self.results["run_values"]["of"].append(of_ratio)

            # Verify it is a single port_number before updating the port number
            if isinstance(self.geometry, OneCircularPort):
                self.results["run_values"]["radius"].append(self.geometry.get_port_radius())
            else:
                # Append a 0 if there it is not a OneCircularPort geometry
                self.results["run_values"]["radius"].append(0)

            # Update the geometry and nozzle
            self.geometry.regress(ox_flow=ox_flow, dt=dt)
            self.nozzle.erode(dt)

            # Update the loop
            k += 1

            # Update the flag for burn-time
            if max_burn_time:
                flag_burn = k * dt <= max_burn_time

        # ------------------------ POST-PROCESS ------------------------

        # Convert to numpy arrays
        self.results['run_values'] = {key: asarray(value) for key, value in self.results["run_values"].items()}

        # Post-process the data
        self.post_process_data(dt=dt)

    def plot_results(self):
        """
        Plot results generates a plot with the variables typically outputted by Maxim's code
        :return: nothing
        """

        # Check the time array is not, if so then raise error, otherwise
        if len(self.results["run_values"]["time"]) != 0:

            # Extract the concerning results
            time = self.results["run_values"]["time"]
            thrust = self.results["run_values"]["thrust"]
            isp = self.results["run_values"]["isp"]
            radius = self.results["run_values"]["radius"]
            of_ratio = self.results["run_values"]["of"]
        else:
            raise ValueError("No values found for time, check results before plotting. \n")

        # Generate the plots
        fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(20, 15), squeeze=True, facecolor='w')
        fig.suptitle('Simulation results', fontsize=14)

        # Thrust-plot
        axs[0].plot(time, thrust, label='Thrust', color='blue')
        axs[0].set_title('')
        axs[0].set_xlabel('time (s)')
        axs[0].set_ylabel('Thrust (N)')
        axs[0].grid(b=True, axis='both')
        axs[0].set_xlim(left=time[0])

        # Port radius-plot
        axs[1].plot(time, radius, linestyle='--', label='Port Radius', color='red')
        axs[1].set_title('')
        axs[1].set_xlabel('time (s)')
        axs[1].set_ylabel('Port Radius (m)')
        axs[1].grid(b=True, axis='both')
        axs[1].set_xlim(left=time[0])

        # Isp-plot
        axs[2].plot(time, isp, label='Isp', color='r')
        axs[2].set_title('')
        axs[2].set_xlabel('time (s)')
        axs[2].set_ylabel('Isp (s)')
        axs[2].grid(b=True, axis='both')
        axs[2].set_xlim(left=time[0])

        # of_ratio-plot
        axs[3].plot(time, of_ratio, linestyle='--', color='green', label='O/F ratio')
        axs[3].set_title('')
        axs[3].set_xlabel('time (s)')
        axs[3].set_ylabel('O/F ratio')
        axs[3].grid(b=True, axis='both')
        axs[3].set_xlim(left=time[0])

        # Show the plot
        plt.show()

    @staticmethod
    def initialize_variables():
        """ initialize the variables that are used during the simulation and later recalculated
        as the simulation advances. These are:
            regression_rate, fuel_flow, total_mass_flow, OF, T_chambre, gamma, masse_molaire_gaz, c_star """
        return 0, 0, 0, 0, 3000, 1.4, 20*10**(-3), 1500
