# Combustion.py nucleates the combustion functionality in order to ease integration of results
# It performs the required integration of the Regression rate in order to obtain the results
# used for modelling of the rocket.
# Author: Maxime Sicat/Jose Felix Zapata Usandivaras
# Date: 29/12/2018

# -------------------------- IMPORT MODULES -------------------------------

from CombustionModule.Geometries import *                                           # Import the geometry classes
from CombustionModule.Nozzle import *                                               # Import the nozzle class
import CombustionModule.Isentropic as iso                                           # Import the isoentropic module
from CombustionModule.Fuel import *                                                 # Import the RocketCEA wrapper
import math                                                                         # Import math
import matplotlib.pyplot as plt                                                     # Import matplotlib to plot results
from matplotlib.ticker import FormatStrFormatter                                    # Import formatter
from numpy import mean, trapz, asarray, nan, nanmean, nanmax, linspace, nanmin      # Import numpy functions
from Libraries.Interpolator import *                                                # Import interpolator
from DataLayer.JsonInterpreter import JsonInterpreter                               # Import the JsonInterpreter
from abc import ABC, abstractmethod                                                 # Import the abstract class mold



# ------------------------- FUNCTION DEFINITIONS ---------------------------

def calculate_flow_speed(cross_section, mass_flow, density):
    """
    calculate_flow_speed determines the average speed of the port flow at the
    exit.
    :param cross_section: float with cross section area of the port [m^2]
    :param mass_flow: mass flow of mixture [kg/s]
    :param density: density of the flow [kg/m^3]
    :return: flow speed in m/s
    """
    return mass_flow / (density * cross_section)

# -------------------------- CLASS DEFINITIONS -----------------------------

class CombustionObject(ABC):
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
        super().__init__()

        # Assert nature of the inputs - JsonInterpreter is a singleton so it is instantiated once
        assert json_interpreter == JsonInterpreter.instance, "Please insert a valid JsonInterpreter. \n"
        assert isinstance(geometry_object, Geometry), " Please insert a valid Geometry instance. \n"
        assert isinstance(nozzle_object, Nozzle), " Please insert a valid Nozzle instance. \n"

        # Allocate the attributes
        self.data_dictionary = json_interpreter.return_combustion_table()
        self.geometry = geometry_object
        self.nozzle = nozzle_object
        self.fuel = Fuel(json_interpreter)
        self.results = {
                        "run_values": {"time": [0],
                                       "thrust": [0],
                                       "isp": [0],
                                       "pressure": [0],
                                       "temperature": [0],
                                       "radius": [self.geometry.get_port_radius() if
                                                  isinstance(self.geometry,
                                                             (OneCircularPort,
                                                              ThreeCircularPorts))
                                                  else nan],
                                       "regression_rate": [0],
                                       "mass_flow": [0],
                                       "of": [0],
                                       "Go": [0],
                                       "nozzle_param": [0],
                                       "c_star": [0],
                                       "chamber_sound_speed": [0],
                                       "chamber_rho": [0],
                                       "chamber_speed": [0],
                                       "hydraulic_port_diameter": [0]},
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
                join(("{name}, {value:.2E}".format(name=name, value=value) for
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
        vars_ = ["t", "m", "gammas", "cstar", "sonvel"]
        # Call out the interpolator with the required variables and get the output
        output = self.interpolator.interpolate_data(o_f_desired_value=desired_of_ratio,
                                                    variables=vars_,
                                                    mole_fractions="")

        # Return the values
        return output['variables']['t']['CHAMBER'], output['variables']['gammas']['CHAMBER'],\
               output['variables']['m']['CHAMBER']/1000, output['variables']['cstar']['THROAT'], \
               output['variables']['sonvel']['CHAMBER'],

    def run_thermochemical_analysis(self, of_ratio):
        """
        run_thermochemical_analysis returns the thermochemical variables at CHAMBER conditions as obtained from
        CEA.
        :param of_ratio: mixture ratio
        :return: tuple with the variables of interest in their right locations ()
        """

        # Determine the expansion ratio from Pamb and  Pc
        # eps = self.data_dictionary['P_chamber_bar'] / pascal2bar(self.data_dictionary['Pa'])
        eps = self.nozzle.expansion_ratio

        t_chamber, gamma, gas_molar_mass, cea_c_star, son_vel, rho_ch = self.fuel.return_combustion_variables(
            Pc=self.data_dictionary['P_chamber_bar'],
            MR=of_ratio,
            eps=eps)

        # Account for combustion efficiency
        eta_comb = self.data_dictionary['combustion_efficiency']            # Extract the combustion efficiency
        t_chamber = eta_comb * t_chamber                                    # Multiply the t_chamber by eta_comb
        rho_ch = rho_ch / eta_comb                                          # Divide the density (of gas)
        son_vel = m.sqrt(eta_comb) * son_vel

        # Return the output
        return t_chamber, gamma, gas_molar_mass, cea_c_star, son_vel, rho_ch

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
               self.data_dictionary['P_chamber_bar'] * (10 ** 5), \
               self.data_dictionary['combustion_efficiency']

    def post_process_data(self, dt):
        """ Compute results averages and total impulse"""

        for key in self.results['run_values'].keys():
            self.results['run_values'][key][0] = self.results['run_values'][key][1]

        self.results["magnitudes"] = {key: nanmean(value) for key, value in self.results["run_values"].items()
                                      if key != "time"}
        self.results["magnitudes"]["impulse"] = trapz(self.results["run_values"]["thrust"],
                                                                        self.results["run_values"]["time"],
                                                                        dt)

        self.results["magnitudes"]["burn_time"] = self.results["run_values"]["time"][-1]

    def return_results(self):
        """ return the results from the simulation """
        return self.results


    @abstractmethod
    def run_simulation_constant_fuel_sliver(self, ox_flow, safety_thickness, dt, max_burn_time):
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
        pass


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
            v_reg = self.results["run_values"]["regression_rate"]
            of_ratio = self.results["run_values"]["of"]
            Go = self.results["run_values"]["Go"]
            chamber_sound_speed = self.results["run_values"]["chamber_sound_speed"]
        else:
            raise ValueError("No values found for time, check results before plotting. \n")

        # Set the font dictionaries (for plot title and axis titles)
        title_font = {'size': '20', 'color': 'black', 'weight': 'normal',
                      'verticalalignment': 'bottom'}  # Bottom vertical alignment for more space
        axis_font = {'size': '16'}

        # Generate the plots
        fig = plt.figure(facecolor='w', figsize=(30, 30))
        fig.suptitle('Combustion Module results', **title_font)
        axs = [plt.subplot2grid((4, 2), (0, 0), rowspan=1, colspan=1),
               plt.subplot2grid((4, 2), (1, 0), rowspan=1, colspan=1),
               plt.subplot2grid((4, 2), (2, 0), rowspan=1, colspan=1),
               plt.subplot2grid((4, 2), (3, 0), rowspan=1, colspan=1),
               plt.subplot2grid((4, 2), (0, 1), rowspan=2, colspan=1),
               plt.subplot2grid((4, 2), (2, 1), rowspan=1, colspan=1),
               plt.subplot2grid((4, 2), (3, 1), rowspan=1, colspan=1)]

        # Set the tick labels font
        for ax in axs:
            for label in (ax.get_xticklabels() + ax.get_yticklabels()):
                label.set_fontname('Arial')
                label.set_fontsize(14)

        # Thrust-plot
        axs[0].plot(time, thrust, label='Thrust', color='blue', linewidth=2.0)
        axs[0].set_title('')
        axs[0].set_ylabel('Thrust (N)', **axis_font)
        axs[0].grid(b=True, axis='both')
        axs[0].set_xlim(left=time[0])
        axs[0].set_ylim(bottom=0, top=max(thrust)*1.5)

        # R-plot
        axs[1].plot(time, radius, linestyle='--', label='Radius', color='red')
        axs[1].set_title('')
        axs[1].set_ylabel('R (m)', **axis_font)
        axs[1].yaxis.set_major_formatter(FormatStrFormatter('%.2E'))
        axs[1].grid(b=True, axis='both')
        axs[1].set_xlim(left=time[0])

        # Isp-plot
        axs[2].plot(time, isp, label='Isp', color='r', linewidth=2.0)
        axs[2].set_title('')
        axs[2].set_ylabel('Isp (s)', **axis_font)
        axs[2].grid(b=True, axis='both')
        axs[2].set_xlim(left=time[0])
        axs[2].set_ylim(bottom=0, top=max(isp)*1.5)

        # of_ratio-plot
        axs[3].plot(time, of_ratio, linestyle='--', color='green', label='O/F ratio', linewidth=2.0)
        axs[3].set_title('')
        axs[3].set_xlabel('Time (s)', **axis_font)
        axs[3].set_ylabel('O/F ratio', **axis_font)
        axs[3].grid(b=True, axis='both')
        axs[3].set_xlim(left=time[0])

        # Regression rate plot
        axs[4].loglog(Go, v_reg, color='black', label='Regression Rate', linewidth=2.0)
        axs[4].set_title('')
        axs[4].set_ylabel('Regression Rate', **axis_font)
        axs[4].grid(True, which='both', ls='-')
        axs[4].set_ylim(bottom=1e-5, top=1e-2)
        axs[4].set_xlim(left=nanmin(Go), right=nanmax(Go))

        # Go-plot
        axs[5].plot(Go, time, linestyle='--', label='Go', color='red')
        axs[5].set_title('')
        axs[5].set_ylabel('Time (s)', **axis_font)
        axs[5].set_xlabel('Go [kg/m^2/sec]', **axis_font)
        axs[5].grid(True, which='both', ls='-')
        axs[5].set_xlim(left=nanmin(Go), right=nanmax(Go))

        # Sonic Speed-plot
        axs[6].plot(time, chamber_sound_speed, ls='-', label='Sound Speed', color='black')
        axs[6].set_title('')
        axs[6].set_ylabel('Chamber c (m/s^2)', **axis_font)
        axs[6].set_xlabel('time (secs)', **axis_font)
        axs[6].grid(True, which='both', ls='-')
        axs[6].set_ylim(bottom=nanmin(chamber_sound_speed), top=nanmax(chamber_sound_speed))


    @staticmethod
    def initialize_variables():
        """ initialize the variables that are used during the simulation and later recalculated
        as the simulation advances. These are:
            regression_rate, fuel_flow, total_mass_flow, OF, T_chambre, gamma, masse_molaire_gaz, c_star """
        return 0, 0, 0, 0, 3000, 1.4, 20*10**(-3), 1500



class CombustionObjectClassic(CombustionObject):

    def __init__(self, json_interpreter, geometry_object, nozzle_object):
        super().__init__(json_interpreter, geometry_object, nozzle_object)


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
        g0, R, pression_atmo, a, n, m, rho_fuel, initial_chamber_pressure, eta_comb = self.unpack_data_dictionary()

        # -------------------- Initialize simulation_variables if required -------------

        pression_chambre = initial_chamber_pressure

        # ---------------------------- MAIN SIMULATION LOOP -----------------------------

        # Set a counter to keep-track of the loop
        k = 1

        # Set the flag for the maximum burn-time
        flag_burn = True

        while self.geometry.min_bloc_thickness() > safety_thickness and flag_burn:

            # Regression rate and mass flow calculations
            fuel_flow = self.geometry.compute_fuel_rate(rho=rho_fuel, ox_flow=ox_flow)
            total_mass_flow = ox_flow + fuel_flow
            of_ratio = ox_flow / fuel_flow
            Go = ox_flow / self.geometry.totalCrossSectionArea()

            # Call CEA process to obtain thermochemical variables
            t_chamber, gamma, gas_molar_mass, cea_c_star, son_vel, rho_ch = self.run_thermochemical_analysis(
                of_ratio)
            r = R / gas_molar_mass  # Specific gaz constant

            # Determine the flow port-speed
            u_ch = calculate_flow_speed(cross_section=self.geometry.totalCrossSectionArea(),
                                        mass_flow=total_mass_flow,
                                        density=rho_ch)

            # Calculate chamber conditions and motor performance
            # Chamber pressure is not recalculated as it is assumed constant
            mach_exit = iso.exit_mach_via_expansion(gamma=gamma, expansion=self.nozzle.get_expansion_ratio(),
                                                    supersonic=True)
            exit_pressure = pression_chambre / iso.pressure_ratio(gamma, mach_exit)
            t_exit = t_chamber / iso.temperature_ratio(gamma, mach_exit)
            v_exit = math.sqrt(gamma * r * t_exit) * mach_exit

            thrust = self.nozzle.get_nozzle_effeciency() * (
            total_mass_flow * v_exit + (exit_pressure - pression_atmo) *
            self.nozzle.get_exit_area())

            # Calculate the nozzle equation validation
            nozzle_p = total_mass_flow * cea_c_star / (self.nozzle.get_throat_area() * initial_chamber_pressure)

            isp = thrust / total_mass_flow / g0

            # Results updates
            self.results['run_values']['time'].append(k * dt)
            self.results["run_values"]["thrust"].append(thrust)
            self.results["run_values"]["isp"].append(isp)
            self.results["run_values"]["pressure"].append(pression_chambre)
            self.results["run_values"]["temperature"].append(t_chamber)
            self.results["run_values"]["of"].append(of_ratio)
            self.results["run_values"]["Go"].append(Go)
            self.results["run_values"]["nozzle_param"].append(nozzle_p)
            self.results["run_values"]["c_star"].append(cea_c_star)
            self.results["run_values"]["chamber_sound_speed"].append(son_vel)
            self.results["run_values"]["chamber_rho"].append(rho_ch)
            self.results["run_values"]["mass_flow"].append(total_mass_flow)
            self.results["run_values"]["chamber_speed"].append(u_ch)
            self.results["run_values"]["hydraulic_port_diameter"].append(self.geometry.get_hydraulic_diameter())

            # TODO: modify the implementation of the radius storage
            # Verify it is a single port_number before updating the port number
            if isinstance(self.geometry, (OneCircularPort,
                                          ThreeCircularPorts)):
                self.results["run_values"]["radius"].append(self.geometry.get_port_radius())

                # Compute regression rate
                delta_r = self.results["run_values"]["radius"][k] - self.results["run_values"]["radius"][k - 1]
                self.results["run_values"]["regression_rate"].append(delta_r / dt)
            else:
                # Append a 0 if there it is not a OneCircularPort geometry
                self.results["run_values"]["radius"].append(nan)
                self.results["run_values"]["regression_rate"].append(nan)

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


class CombustionObjectImage(CombustionObject):

    def __init__(self, json_interpreter, geometry_object, nozzle_object):

        assert type(geometry_object).__name__ == 'SinglePortImageGeometry', "Geometry is not image-based, please " \
                                                                     "use SinglePortGeometry class. \n"

        super().__init__(json_interpreter, geometry_object, nozzle_object)

    def run_simulation_constant_fuel_sliver(self, ox_flow, safety_thickness, dt, max_burn_time):
        """
        Simulate the hybrid rocket burn.
        Every value is expressed in ISU unless specified otherwise. Simulation is only valid for normal ground
        atmosphere conditions (1 bar, 293 K).
        :param ox_flow: value of input oxidizer flow
        :param safety_thickness: minimum allowed thickness for fuel slivers, which conditions burn termination
        :param max_burn_time: maximum burn time (if inputted)
        :return: TBD
        """

        # Check the inputs
        assert ox_flow > 0, "Oxidizer flow not greater than 0, check your inputs to CombustionModule. \n"
        assert safety_thickness > 0, "Safety thickness not greater than 0, check your inputs to CombustionModule. \n"

        # Call for workspace definition methods
        g0, R, pression_atmo, a, n, m, rho_fuel, initial_chamber_pressure, eta_comb = self.unpack_data_dictionary()

        # -------------------- Initialize simulation_variables if required -------------

        # regression_rate, fuel_flow, total_mass_flow, OF, T_chambre, gamma, masse_molaire_gaz, c_star = \
        #     self.initialize_variables()
        pression_chambre = initial_chamber_pressure

        # ---------------------------- MAIN SIMULATION LOOP -----------------------------

        # Set the flag for the maximum burn-time
        flag_burn = True
        time = 0

        # Get the thickness equivalent to 5 pixels
        dr = self.geometry.getMetersPerPixel() * 5


        self.geometry.draw_geometry()

        while self.geometry.min_bloc_thickness() > safety_thickness and time < max_burn_time:

            # print( self.geometry.totalSurfaceArea() / self.geometry.get_length() / 2 / math.pi)

            # Regression rate and mass flow calculations
            fuel_flow = self.geometry.compute_fuel_rate(rho=rho_fuel, ox_flow=ox_flow)
            total_mass_flow = ox_flow + fuel_flow
            of_ratio = ox_flow / fuel_flow
            Go = ox_flow / self.geometry.totalCrossSectionArea()
            delta_t = dr / self.geometry.regressionModel.computeRegressionRate(self.geometry, ox_flow)

            # Call CEA process to obtain thermochemical variables
            t_chamber, gamma, gas_molar_mass, cea_c_star, son_vel, rho_ch = self.run_thermochemical_analysis(
                of_ratio)
            r = R / gas_molar_mass  # Specific gaz constant

            # Determine the flow port-speed
            u_ch = calculate_flow_speed(cross_section=self.geometry.totalCrossSectionArea(),
                                        mass_flow=total_mass_flow,
                                        density=rho_ch)

            # Calculate cea_c_star to consider the combustion efficiency parameter
            # cea_c_star = eta_comb * cea_c_star

            # Calculate chamber conditions and motor performance
            # Chamber pressure is not recalculated as it is assumed constant

            mach_exit = iso.exit_mach_via_expansion(gamma=gamma, expansion=self.nozzle.get_expansion_ratio(),
                                                    supersonic=True)
            exit_pressure = pression_chambre / iso.pressure_ratio(gamma, mach_exit)
            t_exit = t_chamber / iso.temperature_ratio(gamma, mach_exit)
            v_exit = math.sqrt(gamma * r * t_exit) * mach_exit

            thrust = self.nozzle.get_nozzle_effeciency() * (
            total_mass_flow * v_exit + (exit_pressure - pression_atmo) *
            self.nozzle.get_exit_area())

            # Calculate the nozzle equation validation
            nozzle_p = total_mass_flow * cea_c_star / (self.nozzle.get_throat_area() * initial_chamber_pressure)
            # print(nozzle_p)

            isp = thrust / total_mass_flow / g0

            # Results updates
            self.results['run_values']['time'].append(time)
            self.results["run_values"]["thrust"].append(thrust)
            self.results["run_values"]["isp"].append(isp)
            self.results["run_values"]["pressure"].append(pression_chambre)
            self.results["run_values"]["temperature"].append(t_chamber)
            self.results["run_values"]["of"].append(of_ratio)
            self.results["run_values"]["Go"].append(Go)
            self.results["run_values"]["nozzle_param"].append(nozzle_p)
            self.results["run_values"]["c_star"].append(cea_c_star)
            self.results["run_values"]["chamber_sound_speed"].append(son_vel)
            self.results["run_values"]["chamber_rho"].append(rho_ch)
            self.results["run_values"]["mass_flow"].append(total_mass_flow)
            self.results["run_values"]["chamber_speed"].append(u_ch)
            self.results["run_values"]["hydraulic_port_diameter"].append(self.geometry.get_hydraulic_diameter())

            # Verify it is a single port_number before updating the port number
            if isinstance(self.geometry, (OneCircularPort,
                                          ThreeCircularPorts,
                                          MultipleCircularPortsWithCircularCenter)):
                self.results["run_values"]["radius"].append(self.geometry.get_port_radius())

                # Compute regression rate
                self.results["run_values"]["regression_rate"].append(dr / delta_t)
            else:
                # Append a 0 if there it is not a OneCircularPort geometry
                self.results["run_values"]["radius"].append(nan)
                self.results["run_values"]["regression_rate"].append(nan)

            # Update the geometry and nozzle
            self.geometry.regress(ox_flow=ox_flow, dt=delta_t)
            self.nozzle.erode(delta_t)

            time += delta_t

            # Update the flag for burn-time
            flag_burn = time <= max_burn_time

        self.geometry.draw_geometry()

        # ------------------------ POST-PROCESS ------------------------

        # Convert to numpy arrays
        self.results['run_values'] = {key: asarray(value) for key, value in self.results["run_values"].items()}

        # Re-sample the data for dt
        k = 1
        time = []
        thrust = []
        isp = []
        pressure = []
        temperature = []
        of = []
        Go =[]
        nozzle_param = []
        c_star = []
        chamber_sound_speed = []
        chamber_rho = []
        mass_flow = []
        chamber_speed = []
        hydraulic_port_diameter = []
        radius = []
        regression_rate = []

        for i in range(0, int(self.results['run_values']['time'][-1] / dt)):
            t = dt * i
            time.append(t)

            thrust.append(self.results['run_values']['thrust'][k])
            isp.append(self.results['run_values']['isp'][k])
            pressure.append(self.results['run_values']['pressure'][k])
            temperature.append(self.results['run_values']['temperature'][k])
            of.append(self.results['run_values']['of'][k])
            Go.append(self.results['run_values']['Go'][k])
            nozzle_param.append(self.results['run_values']['nozzle_param'][k])
            c_star.append(self.results['run_values']['c_star'][k])
            chamber_sound_speed.append(self.results['run_values']['chamber_sound_speed'][k])
            chamber_rho.append(self.results['run_values']['chamber_rho'][k])
            mass_flow.append(self.results['run_values']['mass_flow'][k])
            chamber_speed.append(self.results['run_values']['chamber_speed'][k])
            hydraulic_port_diameter.append(self.results['run_values']['hydraulic_port_diameter'][k])
            radius.append(self.results['run_values']['radius'][k])
            regression_rate.append(self.results['run_values']['regression_rate'][k])

            if self.results['run_values']['time'][k] < t:
                k += 1

        time.append(time[-1] + dt)
        thrust.append(self.results['run_values']['thrust'][k])
        isp.append(self.results['run_values']['isp'][k])
        pressure.append(self.results['run_values']['pressure'][k])
        temperature.append(self.results['run_values']['temperature'][k])
        of.append(self.results['run_values']['of'][k])
        Go.append(self.results['run_values']['Go'][k])
        nozzle_param.append(self.results['run_values']['nozzle_param'][k])
        c_star.append(self.results['run_values']['c_star'][k])
        chamber_sound_speed.append(self.results['run_values']['chamber_sound_speed'][k])
        chamber_rho.append(self.results['run_values']['chamber_rho'][k])
        mass_flow.append(self.results['run_values']['mass_flow'][k])
        chamber_speed.append(self.results['run_values']['chamber_speed'][k])
        hydraulic_port_diameter.append(self.results['run_values']['hydraulic_port_diameter'][k])
        radius.append(self.results['run_values']['radius'][k])
        regression_rate.append(self.results['run_values']['regression_rate'][k])

        self.results['run_values']['time'] = time
        self.results["run_values"]["thrust"] = thrust
        self.results["run_values"]["isp"] = isp
        self.results["run_values"]["pressure"] = pressure
        self.results["run_values"]["temperature"] = temperature
        self.results["run_values"]["of"] = of
        self.results["run_values"]["Go"] = Go
        self.results["run_values"]["nozzle_param"] = nozzle_param
        self.results["run_values"]["c_star"] = c_star
        self.results["run_values"]["chamber_sound_speed"] = chamber_sound_speed
        self.results["run_values"]["chamber_rho"] = chamber_rho
        self.results["run_values"]["mass_flow"] = mass_flow
        self.results["run_values"]["chamber_speed"] = chamber_speed
        self.results["run_values"]["hydraulic_port_diameter"] = hydraulic_port_diameter
        self.results['run_values']['radius'] = radius
        self.results['run_values']['regression_rate'] = regression_rate

        # Post-process the data
        self.post_process_data(dt=dt)


