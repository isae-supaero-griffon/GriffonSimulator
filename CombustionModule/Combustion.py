# Combustion.py nucleates the combustion functionality in order to ease integration of results
# It performs the required integration of the Regression rate in order to obtain the results
# used for modelling of the rocket.
# Author: Maxime Sicat/Jose Felix Zapata Usandivaras
# Date: 29/12/2018

# -------------------------- IMPORT MODULES -------------------------------

from CombustionModule.Geometries import Geometry                                # Import the geometry abstract class
from CombustionModule.Geometries0D import *                                     # Import the Geometry0D classes
from CombustionModule.Geometries1D import *                                     # Import the Geometry1D classes
from CombustionModule.Nozzle import *                                           # Import the nozzle class
import CombustionModule.Isentropic as iso                                       # Import the isoentropic module
from CombustionModule.Fuel import *                                             # Import the RocketCEA wrapper
from CombustionModule.RegressionModel import *                                  # Import the Regression Model
import math                                                                     # Import math
import matplotlib.pyplot as plt                                                 # Import matplotlib to plot results
from matplotlib.ticker import FormatStrFormatter                                # Import formatter
from numpy import mean, trapz, asarray, nan, nanmean, nanmax, linspace, nanmin  # Import numpy functions
from Libraries.Interpolator import Interpolator                                 # Import interpolator
from DataLayer.JsonInterpreter import JsonInterpreter                           # Import the JsonInterpreter
from abc import ABC, abstractmethod                                             # Import the abstract class mold
from HydraulicModule.Hydraulic import *                                         # Import the Hydraulic Module

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

    # Attributes:
        1. data_dictionary: combustion table that contains all of the static data required for
        object to work properly, (mainly constants stored in the json file and uploaded by the
        JsonInterpreter.
        2. fuel: Fuel object instance, acts as a wrapper of RocketCEA module
        2. nozzle: Nozzle instance to define the parameters of the nozzle and convert regression rate
        into thrust and Isp.
        3. results
    """
    # TODO: include a nested class for simulation management so as to avoid duplicated code fragments

    def __init__(self, json_interpreter, geometry_object, nozzle_object):
        """
        class initializer
        :param json_interpreter: JsonInterpreter instance used to collect the data
        :param geometry_object: Geometry instance used as Geometry
        :param nozzle_object: Nozzle instance

        """
        # Call superclass initializer
        super(CombustionObject, self).__init__()

        # Assert nature of the inputs - JsonInterpreter is a singleton so it is instantiated once
        assert json_interpreter == JsonInterpreter.instance, "Please insert a valid JsonInterpreter. \n"
        assert isinstance(geometry_object, Geometry), "Please insert a valid Geometry instance. \n"
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
                           "mass_flow": [0],
                           "mass_flow_ox": [0],
                           "mass_flow_f": [0],
                           "of": [0],
                           "Go": [0],
                           "nozzle_param": [0],
                           "c_star": [0],
                           "chamber_sound_speed": [0],
                           "chamber_rho": [0],
                           "chamber_speed": [0]},
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
            results_str = "\nAverages values: \n\t" + "\t\n\t". \
                join(("{name}, {value:.2E}".format(name=name, value=value) for
                      name, value in self.results['magnitudes'].items()))
        else:
            results_str = "\nNo results have been yet calculated. \n"

        # Return the string
        return "\n".join([data_str, inputs_str, results_str]) + "\n"

    def _print_last_solution_step(self):
        """
        _print_last_solution_step will print the solution over the last step stored in the results dictionary
        :return: nothing
        """
        # Configure the string
        vars_out = ['time', 'thrust', 'isp', 'of', 'pressure', 'mass_flow_ox']
        vals = [self.results['run_values'][k][-1] for k in vars_out]
        out = ", ".join(("{0:>5s}:{1:8.3f}".format(k, val) for k, val in zip(vars_out, vals)))
        print(out)

    def _run_thermochemical_analysis(self, of_ratio, chamber_pressure):
        """
        _run_thermochemical_analysis returns the thermochemical variables at CHAMBER conditions as obtained from
        CEA.
        :param of_ratio: mixture ratio
        :param chamber_pressure: pressure of the chamber in bars
        :return: tuple with the variables of interest in their right locations ()
        """

        # Determine the expansion ratio from Pamb and  Pc
        eps = self.nozzle.expansion_ratio
        t_chamber, gamma, gas_molar_mass, cea_c_star, son_vel, rho_ch = self.fuel.return_combustion_variables(
            Pc=chamber_pressure,
            MR=of_ratio,
            eps=eps)

        # Account for combustion efficiency
        eta_comb = self.data_dictionary['combustion_efficiency']  # Extract the combustion efficiency
        t_chamber = eta_comb * t_chamber  # Multiply the t_chamber by eta_comb
        rho_ch = rho_ch / eta_comb  # Divide the density (of gas)
        son_vel = m.sqrt(eta_comb) * son_vel

        # Return the output
        return t_chamber, gamma, gas_molar_mass, cea_c_star, son_vel, rho_ch

    def _unpack_data_dictionary(self):
        """ unpack the data dictionary into different variables to reduce code size
        order: (g0, R, pression_atmo, a, n, m, rho_fuel, initial_chamber_pressure) """

        return self.data_dictionary['g0'], \
               self.data_dictionary['R'], \
               self.data_dictionary['Pa'], \
               self.data_dictionary['a'], \
               self.data_dictionary['n'], \
               self.data_dictionary['m'], \
               self.data_dictionary['rho_fuel'], \
               self.data_dictionary['P_chamber_bar'] * (10 ** 5), \
               self.data_dictionary['combustion_efficiency']

    def _post_process_data(self, dt):
        """ Compute results averages and total impulse"""

        for key in self.results['run_values'].keys():
            self.results['run_values'][key][0] = self.results['run_values'][key][1]

        self.results["magnitudes"] = {key: nanmean(value) for key, value in self.results["run_values"].items()
                                      if key != "time"}
        self.results["magnitudes"]["impulse"] = trapz(self.results["run_values"]["thrust"],
                                                      self.results["run_values"]["time"],
                                                      dt)

        self.results["magnitudes"]["burn_time"] = self.results["run_values"]["time"][-1]

        # Calculate global O/F
        total_fuel_mass = trapz(self.results["run_values"]["mass_flow_f"], self.results["run_values"]["time"], dt)
        total_ox_mass = trapz(self.results["run_values"]["mass_flow_ox"], self.results["run_values"]["time"], dt)
        self.results["magnitudes"]["global_of"] = total_ox_mass / total_fuel_mass

    def return_results(self):
        """ return the results from the simulation """
        return self.results

    def plot_results(self):
        """
        Plot results generates a plot with the relevant variables of the code
        :return: nothing
        """
        # Check the time array is not, if so then raise error, otherwise
        if len(self.results["run_values"]["time"]) != 0:
            # Extract the concerning results
            time = self.results["run_values"]["time"]
            thrust = self.results["run_values"]["thrust"]
            isp = self.results["run_values"]["isp"]
            of_ratio = self.results["run_values"]["of"]
            g_ox = self.results["run_values"]["Go"]
        else:
            raise ValueError("No values found for time, check results before plotting. \n")

        # Set the font dictionaries (for plot title and axis titles)
        title_font = {'size': '20', 'color': 'black', 'weight': 'normal',
                      'verticalalignment': 'bottom'}  # Bottom vertical alignment for more space
        axis_font = {'size': '16'}

        # Generate the plots
        fig = plt.figure(facecolor='w', figsize=(30, 30))
        fig.suptitle('Combustion Module results', **title_font)
        axs = [plt.subplot2grid((2, 2), (0, 0), rowspan=1, colspan=1),
               plt.subplot2grid((2, 2), (1, 0), rowspan=1, colspan=1),
               plt.subplot2grid((2, 2), (0, 1), rowspan=1, colspan=1),
               plt.subplot2grid((2, 2), (1, 1), rowspan=1, colspan=1)]

        # Set the tick labels font
        for ax in axs:
            for label in (ax.get_xticklabels() + ax.get_yticklabels()):
                label.set_fontname('Arial')
                label.set_fontsize(14)

        # Thrust-plot
        axs[0].plot(time, thrust, label='Thrust', color='blue', linewidth=4.0)
        axs[0].set_title('')
        axs[0].set_ylabel('Thrust (N)', **axis_font)
        axs[0].grid(b=True, axis='both')
        axs[0].set_xlim(left=time[0])
        axs[0].set_ylim(bottom=0, top=max(thrust) * 1.5)

        # Isp-plot
        axs[1].plot(time, isp, label='Isp', color='r', linewidth=4.0)
        axs[1].set_title('')
        axs[1].set_ylabel('Isp (s)', **axis_font)
        axs[1].grid(b=True, axis='both')
        axs[1].set_xlim(left=time[0])
        axs[1].set_ylim(bottom=0, top=max(isp) * 1.5)

        # of_ratio-plot
        axs[2].plot(time, of_ratio, linestyle='--', color='green', label='O/F ratio', linewidth=4.0)
        axs[2].set_title('')
        axs[2].set_xlabel('Time (s)', **axis_font)
        axs[2].set_ylabel('O/F ratio', **axis_font)
        axs[2].grid(b=True, axis='both')
        axs[2].set_xlim(left=time[0])

        # Go-plot
        axs[3].plot(time, g_ox, linestyle='--', label='Go', color='red', linewidth=4.0)
        axs[3].set_title('')
        axs[3].set_xlabel('Time (s)', **axis_font)
        axs[3].set_ylabel('Go [kg/m^2/sec]', **axis_font)
        axs[3].grid(True, which='both', ls='-')
        axs[3].set_xlim(left=time[0])

    def run_simulation_constant_fuel_sliver(self, ox_flow, safety_thickness, dt, max_burn_time=None, **kwargs):
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
        g0, R, pression_atmo, a, n, m, rho_fuel, initial_chamber_pressure, eta_comb = self._unpack_data_dictionary()

        # -------------------- Initialize simulation_variables if required -------------
        pression_chambre = initial_chamber_pressure

        # ---------------------------- MAIN SIMULATION LOOP -----------------------------
        # Set a counter to keep-track of the loop
        k = 1

        # Set the flag for the maximum burn-time
        flag_burn = True

        while self.geometry.min_bloc_thickness() > safety_thickness and flag_burn:
            # ------------------- Perform the regression of the block ---------------------
            # noinspection PyArgumentList
            solution, mass_flows, of_ratio, mass_fluxes = self.geometry.solve_mass_flux(ox_flow, rho_fuel)
            m_ox, m_fuel = mass_flows
            g_ox, g_f = mass_fluxes
            total_mass_flow = m_ox + m_fuel

            # --------------------- Run Thermochemical Analysis ---------------------------
            t_chamber, gamma, gas_molar_mass, cea_c_star, son_vel, rho_ch = self._run_thermochemical_analysis(
                of_ratio, pression_chambre / 10 ** 5)
            r = R / gas_molar_mass  # Specific gaz constant

            # Determine the flow port exit-speed
            u_ch = calculate_flow_speed(cross_section=self.geometry.mesh.cells[-1].return_area_data(),
                                        mass_flow=total_mass_flow,
                                        density=rho_ch)

            # ------------------------- Determine Propulsion Performances ------------------
            mach_exit = iso.exit_mach_via_expansion(gamma=gamma, expansion=self.nozzle.get_expansion_ratio(),
                                                    supersonic=True)
            exit_pressure = pression_chambre / iso.pressure_ratio(gamma, mach_exit)
            t_exit = t_chamber / iso.temperature_ratio(gamma, mach_exit)
            v_exit = math.sqrt(gamma * r * t_exit) * mach_exit

            thrust = self.nozzle.get_nozzle_effeciency() * (
                    total_mass_flow * v_exit + (exit_pressure - pression_atmo) *
                    self.nozzle.get_exit_area())

            # Calculate the nozzle equation validation
            nozzle_p = total_mass_flow * cea_c_star / (self.nozzle.get_throat_area() * pression_chambre)
            isp = thrust / total_mass_flow / g0

            # ------------------------------ Results updates --------------------------------
            self.results['run_values']['time'].append(k * dt)
            self.results["run_values"]["thrust"].append(thrust)
            self.results["run_values"]["isp"].append(isp)
            self.results["run_values"]["pressure"].append(initial_chamber_pressure)
            self.results["run_values"]["temperature"].append(t_chamber)
            self.results["run_values"]["of"].append(of_ratio)
            self.results["run_values"]["Go"].append(g_ox)
            self.results["run_values"]["nozzle_param"].append(nozzle_p)
            self.results["run_values"]["c_star"].append(cea_c_star)
            self.results["run_values"]["chamber_sound_speed"].append(son_vel)
            self.results["run_values"]["chamber_rho"].append(rho_ch)
            self.results["run_values"]["mass_flow"].append(total_mass_flow)
            self.results["run_values"]["mass_flow_ox"].append(ox_flow)
            self.results["run_values"]["mass_flow_f"].append(m_fuel)
            self.results["run_values"]["chamber_speed"].append(u_ch)

            # ---------------- Perform the regression of the block and erode ----------------
            self.geometry.regress(ox_flow, dt, solution)
            self.nozzle.erode(dt)

            # ------------------------------ Update the loop --------------------------------
            k += 1

            # Update the flag for burn-time
            if max_burn_time:
                flag_burn = k * dt <= max_burn_time

        # ------------------------ POST-PROCESS ------------------------

        # Convert to numpy arrays
        self.results['run_values'] = {key: asarray(value) for key, value in self.results["run_values"].items()}

        # Post-process the data
        self._post_process_data(dt=dt)

    def _solve_coupled_analysis(self, hydraulic_module, chamber_pressure, rho_fuel, tol_press):
        """
        _solve_pressure_loop solves the internal pressure loop for the combustion chamber. Private method
        :param chamber_pressure: initial value of the chamber pressure [Pa]
        :param ox_flow: oxidizer flow [kg/s]
        :param rho_fuel: fuel density [kg/m^3]
        :param tol_press: tolerance on the pressure value
        :return: TBD
        """
        error = 1e3  # Define initial value for the error

        while tol_press < error:
            # ----------- Set the chamber pressure of the Hydraulic Module and Solve ------

            hydraulic_module.set_chamber_pressure(chamber_pressure)
            hydraulic_module.run_simulation()
            m_ox = hydraulic_module.return_exit_flow()

            # --------------------- Solve the equation of the block -----------------------------
            solution, mass_flows, of_ratio, mass_fluxes = self.geometry.solve_mass_flux(m_ox, rho_fuel)

            # --------------------- Run Thermochemical Analysis ---------------------------
            t_chamber, gamma, gas_molar_mass, cea_c_star, son_vel, rho_ch = self._run_thermochemical_analysis(
                of_ratio, chamber_pressure / 10 ** 5)

            m_ox, m_fuel = mass_flows
            g_ox, g_f = mass_fluxes
            total_mass_flow = m_ox + m_fuel

            # Compute by fixed point method the new value of the chamber pressure [Pa]
            old_pressure = chamber_pressure
            chamber_pressure = total_mass_flow * cea_c_star / self.nozzle.get_throat_area()
            error = np.abs(chamber_pressure - old_pressure) / 1e5

        # Return the results of the loop
        # noinspection PyUnboundLocalVariable
        return solution, of_ratio, m_fuel, m_ox, total_mass_flow, g_ox, g_f, chamber_pressure, t_chamber, gamma, \
               gas_molar_mass, cea_c_star, son_vel, rho_ch

    def _solve_pressure_loop(self, chamber_pressure, ox_flow, rho_fuel, tol_press):
        """
        _solve_pressure_loop solves the internal pressure loop for the combustion chamber. Private method
        :param chamber_pressure: initial value of the chamber pressure [Pa]
        :param ox_flow: oxidizer flow [kg/s]
        :param rho_fuel: fuel density [kg/m^3]
        :param tol_press: tolerance on the pressure value
        :return: solution, of_ratio, m_ox, m_fuel, total_mass_flow, g_ox, g_f, chamber_pressure, t_chamber, gamma,
        gas_molar_mass, cea_c_star, son_vel, rho_ch
        """
        error = 1e3  # Define initial value for the error

        while tol_press < error:
            # --------------------- Solve the equation of the block -----------------------------
            solution, mass_flows, of_ratio, mass_fluxes = self.geometry.solve_mass_flux(ox_flow, rho_fuel)

            # --------------------- Run Thermochemical Analysis ---------------------------
            t_chamber, gamma, gas_molar_mass, cea_c_star, son_vel, rho_ch = self._run_thermochemical_analysis(
                of_ratio, chamber_pressure / 10 ** 5)

            m_ox, m_fuel = mass_flows
            g_ox, g_f = mass_fluxes
            total_mass_flow = m_ox + m_fuel

            # Compute by fixed point method the new value of the chamber pressure [Pa]
            old_pressure = chamber_pressure
            chamber_pressure = total_mass_flow * cea_c_star / self.nozzle.get_throat_area()
            error = np.abs(chamber_pressure - old_pressure)

        # Return the results of the loop
        # noinspection PyUnboundLocalVariable
        return solution, of_ratio, m_ox, m_fuel, total_mass_flow, g_ox, g_f, chamber_pressure, t_chamber, gamma, \
               gas_molar_mass, cea_c_star, son_vel, rho_ch

    def run_balanced_nozzle_analysis(self, ox_flow, safety_thickness, dt, max_burn_time=None, tol_press=1e-3, **kwargs):
        """
        run_balanced_nozzle_analysis runs the analysis on the combustion chamber considering the balance of the choked
        nozzle.
        :param ox_flow: value of input oxidizer flow
        :param safety_thickness: minimum allowed thickness for fuel slivers, which conditions burn termination
        :param dt: time increment
        :param max_burn_time: maximum burn time (if inputted)
        :param tol_press: tolerance for the pressure loop solver
        :return: nothing
        """
        # Check the inputs
        assert ox_flow > 0, "Oxidizer flow not greater than 0, check your inputs to CombustionModule. \n"
        assert safety_thickness > 0, "Safety thickness not greater than 0, check your inputs to CombustionModule. \n"
        assert dt > 0, "Simulation time-increment not greater than 0, check your inputs to CombustionModule. \n"

        # Call for workspace definition methods
        g0, R, pression_atmo, a, n, m, rho_fuel, initial_chamber_pressure, eta_comb = self._unpack_data_dictionary()

        # -------------------- Initialize simulation_variables if required -------------
        pression_chambre = initial_chamber_pressure

        # ---------------------------- MAIN SIMULATION LOOP -----------------------------
        # Set a counter to keep-track of the loop
        k = 1

        # Set the flag for the maximum burn-time
        flag_burn = True

        while self.geometry.min_bloc_thickness() > safety_thickness and flag_burn:

            # ------------------- Perform block solution with pressure coupling  ---------------------
            solution, of_ratio, m_ox, m_fuel, total_mass_flow, g_ox, g_f, pression_chambre, t_chamber, gamma,\
                gas_molar_mass, cea_c_star, son_vel, rho_ch = self._solve_pressure_loop(pression_chambre, ox_flow,
                                                                                        rho_fuel, tol_press)
            r = R / gas_molar_mass  # Specific gaz constant

            # Determine the flow port exit-speed
            u_ch = calculate_flow_speed(cross_section=self.geometry.mesh.cells[-1].return_area_data(),
                                        mass_flow=total_mass_flow,
                                        density=rho_ch)

            # ------------------------- Determine Propulsion Performances ------------------
            mach_exit = iso.exit_mach_via_expansion(gamma=gamma, expansion=self.nozzle.get_expansion_ratio(),
                                                    supersonic=True)
            exit_pressure = pression_chambre / iso.pressure_ratio(gamma, mach_exit)
            t_exit = t_chamber / iso.temperature_ratio(gamma, mach_exit)
            v_exit = math.sqrt(gamma * r * t_exit) * mach_exit

            thrust = self.nozzle.get_nozzle_effeciency() * (
                    total_mass_flow * v_exit + (exit_pressure - pression_atmo) *
                    self.nozzle.get_exit_area())

            # Calculate the nozzle equation validation
            nozzle_p = total_mass_flow * cea_c_star / (self.nozzle.get_throat_area() * pression_chambre)
            isp = thrust / total_mass_flow / g0

            # ------------------------------ Results updates --------------------------------
            self.results['run_values']['time'].append(k * dt)
            self.results["run_values"]["thrust"].append(thrust)
            self.results["run_values"]["isp"].append(isp)
            self.results["run_values"]["pressure"].append(pression_chambre)
            self.results["run_values"]["temperature"].append(t_chamber)
            self.results["run_values"]["of"].append(of_ratio)
            self.results["run_values"]["Go"].append(g_ox)
            self.results["run_values"]["nozzle_param"].append(nozzle_p)
            self.results["run_values"]["c_star"].append(cea_c_star)
            self.results["run_values"]["chamber_sound_speed"].append(son_vel)
            self.results["run_values"]["chamber_rho"].append(rho_ch)
            self.results["run_values"]["mass_flow"].append(total_mass_flow)
            self.results["run_values"]["mass_flow_ox"].append(ox_flow)
            self.results["run_values"]["mass_flow_f"].append(m_fuel)
            self.results["run_values"]["chamber_speed"].append(u_ch)

            # ---------------- Perform the regression of the block and erode ----------------
            self.geometry.regress(ox_flow, dt, solution)
            self.nozzle.erode(dt)

            # ------------------------------ Update the loop --------------------------------
            k += 1
            # Update the flag for burn-time
            if max_burn_time:
                flag_burn = k * dt <= max_burn_time

        # ------------------------ POST-PROCESS ------------------------
        self.results['run_values'] = {key: asarray(value) for key, value in self.results["run_values"].items()}

        # Post-process the data
        self._post_process_data(dt=dt)

    def run_full_flow_analysis(self, hydraulic_module, safety_thickness, dt, max_burn_time=None, tol_press=1e-3,
                               **kwargs):
        """
        run_full_flow_analysis executes the multidisciplinary simulation of the hydraulic_module and the
        combustion_module
        :param hydraulic_module: Hydraulic Module instance.
        :param safety_thickness: minimum allowed thickness for fuel slivers, which conditions burn termination
        :param dt: time increment
        :param max_burn_time: maximum allowed burn time [s]
        :param tol_press: tolerance for the pressure loop solver
        :return: TBD
        """
        # Check the inputs
        assert isinstance(hydraulic_module, HydraulicModule), "Hydraulic Module has to be a HydraulicModule class " \
                                                              "instance.\n"
        assert safety_thickness > 0, "Safety thickness not greater than 0, check your inputs to CombustionModule. \n"
        assert dt > 0, "Simulation time-increment not greater than 0, check your inputs to CombustionModule. \n"

        # Call for workspace definition methods
        g0, R, pression_atmo, a, n, m, rho_fuel, initial_chamber_pressure, eta_comb = self._unpack_data_dictionary()

        # -------------------- Initialize simulation_variables if required -------------
        pression_chambre = initial_chamber_pressure

        # ---------------------------- MAIN SIMULATION LOOP -----------------------------

        # Set a counter to keep-track of the loop
        k = 1

        # Set the flag for the maximum burn-time
        flag_burn = True

        while self.geometry.min_bloc_thickness() > safety_thickness and flag_burn:

            # Execute block inside try-except to catch exception
            try:

                # ------------------- Perform block solution with pressure coupling  ---------------------
                solution, of_ratio, m_fuel, m_ox, total_mass_flow, g_ox, g_f, pression_chambre, t_chamber, gamma, \
                    gas_molar_mass, cea_c_star, son_vel, rho_ch = self._solve_coupled_analysis(hydraulic_module,
                                                                                               pression_chambre,
                                                                                               rho_fuel,
                                                                                               tol_press)
                r = R / gas_molar_mass  # Specific gaz constant

                # Determine the flow port-speed
                u_ch = calculate_flow_speed(cross_section=self.geometry.total_cross_section_area(),
                                            mass_flow=total_mass_flow,
                                            density=rho_ch)

                # Calculate chamber conditions and motor performance
                mach_exit = iso.exit_mach_via_expansion(gamma=gamma, expansion=self.nozzle.get_expansion_ratio(),
                                                        supersonic=True)
                exit_pressure = pression_chambre / iso.pressure_ratio(gamma, mach_exit)
                t_exit = t_chamber / iso.temperature_ratio(gamma, mach_exit)
                v_exit = math.sqrt(gamma * r * t_exit) * mach_exit

                thrust = self.nozzle.get_nozzle_effeciency() * (
                        total_mass_flow * v_exit + (exit_pressure - pression_atmo) *
                        self.nozzle.get_exit_area())

                # Calculate the nozzle equation validation
                nozzle_p = total_mass_flow * cea_c_star / (self.nozzle.get_throat_area() * pression_chambre)
                isp = thrust / total_mass_flow / g0

                # Results updates
                self.results['run_values']['time'].append(k * dt)
                self.results["run_values"]["thrust"].append(thrust)
                self.results["run_values"]["isp"].append(isp)
                self.results["run_values"]["pressure"].append(pression_chambre)
                self.results["run_values"]["temperature"].append(t_chamber)
                self.results["run_values"]["of"].append(of_ratio)
                self.results["run_values"]["Go"].append(g_ox)
                self.results["run_values"]["nozzle_param"].append(nozzle_p)
                self.results["run_values"]["c_star"].append(cea_c_star)
                self.results["run_values"]["chamber_sound_speed"].append(son_vel)
                self.results["run_values"]["chamber_rho"].append(rho_ch)
                self.results["run_values"]["mass_flow"].append(total_mass_flow)
                self.results["run_values"]["mass_flow_ox"].append(m_ox)
                self.results["run_values"]["mass_flow_f"].append(m_fuel)
                self.results["run_values"]["chamber_speed"].append(u_ch)

                # ---------------- Perform the regression of the block and erode ----------------
                self.geometry.regress(m_ox, dt, solution)
                self.nozzle.erode(dt)

                # ------------------------------ Update the loop --------------------------------
                k += 1                                      # Update the counter
                hydraulic_module.update(dt)                 # Update the hydraulic module

                # Update the flag for burn-time
                if max_burn_time:
                    flag_burn = k * dt <= max_burn_time

                # Print the solution to the screen
                self._print_last_solution_step()

            # End the burn by setting the flag_burn to False
            except HydraulicModule.RunEndingError:
                flag_burn = False

        # ------------------------ POST-PROCESS ------------------------

        # Convert to numpy arrays
        self.results['run_values'] = {key: asarray(value) for key, value in self.results["run_values"].items()}

        # Post-process the data
        self._post_process_data(dt=dt)


# ------------------------- 0D COMBUSTION DEFINITIONS ----------------------

class CombustionObject0D(CombustionObject, ABC):
    """
    The CombustionObject0D inherits from the CombustionObject class and is in charge
    of working out the simulation at the combustion chamber level in the case of a 0D
    model implementation. Along the attributes of the parent class its other attributes
    are:

    # Attributes:
        1. geometry: Geometry instance used to define the geometry of the port (hybrid rocket)
    """

    def __init__(self, json_interpreter, geometry_object, nozzle_object):
        """
        class initializer
        :param json_interpreter: JsonInterpreter instance used to collect the data
        :param geometry_object: Geometry0D instance
        :param nozzle_object: Nozzle instance
        """
        # Assert nature of the inputs
        assert isinstance(geometry_object, Geometry0D), " Please insert a valid Geometry0D instance. \n"

        # Call superclass initializer
        super(CombustionObject0D, self).__init__(json_interpreter, geometry_object, nozzle_object)


class CombustionObjectClassic(CombustionObject0D):
    """
    CombustionObjectClassic is the CombustionObjectOD implementation of the combustion object
    with regular behavior.
    """

    def __init__(self, json_interpreter, geometry_object, nozzle_object):
        super(CombustionObjectClassic, self).__init__(json_interpreter, geometry_object, nozzle_object)


class CombustionObjectImage(CombustionObject0D):
    """
    CombustionObjectImage implements the behavior of the Geometric profile based
    on an image definition.
    """

    def __init__(self, json_interpreter, geometry_object, nozzle_object):
        """ class initializer """
        # Assert the nature of the input
        assert isinstance(geometry_object, SinglePortImageGeometry), "Geometry is not image-based, please " \
                                                                     "use SinglePortGeometry class. \n"
        # Call superclass initializer
        super(CombustionObjectImage, self).__init__(json_interpreter, geometry_object, nozzle_object)

    def _post_process_data(self, dt):
        """ exclusive post_process data for CombustionObjectImage
        class. Enhance polymorphism
        """
        # Re-map the run results
        self._remap_run_results(dt)

        # Call superclass _post_process_data
        super(CombustionObjectImage, self)._post_process_data(dt)

    def _remap_run_results(self, dt):
        """
        _remap_run_results has for function to remap the results obtained by the Image Geometry
        code.
        :param dt: time-step
        :return: nothing
        """
        # Loop through the variables on run_results and remap by doing interpolation of the results
        # to increase the density of the array
        new_time_array = np.arange(0, self.results['run_values']['time'][-1], dt)
        time_interp = self.results['run_values']['time']

        for key, value in self.results['run_values'].items():
            if key != 'time':
                self.results['run_values'][key] = np.interp(new_time_array, time_interp, value)

        # Set new time array
        self.results['run_values']['time'] = new_time_array

    def run_simulation_constant_fuel_sliver(self, ox_flow, safety_thickness, dt, max_burn_time=None, **kwargs):
        # Check the inputs
        assert ox_flow > 0, "Oxidizer flow not greater than 0, check your inputs to CombustionModule. \n"
        assert safety_thickness > 0, "Safety thickness not greater than 0, check your inputs to CombustionModule. \n"

        # Call for workspace definition methods
        g0, R, pression_atmo, a, n, m, rho_fuel, initial_chamber_pressure, eta_comb = self._unpack_data_dictionary()

        # -------------------- Initialize simulation_variables if required -------------
        pression_chambre = initial_chamber_pressure

        # ---------------------------- MAIN SIMULATION LOOP -----------------------------
        # Set the flag for the maximum burn-time
        flag_burn = True
        time = 0

        # Get the thickness equivalent to 5 pixels
        dr = self.geometry.get_meters_per_pixel() * 5

        while self.geometry.min_bloc_thickness() > safety_thickness and flag_burn:

            # Regression rate and mass flow calculations
            m_fuel = self.geometry.compute_fuel_rate(rho=rho_fuel, ox_flow=ox_flow)
            total_mass_flow = ox_flow + m_fuel
            of_ratio = ox_flow / m_fuel
            g_ox = ox_flow / self.geometry.total_cross_section_area()
            delta_t = dr / self.geometry.regression_model.compute_regression_rate(self.geometry, ox_flow)

            # Call CEA process to obtain thermochemical variables
            t_chamber, gamma, gas_molar_mass, cea_c_star, son_vel, rho_ch = self._run_thermochemical_analysis(
                of_ratio, pression_chambre / 10 ** 5)
            r = R / gas_molar_mass  # Specific gaz constant

            # Determine the flow port-speed
            u_ch = calculate_flow_speed(cross_section=self.geometry.total_cross_section_area(),
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
            nozzle_p = total_mass_flow * cea_c_star / (self.nozzle.get_throat_area() * pression_chambre)
            isp = thrust / total_mass_flow / g0

            # Results updates
            self.results['run_values']['time'].append(time)
            self.results["run_values"]["thrust"].append(thrust)
            self.results["run_values"]["isp"].append(isp)
            self.results["run_values"]["pressure"].append(pression_chambre)
            self.results["run_values"]["temperature"].append(t_chamber)
            self.results["run_values"]["of"].append(of_ratio)
            self.results["run_values"]["Go"].append(g_ox)
            self.results["run_values"]["nozzle_param"].append(nozzle_p)
            self.results["run_values"]["c_star"].append(cea_c_star)
            self.results["run_values"]["chamber_sound_speed"].append(son_vel)
            self.results["run_values"]["chamber_rho"].append(rho_ch)
            self.results["run_values"]["mass_flow"].append(total_mass_flow)
            self.results["run_values"]["mass_flow_ox"].append(ox_flow)
            self.results["run_values"]["mass_flow_f"].append(m_fuel)
            self.results["run_values"]["chamber_speed"].append(u_ch)

            # ---------------- Perform the regression of the block and erode ----------------
            self.geometry.regress(ox_flow=ox_flow, dt=delta_t)
            self.nozzle.erode(delta_t)

            # ------------------------------ Update the loop --------------------------------
            time += delta_t

            # Update the flag for burn-time
            if max_burn_time:
                flag_burn = time <= max_burn_time

        # ------------------------ POST-PROCESS ------------------------

        # Convert to numpy arrays
        self.results['run_values'] = {key: asarray(value) for key, value in self.results["run_values"].items()}

        # Post-process the data
        self._post_process_data(dt=dt)

    def run_balanced_nozzle_analysis(self, ox_flow, safety_thickness, dt, max_burn_time=None, tol_press=1e-3, **kwargs):
        # Check the inputs
        assert ox_flow > 0, "Oxidizer flow not greater than 0, check your inputs to CombustionModule. \n"
        assert safety_thickness > 0, "Safety thickness not greater than 0, check your inputs to CombustionModule. \n"

        # Call for workspace definition methods
        g0, R, pression_atmo, a, n, m, rho_fuel, initial_chamber_pressure, eta_comb = self._unpack_data_dictionary()

        # -------------------- Initialize simulation_variables if required -------------
        pression_chambre = initial_chamber_pressure

        # ---------------------------- MAIN SIMULATION LOOP -----------------------------

        # Set the flag for the maximum burn-time
        flag_burn = True
        time = 0

        # Get the thickness equivalent to 5 pixels
        dr = self.geometry.get_meters_per_pixel() * 5

        while self.geometry.min_bloc_thickness() > safety_thickness and flag_burn:

            # ------------------- Perform block solution with pressure coupling  ---------------------

            solution, of_ratio, m_ox, m_fuel, total_mass_flow, g_ox, g_f, pression_chambre, t_chamber, gamma,\
                gas_molar_mass, cea_c_star, son_vel, rho_ch = self._solve_pressure_loop(pression_chambre, ox_flow,
                                                                                        rho_fuel, tol_press)

            r = R / gas_molar_mass  # Specific gaz constant
            delta_t = dr / self.geometry.regression_model.compute_regression_rate(self.geometry, ox_flow)

            # Determine the flow port-speed
            u_ch = calculate_flow_speed(cross_section=self.geometry.total_cross_section_area(),
                                        mass_flow=total_mass_flow,
                                        density=rho_ch)

            # Calculate chamber conditions and motor performance
            mach_exit = iso.exit_mach_via_expansion(gamma=gamma, expansion=self.nozzle.get_expansion_ratio(),
                                                    supersonic=True)
            exit_pressure = pression_chambre / iso.pressure_ratio(gamma, mach_exit)
            t_exit = t_chamber / iso.temperature_ratio(gamma, mach_exit)
            v_exit = math.sqrt(gamma * r * t_exit) * mach_exit

            thrust = self.nozzle.get_nozzle_effeciency() * (
                    total_mass_flow * v_exit + (exit_pressure - pression_atmo) * self.nozzle.get_exit_area())

            # Calculate the nozzle equation validation
            nozzle_p = total_mass_flow * cea_c_star / (self.nozzle.get_throat_area() * pression_chambre)
            # print(nozzle_p)

            isp = thrust / total_mass_flow / g0

            # Results updates
            self.results['run_values']['time'].append(time)
            self.results["run_values"]["thrust"].append(thrust)
            self.results["run_values"]["isp"].append(isp)
            self.results["run_values"]["pressure"].append(pression_chambre)
            self.results["run_values"]["temperature"].append(t_chamber)
            self.results["run_values"]["of"].append(of_ratio)
            self.results["run_values"]["Go"].append(g_ox)
            self.results["run_values"]["nozzle_param"].append(nozzle_p)
            self.results["run_values"]["c_star"].append(cea_c_star)
            self.results["run_values"]["chamber_sound_speed"].append(son_vel)
            self.results["run_values"]["chamber_rho"].append(rho_ch)
            self.results["run_values"]["mass_flow"].append(total_mass_flow)
            self.results["run_values"]["mass_flow_ox"].append(ox_flow)
            self.results["run_values"]["mass_flow_f"].append(m_fuel)
            self.results["run_values"]["chamber_speed"].append(u_ch)

            # ---------------- Perform the regression of the block and erode ----------------
            self.geometry.regress(ox_flow, delta_t, solution)
            self.nozzle.erode(delta_t)

            # ------------------------------ Update the loop --------------------------------
            time += delta_t

            # Update the flag for burn-time
            if max_burn_time:
                flag_burn = time <= max_burn_time

        # ------------------------ POST-PROCESS ------------------------

        # Convert to numpy arrays
        self.results['run_values'] = {key: asarray(value) for key, value in self.results["run_values"].items()}

        # Post-process the data
        self._post_process_data(dt=dt)

    def run_full_flow_analysis(self, hydraulic_module, safety_thickness, dt, max_burn_time=None, tol_press=1e-3,
                               **kwargs):
        # Check the inputs
        assert isinstance(hydraulic_module, HydraulicModule), "Hydraulic Module has to be a HydraulicModule class " \
                                                              "instance.\n"
        assert safety_thickness > 0, "Safety thickness not greater than 0, check your inputs to CombustionModule. \n"
        assert dt > 0, "Simulation time-increment not greater than 0, check your inputs to CombustionModule. \n"

        # Call for workspace definition methods
        g0, R, pression_atmo, a, n, m, rho_fuel, initial_chamber_pressure, eta_comb = self._unpack_data_dictionary()

        # -------------------- Initialize simulation_variables if required -------------
        pression_chambre = initial_chamber_pressure

        # ---------------------------- MAIN SIMULATION LOOP -----------------------------

        # Set the flag for the maximum burn-time
        flag_burn = True
        time = 0

        # Get the thickness equivalent to 5 pixels
        dr = self.geometry.get_meters_per_pixel() * 5

        while self.geometry.min_bloc_thickness() > safety_thickness and flag_burn:

            # Execute block inside try-except to catch exception
            try:

                # ------------------- Perform block solution with pressure coupling  ---------------------
                solution, of_ratio, m_fuel, m_ox, total_mass_flow, g_ox, g_f, pression_chambre, t_chamber, gamma, \
                gas_molar_mass, cea_c_star, son_vel, rho_ch = self._solve_coupled_analysis(hydraulic_module,
                                                                                           pression_chambre,
                                                                                           rho_fuel,
                                                                                           tol_press)

                r = R / gas_molar_mass  # Specific gaz constant
                delta_t = dr / self.geometry.regression_model.compute_regression_rate(self.geometry, m_ox)

                # Determine the flow port-speed
                u_ch = calculate_flow_speed(cross_section=self.geometry.total_cross_section_area(),
                                            mass_flow=total_mass_flow,
                                            density=rho_ch)

                # Calculate chamber conditions and motor performance
                mach_exit = iso.exit_mach_via_expansion(gamma=gamma, expansion=self.nozzle.get_expansion_ratio(),
                                                        supersonic=True)
                exit_pressure = pression_chambre / iso.pressure_ratio(gamma, mach_exit)
                t_exit = t_chamber / iso.temperature_ratio(gamma, mach_exit)
                v_exit = math.sqrt(gamma * r * t_exit) * mach_exit

                thrust = self.nozzle.get_nozzle_effeciency() * (
                        total_mass_flow * v_exit + (exit_pressure - pression_atmo) *
                        self.nozzle.get_exit_area())

                # Calculate the nozzle equation validation
                nozzle_p = total_mass_flow * cea_c_star / (self.nozzle.get_throat_area() * pression_chambre)

                isp = thrust / total_mass_flow / g0

                # Results updates
                self.results['run_values']['time'].append(time)
                self.results["run_values"]["thrust"].append(thrust)
                self.results["run_values"]["isp"].append(isp)
                self.results["run_values"]["pressure"].append(pression_chambre)
                self.results["run_values"]["temperature"].append(t_chamber)
                self.results["run_values"]["of"].append(of_ratio)
                self.results["run_values"]["Go"].append(g_ox)
                self.results["run_values"]["nozzle_param"].append(nozzle_p)
                self.results["run_values"]["c_star"].append(cea_c_star)
                self.results["run_values"]["chamber_sound_speed"].append(son_vel)
                self.results["run_values"]["chamber_rho"].append(rho_ch)
                self.results["run_values"]["mass_flow"].append(total_mass_flow)
                self.results["run_values"]["mass_flow_ox"].append(m_ox)
                self.results["run_values"]["mass_flow_f"].append(m_fuel)
                self.results["run_values"]["chamber_speed"].append(u_ch)

                # ---------------- Perform the regression of the block and erode ----------------
                self.geometry.regress(m_ox, delta_t, solution)
                self.nozzle.erode(delta_t)

                # ------------------------------ Update the loop --------------------------------
                time += delta_t                             # Calculate the time
                hydraulic_module.update(delta_t)            # Update the hydraulic module

                # Update the flag for burn-time
                if max_burn_time:
                    flag_burn = time <= max_burn_time

                # Print the solution to the screen
                self._print_last_solution_step()

            # End the burn by setting the flag_burn to False
            except HydraulicModule.RunEndingError:
                flag_burn = False

        # ------------------------ POST-PROCESS ------------------------
        # Convert to numpy arrays
        self.results['run_values'] = {key: asarray(value) for key, value in self.results["run_values"].items()}

        # Post-process the data
        self._post_process_data(dt=dt)


# ----------------------------- COMBUSTION OBJECT 1D -----------------------------

class CombustionObject1D(CombustionObject):
    """
    CombustionObject1D is a 1 dimensional implementation of the Combustion Model considering the
    empirical model described in Humble's book. It is implemented as an abstract class

        Attributes:
            1. geometry = Geometry1D class instance

    """

    def __init__(self, json_interpreter, nozzle_object, geometry_object):
        """
        Class initializer
        :param json_interpreter: JsonInterpreter instance used to collect the data
        :param geometry_object: Geometry instance
        :param nozzle_object: Nozzle instance
        """
        # Assert nature of the inputs
        assert isinstance(geometry_object, Geometry1D), " Please insert a valid Geometry instance. \n"

        # Call superclass constructor
        super(CombustionObject1D, self).__init__(json_interpreter, geometry_object, nozzle_object)