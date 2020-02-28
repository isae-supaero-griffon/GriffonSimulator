# Rossiter.py aims to perform the estimation of the Rossiter mode frequency over
# time for the tests performed at ONERA during the month of April 2019.
# @author: Jose Felix Zapata Usandivaras
# Date: 11/05/2019
# ISAE - SUPAERO, Scube GRIFFON, Combustion Team
#

# ---------------------- IMPORT MODULES -------------------------

import numpy as np                                          # Import numpy
import pandas as pd                                         # Import pandas to handle time-series
import matplotlib.pyplot as plt                             # Import matplotlib


# ------------------- FUNCTIONS DEFINITIONS ---------------------

def calculate_flow_speed(port_radius, mass_flow, density):
    """
    calculate_flow_speed determines the average speed of the port flow at the
    exit.
    :param port_radius: float with radius of the port [m]
    :param mass_flow: mass flow of mixture [kg/s]
    :param density: density of the flow [kg/m^3]
    :return: flow speed in m/s
    """
    return mass_flow / (density * np.pi * port_radius**2)


def calculate_rossiter_frequency(u, m_ross, alpha, L, c, k_ross, beta=1):
    """
    calculate_rossiter_frequency determines value of the rossiter frequency
    :param u: flow speed [m/s]
    :param m_ross: integer number associated to the vortex multiplicity
    :param alpha: calibration factor for the rossiter mode
    :param L: length of post-combustion chamber
    :param c: speed of sound in the medium
    :param k_ross: proportionality constant between flow speed and
    :param beta: u speed proportionality constant
    :return: nothing
    """

    # Calculate the effective u
    u_eff = beta * u

    # Calculate the local mach number
    mach = u_eff / c

    # Return the output
    return (u_eff / L) * (m_ross - alpha)/(mach + 1/k_ross)


def calculate_k_ross(m_ross, alpha, L, f0, u0, c0, beta=1):
    """
    calculate_k_ross performs the inverse calculation of the k_ross
    based on initially provided data
    :param m_ross: integer number of rossiter mode associated to multiplicity
    :param alpha: calibration parameter for the rossiter mode
    :param L: rossiter characteristic mode
    :param f0: approximate frequency to match [Hz]
    :param u0: initial flow speed (at time when mode appears)
    :param c0: initial flow speed of sound (obtained from CombustionModule when mode appears)
    :param beta: u speed proportionality constant
    :return: estimated k_ross
    """

    # Calculate the effective u0
    u0_eff = beta * u0

    # Calculate the local Mach Number
    mach = u0_eff / c0

    # Return the output
    return 1 / (u0_eff * (m_ross - alpha) / (L * f0) - mach)


def return_variables_time(data, variables, time):
    """
    return_variables_time returns the set of variables contained in the given dictionary
    at the specified time
    :param data: data dictionary, must contain the time key
    :param variables: tupple of variables from which we want to fetch the dictionary data
    :param time: time of simulation [secs] at which we desire to obtain the data
    :return: dictionary with requested variables values
    """

    assert 'time' in data, "Time key not present in input data dictionary, review dictionary. \n"

    # Extract from the data the required fields
    count = 0
    flag = True

    # Loop until we find the first time value that matches the criteria
    while flag:
        t = data['time'][count]
        if time <= t:
            flag = False
        else:
            count += 1

    # Extract the data
    if not flag:
        # Return the values if a time was found, otherwise return nan
        sub_dict = {key : data[key][count] for key in variables}
    else:
        sub_dict = {key : np.nan for key in variables}

    # Return the sub_dict
    return sub_dict


def calculate_vortex_speed(u, k_ross, beta=1):
    """
    calculate_vortex_speed determines the average vortex speed associated to the rossiter mode
    :param u: chamber speed  [m/s]
    :param k_ross: rossiter constant that relates free-flow speed and vortex speed
    :param beta: u speed proportionality constant
    :return: vortex speed [m/s]
    """
    # Return the output
    return beta * k_ross * u


# --------------------- CLASS DEFINITIONS -----------------------


class Rossiter:
    """
    The Rossiter class intends to perform the calculations required to
    estimate the Rossiter mode and its drift.

    # Attributes:
        1. modes: RossiterMode array with modes of interest
        2. combustion: CombustionModule object already run
        3. test_results: data associated to the tests
        4. test_name: string containing the name of the test to ease identification
    """

    def __init__(self, post_combustion, modes_params, combustion_module, test_file_name, test_name, time_shift):
        """
        class initializer
        :param post_combustion: PostCombustion class instance with post combustion geometric data
        :param modes_params: dict containing the parameters associated to the modes
        :param combustion_module: CombustionModule instance
        :param test_file_name: file name of txt file which contains the data of the test Rossiter modes.
        :param test_name: name of the test
        :param time_shift: time shift between test data and simulation
        """

        # Set the properties
        self.modes = self._initialize_modes(modes_params, post_combustion)
        self.combustion = combustion_module
        self.test_results = self._load_test_data(test_file_name, time_shift)
        self.test_name = test_name

    def calibrate_modes(self):
        """
        calibrate_modes sets the calibration procedure for the different modes
        :return: nothing
        """

        # Use as reference the first mode that appears in the list
        self.modes[0].calibrate_mode(self.combustion.results['run_values'])

        # Loop through the mode and set the k_ross
        for mode in self.modes[1:]:
            mode.k_ross = self.modes[0].k_ross

    def calculate_modes_frequencies(self):
        """
        calculate_modes_frequencies starts the calculations for the modes involved
        :return: nothing
        """

        # Loop through every mode and execute the corresponding method
        for mode in self.modes:
            mode.calculate_frequency_evolution(self.combustion.results['run_values'])

    @staticmethod
    def _initialize_modes(modes_params, post_combustion):
        """
        _initialize_modes creates the RossiterMode objects
        :param modes_params:
        :return:
        """

        # Initialize a list of modes
        modes = []

        # Loop through the input modes and initialize them
        for value in modes_params:
            value['post_combustion'] = post_combustion      # Assign post_combustion to the dictionary
            modes.append(RossiterMode(**value))

        # Return the array of modes
        return modes

    @staticmethod
    def _load_test_data(test_file_name, time_shift):
        """
        _load_test_data loads the data coming from the test exported previously to a txt file
        :param test_file_name: txt file name
        :param time_shift: time shift since the beginning of the test
        :return: dictionary containing the data
        """

        # Use numpy loadtxt function
        times, frequencies = np.loadtxt(test_file_name, skiprows=1, delimiter='\t', unpack=True)

        # Return the dictionary with the test results
        return {'t': times - time_shift, 'f': frequencies}

    def plot_results(self):
        """
        plot_results plots the obtained results from the calculation of the Rossiter modes
        :return: nothing
        """

        # Set the font dictionaries (for plot title and axis titles)
        title_font = {'size': '20', 'color': 'black', 'weight': 'normal',
                      'verticalalignment': 'bottom'}  # Bottom vertical alignment for more space
        axis_font = {'size': '16'}

        # Generate the plots
        fig = plt.figure(facecolor='w', figsize=(20, 20))
        fig.suptitle('Rossiter Modes Analysis', **title_font)

        axs = [plt.subplot2grid((2, 2), (0, 0), rowspan=1, colspan=1),
               plt.subplot2grid((2, 2), (1, 0), rowspan=1, colspan=1),
               plt.subplot2grid((2, 2), (0, 1), rowspan=1, colspan=1),
               plt.subplot2grid((2, 2), (1, 1), rowspan=1, colspan=1)]

        # Set the label properties
        for ax in axs:
            for label in (ax.get_xticklabels() + ax.get_yticklabels()):
                label.set_fontname('Arial')
                label.set_fontsize(10)

        # --------------- ROSSITER PLOT -----------------

        # Plot the data associated to the test
        axs[0].plot(self.test_results['t'], self.test_results['f'], color='lightcoral',
                  label="%s Test" % self.test_name, linestyle='',
                  marker='*', markersize=1.0, markerfacecolor='lightcoral')

        # Plot the modes results
        for mode in self.modes:
            axs[0].plot(mode.results['t'], mode.results['f'], label='m = %d'%mode.m_ross,
                        linestyle='--', linewidth=1.5)

        # Set remaining axes properties
        axs[0].set_title('Rossiter Modes frequencies')
        axs[0].set_ylabel('Frequency [Hz]', **axis_font)
        axs[0].set_xlabel('Time [s]', **axis_font)
        axs[0].grid(b=True, axis='both')

        # Show the legends
        axs[0].legend(loc=0, fancybox=False, frameon=False, fontsize=12)

        # --------------- VORTEX SPEED PLOT -------------

        # Extract the data
        x_data, y_data = self.modes[0].results['t'], self.modes[0].results['u_vor']


        # Plot the modes results
        axs[1].plot(x_data, y_data, label='u_vor',
                    linestyle='-', linewidth=1.5, color='red')

        # Set remaining axes properties
        axs[1].set_title('Vortex Speeds')
        axs[1].set_ylabel('u_vor [m/s]', **axis_font)
        axs[1].set_xlabel('Time [s]', **axis_font)
        axs[1].grid(b=True, axis='both')

        # Show the legends
        axs[1].legend(loc=0, fancybox=False, frameon=False, fontsize=12)

        # ------------------ FLOW SPEED ------------------

        ax2 = axs[1].twinx()

        # Extract the data
        x_data, y_data = self.combustion.results['run_values']['time'],\
                         self.combustion.results['run_values']['chamber_speed']

        # Plot the data
        ax2.tick_params(axis='y')
        ax2.plot(x_data, y_data,  label='u',
                    linestyle='-', linewidth=1.5, color='blue')

        # Set remaining axes properties
        ax2.set_ylabel('u [m/s]', **axis_font)

        # Show the legends
        ax2.legend(loc=0, fancybox=False, frameon=False, fontsize=12)

        # --------------- DIAMETER GROWTH ----------------

        # Extract the data
        x_data, y_data = self.combustion.results['run_values']['time'],\
                         self.combustion.results['run_values']['hydraulic_port_diameter']
        axs[2].plot(x_data, 1000 * y_data, label='Hydraulic Diameter', linestyle='-', linewidth=1.5, color='blue')

        # Set remaining axes properties
        axs[2].set_title('Port Hydraulic Diameter')
        axs[2].set_ylabel('Diameter [mm]', **axis_font)
        axs[2].set_xlabel('Time [s]', **axis_font)
        axs[2].grid(b=True, axis='both')

        # Show the legends
        axs[2].legend(loc=0, fancybox=False, frameon=False, fontsize=12)

        # ------------- O/F RATIO PLOT ---------------------

        # Extract the data
        y_data = self.combustion.results['run_values']['of']
        axs[3].plot(x_data, y_data, label='O/F Ratio', linestyle='-', linewidth=1.5, color='green')

        # Set the remaining axes properties
        axs[3].set_title('O/F Ratio')
        axs[3].set_ylabel('O/F', **axis_font)
        axs[3].set_xlabel('Time [s]', **axis_font)
        axs[3].set_ylim(bottom=0)
        axs[3].grid(b=True, axis='both')

        # Show the legends
        axs[3].legend(loc=0, fancybox=False, frameon=False, fontsize=12)


class RossiterMode:
    """
    RossiterMode is a class that is to be used to calculate the frequencies
    associated to a particular mode.

    # Attributes:
        0. post_combustion: PostCombustion class instance with post combustion
        geometric data
        1. m_ross: integer indicating the mode's vortex multiplicity
        2. k_ross: multiplier of velocity to obtain the mode's vortex
        3. alpha: float rossiter fitting coefficient
        characteristic speed.
        4. results: dict containing the results of the frequency associated to
        the mode.
    """

    def __init__(self, m_ross, alpha, start_time, f0, post_combustion):
        """
        class initializer
        :param m_ross: integer indicating the mode's vortex multiplicity
        :param alpha: float for rossiter fitting coefficient
        :param start_time: time at which the mode appears in the simulation
        :param f0: initial frequency [Hz] at which the mode appears
        :param post_combustion: PostCombustion class instance with post combustion
        geometric data.
        """

        # Set the properties values
        self.m_ross = m_ross
        self.k_ross = np.nan                                # Set k_ross to nan since later it will be calibrated
        self.alpha = alpha
        self.t0 = start_time
        self.f0 = f0
        self.post_combustion = post_combustion
        self.results = {'t': [], 'f': [], 'u_vor': []}       # Initialize an empty dictionary for the results

    def calibrate_mode(self, combustion_data):
        """
        calibrate_mode determines the value of the k_ross quantity to match
        :param combustion_data: dictionary containing the data obtained from combustion simulation
        :return: nothing
        """

        # Extract parameters of relevance
        variables = ('chamber_sound_speed', 'hydraulic_port_diameter', 'chamber_speed')
        sub_dict = return_variables_time(combustion_data, variables, self.t0)

        # Determine the beta coefficient
        beta = self.post_combustion.calculate_beta(sub_dict['hydraulic_port_diameter'] / 2)

        # Calculate the k_ross
        self.k_ross = calculate_k_ross(m_ross=self.m_ross, alpha=self.alpha,
                                       L=self.post_combustion.L, f0=self.f0, u0=sub_dict['chamber_speed'],
                                       c0=sub_dict['chamber_sound_speed'], beta=beta)

    def calculate_frequency_evolution(self, combustion_results):
        """
        calculate_frequency_evolution estimates the frequency evolution over time
        of the burn.
        :param combustion_results:
        :return: nothing
        """

        # Extract combustion_results that are of interest
        vars_ = ('time', 'hydraulic_port_diameter', 'chamber_sound_speed', 'chamber_speed')
        sub_dict = {k : combustion_results.get(k, None) for k in vars_}

        # Reshape the dictionary
        combustion_data = pd.DataFrame(sub_dict).to_dict('records')

        # Loop through the records
        for record in combustion_data:
            t = record['time']

            # Check for the initial time in which the mode appears
            if self.t0 <= t:
                # Extract the required data
                r, c, u = record['hydraulic_port_diameter']/2, record['chamber_sound_speed'], record['chamber_speed']

                # Calculate the beta-coefficient
                beta = self.post_combustion.calculate_beta(r_port=r)

                # Determine the frequency
                f = calculate_rossiter_frequency(u, self.m_ross, self.alpha, self.post_combustion.L,
                                                 c, self.k_ross, beta)

                # Determine the vortex speed
                u_vor = calculate_vortex_speed(u, self.k_ross, beta)

                # Append result to dictionary
                self.results['t'].append(t)
                self.results['f'].append(f)
                self.results['u_vor'].append(u_vor)

        # Reconvert the results to numpy arrays
        self.results = {key: np.asarray(value) for key, value in self.results.items()}


class PostCombustion:
    """
    PostCombustion is a class that is to be used to help in the calculation of
    the hydrodinamic instabilities.

    # Attributes:
        1. length: float length of post-combustion chamber
        2. D: float diameter of post-combustion chamber
    """

    def __init__(self, L, D):
        """
        class initializer
        :param L: float length of post-combustion chamber
        :param D: float diameter of post-combustion chamber
        """

        # Set attributes
        self.L = L
        self.D = D

    def calculate_beta(self, r_port):
        """
        calculate the beta expansion coefficient
        :param r_port: radius of the port
        :return: beta coefficient
        """
        return (2 * r_port / self.D) ** 2
