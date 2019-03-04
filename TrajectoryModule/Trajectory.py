# -*- coding: utf-8 -*-
"""
Rocket parabolic flight calculator, from Thrust Curve - based on Eric Brunner code.
Author: Iñigo Reizabal Arregui /  Jose Felix Zapata Usandivaras
Date: 14 Dec 2018
"""

# ------------------- IMPORT MODULES ----------------------

from TrajectoryModule.Drag import *                         # Import the drag module
from TrajectoryModule.Density import *                      # Import the density module
import numpy as np                                          # Import numpy
import matplotlib.pyplot as plt                             # Import matplotlib for plot generation


# ----------------- FUNCTION DEFINITIONS ------------------


def mass_calc(thrust_vector, isp, g0, m0):
    """
    Calculation of the mass in the point where the thrust curve is defined
    for a constant Isp
    Necessary variables are
        -Thrust vector, defined as 2 line n column matrix where
         first line is the thrust in a discrete time and second line is
         the time
        -m0 is the mass of the rocket at t=0
        -g0=gravity (9,81m/s^2)
        -Isp is the specific impulse vector
    It returns a vector with 3 lines and n columns where
        -First line is the mass at each time
        -Second line is the time
        -Third line is the dm(kg/s) at each time

    """

    # Initialize the values
    dm = 0
    pos = 0
    mass_vector = np.zeros((3, thrust_vector.shape[1]))
    mass_vector[0, pos] = m0
    mass_vector[1, pos] = thrust_vector[1, pos]

    # Run the loop
    for _ in thrust_vector[1, 1:]:
        pos = pos + 1
        mass_vector[0, pos] = mass_vector[0, pos - 1] - dm
        mass_vector[1, pos] = thrust_vector[1, pos]
        mass_vector[2, pos - 1] = dm
        b = thrust_vector[0, pos]
        t1 = thrust_vector[1, pos]
        t0 = thrust_vector[1, pos - 1]
        isp1 = isp[pos]
        isp0 = isp[pos - 1]
        if not np.isnan(isp1):
            # TODO: Check the ISP implementation
            dm = b * t1 / (isp1 * g0) - b * t0 / (isp0 * g0)
        else:
            dm = 0

    mass_vector[0, pos] = mass_vector[0, pos - 1]
    mass_vector[1, pos] = thrust_vector[1, pos]
    mass_vector[2, pos] = dm

    # Return output
    return mass_vector


def interpol(vector, line1, line2, pos, x):
    """
    Linear interpolation for a vector where:
        -line1 are the points where "y" is defined
        -line2 are the points where "x" is defined
        -pos is the position the n value in the vector
        -x is the variable that has to be interpolated to obtain y
        y(x)=y(x1)+(y(x2)-y(x1))/(x2-x1)*(x-x1)
    """
    m = (vector[line1, pos + 1] - vector[line1, pos]) / (vector[line2, pos + 1] - vector[line2, pos])
    y = vector[line1, pos] + m * (x - vector[line2, pos])

    # Return output
    return y


def newton(thrust_vector, mass_vector, drag, pos, t, h, v, dv):
    """
    Calculation of the acceleration at discrete time=t
        -Thrust vector, defined as 2 line n column matrix where
         first line is the thrust(N) in a discrete time and second line is
         the time(s)
        -Mass vector, defined as 3 line n column matrix where
         first line is the mass(kg) at each time, second line is the time(s)
         and third line is the dm(kg/s) at each time
        -drag: drag object used to compute the drag exerted by the fluid
        -pos is the position of the time in the vector
        -t(s) is the time for which we calculate the acc.
        -h is the altitude
        -v is the velocity
        -dv is the differential of the velocity of an intermediate step

    Hypothesis:
        -for each step, rho(n+1)=rho(n)
    """

    thrust = interpol(thrust_vector, 0, 1, pos, t)

    # Set T to 0 if it yields a value less than 0
    if thrust < 0:
        thrust = 0

    # Perform newton's method
    drag_force = drag.compute_drag_force(speed=v + dv, altitude=h)
    weight = -interpol(mass_vector, 0, 1, pos, t) * drag.density.g0
    momentum = -(v + dv) * interpol(mass_vector, 2, 1, pos, t)
    acceleration = (thrust + drag_force + weight + momentum) / interpol(mass_vector, 0, 1, pos, t)

    # Return the acceleration
    return acceleration


def runge_kutta_4(thrust_vector, mass_vector, initial_conditions, drag):
    """
    -Thrust vector, defined as 2 line n column matrix where
         first line is the thrust(N) in a discrete time and second line is
         the time(s)
        -Mass vector, defined as 3 line n column matrix where
         first line is the mass(kg) at each time, second line is the time(s)
         and third line is the dm(kg/s) at each time
        -initial_position: initial conditions dictionary
        -drag: drag object used to calculate the drag exerted by the fluid
        -g0 is gravity (9,81m/s^2)
    It returns a vector with 4 lines and n columns where
        -First line is the velocity(m/s) at each time (line 0)
        -Second line is the time (line 1)
        -Third line is the acceleration(m/s^2 at each time (line 2)
        -Forth line is the altitude(m) at each time (line 3)
    """

    # Initialize the loop
    output = np.zeros((4, thrust_vector.shape[1]))
    pos = 0
    h = initial_conditions['h0']
    v = initial_conditions['v0']
    output[0, pos] = v
    output[1, :] = thrust_vector[1, :]

    # Run the loop within a try-catch block
    for time in thrust_vector[1, :][:-1]:

        try:
            v1 = v
            # Perform integration in time
            delta_t = thrust_vector[1, pos + 1] - time
            dv1 = delta_t * newton(thrust_vector, mass_vector, drag, pos, time, h, v, 0)
            dv2 = delta_t * newton(thrust_vector, mass_vector, drag, pos, time + delta_t / 2, h, v, dv1 / 2)
            dv3 = delta_t * newton(thrust_vector, mass_vector, drag, pos, time + delta_t / 2, h, v, dv2 / 2)
            dv4 = delta_t * newton(thrust_vector, mass_vector, drag, pos, time + delta_t / 2, h, v, dv3)
            v = v + (dv1 + 2 * dv2 + 2 * dv3 + dv4) / 6
            h = h + v * delta_t
            pos = pos + 1
            output[0, pos] = v  # store velocity
            output[2, pos] = (v - v1) / delta_t  # store acceleration
            output[3, pos] = h  # store altitude

            # Break the cycle when the rocket returns to ground
            if h < 0:
                break

        except (ArithmeticError, ValueError, ZeroDivisionError, IndexError):
            print("Error encountered while running RK4 loop. \n")

    return output


# -------------------- CLASS DEFINITIONS --------------------

class TrajectoryObject:
    """
    TrajectoryObject is in charge of handling all of the altitude computations for the rocket
    based on the Isp and Thrust values. The current class acts as a wrapper to ease integration
    into the simulation code and make it object oriented.

    Attributes:
        0.data_dictionary: dictionary containing all of the relevant data to the Trajectory
          calculator.
        1. drag: Drag instance object which determines the object drag (can be different from constant)
        2. density: Density instance which determines the local density at a given altitude.
        3. results: dictionary containing the results of the calculation made by the program
    """

    def __init__(self, density_obj, drag_obj):
        """
        class initializer
        :param density_obj: DensityLaw instance with the density calculation
        :param drag_obj: Drag instance
        """

        # First perform check on inputs
        assert isinstance(density_obj, DensityLaw), "Please insert a valid DensityLaw. \n"
        assert isinstance(drag_obj, Drag), "Please insert a valid Drag model. \n"

        # Allocate the attributes
        self.density = density_obj
        self.drag = drag_obj
        self.results = {'run_values': {'time': np.nan, 'altitude': np.nan, 'velocity': np.nan, 'acceleration': np.nan},
                        'magnitudes': {}}

    def __str__(self):
        """ redefine the str method for the TrajectoryObject """

        results_str = "\n".join(("{0}, {1:3.4f}".format(key, value)
                                 for key, value in self.results['magnitudes'].items()))

        # Return the output
        return "\n".join(["\nTrajectory,", str(self.drag), "\nResults,", results_str])

    def post_process_results(self):
        """
        post_process_results determines the magnitudes of interest that are to be outputted from the TrajectoryModule
        :return: nothing
        """

        self.results['magnitudes'] = {"max_{key}".format(key=key): np.nanmax(value) for key, value in
                                      self.results['run_values'].items() if key != 'time'}

    def run_simulation_on_trajectory(self, time, thrust, isp, initial_conditions):
        """
        run_simulation_on_trajectory performs the calculation of the rocket performance based on a given
        thrust and isp curve and initial_conditions.
        :param time: time array
        :param thrust: thrust array - mapping with time has to be one-2-one
        :param isp: isp array/value
        :param initial_conditions: dictionary with initial conditions {'h0', 'v0', 'm0'}
        :return: nothing
        """

        # Perform basic checks on the inputs
        assert len(time) == len(thrust), "Time and Thrust arrays do not have same length, please check inputs. \n"
        assert 'h0' in initial_conditions, "Please make sure initial conditions are properly set. \n"
        assert 'v0' in initial_conditions, "Please make sure initial conditions are properly set. \n"
        assert 'm0' in initial_conditions, "Please make sure initial conditions are properly set. \n"

        # Prepare inputs and call the simulation methods
        n_points = len(time)
        thrust_vector = np.empty((2, n_points))
        thrust_vector[0, :] = thrust
        thrust_vector[1, :] = time

        # Calculate the mass over time for the rocket
        mass_vector = mass_calc(thrust_vector, isp, self.drag.density.g0, initial_conditions['m0'])

        # Calculate acceleration, speed and altitude
        output = runge_kutta_4(thrust_vector, mass_vector, initial_conditions, self.drag)

        # Unpack output and store results
        self.results['run_values'] = {'time': output[1, :],
                                      'altitude': output[3, :],
                                      'velocity': output[0, :],
                                      'acceleration': output[2, :]}

        # Post-process results
        self.post_process_results()

    def plot_results(self):
        """
        plot_results plots the results from the simulation into a matplotlib figure
        :return: nothing
        """

        # Check the time array is not, if so then raise error, otherwise
        if self.results['run_values']["time"].any != np.nan:

            # Extract the concerning results
            time = self.results['run_values']["time"]
            altitude = self.results['run_values']["altitude"]
            velocity = self.results['run_values']["velocity"]
            acceleration = self.results['run_values']["acceleration"]
        else:
            raise ValueError("No values found for time, check results before plotting. \n")

        # Set the font dictionaries (for plot title and axis titles)
        title_font = {'size': '20', 'color': 'black', 'weight': 'normal',
                      'verticalalignment': 'bottom'}  # Bottom vertical alignment for more space
        axis_font = {'size': '16'}

        # Generate the plots
        fig = plt.figure(facecolor='w', figsize=(30, 30))
        fig.suptitle('Simulation results', **title_font)
        axs = [plt.subplot2grid((3, 2), (0, 0), rowspan=1, colspan=1),
               plt.subplot2grid((3, 2), (1, 0), rowspan=1, colspan=1),
               plt.subplot2grid((3, 2), (2, 0), rowspan=1, colspan=1),
               plt.subplot2grid((3, 2), (0, 1), rowspan=3, colspan=1)]

        # Set the tick labels font
        for ax in axs:
            for label in (ax.get_xticklabels() + ax.get_yticklabels()):
                label.set_fontname('Arial')
                label.set_fontsize(14)

        # Altitude-plot
        axs[0].plot(time, altitude, label='Altitude', color='green', linewidth=2.0)
        axs[0].set_title('')
        # axs[0].set_xlabel('time (s)')
        axs[0].set_ylabel('Altitude (m)', **axis_font)
        axs[0].grid(b=True, axis='both')
        axs[0].set_xlim(left=time[0])
        axs[0].set_ylim(bottom=altitude[0])

        # Port Velocity-plot
        axs[1].plot(time, velocity, linestyle='--', label='Velocity', color='red', linewidth=2.0)
        axs[1].set_title('')
        # axs[1].set_xlabel('time (s)')
        axs[1].set_ylabel('Vertical-Velocity (m/seg)', **axis_font)
        axs[1].grid(b=True, axis='both')
        axs[1].set_xlim(left=time[0])

        # Acceleration-plot
        axs[2].plot(time, acceleration, label='Acceleration', color='blue', linewidth=2.0)
        axs[2].set_title('')
        axs[2].set_xlabel('time (s)', **axis_font)
        axs[2].set_ylabel('Acceleration (m/s^2)', **axis_font)
        axs[2].grid(b=True, axis='both')
        axs[2].set_xlim(left=time[0])

        # Lift off-plot
        axs[3].plot(altitude, velocity, label='Lift off', color='black', linewidth=2.0)
        axs[3].set_title('')
        axs[3].set_xlabel('Altitude (m)', **axis_font)
        axs[3].set_ylabel('Speed (m/s)', **axis_font)
        axs[3].grid(b=True, axis='both')
        axs[3].set_xlim(left=0, right=10)
        axs[3].set_ylim(bottom=velocity[0], top=100)
        axs[3].set_xticks(np.linspace(0, 10, 11))
        axs[3].set_yticks(np.linspace(0, 100, 11))

    def return_results(self):
        """ return the results from the simulation """
        return self.results
