# testsOnera.py is a script specifically written to test different configurations to be later
# sent to the test-bench to obtain sensible data for Griffon.
# Author: Jose Felix Zapata Usandivaras
# Date: 5/3/2019
# ISAE - SUPAERO / Griffon Project - Scube

# ------------------------ IMPORT MODULES ------------------------------

from DataLayer.JsonInterpreter import JsonInterpreter           # Import the json interpreter
import CombustionModule.RegressionModel as Reg                  # Import the RegressionModel module
import CombustionModule.Geometries as Geom                      # Import the Geometry module
import CombustionModule.Nozzle as Noz                           # Import the Nozzle module
from CombustionModule.Combustion import CombustionObject        # Import the CombustionObject
from MassEstimationModule.system import System                  # Import the system class
from TrajectoryModule.Drag import *                             # Import the Drag library
from TrajectoryModule.Atmosphere import DensityLaw                 # Import the density-law library
from TrajectoryModule.Trajectory import TrajectoryObject        # Import the trajectory object
import numpy as np                                              # Import numpy
import matplotlib.pyplot as plt                                 # Import matplotlib

# ------------------------ FUNCTION DEFINITIONS ------------------------


def generate_data_layer(data_file="Thermodynamic Data 36 bar OF 0,1 to 8,0 H2O2 87,5.json"):
    """ generate_data_layer instantiates all of the objects that
    are required as data layer for the program to run properly.
    :param data_file: name of the json file to be used as reference for the simulation (str)
    :return JsonInterpreter """

    # Pass the file name
    data_directory = "../data"

    # Return the output
    return JsonInterpreter(file_name="/".join([data_directory, data_file]))


def test_onera_physical_test_1():
    """ Study the potential evolution of the geometry associated to the 1st physical test
        to be performed at ONERA. The characteristics of the Test are:
        1. Single Port Geometry
        2. Go range: [100 - etc] kg/m^2/sec
        3. Chamber pressure: 41 bar, closest possible to the target pressure of Griffon.
    """

    # ------------ Generate the data layer:

    json_interpreter = generate_data_layer(data_file="Thermodynamic Data Onera 41 bar H2O2 87_5.json")
    # Use same data-layer used for Griffon (same pressure)
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'length': 0.157,
                        'r_int_initial': 0.02,
                        'r_ext': 0.041,
                        'regression_model': Reg.MarxmanAndConstantFloodingRegimeModel(**combustion_table)}

    nozzle_params = {'At': 0.000038, 'expansion': 6.3, 'lambda_e': 0.98, 'erosion': 0}

    simulation_params = {'ox_flow': 0.08, 'safety_thickness': 0.005, 'dt': 0.01, 'max_burn_time': 8}

    # ------------- Generate objects:

    geometry_obj = Geom.OneCircularPort(**geometric_params)
    nozzle_obj = Noz.Nozzle(**nozzle_params)

    # Instantiate the combustion module
    combustion_obj = CombustionObject(json_interpreter=json_interpreter,
                                      geometry_object=geometry_obj,
                                      nozzle_object=nozzle_obj)

    # -------------- Run simulation & Plot:

    combustion_obj.run_simulation_constant_fuel_sliver(**simulation_params)

    # Print the module
    print(combustion_obj)

    # Plot the results
    combustion_obj.plot_results()


def test_onera_physical_test_1_2():
    """ Study the potential evolution of the geometry associated to the 1st physical test
        to be performed at ONERA (version at 36 bar). The characteristics of the Test are:
        1. Single Port Geometry
        2. Go range: [100 - etc] kg/m^2/sec
        3. Chamber pressure: TBD, closest possible to the target pressure of Griffon.
    """

    # ------------ Generate the data layer:

    json_interpreter = generate_data_layer()
    # Use same data-layer used for Griffon (same pressure)
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'length': 0.130,
                        'r_int_initial': 0.015,
                        'r_ext': 94.5e-3 / 2,
                        'regression_model': Reg.MarxmanAndConstantFloodingRegimeModel(**combustion_table)}

    nozzle_params = {'At': 0.000038, 'expansion': 6.3, 'lambda_e': 0.98, 'erosion': 0}

    simulation_params = {'ox_flow': 0.07, 'safety_thickness': 0.005, 'dt': 0.01, 'max_burn_time': 8}

    # ------------- Generate objects:

    geometry_obj = Geom.OneCircularPort(**geometric_params)
    nozzle_obj = Noz.Nozzle(**nozzle_params)

    # Instantiate the combustion module
    combustion_obj = CombustionObject(json_interpreter=json_interpreter,
                                      geometry_object=geometry_obj,
                                      nozzle_object=nozzle_obj)

    # -------------- Run simulation & Plot:

    combustion_obj.run_simulation_constant_fuel_sliver(**simulation_params)

    # Print the module
    print(combustion_obj)

    # Plot the results
    combustion_obj.plot_results()


def test_onera_physical_test_2():
    """ Study the potential evolution of the geometry associated to the 2nd physical test
        to be performed at ONERA. The characteristics of the Test are:
        1. Single Port Geometry
        2. Go range: [100, 300] kg/m^2/sec
        3. Chamber pressure: 36 bar, closest possible to the target pressure of Griffon.
    """

    # ------------ Generate the data layer:

    json_interpreter = generate_data_layer()                    # Use same data-layer used for Griffon (same pressure)
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'length': 0.130,
                        'r_int_initial': 17.23e-3 / 2,
                        'r_ext': 94.5e-3 / 2,
                        'regression_model': Reg.MarxmanAndConstantFloodingRegimeModel(**combustion_table)}

    nozzle_params = {'At': 0.000038, 'expansion': 6.3, 'lambda_e': 0.98, 'erosion': 0}

    simulation_params = {'ox_flow': 0.07, 'safety_thickness': 0.005, 'dt': 0.01, 'max_burn_time': 8}

    # ------------- Generate objects:

    geometry_obj = Geom.OneCircularPort(**geometric_params)
    nozzle_obj = Noz.Nozzle(**nozzle_params)

    # Instantiate the combustion module
    combustion_obj = CombustionObject(json_interpreter=json_interpreter,
                                      geometry_object=geometry_obj,
                                      nozzle_object=nozzle_obj)

    # -------------- Run simulation & Plot:

    combustion_obj.run_simulation_constant_fuel_sliver(**simulation_params)

    # Print the module
    print(combustion_obj)

    # Plot the results
    combustion_obj.plot_results()


def test_onera_physical_test_3():
    """ Study the potential evolution of the geometry associated to the 2nd physical test
        to be performed at ONERA. The characteristics of the Test are:
        1. Single Port Geometry
        2. Go range: [400, 500] kg/m^2/sec
        3. Chamber pressure: 36 bar, closest possible to the target pressure of Griffon.
    """

    # ------------ Generate the data layer:

    json_interpreter = generate_data_layer()                    # Use same data-layer used for Griffon (same pressure)
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'length': 0.130,
                        'r_int_initial': 13.33e-3 / 2,
                        'r_ext': 94.5e-3 / 2,
                        'regression_model': Reg.MarxmanAndConstantFloodingRegimeModel(**combustion_table)}

    nozzle_params = {'At': 0.000038, 'expansion': 6.3, 'lambda_e': 0.98, 'erosion': 0}

    simulation_params = {'ox_flow': 0.07, 'safety_thickness': 0.005, 'dt': 0.01, 'max_burn_time': 7}

    # ------------- Generate objects:

    geometry_obj = Geom.OneCircularPort(**geometric_params)
    nozzle_obj = Noz.Nozzle(**nozzle_params)
    json_interpreter = generate_data_layer(data_file="Thermodynamic Data Onera 41 bar H2O2 87_5.json")

    # Instantiate the combustion module
    combustion_obj = CombustionObject(json_interpreter=json_interpreter,
                                      geometry_object=geometry_obj,
                                      nozzle_object=nozzle_obj)

    # -------------- Run simulation & Plot:

    combustion_obj.run_simulation_constant_fuel_sliver(**simulation_params)

    # Print the module
    print(combustion_obj)

    # Plot the results
    combustion_obj.plot_results()


def test_onera_three_port_geometry_with_center_port():
    """ Study the potential evolution of the geometry associated to a Three Port Geometry
        to be performed at ONERA. The characteristics of the Test are:
        1. Three Port Geometry
        2. Go range: TBD kg/m^2/sec
        3. ox flow: same
        4. Chamber pressure: 36 bar, closest possible to the target pressure of Griffon.
    """

    # ------------ Generate the data layer:

    json_interpreter = generate_data_layer()  # Use same data-layer used for Griffon (same pressure)
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'nodes_number': 4,
                        'length': 0.130,
                        'ring_ports_initial_radius': 9e-3 / 2,
                        'central_ports_initial_radius': 9e-3 / 2,
                        'r_ext': 94.5e-3 / 2,
                        'regression_model': Reg.MarxmanAndConstantFloodingRegimeModel(**combustion_table)}

    nozzle_params = {'At': 0.000038, 'expansion': 6.3, 'lambda_e': 0.98, 'erosion': 0}

    simulation_params = {'ox_flow': 0.065, 'safety_thickness': 0.005, 'dt': 0.01, 'max_burn_time': 8}

    # ------------- Generate objects:

    geometry_obj = Geom.MultipleCircularPortsWithCircularCenter(**geometric_params)
    geometry_obj.draw_geometry()
    nozzle_obj = Noz.Nozzle(**nozzle_params)
    json_interpreter = generate_data_layer(data_file="Thermodynamic Data Onera 41 bar H2O2 87_5.json")

    # Instantiate the combustion module
    combustion_obj = CombustionObject(json_interpreter=json_interpreter,
                                      geometry_object=geometry_obj,
                                      nozzle_object=nozzle_obj)

    # -------------- Run simulation & Plot:

    combustion_obj.run_simulation_constant_fuel_sliver(**simulation_params)

    # Print the module
    print(combustion_obj)

    # Plot the results
    combustion_obj.plot_results()


def compute_sampling_frequency_for_hycomm():
    """
    compute_sampling_frequency_for_hycomm determines the sampling frequency required based
    on the sonic speed expected at 36 bar at different O/F ratios.
    :return: nothing
    """

    # Define the parameters for the engine
    L = 240e-3                          # Length of the engine in m
    n = 4                               # Coefficient associated to the acoustic mode to be observe
    uncertainty = 0.30                  # Uncertainty coefficient
    sampling_coefficient = 9           # Ratio between the sampling frequency and the nyquist frequency

    # Collect the data from the data layer generated for 36 bars
    json_interpreter = generate_data_layer()
    cea_data = json_interpreter.return_cea_data(CombustionObject.simulation_type)
    of_data = np.asarray([record['state']['o/f'] for record in cea_data])
    sonic_velocity = np.asarray([record['variables']['sonvel']['CHAMBER'] for record in cea_data])
    c_star = np.asarray([record['variables']['cstar']['THROAT'] for record in cea_data])

    # Get the max value of the sonic velocity
    max_c = max(sonic_velocity)

    # Determine the max frequency that we'll have
    max_freq = (1 + uncertainty) * max_c / (n * L)
    nyquist_frequency = 2*max_freq
    sampling_frequency = sampling_coefficient * nyquist_frequency

    # Print the results
    print("The maximum Speed of Sound observed is {0} m/s".format(max_c))
    print("The maximum estimated frequency for the mode is {0:.2E} Hz".format(max_freq))
    print("The Nyquist Frequency for the mode is {0:.2E} Hz".format(nyquist_frequency))
    print("The Sampling Frequency Suggested is {0:.2E} Hz".format(sampling_frequency))

    # Plot the of and the sonic velocity
    fig = plt.figure(facecolor='w', figsize=(10, 10))
    fig.suptitle('OF vs. C')
    ax = plt.axes()
    ax.plot(of_data, sonic_velocity, label='Sonic Velocity', color='black', linewidth=2.0)
    ax.set_title('')
    ax.set_ylabel('Sonic Velocity (m/s)')
    ax.set_xlabel('O/F ratio')
    ax.grid(b=True, axis='both')


# ---------------------------- MAIN  ---------------------------------

if __name__ == '__main__':

    plt.close()
    # Call on desired test method
    # test_onera_physical_test_1()
    # test_onera_physical_test_1_2()
    # test_onera_physical_test_2()
    # test_onera_physical_test_3()
    # test_onera_three_port_geometry_with_center_port()
    compute_sampling_frequency_for_hycomm()
    # Show any plots
    plt.show()
