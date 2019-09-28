# InstabilityAnalysisScript.py intends to perform the instability analysis
# associated to the Rossiter modes in order to possibly estimate the behavior
# of the grain.
# @author: Jose Felix Zapata Usandivaras
# Date: 12/05/2019
# ISAE - SUPAERO, Scube GRIFFON, Combustion Team
#

# ---------------------------- IMPORT MODULES -------------------------------

from DataLayer.JsonInterpreter import JsonInterpreter           # Import the json interpreter
from CombustionModule.Combustion import *
import CombustionModule.RegressionModel as Reg                  # Import the RegressionModel module
import CombustionModule.Geometries as Geom                      # Import the Geometry module
import CombustionModule.Nozzle as Noz                           # Import the Nozzle module
from CombustionModule.Combustion import CombustionObject        # Import the CombustionObject
import numpy as np                                              # Import numpy
import matplotlib.pyplot as plt                                 # Import matplotlib
from CombustionModule.Rossiter import *                         # Import the Rossiter values


# -------------------------- FUNCTION DEFINITIONS ---------------------------

def generate_data_layer(data_file):
    """ generate_data_layer instantiates all of the objects that
    are required as data layer for the program to run properly.
    :param data_file: name of the json file to be used as reference for the simulation (str)
    :return JsonInterpreter """

    # Pass the file name
    data_directory = "../data"

    # Return the output
    return JsonInterpreter(file_name="/".join([data_directory, data_file]))


def define_onera_physical_hycom_13_test_combustion_object():
    """ Define the ONERA physical test CombustionObject for Hycom 13 test
    :return combustion object
    """

    # ------------ Generate the data layer:

    json_interpreter = generate_data_layer("DataLayerONERATests_Hycom13.json")
    # Use same data-layer used for Griffon (same pressure)
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'L': 0.130,
                        'rintInitial': 0.015,
                        'rext0': 94.5e-3 / 2,
                        'regressionModel': Reg.TwoRegimesMarxmanAndFloodedModel(**combustion_table)}

    nozzle_params = {'At': 4.03764e-5, 'expansion': 6.3, 'lambda_e': 0.98, 'erosion': 0}

    simulation_params = {'ox_flow': 0.079, 'safety_thickness': 0.005, 'dt': 0.01, 'max_burn_time': 8.0  }

    # ------------- Generate objects:

    geometry_obj = Geom.OneCircularPort(**geometric_params)
    nozzle_obj = Noz.Nozzle(**nozzle_params)

    # Instantiate the combustion module
    combustion_object = CombustionObjectClassic(json_interpreter=json_interpreter,
                                                geometry_object=geometry_obj,
                                                nozzle_object=nozzle_obj)

    # -------------- Run simulation & Plot:

    combustion_object.run_simulation_constant_fuel_sliver(**simulation_params)

    # Return the combustion object
    return combustion_object


def define_rossiter_modes_hycom_13_test(combustion_module):
    """
    define_rossiter_modes generates de objects associated to the calculation of the rossiter modes
    :param combustion_module: CombustionObject instance with already run_values populated
    :return: nothing
    """

    # -------------------- Generate Post-Combustion Object:

    post_combustion_params = {'L': 0.177,        # Length of post-combustion chamber
                              'D': 0.047}        # Diameter of post-combustion chamber

    post_combustion_obj = PostCombustion(**post_combustion_params)


    # --------------------- Define a mode start_time for firing synchronization:

    start_time = 3.0            # Mode appears at 2.5 seconds
    time_shift = 2.5            # Time shift between test data and simulation
    f0 = 820                    # Initial frequency of the mode [Hz]

    # --------------------- Define modes to be instantiated:

    alpha = 0                                     # Define the constant alpha for all modes
    modes_params = [{'m_ross': 20, 'alpha': alpha, 'start_time': 3.0, 'f0': 820},
                    {'m_ross': 19, 'alpha': alpha, 'start_time': 0.75, 'f0': 775},
                    {'m_ross': 18, 'alpha': alpha, 'start_time': start_time, 'f0': f0},
                    {'m_ross': 15, 'alpha': alpha, 'start_time': start_time, 'f0': f0},
                    {'m_ross': 10, 'alpha': alpha, 'start_time': start_time, 'f0': f0}]


    # --------------------- Define the test file:

    data_directory = "../data"
    data_file = "Hycom_13_Rossiter_extract.txt"
    test_file_name_complete = "/".join([data_directory, data_file])
    test_name = "Hycom_13"


    # ---------------------- Instantiate the objects:

    rossiter_obj = Rossiter(post_combustion=post_combustion_obj,
                            modes_params=modes_params,
                            combustion_module=combustion_module,
                            test_file_name=test_file_name_complete,
                            test_name=test_name,
                            time_shift=time_shift)

    # ----------------------- Perform the modes calibration:

    rossiter_obj.calibrate_modes()


    # ----------------------- Calculate the frequencies:

    rossiter_obj.calculate_modes_frequencies()


    # ----------------------- Plot results:

    rossiter_obj.plot_results()


def define_onera_physical_hycom_16_test_combustion_object():
    """ Define the ONERA physical test CombustionObject for Hycom 16 test
    :return combustion object
    """

    # ------------ Generate the data layer:

    json_interpreter = generate_data_layer("DataLayerONERATests_Hycom16.json")
    # Use same data-layer used for Griffon (same pressure)
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    shape_params = {'polynom': [-0.2701822916666667, 0.15625, 0],
                    'baseRadius': 0.02,
                    'branches': 5,
                    'n': 50}

    geometric_params = {'L': 0.130,
                        'externalRadius': 74.5 / 2000,
                        'imagePixelSize': 2048,
                        'imageMeterSize': 0.09,
                        'regressionModel': Reg.TwoRegimesMarxmanAndFloodedModel(**combustion_table)}

    nozzle_params = {'At': 4.243e-5, 'expansion': 6.3, 'lambda_e': 0.98, 'erosion': 0}

    simulation_params = {'ox_flow': 0.0975, 'safety_thickness': 0.005, 'max_burn_time': 7.5}

    # ------------- Generate objects:

    geometry_obj = Geom.SinglePortImageGeometry(**geometric_params)
    geometry_obj.generatePolynom(**shape_params)
    nozzle_obj = Noz.Nozzle(**nozzle_params)

    # Instantiate the combustion module
    combustion_object = CombustionObjectClassic(json_interpreter=json_interpreter,
                                                geometry_object=geometry_obj,
                                                nozzle_object=nozzle_obj)

    # -------------- Run simulation & Plot:

    combustion_object.run_simulation_constant_fuel_sliver_image_geometry(**simulation_params)

    # Return the combustion object
    return combustion_object


def define_rossiter_modes_hycom_16_test(combustion_module):
    """
    define_rossiter_modes generates de objects associated to the calculation of the rossiter modes
    :param combustion_module: CombustionObject instance with already run_values populated
    :return: nothing
    """

    # -------------------- Generate Post-Combustion Object:

    post_combustion_params = {'L': 0.135,        # Length of post-combustion chamber
                              'D': 0.09}        # Diameter of post-combustion chamber

    post_combustion_obj = PostCombustion(**post_combustion_params)


    # --------------------- Define a mode start_time for firing synchronization:

    start_time = 3.0            # Mode appears at 2.5 seconds
    time_shift = 2.5            # Time shift between test data and simulation
    f0 = 930                    # Initial frequency of the mode [Hz]

    # --------------------- Define modes to be instantiated:

    alpha = 0                                     # Define the constant alpha for all modes
    modes_params = [{'m_ross': 20, 'alpha': alpha, 'start_time': start_time, 'f0': f0},
                    {'m_ross': 19, 'alpha': alpha, 'start_time': start_time, 'f0': f0},
                    {'m_ross': 18, 'alpha': alpha, 'start_time': start_time, 'f0': f0},
                    {'m_ross': 15, 'alpha': alpha, 'start_time': start_time, 'f0': f0},
                    {'m_ross': 10, 'alpha': alpha, 'start_time': start_time, 'f0': f0}]


    # --------------------- Define the test file:

    data_directory = "../data"
    data_file = "Hycom_16_Rossiter_extract.txt"
    test_file_name_complete = "/".join([data_directory, data_file])
    test_name = "Hycom_16"


    # ---------------------- Instantiate the objects:

    rossiter_obj = Rossiter(post_combustion=post_combustion_obj,
                            modes_params=modes_params,
                            combustion_module=combustion_module,
                            test_file_name=test_file_name_complete,
                            test_name=test_name,
                            time_shift=time_shift)

    # ----------------------- Perform the modes calibration:

    rossiter_obj.calibrate_modes()


    # ----------------------- Calculate the frequencies:

    rossiter_obj.calculate_modes_frequencies()


    # ----------------------- Plot results:

    rossiter_obj.plot_results()


def define_onera_physical_hycom_12_test_combustion_object():
    """ Define the ONERA physical test CombustionObject for Hycom 12 test
    :return combustion object
    """

    # ------------ Generate the data layer:

    json_interpreter = generate_data_layer("DataLayerONERATests_Hycom12.json")
    # Use same data-layer used for Griffon (same pressure)
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'L': 0.130,
                        'rintInitial': 17.23e-3 / 2,
                        'rext0': 94.5e-3 / 2,
                        'regressionModel': Reg.TwoRegimesMarxmanAndFloodedModel(**combustion_table)}

    nozzle_params = {'At': (np.pi / 4) * 7.1e-3 ** 2, 'expansion': 6.3, 'lambda_e': 0.98, 'erosion': 0}

    simulation_params = {'ox_flow': 0.077, 'safety_thickness': 0.005, 'dt': 0.01, 'max_burn_time': 8.0  }

    # ------------- Generate objects:

    geometry_obj = Geom.OneCircularPort(**geometric_params)
    nozzle_obj = Noz.Nozzle(**nozzle_params)

    # Instantiate the combustion module
    combustion_object = CombustionObject(json_interpreter=json_interpreter,
                                         geometry_object=geometry_obj,
                                          nozzle_object=nozzle_obj)

    # -------------- Run simulation & Plot:

    combustion_object.run_simulation_constant_fuel_sliver(**simulation_params)

    # Return the combustion object
    return combustion_object


def define_rossiter_modes_hycom_12_test(combustion_module):
    """
    define_rossiter_modes generates de objects associated to the calculation of the rossiter modes
    :param combustion_module: CombustionObject instance with already run_values populated
    :return: nothing
    """

    # -------------------- Generate Post-Combustion Object:

    post_combustion_params = {'L': 0.135,        # Length of post-combustion chamber
                              'D': 0.09}        # Diameter of post-combustion chamber

    post_combustion_obj = PostCombustion(**post_combustion_params)


    # --------------------- Define a mode start_time for firing synchronization:

    start_time = 0.75            # Mode appears at 2.5 seconds
    time_shift = 2.5            # Time shift between test data and simulation
    f0 = 400                    # Initial frequency of the mode [Hz]

    # --------------------- Define modes to be instantiated:

    alpha = 0                                     # Define the constant alpha for all modes
    modes_params = [{'m_ross': 7, 'alpha': alpha, 'start_time': start_time, 'f0': f0},
                    {'m_ross': 8, 'alpha': alpha, 'start_time': start_time, 'f0': f0},
                    {'m_ross': 9, 'alpha': alpha, 'start_time': start_time, 'f0': f0},
                    {'m_ross': 10, 'alpha': alpha, 'start_time': start_time, 'f0': f0},
                    {'m_ross': 11, 'alpha': alpha, 'start_time': start_time, 'f0': f0}]


    # --------------------- Define the test file:

    data_directory = "../data"
    data_file = "Hycom_12_Rossiter_extract.txt"
    test_file_name_complete = "/".join([data_directory, data_file])
    test_name = "Hycom_12"


    # ---------------------- Instantiate the objects:

    rossiter_obj = Rossiter(post_combustion=post_combustion_obj,
                            modes_params=modes_params,
                            combustion_module=combustion_module,
                            test_file_name=test_file_name_complete,
                            test_name=test_name,
                            time_shift=time_shift)

    # ----------------------- Perform the modes calibration:

    rossiter_obj.calibrate_modes()


    # ----------------------- Calculate the frequencies:

    rossiter_obj.calculate_modes_frequencies()


    # ----------------------- Plot results:

    rossiter_obj.plot_results()


def define_onera_physical_hycom_14_test_combustion_object():
    """ Define the ONERA physical test CombustionObject for Hycom 14 test
    :return combustion object
    """

    # ------------ Generate the data layer:

    json_interpreter = generate_data_layer("DataLayerONERATests_Hycom14.json")
    # Use same data-layer used for Griffon (same pressure)
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'L': 0.130,
                        'rintInitial': 13.30e-3 / 2,
                        'rext0': 94.5e-3 / 2,
                        'regressionModel': Reg.TwoRegimesMarxmanAndFloodedModel(**combustion_table)}

    nozzle_params = {'At': (np.pi / 4) * 7.25e-3 ** 2, 'expansion': 6.3, 'lambda_e': 0.98, 'erosion': 0}

    simulation_params = {'ox_flow': 0.080, 'safety_thickness': 0.005, 'dt': 0.01, 'max_burn_time': 8.5  }

    # ------------- Generate objects:

    geometry_obj = Geom.OneCircularPort(**geometric_params)
    nozzle_obj = Noz.Nozzle(**nozzle_params)

    # Instantiate the combustion module
    combustion_object = CombustionObject(json_interpreter=json_interpreter,
                                         geometry_object=geometry_obj,
                                          nozzle_object=nozzle_obj)

    # -------------- Run simulation & Plot:

    combustion_object.run_simulation_constant_fuel_sliver(**simulation_params)

    # Return the combustion object
    return combustion_object



# ------------------------------ MAIN ----------------------------------

if __name__ == '__main__':

    # Define the test name
    test_name = 'hycom_13'
    # Call the methods
    if test_name == 'hycom_12':
        combustion_obj = define_onera_physical_hycom_12_test_combustion_object()
        define_rossiter_modes_hycom_12_test(combustion_obj)
    elif test_name == 'hycom_13':
        combustion_obj = define_onera_physical_hycom_13_test_combustion_object()
        define_rossiter_modes_hycom_13_test(combustion_obj)
    elif test_name == 'hycom_14':
        combustion_obj = define_onera_physical_hycom_14_test_combustion_object()
    elif test_name == 'hycom_16':
        combustion_obj = define_onera_physical_hycom_16_test_combustion_object()
        define_rossiter_modes_hycom_16_test(combustion_obj)

    # Plot the combustion obj
    combustion_obj.plot_results()
    print(combustion_obj)

    # Show any plots
    plt.show()
