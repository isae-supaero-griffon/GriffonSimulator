# test1D_implementation.py script has for aim test the 1D CombustionModule implementation.
# @author: Jose Felix Zapata Usandivaras
# Date: 10/10/2019
# ISAE - SUPAERO Space Section, Project Griffon.

# ------------------------- IMPORT MODULES --------------------------- #

from DataLayer.JsonInterpreter import JsonInterpreter         # Import the Data Layer
import CombustionModule.Combustion as Combustion              # Import only the CombustionObject1D
import CombustionModule.Geometries as Geom                    # Import only the SingleCircularPort1D Geometry
import CombustionModule.RegressionModel as Reg                # Import all the possible regression models
import CombustionModule.Nozzle as Noz                         # Import the Nozzle module
from HydraulicModule.Hydraulic import HydraulicModule         # Import the Hydraulic module class
import matplotlib.pyplot as plt                               # Import matplotlib
import os.path                                                # Import os.path
import pprint                                                 # Import pprint
import json                                                   # Import json
import numpy as np                                            # Import numpy library
from copy import deepcopy                                     # Import the deepcopy method


# ----------------------- FUNCTION DEFINITIONS ----------------------- #

def generate_data_layer(data_file="Griffon Data - ABS - H2O2 - 36 bar.json"):
    """ generate_data_layer instantiates all of the objects that
    are required as data layer for the program to run properly.
    :param data_file: name of the json file to be used as reference for the simulation (str)
    :return JsonInterpreter """

    # Pass the file name
    data_directory = "../data"

    # Return the output
    return JsonInterpreter(file_name="/".join([data_directory, data_file]))


def test_1d_implementation_single_circular_port():
    """ 1D test of the Combustion Module 1D implementation """

    # ------------ Generate the data-layer:

    json_interpreter = generate_data_layer("Griffon Data - Mock.json")
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'L': 4.572,
                        'r_init': 0.152/2,
                        'r_ext': 0.5,
                        'N': 12}

    nozzle_params = {'At': 0.000589, 'expansion': 5.7, 'lambda_e': 0.98, 'erosion': 0}

    simulation_params = {'ox_flow': 7.95, 'safety_thickness': 0.005, 'dt': 0.1, 'max_burn_time': 70,
                         'file_name': '../data/data_tests/myGeometry.txt',
                         'times': [0.1, 20, 60]}

    # ------------- Generate objects:
    geometric_params['regressionModel'] = Reg.TwoRegimesMarxmanAndFloodedModel(**combustion_table)
    geometry_obj = Geom.SingleCircularPort1D(**geometric_params)
    nozzle_obj = Noz.Nozzle(**nozzle_params)

    # Instantiate the combustion module
    combustion_obj = Combustion.CombustionObject1D(json_interpreter=json_interpreter,
                                                   nozzle_object=nozzle_obj,
                                                   geometry_object=geometry_obj)

    # -------------- Run simulation & Plot:

    combustion_obj.run_simulation_constant_fuel_sliver(**simulation_params)

    # Print the module
    print(combustion_obj)

    # Plot the results
    combustion_obj.plot_results()

    # Show plots
    plt.show()


def test_1d_implementation_single_circular_port_balanced_nozzle():
    """ test_1d_implementation_single_circular_port_balanced_nozzle intends
    to test the implementation of the balanced nozzle based on the fixed point
    method. """

    # ------------ Generate the data-layer:

    json_interpreter = generate_data_layer("Griffon Data - Mock.json")
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'L': 0.4,
                        'r_init': 0.052/2,
                        'r_final': 0.094/2,
                        'r_exit': 0.08/2,
                        'exit_depth': 0.02,
                        'depth': 0.3,
                        'r_ext': 0.05,
                        'N': 200}

    nozzle_params = {'At': 0.000589, 'expansion': 5.7, 'lambda_e': 0.98, 'erosion': 0}

    simulation_params = {'ox_flow': 1.08, 'safety_thickness': 0.004, 'dt': 0.05,
                         'max_burn_time': 4.6, 'tol_press': 1e-3}

    # ------------- Generate objects:
    geometric_params['regressionModel'] = Reg.TwoRegimesMarxmanAndFloodedModel(**combustion_table)
    geometry_obj = Geom.ConicCircularPort1D(**geometric_params)
    nozzle_obj = Noz.Nozzle(**nozzle_params)

    # -------------- Generate deep copy of geometry object:

    geometry_object_original = deepcopy(geometry_obj)

    # Instantiate the combustion module
    combustion_obj = Combustion.CombustionObject1D(json_interpreter=json_interpreter,
                                                   nozzle_object=nozzle_obj,
                                                   geometry_object=geometry_obj)

    # -------------- Run simulation & Plot:

    combustion_obj.run_balanced_nozzle_analysis(**simulation_params)

    # -------------- Export the results to File:

    geometry_obj.export_geometry(file_name="C:/Users/Felix Zapata/Google Drive/MSc in AE ISAE - SUPAERO/"
                                                       "Griffon Project/Design Files/1D March Test Design/"
                                                       "geometry_after_burn_1D.txt")
    geometry_object_original.export_geometry(file_name="C:/Users/Felix Zapata/Google Drive/MSc in AE ISAE - SUPAERO/"
                                                       "Griffon Project/Design Files/1D March Test Design/"
                                                       "original_geometry_1D.txt")

    # Print the module
    print(combustion_obj)

    # Plot the results
    combustion_obj.plot_results()

    # Show plots
    plt.show()

# ------------------------- MAIN ------------------------ #


if __name__ == '__main__':

    # Call the test method
    # test_1d_implementation_single_circular_port()
    test_1d_implementation_single_circular_port_balanced_nozzle()
