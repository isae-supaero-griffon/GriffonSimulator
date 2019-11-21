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
                        'r_init':0.152/2,
                        'r_ext': 0.5,
                        'N': 12}

    nozzle_params = {'At': 0.000589, 'expansion': 5.7, 'lambda_e': 0.98, 'erosion': 0}

    simulation_params = {'ox_flow': 7.95, 'safety_thickness': 0.005, 'dt': 0.1, 'max_burn_time': 70,
                         'file_name': '../data/data_tests/myGeometry.txt',
                         'times': [0.1, 20, 60]}

    # ------------- Generate objects:

    geometry_obj = Geom.SingleCircularPort1D(**geometric_params)
    nozzle_obj = Noz.Nozzle(**nozzle_params)
    regression_mod = Reg.SingleRegimeMarxmanModel(**combustion_table)

    # Instantiate the combustion module
    combustion_obj = Combustion.CombustionObject1D(json_interpreter=json_interpreter,
                                                   nozzle_object=nozzle_obj,
                                                   geometry_object=geometry_obj,
                                                   regression_model=regression_mod)

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

    json_interpreter = generate_data_layer("Griffon Data - ABS - H2O2 - 36 bar.json")
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'L': 0.4,
                        'r_init':0.03,
                        'r_ext': 0.05,
                        'N': 200}

    nozzle_params = {'At': 0.000589, 'expansion': 5.7, 'lambda_e': 0.98, 'erosion': 0}

    simulation_params = {'ox_flow': 1.43, 'safety_thickness': 0.005, 'dt': 0.01,
                         'max_burn_time': None, 'tol_press': 1e-1}

    # ------------- Generate objects:

    geometry_obj = Geom.SingleCircularPort1D(**geometric_params)
    nozzle_obj = Noz.Nozzle(**nozzle_params)
    regression_mod = Reg.MarxmanAndConstantFloodingRegimeModel(**combustion_table)

    # Instantiate the combustion module
    combustion_obj = Combustion.CombustionObject1D(json_interpreter=json_interpreter,
                                                   nozzle_object=nozzle_obj,
                                                   geometry_object=geometry_obj,
                                                   regression_model=regression_mod)

    # -------------- Run simulation & Plot:

    combustion_obj.run_balanced_nozzle_analysis(**simulation_params)

    # Print the module
    print(combustion_obj)

    # Plot the results
    combustion_obj.plot_results()

    # Show plots
    plt.show()


def test_hydraulic_table_generation():
    """ test the hydraulic table extraction by JsonInterpreter """

    # ------------ Generate the data-layer:

    json_interpreter = generate_data_layer("Griffon Data - Mock.json")
    hydraulic_table = json_interpreter.return_hydraulic_table()

    # Print the dictionary
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(hydraulic_table)


def test_hydraulic_module_initiation():
    """ test if the hydraulic module is assembled accordingly and in the right
    way from the hydraulic_table declaration """

    # ------------ Generate the data-layer:

    json_interpreter = generate_data_layer("Griffon Data - Mock.json")
    hydraulic_table = json_interpreter.return_hydraulic_table()

    # Generate the hydraulic module
    my_module = HydraulicModule(hydraulic_table)

    # Set the exit chamber pressure
    my_module.set_chamber_pressure(value=3600000)

    # Initialize the pressure dofs
    my_module.initialize_pressure_dofs()

    # print the module dofs
    print(my_module.dof_print())

    # print the module nodes
    print(my_module.node_print())

    # print the module
    print(my_module)


def test_hydraulic_module_network_solution():
    """ test if the hydraulic module solves properly for the network """

    # ------------ Generate the data-layer:

    json_interpreter = generate_data_layer("Griffon Data - Mock.json")
    hydraulic_table = json_interpreter.return_hydraulic_table()

    # Generate the hydraulic module
    my_module = HydraulicModule(hydraulic_table)

    # Set the pressure of the oxidizer tank
    my_dof = my_module.dofs.return_element(criteria=6)
    my_dof.fix()
    my_dof.set_value(56e5)

    # Set the exit chamber pressure
    my_module.set_chamber_pressure(value=3600000)

    # Initialize the pressure dofs
    my_module.initialize_pressure_dofs()

    # Solve the module
    sol = my_module.run_simulation(tol=1e-6, maxiter=500, scale=1e6)

    # print the module dofs
    print(my_module.dof_print())

    # print the module nodes
    print(my_module.node_print())

    # print the module
    print(my_module)


    # Print sol
    # print(sol)


def generate_injector_json_table():
    """ generate the table associated to the injector in a json file """

    file_name = "../data/Config Tables/BIM49-Curve.txt"
    values = {'table': []}

    with open(file_name, 'r') as f:
        for line in f:
            v = tuple(float(elem) for elem in line.split('\t'))
            values['table'].append(v)

    with open(file_name[:-4] + ".json", 'w') as f:
        json.dump(values, f)

    pprint.pprint(values)



# ------------------------- MAIN ------------------------ #

if __name__ == '__main__':

    # Call the test method
    test_1d_implementation_single_circular_port()
    # test_1d_implementation_single_circular_port_balanced_nozzle()
    # test_hydraulic_table_generation()
    # test_hydraulic_module_initiation()
    # test_hydraulic_module_network_solution()
    # generate_injector_json_table()