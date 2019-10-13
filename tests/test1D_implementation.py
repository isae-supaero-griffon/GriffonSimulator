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
import matplotlib.pyplot as plt                               # Import matplotlib

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

    json_interpreter = generate_data_layer("Griffon Data - ABS - H2O2 - 36 bar.json")
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'L': 0.4,
                        'r_init':0.03,
                        'r_ext': 0.05,
                        'N': 200}

    nozzle_params = {'At': 0.000589, 'expansion': 5.7, 'lambda_e': 0.98, 'erosion': 0}

    simulation_params = {'ox_flow': 1.01, 'safety_thickness': 0.005, 'dt': 0.01, 'max_burn_time': 5.0  }

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

    combustion_obj.run_simulation_constant_fuel_sliver(**simulation_params)

    # Print the module
    print(combustion_obj)

    # Plot the results
    combustion_obj.plot_results()

    # Show plots
    plt.show()


# ------------------------- MAIN ------------------------ #

if __name__ == '__main__':

    # Call the test method
    test_1d_implementation_single_circular_port()