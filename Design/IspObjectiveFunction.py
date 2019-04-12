# IspObjectiveFunction.py is meant to work with the Griffon genetic algorithm,
# by providing the engine Isp as an objective function for selection
# Author: Maxime Sicat
# Date: 23/03/2019
# ISAE-SUPAERO Space Section / Griffon Project.

# --------------------------- IMPORT MODULES ---------------------------

from DataLayer.JsonInterpreter import JsonInterpreter           # Import the json interpreter
from tests.testModules import generate_data_layer
import CombustionModule.RegressionModel as Reg                  # Import the RegressionModel module
import CombustionModule.Geometries as Geom                      # Import the Geometry module
import CombustionModule.Nozzle as Noz                           # Import the Nozzle module
from CombustionModule.Combustion import CombustionObject        # Import the CombustionObject

# ------------------------ FUNCTION DEFINITIONS ------------------------



def test_combustion_image_geometry_polynom(polynom, baseRadius, branches, ox_flow):
    """ Give the geometry mean Isp as an objective function """

    # ------------ Generate Data Layer:

    json_interpreter = generate_data_layer()
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'L': 0.23,
                        'externalRadius': 88.5 / 2000,
                        'imagePixelSize': 2048,
                        'imageMeterSize': 0.09,
                        'regressionModel': Reg.MarxmanAndConstantFloodingRegimeModel(**combustion_table)}

    shape_params = {'polynom': polynom,
                    'baseRadius': baseRadius,
                    'branches': branches,
                    'n': 20}

    nozzle_params = {'At': 0.000038, 'expansion': 6.3, 'lambda_e': 0.98, 'erosion': 0}

    simulation_params = {'ox_flow': ox_flow, 'safety_thickness': 0.005}

    # ------------- Generate objects:

    geometry_obj = Geom.SinglePortImageGeometry(**geometric_params)
    geometry_obj.generatePolynom(**shape_params)
    #geometry_obj.draw_geometry()
    nozzle_obj = Noz.Nozzle(**nozzle_params)
    json_interpreter = generate_data_layer()

    # Instantiate the combustion module
    combustion_obj = CombustionObject(json_interpreter=json_interpreter,
                                      geometry_object=geometry_obj,
                                      nozzle_object=nozzle_obj)

    # -------------- Run simulation

    combustion_obj.run_simulation_constant_fuel_sliver_image_geometry(**simulation_params)

    mean_isp = 0


    for k in range(1, len(combustion_obj.results["run_values"]['isp'])):

        mean_isp += combustion_obj.results["run_values"]['isp'][k] * (combustion_obj.results["run_values"]['time'][k] - combustion_obj.results["run_values"]['time'][k-1])

    mean_isp = mean_isp / combustion_obj.results["run_values"]['time'][-1]

    return mean_isp


def test_combustion_image_geometry_fourier(a, b, baseRadius, branches, ox_flow):
    """ Give the geometry mean Isp as an objective function """

    # ------------ Generate Data Layer:

    json_interpreter = generate_data_layer()
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'L': 0.13,
                        'externalRadius': 88.5 / 2000,
                        'imagePixelSize': 2048,
                        'imageMeterSize': 0.09,
                        'regressionModel': Reg.MarxmanAndConstantFloodingRegimeModel(**combustion_table)}

    shape_params = {'a': a,
                    'b': b,
                    'impact' : 0.5,
                    'branches': branches,
                    'baseRadius': baseRadius,
                    'n': 20}

    nozzle_params = {'At': 0.000038, 'expansion': 6.3, 'lambda_e': 0.98, 'erosion': 0}

    simulation_params = {'ox_flow': ox_flow, 'safety_thickness': 0.005}

    # ------------- Generate objects:

    geometry_obj = Geom.SinglePortImageGeometry(**geometric_params)
    geometry_obj.generateFourier(**shape_params)
    #geometry_obj.draw_geometry()
    nozzle_obj = Noz.Nozzle(**nozzle_params)
    json_interpreter = generate_data_layer()

    # Instantiate the combustion module
    combustion_obj = CombustionObject(json_interpreter=json_interpreter,
                                      geometry_object=geometry_obj,
                                      nozzle_object=nozzle_obj)

    # -------------- Run simulation

    combustion_obj.run_simulation_constant_fuel_sliver_image_geometry(**simulation_params)

    mean_isp = 0

    if(len(combustion_obj.results["run_values"]['isp']) > 5):

        for k in range(1, len(combustion_obj.results["run_values"]['isp'])):

            mean_isp += combustion_obj.results["run_values"]['isp'][k] * (combustion_obj.results["run_values"]['time'][k] - combustion_obj.results["run_values"]['time'][k-1])

        mean_isp = mean_isp / combustion_obj.results["run_values"]['time'][-1]

    return mean_isp


def test_geom_validity(polynom, baseRadius, branches):

    json_interpreter = generate_data_layer()
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'L': 0.13,
                        'externalRadius': 88.5 / 2000,
                        'imagePixelSize': 2048,
                        'imageMeterSize': 0.09,
                        'regressionModel': Reg.MarxmanAndConstantFloodingRegimeModel(**combustion_table)}

    shape_params = {'polynom': polynom,
                    'baseRadius': baseRadius,
                    'branches': branches,
                    'n': 20}

    # ------------- Generate objects:

    geometry_obj = Geom.SinglePortImageGeometry(**geometric_params)
    geometry_obj.generatePolynom(**shape_params)

    return geometry_obj.min_bloc_thickness() > 0.01


def test_geom_validity_fourier(a, b, baseRadius, branches):

    json_interpreter = generate_data_layer()
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'L': 0.13,
                        'externalRadius': 88.5 / 2000,
                        'imagePixelSize': 2048,
                        'imageMeterSize': 0.09,
                        'regressionModel': Reg.MarxmanAndConstantFloodingRegimeModel(**combustion_table)}

    shape_params = {'a': a,
                    'b': b,
                    'impact': 0.5,
                    'branches': branches,
                    'baseRadius': baseRadius,
                    'n': 20}

    # ------------- Generate objects:

    geometry_obj = Geom.SinglePortImageGeometry(**geometric_params)
    geometry_obj.generatePolynom(**shape_params)

    return geometry_obj.min_bloc_thickness() > 0.01

# ---------------------------- MAIN  ---------------------------------

if __name__ == '__main__':

    print(test_combustion_image_geometry_fourier([0.6145833333333334, 0.005208333333333333, -0.3333333333333333], [-0.06770833333333333, -0.5677083333333334, 0.3697916666666667], 0.02828125, 6, 0.0921875))