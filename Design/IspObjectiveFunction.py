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

    for k in range(1, len(combustion_obj.results['isp'])):

        mean_isp += combustion_obj.results['isp'][k] * (combustion_obj.results['time'][k] - combustion_obj.results['time'][k-1])

    mean_isp = mean_isp / combustion_obj.results['time'][-1]

    return mean_isp


def test_combustion_image_geometry_fourier(a, b, baseRadius, branches, ox_flow):
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

    shape_params = {'a': polynom,
                    'b': baseRadius,
                    'impact' : 0.5,
                    'branches': branches,
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

    for k in range(1, len(combustion_obj.results['isp'])):

        mean_isp += combustion_obj.results['isp'][k] * (combustion_obj.results['time'][k] - combustion_obj.results['time'][k-1])

    mean_isp = mean_isp / combustion_obj.results['time'][-1]

    return mean_isp


def test_geom_validity(polynomn, baseRadius, branches):

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

    # ------------- Generate objects:

    geometry_obj = Geom.SinglePortImageGeometry(**geometric_params)
    geometry_obj.generatePolynom(**shape_params)

    return geometry_obj.min_bloc_thickness() > 0.01

# ---------------------------- MAIN  ---------------------------------

if __name__ == '__main__':

    test_combustion_image_geometry([0.1,0.1,0.1], 0.02, 4, 1)