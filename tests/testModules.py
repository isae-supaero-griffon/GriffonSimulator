# Test the different modules separately with a basic test case as part of the integration of the modules.
# Author/s: Jose Felix Zapata Usandivaras
# Date: 30/12/2018

# --------------------------- IMPORT MODULES ---------------------------

from DataLayer.JsonInterpreter import JsonInterpreter           # Import the json interpreter
import CombustionModule.RegressionModel as Reg                  # Import the RegressionModel module
import CombustionModule.Geometries as Geom                      # Import the Geometry module
import CombustionModule.Nozzle as Noz                           # Import the Nozzle module
import CombustionModule.Combustion as Comb                      # Import the Combustion module
from CombustionModule.Fuel import *                             # Import the Fuel class
from MassEstimationModule.system import System                  # Import the system class
from TrajectoryModule.Drag import *                             # Import the Drag library
from TrajectoryModule.Density import DensityLaw                 # Import the density-law library
from TrajectoryModule.Trajectory import TrajectoryObject        # Import the trajectory object
import numpy as np                                              # Import numpy
import matplotlib.pyplot as plt                                 # Import matplotlib

# ------------------------ FUNCTION DEFINITIONS ------------------------


def generate_data_layer(data_file="DataLayerONERATests_Hycom16.json"):
    """ generate_data_layer instantiates all of the objects that
    are required as data layer for the program to run properly.
    :param data_file: name of the json file to be used as reference for the simulation (str)
    :return JsonInterpreter """

    # Pass the file name
    data_directory = "../data"

    # Return the output
    return JsonInterpreter(file_name="/".join([data_directory, data_file]))


def test_combustion():
    """ perform the test over the combustion module """

    # ------------ Generate Data Layer:

    json_interpreter = generate_data_layer()
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'L': 0.4,
                        'rintInitial': 0.05,
                        'rext0': 0.07,
                        'regressionModel': Reg.MarxmanAndConstantFloodingRegimeModel(**combustion_table)}

    nozzle_params = {'At': 0.000589, 'expansion': 5.7, 'lambda_e': 0.98, 'erosion': 0}

    design_params = {'gamma': 1.27, 'p_chamber': 3200000, 'p_exit': 100000,
                     'c_star': 1500, 'ox_flow': 1.2, 'OF': 5}

    simulation_params = {'ox_flow': 1.3, 'safety_thickness': 0.005, 'dt': 0.05}

    # ------------- Generate objects:

    geometry_obj = Geom.OneCircularPort(**geometric_params)
    nozzle_obj = Noz.Nozzle(**nozzle_params)
    nozzle_obj.set_design(**design_params)
    json_interpreter = generate_data_layer()

    # Instantiate the combustion module
    combustion_obj = CombustionObjectClassic(json_interpreter=json_interpreter,
                                      geometry_object=geometry_obj,
                                      nozzle_object=nozzle_obj)

    # -------------- Run simulation & Plot:

    combustion_obj.run_simulation_constant_fuel_sliver(**simulation_params)

    # Print the module
    print(combustion_obj)

    # Plot the results
    combustion_obj.plot_results()


def test_combustion_image_geometry():
    """ perform the test over the combustion module """

    # ------------ Generate Data Layer:

    json_interpreter = generate_data_layer()
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'L': 0.4,
                        'externalRadius': 0.45,
                        'imagePixelSize': 2048,
                        'imageMeterSize': 0.1,
                        'regressionModel': Reg.MarxmanAndConstantFloodingRegimeModel(**combustion_table)}

    shape_params = {'a' : [0.07552083333333333, -0.22916666666666666, 0.3151041666666667],
                    'b' : [],
                    'baseRadius' : 0.036,
                    'branches' : 8,
                    'impact' : 0.3,
                    'n' : 40}

    nozzle_params = {'At': 0.000589, 'expansion': 5.7, 'lambda_e': 0.98, 'erosion': 0}

    design_params = {'gamma': 1.27, 'p_chamber': 3200000, 'p_exit': 100000,
                     'c_star': 1500, 'ox_flow': 1., 'OF': 9}

    simulation_params = {'ox_flow': 0.96875, 'safety_thickness': 0.000000001, 'max_burn_time' : 5, 'dt': 0.05}

    # ------------- Generate objects:

    geometry_obj = Geom.SinglePortImageGeometry(**geometric_params)
    geometry_obj.generateFourier(**shape_params)
    #geometry_obj.draw_geometry()
    nozzle_obj = Noz.Nozzle(**nozzle_params)
    nozzle_obj.set_design(**design_params)
    json_interpreter = generate_data_layer()

    # Instantiate the combustion module
    combustion_obj = Comb.CombustionObjectImage(json_interpreter=json_interpreter,
                                      geometry_object=geometry_obj,
                                      nozzle_object=nozzle_obj)

    # -------------- Run simulation & Plot:

    combustion_obj.run_simulation_constant_fuel_sliver(**simulation_params)

    # Print the module
    print(combustion_obj)

    # Plot the results
    combustion_obj.plot_results()


def test_combustion_three_port_geometry():
    """ perform the test over the combustion module """

    # ------------ Define parameters:

    geometric_params = {'L': 0.6, 'portsIntialRadius': 0.01, 'r_ext': 0.07}

    nozzle_params = {'At': 0.000589, 'expansion': 5.7, 'lambda_e': 0.98, 'erosion': 0}

    design_params = {'gamma': 1.27, 'p_chamber': 4000000, 'p_exit': 100000,
                     'c_star': 1500, 'isp': 230, 'thrust': 30000}

    simulation_params = {'ox_flow': 1.0, 'safety_thickness': 0.005, 'dt': 0.05}

    # ------------- Generate objects:

    geometry_obj = Geom.ThreeCircularPorts(**geometric_params)
    nozzle_obj = Noz.Nozzle(**nozzle_params)
    nozzle_obj.set_design(**design_params)
    json_interpreter = generate_data_layer()

    # Instantiate the combustion module
    combustion_obj = CombustionObjectClassic(json_interpreter=json_interpreter,
                                      geometry_object=geometry_obj,
                                      nozzle_object=nozzle_obj)

    # -------------- Run simulation & Plot:

    combustion_obj.run_simulation_constant_fuel_sliver(**simulation_params)

    # Print the module
    print(combustion_obj)

    # Plot the results
    combustion_obj.plot_results()


def test_combustion_onera_data():
    """ perform the test over the combustion module but with Onera Test Data """

    # ------------ Generate the data layer:

    json_interpreter = generate_data_layer(data_file="Thermodynamic Data Onera 41 bar H2O2 87_5.json")
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'L': 0.235,
                        'rintInitial': 0.0186658954,
                        'rext0': 0.2,
                        'regressionModel': Reg.MarxmanAndConstantFloodingRegimeModel(**combustion_table)}

    nozzle_params = {'At': 0.000038, 'expansion': 6.3, 'lambda_e': 0.98, 'erosion': 0}

    simulation_params = {'ox_flow': 0.0855, 'safety_thickness': 0.005, 'dt': 0.05, 'max_burn_time': 8.25}

    # ------------- Generate objects:

    geometry_obj = Geom.OneCircularPort(**geometric_params)
    nozzle_obj = Noz.Nozzle(**nozzle_params)
    json_interpreter = generate_data_layer(data_file="Thermodynamic Data Onera 41 bar H2O2 87_5.json")

    # Instantiate the combustion module
    combustion_obj = CombustionObjectClassic(json_interpreter=json_interpreter,
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
        2. Go range: [100, 500] kg/m^2/sec
        3. Chamber pressure: 36 bar, closest possible to the target pressure of Griffon.
    """

    # ------------ Generate the data layer:

    json_interpreter = generate_data_layer()                    # Use same data-layer used for Griffon (same pressure)
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Define parameters:

    geometric_params = {'L': 0.157,
                        'rintInitial': 0.007465,
                        'rext0': 0.041,
                        'regressionModel': Reg.TwoRegimesMarxmanAndFloodedModel(**combustion_table)}

    nozzle_params = {'At': 0.000038, 'expansion': 6.3, 'lambda_e': 0.98, 'erosion': 0}

    simulation_params = {'ox_flow': 0.0876, 'safety_thickness': 0.005, 'dt': 0.01, 'max_burn_time': 8}

    # ------------- Generate objects:

    geometry_obj = Geom.OneCircularPort(**geometric_params)
    nozzle_obj = Noz.Nozzle(**nozzle_params)
    json_interpreter = generate_data_layer(data_file="Thermodynamic Data Onera 41 bar H2O2 87_5.json")

    # Instantiate the combustion module
    combustion_obj = CombustionObjectClassic(json_interpreter=json_interpreter,
                                      geometry_object=geometry_obj,
                                      nozzle_object=nozzle_obj)

    # -------------- Run simulation & Plot:

    combustion_obj.run_simulation_constant_fuel_sliver(**simulation_params)

    # Print the module
    print(combustion_obj)

    # Plot the results
    combustion_obj.plot_results()


def test_mass_simulator():
    """ test the mass simulator module with a simple test-case """

    # ----------- Generate objects:

    json_interpreter = generate_data_layer()
    mass_simulator_obj = System(json_interpreter.return_mass_simulator_table())

    # Print the results
    print(mass_simulator_obj)


def test_trajectory():
    """ perform the test over the trajectory module """

    # -------------- Define parameters:

    # thrust & time
    delta_t = 0.01                                          # delta-time in seconds
    simulation_time = 60                                    # simulation-time
    n_points = int(simulation_time / delta_t) + 1           # total number of points
    thrust = []                                             # Initiate the thrust array
    isp = []                                                # Initiate the isp array
    time = []                                               # Initiate the time array
    burn_time = 4.5
    # Burn time in seconds
    constant_thrust = 2600                         # Thrust value in newtons
    constant_isp = 210                              # Isp value in seconds

    for i in range(0, n_points):
        t = delta_t*i
        time.append(t)
        if t < burn_time:
            thrust.append(constant_thrust)
            isp.append(constant_isp)
        else:
            thrust.append(0)
            isp.append(np.nan)

    # isp, area_ref, initial conditions
    initial_conditions = {'h0': 0, 'v0': 0, 'm0': 35}       # Initial mass in kg

    # ------------- Generate objects:

    json_interpreter = generate_data_layer()
    trajectory_data = json_interpreter.return_trajectory_table()
    density_obj = DensityLaw(trajectory_data['density'])

    drag_parameters = {'drag_coefficient': trajectory_data['drag']['Cd'],
                       'area_ref': trajectory_data['drag']['area_ref'],
                       'density': density_obj}

    drag_obj = DragCatalogue.return_drag_object(trajectory_data['drag']['type'], drag_parameters)
    trajectory_obj = TrajectoryObject(density_obj=density_obj, drag_obj=drag_obj)

    # -------------- Run simulation:

    trajectory_obj.run_simulation_on_trajectory(time=np.asarray(time),
                                                thrust=np.asarray(thrust),
                                                isp=np.asarray(isp),
                                                initial_conditions=initial_conditions)

    # Print the results
    print(trajectory_obj)

    # Plot the results
    trajectory_obj.plot_results()


def test_fuel():
    """ Test the Fuel Module first definition """

    # ------------ Generate Data Layer:

    json_interpreter = generate_data_layer()

    # ------------ Instantiate Fuel

    my_fuel = Fuel(json_interpreter)

    # ------------ Test the object

    # Obtain the output
    params = {'Pc': bar2psia(36), 'eps': 36, 'MR': 5.5, 'short_output': 1}
    s = my_fuel.ispObj.get_full_cea_output(**params)

    # Print the result
    print(s)


# ---------------------------- MAIN  ---------------------------------

if __name__ == '__main__':

    # Call on test_combustion method
    # test_combustion()
    # test_combustion_three_port_geometry()
    # test_mass_simulator()
    #test_trajectory()
    #test_combustion_onera_data()
    # test_onera_physical_test_2()
    test_combustion_image_geometry()
    # test_fuel()
    # Show any plots
    plt.show()
