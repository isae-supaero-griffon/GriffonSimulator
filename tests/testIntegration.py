# testIntegration.py script is in charge of testing the integration of the different modules in
# order to allow the interaction between them.
# Author: Jose Felix Zapata Usandivaras
# Date: 2/1/2018.

# ---------------------- IMPORT MODULES -----------------------

from Initializer.Initializer import *                               # Import the Initializer object
from CombustionModule.Combustion import *                           # Import the Combustion module
from IntegrationModule.SimulationObject import SimulationObject     # Import the SimulationObject


# -------------------- FUNCTIONS DEFINITIONS ------------------


def generate_fourier_coefficients(n_coefs, period, fun, *args):
    """ generate_fourier_coefficients for given configuration """

    # Generate the x coordinates
    x = np.linspace(0, 1, 1000)

    # Preallocate results as numpy arrays
    a_s = np.empty(shape=(n_coefs,1))
    b_s = np.empty(shape=(n_coefs,1))

    for i in range(1, n_coefs+1):
        sin_fourier = lambda x_: m.sin(2 * m.pi * i * x_ / period) * fun(x_, *args)
        cos_fourier = lambda x_: m.cos(2 * m.pi * i * x_ / period) * fun(x_, *args)
        sin_fourier = np.vectorize(sin_fourier)
        cos_fourier = np.vectorize(cos_fourier)
        a_s[i-1] = 2 / period * trapz(cos_fourier(x), x)
        b_s[i-1] = 2 / period * trapz(sin_fourier(x), x)

    # Return the results
    return a_s, b_s


def my_fun_2(x, *args):
    delta1, delta2 = args[0], args[1]
    if 0<=x<delta1:
        return 0
    elif delta1<=x<delta2:
        return 0.5
    elif delta2<=x<1-delta2:
        return 1 + 0.12*m.fabs(x - 0.5)
    elif 1-delta2<=x<1-delta1:
        return 0.5
    elif 1-delta1<=x<=1:
        return 0


def my_fun(x, *args):
    """
    my_fun generates the desired profile for the branch
    :param x: float with x-coordinate
    :param br: float with branch radius position
    :param r: radius
    :return: value of function
    """
    br, r, delta = args[0], args[1], args[2]

    # Check the inputs
    assert r < 0.5, 'The radius has to be less than 0.5'
    assert 0<= x <=1, 'x has to be contained between 0 and 1'
    assert r < br, 'br has to be greater than r'

    if delta <= x <= 1-delta:
        # Scale x
        x = (x - delta)/(1 - 2*delta)
        # Evaluate the function
        if 0<=x <br-r:
            return (1-r)/(br-r)*x
        elif br-r<=x<br:
            return m.sqrt(r**2 - (x - br)**2) + (1 - r)
        elif br<=x<1-br:
            return 1
        elif 1-br<=x<1-br+r:
            return sqrt(r**2 - (x - 1 + br)**2) + (1 - r)
        elif 1-br+r<=x<=1:
            return (1-r)/(br-r)*(1-x)
    else:
        return 0


def generate_data_layer(data_file):
    """ generate_data_layer instantiates all of the objects that
    are required as data layer for the program to run properly.
    :param data_file: name of the json file to be used as reference for the simulation (str)
    :return JsonInterpreter """

    # Pass the file name
    data_directory = "../data"

    # Return the output
    return JsonInterpreter(file_name="/".join([data_directory, data_file]))


def test_simulation_initializer():
    """ the method is in charge of generating the inputs for the initializer and
    instantiating it. """

    # ------------ Generate the data-layer:

    json_interpreter = generate_data_layer("DataLayerONERATests_Hycom14.json")

    # ---------- Pack the inputs:

    init_parameters = {
                        'combustion': {
                                       'geometric_params': {'type': SinglePortImageGeometry, 'length': 0.4,
                                                            'r_ext': 0.05, 'image_pixel_size': 1024,
                                                            'image_meter_size': 0.1},

                                       'shape_params': {'a': [0, 1], 'b': [1, 0, 0],
                                                            'base_radius': 0.036, 'branches': 10, 'impact': 0.1,
                                                            'n': 40},

                                       'nozzle_params': {'At': 0.000589, 'expansion': 5.7, 'lambda_e': 0.98,
                                                         'erosion': 0},

                                       'set_nozzle_design': True,

                                       'design_params': {'gamma': 1.27, 'p_chamber': 3200000, 'p_exit': 100000,
                                                         'c_star': 1500, 'ox_flow': 1.2, 'OF': 5},
                                      },


                      }
    simulation_parameters = {
                              'combustion': {'ox_flow': 1, 'safety_thickness': 0.005, 'dt': 0.05,
                                             'max_burn_time': None
                                             },

                              'mass_simulator': {'ox_flow': 1, 'burn_time': 'TBD', 'extra_filling': 0.05,
                                                 'injection_loss': 0.5, 'area_injection': 0.000105,
                                                 'system': SystemDynamic
                                                 },

                              'trajectory': {'initial_conditions': {'h0': 0, 'v0': 0, 'm0': 'TBD'},
                                             'simulation_time': 100
                                             }
                            }

    # -------------- Generate the initializer:

    init_obj = Initializer(init_parameters=init_parameters,
                           simulation_parameters=simulation_parameters,
                           json_interpreter=json_interpreter)

    # -------------- Generate the simulation object:

    simulation_object = SimulationObject(initializer_collection=init_obj)

    # --------------- Run the simulation:

    simulation_object.run_simulation_in_batch()

    # Print the total mass
    print("\nRockets Total Mass: {0} kgs".format(simulation_object.mass_simulator_module.get_mass()))

    # Print the splitted masses
    print(simulation_object.mass_simulator_module)

    # --------------- Export results to csv files:

    # data directory
    data_directory = "../data/data_tests"

    file_name_expression = "Griffon Output Test {number}.csv"

    simulation_object.export_results_to_file(file_name_expression="/".join([data_directory,
                                                                            file_name_expression]))

    # --------------- Plot the results

    simulation_object.results_collection.elements_list[0].combustion.plot_results()
    simulation_object.results_collection.elements_list[0].trajectory.plot_results()


def test_hydraulic_module_integration_with_combustion():
    """ The aim of this function is to integrate and test such integration between the
    Hydraulic module and the Combustion module """

    # ------------ Generate the data-layer:
    json_interpreter = generate_data_layer("Griffon Data - Mock.json")
    combustion_table = json_interpreter.return_combustion_table()

    # ------------ Generate the Fourier Coefficients:
    r, br, delta, n_coefs, period = 0.2, 0.3, 0.17, 100, 1
    a_s, b_s = generate_fourier_coefficients(n_coefs, period, my_fun, br, r, delta)

    # ---------- Pack the inputs:
    init_parameters = {
        'combustion': {
            'geometric_params': {'length': 0.40,
                                 'regression_model': TwoRegimesMarxmanAndFloodedModel(**combustion_table),
                                 'r_ext': 0.05,
                                 'image_pixel_size': 2048*4,
                                 'image_meter_size': 0.1},

            'shape_params': {'a': a_s, 'b': b_s,
                             'base_radius': 0.032, 'branches': 12, 'impact': 0.8,
                             'n': 50},

            'nozzle_params': {'At': 0.000589, 'expansion': 5.7, 'lambda_e': 0.98,
                              'erosion': 0},

            'set_nozzle_design': False
        },
        'hydraulic': {
                      'chamber_pressure': combustion_table['P_chamber_bar']*1e5,
                      'ox_flow': 1.08,
                      'pressurizer_flow': 0.01137
        }
    }
    simulation_parameters = {'hydraulic_module': None, 'safety_thickness': 0.004, 'dt': 0.05,
                             'max_burn_time': None, 'tol_press': 1e-3}

    # ---------- Generate objects:

    geometry_obj = SinglePortImageGeometry(**init_parameters['combustion']['geometric_params'])
    geometry_obj.generate_fourier(**init_parameters['combustion']['shape_params'])
    nozzle_obj = Nozzle(**init_parameters['combustion']['nozzle_params'])

    # Generate the Combustion Module:
    combustion_obj = CombustionObjectImage(json_interpreter=json_interpreter,
                                           geometry_object=geometry_obj,
                                           nozzle_object=nozzle_obj)

    # Generate the Hydraulic Module:
    hydraulic_obj = HydraulicModule(json_interpreter.return_hydraulic_table())
    hydraulic_obj.initialize(init_parameters['hydraulic'])
    simulation_parameters['hydraulic_module'] = hydraulic_obj

    # Generate the plot before it runs
    geometry_obj.draw_geometry()

    # ---------- Run simulation and & Plot:

    combustion_obj.run_full_flow_analysis(**simulation_parameters)

    # Print the module
    print(combustion_obj)

    # Plot the results
    combustion_obj.plot_results()
    geometry_obj.draw_geometry()
    plt.show()


def test_hydraulic_module_integration_with_combustion_1d():

    # ------------ Generate the data-layer:
    json_interpreter = generate_data_layer("Griffon Data - Mock.json")
    combustion_table = json_interpreter.return_combustion_table()

    # ---------- Pack the inputs:

    # Define the oxidizer flow
    init_parameters = {
        'combustion': {
            'geometric_params': {
                                 'length': 0.4,
                                 'regression_model': TwoRegimesMarxmanAndFloodedModel(**combustion_table),
                                 'r_init': 0.052 / 2,
                                 'r_final': 0.094 / 2,
                                 'r_exit': 0.08 / 2,
                                 'exit_depth': 0.02,
                                 'depth': 0.3,
                                 'r_ext': 0.05,
                                 'nodes_number': 200
                                },

            'nozzle_params': {'At': 0.000589, 'expansion': 5.7, 'lambda_e': 0.98,
                              'erosion': 0},

            'set_nozzle_design': False
        },
        'hydraulic': {
            'chamber_pressure': combustion_table['P_chamber_bar']*1e5,
            'ox_flow': 1.08,
            'pressurizer_flow': 0.01137
        }
    }

    simulation_parameters = {'hydraulic_module': None, 'safety_thickness': 0.004, 'dt': 0.05,
                             'max_burn_time': None, 'tol_press': 1e-3}

    # ---------- Generate objects:

    geometry_obj = ConicCircularPort1D(**init_parameters['combustion']['geometric_params'])
    nozzle_obj = Nozzle(**init_parameters['combustion']['nozzle_params'])

    # Generate the Combustion Module:
    combustion_obj = CombustionObject1D(json_interpreter=json_interpreter,
                                        geometry_object=geometry_obj,
                                        nozzle_object=nozzle_obj)

    # Generate the Hydraulic Module:
    hydraulic_obj = HydraulicModule(json_interpreter.return_hydraulic_table())
    hydraulic_obj.initialize(init_parameters['hydraulic'])
    simulation_parameters['hydraulic_module'] = hydraulic_obj

    # ---------- Run simulation and & Plot:

    combustion_obj.run_full_flow_analysis(**simulation_parameters)

    # Print the module
    print(combustion_obj)

    # Plot the results
    combustion_obj.plot_results()
    plt.show()

# ----------------------------- MAIN ------------------------------


if __name__ == '__main__':

    # Call the test functions
    # test_simulation_initializer()
    # test_hydraulic_module_integration_with_combustion()
    test_hydraulic_module_integration_with_combustion_1d()
