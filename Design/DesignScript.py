# DesignScript.py is a script intended to perform iterations on the design of the Griffon Combustion Chamber based on
# the expressed requirements.
# Author: Jose Felix Zapata Usandivaras
# Date: 6/02/2019
# ISAE-SUPAERO Space Section / Griffon Project.

# ---------------------- IMPORT MODULES -----------------------

from Initializer.Initializer import *                               # Import the Initializer object
from CombustionModule.Combustion import *                           # Import the Combustion Models
from IntegrationModule.SimulationObject import SimulationObject     # Import the SimulationObject
import numpy as np                                                  # Import numpy
import itertools                                                    # Import itertools
import CombustionModule.RegressionModel as Reg                      # Import the RegressionModel class
import json                                                         # Import the json library

# -------------------- FUNCTIONS DEFINITIONS ------------------


def generate_data_layer(data_file="Griffon Data - ABS - H2O2 - 36 bar.json"):
    """ generate_data_layer instantiates all of the objects that
    are required as data layer for the program to run properly.
    :param data_file: name of the json file to be used as reference for the simulation (str)
    :return JsonInterpreter """

    # Pass the file name
    data_directory = "../data"

    # Return the output
    return JsonInterpreter(file_name="/".join([data_directory, data_file]))


def single_case_analysis_one_circular_port():
    """
    Study the behavior of the rocket for a single case with OneCircularPort
    :return: nothing
    """

    # ------------ Generate the data-layer:

    json_interpreter = generate_data_layer("Griffon Data - ABS - H2O2 - 36 bar - Mock.json")

    # ---------- Pack the inputs:

    ox_flow = 1.09

    init_parameters = {
                        'combustion': {
                                       'geometric_params': {'type': OneCircularPort,
                                                            'L': 0.325,
                                                            'rintInitial': 0.03,
                                                            'rext0': 0.05,
                                                            'regressionModel': Reg.TwoRegimesMarxmanAndFloodedModel},

                                       'nozzle_params': {'At': 0.000545, 'expansion': 4.96559, 'lambda_e': 0.98,
                                                         'erosion': 0},

                                       'set_nozzle_design': False,

                                       'design_params': {'gamma': 1.27, 'p_chamber': 3600000, 'p_exit': 100000,
                                                         'c_star': 1430, 'ox_flow': ox_flow, 'OF': 12.5},
                                      },


                      }
    simulation_parameters = {
                              'CombustionModel': CombustionObjectClassic,

                              'combustion': {
                                             'ox_flow': ox_flow,
                                             'safety_thickness': 0.005,
                                             'dt': 0.01,
                                             'max_burn_time': 4.5},

                              'mass_simulator': {'ox_flow': ox_flow, 'burn_time': 'TBD',
                                                 'system': SystemDynamic, 'extra_filling': 0.1,
                                                 'injection_loss': 0.3, 'area_injection': 0.000105},

                              'trajectory': {'initial_conditions': {'h0': 0, 'v0': 0, 'm0': 'TBD'},
                                             'simulation_time': 60}
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

    # Print combustion results
    print(simulation_object.combustion_module)

    # --------------- Plot the results

    simulation_object.results_collection.elements_list[0].combustion.plot_results()
    simulation_object.results_collection.elements_list[0].trajectory.plot_results()

    # # data directory
    # data_directory = "../data/data_tests"
    #
    # file_name_expression = "Tentative Design Single Port {number}.csv"
    #
    # simulation_object.export_results_to_file(file_name_expression="/".join([data_directory,
    #                                                                         file_name_expression]))
    #
    # # Save to json the mass simulator dict
    # output_file = "Tentative Design Single Port 1.json"
    # with open("/".join([data_directory, output_file]), 'w') as f:
    #     json.dump(simulation_object.mass_simulator_module.dict, f)

    # Show any plots
    plt.show()


def single_case_analysis_one_port_review():
    """
    Study the behavior of the rocket for a single case with OneCircularPort
    :return: nothing
    """

    # ------------ Generate the data-layer:

    json_interpreter = generate_data_layer("Thermodynamic Data 36 bar OF 0,1 to 8,0 H2O2 87,5.json")

    # ---------- Pack the inputs:

    ox_flow = 1.09

    init_parameters = {
                        'combustion': {
                                       'geometric_params': {'type': OneCircularPort,
                                                            'L': 0.325,
                                                            'rintInitial': 0.03,
                                                            'rext0': 0.05,
                                                            'regressionModel': Reg.MarxmanAndConstantFloodingRegimeModel},

                                       'nozzle_params': {'At': 0.000589, 'expansion': 5.7, 'lambda_e': 0.98,
                                                         'erosion': 0},

                                       'set_nozzle_design': True,

                                       'design_params': {'gamma': 1.27, 'p_chamber': 3600000, 'p_exit': 100000,
                                                         'c_star': 1580, 'ox_flow': ox_flow, 'OF': 5.5},
                                      },


                      }
    simulation_parameters = {
                              'combustion': {'ox_flow': ox_flow, 'safety_thickness': 0.005, 'dt': 0.01,
                                             'max_burn_time': 4.75},

                              'mass_simulator': {'ox_flow': ox_flow, 'burn_time': 'TBD',
                                                 'system': SystemStatic,
                                                 'system_dict': {'dry_mass': 12.35 + 20 + 1.3,
                                                                 'pressurizer_mass': 0.17,
                                                                 'oxidizer_mass': 'TBD',
                                                                 'fuel_mass': 'TBD'},
                                                 'extra_filling': 0.1,
                                                 'injection_loss': 0.5,
                                                 'area_injection': 0.000105},

                              'trajectory': {'initial_conditions': {'h0': 0, 'v0': 0, 'm0': 'TBD'},
                                             'simulation_time': 60}
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

    # Print combustion results
    print(simulation_object.combustion_module)

    # Print Trajectory results
    print(simulation_object.trajectory_module)

    # --------------- Plot the results

    simulation_object.results_collection.elements_list[0].combustion.plot_results()
    simulation_object.results_collection.elements_list[0].trajectory.plot_results()

    # # data directory
    # data_directory = "../Design/Design Files/Proposed Designs/Single Port Geometry/"
    #
    # file_name_expression = "Tentative Design Single Port Reviewed {number}.csv"
    # 
    # simulation_object.export_results_to_file(file_name_expression="/".join([data_directory,
    #                                                                         file_name_expression]))
    #
    # # Save to json the mass simulator dict
    # output_file = "Tentative Design Single Port Reviewed 1.json"
    # with open("/".join([data_directory, output_file]), 'w') as f:
    #     json.dump(simulation_object.mass_simulator_module.dict, f)

    # Show any plots
    plt.show()


def single_case_analysis_one_port_image_geometry():
    """ case study for image geometry design """

    # ------------ Generate the data-layer:

    json_interpreter = generate_data_layer()

    # ---------- Pack the inputs:

    init_parameters = {
                        'combustion': {
                                       'geometric_params': {'type': SinglePortImageGeometry, 'L': 0.4,
                                                            'externalRadius': 0.05, 'imagePixelSize': 1024,
                                                            'imageMeterSize': 0.1},

                                       'shape_params': {'a': [0,1], 'b': [1, 0, 0],
                                                            'baseRadius': 0.036, 'branches': 10, 'impact': 0.1,
                                                            'n': 40},

                                       'nozzle_params': {'At': 0.000589, 'expansion': 5.7, 'lambda_e': 0.98,
                                                         'erosion': 0},

                                       'set_nozzle_design': True,

                                       'design_params': {'gamma': 1.27, 'p_chamber': 3200000, 'p_exit': 100000,
                                                         'c_star': 1500, 'ox_flow': 1.2, 'OF': 5},
                                      },


                      }
    simulation_parameters = {
                              'CombustionModel': CombustionObjectImage,

                              'combustion': {'ox_flow': 1, 'safety_thickness': 0.005, 'dt': 0.05,
                                             'max_burn_time': None},

                              'mass_simulator': {'ox_flow': 1, 'burn_time': 'TBD', 'extra_filling': 0.05,
                                                 'injection_loss': 0.5, 'area_injection': 0.000105, 'system' : SystemDynamic},

                              'trajectory': {'initial_conditions': {'h0': 0, 'v0': 0, 'm0': 'TBD'},
                                             'simulation_time': 100}
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
    #
    # # data directory
    # data_directory = "../data/data_tests"
    #
    # file_name_expression = "Griffon Output Test {number}.csv"
    #
    # simulation_object.export_results_to_file(file_name_expression="/".join([data_directory,
    #                                                                         file_name_expression]))

    # --------------- Plot the results

    simulation_object.results_collection.elements_list[0].combustion.plot_results()
    simulation_object.results_collection.elements_list[0].trajectory.plot_results()


def single_case_analysis_three_circular_ports():
    """
    Study the behavior of the rocket for a single case with ThreeCircularPort
    :return: nothing
    """

    # ------------ Generate the data-layer:

    json_interpreter = generate_data_layer("Thermodynamic Data 36 bar OF 0,1 to 8,0 H2O2 87,5.json")

    # ---------- Pack the inputs:

    # Set the ox_flow
    ox_flow = 1.12

    init_parameters = {
                        'combustion': {
                                       'geometric_params': {'type': ThreeCircularPorts, 'L': 0.325,
                                                                'portsIntialRadius': 0.016,
                                                                'r_ext': 0.07,
                                                                'regressionModel': Reg.TwoRegimesMarxmanAndFloodedModel},

                                       'nozzle_params': {'At': 0.000589, 'expansion': 5.7, 'lambda_e': 0.98,
                                                         'erosion': 0},

                                       'set_nozzle_design': True,

                                       'design_params': {'gamma': 1.27, 'p_chamber': 3600000, 'p_exit': 100000,
                                                         'c_star': 1500, 'ox_flow': ox_flow, 'OF': 5},
                                      },


                      }
    simulation_parameters = {
                              'combustion': {'ox_flow': ox_flow, 'safety_thickness': 0.005, 'dt': 0.01,
                                             'max_burn_time': 5},

                              'mass_simulator': {'ox_flow': ox_flow, 'burn_time': 'TBD',
                                                 'system': SystemDynamic,'extra_filling': 0.1,
                                                 'injection_loss': 0.5, 'area_injection': 0.000105},

                              'trajectory': {'initial_conditions': {'h0': 0, 'v0': 0, 'm0': 'TBD'},
                                             'simulation_time': 60}
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

    # # data directory
    # data_directory = "../data/data_tests"
    #
    # file_name_expression = "Tentative Design {number}.csv"
    #
    # simulation_object.export_results_to_file(file_name_expression="/".join([data_directory,
    #                                                                         file_name_expression]))
    #
    # # Save to json the mass simulator dict
    # output_file = "Tentative Design 1.json"
    # with open("/".join([data_directory, output_file]), 'w') as f:
    #     json.dump(simulation_object.mass_simulator_module.dict, f)

    # --------------- Plot the results

    simulation_object.results_collection.elements_list[0].combustion.plot_results()
    simulation_object.results_collection.elements_list[0].trajectory.plot_results()

    # Show any plots
    plt.show()


def generate_analysis_cases_three_port_geometry():
    """ generate analysis cases generates the associated dictionaries that will be used in the
        run in batch method.
        :return InitializerCollection"""

    # ------------ Generate the data-layer:

    json_interpreter = generate_data_layer("Thermodynamic Data 36 bar OF 0,1 to 8,0 H2O2 87,5.json")

    # Set the parameters over which we want to iterate
    burn_time = None
    external_radius = 0.01*np.linspace(5, 8, 10)
    chamber_length = 0.4
    internal_radius = 0.015
    ox_flows = 0.001*np.linspace(50, 1000, 20)

    # Get the Initializer Collection
    collection = InitializerCollection([])
    combinations = list(itertools.product(external_radius, ox_flows))

    for r_ext, ox_flow in combinations:

        # ---------- Pack the inputs:

        init_parameters = {
                            'combustion': {
                                           'geometric_params': {'type': ThreeCircularPorts, 'L': chamber_length,
                                                                'portsIntialRadius': internal_radius,
                                                                'r_ext': r_ext,
                                                                'regressionModel': Reg.MarxmanAndConstantFloodingRegimeModel},

                                           'nozzle_params': {'At': 0.000589,
                                                             'expansion': 5.7,
                                                             'lambda_e': 0.98,
                                                             'erosion': 0},

                                           'set_nozzle_design': True,

                                           'design_params': {'gamma': 1.27,
                                                             'p_chamber': 3600000,
                                                             'p_exit': 100000,
                                                             'c_star': 1500,
                                                             'ox_flow': ox_flow,
                                                             'OF': 5},
                                          },


                            }

        simulation_parameters = {
                                  'combustion': {'ox_flow': ox_flow, 'safety_thickness': 0.0025, 'dt': 0.05,
                                                 'max_burn_time': burn_time},

                                  'mass_simulator': {'ox_flow': ox_flow, 'burn_time': 'TBD',
                                                     'system': SystemDynamic, 'extra_filling': 0.1,
                                                     'injection_loss': 0.5, 'area_injection': 0.000105},

                                  'trajectory': {'initial_conditions': {'h0': 0, 'v0': 0, 'm0': 'TBD'},
                                                 'simulation_time': 60}
                                }

        # Add initializer to the collection
        new = Initializer(init_parameters=init_parameters,
                          simulation_parameters=simulation_parameters,
                          json_interpreter=json_interpreter)
        collection.add_element(new)

    # Return the collection
    return collection


def generate_analysis_cases_port_geometry():
    """ generate analysis cases generates the associated dictionaries that will be used in the
        run in batch method.
        :return InitializerCollection"""

    # ------------ Generate the data-layer:

    json_interpreter = generate_data_layer("Thermodynamic Data 36 bar OF 0,1 to 8,0 H2O2 87,5.json")

    # Set the parameters over which we want to iterate
    chamber_length = 0.3
    external_radius = 0.05
    internal_radius = 0.01*np.linspace(1, 3, 10)
    ox_flows = 0.001*np.linspace(50, 200, 10)

    # Get the Initializer Collection
    collection = InitializerCollection([])
    combinations = list(itertools.product(internal_radius, ox_flows))

    for r_int, ox_flow in combinations:

        # ---------- Pack the inputs:

        init_parameters = {
                            'combustion': {
                                           'geometric_params': {'type': OneCircularPort, 'L': chamber_length,
                                                                'rintInitial': r_int,
                                                                'rext0': external_radius,
                                                                'regressionModel': Reg.MarxmanAndConstantFloodingRegimeModel},

                                           'nozzle_params': {'At': 0.000589,
                                                             'expansion': 5.7,
                                                             'lambda_e': 0.98,
                                                             'erosion': 0},

                                           'set_nozzle_design': True,

                                           'design_params': {'gamma': 1.27,
                                                             'p_chamber': 3600000,
                                                             'p_exit': 100000,
                                                             'c_star': 1500,
                                                             'ox_flow': ox_flow,
                                                             'OF': 5},
                                          },


                            }

        simulation_parameters = {
                                  'combustion': {'ox_flow': ox_flow, 'safety_thickness': 0.005, 'dt': 0.05,
                                                 'max_burn_time': 8},

                                  'mass_simulator': {'ox_flow': ox_flow, 'burn_time': 'TBD',
                                                     'system': SystemDynamic, 'extra_filling': 0.05,
                                                     'injection_loss': 0.5, 'area_injection': 0.000105},

                                  'trajectory': {'initial_conditions': {'h0': 0, 'v0': 0, 'm0': 'TBD'},
                                                 'simulation_time': 60}
                                }

        # Add initializer to the collection
        new = Initializer(init_parameters=init_parameters,
                          simulation_parameters=simulation_parameters,
                          json_interpreter=json_interpreter)
        collection.add_element(new)

    # Return the collection
    return collection


def generate_analysis_cases_multi_port_geometry():
    """ generate analysis cases generates the associated dictionaries that will be used in the
        run in batch method.
        :return InitializerCollection"""

    # ------------ Generate the data-layer:

    json_interpreter = generate_data_layer("Thermodynamic Data 36 bar OF 0,1 to 8,0 H2O2 87,5.json")

    # Set the parameters over which we want to iterate
    ports_number = 4
    burn_time = 5
    external_radius = 0.08
    chamber_length = 0.4
    internal_radius = 0.01*np.linspace(0.5, 2, 10)
    r_center = 0.01
    ox_flows = 0.001*np.linspace(50, 400, 10)

    # Get the Initializer Collection
    collection = InitializerCollection([])
    combinations = list(itertools.product(internal_radius, ox_flows))

    for r_int, ox_flow in combinations:

        # ---------- Pack the inputs:

        init_parameters = {
                            'combustion': {
                                           'geometric_params': {'type': MultipleCircularPortsWithCircularCenter,
                                                                'L': chamber_length,
                                                                'N': ports_number,
                                                                'ringPortsIntialRadius': r_int,
                                                                'centralPortInitialRadius': r_center,
                                                                'r_ext': external_radius,
                                                                'regressionModel': Reg.MarxmanAndConstantFloodingRegimeModel},

                                           'nozzle_params': {'At': 0.000589,
                                                             'expansion': 5.7,
                                                             'lambda_e': 0.98,
                                                             'erosion': 0},

                                           'set_nozzle_design': True,

                                           'design_params': {'gamma': 1.27,
                                                             'p_chamber': 3600000,
                                                             'p_exit': 100000,
                                                             'c_star': 1500,
                                                             'ox_flow': ox_flow,
                                                             'OF': 5},
                                          },


                            }

        simulation_parameters = {
                                  'combustion': {'ox_flow': ox_flow, 'safety_thickness': 0.0025, 'dt': 0.05,
                                                 'max_burn_time': burn_time},

                                  'mass_simulator': {'ox_flow': ox_flow, 'burn_time': 'TBD',
                                                     'system': SystemDynamic, 'extra_filling': 0.05,
                                                     'injection_loss': 0.5, 'area_injection': 0.000105},

                                  'trajectory': {'initial_conditions': {'h0': 0, 'v0': 0, 'm0': 'TBD'},
                                                 'simulation_time': 60}
                                }

        # Add initializer to the collection
        new = Initializer(init_parameters=init_parameters,
                          simulation_parameters=simulation_parameters,
                          json_interpreter=json_interpreter)
        collection.add_element(new)

    # Return the collection
    return collection


def run_design_cases():
    """
    run_design_cases performs the run in batch of the design cases defined in previous function
    """

    # ----------------- Get the Initializer Collection:
    cases = generate_analysis_cases_three_port_geometry()   # Change the function for different types of geometries

    # -------------- Generate the simulation object:

    simulation_object = SimulationObject(initializer_collection=cases)

    # --------------- Run the simulation:

    simulation_object.run_simulation_in_batch()

    # --------------- Export results to csv files:

    # data directory
    data_directory = "../Design/Design Files/ThreePortGeometry_OXFlow"

    file_name_expression = "Griffon Output Multi-Port Geometry {number}.csv"

    simulation_object.export_results_to_file(file_name_expression="/".join([data_directory,
                                                                            file_name_expression]))


# ----------------------------- MAIN ------------------------------


if __name__ == '__main__':

    # Execute the program
    # run_design_cases()
    # single_case_analysis_three_circular_ports()
    # single_case_analysis_one_port_review()
    # single_case_analysis_one_port_image_geometry()
    single_case_analysis_one_circular_port()
