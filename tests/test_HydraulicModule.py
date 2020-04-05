# test_HydraulicModule.py implements the different tests run to check on Hydraulic Module
# functionality.
# @author: Jose Felix Zapata Usandivaras
# Date: 1/6/2020
# ISAE-SUPAERO Space Section, Project "Griffon", Combustion Team

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
import cProfile                                               # Import the cProfile library
import pstats                                                 # Import pstats
from pstats import SortKey                                    # Import the SortKey
import pandas as pd                                           # Import pandas to handle time-series
import datetime                                               # Import the datetime to record the date


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

    # Generate the hydraulic module/
    my_module = HydraulicModule(hydraulic_table)

    # Set the exit chamber pressure
    my_module.set_chamber_pressure(value=3600000)

    # Initialize the pressure dofs
    my_module.initialize_pressure_dofs()

    # Initialize the mass flow
    ox_flow, pressurizer_flow = 1, 0.001
    oxidizer_tank = my_module.checkout_component("oxidizer_tank")
    oxidizer_tank.mass_node[1].dof.set_value(value=ox_flow)
    pressurizer_tank = my_module.checkout_component("pressurizer_tank")
    pressurizer_tank.mass_node[1].dof.set_value(value=pressurizer_flow)

    # Plot the figure of the pressure regulator
    # pressure_regulator = my_module.checkout_component("pressure_regulator")
    # pressure_regulator.print_curve(inlet_p=200, flow_range=(0, 1000))
    # plt.show()

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
    combustion_table = json_interpreter.return_combustion_table()

    # Generate the hydraulic module
    my_module = HydraulicModule(hydraulic_table)

    # Set the exit chamber pressure
    chamber_pressure = combustion_table['P_chamber_bar']*1e5
    my_module.set_chamber_pressure(value=chamber_pressure)

    # Initialize the pressure dofs
    my_module.initialize_pressure_dofs()

    # Initialize the mass flow
    ox_flow, pressurizer_flow = 1.11304, 0.03429
    oxidizer_tank = my_module.checkout_component("oxidizer_tank")
    oxidizer_tank.mass_node[1].dof.set_value(value=ox_flow)
    pressurizer_tank = my_module.checkout_component("pressurizer_tank")
    pressurizer_tank.mass_node[1].dof.set_value(value=pressurizer_flow)

    # Solve the module
    _, res = my_module.run_simulation()

    print(res/my_module.solver_params['scale'])

    # print the module dofs
    print(my_module.dof_print())

    # print the module nodes
    print(my_module.node_print())

    # print the module
    print(my_module)


def test_hydraulic_module_network_solution_transient():
    """ test the capability to operate in transient mode of the hydraulic module """

    # ------------ Generate the data-layer:
    json_interpreter = generate_data_layer("Griffon Data - Mock.json")
    hydraulic_table = json_interpreter.return_hydraulic_table()
    combustion_table = json_interpreter.return_combustion_table()

    # Generate the hydraulic module
    my_module = HydraulicModule(hydraulic_table)

    # Set the exit chamber pressure
    chamber_pressure = combustion_table['P_chamber_bar']*1e5
    my_module.set_chamber_pressure(value=chamber_pressure)

    # Initialize the pressure dofs
    my_module.initialize_pressure_dofs()

    # Initialize the mass flow
    ox_flow, pressurizer_flow = 1.17484, 0.01137
    oxidizer_tank = my_module.checkout_component("oxidizer_tank")
    oxidizer_tank.mass_node[1].dof.set_value(value=ox_flow)
    pressurizer_tank = my_module.checkout_component("pressurizer_tank")
    pressurizer_tank.mass_node[1].dof.set_value(value=pressurizer_flow)

    # Define delta T
    dt = 0.02
    results = {dof.number: [] for dof in my_module.dofs.elements_list}
    results['time'] = []
    time_count = 0

    # Run a transient loop
    while True:
        try:
            _, _ = my_module.run_simulation()
            for dof in my_module.dofs.elements_list:
                results[dof.number].append(dof.get_value())
            results['time'].append(time_count)
            time_count += dt
            my_module.update(dt)

        except my_module.RunEndingError:
            print("Transient Solution Finished")
            break

    # Write the results to a file
    write_results_2_file('transient_analysis_HM.txt', results,
                         [(dof.number, dof.type) for dof in my_module.dofs.elements_list])

    # Extract the oxidizer flow history
    time, ox_flow_history, ox_pressure, pressurizer_tank_p = np.array(results['time']),\
                                         np.array(results[oxidizer_tank.mass_node[1].dof.number]), \
                                         np.array(results[oxidizer_tank.pressure_nodes[1].dof.number]), \
                                         np.array(results[pressurizer_tank.pressure_nodes[1].dof.number])

    # Set the font dictionaries (for plot title and axis titles)
    title_font = {'size': '24', 'color': 'black', 'weight': 'normal',
                  'verticalalignment': 'bottom'}  # Bottom vertical alignment for more space
    axis_font = {'size': '20'}

    # plot the oxidizer flow history
    fig = plt.figure(facecolor='w', figsize=(30, 20))
    fig.suptitle('Hydraulic Module results', **title_font)
    ax = plt.axes()
    ax1 = ax.twinx()
    ax.plot(time, ox_flow_history, label='Ox Flow History', color='blue', linewidth=4.0)
    ax.set_title('')
    ax.set_ylabel('Ox Flow [kg/sec]', **axis_font)
    ax.set_xlabel('Time [sec]', **axis_font)
    ax.grid(b=True, axis='both')
    ax.set_xlim(left=time[0])
    ax.set_ylim(bottom=0.0, top=1.5)
    ax.tick_params(axis='both', labelsize=16)

    # Plot the pressure of the tanks
    ax1.plot(time, ox_pressure/1e5, color='green', label='Oxidizer Pressure', linewidth=4.0)
    ax1.plot(time, pressurizer_tank_p/1e5, color='red', label='Pressurizer Tank Pressure', linewidth=4.0)
    ax1.set_ylabel('Pressure [bar]', **axis_font)
    ax1.tick_params(axis='y', labelsize=16)

    ax.legend()
    ax1.legend()

    # Show plot
    plt.show()


def write_results_2_file(file_name, results, dofs_list):
    """
    write_results_2_file writes the results obtained from the transient simulation to a
    txt file.
    :param file_name: name of the txt file
    :param results: results dictionary
    :param dofs_list: dofs list
    :return: nothing
    """
    file_name = '../data/data_tests/' + file_name
    date_time = datetime.datetime.now().strftime("%I:%M%p on %B %d - %Y")
    header = 'Results Transitory Simulation - Date: {0} \n'.format(date_time)
    column_headers = '\n{0:<19s}\t'.format('time') + '\t'.join(('DOF_{id:>4d}_{t:<10s}'.format(id=num, t=d)
                                                                for num, d in dofs_list))
    records = pd.DataFrame(results).to_dict('records')
    # Write to the file
    with open(file_name, 'w') as f:
        f.write(header)                                     # Write the header
        f.write(column_headers)                             # Write the columns headers
        # Write the lines
        for r in records:
            line = '\n{0:>19.3f}\t'.format(r['time']) + \
                   '\t'.join(('{val:>19.3f}'.format(val=r[id_]) for id_, _ in dofs_list))
            f.write(line)


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
    # test_hydraulic_table_generation()
    # cProfile.run('test_hydraulic_table_generation')
    # p = pstats.Stats('../data/data_tests/my_code_stats')
    # p.strip_dirs().sort_stats(-1).print_stats()
    # test_hydraulic_module_initiation()
    test_hydraulic_module_network_solution()
    # test_hydraulic_module_network_solution_transient()
    # test_hydraulic_module_valve()
    # generate_injector_json_table()
