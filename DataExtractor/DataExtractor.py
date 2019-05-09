# Griffon Project: Program to extract the Chemical Data from the Analysis file
# Author: Jose Felix Zapata Usandivaras
# Date: 29/11/2018


# ------------------------- IMPORT MODULES --------------------------

import json                                     # Import json module
from MassEstimationModule.parameters import *   # Import the system_dict
from math import pi                             # Import pi from math library

# ------------------------- FUNCTION DEFINITIONS --------------------


def return_float(value):
    """is_float checks if the string inputted is a float or not"""
    try:
        new_value = float(value)
        return new_value
    except ValueError:
        return None


def extract_data_from_line(zones, zones_data_line, record_dict):
    """ extract the data out of the line of file and store it in records_dict """

    # Instantiate a new variable_dict
    variable_dict = {}

    # Get the data, split the line and obtain the variable name, units and associated numbers
    zone_data_variable = zones_data_line[:17]
    zones_data_numbers = [zones_data_line[17:25], zones_data_line[26:34], zones_data_line[35:43]]
    zones_data_numbers = [return_float(number) for number in zones_data_numbers]
    zones_data_variable = zone_data_variable.replace(' ', '')
    zones_data_variable = zones_data_variable.split(',')

    # Set the variable name and the variable unit (if there is)
    variable_name = zones_data_variable[0]
    variable_unit = zones_data_variable[1] if len(zones_data_variable) > 1 else ''

    # TODO: Rethink how to deal with RHO values line
    if len(zones_data_numbers) == len(zones):
        variable_dict['name'] = variable_name.lower()
        variable_dict['unit'] = variable_unit.lower()
        variable_dict.update({zone_name: zone_data for zone_name, zone_data
                              in zip(zones, zones_data_numbers)})

        # Append the zones_data_line to the variable_dict
        record_dict['variables'][variable_dict['name']] = variable_dict


def extract_mole_fractions_equilibrium(zones, zones_data_line, record_dict):
    """ extract the data out of the line of file and store it in records_dict """

    # Instantiate a new variable_dict
    variable_dict = {}

    # Get the data, split the line and obtain the variable name, units and associated numbers
    zone_data_variable = zones_data_line[:17]
    zones_data_numbers = [zones_data_line[17:25], zones_data_line[26:34], zones_data_line[35:43]]
    zones_data_numbers = [return_float(number) for number in zones_data_numbers]
    variable_name = zone_data_variable.replace(' ', '')

    if len(zones_data_numbers) == len(zones):
        variable_dict.update({zone_name: zone_data for zone_name, zone_data
                              in zip(zones, zones_data_numbers)})

        # Append the zones_data_line to the variable_dict
        record_dict['mole_fractions'][variable_name] = variable_dict


def extract_mole_fractions_frozen(data_line, record_dict):
    """ extract the data out of the line of file and store it in records_dict for frozen composition """

    # Get the data, split the line and obtain the composition name and associated numbers
    data_line = data_line.split(' ')
    data_line = [value for value in data_line if value != '']

    # Populate the dictionary
    for i in range(0, len(data_line)-1, 2):
        record_dict['mole_fractions'][data_line[i]] = return_float(data_line[i+1])


def generate_cea_table(cea_txt_file):
    """
    generate_cea_table generates the dictionary which contains the CEA_file
    :param cea_txt_file: .txt file name with CEA data
    :return: dict object with the requested table
    """

    # Read the file and extract the data from the lines
    with open(cea_txt_file, 'r') as f:
        text = f.readlines()

    values_dict = {}
    file_length = len(text)
    count1 = 0
    while count1 < file_length:
        line = text[count1].replace('\n', '')

        # Extract the simulation type
        if 'THEORETICAL' in line:
            simulation_type = line.strip()
            # Add the simulation type if not existent in the values_dict
            if simulation_type not in values_dict:
                values_dict[simulation_type] = []

        if 'Pin = ' in line:
            data_dict = {'Pin': return_float(line[9:15])/14.5038}

        record_dict = {'state': None, 'variables': {}, 'mole_fractions': {}}
        # Start to parse the table
        if 'O/F' in line:

            # Extract first line of the table
            data = line.replace("=", " ")
            data = data.split(" ")
            data = [value for value in data if value != '']
            for j in range(0, len(data)-1, 2):
                data_dict[data[j].lower()] = float(data[j+1])

            # Set the state of the record_dict
            record_dict['state'] = data_dict
            #  Get the chamber/throat/exit data
            count2 = 2
            zones = text[count1+count2].replace('\n', '')
            zones = zones.split(' ')
            zones = [zone for zone in zones if zone != '']

            flag_performance = True
            while flag_performance:
                count2 += 1                                                 # Add 1 to the second count
                zones_data_line = text[count1+count2].replace('\n', '')     # Get the line of interest

                if zones_data_line != '':

                    # Check if we've reached the performance table
                    if 'PERFORMANCE' not in zones_data_line:
                        # Extract the data out of the line
                        extract_data_from_line(zones, zones_data_line, record_dict)
                    else:
                        flag_performance = False                                # Set the flag to False to stop the loop

            flag_fractions = True
            while flag_fractions:
                count2 += 1                                                 # Add 1 to the second count
                zones_data_line = text[count1+count2].replace('\n', '')     # Get the line of interest

                if zones_data_line != '':

                    # Check if we've reached the performance table
                    if 'MOLE FRACTIONS' not in zones_data_line:
                        # Extract the data out of the line
                        extract_data_from_line(zones, zones_data_line, record_dict)
                    else:
                        flag_fractions = False                                # Set the flag to False to stop the loop

            flag_thermal = True
            while flag_thermal:
                count2 += 1
                zones_data_line = text[count1+count2].replace('\n', '')         # Get the line of interest

                if zones_data_line != '':

                    # Check if we've reached the performance table
                    if 'THERMODYNAMIC' not in zones_data_line:
                        # Extract the data out of the line
                        if 'EQUILIBRIUM' in simulation_type:
                            extract_mole_fractions_equilibrium(zones, zones_data_line, record_dict)
                        elif 'FROZEN' in simulation_type:
                            extract_mole_fractions_frozen(zones_data_line, record_dict)

                    else:
                        flag_thermal = False                                 # Set the flag to False to stop the loop

            # Append de record_dict to the values_dict
            values_dict[simulation_type].append(record_dict)

        else:

            # Set the count2 if no line of 'O\F' has been found
            count2 = 1

        # Add count2 to count1
        count1 += count2

    # Return the dictionary of values
    return values_dict


def generate_combustion_table():
    """
    generate_combustion_table generates the dictionary containing all of the combustion's module
    static_data
    :return: dict with combustion table properties
    """

    # Introduce the dictionary by hand and output it later
    table = {'g0': 9.81,
             'R': 8.314,
             'Pa': 1.013*10**5,
             'a': 3.27*10**-5,
             'n': 0.84,
             'm': -0.2,
             'rho_fuel': 1055,
             'P_chamber_bar': 36,
             'rho_ox_pure': 1.45,
             'ox_purity': 0.875
             }

    # Return the dict
    return table


def generate_mass_simulator_table():
    """
    generate_mass_simulator_table generates the dictionary containing all of the mass simulator module
    static_data
    :return: dict with mass simulator table properties
    """
    # Return the dict
    return system_dict


def generate_trajectory_table():
    """
    generate_trajectory_table generates the dictionary containing all of the trajectory module
    static_data
    :return: dict with trajectory table properties
    """

    # Introduce the dictionary by hand and output it later (all in SI Units)
    table = {
              'density':    {
                             'g0': 9.81,                    # Gravity acceleration used for density calculation
                             'R': 8.314,                    # Gas Constant
                             'M': 0.02896,                  # Air molar mass
                             'rho_0': 1.225,                # Density of air
                             'T_0': 298,                    # Air Reference Temperature
                             'h_0': 0                       # Departure height
                            },
              'drag':       {
                             'type': 'constant',
                             'Cd': 0.4,
                             'area_ref': pi * (0.2 * 0.2) / 4
                            }
            }                    # Reference drag coefficient

    # Return the dict
    return table


def generate_propellants_table():
    """
    generate_propellants_table generates the table associated to the propellant that will be input
    into the RocketCEA wrapper
    :return: dict with propellants table
    """

    # Introduce the dictionary by hand and output it later (all in SI Units)
    table = {
                # Define the oxidizer/s table
                'oxidizer': {
                                'name': "H2O2GriffonOxidizer",
                                'components':   [{
                                                    'oxid': "H2O2(L)",                                # Name of the oxidizer
                                                    'wt': 87.5,                                       # Mass weight percentage of the Oxidizer
                                                    't': {'value': 298, 'units': '(k)'}               # Reference Temperature in Kelvin
                                                  },
                                                {
                                                    'oxid': "H2O(L)",                                 # Name of the oxidizer
                                                    'wt': 12.5,                                       # Mass weight percentage of the Oxidiser
                                                    't': {'value': 298, 'units': '(k)'},              # Reference Temperature in Kelvin
                                                }]},

                # Define the fuel/s table
                'fuel': {
                            'name': "3DPrintedABSGriffon",
                            'components':

                                                [{
                                                     'fuel': "Air",  # Name of the fuel
                                                     'wt': 0.01,  # Mass weight percentage of the Fuel
                                                     't': {'value': 298, 'units': '(k)'}
                                                 # Reference Temperature in Kelvin
                                                 },

                                                 {
                                                     'fuel': "ABS",  # Name of the fuel
                                                     'wt': 99.9,  # Mass weight percentage of the Fuel
                                                     't': {'value': 298, 'units': '(k)'},
                                                 # Reference Temperature in Kelvin
                                                     'h': {'value': 62.63, 'units': 'kj/mol'},
                                                 # Reference Enthalpy in kj/mol
                                                     'composition': {'C': 3.85, 'H': 4.85, 'N': 0.43}
                                                 # Composition of the fuel
                                                 }]}
    }

    # Return the output
    return table


def generate_file_from_scratch():
    """ Generate the json file from scratch """

    # Generate the combustion_table
    combustion_table = generate_combustion_table()

    # Generate the mass_simulator_values
    mass_simulator_table = generate_mass_simulator_table()

    # Generate the trajectory_table
    trajectory_table = generate_trajectory_table()

    # Generate the propellant_table
    propellant_table = generate_propellants_table()

    # Get the data repository
    data_directory = "../data"
    # Generate the cea_table
    cea_file_name = "Analysis Data 36 bar OF 0,1 to 8,0 H2O2 87,5.txt"
    CEA_Data = generate_cea_table(cea_txt_file="/".join([data_directory, cea_file_name]))

    # instantiate the global dictionary
    values_dictionary = {'combustion_table': combustion_table,
                         'mass_simulator_table': mass_simulator_table,
                         'trajectory_table': trajectory_table,
                         'propellant_table': propellant_table,
                         'CEA_Data': CEA_Data}

    # Store the generated dictionary to a json file
    output_file = "Thermodynamic Data 36 bar OF 0,1 to 8,0 H2O2 87,5 V2.json"
    with open("/".join([data_directory, output_file]), 'w') as f:
        json.dump(values_dictionary, f)


def append_table_to_file():
    """ Append the propellants table to the already existing file"""

    # Generate the table to append
    new_table = generate_propellants_table()

    # File name & Data Directory
    data_directory = "../data"
    file_name = "Thermodynamic Data 36 bar OF 0,1 to 8,0 H2O2 87,5.json"

    # Read the dictionary from the json file
    with open("/".join([data_directory, file_name]), 'r') as f:
        # Read the file
        json_str = f.read()
        json_data = json.loads(json_str)

    # Add the new table to the extracted dictionary
    json_data['propellant_table'] = new_table

    # Overwrite the file
    with open("/".join([data_directory, file_name]), 'w') as f:
        # Dump the data into the json file
        json.dump(json_data, f)


# ---------------------------- MAIN ---------------------------------

if __name__ == '__main__':

    # Call the corresponding method
    append_table_to_file()
