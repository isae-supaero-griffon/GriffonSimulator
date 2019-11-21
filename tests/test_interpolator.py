# Test Interpolator file generates the tests for the Interpolator Class
# Author: Jose Felix Zapata Usandivaras
# Date: 5/12/2018

# ------------------- IMPORT MODULES -----------------------

from random import uniform
from Libraries.Interpolator_old import *
from DataLayer.JsonInterpreter import *

# ------------------ FUNCTIONS DEFINITIONS -----------------


def generate_data_layer(data_file="Thermodynamic Data 36 bar OF 0,1 to 8,0.json"):
    """ generate_data_layer instantiates all of the objects that
    are required as data layer for the program to run properly.
    :param data_file: name of the json file to be used as reference for the simulation (str)
    :return JsonInterpreter """

    # Pass the file name
    data_directory = "../data"

    # Return the output
    return JsonInterpreter(file_name="/".join([data_directory, data_file]))


def save_to_json(dictionary, name):
    """
    store the result to a json file
    :param dictionary: dictionary we want to store as a json file
    :param name: name of the destination file
    :return: nothing
    """
    # Check the inputs
    assert isinstance(dictionary, dict), "Dictionary is not a dict instance"
    assert isinstance(name, str), "Insert a valid str as file_name"

    # Pass the file name
    data_directory = "../data/data_tests"

    with open("/".join([data_directory, name]), 'w') as f:
        json.dump(dictionary, f)


def test_1():
    """ run the first test for the Interpolator - with internal values """

    # Set type, data_file, and o_f_range
    type_ = "THEORETICAL ROCKET PERFORMANCE ASSUMING EQUILIBRIUM"
    data_file = generate_data_layer()
    json_interpreter = JsonInterpreter(data_file)
    interp_obj = Interpolator(simulation_type=type_, json_interpreter=json_interpreter)
    o_f_range = [interp_obj.data_dictionary[0]['state']['o/f'],
                 interp_obj.data_dictionary[-1]['state']['o/f']]
    o_f_desired_value = uniform(o_f_range[0], o_f_range[1])
    vars_ = ["t", "m", "gammas", "cstar"]
    output = interp_obj.interpolate_data(o_f_desired_value, variables=vars_, mole_fractions='')

    # Save to json
    save_to_json(output, "test_1of_{0:2.2f}.json".format(o_f_desired_value))


def test_2():
    """ run the first test for the Interpolator - with internal values """

    # Set type, data_file, and o_f_range
    type_ = "THEORETICAL ROCKET PERFORMANCE ASSUMING FROZEN COMPOSITION"
    data_file = "Thermodynamic Data.json"
    json_interpreter = JsonInterpreter(data_file)
    interp_obj = Interpolator(simulation_type=type_, json_interpreter=json_interpreter)
    o_f_range = [interp_obj.data_dictionary[0]['state']['o/f'],
                 interp_obj.data_dictionary[-1]['state']['o/f']]
    o_f_desired_value = uniform(o_f_range[0], o_f_range[1])
    output = interp_obj.interpolate_data(o_f_desired_value)

    # Save to json
    save_to_json(output, "test_2_of_{0:2.2f}.json".format(o_f_desired_value))


def test_3():
    """ run the test for the interpolator for external values (extrapolation)"""

    # Set type, data_file, and o_f_range
    type_ = "THEORETICAL ROCKET PERFORMANCE ASSUMING FROZEN COMPOSITION"
    data_file = "Thermodynamic Data.json"
    json_interpreter = JsonInterpreter(data_file)
    interp_obj = Interpolator(simulation_type=type_, json_interpreter=json_interpreter)
    o_f_range = [0.9*interp_obj.data_dictionary[0]['state']['o/f'],
                 interp_obj.data_dictionary[0]['state']['o/f']]
    o_f_desired_value = uniform(o_f_range[0], o_f_range[1])
    output = interp_obj.interpolate_data(o_f_desired_value)

    # Save to json
    save_to_json(output, "test_3_of_{0:2.2f}.json".format(o_f_desired_value))


def test_4():
    """ run the test for the interpolator for external values (extrapolation)"""

    # Set type, data_file, and o_f_range
    type_ = "THEORETICAL ROCKET PERFORMANCE ASSUMING FROZEN COMPOSITION"
    data_file = "Thermodynamic Data.json"
    json_interpreter = JsonInterpreter(data_file)
    interp_obj = Interpolator(simulation_type=type_, json_interpreter=json_interpreter)
    o_f_range = [interp_obj.data_dictionary[-1]['state']['o/f'],
                 1.1*interp_obj.data_dictionary[-1]['state']['o/f']]
    o_f_desired_value = uniform(o_f_range[0], o_f_range[1])
    output = interp_obj.interpolate_data(o_f_desired_value)

    # Save to json
    save_to_json(output, "test_4_of_{0:2.2f}.json".format(o_f_desired_value))


def test_interpreter():
    """ test the interpreter """
    data_file = "Thermodynamic Data.json"
    json_interpreter = JsonInterpreter(data_file)
    combustion_table = json_interpreter.return_combustion_table()
    print(combustion_table)


# -------------------------- MAIN -------------------------

if __name__ == '__main__':

    # Call the test functions

    test_1()
    # test_2()
    # test_3()
    # test_4()
    # test_interpreter()
