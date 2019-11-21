# testRocketCEA.py file has as purpose to test the behavior of the RocketCEA python library available through pip.
# For documentation on the library's resources refer to: https://rocketcea.readthedocs.io/en/latest/
# @author: Jose Felix Zapata Usandivaras
# Date: 26/04/2019
# ISAE - SUPAERO SCUBE Project Griffon
#

# --------------------------- IMPORT MODULES -------------------------

import rocketcea.cea_obj as cea_obj  # Import the RocketCEA module
import matplotlib.pyplot as plt      # Import Matplotlib
import numpy as np                   # Import Numpy

# ------------------------ FUNCTION DEFINITIONS -----------------------


def bar2psia(a):
    """
    bar2psia converts from bar to psia
    :param a: float with pressure in bar
    :return: pressure in psia
    """
    return (a + 1.013) * 14.50377


def define_propellants():
    """ define the string that characterizes the fuel """

    # Introduce the fuel definition
    fuel_card = "fuel=Air  wt=0.01  t(k)=298 " \
                "fuel=ABS  wt=99.9  t(k)=298 " \
                "h,kj/mol=62.63  C 3.85 H 4.85 N 0.43"

    # Introduce the oxidizer definition
    oxidizer_card = "oxid=H2O2(L) wt=87.5  t(k)=298 " \
                    "oxid=H2O(L) wt=12.5  t(k)=298"

    # Add the new fuel and oxidizer
    cea_obj.add_new_fuel("3DPrinted_ABS", fuel_card)
    cea_obj.add_new_oxidizer("GriffonOxydizer_H2O2", oxidizer_card)


def run_cea():
    """ Run CEA returns runs the code of CEA for the inputted variables """

    # Define the parameters of the run
    params = {'Pc': bar2psia(36), 'eps': 5.5072, 'MR': 8.3, 'short_output': 1}
    propellant = {'oxName': 'GriffonOxydizer_H2O2', 'fuelName': '3DPrinted_ABS'}

    # Instantiate the CEA object
    ispObj = cea_obj.CEA_Obj(**propellant)

    # Obtain the output
    s = ispObj.get_full_cea_output(**params)

    # Print the result
    my_file = 'CEANASA_Results_Equilibrium.txt'
    with open(my_file, 'w') as f:
        f.write(s)

    print(s)


# --------------------------- MAIN -----------------------------

if __name__ == '__main__':

    # Define propellants and then run cea
    define_propellants()
    run_cea()
