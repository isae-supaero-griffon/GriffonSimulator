# Interpolator Class file for extracting Thermo-chemical data
# from the Thermodynamic Data.json file.
# Author: Jose Felix Zapata Usandivaras
# 30/11/2018

# ------------------------- IMPORT MODULES ----------------------


# ------------------------ FUNCTION DEFINITIONS -----------------


def interpolate(y2, y1, x2, x1, x):
    """
    interpolate function performs the interpolation given the values of y and x for a desired x value
    :param y2: value of y at x2
    :param y1: value of y at x1
    :param x2: right value of interpolation interval
    :param x1: left value of interpolation interval
    :param x: desired coordinate for which we perform the interpolation
    :return: interpolated value
    """
    # Check the inputs
    flag = isinstance(y2, float) and isinstance(y1, float) and isinstance(x2, float) and \
           isinstance(x1, float) and isinstance(x, float)

    # Check on the flag value, otherwise raise and exception
    if flag:
        # Perform the interpolation
        y = y1 + (y2 - y1) * (x - x1) / (x2 - x1)
    else:
        raise TypeError("Invalid type, not usefull for interpolation")

    # Return the interpolated value
    return y


# -------------------------- CLASS DEFINITIONS ------------------

class Interpolator:
    """
    Interpolator class performs the interpolation from the thermo-chemical
    data provided in Thermodynamic Data.json file generated from CEA Nasa
    runs. The json has been generated from the tables obtained in the .out file.
    To simplify data treatment all the data has been arranged into a python dictionary

    Attributes:
        data_dictionary: dictionary containing the data
    """

    def __init__(self, simulation_type, json_interpreter):
        """
        class initializer
        :param simulation_type: type of simulation (EQUILIBRIUM/FROZEN)
        :param json_interpreter: JsonInterpreter instance
        """

        # Select the simulation type and set the corresponding attribute
        self.data_dictionary = json_interpreter.return_cea_data(simulation_type)

        # Check the dictionary is not empty
        if not self.data_dictionary:
            raise ValueError('Data dictionary of Interpolator is empty, please check simulation_type')

    def interpolate_data(self, o_f_desired_value, variables='all', mole_fractions='all'):
        """
        Interpolate data performs the interpolation of the required variables given the O/F ratio provided
        :param o_f_desired_value: float value of O/F ratio desired
        :param variables: we are interested in getting out of the program
        :param mole_fractions: mole fractions specified by the user to be included in the output
        :return: dictionary with results
        """

        # Perform a search on the data dictionary and find the nearest O/F values
        o_f_value = self.data_dictionary[0]['state']['o/f']
        count = 0
        # Set a flag for the loop
        top_o_f, bottom_o_f = None, None

        # Use a while-loop to set the desired value
        while o_f_value < o_f_desired_value and count < len(self.data_dictionary) - 1:
            # Get the new of_value
            next_o_f_value = self.data_dictionary[count + 1]['state']['o/f']
            if o_f_desired_value < next_o_f_value:
                top_o_f, bottom_o_f = count + 1, count
                break
            else:
                # Perform an increment in the count
                count += 1
                bottom_o_f = count
                o_f_value = next_o_f_value

        # Perform additional checks before going with the interpolation
        if bottom_o_f is None:
            top_o_f, bottom_o_f = 1, 0
        elif top_o_f is None:
            top_o_f = len(self.data_dictionary) - 1
            bottom_o_f = top_o_f - 1

        # Call the interpolate private method to get the output
        return self._interpolate(top_o_f, bottom_o_f, o_f_desired_value, variables, mole_fractions)

    def _interpolate(self, top_o_f, bottom_o_f, o_f_desired_value, vars_list, m_fracs):
        """
        _interpolate is a private method used to perform the interpolation and get the output dictionary
        :param top_o_f: value of the top
        :param bottom_o_f: value of the bottom
        :param vars_list: variables to be included into the output
        :param m_fracs: mole fractions to be included in the output
        :return: output dictionary containing the data
        """

        # Instantiate an output dictionary
        output = {'state': {}, 'variables': {}, 'mole_fractions': {}}
        # Initialize x values
        x1 = self.data_dictionary[bottom_o_f]['state']['o/f']
        x2 = self.data_dictionary[top_o_f]['state']['o/f']

        # -------------------- STATE INTERPOLATION --------------------
        for key in self.data_dictionary[0]['state']:
            # Extract y2 and y1
            y2 = self.data_dictionary[top_o_f]['state'][key]
            y1 = self.data_dictionary[bottom_o_f]['state'][key]

            # Run interpolation in try-catch block
            try:
                # Interpolate
                y = interpolate(y2, y1, x2, x1, o_f_desired_value)

            except TypeError:
                # If type error encountered (not a float value to be interpolated) store the original output of
                # bottom_o_f
                y = y1

            # Allocate result into dictionary
            output['state'][key] = y

        # ------------------- VARIABLES INTERPOLATION -----------------

        # Check if the inputted value of vars has been set to all, in such case
        # update vars to the entire list of variables
        if vars_list == 'all':
            vars_list = self.data_dictionary[0]['variables'].keys()

        for var_name in vars_list:
            # Get the values dictionary
            values = {}
            # Go through the keys of the variable dictionary
            for key in self.data_dictionary[0]['variables'][var_name].keys():
                # Extract y2 and y1
                y2 = self.data_dictionary[top_o_f]['variables'][var_name][key]
                y1 = self.data_dictionary[bottom_o_f]['variables'][var_name][key]

                # Run interpolation in a try-catch block
                try:
                    # Interpolate
                    y = interpolate(y2, y1, x2, x1, o_f_desired_value)

                except TypeError:
                    # If type error encountered (not a float value to be interpolated) store the original output of
                    # bottom_o_f
                    y = y1

                # Store output
                values[key] = y

            # Allocate result of interpolation into output dictionary
            output['variables'][var_name] = values

        # ------------------- SPECIES INTERPOLATION --------------------
        for species_name, species_data in self.data_dictionary[bottom_o_f]['mole_fractions'].items():
            # First check if we have a dictionary or not for the species data
            values = {}

            # Check if the mole fraction requested is given
            if m_fracs == 'all' or species_name in m_fracs:

                # Check the species name is included into the list of species
                if isinstance(species_data, dict):
                    for key, y1 in species_data.items():
                        # Extract other data point
                        y2 = self.data_dictionary[top_o_f]['mole_fractions'][species_name][key]

                        # Run interpolation in a try-catch block
                        try:
                            # Interpolate
                            y = interpolate(y2, y1, x2, x1, o_f_desired_value)

                        except TypeError:
                            # If type error encountered (not a float value to be interpolated) store the original
                            # output of bottom_o_f

                            y = y1

                        # Store the output
                        values[key] = y

                        # Allocate result of interpolation into output dictionary
                        output['mole_fractions'][species_name] = values

                else:
                    # Extract the values for interpolation
                    y1 = species_data
                    y2 = self.data_dictionary[top_o_f]['mole_fractions'][species_name]

                    # Run interpolation in a try-catch block
                    try:
                        # Interpolate
                        y = interpolate(y2, y1, x2, x1, o_f_desired_value)

                    except TypeError:
                        # If type error encountered (not a float value to be interpolated) store the original output of
                        # bottom_o_f
                        y = y1

                    # Allocate the result of interpolation into output dictionary
                    output['mole_fractions'][species_name] = y

        # Return the output dictionary
        return output
