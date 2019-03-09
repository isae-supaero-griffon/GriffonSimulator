from MassEstimationModule.subsystem import *
from MassEstimationModule.material import Material
from abc import ABC, abstractmethod


# Define the abstract class for the system
class System(ABC):
    """
    System abstract class has by objective to define a class for the mass simulator module with
    the most important methods to be implemented.
    """
    def __init__(self, system_dict):
        """ class initializer """

        # Define a total mass property if requested from the system
        self.dict = system_dict

    @abstractmethod
    def get_mass(self):
        """ get mass is an abstract method """


class SystemStatic(System):
    """
    SystemStatic is a class which inherits from System and implements a system with
    a given dry mass
    """
    def __init__(self, system_dict):
        """
        class initializer
        :param system_dict: dictionary containing the detailed description of the system
        The dictionary must contain:
            0. system dry mass
            1. mass of propellant required (fuel)
            2. mass of oxidizer
            3. mass of pressurizer fluid
        """
        # Call superclass constructor
        super().__init__(system_dict)

    def get_mass(self):
        """ calculate the wet initial mass of the system
        :return initial wet mass"""
        return sum((value for value in self.dict.values()))


class SystemDynamic(System):
    def __init__(self, system_dict):
        super().__init__(system_dict)
        self.materials_dict = {}
        self.subsystems_dict = {}
        self.safety_factor = self.dict["safety_factor"]
        for mat in self.dict["materials"].values():
            new_material = Material(**mat)
            self.add_material(new_material)
        for subsystem_values in self.dict["subsystems"].values():
            new_subsystem = Subsystem(subsystem_values, self)
            self.add_subsystem(new_subsystem)

    def add_material(self, material):
        self.materials_dict[material.name] = material

    def add_subsystem(self, subsystem):
        self.subsystems_dict[subsystem.name] = subsystem

    def get_mass(self):
        return sum(ss.get_mass() for ss in self.subsystems_dict.values())

    def get_part_link(self, part_subsystem, part_name):
        """ 
        :param part_subsystem: name of the subsystem of the desired part
        :param part_name: name of the desired part
        :return: part object with given name
        """
        subsystem = self.subsystems_dict[part_subsystem]
        part_link = subsystem.get_part(part_name)
        return part_link

    def __str__(self):

        output_str = " \n".join(("\n {0} \n\t {1}".format(key, value) for key, value in
                                 sorted(self.subsystems_dict.items())))

        total_mass = " \n\n Rockets Total Mass, {0} kgs".format(self.get_mass())

        return output_str + total_mass

