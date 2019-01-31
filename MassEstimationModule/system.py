from MassEstimationModule.subsystem import *
from MassEstimationModule.material import Material


class System:
    def __init__(self, system_dict):
        self.dict = system_dict
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
        return output_str
