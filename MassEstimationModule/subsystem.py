from MassEstimationModule.part import *


class Subsystem:
    def __init__(self, subsystem_dict, sender):
        self.parent = sender
        self.dict = subsystem_dict
        self.name = self.dict["identifier"]
        self.parts_dict = {}
        for part_dict in self.dict["parts"].values():
            self.add_part(part_dict)

    def add_part(self, part_dict):
        part_name = part_dict['identifier']
        part_type = PartsCatalogue.search_for_part(part_name=part_dict['type'])
        part = part_type(sender=self, **part_dict)
        self.parts_dict[part_name] = part

    def get_mass(self):
        return sum(p.get_mass() for p in self.parts_dict.values())

    @staticmethod
    def return_part_property(self, part_name, property_name):
        """ return a part property"""

    def get_part(self, part_name):
        part = self.parts_dict[part_name]
        return part

    def __str__(self):
        output_str = "".join(("\n\t {0}, {1:5.3} kg".format(key, value.get_mass()) for key, value in
                                 sorted(self.parts_dict.items())))
        return output_str
