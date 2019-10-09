from abc import ABC, abstractmethod
import math


class Part(ABC):
    """Abstract base class for a generic part.

    Attributes:
        subsystem (subsystem): Parent subsystem of part.

    """

    def __init__(self, sender, **kwargs):
        """The constructor for a generic part.

        :param sender: subsystem which is generating the part.
        :param kwargs: dictionary of additional part parameters.

        """
        self.subsystem = sender

    def get_part_link(self, subsystem_ref, part_ref):
        """Get a reference to a specific part in the system.

        :param subsystem_ref: identifier of subsystem to which the part belongs.
        :param part_ref: identifier of part.
        :return: reference to part (subsystem_ref, part_ref).

        """
        return self.subsystem.parent.get_part_link(subsystem_ref, part_ref)

    @abstractmethod
    def get_mass(self):
        """Abstract method for calculating mass of part.

        This method must be implemented for each subclass of Part.

        :return: part mass

        """
        pass


class StaticPart(Part):
    """Part with static mass value

    This type of part only has a single value attributed to its mass which is unaffected by other design parameters.

    Attributes:
        mass (real): Fixed mass of part.

    """
    def __init__(self, mass, sender, **kwargs):
        """The constructor for a static part.

        :param mass: fixed mass value of part.
        :param sender: subsystem which is generating the part.
        :param kwargs: dictionary of additional part parameters.

        """
        super().__init__(sender, **kwargs)
        self.mass = mass

    def get_mass(self):
        """Method for getting mass of part.

        Overrides superclass method.

        :return: part mass.
        """

        return self.mass


class DynamicPart(Part):
    """Abstract class for a generic dynamic part.

    Subclasses of this class should define a method for determining part mass based on its parameters.

    """
    @abstractmethod
    def get_mass(self):
        """Abstract method for calculating mass of part.

        This method must be implemented for each subclass of Part.

        :return: part mass

        """
        pass


class Tank(DynamicPart):
    """Generic tank class.

    Defines a type of dynamic part whose geometry is assumed to be a cylinder with circular end caps.

    Attributes:
        radius (real): radius of the tank.
        height (real): height of the tank.

    """
    def __init__(self, tank_radius, tank_height, sender, **kwargs):
        """The constructor for a tank part.

        :param tank_radius: radius of the tank
        :param tank_height: height of the tank.
        :param sender: subsystem which is generating the part.
        :param kwargs: dictionary of additional part parameters.

        """
        super().__init__(sender, **kwargs)
        self.radius = tank_radius
        self.height = tank_height

    @abstractmethod
    def get_mass(self):
        """Abstract method for calculating mass of part.

        This method must be implemented for each subclass of Part.

        :return: part mass

        """
        pass

    def get_internal_volume(self):
        """Method for calculating material volume of tank part.

        Calculates the volume of a cylinder with tank dimensions.

        :return: internal volume of tank part.

        """
        return math.pi * self.radius**2 * self.height


class DimensionedTank(Tank):
    """Class for tanks with wall thickness affected by its parameters.

    Attributes:
        pressure (real): chamber pressure in tank.
        propellant_mass (real): mass of propellant stored in tank.
        safety_factor (real): safety factor on wall thickness.
        material (material_id): id of tank material

    """
    def __init__(self, radius, height, material_id, pressure, propellant_mass, sender, safety_factor=None, **kwargs):
        super().__init__(radius, height, sender, **kwargs)
        self.pressure = pressure
        self.propellant_mass = propellant_mass
        self.safety_factor = safety_factor
        self.material = self.subsystem.parent.materials_dict[material_id]
        self.safety_factor = self.subsystem.parent.safety_factor if not safety_factor else safety_factor

    def get_thickness(self):
        return self.pressure * self.radius / self.material.yield_stress * self.safety_factor

    def get_shell_volume(self):
        return (2 * math.pi * self.radius**2 + 2 * math.pi * self.radius * self.height) * self.get_thickness()

    def get_mass(self):
        return self.get_shell_volume() * self.material.density + self.propellant_mass
    

class OxidantTankStatic(StaticPart, Tank):
    """Class for tanks with already defined mass and dimensions. Inherits from both
       StaticPart and Tank. 
    """
    
    def __init__(self, tank_radius, tank_height, mass, propellant_mass, sender, **kwargs):
        """ class initializer """
        # Call superclass initializer
        StaticPart.__init__(mass, sender, **kwargs)
        Tank.__init__(tank_radius, tank_height, sender, **kwargs)

        # Set additional attribute
        self.propellant_mass = propellant_mass


    def get_mass(self):
        return self.get_mass() + self.propellant_mass


class CombustionChamber(DimensionedTank):

    # TODO: implement the class
    pass


class OxidantTank(DimensionedTank):

    # TODO: implement the class
    pass


class PressuriserTank(Tank):
    def __init__(self, radius, height, structure_mass, oxidiser_ref,
                 pressuriser_gas_constant, pressuriser_gamma, ambient_temperature, sender, **kwargs):
        super().__init__(radius, height, sender, **kwargs)
        self.structure_mass = structure_mass
        self.oxidiser_subsystem_ref = oxidiser_ref["subsystem_id"]
        self.oxidiser_ref = oxidiser_ref["part_id"]
        self.oxidiser = None
        self.gas_constant = pressuriser_gas_constant
        self.gamma = pressuriser_gamma
        self.temperature = ambient_temperature

    def get_mass(self):
        if not self.oxidiser:
            self.link_oxidiser()
        return self.structure_mass + self.get_pressuriser_mass()

    def get_pressuriser_mass(self):
        return self.oxidiser.pressure / (self.gas_constant * self.temperature)\
               * (1 + self.oxidiser.get_internal_volume()/self.get_internal_volume())**self.gamma\
               * self.get_internal_volume()

    def link_oxidiser(self):
        self.oxidiser = self.get_part_link(self.oxidiser_subsystem_ref, self.oxidiser_ref)


class ThermalLiner(DynamicPart):
    def __init__(self, chamber_ref, thickness, material_id, sender, **kwargs):
        super().__init__(sender, **kwargs)
        self.combustion_subsystem_ref = chamber_ref["subsystem_id"]
        self.chamber_ref = chamber_ref["part_id"]
        self.chamber = None
        self.thickness = thickness
        self.material = self.subsystem.parent.materials_dict[material_id]

    def get_mass(self):
        if not self.chamber:
            self.link_chamber()
        return self.get_volume() * self.material.density

    def get_volume(self):
        return math.pi * ((self.chamber.radius + self.thickness)**2 - self.chamber.radius**2)\
               * self.chamber.height

    def link_chamber(self):
        self.chamber = self.get_part_link(self.combustion_subsystem_ref, self.chamber_ref)


class PartsCatalogue:
    """PartsCatalogue is a static class containing a catalogue of the different parts defined
       It is used to return the desired object by passing the string name of the class
    """

    # Define a static attribute with a dictionary containing the constructor
    catalogue = {'StaticPart': StaticPart,
                 'DynamicPart': DynamicPart,
                 'Tank': Tank,
                 'DimensionedTank': DimensionedTank,
                 'CombustionChamber': CombustionChamber,
                 'OxidantTank': OxidantTank,
                 'PressuriserTank': PressuriserTank,
                 'ThermalLiner': ThermalLiner,
                 'OxidantTankStatic': OxidantTankStatic
                 }

    @staticmethod
    def search_for_part(part_name):
        return PartsCatalogue.catalogue[part_name]
