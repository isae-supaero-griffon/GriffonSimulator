# Drag.py contains the declaration of the Drag class alongside all of the other
# classes associated to it.
# Author: Jose Felix Zapata Usandivaras
# Date: 30/12/2018

# ---------------------- IMPORT MODULES -----------------------

from abc import ABC, abstractmethod     # Import ABC for abstract methods definition
from TrajectoryModule.Density import *  # Import the density module

# --------------------- CLASS DEFINITIONS ---------------------


class Drag(ABC):
    """
    Drag class performs the calculation of the Drag force present on the rocket.
    It is an abstract class with no attributes except for the methods. The default drag
    calculation is by means of a Drag coefficient.
    """

    def __init__(self, constants, density):
        """
        class initializer
        :param constants: dictionary containing the relevant constants for the calculation of drag
        """

        # Call superclass initializer
        super().__init__()

        # Set the attributes
        self.constants = constants
        self.density = density

    def __str__(self):
        """ redefine the str method for the Drag class to print the constants """
        return "\n".join(["\nDrag,",
                          "\n".join(("{0}, {1:3.4f}".format(key, value) for key, value in self.constants.items()))])

    @abstractmethod
    def compute_drag_coefficient(self, *args, **kwargs):
        """ compute_drag_coefficient computes de drag coefficient
        By default it takes a variable number of inputs. To account for different cases """

        # Return 0 to comply with code-checking
        return 0

    @abstractmethod
    def compute_area_reference(self, *args, **kwargs):
        """ compute_area_reference computes the reference area
        By default it takes a variable number of inputs to account for the different cases. """

        # Return 0 to comply with code-checking
        return 0

    def compute_drag_force(self, speed, altitude):
        """
        compute_drag_force computes the drag force present based on the
        :param speed: float- speed [m/sec]
        :param altitude: float-altitude [m]
        :return: float-drag force [nodes_number]
        """

        # Compute the drag coefficient
        rho_ = self.density.compute_density(altitude)
        cd = self.compute_drag_coefficient()
        s = self.compute_area_reference()

        # Return overall drag force
        return -0.5*rho_*speed*abs(speed)*cd*s


class ConstantDrag(Drag):
    """ ConstantDrag is used to calculate Drag force when a constant
        drag-coefficient is used. """

    def __init__(self, drag_coefficient, area_ref, density):
        """
        class initializer
        """
        super().__init__({'cd': drag_coefficient, 'S': area_ref}, density)

    def compute_drag_coefficient(self):
        """ compute the drag coefficient (return the constant) """
        return self.constants['cd']

    def compute_area_reference(self):
        """ compute the reference area (return the constant) """
        return self.constants['S']


class DragCatalogue:
    """ DragCatalogue is a static class which helps on choosing the type of drag to
    use for the calculations """

    drag_collection = {'constant': ConstantDrag}

    @staticmethod
    def return_drag_object(obj_type, input_parameters):
        """
        return the object initializer
        :param obj_type: str containing the objects name
        :param input_parameters: input parameters for the object to be created
        :return: object initializer
        """
        return DragCatalogue.drag_collection[obj_type](**input_parameters)
