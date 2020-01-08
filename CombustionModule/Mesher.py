# The Mesher.py file includes the class Mesh, used to define the mesh. The
# mesh is important to define the domain of discretization and the specific
# locations at which the problem is to be solved. By default the mesh strategy
# is to define an uniformly spaced mesh. Further developments can enhance this
# capability from a programming standpoint, however it is to be noticed that
# the Hydres model is derived for a uniformly spaced mesh.
# In addition the class cell is defined in this document.
# @author: Jose Felix Zapata Usandivaras
# @date: 19/08/2019
# ISAE - SUPAERO, MAE1 Research Project
# Reference: Messineo, Jérôme. Modélisation des instabilités hydrodynamiques dans les moteurs-fusées hybrides.
# Diss. Toulouse, ISAE, 2016.

# ------------------------- IMPORT MODULES --------------------------- #

from abc import ABC, abstractmethod                     # Import abstract class module
from Libraries.Interpolator import Interpolator         # Import Interpolator
import CombustionModule.Geometries as Geom              # Import Geometries object
import math as m                                        # Import math module
import numpy as np                                      # Import numpy



# ----------------------- FUNCTION DEFINITIONS ----------------------- #




# ------------------------- CLASS DEFINITIONS ------------------------ #

class Mesh(ABC):
    """
    The class Mesh contains the cells assembly where it is to be calculated
    the magnitudes of the problem.
    Its attributes are:
        0. name: string containing the name of the mesh.
        1. cells: array of cells
    """

    def __init__(self, name, geometry_obj):
        """ class initializer
            :param  name: string containing the name of the mesh.
            :param geometry_obj: Geometry class instance which defined the domain of the problem.
        """
        # Call superclass initializer
        super(Mesh, self).__init__()

        # Check the inputs
        assert isinstance(name, str), "Failed assertion: name of mesh has to be of str type.\n"
        assert isinstance(geometry_obj, Geom.Geometry1D), "Failed assertion: geometry_obj has to be of Geometry type.\n"

        # Initialize attributes
        self.name = name
        self.geometry = geometry_obj
        # TODO: find a more convenient array type with iterator protocol to operate with cells
        self.cells = []

    @abstractmethod
    def _generate_cells(self, cell_factory):
        pass

    def return_x_cor(self):
        """
        return_x_corr returns a numpy array with the x_cor of the mesh
        :return: np.array of x_cor
        """
        return np.array([my_cell.x_cor for my_cell in self.cells])

    def return_profile_data(self):
        """
        return_profile_data returns a numpy array with the profile characteristic of the mesh
        for a circular port its the radius. For others is to be defined
        :return: np.array of profile data
        """
        return np.array([my_cell.return_profile_data() for my_cell in self.cells])

    def return_data(self):
        """
        return_data provides the data to run the interpolation
        :return: x, area array, perimeter array
        """
        x, areas, perimeters = zip(*[(my_cell.x_cor, my_cell.return_area_data(), my_cell.return_perimeter_data()) for
                                    my_cell in self.cells])
        return np.array(x), np.array(areas), np.array(perimeters)



class UniformlySpacedMesh(Mesh):
    """
    The UniformlySpacedMesh inherits from the Mesh class.
    """

    def __init__(self, name, geometry_obj, n_el):
        """
        class initializer
        :param name: string defining the name of the Mesh object
        :param geometry_obj: Geometry class instance which defined the domain of the problem.
        :param n_el: number of cells integer.
        """
        # Call superclass constructor
        super(UniformlySpacedMesh, self).__init__(name, geometry_obj)

        # Check the input
        assert isinstance(n_el, int), "Failed assertion: check n_el is an int.\n"

        # Attributes definition
        self.n_el = n_el
        self.dx = geometry_obj.L / n_el
        self._generate_cells(geometry_obj.my_cell_factory)

    def _generate_cells(self, cell_factory):
        """
        Generate the cells according to the law proposed.
        :param cell_factory: cell factory method associated to the geometry object
        :return: nothing
        """

        for i in range(0, self.n_el):
            new_cell = cell_factory(i, (i+1)*self.dx)
            self.cells.append(new_cell)


class GeometricMesh(Mesh):
    """
    The Geometric Mesh inherits from the mesh class and it implements a geometric
    series to define the mesh.
    """

    def __init__(self, name, geometry_obj, n_el, bias):
        """
        class initializer
        :param name: string defining the name of the Mesh object
        :param geometry_obj: Geometry class instance which defined the domain of the problem.
        :param n_el: number of cells integer.
        :param bias: bias to apply on the geometric definition of the cells
        """

        # Call superclass constructor
        super(GeometricMesh, self).__init__(name, geometry_obj)

        # Check the input
        assert isinstance(n_el, int), "Failed assertion: check n_el is an int.\n"
        assert 0 < bias, "Bias has to be strictly greater than 0. \n"

        # Attributes definition
        self.n_el = n_el
        self.bias = bias
        self._generate_cells(geometry_obj.my_cell_factory)

    def _generate_cells(self, cell_factory):
        """
        Generate the cells according to the law proposed.
        :param cell_factory: cell factory method associated to the geometry object
        :return: nothing
        """

        # Calculate the initial value displacement
        initial_dx = self.geometry.L * (1 - self.bias) / (1 - self.bias ** self.n_el)
        accumulated_distance = float(0)

        for i in range(0, self.n_el+1):
            new_cell = cell_factory(i, accumulated_distance)
            self.cells.append(new_cell)
            accumulated_distance += initial_dx * (self.bias ** i)


class Cell(ABC):
    """
    The Cell class is in charge of storing the local magnitudes associated to a given
    position within the discretized domain of the rocket.
    Its attributes are:
        0. number: integer starting from 0 defining the position on the mesh array
        1. x_cor: float defining the x coordinate of the cell. First cell starts at 0 and last
                  finishes at L equivalent to the length of the domain.
        2. profile: associated profile to the cell depends on the type of geometry
        3. r_ext: external radius of the cell
    """

    def __init__(self, num, x, profile, r_ext):
        """
        class initializer
        :param num: integer indicating the cell number.
        :param x: float which indicates the cell position.
        :param profile: profile which characterizes the cell
        :param r_ext: external radius of the cell
        """
        # Call superclass initializers
        super(Cell, self).__init__()

        # Initialize attributes
        assert isinstance(num, int), "Failed assertion: cell number is not an int instance.\n"
        assert isinstance(x, float), "Failed assertion: cell position is not a float instance.\n"

        # Set the attributes
        self.number = num
        self.x_cor = x
        self.profile = profile
        self.r_ext = r_ext
        self.min_thickness = None

    def set_profile(self, my_profile):
        """
        set_profile sets the profile property
        :param my_profile: depends on the type geometry but it characterizes the geometry
        description
        :return: nothing
        """
        self.profile = my_profile

    @abstractmethod
    def total_cross_section_area(self):
        """ return total cross section area"""
        pass

    @abstractmethod
    def cross_section_perimeter(self):
        """ return the cross section perimeter """
        pass

    @abstractmethod
    def return_area_data(self):
        """ return area data
        :return single element or tupple if its multi-dimensional"""
        pass

    @abstractmethod
    def return_perimeter_data(self):
        """ return perimeter data
        :return single element or tupple if its multi-dimensional """
        pass

    # TODO: get a proper definition for profile data returned by return_profile_data of Cell class
    @abstractmethod
    def return_profile_data(self):
        """
        return_profile_data returns a number associated to the profile that is representative of the cell
        :return: single number
        """
        pass


    @abstractmethod
    def update_min_thickness(self):
        """
        update the value of the minimum thickness of the cell
        :return: nothing
        """
        pass

    @abstractmethod
    def regress(self, local_r_dot, dt):
        """
        regress performs the regression of the cell
        :param local_r_dot: local regression rate
        :param dt: time-step
        :return: nothing
        """
        pass


class CircularPortCell(Cell):
    """
    CircularPortCell is a concrete class which inherits from Cell, only to be instantiated
    by the proper factory method
    """

    def __init__(self, num, x, profile, r_ext):
        """ class initializer """

        # Call superclass constructor
        super(CircularPortCell, self).__init__(num, x, profile, r_ext)

        # Set minimum thickness
        self.min_thickness = self.r_ext - self.profile

    def update_min_thickness(self):
        self.min_thickness = self.r_ext - self.profile

        # If thickness is smaller than 0 raise exception
        if self.min_thickness < 0:
            raise ValueError("Minimum Thickness cannot be smaller than 0 \n")

    def total_cross_section_area(self):
        """ Calculate cross section area for circular port cell """
        return m.pi * self.profile ** 2

    def cross_section_perimeter(self):
        """ Calculate the cross section perimeter for a circular port cell"""
        return 2 * m.pi * self.profile

    def return_area_data(self):
        return self.total_cross_section_area()

    def return_perimeter_data(self):
        return self.cross_section_perimeter()

    def return_profile_data(self):
        return self.profile

    def regress(self, local_r_dot, dt):
        # Subtract to the radius the regressed value
        self.profile += local_r_dot * dt
        # Update the minimum thickness
        self.update_min_thickness()
