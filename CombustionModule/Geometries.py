# Geometry abstract class and its various implementations for hybrid rocket burn calculations
# Author: Maxime Sicat
# 28/12/2018

# ------------------------- IMPORT MODULES ----------------------

from abc import ABC, abstractmethod                 # Import abc module for abstract classes definitions
import matplotlib.pyplot as plt                     # Import matplotlib
import math as m                                    # Import the math module
import numpy as np                                  # Import numpy
import cv2                                          # Import cv2
from CombustionModule.Mesher import *               # Import the different kinds of meshes
from CombustionModule.RegressionModel import *      # Import the regression models
from scipy.integrate import solve_ivp               # Import solve_ivp method from scipy
import os.path                                      # Import the OS library
from Libraries.Interpolator import Interpolator     # Import the interpolator

# ------------------------ FUNCTION DEFINITIONS -----------------


def cylinder_volume(length, radius):
    return length * m.pi * (radius**2)


def n_star_volume(length, base_radius, int_radius, number_branches) :
    return 2 * length * number_branches * base_radius * int_radius * sin (pi / number_branches)


def compute_regression_rate(geometry, ox_flow, max_regression_rate, a, n, m):
    """
    Return the instantaneous fuel regression rate of the geometry given flow properties.
    :param geometry: Port geometry from which the regression rate should be calculated. Port number must be 1.
    :param ox_flow: instantaneous oxidizer flow
    :param max_regression_rate: maximum allowed regression rate until it reaches saturation
    :param a: classical regression rate factor (depends on fuel choice)
    :param n: classical regression rate exponent for oxidizer flux (depends on fuel choice)
    :param m: classical regression rate exponent for fuel length(usually set to -0.2)
    :raise ValueError: if port number is not 1, as regression rate has a physical meaning only
        if calculated in a single port.
    """

    if geometry.port_number != 1:

        raise ValueError("Geometry must be single port")

    Go = ox_flow / geometry.totalCrossSectionArea()

    # Compute the regression rate following the exponential regression law.
    regression_rate = a * (geometry.length ** m) * (Go ** n)

    # If the regression rate is higher than the maximum return the maximum, else use the exponential law
    if max_regression_rate < regression_rate:
        print("Warning: Exceeding Max allowed Regression Rate, Actual: {0:.2E} m/seg, Setting Cap "
              "to: {1:.2E} m/seg".format(regression_rate, max_regression_rate))
        regression_rate = max_regression_rate

    # Return the output
    return regression_rate


def evaluate_fourier(a, b, x):
    """
    evaluate_fourier gets the value of the fourier series for x, given the
    coefficients of the series.
    :param a: array of coefficients which multiply cos(x)
    :param b: array of coefficients which multiply sin(x)
    :param x: point where we want to evaluate the series
    :return: value of the sum
    """
    assert 0 <= x <= 1, "x variable not contained in [0, 1] \n"

    # Generate the sum for both arrays
    sum_a = sum((a[k] * m.cos(2 * m.pi * (k+1) * x) for k in range(len(a))))
    sum_b = sum((b[k] * m.sin(2 * m.pi * (k+1) * x) for k in range(len(b))))

    # Return the output
    return sum_a + sum_b


def draw_circular_port(ax, center, port):
    """
    draw_circular_port plots the circular port into an axes
    :param ax: axes object where to plot the port image
    :param center: coordinates of the center of the port
    :param port: Geometry instance associated to the port
    :return: nothing
    """
    # Plot the combustion port
    ax.add_patch(plt.Circle(center, radius=port.r_int * 1000, color='w'))


def draw_n_star_branch(p1_0, p2_0, p3_0, rotation_angle):
    """
    draw_n_star_branch draws a star
    :param p1_0: first vertices coordinate
    :param p2_0: 2nd vertices coordinate
    :param p3_0: 3rd vertices coordinate
    :param rotation_angle: angle we want to rotate the branch
    :return: nothing
    """

    # Define the rotation matrix
    rot = [[cos(rotation_angle), - sin(rotation_angle)],
           [sin(rotation_angle), cos(rotation_angle)]]

    # Rotate each of the points
    p1 = p1_0 * rot
    p2 = p2_0 * rot
    p3 = p3_0 * rot

    # Plot the combustion port
    plt.plot(p1, p3, 'w-', lw=2)
    plt.plot(p2, p3, 'w - ', lw=2)


def calculate_hydraulic_diameter(area, perimeter):
    """
    calculate_hydraulic_diameter implements the equation of the hydraulic for
    a given section
    :param area: area of the section of interest
    :param perimeter: perimeter of the section of interest
    :return: hydraulic diameter
    """
    return 4 * area / perimeter


# -------------------------- CLASS DEFINITIONS ------------------

class Geometry(ABC):
    """
    Geometry abstract class sets the common methods and implementation
    for all geometries meant to be tested in the Griffon hybrid burn simulation
    code.

    Attributes:
        1. length: Length of the fuel slab
        2. port_number: number of combustion ports for this geometry
        3. regressionModel: RegressionModel object for regression
        rate calculation
    """

    def __init__(self, L, N, regressionModel):
        """
        class initializer
        """

        super().__init__()
        self.length = L
        self.port_number = N
        self.regressionModel = regressionModel

    def __str__(self):
        """ return a string of the objects representation """

        # Return the objects type
        class_name = self.__class__.__name__
        # Loop through the objects properties and values
        return "Geometry" + "\t\n\ttype, {0}\t\n\t".format(class_name) + \
               "\t\n\t".join(("{prop}, {val}".format(prop=prop, val=value) for prop, value in vars(self).items()))

    def get_length(self):
        """
        Return the fuel slab length.
        """

        return self.length

    def number_of_ports(self):
        """
        Return the number of ports for this geometry
        """

        return self.port_number

    @abstractmethod
    def compute_fuel_rate(self, rho, ox_flow):
        """
        Return the instantaneous fuel mass flow rate for the geometry
        :param rho: fuel density
        """
        pass

    @abstractmethod
    def totalCrossSectionArea(self):
        """
        Return the current cross section area of the geometry
        """
        pass

    @abstractmethod
    def totalSurfaceArea(self):
        """
        Return the current total lateral surface area of the geometry
        """
        pass

    @abstractmethod
    def regress(self, ox_flow, dt):
        """
        Apply regression to the geometry
        :param ox_flow: instantaneous oxidizer flow
        :param max_regression_rate: maximum allowed regression rate for the fuel under consideration
        :param a: classical regression rate factor (depends on fuel choice)
        :param n: classical regression rate exponent for oxidizer flux (depends on fuel choice)
        :param m: classical regression rate exponent for fuel length(usually set to -0.2)
        :param dt: time increment
        """
        pass

    @abstractmethod
    def min_bloc_thickness(self):
        """
        Return the smallest fuel space between ports or fuel edge.
        Useful for setting an eburn termination safety margin.
        """
        pass

    @abstractmethod
    def draw_geometry(self):
        """
        Plot and show the geometry cross-section in pyplot
        """

        pass

    @abstractmethod
    def get_fuel_mass(self, fuel_density):
        """
        Return the initial fuel mass based on fuel_density
        :param fuel_density: density of the fuel [kg/m^3]
        """
        pass

    @abstractmethod
    def return_external_radius(self):
        """
        Return external radius of the geometry
        :return: external radius of the geometry
        """
        pass

    @abstractmethod
    def get_hydraulic_diameter(self):
        """
        Return the hydraulic diameter calculation of the geometry
        :return: hydraulic diameter of the geometry
        """
        pass


class OneCircularPort(Geometry):
    """
    Implementation of the Geometry abstract class for a simple, one centered
    circular combustion port geometry.

    Additional attributes:
    r_int: circular port radius
    r_ext: total fuel slab radius
    """

    def __init__(self, L, regressionModel, rintInitial, rext0):
        """
        class initializer
        """
        super().__init__(L, 1, regressionModel)
        # N parameter is set as 1 in Geometry initializer as there is a single combustion port
        self.initial_r = rintInitial
        self.r_int = rintInitial
        self.r_ext = rext0

    # Methods specific to this geometry

    def get_port_radius(self):
        """
        Return the current combustion port radius
        """
        return self.r_int

    # Abstract methods implementation

    def compute_fuel_rate(self, rho, ox_flow):
        """
        Return the instantaneous fuel regression rate of the geometry given flow properties.
        :param rho: solid fuel density
        :param ox_flow: instantaneous oxidizer flow
        """

        return rho * self.regressionModel.computeRegressionRate(self, ox_flow) * self.totalSurfaceArea()

    def totalCrossSectionArea(self):

        return m.pi * (self.r_int**2)

    def totalSurfaceArea(self):

        return self.length * 2 * m.pi * self.r_int

    def regress(self, ox_flow, dt):

        self.r_int += self.regressionModel.computeRegressionRate(self, ox_flow) * dt

        # print("Port radius : " + str(self.r_int*1000) + " mm")

    def min_bloc_thickness(self):

        return self.r_ext - self.r_int

    def draw_geometry(self):

        ax = plt.gca()

        # Plot the outer shape
        ax.add_patch(plt.Circle((0,0), radius=self.r_ext * 1000, color='g'))

        # plot the combustion port
        draw_circular_port(ax, (0,0), self)

        # Adjust axis and show
        plt.axis("scaled")
        plt.show()

    def get_fuel_mass(self, fuel_density):

        return fuel_density * (cylinder_volume(self.length, self.r_ext) - cylinder_volume(self.length, self.initial_r))

    def return_external_radius(self):
        return self.r_ext

    def get_hydraulic_diameter(self):
        return 2 * self.r_int


class MultipleCircularPortsWithCircularCenter(Geometry):
    """
    Implementation of the Geometry abstract class for a more complex, multiple
    circular combustion ports surrounding one central circular port.
    Total port number must be >= 4

    Additional attributes:
    ring_ports: OneCircularGeometry object representing the outside ring ports
    central_port: OneCircularGeometry object representing the central port

    Those two attributes combined with port number encapsulate all geometrical aspects such as total radius.
    """

    def __init__(self, N, L, regressionModel, ringPortsIntialRadius, centralPortInitialRadius, r_ext):
        """
        class initializer
        """

        if N < 4:

            raise ValueError("Geometry must have at least 3 ring ports")

        super().__init__(L, N, regressionModel)

        c = m.sin(m.pi / (N-1))  # Shape parameter linked to the number of ring ports (= N-1)

        self.r_ext = r_ext
        self.ring_ports = OneCircularPort(L, regressionModel, ringPortsIntialRadius, r_ext * c / (1+c))
        self.central_port = OneCircularPort(L, regressionModel, centralPortInitialRadius, r_ext - self.ring_ports.r_ext - ringPortsIntialRadius)

    # Methods specific to this geometry

    def get_total_outer_radius(self):
        """
        Return the total outer radius of the fuel slab
        """
        return self.ring_ports.r_int + self.ring_ports.r_ext + self.central_port.r_ext

    # Abstract methods implementation

    def compute_fuel_rate(self, rho, ox_flow):
        """
        Return the instantaneous fuel regression rate of the geometry given flow properties.
        :param ox_flow: instantaneous oxidizer flow
        :param rho: density of the fuel [kg/m^3]
        """

        # Compute oxidizer flow in each port
        ring_individual_port_ox_flow = ox_flow * self.ring_ports.totalCrossSectionArea() / self.totalCrossSectionArea()
        center_port_ox_flow = ox_flow * self.central_port.totalCrossSectionArea() / self.totalCrossSectionArea()

        return self.central_port.compute_fuel_rate(rho, center_port_ox_flow) + \
               (self.port_number - 1) * self.ring_ports.compute_fuel_rate(rho, ring_individual_port_ox_flow)

    def totalCrossSectionArea(self):

        return self.central_port.totalCrossSectionArea() + (self.port_number - 1) * self.ring_ports.totalCrossSectionArea()

    def totalSurfaceArea(self):

        return self.central_port.totalSurfaceArea() + (self.port_number - 1) * self.ring_ports.totalSurfaceArea()

    def regress(self, ox_flow, dt):

        # Compute oxidizer flow in each port
        ring_individual_port_ox_flow = ox_flow * self.ring_ports.totalCrossSectionArea() / self.totalCrossSectionArea()
        center_port_ox_flow = ox_flow * self.central_port.totalCrossSectionArea() / self.totalCrossSectionArea()

        initial_ring_radius = self.ring_ports.get_port_radius()

        self.central_port.regress(center_port_ox_flow, dt)
        self.ring_ports.regress(ring_individual_port_ox_flow, dt)

        self.central_port.r_ext += initial_ring_radius - self.ring_ports.get_port_radius()

    def min_bloc_thickness(self):

        return min(self.central_port.min_bloc_thickness(), self.ring_ports.min_bloc_thickness())

    def draw_geometry(self):

        ax = plt.gca()

        # Plot the outer shape
        ax.add_patch(plt.Circle((0, 0), radius=self.get_total_outer_radius() * 1000, color='g'))

        # Plot the central port
        draw_circular_port(ax, (0,0), self.central_port)

        # Plot the ring ports

        n = self.port_number - 1

        R = self.get_total_outer_radius() / (1 + m.sin(m.pi / n)) * 1000

        for k in range(n):

            draw_circular_port(ax, (R*m.cos(2*k*m.pi/n), R*m.sin(2*k*m.pi/n)), self.ring_ports)

        # Adjust axis and show
        plt.axis("scaled")
        plt.show()

    def get_fuel_mass(self, fuel_density):

        return fuel_density * (cylinder_volume(self.length, self.get_total_outer_radius()) - (self.port_number - 1) *
                               cylinder_volume(self.length, self.ring_ports.initial_r) -
                               cylinder_volume(self.length, self.central_port.initial_r))

    def return_external_radius(self):
        return self.r_ext

    def get_hydraulic_diameter(self):
        return 2 * self.central_port.get_port_radius()


class ThreeCircularPorts(Geometry):
    """
    Implementation of the Geometry abstract class for a three circular
    combustion ports configuration.

    Additional attributes:
    ring_ports: OneCircularGeometry object representing the three identical ports.
    central_port: OneCircularGeometry object representing the central port.
    """

    def __init__(self, L, regressionModel, portsIntialRadius, r_ext):
        """
        class initializer
        """

        super().__init__(L, 3, regressionModel)

        self.r_ext = r_ext
        self.initial_r = portsIntialRadius
        self.ports = OneCircularPort(L, regressionModel, portsIntialRadius, r_ext / (1 + 2/m.sqrt(3)))

    # Methods specific to this geometry

    def get_port_radius(self):
        """
        Return the current combustion port radius
        """
        return self.ports.get_port_radius()

    def get_total_outer_radius(self):
        """
        Return the total outer radius of the fuel slab
        """
        return self.ports.r_ext * (1 + 2/m.sqrt(3))

    # Abstract methods implementation

    def compute_fuel_rate(self, rho, ox_flow):
        """
        Return the instantaneous fuel regression rate of the geometry given flow properties.
        :param ox_flow: instantaneous oxidizer flow
        """

        return 3 * self.ports.compute_fuel_rate(rho, ox_flow / 3)

    def totalCrossSectionArea(self):

        return 3 * self.ports.totalCrossSectionArea()

    def totalSurfaceArea(self):

        return 3 * self.ports.totalSurfaceArea()

    def regress(self, ox_flow, dt):

        self.ports.regress(ox_flow / 3, dt)

    def min_bloc_thickness(self):

        return self.ports.min_bloc_thickness()

    def draw_geometry(self):

        ax = plt.gca()

        # Plot the outer shape
        ax.add_patch(plt.Circle((0, 0), radius=self.get_total_outer_radius() * 1000, color='g'))

        R = self.ports.r_ext * 2 / m.sqrt(3) * 1000

        for k in range(3):

            draw_circular_port(ax, (R*m.cos(2*k*m.pi/3), R*m.sin(2*k*m.pi/3)), self.ports)

        # Adjust axis and show
        plt.axis("scaled")
        plt.show()

    def get_fuel_mass(self, fuel_density):
        return fuel_density * (cylinder_volume(self.length, self.get_total_outer_radius()) -
                               3 * cylinder_volume(self.length, self.ports.initial_r))

    def return_external_radius(self):
        return self.r_ext

    def get_hydraulic_diameter(self):
        return 2 * self.get_port_radius()


class NBranchStarPort(Geometry):
    """
    Implementation of the Geometry abstract class for a more simple, one centered
    N branches star.
    Total branches number must be between 3 and 8

    Additional attributes:
    r_b: base circonscrit circle radius
    r_int :  complete star circonscrit circle radius
    r_ext: total fuel slab radius
    n : number of branches of the star geometry
    """

    def __init__(self, L, regressionModel,rint0, rext0, rb0, n0):
        """
        class initializer
        :param L: length of the block in [mm]
        :param regressionModel: RegressionModel instance to account for the regression law
        of the fuel
        :param rint0: internal radius of the star [mm]
        :param rext0: external redius of the geometry [mm]
        :param rb0: external radius which circunscribes de branches of the stars
        :param n0: number of branches in the star
        """

        if n0 < 3:
            raise ValueError("Geometry must have at least 3 branches")
        elif n0 > 8:
            raise ValueError("Geometry must have at most 8 branches")

        # N parameter is set as 1 in Geometry initializer as there is a single combustion port
        super().__init__(L, 1, regressionModel)
        self.r_b = rb0
        self.r_ext = rext0
        self.r_int = rint0
        self.n = n0

        # ------------------ Methods specific to this geometry ----------------------

    def get_port_base_radius(self):
        """
        Return the current combustion port base radius
        """
        return self.r_b

    def get_port_int_radius(self):
        """
        Return the current combustion port int radius
        """
        return self.r_int

    def get_port_ext_radius(self):
        """
        Return the current combustion port ext radius
        """
        return self.r_ext

    # -------------------- Abstract methods implementation ---------------------

    def compute_fuel_rate(self, rho, ox_flow):
        """
        Return the instantaneous fuel regression rate of the geometry given flow properties.
        :param ox_flow: instantaneous oxidizer flow
        :param rho: density of the fuel [kg/m^3]
        :return calculated fuel rate from regression law [kg/s]
        """
        return rho * self.regressionModel.computeRegressionRate(ox_flow) * self.totalSurfaceArea()

    def totalCrossSectionArea(self):
        return 2 * self.n * self.r_b * self.r_int * sin(pi / self.n)

    def totalSurfaceArea(self):
        return self.length * 2 * self.n * sqrt( (self.r_int - self.r_b * cos(pi / self.n)) ** 2 + (self.r_b * sin(pi / self.n)) ** 2)

    def regress(self, ox_flow, dt):
        self.r_int += self.regressionModel.computeRegressionRate(ox_flow) * dt
        self.r_b += self.regressionModel.computeRegressionRate(ox_flow) * dt


    def min_bloc_thickness(self):
        return self.r_ext - self.r_int

    def draw_geometry(self):
        p1_0 = [self.r_b / 2 * sin(pi / self.n), self.r_b * cos(pi / self.n)]
        p2_0 = [- self.r_b / 2 * sin(pi / self.n), - self.r_b * cos(pi / self.n)]
        p3_0 = [0, - self.r_int]

        for i in range(0, self.n):
            rotation_angle = pi / self.n * i
            draw_n_star_branch(p1_0, p2_0, p3_0, rotation_angle)

        # Ajust axis and show
        plt.axis("scaled")
        plt.show()


    def get_fuel_mass(self, fuel_density):

        return fuel_density * (cylinder_volume(self.length, self.r_ext) - n_star_volume(self.length, self.r_b, self.r_int, self.n))

    def get_hydraulic_diameter(self):
        return np.nan


class NBranchRectangleStarPort(Geometry):
    """
    Implementation of the Geometry abstract class for a more simple, one centered
    N branches star with rectangular branches.
    Total branches number must be between 3 and 8

    Additional attributes:
    r_b: base circonscrit circle radius
    r_int :  complete star circonscrit circle radius
    r_ext: total fuel slab radius
    n : number of branches of the star geometry
    """

    def __init__(self, L, regressionModel, rint0, rext0, rb0, n0):
        """
        class initializer
        """

        if n0 < 3:

            raise ValueError("Geometry must have at least 3 branches")

        elif n0 > 8:

            raise ValueError("Geometry must have at most 8 branches")

        super().__init__(L, regressionModel, 1)  # N parameter is set as 1 in Geometry initializer as there is a single combustion port
        self.r_b = rb0
        self.r_ext = rext0
        self.r_int = rint0
        self.n = n0

        # Methods specific to this geometry

    def get_port_base_radius(self):
        """
        Return the current combustion port base radius
        """
        return self.r_b

    def get_port_int_radius(self):
        """
        Return the current combustion port int radius
        """
        return self.r_int

    def get_port_ext_radius(self):
        """
        Return the current combustion port ext radius
        """
        return self.r_ext

    # Abstract methods implementation

    def compute_fuel_rate(self, rho, ox_flow):
        """
        Return the instantaneous fuel regression rate of the geometry given flow properties.
        :param ox_flow: instantaneous oxidizer flow
        :param rho: density of the fuel [kg/m^3]
        :return fuel mass flow that's regressed [kg/s]
        """

        return rho * self.regressionModel.computeRegressionRate(ox_flow) * self.totalSurfaceArea()

    def totalCrossSectionArea(self):
        return 2 * self.n * self.r_b * self.r_int * sin(pi / self.n) * sqrt (1 - (self.r_b / self.r_int * sin(pi / self.n)) ** 2)

    def totalSurfaceArea(self):
        return self.length * 2 * self.n * self.r_b * ( ((self.r_int / self.r_b) ** 2 - sin(pi / self.n) ** 2)**0.5 + sin(pi / self.n) - cos(pi / self.n))

    def regress(self, ox_flow, dt):
        self.r_int += self.regressionModel.computeRegressionRate(ox_flow) * dt
        self.r_b += self.regressionModel.computeRegressionRate(ox_flow) * dt

    def min_bloc_thickness(self):

        return self.r_ext - self.r_int

    def draw_geometry(self):

        p1_0 = [self.r_b / 2 * sin(pi / self.n), self.r_b * cos(pi / self.n)]
        p2_0 = [- self.r_b / 2 * sin(pi / self.n), - self.r_b * cos(pi / self.n)]
        p3_0 = [0, - self.r_int]

        for i in range(0, self.n):
            rotation_angle = pi / self.n * i
            draw_n_star_branch(p1_0, p2_0, p3_0, rotation_angle)

        # Ajust axis and show
        plt.axis("scaled")
        plt.show()

    def get_fuel_mass(self, fuel_density):
        return fuel_density * (cylinder_volume(self.length, self.r_ext) - n_star_volume(self.length, self.r_b, self.r_int, self.n) )


class SinglePortImageGeometry(Geometry):

    def __init__(self, L, externalRadius, regressionModel, imagePixelSize, imageMeterSize):
        """
        class initializer
        :param externalRadius: external radius of the geometry [mm]
        :param regressionModel: RegressionModel instance that characterizes the regression
        of the geometry
        :param imagePixelSize: pixel size of the image
        :param imageMeterSize: size in meters of the image
        """

        # Call superclass initializer
        super().__init__(L, 1, regressionModel)

        self.externalRadius = externalRadius
        self.imagePixelSize = imagePixelSize
        self.imageMeterSize = imageMeterSize
        self.image = np.zeros((imagePixelSize,imagePixelSize), np.uint8)
        self.portGeometryIsGenerated = False
        self.stepPixelRegression = 5

        # kernel is tailored for 5 pixels/step regression
        size = 12
        R = size // 2
        kernel = np.zeros((size, size))
        for i in range(size):
            for j in range(size):
                x = i - R
                y = j - R
                r = m.sqrt(x ** 2 + y ** 2)
                kernel[i, j] = max(1 - m.sqrt(r / R), 0)

        self.kernel = kernel

    def __str__(self):
        """ return a string of the objects representation """

        # Return the objects type
        class_name = self.__class__.__name__
        # Variables to avoid printing
        var_print = ('image', 'portGeometryIsGenerated', 'kernel')
        # Loop through the objects properties and values
        return "Geometry" + "\t\n\ttype, {0}\t\n\t".format(class_name) + \
               "\t\n\t".join(("{prop}, {val}".format(prop=prop, val=value) for prop, value in vars(self).items()
                              if prop not in var_print))

    def getMetersPerPixel(self):
        return self.imageMeterSize / self.imagePixelSize

    def totalCrossSectionArea(self):
        if self.portGeometryIsGenerated:
            return cv2.countNonZero(self.image) * (self.getMetersPerPixel() **2)
        else:
            raise ValueError("Image is black : please generate geometry")

    def totalSurfaceArea(self):
        return self.get_length() * self.get_perimeter()

    def get_perimeter(self):
        """ return the perimeter of the image """
        if self.portGeometryIsGenerated:
            perimeter = 0
            proxy_image = self.image
            contours, hierarchy = cv2.findContours(proxy_image, 0, 2)
            for cnt in contours :
                perimeter += cv2.arcLength(cnt, True) * self.getMetersPerPixel()
            return perimeter
        else:
            raise ValueError("Image is black : please generate geometry")

    def get_hydraulic_diameter(self):

        # Calculate the perimeter and the cross-section, then compute hydraulic diameter
        return 4 * self.totalCrossSectionArea() / self.get_perimeter()

    def compute_fuel_rate(self, rho, ox_flow):
        """
        Return the instantaneous fuel mass flow rate for the geometry
        :param ox_flow: total oxidizer flow [kg/s]
        :param rho: fuel density
        :return mass flow of fuel [kg/s]
        """
        regression_rate = self.regressionModel.computeRegressionRate(self, ox_flow)
        return rho * regression_rate * self.totalSurfaceArea()

    def regress(self, ox_flow, dt):
        """
        Apply regression to the geometry
        :param ox_flow: instantaneous oxidizer flow, ignored
        :param dt: time increment, ignored
        """
        # A 5 pixel radial regression is applied,
        # regardless of dt

        if self.portGeometryIsGenerated:
            self.image = cv2.filter2D(self.image, -1, self.kernel)
            ret, self.image = cv2.threshold(self.image, 20, 255, cv2.THRESH_BINARY)
        else:
            raise ValueError("Image is black : please generate geometry")

    def get_fuel_mass(self, fuel_density):
        """
        Return the initial fuel mass based on fuel_density
        :param fuel_density: density of the fuel [kg/m^3]
        """
        volume = self.get_length() * (m.pi * (self.externalRadius**2) - self.totalCrossSectionArea())
        return fuel_density * volume

    def return_external_radius(self):
        """
        Return external radius of the geometry
        :return: external radius of the geometry
        """
        return self.externalRadius

    def min_bloc_thickness(self):
        """
        Return the smallest fuel space between ports or fuel edge.
        Useful for setting an eburn termination safety margin.
        """

        if self.portGeometryIsGenerated:

            proxy_image = self.image
            contours, hierarchy = cv2.findContours(proxy_image, 1, 2)
            cnt = contours[0]
            (x, y), radius = cv2.minEnclosingCircle(cnt)

            return self.externalRadius - radius*self.getMetersPerPixel()

        else:
            raise ValueError("Image is black : please generate geometry")

    def export_geometry(self, file_name='../data/myGeometry.txt'):
        """
        export_geometry generates an array of points x,y to describe the image of the geometry.
        :param file_name: name of the txt file to which we want to export the curve
        :return: nothing
        """

        if self.portGeometryIsGenerated:

            # Extract the first contour of the image
            proxy_image = self.image
            contours, _ = cv2.findContours(proxy_image, 1, 2)
            cnt = contours[0]

            # Write the file
            with open(file_name, 'w') as f:
                f.write("{x:10s} \t {y:10s} \n".format(x='X [mm]', y='Y [mm]'))
                for counter in range(0, len(cnt), 2):
                    couple = cnt[counter]
                    f.write("{x:5.5f} \t {y:5.5f}\n".format(x=1000*(couple[0][0]*self.getMetersPerPixel() -
                                                                    self.imageMeterSize/2),
                                                             y=1000*(couple[0][1]*self.getMetersPerPixel() -
                                                                     self.imageMeterSize/2)))
        else:
            raise ValueError("Image is black : please generate geometry")


    def draw_geometry(self):

        # Generate the plots
        fig = plt.figure(facecolor='w', figsize=(15, 15))
        fig.suptitle('Geometry drawing')
        plt.imshow(self.image, cmap='gray')

    def generatePolynom(self, polynom, baseRadius, branches, n):

        polynom[-1] = 0 # Make sure there is no term of order 0 in polynom
        pixelBaseRadius = m.floor( baseRadius / self.getMetersPerPixel() ) # Get the base radius in pixels

        # Build the polygon

        points = [ [self.imagePixelSize//2, self.imagePixelSize//2 + pixelBaseRadius] ]

        for i in range(branches):

            for k in range(n+1):
                modifier = 1 + np.polyval(polynom, k / n)
                points.append([self.imagePixelSize // 2 + pixelBaseRadius * m.sin(2*m.pi/branches *(i + k/n)) * modifier,
                               self.imagePixelSize // 2 + pixelBaseRadius * m.cos(2*m.pi/branches *(i + k/n)) * modifier])

        points = np.floor(np.array(points))

        # Draw the shape
        cv2.fillPoly(self.image, np.int32([points]), 1, 255)

        # Indicate that a shape has been generated
        self.portGeometryIsGenerated = True

    def generateFourier(self, a, b, baseRadius, branches, impact, n):

        pixelBaseRadius = m.floor(baseRadius / self.getMetersPerPixel())  # Get the base radius in pixels

        # Build the polygon
        points = []

        for i in range(branches):

            for k in range(n):
                modifier = 1 + impact * evaluate_fourier(a, b, k / n)
                points.append([self.imagePixelSize // 2 + pixelBaseRadius * m.sin(2*m.pi/branches *(i + k/n)) * modifier,
                               self.imagePixelSize // 2 + pixelBaseRadius * m.cos(2*m.pi/branches *(i + k/n)) * modifier])

        points = np.floor(np.array(points))

        # Draw the shape
        cv2.fillPoly(self.image, np.int32([points]), 1, 255)

        # Indicate that a shape has been generated
        self.portGeometryIsGenerated = True

    def generateMultipleCircular(self, n, ringPortsInitialRadius, centralPortInitialRadius):

        # Draw the central port
        cv2.circle(self.image, (int(self.imagePixelSize//2), int(self.imagePixelSize//2)), int(centralPortInitialRadius // self.getMetersPerPixel()), 255, -1)

        # Plot the ring ports
        if n > 0:
            R = self.externalRadius / (1 + m.sin(m.pi / n)) // self.getMetersPerPixel()

            for k in range(n):

                cv2.circle(self.image, (int(R * m.cos(2 * k * m.pi / n) + self.imagePixelSize//2),
                                        int(R * m.sin(2 * k * m.pi / n)) + self.imagePixelSize//2), int(ringPortsInitialRadius // self.getMetersPerPixel()), 255, -1)

        # Indicate that a shape has been generated
        self.portGeometryIsGenerated = True


# ---------------------- 1D Geometry definition

class Geometry1D(ABC):
    """
    The Geometry1D class is an abstract class aimed at providing 1D simulation
    capabilities to the Griffon Simulator, based on the Haltman's book.

    Attributes:
        1. Mesh: Mesh class instance
        2. r_ext: External cylinder radius (all grains are assumed to be cylindrical
        3. min_thickness: min material thickness with the external part of the cylinder
        3. interpolator: method used to aid in the interpolation of the geometries
    """

    #TODO: associate the Geometry1D model with the Geometry - 0D models in a redesign of the code

    def __init__(self, L, N, r_ext):
        """
        Class initializer
        :param L: float containing the grain-length
        :param r_ext: float indicating the external radius of the geometry
        :param N: integer defining the number of cells in the mesh
        """
        # Call superclass initializer
        super(Geometry1D, self).__init__()

        # Check the inputs
        assert  L > 0, "Grain length (L) must be greater than 0 \n"
        assert  r_ext > 0, "Grain external radius (r_ext) must be greater than 0 \n"
        assert  N > 0, "Number of nodes (N) must be greater than 0 \n"
        assert  isinstance(N, int), "Number of nodes (N) must be an integer \n"

        # Set attributes
        self.L = L
        self.r_ext = r_ext
        self.interpolator = None

    def _generate_mesh(self, N):
        """
        _generate_mesh is a private method which generates a uniformly spaced mesh
        :param N: integer indicating the number of cells of the mesh
        :return: UniformSpacedMesh instance
        """
        return UniformlySpacedMesh("Mesh-1D", self, N)
        # return GeometricMesh("Mesh-1D", self, N, bias=1.01)

    @abstractmethod
    def _generate_interpolator(self):
        """
        _generate_interpolator handles the geneation of the interpolator object.
        :return: interpolator object or collection of interpolators
        """

    @abstractmethod
    def my_cell_factory(self, i, x):
        """
        my_cell_factory is a method which helps instantiate the cells associated to the geometry
        :param i: integer indicating cell number
        :param x: cell position
        :return: Cell instance of interest
        """
        pass

    @abstractmethod
    def _update_interpolator(self):
        """ _update_interpolator updates the values present in the interpolator to the
        values present in the cells.
        """
        pass

    @abstractmethod
    def print_geometry_to_file(self, file_name, time=0):
        """ print_geometry_to_file prints the local current state of the geometry
        :param file_name: string containing name of file towards which to print the geometry
        :param time: simulation time at which to print the geometry
        """
        pass

    @abstractmethod
    def solve_mass_flux(self, regression_model, m_ox, rho_f):
        """
        solve_mass_flux solves the equation by using the RK4 method
        :param regression_model: RegressionModel instance.
        :param m_ox: mass flow of oxidizer
        :param rho_f: fuel density
        :return: solution_vector, mass_flows tupple, mass_fluxes tupple, O/F ratio at the end of port,
        """
        pass

    @abstractmethod
    def regress(self, solution, regression_model, m_ox, dt):
        """
        regress method performs the regression of the fuel grain. Updates the geometry
        along the port by updating the cells on the mesh from the solution obtained in solve_mass_flux
        :param solution: vector of fuel mass fluxes obtained from solve_mass_flux
        :param regression_model: RegressionModel instance.
        :param m_ox: mass flow of oxidizer
        :param dt: time of simulation
        :return: nothing
        """
        pass

    def min_bloc_thickness(self):
        """
        min_block_thickness returns the minimum thickness of the block
        :return: min_block thickness
        """
        min_block_thickness = self.r_ext
        for my_cell in self.mesh.cells:
            if my_cell.min_thickness < min_block_thickness:
                min_block_thickness = my_cell.min_thickness

        # Return the min block thickness
        return min_block_thickness



class SingleCircularPort1D(Geometry1D):
    """
    SingleCircularPort1D class inherits from the Geometry1D class and implements
    the case of a single circular port. For such case then, the profile expected for the
    port is a circle of constant radius.

    Attributes:
        1. r_int: float indicating the internal radius
    """

    def __init__(self, L, r_ext, N, r_init):
        """
        class initializer
        :param r_init: float indicating the internal radius.
        """

        # Call parent class constructor
        super(SingleCircularPort1D, self).__init__(L, N, r_ext)

        # Check the inputs
        assert  r_init > 0, "Internal radius (r_init) has to be greater than 0 \n"

        # Set the properties
        self.r_init = r_init
        self.mesh = self._generate_mesh(N)
        self.interpolator = self._generate_interpolator()

    def my_cell_factory(self, i, x):
        """ implement cell factory for single circular port """
        return CircularPortCell(i, x, self.r_init, self.r_ext)

    def _generate_interpolator(self):
        """ generate the interpolator """

        # Collect the area data
        x, areas, perimeters = self.mesh.return_data()

        # Return the output
        # return [Interpolator(0, self.L, x, areas, perimeters)]
        return [Interpolator(x, areas), Interpolator(x, perimeters)]

    def _update_interpolator(self):
        # Collect the area data
        x, updated_areas, updated_perimeters = self.mesh.return_data()

        # Set the new areas
        self.interpolator[0].set_coordinates(x, updated_areas)
        self.interpolator[1].set_coordinates(x, updated_perimeters)

    def print_geometry_to_file(self, file_name, time=0):
        # Append the data to the file
        if not os.path.isfile(file_name):
            data_set = [(my_cell.x_cor, my_cell.return_profile_data()) for my_cell in self.mesh.cells]
            x, profile_data = zip(*data_set)
            x = [-1] + list(x)
            profile_data = [time] + list(profile_data)
            data_set = np.column_stack((np.array(x), np.array(profile_data)))
        else:
            data_set = np.genfromtxt(file_name)
            profile_data = np.array([time]+[my_cell.return_profile_data() for my_cell in self.mesh.cells])
            data_set = np.hstack((data_set, np.reshape(profile_data, (len(profile_data), 1))))

        # Write the file
        np.savetxt(file_name, data_set, delimiter=" ", fmt="%15.12f")

    def solve_mass_flux(self, regression_model, m_ox, rho_f):
        """ perform the solution of the equation of regression of the fuel geometry """

        # Define the Runge-Kutta function
        def mass_flux_gradient(x, mass_flux_fuel):
            """
            mass_flux_gradient is the function that is inputted into the solver
            :param x: x_coordinate
            :param mass_flux_fuel: mass flux of fuel at coordinate x
            :return: mass flux gradient
            """

            # Interpolate the area and perimeter from the mesh
            area, perimeter = self.interpolator[0].interpolate(x), self.interpolator[1].interpolate(x)

            # Calculate fluxes and hydraulic_diameter
            mass_flux_ox = m_ox / area
            hydraulic_diameter = calculate_hydraulic_diameter(area, perimeter)
            total_mass_flux = mass_flux_fuel + mass_flux_ox

            # Return the output of the fuel flux
            return 4 * rho_f * regression_model.compute_regression_rate_haltman(x, total_mass_flux) / hydraulic_diameter

        # Implement the Runge-Kutta solver for the cells coordinates
        x_cor = self.mesh.return_x_cor()

        # Problem defined as follows
        solution = solve_ivp(fun=mass_flux_gradient, t_span=(x_cor[0], self.L), y0=[0], method='RK45', t_eval=x_cor,
                             rtol=1e-3)

        # Compute the mass flow rates and mass fluxes and OF at the end of the port
        end_cell = self.mesh.cells[-1]
        end_area = end_cell.return_area_data()
        g_f = solution.y[0][end_cell.number]
        g_ox = m_ox / end_area
        m_f = g_f * end_area

        # Return the results
        return solution, (m_ox, m_f), m_ox / m_f, (g_ox, g_f)

    def regress(self, solution, regression_model, m_ox, dt):
        # Update the profile by regressing it
        for my_cell in self.mesh.cells:
            g_ox = m_ox / my_cell.return_area_data()                                           # Get ox flux
            g_f = solution.y[0][my_cell.number]                                                # Get fuel flux
            g_tot = g_f  +  g_ox                                                               # Get the total mass flux
            r_dot = regression_model.compute_regression_rate_haltman(my_cell.x_cor, g_tot)     # Calculate regression r
            my_cell.regress(r_dot, dt)

        # Update the interpolator to account for the changes in the geometry
        self._update_interpolator()















