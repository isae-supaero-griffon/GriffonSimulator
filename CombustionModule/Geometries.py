# Geometry abstract class and its various implementations for hybrid rocket burn calculations
# Author: Maxime Sicat
# 28/12/2018

# ------------------------- IMPORT MODULES ----------------------

from abc import ABC, abstractmethod
import matplotlib.pyplot as plt
import math as m
import numpy as np
import cv2

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

    return regression_rate

def evaluate_fourier(a, b, x):

    assert(0 <= x <= 1)

    sum = 0

    for k in range(len(a)):

        sum += a[k] * m.cos(x * (k+1) * 2 * m.pi)

    for k in range(len(b)):

        sum += b[k] * m.sin(x * (k+1) * 2 * m.pi)

    return sum


def draw_circular_port(ax, center, port):

    # Plot the combustion port
    ax.add_patch(plt.Circle(center, radius=port.r_int * 1000, color='w'))


def draw_n_star_branch(p1_0, p2_0, p3_0, rotation_angle):

    rot = [[cos(rotation_angle), - sin(rotation_angle)],
           [sin(rotation_angle), cos(rotation_angle)]]

    p1 = p1_0 * rot
    p2 = p2_0 * rot
    p3 = p3_0 * rot

    # Plot the combustion port
    plt.plot(p1, p3, 'w-', lw=2)
    plt.plot(p2, p3, 'w - ', lw=2)

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

    def __init__(self, L, rint0, rext0, rb0, n0):
        """
        class initializer
        """

        if n0 < 3:

            raise ValueError("Geometry must have at least 3 branches")

        elif n0 > 8:

            raise ValueError("Geometry must have at most 8 branches")


        super().__init__(L, 1)  # N parameter is set as 1 in Geometry initializer as there is a single combustion port
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

    def compute_fuel_rate(self, rho, ox_flow, a, n, m):
        """
        Return the instantaneous fuel regression rate of the geometry given flow properties.
        :param ox_flow: instantaneous oxidizer flow
        :param a: classical regression rate factor (depends on fuel choice)
        :param n: classical regression rate exponent for oxidizer flux (depends on fuel choice)
        :param m: classical regression rate exponent for fuel length(usually set to -0.2)
        """

        return rho * compute_regression_rate(self, ox_flow, a, n, m) * self.totalSurfaceArea()

    def totalCrossSectionArea(self):

        return 2 * self.n * self.r_b * self.r_int * sin(pi / self.n)

    def totalSurfaceArea(self):

        return self.length * 2 * self.n * sqrt( (self.r_int - self.r_b * cos(pi / self.n)) ** 2 + (self.r_b * sin(pi / self.n)) ** 2)

    def regress(self, ox_flow, a, n, m, dt):

        self.r_int += compute_regression_rate(self, ox_flow, a, n, m) * dt
        self.r_b += compute_regression_rate(self, ox_flow, a, n, m) * dt


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

    def __init__(self, L, rint0, rext0, rb0, n0):
        """
        class initializer
        """

        if n0 < 3:

            raise ValueError("Geometry must have at least 3 branches")

        elif n0 > 8:

            raise ValueError("Geometry must have at most 8 branches")

        super().__init__(L, 1)  # N parameter is set as 1 in Geometry initializer as there is a single combustion port
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

    def compute_fuel_rate(self, rho, ox_flow, a, n, m):
        """
        Return the instantaneous fuel regression rate of the geometry given flow properties.
        :param ox_flow: instantaneous oxidizer flow
        :param a: classical regression rate factor (depends on fuel choice)
        :param n: classical regression rate exponent for oxidizer flux (depends on fuel choice)
        :param m: classical regression rate exponent for fuel length(usually set to -0.2)
        """

        return rho * compute_regression_rate(self, ox_flow, a, n, m) * self.totalSurfaceArea()

    def totalCrossSectionArea(self):

        return 2 * self.n * self.r_b * self.r_int * sin(pi / self.n) * sqrt (1 - (self.r_b / self.r_int * sin(pi / self.n)) ** 2)

    def totalSurfaceArea(self):

        return self.length * 2 * self.n * self.r_b * ( ((self.r_int / self.r_b) ** 2 - sin(pi / self.n) ** 2)**0.5 + sin(pi / self.n) - cos(pi / self.n))

    def regress(self, ox_flow, a, n, m, dt):

        self.r_int += compute_regression_rate(self, ox_flow, a, n, m) * dt
        self.r_b += compute_regression_rate(self, ox_flow, a, n, m) * dt


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
        """

        super().__init__(L, 1, regressionModel)

        self.externalRadius = externalRadius
        self.imagePixelSize = imagePixelSize
        self.imageMeterSize = imageMeterSize
        self.image = img = np.zeros((imagePixelSize,imagePixelSize), np.uint8)
        self.portGeometryIsGenerated = False
        self.stepPixelRegression = 5

        # kernel is tailored for 5 pixels/step regression
        size = 12
        R = size // 2
        print("R = ", R)
        kernel = np.zeros((size, size))
        for i in range(size):
            for j in range(size):
                x = i - R
                y = j - R
                r = m.sqrt(x ** 2 + y ** 2)
                kernel[i, j] = max(1 - m.sqrt(r / R), 0)

        self.kernel = kernel

    def getMetersPerPixel(self):

        return self.imageMeterSize / self.imagePixelSize

    def totalCrossSectionArea(self):

        if (self.portGeometryIsGenerated):

            proxy_image = self.image

            return cv2.countNonZero(self.image) * (self.getMetersPerPixel() **2)

        else:
            raise ValueError("Image is black : please generate geometry")

    def totalSurfaceArea(self):

        if (self.portGeometryIsGenerated):

            proxy_image = self.image

            contours, hierarchy = cv2.findContours(proxy_image, 0, 2)
            cnt = contours[0]

            perimeter = cv2.arcLength(cnt, True) * self.getMetersPerPixel()

            return perimeter * self.get_length()

        else:
            raise ValueError("Image is black : please generate geometry")


    def compute_fuel_rate(self, rho, ox_flow):
        """
        Return the instantaneous fuel mass flow rate for the geometry
        :param rho: fuel density
        """
        regressionRate = self.regressionModel.computeRegressionRate(self, ox_flow)

        return rho * regressionRate * self.totalSurfaceArea()

    def regress(self, ox_flow, dt):
        """
        Apply regression to the geometry
        :param ox_flow: instantaneous oxidizer flow, ignored
        :param dt: time increment, ignored
        """
        # A 5 pixel radial regression is applied,
        # regardless of dt

        if (self.portGeometryIsGenerated):

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

        if (self.portGeometryIsGenerated):

            proxy_image = self.image
            contours, hierarchy = cv2.findContours(proxy_image, 1, 2)
            cnt = contours[0]
            (x, y), radius = cv2.minEnclosingCircle(cnt)

            return self.externalRadius - radius*self.getMetersPerPixel()

        else:
            raise ValueError("Image is black : please generate geometry")

    def draw_geometry(self):

        plt.imshow(self.image, cmap='gray')
        plt.show()

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