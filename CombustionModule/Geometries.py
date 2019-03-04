# Geometry abstract class and its various implementations for hybrid rocket burn calculations
# Author: Maxime Sicat
# 28/12/2018

# ------------------------- IMPORT MODULES ----------------------

from abc import ABC, abstractmethod
import matplotlib.pyplot as plt
import math as m

# ------------------------ FUNCTION DEFINITIONS -----------------


def cylinder_volume(length, radius):

    return length * m.pi * (radius**2)


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


def draw_circular_port(ax, center, port):

    # Plot the combustion port
    ax.add_patch(plt.Circle(center, radius=port.r_int * 1000, color='w'))


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

        return fuel_density * (cylinder_volume(self.length, self.r_ext) - cylinder_volume(self.length, self.r_int))

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
                               cylinder_volume(self.length, self.ring_ports.r_int) -
                               cylinder_volume(self.length, self.central_port.r_int))

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
                                3 * cylinder_volume(self.length, self.ports.r_int))

    def return_external_radius(self):
        return self.r_ext


class SinglePortImageGeometry(Geometry):

    def __init__(self, L, externalRadius, regressionModel, imagePixelSize, imageMeterSize):
        """
        class initializer
        """

        super().__init__(L, 1, regressionModel)

        self.externalRadius = externalRadius
        self.imagePixelSize = imagePixelSize
        self.imageMeterSize = imageMeterSize
        self.image = np.zeros((imagePixelSize, imagePixelSize, 1), np.uint8)
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

            contours, hierarchy = cv2.findContours(self.image, 1, 2)
            cnt = contours[0]
            return cv2.contourArea(cnt)

        else:
            raise ValueError("Image is black : please generate geometry")

    def totalSurfaceArea(self):

        if (self.portGeometryIsGenerated):

            contours, hierarchy = cv2.findContours(self.image, 1, 2)
            cnt = contours[0]
            perimeter = cv2.arcLength(cnt, True)

            return perimeter * 2 * m.pi() * self.get_length()

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
        :param ox_flow: instantaneous oxidizer flow
        :param dt: time increment, ignored
        """
        # A 5 pixel radial regression is applied,
        # regardless of dt

        if (self.portGeometryIsGenerated):

            blurred = cv2.filter2D(self.image, -1, self.kernel)
            ret, self.image = cv2.threshold(blurred, 30, 255, cv2.THRESH_BINARY)

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

            contours, hierarchy = cv2.findContours(self.image, 1, 2)
            cnt = contours[0]
            (x, y), radius = cv2.minEnclosingCircle(cnt)

            return self.externalRadius - radius

        else:
            raise ValueError("Image is black : please generate geometry")

