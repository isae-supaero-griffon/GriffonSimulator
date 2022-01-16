# Geometries0D.py includes the definition of the 0D geometries solved by the GriffonSimulator combustion code.
# This includes the single-port as well as multi-port geometries.
# @author Jose Felix Zapata Usandivaras
# Date: 24/02/2020
# ISAE-SUPAERO Space Section, Project Griffon, Combustion team

# ------------------------- IMPORT MODULES ----------------------

from CombustionModule.Geometries import *               # Import all of the components from Geometries.py

# ------------------------ FUNCTION DEFINITIONS -----------------


def cylinder_volume(length, radius):
    return length * m.pi * (radius**2)


def n_star_volume(length, base_radius, int_radius, number_branches):
    return 2 * length * number_branches * base_radius * int_radius * sin(pi / number_branches)


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

    # Equation 13 used from https://mathworld.wolfram.com/FourierSeries.html
    # Generate the sum for both arrays
    sum_a = sum((a[k] * m.cos(2 * m.pi * k * x) for k in range(len(a))))
    sum_b = sum((b[k] * m.sin(2 * m.pi * k * x) for k in range(len(b))))

    # Return the output
    return sum_a + sum_b - a[0] / 2


def compute_point_image_position(image_pixel_size, pixel_base_radius, branches,
                                 i, k, n, modifier):
    return [image_pixel_size // 2 + pixel_base_radius * m.sin(2 * m.pi / branches * (i + k / n)) * modifier,
            image_pixel_size // 2 + pixel_base_radius * m.cos(2 * m.pi / branches * (i + k / n)) * modifier]

# -------------------------- CLASS DEFINITIONS ------------------


class Geometry0D(Geometry, ABC):
    """
    Geometry0D is an abstract class which inherits from the geometry
    class and implements the necessary methods for the 0D definition of it.

        Additional Attributes:
            1. port_number: number of ports
    """
    def __init__(self, length, regression_model, r_ext, port_number):
        """
        class initializer
        :param length: float with length of the grain [m]
        :param regression_model: RegressionModel instance
        :param r_ext: float with external radius [m]
        :param port_number: integer with number of ports in the grain
        """
        # Call superclass initializer
        super(Geometry0D, self).__init__(length, regression_model, r_ext)
        self.port_number = port_number

    def get_ports_number(self):
        """
        Return the number of ports for this geometry
        """
        return self.port_number

    def compute_fuel_rate(self, rho, ox_flow):
        """
        Return the instantaneous fuel regression rate of the geometry given flow properties.
        :param rho: solid fuel density
        :param ox_flow: instantaneous oxidizer flow
        """
        return rho * self.regression_model.compute_regression_rate(self, ox_flow) * self.total_surface_area()

    def solve_mass_flux(self, m_ox, rho_f):
        """
        solve_mass_flux solves for the mass flux of oxidizer and fuel
        :param m_ox: mass flow of oxidizer [kg/s]
        :param rho_f: fuel density in [kg/m^3]
        :return: solution_vector, mass_flows tupple, of_ratio, mass_fluxes tupple
        """
        m_f = self.compute_fuel_rate(rho_f, m_ox)               # Compute the fuel rate
        area = self.total_cross_section_area()                  # Compute the total cross section area
        g_ox, g_f = m_ox / area, m_f / area                     # Compute the mass fluxes
        solution = np.nan                                       # Set the solution (np.nan for 0D code)

        # Return output
        return solution, (m_ox, m_f), m_ox / m_f, (g_ox, g_f)

    @abstractmethod
    def total_cross_section_area(self):
        """
        Return the current cross section area of the geometry
        """
        pass

    @abstractmethod
    def total_surface_area(self):
        """
        Return the current total lateral surface area of the geometry
        """
        pass

    @abstractmethod
    def get_hydraulic_diameter(self):
        """
        Return the hydraulic diameter calculation of the geometry
        :return: hydraulic diameter of the geometry
        """
        pass

    @abstractmethod
    def draw_geometry(self):
        """
        Plot and show the geometry cross-section in pyplot
        """
        pass


class OneCircularPort(Geometry0D):
    """
    Implementation of the Geometry abstract class for a simple, one centered
    circular combustion port geometry.

        Additional attributes:
            1. r_int: circular port radius
            2. r_ext: total fuel slab radius
    """

    def __init__(self, length, regression_model, r_ext, r_int_initial):
        """
        class initializer
        :param length: float with length of the grain [m]
        :param regression_model: RegressionModel instance
        :param r_ext: float with external radius [m]
        :param r_int_initial: float with initial internal radius [m]
        """
        # Call superclass initializer
        super(OneCircularPort, self).__init__(length, regression_model, r_ext, 1)      # Set the port_number to 1
        self.initial_r = r_int_initial
        self.r_int = rintInitial

    def get_port_radius(self):
        """
        Return the current combustion port radius
        """
        return self.r_int

    def total_cross_section_area(self):

        return m.pi * (self.r_int**2)

    def total_surface_area(self):

        return self.length * 2 * m.pi * self.r_int

    def regress(self, ox_flow, dt, *args):
        self.r_int += self.regression_model.compute_regression_rate(self, ox_flow) * dt

    def min_bloc_thickness(self):

        return self.r_ext - self.r_int

    def draw_geometry(self):
        ax = plt.gca()                                                            # Generate an axes
        ax.add_patch(plt.Circle((0, 0), radius=self.r_ext * 1000, color='g'))     # Plot the outer shape
        draw_circular_port(ax, (0, 0), self)                                      # plot the combustion port
        plt.axis("scaled")                                                        # Adjust axis and show
        plt.show()

    def get_fuel_mass(self, fuel_density):
        return fuel_density * (cylinder_volume(self.length, self.r_ext) -
                               cylinder_volume(self.length, self.initial_r))

    def get_hydraulic_diameter(self):
        return 2 * self.r_int


class MultipleCircularPortsWithCircularCenter(Geometry0D):
    """
    Implementation of the Geometry0D abstract class for a more complex, multiple
    circular combustion ports surrounding one central circular port.
    Total port number must be >= 4

    Additional attributes:
        1. ring_ports: OneCircularGeometry object representing the outside ring ports
        2. central_port: OneCircularGeometry object representing the central port

    Those two attributes combined with port number encapsulate all geometrical aspects such as total radius.
    """

    def __init__(self, length, regression_model, r_ext, ports_number,
                 ring_ports_initial_radius, central_port_initial_radius):
        """
        class initializer
        :param length: float with length of grain [m]
        :param regression_model: RegressionModel instance
        :param r_ext: float with external radius [m]
        :param ports_number: integer with number of ports
        :param ring_ports_initial_radius: float with initial internal radius of ring ports [m]
        :param central_port_initial_radius: float with initial internal radius of central port [m]
        """
        # Check the attributes
        if N < 4:
            raise ValueError("Geometry must have at least 3 ring ports")

        # Call superclass initializer
        super(MultipleCircularPortsWithCircularCenter, self).__init__(length, regression_model,
                                                                      r_ext, ports_number)

        # Set the remaining attributes
        c = m.sin(m.pi / (N-1))  # Shape parameter linked to the number of ring ports (= nodes_number-1)
        self.ring_ports = OneCircularPort(length, regression_model, ring_ports_initial_radius, r_ext * c / (1+c))
        self.central_port = OneCircularPort(length, regression_model, central_port_initial_radius, r_ext -
                                            self.ring_ports.r_ext - ringPortsIntialRadius)

    def get_total_outer_radius(self):
        """
        Return the total outer radius of the fuel slab
        """
        return self.ring_ports.r_int + self.ring_ports.r_ext + self.central_port.r_ext

    def compute_fuel_rate(self, rho, ox_flow):
        """
        Return the instantaneous fuel regression rate of the geometry given flow properties.
        :param ox_flow: instantaneous oxidizer flow
        :param rho: density of the fuel [kg/m^3]
        """

        # Compute oxidizer flow in each port
        ring_individual_port_ox_flow = ox_flow * self.ring_ports.total_cross_section_area() / \
                                       self.total_cross_section_area()
        center_port_ox_flow = ox_flow * self.central_port.total_cross_section_area() / \
                              self.total_cross_section_area()

        return self.central_port.compute_fuel_rate(rho, center_port_ox_flow) + \
               (self.port_number - 1) * self.ring_ports.compute_fuel_rate(rho, ring_individual_port_ox_flow)

    def total_cross_section_area(self):
        return self.central_port.total_cross_section_area() + (self.port_number - 1) * \
               self.ring_ports.total_cross_section_area()

    def total_surface_area(self):

        return self.central_port.total_surface_area() + (self.port_number - 1) * self.ring_ports.total_surface_area()

    def regress(self, ox_flow, dt, *args):
        # Compute oxidizer flow in each port
        ring_individual_port_ox_flow = ox_flow * self.ring_ports.total_cross_section_area() / \
                                       self.total_cross_section_area()
        center_port_ox_flow = ox_flow * self.central_port.total_cross_section_area() / \
                              self.total_cross_section_area()

        initial_ring_radius = self.ring_ports.get_port_radius()

        # Proceed with regression
        self.central_port.regress(center_port_ox_flow, dt)
        self.ring_ports.regress(ring_individual_port_ox_flow, dt)
        self.central_port.r_ext += initial_ring_radius - self.ring_ports.get_port_radius()

    def min_bloc_thickness(self):

        return min(self.central_port.min_bloc_thickness(), self.ring_ports.min_bloc_thickness())

    def draw_geometry(self):
        ax = plt.gca()                                                           # Generate the axes
        ax.add_patch(plt.Circle((0, 0),
                                radius=self.get_total_outer_radius() * 1000,
                                color='g'))                                      # Plot the outer shape
        draw_circular_port(ax, (0, 0), self.central_port)                        # Plot the central port
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

    def get_hydraulic_diameter(self):
        return 2 * self.central_port.get_port_radius()


class ThreeCircularPorts(Geometry0D):
    """
    Implementation of the Geometry0D abstract class for a three circular
    combustion ports configuration.

    Additional attributes:
        1. ring_ports: OneCircularGeometry object representing the three identical ports.
        2. central_port: OneCircularGeometry object representing the central port.
    """

    def __init__(self, length, regression_model, r_ext, ports_initial_radius):
        """
        class initializer
        :param length: float with length of grain [m]
        :param regression_model: RegressionModel instance
        :param r_ext: float with external radius [m]
        :param ports_initial_radius: float with ports initial radius [m]
        """
        # Call superclass initializer
        super(ThreeCircularPorts, self).__init__(length, regression_model, r_ext, 3)

        # Set remaining attributes
        self.initial_r = ports_initial_radius
        self.ports = OneCircularPort(length, regression_model, r_ext / (1 + 2/m.sqrt(3)), portsIntialRadius)

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

    def compute_fuel_rate(self, rho, ox_flow):
        """
        Return the instantaneous fuel regression rate of the geometry given flow properties.
        :param rho: density of the grain [kg/m^3]
        :param ox_flow: instantaneous oxidizer flow [kg/s]
        """
        return 3 * self.ports.compute_fuel_rate(rho, ox_flow / 3)

    def total_cross_section_area(self):
        return 3 * self.ports.total_cross_section_area()

    def total_surface_area(self):
        return 3 * self.ports.total_surface_area()

    def regress(self, ox_flow, dt, *args):
        self.ports.regress(ox_flow / 3, dt)

    def min_bloc_thickness(self):
        return self.ports.min_bloc_thickness()

    def draw_geometry(self):
        ax = plt.gca()                                  # Generate the axes
        ax.add_patch(plt.Circle((0, 0),
                                radius=self.get_total_outer_radius() * 1000,
                                color='g'))             # Plot the outer shape

        # Plot the outer rings
        R = self.ports.r_ext * 2 / m.sqrt(3) * 1000
        for k in range(3):
            draw_circular_port(ax, (R*m.cos(2*k*m.pi/3), R*m.sin(2*k*m.pi/3)), self.ports)

        # Adjust axis and show
        plt.axis("scaled")
        plt.show()

    def get_fuel_mass(self, fuel_density):
        return fuel_density * (cylinder_volume(self.length, self.get_total_outer_radius()) -
                               3 * cylinder_volume(self.length, self.ports.initial_r))

    def get_hydraulic_diameter(self):
        return 2 * self.get_port_radius()


class SinglePortImageGeometry(Geometry0D):
    """
    SinglePortImageGeometry class inherits from the Geometry0D abstract class
    however using openCV library to estimate the regression of the block through parallel
    layers.
    """

    def __init__(self, length, regression_model,
                 r_ext, image_pixel_size,
                 image_meter_size):
        """
        class initializer
        :param length: float with length of grain [m]
        :param regression_model: RegressionModel instance
        :param r_ext: float external radius of the geometry [mm]
        :param image_pixel_size: pixel size of the image (the number of pixels being image_pixel_size^2)
        :param image_meter_size: size in meters of the image
        """
        # Call superclass initializer
        super(SinglePortImageGeometry, self).__init__(length, regression_model, r_ext, 1)

        # Set the remaining attributes
        self.image_pixel_size = image_pixel_size
        self.image_meter_size = image_meter_size
        self.image = np.zeros((image_pixel_size, image_pixel_size), np.uint8)
        self.is_port_geometry = False
        self.step_pixel_regression = 5
        self.kernel = self._set_kernel()

    def __str__(self):
        """ return a string of the objects representation """

        # Return the objects type
        class_name = self.__class__.__name__
        # Variables to avoid printing
        var_print = ('image', 'is_port_geometry', 'kernel')
        # Loop through the objects properties and values
        return "Geometry" + "\t\n\ttype, {0}\t\n\t".format(class_name) + \
               "\t\n\t".join(("{prop}, {val}".format(prop=prop, val=value) for prop, value in vars(self).items()
                              if prop not in var_print))

    @staticmethod
    def _set_kernel():
        """
        _set_kernel is a private static method which is used to initialize the kernel
        :return: kernel matrix
        """
        # Remark: kernel is tailored for 5 pixels/step regression
        size = 12                           # Define a size for the kernel
        R = size // 2                       # Get the radius of influence
        kernel = np.zeros((size, size))     # Initialize the matrix
        for i in range(size):
            for j in range(size):
                x = i - R
                y = j - R
                r = m.sqrt(x ** 2 + y ** 2)
                kernel[i, j] = max(1 - m.sqrt(r / R), 0)

        # Return the kernel
        return kernel

    def get_meters_per_pixel(self):
        return self.image_meter_size / self.image_pixel_size

    def total_cross_section_area(self):
        if self.is_port_geometry:
            return cv2.countNonZero(self.image) * (self.get_meters_per_pixel() ** 2)
        else:
            raise ValueError("Image is black : please generate geometry")

    def total_surface_area(self):
        return self.get_length() * self.get_perimeter()

    def get_perimeter(self):
        """ return the perimeter of the image """
        if self.is_port_geometry:
            perimeter = 0
            proxy_image = self.image
            contours, hierarchy = cv2.findContours(proxy_image, 0, 2)
            for cnt in contours:
                perimeter += cv2.arcLength(cnt, True) * self.get_meters_per_pixel()
            return perimeter
        else:
            raise ValueError("Image is black : please generate geometry")

    def get_hydraulic_diameter(self):
        # Calculate the perimeter and the cross-section, then compute hydraulic diameter
        return 4 * self.total_cross_section_area() / self.get_perimeter()

    def regress(self, ox_flow, dt, *args):
        """
        Apply regression to the geometry
        :param ox_flow: instantaneous oxidizer flow, ignored
        :param dt: time increment, ignored
        """
        # A 5 pixel radial regression is applied,
        # regardless of dt

        if self.is_port_geometry:
            self.image = cv2.filter2D(self.image, -1, self.kernel)
            ret, self.image = cv2.threshold(self.image, 20, 255, cv2.THRESH_BINARY)
        else:
            raise ValueError("Image is black : please generate geometry")

    def get_fuel_mass(self, fuel_density):
        """
        Return the initial fuel mass based on fuel_density
        :param fuel_density: density of the fuel [kg/m^3]
        """
        volume = self.get_length() * (m.pi * (self.r_ext**2) - self.total_cross_section_area())
        return fuel_density * volume

    def min_bloc_thickness(self):
        """
        Return the smallest fuel space between ports or fuel edge.
        Useful for setting an eburn termination safety margin.
        """
        if self.is_port_geometry:
            proxy_image = self.image
            contours, hierarchy = cv2.findContours(proxy_image, 1, 2)
            cnt = contours[0]
            _, radius = cv2.minEnclosingCircle(cnt)
            return self.r_ext - radius*self.get_meters_per_pixel()

        else:
            raise ValueError("Image is black : please generate geometry")

    def export_geometry(self, file_name='../data/myGeometry.txt'):
        """
        export_geometry generates an array of points x,y to describe the image of the geometry.
        :param file_name: name of the txt file to which we want to export the curve
        :return: nothing
        """
        def write_line_to_file(f_, couple_):
            """
            write a line to the given file identifier
            :param f_: file identifier
            :param couple_: parameters to print
            :return: nothing
            """
            f_.write("{x:5.5f} \t {y:5.5f}\n".format(x=1000 * (couple_[0][0] * self.get_meters_per_pixel() -
                                                               self.image_meter_size / 2),
                                                     y=1000 * (couple_[0][1] * self.get_meters_per_pixel() -
                                                               self.image_meter_size / 2)))

        if self.is_port_geometry:
            # Extract the first contour of the image
            proxy_image = self.image
            contours, _ = cv2.findContours(proxy_image, 1, 2)
            cnt = max(contours, key=lambda x: x.shape[0])

            # Write the file
            with open(file_name, 'w') as f:
                f.write("{x:10s} \t {y:10s} \n".format(x='X [mm]', y='Y [mm]'))
                for counter in range(0, len(cnt), 2):
                    couple = cnt[counter]
                    write_line_to_file(f, couple)
        else:
            raise ValueError("Image is black : please generate geometry")

    def draw_geometry(self):
        # Generate the plots
        fig = plt.figure(facecolor='w', figsize=(15, 15))
        fig.suptitle('Geometry drawing')
        plt.imshow(self.image, cmap='gray')

    def generate_polynom(self, polynom, base_radius, branches, n):

        polynom[-1] = 0                                             # Make sure there is no term of order 0 in polynom
        pixel_base_radius = m.floor(base_radius /
                                    self.get_meters_per_pixel())    # Get the base radius in pixels

        # Build the polygon
        points = [[self.image_pixel_size // 2, self.image_pixel_size // 2 + pixel_base_radius]]
        for i in range(branches):
            for k in range(n+1):
                modifier = 1 + np.polyval(polynom, k / n)
                points.append(compute_point_image_position(self.image_pixel_size, pixel_base_radius,
                                                           branches, i, k, n, modifier))

        points = np.floor(np.array(points))

        # Draw the shape
        cv2.fillPoly(self.image, np.int32([points]), 1, 255)

        # Indicate that a shape has been generated
        self.is_port_geometry = True

    def generate_fourier(self, a, b, base_radius, branches, impact, n):
        pixel_base_radius = m.floor(base_radius / self.get_meters_per_pixel())  # Get the base radius in pixels

        # Build the polygon
        points = []
        for i in range(branches):
            for k in range(n):
                modifier = 1 + impact * evaluate_fourier(a, b, k / n)
                points.append(compute_point_image_position(self.image_pixel_size, pixel_base_radius,
                                                           branches, i, k, n, modifier))
        points = np.floor(np.array(points))

        # Draw the shape
        cv2.fillPoly(self.image, np.int32([points]), 1, 255)

        # Indicate that a shape has been generated
        self.is_port_geometry = True

    def generate_multiple_circular(self, n, ring_ports_initial_radius, central_ports_initial_radius):
        # Draw the central port
        cv2.circle(self.image, (int(self.image_pixel_size // 2),
                                int(self.image_pixel_size // 2)),
                   int(central_ports_initial_radius // self.get_meters_per_pixel()), 255, -1)

        # Plot the ring ports
        if n > 0:
            R = self.r_ext / (1 + m.sin(m.pi / n)) // self.get_meters_per_pixel()
            for k in range(n):
                cv2.circle(self.image, (int(R * m.cos(2 * k * m.pi / n) + self.image_pixel_size // 2),
                                        int(R * m.sin(2 * k * m.pi / n)) + self.image_pixel_size // 2),
                           int(ring_ports_initial_radius // self.get_meters_per_pixel()), 255, -1)

        # Indicate that a shape has been generated
        self.is_port_geometry = True
