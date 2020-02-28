# Geometries1D.py includes the definition of the mono-dimensional geometry codes used in GriffonSimulator.py
# @author Jose Felix Zapata Usandivaras
# Date: 24/02/2020
# ISAE-SUPAERO Space Section, Project Griffon, Combustion Team

# ------------------------- IMPORT MODULES ----------------------

from CombustionModule.Geometries import *           # Import all of the components from Geometries
from CombustionModule.Mesher import *               # Import the different kinds of meshes

# ------------------------ FUNCTION DEFINITIONS -----------------


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


class Geometry1D(Geometry, ABC):
    """
    The Geometry1D class is an abstract class aimed at providing 1D simulation
    capabilities to the Griffon Simulator, based on the Haltman's book.

    Attributes:
        1. Mesh: Mesh class instance
        2. r_ext: External cylinder radius (all grains are assumed to be cylindrical
        3. min_thickness: min material thickness with the external part of the cylinder
        3. interpolator: method used to aid in the interpolation of the geometries
    """

    def __init__(self, length, nodes_number, r_ext, regression_model):
        """
        Class initializer
        :param length: float containing the grain-length
        :param r_ext: float indicating the external radius of the geometry
        :param nodes_number: integer defining the number of cells in the mesh
        """
        # Check the inputs
        assert nodes_number > 0, "Number of nodes (nodes_number) must be greater than 0 \n"
        assert isinstance(nodes_number, int), "Number of nodes (nodes_number) must be an integer \n"

        # Call superclass initializer
        super(Geometry1D, self).__init__(length, regression_model, r_ext)

        # Set attributes
        self.interpolator = None

    def _generate_mesh(self, nodes_number):
        """
        _generate_mesh is a private method which generates a uniformly spaced mesh
        :param nodes_number: integer indicating the number of cells of the mesh
        :return: UniformSpacedMesh instance
        """
        return UniformlySpacedMesh("Mesh-1D", self, nodes_number)

    def _generate_interpolator(self):
        """
        _generate_interpolator handles the geneation of the interpolator object.
        :return: interpolator object or collection of interpolators
        """
        # Collect the area data
        x, areas, perimeters = self.mesh.return_data()

        # Return the output
        return [Interpolator(x, areas), Interpolator(x, perimeters)]

    @abstractmethod
    def my_cell_factory(self, i, x):
        """
        my_cell_factory is a method which helps instantiate the cells associated to the geometry
        :param i: integer indicating cell number
        :param x: cell position
        :return: Cell instance of interest
        """
        pass

    def _update_interpolator(self):
        """ _update_interpolator updates the values present in the interpolator to the
        values present in the cells.
        """
        # Collect the area data
        x, updated_areas, updated_perimeters = self.mesh.return_data()

        # Set the new areas
        self.interpolator[0].set_coordinates(x, updated_areas)
        self.interpolator[1].set_coordinates(x, updated_perimeters)

    def print_geometry_to_file(self, file_name, time=0):
        """ print_geometry_to_file prints the local current state of the geometry
        :param file_name: string containing name of file towards which to print the geometry
        :param time: simulation time at which to print the geometry
        """
        # Append the data to the file
        if not os.path.isfile(file_name):
            data_set = [(my_cell.x_cor, my_cell.return_profile_data()) for my_cell in self.mesh.cells]
            x, profile_data = zip(*data_set)
            x = [-1] + list(x)
            profile_data = [time] + list(profile_data)
            data_set = np.column_stack((np.array(x), np.array(profile_data)))
        else:
            data_set = np.genfromtxt(file_name)
            profile_data = np.array([time] + [my_cell.return_profile_data() for my_cell in self.mesh.cells])
            data_set = np.hstack((data_set, np.reshape(profile_data, (len(profile_data), 1))))

        # Write the file
        np.savetxt(file_name, data_set, delimiter=" ", fmt="%15.12f")

    def export_geometry(self, file_name):
        # Extract the data
        data_set = [(my_cell.x_cor, my_cell.return_profile_data()) for my_cell in self.mesh.cells]

        # Write the file
        with open(file_name, 'w') as f:
            f.write("{x:10s} \t {y:10s} \n".format(x='X [mm]', y='Y [mm]'))
            for x, data in data_set:
                f.write("{x:5.5f} \t {y:5.5f}\n".format(x=1000*x, y=1000*data))

    def get_fuel_mass(self, fuel_density):
        # Proceed to integrate the area to obtain the internal volume
        x, areas, _ = self.mesh.return_data()
        internal_volume = np.trapz(areas, x)
        return fuel_density * (np.pi * self.length * self.r_ext ** 2 - internal_volume)

    def compute_fuel_rate(self, rho, ox_flow):
        """
        compute_fuel_rate determines the rate of fuel that is regressed from the
        grain. This method is obsolete.
        :param rho: density of the fuel [kg/m^3]
        :param ox_flow: oxidizer flow at the beginning of the port
        :return: fuel rate [kg/s]
        """
        _, m_flows, _, _ = self.solve_mass_flux(ox_flow, rho)
        return m_flows[1]

    def solve_mass_flux(self, m_ox, rho_f):
        """
        solve_mass_flux solves the equation by using the RK4 method
        :param m_ox: mass flow of oxidizer
        :param rho_f: fuel density
        :return: solution_vector, mass_flows tupple, of_ratio, mass_fluxes tupple
        """
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
            return 4 * rho_f * self.regression_model.compute_regression_rate_haltman(x,
                                                                                     total_mass_flux) / \
                   hydraulic_diameter

        # Implement the Runge-Kutta solver for the cells coordinates
        x_cor = self.mesh.return_x_cor()

        # Problem defined as follows
        solution = solve_ivp(fun=mass_flux_gradient, t_span=(x_cor[0], self.length), y0=[0], method='RK45', t_eval=x_cor,
                             rtol=1e-3)

        # Compute the mass flow rates and mass fluxes and OF at the end of the port
        end_cell = self.mesh.cells[-1]
        end_area = end_cell.return_area_data()
        g_f = solution.y[0][end_cell.number]
        g_ox = m_ox / end_area
        m_f = g_f * end_area

        # Return the results
        return solution, (m_ox, m_f), m_ox / m_f, (g_ox, g_f)

    def regress(self, ox_flow, dt, *args):
        """
        Apply regression to the geometry
        :param ox_flow: instantaneous oxidizer flow
        :param dt: time increment
        :param args: additional arguments (depending on the class) - in this case it has to be the
        solution of the regression
        """
        # Update the profile by regressing it
        for my_cell in self.mesh.cells:
            g_ox = ox_flow / my_cell.return_area_data()                                         # Get ox flux
            g_f = args[0].y[0][my_cell.number]                                                  # Get fuel flux
            g_tot = g_f + g_ox                                                                  # Get the total flux
            r_dot = self.regression_model.compute_regression_rate_haltman(my_cell.x_cor,
                                                                          g_tot)                # Calculate regression r
            my_cell.regress(r_dot, dt)

        # Update the interpolator to account for the changes in the geometry
        self._update_interpolator()

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

    def get_grain_profile(self):
        """
        get_grain_profile returns two numpy arrays which have the profile distribution of the Geometry1D
        :return: x coordinates numpy array, profile data numpy array
        """
        return self.mesh.return_x_cor(), self.mesh.return_profile_data()

    def total_cross_section_area(self):
        """ return the cross section area of the last cell """
        return self.interpolator[0].y_cor[-1]


class SingleCircularPort1D(Geometry1D):
    """
    SingleCircularPort1D class inherits from the Geometry1D class and implements
    the case of a single circular port. For such case then, the profile expected for the
    port is a circle of constant radius.

    Attributes:
        1. r_int: float indicating the internal radius
    """

    def __init__(self, length, r_ext, nodes_number, regression_model, r_init):
        """
        class initializer
        :param r_init: float indicating the internal radius.
        """

        # Call parent class constructor
        super(SingleCircularPort1D, self).__init__(length, nodes_number, r_ext, regression_model)

        # Check the inputs
        assert r_ext > r_init > 0, "Internal radius (r_init) has to be greater than 0 and smaller than External " \
                                   "radius (r_ext).\n"

        # Set the properties
        self.r_init = r_init
        self.mesh = self._generate_mesh(nodes_number)
        self.interpolator = self._generate_interpolator()

    def my_cell_factory(self, i, x):
        """ implement cell factory for single circular port """
        return CircularPortCell(i, x, self.r_init, self.r_ext)


class ConicCircularPort1D(Geometry1D):
    """
    ConicCircularPort1D inherits from SingleCircularPort1D and implements a conic port
    defined by its initial and final port radius.

        Attributes:
            1. r_init: initial radius of the port origin [m]
            2. r_final: radius at the end of the port [m]
    """

    def __init__(self, length, r_ext, nodes_number, regression_model, r_init, r_final, depth, **kwargs):
        """
        class initializer
        :param r_init: initial radius of the port 0 [m]
        :param r_final: radius at the end of the port [m]
        """

        # Call the parent class constructor
        super(ConicCircularPort1D, self).__init__(length, nodes_number, r_ext, regression_model)

        # Check the inputs
        assert r_ext > r_init > 0, "Origin internal radius (r_init) has to be greater than 0 and smaller than " \
                                   "External radius (r_ext).\n"
        assert r_ext > r_final > 0, "Final internal radius (r_init) has to be greater than 0 and smaller than " \
                                    "External radius (r_ext).\n"

        # Set the properties
        self.r_init, self.r_final, self.depth = r_init, r_final, depth
        for key, val in kwargs.items():
            setattr(self, key, val)
        self.mesh = self._generate_mesh(nodes_number)
        self.interpolator = self._generate_interpolator()

    def my_cell_factory(self, i, x):
        """ implement the cell factory for the class """
        return CircularPortCell(i, x, self.my_function(x), self.r_ext)

    def my_function(self, x):
        delta_h = self.r_final - self.r_init
        epsilon = 1e-4
        c = np.log(delta_h / epsilon - 1) / self.depth
        if x < self.length - self.exit_depth:
            r = self.r_init + delta_h / (1 + np.exp(-1 * c * (x - self.depth)))
        else:
            r_i = self.r_init + delta_h / (1 + np.exp(-1 * c * (self.length - self.exit_depth - self.depth)))
            r = self.r_exit - (self.r_exit - r_i) * (self.length - x) / self.exit_depth
        return r
