# Hydraulic.py file contains the Hydraulic Module objects definitions.
# These includes the definition of the different hydraulic components included
# in the definition of the rockets hydraulic circuit, such as pressurizer tank,
# pressure regulator, oxidizer tank, valves, etc.
# The problem is considered as pseudo-stationary, thus local components equations
# are considered in their stationary form.
# @author: Jose Felix Zapata Usandivaras
# Date: 28/10/2019
# ISAE - SUPAERO Space Section, Project "Griffon"
#


# -------------------------- IMPORT MODULES -------------------------------

from abc import ABC, abstractmethod                                                 # Import the abstract class mold
import numpy as np                                                                  # Import  numpy
import math as m                                                                    # Import math library
from Libraries.Collections import Collections, UniqueCollections                    # Import the collections class
from HydraulicModule.Fluid import *                                                 # Import the different fluids
from HydraulicModule.Network import *                                               # Import the Network classes
from HydraulicModule.Components import *                                            # Import the Hydraulic Components
import scipy.interpolate as interp                                                  # Import the interpolate module
import scipy.optimize as opt                                                        # Import the optimization library

# ------------------------- FUNCTION DEFINITIONS ---------------------------.


def compute_error(old_vector, new_vector):
    """
    compute_error calculates the relative error on the norm of two different vectors.
    :param old_vector: old numpy array
    :param new_vector: new numpy array (same size of old vector)
    :return: relative error
    """
    old_norm = np.linalg.norm(old_vector)
    if old_norm > 1e-12:
        return np.linalg.norm(new_vector - old_vector) / old_norm
    else:
        raise ZeroDivisionError("Norm of vector already 0 \n")


def get_x_activity(x):
    return x.is_active

# -------------------------- CLASS DEFINITIONS -----------------------------


class HydraulicModule:
    """
    HydraulicModule class is the master class in charge of commanding the hydraulic
    module physics. Its functions include, creation of the different components, initialization
    and later commanding the network solving.

        Attributes:
            1. fluid: list of fluids present on the hydraulic module
            2. components: list of components that integrate the hydraulic network
            3. nodes: list of nodes that integrate the hydraulic network
            4. dofs: list of dofs that integrate the hydraulic network
    """

    # Define useful functions to be generated once during the simulation
    activity_function = np.vectorize(get_x_activity)

    class RunEndingError(Exception):
        """
        RunEndingException is a nested class of the HydraulicModule. It shall be explicitly used
        to indicate that a run has finished do to the cease in activity of one of the components
        """
        pass

    def __init__(self, hydraulic_table):
        """
        class initializer
        :param hydraulic_table: dictionary defining the data necessary for the hydraulic module
        """
        # Initialize attributes as lists
        self._dimension = 0
        self.fluids = {}
        self.components = []
        self.nodes = NodesCollection([])
        self.dofs = DofCollection([])
        self.is_initialized = False
        self.is_active = True
        self.solver_params = {'name': 'Nelder-Mead', 'tol': 1e-3, 'maxiter': 100, 'scale': 1e5}

        # Initialize the network
        self._initialize(hydraulic_table)

    def __str__(self):
        presentation_text = "Hydraulic Module:: \n " + self.solver_print() + "\n".join((str(obj)
                                                                                        for obj in self.components))
        return presentation_text

    def solver_print(self):
        return "\n Solver:: " + ", ".join(("{0}: {1}".format(key, value) for key, value in self.solver_params.items()))

    def dof_print(self):
        return "Hydraulic Module:: \n" + str(self.dofs)

    def node_print(self):
        return "Hydraulic Module:: \n" + str(self.nodes)

    def update(self, dt):
        """
        update methods runs the update method on the components that make the hydraulic module
        :param dt: time-step [sec]
        :return: nothing
        """
        for component in self.components:
            component.update(dt)

    def checkout_component(self, name):
        """
        return_component returns the component specified by the name inputted
        :param name: string with the name of the component
        :return: component instance
        """
        return next((elem for elem in self.components if elem.name == name), None)

    def _initialize(self, hydraulic_table):
        """
        _initialize is a private method of he HydraulicModule class
        :param hydraulic_table: dictionary defining the data necessary for the hydraulic module
        :return: nothing
        """
        # Generate the fluids
        self.fluids = {r['name']: r for r in hydraulic_table['fluids']}
        # Instantiate the components
        for component_dict in hydraulic_table['components']:
            self._create_component(component_dict)
        # Set the solver params
        if 'solver' in hydraulic_table:
            self.solver_params = hydraulic_table['solver']
        # Check activity
        self.check_activity()

    def _create_dof(self, number, type):
        """
        _create_dof instantiates a new Dof object and appends it to the dof list
        :param number: dof identifier
        :param type: objects variable type
        :return: Dof instance of interest
        """

        # Define the scale for the degree of freedom
        if type == "pressure":
            s = self.solver_params['scale']
        else:
            s = 1
        self.dofs.add_element(Dof(number, type, s))
        return self.dofs.return_element(number)

    def _create_node(self, node_dict, link_iterator=None):
        """
        _create_node instantiates a new Node object and appends it to the nodes list
        :param node_dict: dictionary defining the node keys: ["id", "type"}
        :param link_iterator: iterator of components used for instantiation of the mass dofs
        :return: Node instance
        """

        # Extract the node data
        node_number, dof_type = node_dict['identifier'], node_dict['type']

        # Create associated Dofs and DofCollection
        if dof_type == "mass":
            link = next(link_iterator)
            my_dof = link.request_link() if link is not None else self._create_dof(node_number, dof_type)
        else:
            my_dof = self._create_dof(node_number, dof_type)

        # Instantiate the node and add it to the NodesCollection
        self.nodes.add_element(Node(node_number, my_dof))
        return self.nodes.return_element(node_number)

    def _create_fluid(self, name, T):
        """
        _create_fluid is a private method who's in charge of generating the fluid
        :param name: name of the fluid of interest
        :param`T: temperature of initialization
        :return: Fluid instance
        """
        # Check the inputs
        assert name in self.fluids, "Fluid requested: {0} not present in fluids list. \n".format(name)
        assert T > 0, "Temperature declared is in Kelvin \n"

        # Extract the fluid dictionary
        my_dict = self.fluids[name]                                     # Extract the fluid dictionary
        my_init = FluidCatalogue.return_initializer(my_dict['fluid'])   # Extract the initializer from the catalogue
        my_obj = my_init(**my_dict)                                     # Create the object
        my_obj.set_temperature(T)                                       # Set the temperature of the fluid

        # Return the object
        return my_obj

    def return_component(self, number):
        """
        return_component returns the requested component according to its number
        :param number: identifier of the component
        :return: Component instance
        """
        return next((elem for elem in self.components if elem.number == number), None)

    def _create_component(self, component_dictionary):
        """
        _create_component instantiates the components
        :param component_dictionary: dictionary associated to the component definition
        :return: nothing
        """
        my_dict = dict(component_dictionary)                                            # Make deep copy of dict
        my_dict['link'] = [self.return_component(value) for value in my_dict['link']]   # Replace the link by the object
        my_iter = iter(my_dict['link'])
        my_dict['nodes'] = [self._create_node(node_dict, my_iter)
                            for node_dict in my_dict['nodes']]                          # Create the nodes of the system
        my_dict['fluid'] = [self._create_fluid(**r) for r in my_dict['fluid']]          # Store the fluid
        if len(my_dict['fluid']) == 1:
            my_dict['fluid'] = my_dict['fluid'][0]

        # Create the component
        my_init = ComponentsCatalogue.return_component(my_dict.pop('type'))
        new_component = my_init(**my_dict)
        new_component.set_queue(range(self._dimension, self._dimension + new_component.dimension))
        self.components.append(new_component)
        self._dimension += new_component.dimension

    def set_chamber_pressure(self, value):
        """
        set_chamber_pressure sets the value of the chamber pressure at the final node in the network
        :param value: value of the chamber pressure
        :return: nothing
        """
        self.dofs.elements_list[-1].fix()
        self.dofs.elements_list[-1].set_value(value)

    def initialize_pressure_dofs(self):
        """
        initialize_pressure_dofs initializes through linear interpolation the
        values of the non-fixed pressure dofs
        :return:
        """
        # Get pressure dofs and interpolate the values from the fixed ones
        pressure_dofs = self.dofs.search(prop="type", value="pressure")
        x, y = zip(*[(dof.number, dof.get_value()) for dof in pressure_dofs if dof.isFixed])

        # Check that there are fixed pressure dofs
        if not x:
            raise ValueError("No Fixed pressure nodes are present \n")

        # Interpolate the values
        f_interp = interp.interp1d(x, y, bounds_error=False, fill_value='extrapolate')
        for dof in pressure_dofs:
            dof.set_value(f_interp(dof.number))

        # Set isInitialized
        self.is_initialized = True

    def check_activity(self):
        """
        check_activity determines is the network is still operating, if one element of the
        network is down it shuts-down the object
        :return: nothing
        """
        self.is_active = bool(np.prod(self.activity_function(self.components)))

    def return_exit_flow(self):
        """
        return_exit_flow returns the flow at the end of the hydraulic module. Typically an injector
        :return: exit flow in [kg/s]
        """
        return self.components[-1].mass_dof.dof.get_value()

    def run_simulation(self):
        """
        run_simulation method solves for the current state of the hydraulic network provided that boundary
        conditions have already been set.
        :return: nothing
        """

        # Search dofs that are not fixed
        scale = self.solver_params['scale']
        non_fixed_dofs = self.dofs.search("isFixed", False)
        res = np.empty(shape=(self._dimension, 1))

        def residual(x):
            """
            residual function implements the function which will be later passed into the solver
            :param x: non-fixed dof values in iterable format according to the right order
            :return: residual value
            """
            # Set their value - only for pressure dofs
            for value, dof in zip(x, non_fixed_dofs):
                dof.set_scaled_value(value)

            # Run the methods from the component
            for comp in self.components:
                comp.my_method(res)

            # Get the residual
            return res

        def optimize_fun(x):
            """
            optimize_fun defines the optimization criteria for the function to be minimized
            :param x: on-fixed dof values in iterable format according to the right order
            :return: objective function (norm of residual ---> solution tends to 0)
            """
            return np.linalg.norm(residual(x))/scale

        # Define the optimization bounds
        # TODO: implement a solver for root finding, investigate if its better than optimize - bounds definition
        x0_fun = np.vectorize(lambda x: x.get_value()/scale if x.type == "pressure" else x.get_value())
        x0 = x0_fun(non_fixed_dofs)

        # Check the activity of the system
        self.check_activity()
        if not self.is_active:
            raise self.RunEndingError("\n At least one of the components is not active \n")

        # Solve the problem:
        sol = opt.minimize(optimize_fun, x0, method=self.solver_params['name'], jac=None,
                           tol=self.solver_params['tol'],
                           options=self.solver_params['options'])
        return sol, res
