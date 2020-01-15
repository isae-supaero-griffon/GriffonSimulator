# Components.py file includes the declaration of the different objects
# associated to the hydraulic network at a more concrete level. For instance,
# valves, injectors, tanks (PressurizerTank and OxidizerTank) among others.
# Each of this elements interacts with their nodes in a particular way. The physical
# description of this interaction is given by the objects implementation.
# @author: Jose Felix Zapata Usandivaras
# Date: 13/11/2019
# ISAE - SUPAERO, Space Section, Project "Griffon", Combustion Team

# -------------------------- IMPORT MODULES -------------------------------

from abc import ABC, abstractmethod                                                 # Import the abstract class mold
import numpy as np                                                                  # Import  numpy
from scipy.interpolate import interp2d                                              # Import interp2d from scipy
import math as m                                                                    # Import math library
from HydraulicModule.Fluid import *                                                 # Import the Fluids
from HydraulicModule.Network import *                                               # Import the Network
from Libraries.Reynolds import Reynolds                                             # Import the Reynolds class
from Libraries.Interpolator import Interpolator                                     # Import the Interpolator class

# ------------------------- FUNCTION DEFINITIONS ---------------------------

def liters_per_minute_2_m3_per_second(a):
    """
    liters_per_minute_2_m3_per_second converts from lt/min to m^3/sec
    :param a: float in lt/min
    :return: float in m^3/sec
    """
    return 1.66666667e-5 * a

def bar_2_pascals(a):
    """
    bar_2_pascals converts from bars to pascals
    :param a: float in bar
    :return: float in pascals
    """
    return 1e5 * a

def pascals_2_bar(a):
    """
    pascals_2_bar converts from pascals to bars.
    :param a: float in pascal
    :return: float in bars
    """
    return a/1e5

def psig_2_bar(a):
    """
    psig_2_bar converts from psig to bars
    :param a: float in psig
    :return: float in bars
    """
    return a*0.0689476

def circle_area(diameter):
    """
    circle_area determines the area of the circle  given its diameter
    :param diameter: diameter of the circle
    :return: area of the circle
    """
    return np.pi * diameter ** 2 / 4

def cv_2_kv(a):
    """
    cv_2_kv converts from the US units flow coefficient to the flow factor
    :param a: flow coefficient in US units
    :return: flow factor in SI units [m^3/h*sqrt(bar)]
    """
    return 0.865*a

# -------------------------- CLASS DEFINITIONS -----------------------------


class Component(ABC):
    """
    Component is the top class in the Hydraulic Components definition, except for Tanks.
    It is implemented as an abstract class.

        Attributes:
            1. name: string defining the name of the component
            2. number: integer which identifies the component ordering
            3. link: associated component of the hydraulic network
            4. is_active: boolean variable to indicate the activity of the component (by default is True)
    """

    def __init__(self, name, identifier, link):
        """
        class initializer
        :param name: string which defines the name of the component
        :param identifier: integer which defines the identifier of the component
        :param link: associated component of the hydraulic network
        """

        # Set the attritbute of component
        self.name = name
        self.number = identifier
        self.is_active = True
        self.pressure_nodes = None
        self.mass_node = None
        self.link = link

    @abstractmethod
    def _set_nodes(self, nodes):
        """
        _set_nodes is a private method that aims to aid in the setting of the nodes in the component
        :param nodes: list of nodes (already instantiated, associated to the component)
        :return: nothing
        """
        pass

    @abstractmethod
    def request_link(self):
        """
        request_link returns the mass dof through which the
        :return: mass node assocaited to the link
        """
        pass

    def request_link_temperature(self):
        """
        request_link_temperature returns the stored temperature of the link
        :return: link temperature, provided by its fluid
        """
        return self.link[0].fluid.get_temperature()

    @abstractmethod
    def calculate_delta_p(self):
        """
        calculate_delta_p determines the pressure loss across the component
        :return: float with delta P [Pa]
        """
        pass

    @abstractmethod
    def my_method(self, *args, **kwargs):
        """
        my_method is in charge of calculating the relationship between the dofs which are affected
        by the component.
        :param args: variable arguments expected by the method, time could be a factor.
        :param kwargs:
        :return: nothing.
        """
        pass

    @abstractmethod
    def update(self, dt):
        """
        update the object definition after one solver iteration if necessary (possible time step)
        :param dt: time-step
        :return: nothing
        """
        pass

    @abstractmethod
    def check_activity(self):
        """
        check_activity helps determine if the operating state of the component allows it
        to continue working or not
        :return: nothing
        """
        pass


class FlowPassage(Component):
    """
    FlowPassage class inherits from the Component class and aids in the implementation
    of the different components and differentiate from the Accumulator class which has
    a different behavior.

        Attributes:
            1. fluid: Fluid that goes through the component and whose properties will be used for calculation
            2. pressure_nodes: list of nodes, it can be from 1 to * depending on the component
            3. mass_nodes: node associated to the mass of the component
    """

    def __init__(self, name, identifier, link, fluid, nodes):
        """
        class initializer
        :param name: string which defines the name of the component
        :param identifier: integer which defines the identifier of the component
        :param link: associated component of the hydraulic network
        :param fluid: Fluid instance defining the fluid type that is present in the component
        :param nodes: list of nodes, it can be from 1 to * depending on the component
        """
        # Check the inputs
        assert isinstance(fluid, Fluid), "Input must be of Fluid type \n"
        assert len(nodes) == 3, "A FlowPassage object requires 3 nodes to work properly.\n"

        # Call superclass initializer
        super(FlowPassage, self).__init__(name, identifier, link)

        # Set the attributes / FlowPassage can only handle one fluid
        self.fluid = fluid

        # Set the nodes
        self._set_nodes(nodes)

    def __str__(self):
        nodes_text = "Nodes:: \n\t" + "\n\t".join([str(node) for node in self.pressure_nodes + [self.mass_node]])
        link_text = "Link:" + "\t".join([str(link.number) for link in self.link if link is not None])
        presentation_text = "\nComponent:: \n\t ID: {id} \n\t Name: {name} \n\t Fluid: {fluid} \n\t isActive: {act} " \
                            "\n\t {link} \n\t {nod}".format(id=self.number,
                                                            name=self.name,
                                                            fluid=self.fluid.name,
                                                            act=self.is_active,
                                                            link=link_text,
                                                            nod=nodes_text)
        return presentation_text

    def _set_nodes(self, nodes):
        self.pressure_nodes = []        # Define the pressure_nodes as a list
        for node in nodes:
            if node.dof.type == "pressure":
                self.pressure_nodes.append(node)
            elif node.dof.type == "mass":
                self.mass_node = node
            else:
                raise ValueError("Unrecognised node type {0}, check data file \n".format(node.type))

    def request_link(self):
        """ return the mass node dof """
        return self.mass_node.get_dof()

    def determine_mean_density(self):
        """
        determine_mean_density gives the average density of the fluid as a function of the pressure
        :return: average density of the fluid given the pressures
        """
        # First compute the density at the components mean pressure
        mean_pressure = (self.pressure_nodes[0].dof.get_value() +
                         self.pressure_nodes[1].dof.get_value()) / 2     # Compute the mean pressure at the component
        # Compute the density of the fluid
        return self.fluid.compute_density(mean_pressure, self.request_link_temperature())

    def my_method(self, *args, **kwargs):
        residual = self.pressure_nodes[0].dof.get_value() - self.pressure_nodes[1].dof.get_value() - \
                   self.calculate_delta_p()
        return np.array([pascals_2_bar(residual)])

    def update(self, dt):
        pass

    def check_activity(self):
        pass


class Accumulator(Component):
    """
    Accumulator class is an abstract class which inherits from Components and
    is designed to mimic the behavior of an accumulator by considering two mass nodes.
    It must be noted that an accumulator may be single phase (with a single fluid)
    or multiphase.

        Attributes:
            1. pressure_nodes: there are 3 pressure nodes in an accumulator
            2. mass_node: there are 2 mass node in an accumulator

    """
    def __init__(self, name, identifier, link, nodes):
        """
        class initializer
        :param name: string which defines the name of the component
        :param identifier: integer which defines the identifier of the component
        :param link: associated component of the hydraulic network
        :param nodes: list of nodes, it can be from 1 to * depending on the component
        """
        # Check the inputs
        assert len(nodes) == 5, "An Accumulator Requires 5 nodes to work. \n"
        # Call superclass initializer
        super(Accumulator, self).__init__(name, identifier, link)

        # Set the attributes
        self._set_nodes(nodes)

    def __str__(self):
        nodes_text = "Nodes:: \n\t" + "\n\t".join([str(node) for node in self.pressure_nodes + self.mass_node])
        link_text = "Link:" + "\t".join([str(link.number) for link in self.link if link is not None])
        presentation_text = "\nComponent:: \n\t ID: {id} \n\t Name: {name} \n\t isActive: {act} " \
                            "\n\t {link} \n\t {nod}".format(id=self.number,
                                                            name=self.name,
                                                            act=self.is_active,
                                                            link=link_text,
                                                            nod=nodes_text)
        return presentation_text

    def _set_nodes(self, nodes):
        self.pressure_nodes = []        # Set the pressure_nodes list
        self.mass_node = []             # Set the mass_node list
        # Allocate the nodes
        for node in nodes:
            if node.dof.type == "pressure":
                self.pressure_nodes.append(node)
            elif node.dof.type == "mass":
                self.mass_node.append(node)
            else:
                raise ValueError("Unrecognised node type {0}, check data file \n".format(node.type))

    def request_link(self):
        """ return the downstream mass dof """
        return self.mass_node[1].get_dof()

    @abstractmethod
    def calculate_delta_p_inlet(self):
        """
        calculate_delta_p_inlet determines the delta p at the inlet of the accumulator
        :return: delta_p at the inlet [Pa]
        """
        pass

    @abstractmethod
    def calculate_delta_p_outlet(self):
        """
        calculate_delta_p_outlet determines the delta p at the outlet of the accumulator
        :return: delta_p at the outlet [Pa]
        """
        pass

    def calculate_delta_p(self):
        return np.array([self.calculate_delta_p_inlet(),
                         self.calculate_delta_p_outlet()])

    def my_method(self, *args, **kwargs):
        delta_p_1 = self.pressure_nodes[0].dof.get_value() - self.pressure_nodes[1].dof.get_value()
        delta_p_2 = self.pressure_nodes[1].dof.get_value() - self.pressure_nodes[2].dof.get_value()
        vf = np.vectorize(pascals_2_bar)
        residual = np.subtract(np.array([delta_p_1, delta_p_2]), self.calculate_delta_p())
        return vf(residual)


class SinglePhaseTank(Accumulator):
    """
    SinglePhaseTank is a class which inherits from the abstract class Accumulator. It is in
    charge of implementing the behavior for a tank which a single phase behavior.

        Attributes:
            1. volume: tank volume in m^3
            2. fluid: Fluid class instance
            3. flow_factors: pressure loss at the constriction
            4. inlet_area: tank inlet area [m^2]
            5. outlet_area: tank outlet area [m^2]
            6. fluid_mass: mass of fluid contained in the tank
    """
    def __init__(self, name, identifier, link, nodes, fluid, volume, flow_factors, inlet_diameter, outlet_diameter,
                 initial_mass):
        """
        class initializer
        :param name: string which defines the name of the component
        :param identifier: integer which defines the identifier of the component
        :param link: associated component of the hydraulic network
        :param nodes: list of nodes, it can be from 1 to * depending on the component
        :param fluid: Fluid class instance
        :param volume: volume of the tank in m^3
        :param flow_factors: flow factors associated to the flow
        :param inlet_diameter: tank inlet diameter [m]
        :param outlet_diameter: tank outlet diameter [m]
        :param initial_mass: initial mass of fluid [kg]
        """
        # Call superclass initializer
        super(SinglePhaseTank, self).__init__(name, identifier, link, nodes)

        # Set the remaining attributes / SinglePhaseTank can only handle one fluid
        self.fluid = fluid
        self.volume = volume
        self.flow_factors = flow_factors
        self.inlet_area = circle_area(inlet_diameter)
        self.outlet_area = circle_area(outlet_diameter)
        self.fluid_mass = initial_mass

    def set_volume(self, new_volume):
        self.volume = new_volume

    def get_volume(self):
        return self.volume

    def determine_density_at_inlet(self):
        mean_pressure = (self.pressure_nodes[0].dof.get_value() +
                         self.pressure_nodes[1].dof.get_value()) / 2       # Compute the mean pressure at the component
        # Compute the density of the fluid
        return self.fluid.compute_density(mean_pressure, self.request_link_temperature())

    def determine_density_at_outlet(self):
        mean_pressure = (self.pressure_nodes[1].dof.get_value() +
                         self.pressure_nodes[2].dof.get_value()) / 2       # Compute the mean pressure at the component
        # Compute the density of the fluid with own temperature
        return self.fluid.compute_density(mean_pressure, self.fluid.get_temperature())

    def calculate_delta_p_inlet(self):
        rho_ = self.determine_density_at_inlet()
        return 0.5 * self.flow_factors[0] * (self.mass_node[0].dof.get_value() / self.inlet_area) ** 2 / rho_

    def calculate_delta_p_outlet(self):
        rho_ = self.determine_density_at_inlet()
        return 0.5 * self.flow_factors[1] * (self.mass_node[1].dof.get_value() / self.outlet_area) ** 2 / rho_

    def update(self, dt):
        """update the mass of the system"""
        self.fluid_mass += (self.mass_node[0].dof.get_value() - self.mass_node[1].dof.get_value())*dt

    def check_activity(self):
        if self.fluid_mass < 1e-6: self.is_active = False


class LiquidTank(SinglePhaseTank):
    """
    LiquidTank is an abstract class which is implemented to add functionality to
    SinglePhaseTank but it is never to be instantiated (hence the abstract class definition)
    since a LiquidTank cannot exist by itself due to being a DualPhaseTank. The only possibility of
    it existing is when it is filled.
    """
    def __init__(self, name, identifier, link, nodes, fluid, volume, flow_factors, inlet_diameter, outlet_diameter,
                 filling):
        """
        class initializer
        :param name: string which defines the name of the component
        :param identifier: integer which defines the identifier of the component
        :param link: associated component of the hydraulic network
        :param nodes: list of nodes, it can be from 1 to * depending on the component
        :param fluid: Fluid class instance
        :param volume: volume of the tank in m^3
        :param flow_factors: flow factors associated to the flow
        :param inlet_diameter: tank inlet diameter [m]
        :param outlet_diameter: tank outlet diameter [m]
        :param filling: filling percentage of the tank volume
        """
        # Check the inputs
        assert isinstance(fluid, Liquid), "For a LiquidTank instance the fluid provided must be a Liquid instance.\n"
        # Call superclass initializer
        initial_mass = filling * fluid.get_density() * volume
        super(LiquidTank, self).__init__(name, identifier, link, nodes, fluid, volume, flow_factors,
                                         inlet_diameter, outlet_diameter, initial_mass)

        # Finish initialization
        self.check_activity()                      # Check the activity of the component to the right value

    def get_empty_volume(self):
        """
        get_empty_volume returns the empty volume of the tank
        :return: empty volume of the tank
        """
        return self.volume - self.fluid_mass / self.fluid.get_density()


class GasTank(SinglePhaseTank):
    """
    PressurizerTank is a class that implements the behavior of a gas tank in a fully SinglePhaseTank

        Attributes:
            1. initial_pressure: initial_pressure of the tank [Pa]
            2. initial_temperature: initial temperature of the tank [K]
            3. discharge_constant: polytropic discharge constant value
            4. gas_mass: mass of gas contained in the tank [kg]
            5. discharge_constant: polytropic discharge constant value
    """

    def __init__(self, name, identifier, link, nodes, fluid, volume, flow_factors, inlet_diameter, outlet_diameter,
                 initial_pressure):
        """
        class initializer
        :param name: string which defines the name of the component
        :param identifier: integer which defines the identifier of the component
        :param link: associated component of the hydraulic network
        :param nodes: list of nodes, it can be from 1 to * depending on the component
        :param fluid: Fluid class instance
        :param volume: volume of the tank in m^3
        :param flow_factors: flow factors associated to the flow
        :param inlet_diameter: tank inlet diameter [m]
        :param outlet_diameter: tank outlet diameter [m]
        :param initial_pressure: initial pressure of the PressurizerTank in [Pa]
        """
        # Check the input
        assert isinstance(fluid, Gas), "The Fluid instance must be of Gas type.\n"
        # Call superclass initializer
        super(GasTank, self).__init__(name, identifier, link, nodes, fluid, volume, flow_factors, inlet_diameter,
                                      outlet_diameter, None)

        # Set the attributes
        self.initial_pressure = initial_pressure
        self.discharge_constant = None

        # Initialize the variables
        self._initialize_pressure()                 # Set the initial pressure at the tank intermediate node
        self._initialize_fluid()                    # Compute the density of the gas
        self._initialize_mass()                     # Set the initial mass of gas at the pressurizer tank
        self.set_discharge_constant()               # Set the discharge constant
        self.check_activity()                      # Check the activity of the component to the right value

    def _compute_local_gas_density(self):
        """
        _compute_local_gas_density updates the density of the fluid managed by this tank
        :return: nothing
        """
        self.fluid.set_density(self.fluid_mass / self.volume)

    def _compute_local_gas_temperature(self):
        """
        _compute_local_gas_temperature computes the local temperature of the gas given the current
        mass contained in the system
        :return: nothing
        """
        self.fluid.compute_local_gas_temperature(self.pressure_nodes[1].dof.get_value())

    def _compute_local_gas_pressure(self):
        """
        _compute_local_gas_pressure computes the local pressure of the gas given the current mass
        contained in the system
        :return: nothing
        """
        self.pressure_nodes[1].dof.set_value(self.discharge_constant * self.fluid.get_density() ** self.fluid.gamma)

    def set_discharge_constant(self):
        """ set_discharge_constant computes the value of the polytropic expansion
        constant associated to the discharge """

        # Compute mass of gas present in tank at discharge initiation
        self.discharge_constant = self.initial_pressure * (1 / self.fluid.get_density()) ** self.fluid.gamma

    def _initialize_pressure(self):
        """ _initialize_dof sets the pressure of the dof to the right value during
        creation of the object. It is a private method"""
        self.pressure_nodes[1].dof.fix()
        self.pressure_nodes[1].dof.set_value(self.initial_pressure)

    def _initialize_fluid(self):
        """ _initialize_fluid sets the right density of the fluid object given the
        initial pressure and temperature """
        rho_ = self.fluid.compute_density(self.initial_pressure, self.fluid.get_temperature())
        self.fluid.set_density(rho_)

    def _initialize_mass(self):
        """
        _initialize the mass value for the gas at the given pressure """
        self.fluid_mass = self.fluid.get_density() * self.volume

    def update(self, dt):
        super(GasTank, self).update(dt)                                        # Call superclass method
        self._compute_local_gas_density()                                      # Update the local gas density
        self._compute_local_gas_pressure()                                     # Compute the new local gas pressure
        self._compute_local_gas_temperature()                                  # Update the local gas temperature
        self.check_activity()                                                  # Set the component activity


class PressurizerTank(GasTank):
    """
    PressurizerTank is a particular implementation of the GasTank to account for the particular
    behavior of the PressurizerTank. It differences itself in the fact that it doesn't have an
    inlet so the GasTank is closed on one side. This is customaly forced in the code.
    """
    def __init__(self, name, identifier, link, nodes, fluid, volume, flow_factor, outlet_diameter,
                 initial_pressure):
        """
        class initializer
        :param name: string which defines the name of the component
        :param identifier: integer which defines the identifier of the component
        :param link: associated component of the hydraulic network
        :param nodes: list of nodes, it can be from 1 to * depending on the component
        :param fluid: Fluid class instance
        :param volume: volume of the tank in m^3
        :param flow_factor: flow factors associated to the outlet flow
        :param outlet_diameter: tank outlet diameter [m]
        :param initial_pressure: initial pressure of the PressurizerTank in [Pa]
        """
        # Check the inputs
        assert len(nodes) == 3, "Incorrect number of nodes in PressurizerTank.\n"

        # Create the fake nodes of the inlet (fixed value).
        mass_inlet, pressure_inlet = Node(0, Dof(0, "mass", 0)), Node(0, Dof(0, "pressure", initial_pressure))
        mass_inlet.dof.fix()
        pressure_inlet.dof.fix()

        # Append nodes
        nodes = [pressure_inlet, mass_inlet] + nodes

        # Call superclass initializer
        super(PressurizerTank, self).__init__(name, identifier, link, nodes, fluid, volume, [np.nan, flow_factor],
                                              np.nan, outlet_diameter, initial_pressure)

    def calculate_delta_p(self):
        return self.calculate_delta_p_outlet()

    def my_method(self, *args, **kwargs):
        residual = self.pressure_nodes[2].dof.get_value() - self.pressure_nodes[1].dof.get_value() -\
                   self.calculate_delta_p()
        return np.array([pascals_2_bar(residual)])


class DualPhaseTank(Accumulator):
    """
    OxidizerTank implements the behavior of the Oxidizer Tank by modelling the behavior between
    two phases (one gas and one liquid). Both phases are separated system which exchange volume
    and an interface which transmits pressure. The pressure inside the tank is determined by
    the gaseous phase.

        Attributes:
            1. phases_tanks: SinglePhaseTank implemented as objects with no relationship to the network
            2.
    """
    def __init__(self, name, identifier, link, nodes, fluid, volume, flow_factors, inlet_diameter, outlet_diameter,
                 initial_pressure, filling):
        """
        class initializer
        :param name: string which defines the name of the component
        :param identifier: integer which defines the identifier of the component
        :param link: associated component of the hydraulic network
        :param nodes: list of nodes, it can be from 1 to * depending on the component
        :param fluid: Fluid class instance
        :param volume: volume of the tank in m^3
        :param flow_factors: flow factors associated to the flow
        :param inlet_diameter: tank inlet diameter [m]
        :param outlet_diameter: tank outlet diameter [m]
        :param initial_pressure: initial pressure of the PressurizerTank in [Pa]
        :param filling: percentage filling of the oxidizer tank
        """
        # Check the inputs
        assert len(fluid)==2 , "Two different fluids are required.\n"

        # Call superclass initializer
        super(DualPhaseTank, self).__init__(name, identifier, link, nodes)

        # Set the attributes
        self.phases_tanks = self._create_internal_tanks(link, fluid, volume, flow_factors, inlet_diameter,
                                                        outlet_diameter, initial_pressure, filling)

    def _create_internal_tanks(self, link, fluid, volume, flow_factors, inlet_diameter, outlet_diameter,
                               initial_pressure, filling):
        """
        create_internal_tanks creates the tanks that are contained in the OxidizerTank object
        :return: internal tanks list
        """
        # Define the fake mass-node
        fake_node = Node(0, Dof(0, 'mass'))
        fake_node.dof.fix()

        # Create the first gas tank
        tank_1 = GasTank('FakeTank0',0, link[0], [self.pressure_nodes[0], self.pressure_nodes[1], self.pressure_nodes[1]]
                         + [self.mass_node[0], fake_node], fluid[0], (1-filling)*volume, [flow_factors[0], np.nan],
                         inlet_diameter, np.nan, initial_pressure)
        # Create the second liquid tank
        tank_2 = LiquidTank('FakeTank1', 1, link[1], [self.pressure_nodes[1], self.pressure_nodes[1],
                            self.pressure_nodes[2]] + [fake_node, self.mass_node[1]], fluid[1], volume,
                            [np.nan, flow_factors[1]], np.nan, outlet_diameter, filling)
        # Return the tanks
        return [tank_1, tank_2]

    def calculate_delta_p_inlet(self):
        return self.phases_tanks[0].calculate_delta_p_inlet()

    def calculate_delta_p_outlet(self):
        return self.phases_tanks[1].calculate_delta_p_outlet()

    def update(self, dt):
        # When updating the object the volume of the gas tank must be updated. At initialization
        # the pressure is already fixed.

        self.phases_tanks[1].update(dt)                                             # Update the liquid tank first
        self.phases_tanks[0].set_volume(self.phases_tanks[1].get_empty_volume())    # Determine the empty volume
        self.phases_tanks[0].update(dt)                                             # Update the gas tank

    def check_activity(self):
        #TODO: check this warning
        return self.phases_tanks[0].check_activity() and self.phases_tanks[1].check_activity()


class Valve(FlowPassage):
    """
    Valve class is an abstract class which serves as a parent class for Valves.

        Attributes:
            1. flow_factor: flow factor of the valve (constant)
    """

    def __init__(self, name, identifier, nodes, fluid, link, flow_factor):
        """
        class initializer
        :param name: string which defines the name of the component
        :param identifier: integer which defines the identifier of the component
        :param nodes: nodes associated to the valve
        :param fluid: Fluid instance defining the fluid type that is present in the component
        :param link: associated component of the hydraulic network
        :param flow_factor: flow factor of the valve (constant) in [m^3 h^-1 bar^-0.5]
        """
        # Call superclass initializer
        super(Valve, self).__init__(name, identifier, link, fluid, nodes)

        # Set other attributes
        self.flow_factor = cv_2_kv(flow_factor)


class LiquidValve(Valve):
    """
    LiquidValve implements the Delta P formula applied to a liquid valve
    """
    def __init__(self, name, identifier, nodes, fluid, link, flow_factor):
        # Check the input
        assert isinstance(fluid, Liquid), "The Fluid must be a Liquid. \n"
        # Call superclass initializer
        super(LiquidValve, self).__init__(name, identifier, nodes, fluid, link, flow_factor)

    def calculate_delta_p(self):
        rho_ = self.determine_mean_density()                        # Get the density of the fluid
        sp_g = self.fluid.sp                                        # Get the specific gravity
        volume_flow = 3600 *self.mass_node.dot.get_value() / rho_   # Get the volume flow in m^3/h
        return 1e5 * sp_g * (volume_flow / self.flow_factor) ** 2


class GasValve(Valve):
    """
    GasValve implements the Delta P formula applied to a Gas Valve considering the constants
    are defined for standard flow.
    """
    def __init__(self, name, identifier, nodes, fluid, link, flow_factor):
        # Check the input
        assert isinstance(fluid, Gas), "The Fluid must be a Gas. \n"
        # Call superclass initializer
        super(GasValve, self).__init__(name, identifier, nodes, fluid, link, flow_factor)

    def calculate_delta_p(self):
        volume_flow = 3600 *self.mass_node.dof.get_value() / self.fluid.std_density     # Get the standard flow [m^3/h]
        sp_g = self.fluid.sp                                                            # Get the specific gravity
        return 1e5 * sp_g * (volume_flow / self.flow_factor) ** 2


class Pipe(FlowPassage):
    """
    Pipe class implements a standard pipe with constants given

        Attributes:
            1. length: length of the pipe
            2. reynolds: Reynolds object instance
    """

    def __init__(self, name, identifier, nodes, fluid, link, length, diameter, roughness):
        """
        class initializer
        :param name: string which defines the name of the component
        :param identifier: integer which defines the identifier of the component
        :param nodes: nodes associated to the valve
        :param fluid: Fluid instance defining the fluid type that is present in the component
        :param link: associated component of the hydraulic network
        :param length: equivalent length of the valve
        :param diameter: diameter of piping
        :param density: density of the fluid
        :param roughness: relative roughness of the pipe
        """
        # Call superclass initializer
        super(Pipe, self).__init__(name, identifier, link, fluid, nodes)

        # Set other attributes
        self.length = length
        self.reynolds = Reynolds(roughness, diameter)

    def calculate_delta_p(self):
        rho_ = self.determine_mean_density()
        f = self.reynolds.solve_for_friction_factor(rho_, fluid.get_viscosity(), self.mass_node.dof.get_value())
        return f * self.length * self.mass_node.dof.get_value() ** 2 / \
               (2 * self.rho_ * self.reynolds.hydraulic_diameter * self.reynolds.area ** 2)


class Fitting(FlowPassage):
    """
    Fitting class implements a standard fitting with classic resistance coefficient

        Attributes:
            1. coefficient: resistance coefficient (non-dimensional)
    """

    def __init__(self, name, identifier, nodes, fluid, link, coefficient, diameter):
        """
        class initializer
        :param name: string which defines the name of the component
        :param identifier: integer which defines the identifier of the component
        :param nodes: nodes associated to the valve
        :param fluid: Fluid instance defining the fluid type that is present in the component
        :param link: associated component of the hydraulic network
        :param coefficient: resistance coefficient of the fitting
        :param diameter: diameter of piping
        """
        # Check the inputs
        assert coefficient > 0, "Resistance coefficient must be greater than 0. \n"
        assert diameter > 0, "Diameter of pipe must be greater than 0. \n"
        # Call superclass initializer
        super(Fitting, self).__init__(name, identifier, link, fluid, nodes)

        # Set other attributes
        self.coefficient = coefficient
        self.area = circle_area(diameter)

    def calculate_delta_p(self):
        rho_ = self.determine_mean_density()
        return 0.5 * self.coefficient * (self.mass_node.dof.get_value() / self.area) ** 2 / rho_


class Injector(FlowPassage):
    """
    Injector class implements a standard injector with tabulated data.

        # Attributes:
            1. Interpolator: interpolator object used to obtain the data.
    """

    def __init__(self, name, identifier, nodes, fluid, link, table):
        """
        class initializer
        :param name: string which defines the name of the component
        :param identifier: integer which defines the identifier of the component
        :param nodes: nodes associated to the valve
        :param fluid: Fluid instance defining the fluid type that is present in the component
        :param link: associated component of the hydraulic network
        :param table: json table to be implemented in the interpolator
        """
        # Check the fluid is correct
        assert isinstance(fluid, Liquid), "Fluid instance at injector must be of Liquid type. \n"
        # Call superclass initializer
        super(Injector, self).__init__(name, identifier, link, fluid, nodes)

        # Set the other attributes
        delta_ps, flows = map(np.array, zip(*table))
        flow_fun = np.vectorize(liters_per_minute_2_m3_per_second)
        dps_fun = np.vectorize(bar_2_pascals)
        self.interpolator = Interpolator(flow_fun(flows), dps_fun(delta_ps), kind='linear')

    def calculate_delta_p(self):
        rho_ = self.determine_mean_density()
        volume_flow = self.mass_node.dof.get_value() / rho_
        delta_p = self.interpolator.extrapolate(volume_flow)
        return delta_p


class PressureRegulator(FlowPassage):
    """
    PressureRegulator implements the behavior of the regulator of pressure which takes an
    input pressure and regulates it at the given setting pressure.
        Attributes:
            1. interpolator: 2d interpolator [x: inlet pressure [psig], y: N2_flow [slpm], z: outlet pressure [psig]
            2. threshold: minimum DeltaP above which the instrument can regulate the pressure to the output given
    """

    def __init__(self, name, identifier, nodes, fluid, link, table):
        """
        class initializer
        :param name: string which defines the name of the component
        :param identifier: integer which defines the identifier of the component
        :param nodes: nodes associated to the valve
        :param fluid: Fluid instance defining the fluid type that is present in the component
        :param link: associated component of the hydraulic network
        :param table: json table to be implemented in the interpolator
        """
        # Call superclass initializer
        super(PressureRegulator, self).__init__(name, identifier, link, fluid, nodes)

        # Set the other attributes
        curves = [np.append(record['pressure'] * np.ones((len(record['curve']), 1)),
                            np.array(record['curve']), axis=1) for record in table]
        curves = np.concatenate(curves, axis=0)
        f = np.vectorize(psig_2_bar)
        x, y, z = f(curves[:, 0]), curves[:, 1], f(curves[:, 2])

        # Set the remaining attributes
        self.interpolator = interp2d(x, y, z, kind='linear')

    @staticmethod
    def _recondition_data(record):
        """
        _recondition_data takes the record from the json file and outputs the x,y,z vectors
        :param record:
        :return:
        """

    def determine_mean_density(self):
        """ determine the mean density at the inlet """
        # First compute the density at the components mean pressure
        inlet_pressure = self.pressure_nodes[0].dof.get_value()                                 # Get the pressure at the inlet
        return self.fluid.compute_density(inlet_pressure, self.request_link_temperature())      # Compute the density of the fluid

    def calculate_delta_p(self):
        rho_ = self.determine_mean_density()                                                    # Calculate the density
        standard_rho = rho_ * 1.013 / pascals_2_bar(self.pressure_nodes[0].dof.get_value())     # Calculate std density
        standard_flow = self.mass_node.dof.get_value() / standard_rho                           # Calculate std flow
        inlet_pressure = self.pressure_nodes[0].dof.get_value()                                 # Extract internal p

        # Interpolate the delta P
        delta_p = inlet_pressure - bar_2_pascals(self.interpolator(pascals_2_bar(inlet_pressure), standard_flow))
        return delta_p


class ComponentsCatalogue:
    """ ComponentsCatalogue is a static class which serves as a
    reservoir of components constructors. """

    components = {'PressurizerTank': PressurizerTank,
                  'OxidizerTank': DualPhaseTank,
                  'PressureRegulator': PressureRegulator,
                  'LiquidValve': LiquidValve,
                  'GasValve': GasValve,
                  'Pipe': Pipe,
                  'Fitting': Fitting,
                  'Injector': Injector}

    @staticmethod
    def return_component(name):
        return ComponentsCatalogue.components[name]