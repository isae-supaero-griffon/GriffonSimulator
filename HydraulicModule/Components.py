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

# -------------------------- CLASS DEFINITIONS -----------------------------


class Component(ABC):
    """
    Component is the top class in the Hydraulic Components definition. It is
    implemented as an abstract class. Only nodes attribute is considered for the
    present moment, additional data is added in the subsequent classes.

        Attributes:
            1. name: string defining the name of the component
            2. number: integer which identifies the component ordering
            3. pressure_nodes: list of nodes, it can be from 1 to * depending on the component
            4. mass_nodes: node associated to the mass of the component
            4. fluid: Fluid that goes through the component and whose properties will be used for calculation
            5. link: associated component of the hydraulic network
            5. is_active: boolean variable to indicate the activity of the component (by default is True)
    """

    def __init__(self, name, identifier, nodes,  fluid, link):
        """
        class initializer
        :param name: string which defines the name of the component
        :param identifier: integer which defines the identifier of the component
        :param nodes: list of nodes, it can be from 1 to * depending on the component
        :param fluid: Fluid instance defining the fluid type that is present in the component
        :param link: associated component of the hydraulic network
        """
        # Check the inputs
        assert isinstance(fluid, Fluid), "Input must be of Fluid type \n"
        assert isinstance(link, Component) or link is None, "The link must be of Component type or None \n"

        # Set the attritbute of component
        self.name = name
        self.number = identifier
        self.fluid = fluid
        self.is_active = True
        self.pressure_nodes = []
        self.mass_node = None
        self.link = link

        # Set the nodes
        self._set_nodes(nodes)

    def __str__(self):
        nodes_text = "Nodes:: \n\t" + "\n\t".join([str(node) for node in self.pressure_nodes+ [self.mass_node]])
        link_text = "Link: {link} \n".format(link=self.link.number if self.link is not None else None)
        presentation_text = "\nComponent:: \n\t ID: {id} \n\t Name: {name} \n\t Fluid: {fluid} \n\t isActive: {act} " \
                            "\n\t {link} \n\t {nod}".format(id=self.number,
                                                            name=self.name,
                                                            fluid=self.fluid.name,
                                                            act=self.is_active,
                                                            link=link_text,
                                                            nod=nodes_text)
        return presentation_text

    def _set_nodes(self, nodes):
        """
        _set_nodes is a private method that aims to aid in the setting of the nodes in the component
        :param nodes: list of nodes (already instantiated, associated to the component)
        :return: nothing
        """
        for node in nodes:
            if node.dof.type == "pressure":
                self.pressure_nodes.append(node)
            elif node.dof.type == "mass":
                self.mass_node = node
            else:
                raise ValueError("Unrecognised node type {0}, check data file \n".format(node.type))

    @abstractmethod
    def calculate_delta_p(self):
        """
        calculate_delta_p determines the pressure loss across the component
        :return: float with delta P [Pa]
        """
        pass

    def my_method(self, *args, **kwargs):
        """
        my_method is in charge of calculating the relationship between the dofs which are affected
        by the component.
        :param args: variable arguments expected by the method, time could be a factor.
        :param kwargs:
        :return: nothing.
        """
        return self.pressure_nodes[0].dof.get_value() - self.calculate_delta_p() -\
               self.pressure_nodes[1].dof.get_value()

    @abstractmethod
    def update(self, dt):
        """
        update the object definition after one solver iteration if necessary (possible time step)
        :param dt: time-step
        :return: nothing
        """
        pass

    @abstractmethod
    def _check_activity(self):
        """
        _check_activity helps determine if the operating state of the component allows it
        to continue working or not
        :return: nothing
        """
        pass



class PressurizerTank(Component):
    """
    PressurizerTank is in charge of describing the physics of the pressurizer tank.
    Since it is a discharge tank it only accepts one node in his NodesCollection.

        Attributes:
            1. volume: internal volume of the tank
            2. gamma: ratio of specific heats
            3. initial_pressure: initial pressure of the tank [Pa]
            4. temperature: initial temperature of the tank [K]
            5. gas_constant: gas constant [J/KgK]
            6. gas_mass: initial mass of gas present
            7. discharge_constant: polytropic discharge constant value
            8. exit_area: exit area [m^2] of the tank open end
            8. cd: discharge coefficient: discharge coefficient associated to the flow of the gas
            out of the reservoir.
    """

    def __init__(self, name, identifier, nodes, fluid, link, volume, initial_pressure, temperature, exit_area, cd):
        """
        class initializer
        :param name: string which defines the name of the component
        :param identifier: integer which defines the identifier of the component
        :param nodes: nodes collection associated to the tank. Must be of length 1.
        :param link: associated component of the hydraulic network
        :param fluid:  Fluid instance defining the fluid type that is present in the component
        :param volume: internal volume of the tank [m^3]
        :param initial_pressure: initial setting pressure
        :param temperature: initial temperature of the gas [K]
        :param exit_area: exit area [m^2] of the tank open end
        :param cd: discharge coefficient of the tank
        """
        # Call superclass constructor
        super(PressurizerTank, self).__init__(name, identifier, nodes, fluid, link)

        # Set the remaining attributes
        self.volume = volume
        self.initial_pressure = initial_pressure
        self.temperature = temperature
        self.gas_mass = None
        self.discharge_constant = None
        self.exit_area = exit_area
        self.cd = cd

        # Initialize the variables
        self._initialize_pressure()                 # Set the initial pressure at the tank node
        self._initialize_fluid()                    # Compute the density of the gas
        self._initialize_mass()                     # Set the initial mass of gas at the pressurizer tank
        self.set_discharge_constant()               # Set the discharge constant
        self._check_activity()                      # Check the activity of the component to the right value

    def _compute_local_gas_density(self):
        """
        _compute_local_gas_density updates the density of the fluid managed by this tank
        :return: nothing
        """
        self.fluid.set_density(self.gas_mass / self.volume)

    def _compute_local_gas_temperature(self):
        """
        _compute_local_gas_temperature computes the local temperature of the gas given the current
        mass contained in the system
        :return: nothing
        """
        self.temperature = self.fluid.compute_local_gas_temperature(
            self.nodes.elements_list[0].dofs.elements_list[0].get_value())

    def _compute_local_gas_pressure(self):
        """
        _compute_local_gas_pressure computes the local pressure of the gas given the current mass
        contained in the system
        :return: nothing
        """
        self.pressure_nodes[0].dof.set_value(self.discharge_constant * self.fluid.get_density() ** self.fluid.gamma)

    def _check_activity(self):
        if self.gas_mass < 1e-6: self.is_active = False

    def update(self, dt):
        self.gas_mass -= self.mass_node.dof.get_value() * dt                 # Update mass contained in the reservoir
        self._compute_local_gas_density()                                    # Update the local gas density
        self._compute_local_gas_temperature()                                # Update the local gas temperature
        self._compute_local_gas_pressure()                                   # Compute the new local gas pressure
        self._check_activity()                                               # Set the component activity

    def _initialize_pressure(self):
        """ _initialize_dof sets the pressure of the dof to the right value during
        creation of the object. It is a private method"""
        self.pressure_nodes[0].dof.fix()
        self.pressure_nodes[0].dof.set_value(self.initial_pressure)

    def _initialize_fluid(self):
        """ _initialize_fluid sets the right density of the fluid object given the
        initial pressure and temperature """
        self.fluid.compute_density(self.initial_pressure, self.temperature)

    def _initialize_mass(self):
        """
        _initialize the mass value for the gas at the given pressure """
        self.gas_mass = self.fluid.get_density() * self.volume

    def set_discharge_constant(self):
        """ set_discharge_constant computes the value of the polytropic expansion
        constant associated to the discharge """

        # Compute mass of gas present in tank at discharge initiation
        self.discharge_constant = self.initial_pressure * (1 / self.fluid.get_density()) ** self.fluid.gamma

    def calculate_delta_p(self):
        # Comments: The delta P is calculated using the discharge factor of the tank
        # TODO: review if the tanks discharge factor variate with the flow or not
        return self.mass_node.dof.get_value() ** 2 / (2 * self.fluid*get_density() * self.cd ** 2 * self.exit_area ** 2)


class OxidizerTank(Component):
    """
    OxidizerTank class implements the behaviour of the oxidizer tank.
        Attributes:
            1. mass: mass of oxidizer [kg] left in the tank
            2. exit_area: discharge area of the tank [m^2]
            3. cd: discharge coefficient of the tank at its exit
    """

    def __init__(self, name, identifier, nodes, fluid, link, volume, filling, exit_area,cd):
        """
        class initializer
        :param name: string which defines the name of the component
        :param identifier: integer which defines the identifier of the component
        :param nodes: nodes associated to the hydraulic module
        :param fluid: Fluid instance defining the fluid type that is present in the component
        :param link: associated component of the hydraulic network
        :param volume: internal volume of the tank [m^3]
        :param filling: percentage filling of the tank
        :param exit_area: discharge area of the tank [m^2]
        :param cd: discharge coefficient
        """
        # Call super class constructor
        super(OxidizerTank, self).__init__(name, identifier, nodes, fluid, link)

        # Set remaining attributes
        self.mass = volume * filling * self.fluid.get_density() / 1000
        self.exit_area = exit_area
        self.cd = cd

        # Set the activity
        self._check_activity()

    def update(self, dt):
        self.mass -= self.mass_node.dof.get_value() * dt
        self._check_activity()

    def _check_activity(self):
        if self.mass < 1e-6: self.is_active = False

    def calculate_delta_p(self):
        # Comments: the delta P is based on the discharge coefficient of the tank
        # TODO: check if the discharge coefficient of the tank varies with the another parameter
        return self.mass_node.dof.get_value() ** 2 / (2 * self.fluid.get_density() * self.cd ** 2 * self.exit_area ** 2)


class Valve(Component):
    """
    Valve class implements a standard valve with constants given.

        Attributes:
            1. length: equivalent length of valve [m]
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
        super(Valve, self).__init__(name, identifier, nodes, fluid, link)

        # Set other attributes
        self.length = length
        self.reynolds = Reynolds(roughness, diameter)

    def calculate_delta_p(self):
        f = self.reynolds.solve_for_friction_factor(self.fluid, self.mass_node.dof.get_value())
        return f * self.length * self.mass_node.dof.get_value() ** 2 / \
               (2 * self.fluid.get_density() * self.reynolds.hydraulic_diameter * self.reynolds.area ** 2)

    def update(self, dt):
        pass

    def _check_activity(self):
        pass


class Injector(Component):
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
        # Call superclass initializer
        super(Injector, self).__init__(name, identifier, nodes, fluid, link)

        # Set the other attributes
        delta_ps, flows = map(np.array, zip(*table))
        flow_fun = np.vectorize(liters_per_minute_2_m3_per_second)
        dps_fun = np.vectorize(bar_2_pascals)
        self.interpolator = Interpolator(flow_fun(flows), dps_fun(delta_ps), kind='linear')

    def calculate_delta_p(self):
        # Comment: calculate the delta p from extrapolating from the curve
        volume_flow = self.mass_node.dof.get_value() / self.fluid.get_density()
        delta_p = self.interpolator.extrapolate(volume_flow)
        return delta_p

    def _check_activity(self):
        pass

    def update(self, dt):
        pass


class PressureRegulator(Component):
    """
    PressureRegulator implements the behavior of the regulator of pressure which takes an
    input pressure and regulates it at the given setting pressure.
        Attributes:
            1. output_pressure: setting pressure at which the mechanism operates
            2. threshold: minimum DeltaP above which the instrument can regulate the pressure to the output given
    """

    def __init__(self, name, identifier, nodes, fluid, link, set_point, threshold):
        """
        class initializer
        :param name: string which defines the name of the component
        :param identifier: integer which defines the identifier of the component
        :param nodes: nodes associated to the valve
        :param fluid: Fluid instance defining the fluid type that is present in the component
        :param link: associated component of the hydraulic network
        :param set_point: set_point at which the pressure regulator is set
        :param threshold: minimum deltaP required for the device to work
        """
        # Call superclass initializer
        super(PressureRegulator, self).__init__(name, identifier, nodes, fluid, link)
        assert set_point > 0, "set_point of Pressure Regulator must be positive value \n"


        # Set remaining attributes
        self.set_point = set_point
        self.threshold = threshold

        # Initialize the output pressure
        self._initialize_output_pressure()

    def _initialize_output_pressure(self):
        """
        _initialize_output_pressure sets the value of the output pressure given the set_point provided
        :return: nothing
        """
        self.pressure_nodes[1].dof.fix()
        self.pressure_nodes[1].dof.set_value(self.set_point)

    def _check_activity(self):
        if self.pressure_nodes[0].dof.get_value() - self.pressure_nodes[1].dof.get_value() < self.threshold:
            self.is_active = False

    def calculate_delta_p(self):
        return self.pressure_nodes[0].dof.get_value() - self.pressure_nodes[1].dof.get_value()

    def my_method(self, *args, **kwargs):
        # Comment: calculate the mass flow associated to the equipment
        pass

    def update(self, dt):
        self._check_activity()


class ComponentsCatalogue:
    """ ComponentsCatalogue is a static class which serves as a
    reservoir of components constructors. """

    components = {'PressurizerTank': PressurizerTank,
                  'OxidizerTank': OxidizerTank,
                  'PressureRegulator': PressureRegulator,
                  'Valve': Valve,
                  'Injector': Injector}

    @staticmethod
    def return_component(name):
        return ComponentsCatalogue.components[name]