# Network.py file contains the declarations of the classes associated to the network
# definition at a mathematical level, seen as a combination of linked nodes and degrees
# of freedom (Dof). Edges in the network and their corresponding behavior are detailed in
# the Components.py file.
# @author: Jose Felix Zapata Usandivaras
# Date: 13/11/2019
# ISAE - SUPAERO, Space Section, Project "Griffon", Combustion team.

# -------------------------- IMPORT MODULES -------------------------------

from abc import ABC, abstractmethod                                                 # Import the abstract class mold
import numpy as np                                                                  # Import  numpy
import math as m                                                                    # Import math library
from Libraries.Collections import Collections, UniqueCollections                    # Import the collections class

# -------------------------- CLASS DEFINITIONS -----------------------------


class Dof:
    """
    Dof is a class associated to the degree of freedom of the network.

        Attributes:
            1. number: integer indicating dof number
            2. type: string variable the dof is associated with
            2. value: float indicating the value of the degree of freedom
    """

    def __init__(self, number, type, scale=1, save=False, value=0):
        """
        class constructor
        :param number: degree of freedom number
        :param type: string variable the dof is associated with
        :param scale: scale to which responds the degree of freedom
        :param save: boolean which indicates if the dof should be saved or not
        :param value: float containing the value of the degree of freedom.
        By default it is initialized to 0.
        """
        # Check the inputs and set the attributes
        assert isinstance(number, int), "Dof number has to be an integer \n"
        assert scale > 0, "The scaling factor for the degrees of freedom must be greater than 0.\n"

        self.number = number
        self.type = type
        self.isFixed = False
        self._scale = scale
        self.value = value
        self.save = save

    def __str__(self):
        return "Dof:: Type: {type:10s}, ID: {id:5d}, isFixed: {fix:6b},Value: {val:>10.5f}".format(type=self.type,
                                                                                                   id=self.number,
                                                                                                   val=self.value,
                                                                                                   fix=self.isFixed)

    def set_value(self, value):
        self.value = value

    def set_scale(self, value):
        self._scale = value

    def get_scale(self):
        return self._scale

    def set_scaled_value(self, value):
        self.value = self._scale * value

    def get_scaled_value(self):
        return self.value / self._scale

    def get_value(self):
        return self.value

    def fix(self):
        self.isFixed = True


class DofCollection(UniqueCollections):
    """ Container of Dofs """
    def __init__(self, dofs):
        """
        class initializer
        :param dofs: dof object, iterable of dof, empty list
        """
        # Call superclass method
        super(DofCollection, self).__init__(dofs, Dof, lambda x: x.number == element.number)

    def add_element(self, element, *args):
        """ Override method from parent class"""
        super(DofCollection, self).add_element(element, Dof, lambda x: x.number == element.number)

    def sort(self):
        self.elements_list.sort(key=lambda x: x.number, reverse=False)

    def return_element(self, criteria):
        return next((elem for elem in self.elements_list if elem.number == criteria), None)

    def search(self, prop, value):
        """
        search an element by its property
        :param prop: attribute in string
        :param value: value of the attribute
        :return: list of found elements
        """
        return [elem for elem in self.elements_list if getattr(elem, prop) == value]

    def set_values(self, values_array):
        """
        set_values is a method that sets the values of the dofs collection based on a numpy values array
        :param values_array: numpy array
        :return: nothing
        """
        for elem, value in zip(self.elements_list, values_array):
            elem.set_value(value)


class Node:
    """
    Node is an abstract class in charge of storing the data at the nodes
    of the hydraulic circuit. That is, places in which we are interested in
    knowing a given degree of freedom.

        Attributes:
            1. number: integer indicating node number
            2. dof: degree of freedom of interest of the node
    """

    def __init__(self, number, dof):
        """
        class initializer
        :param number: integer indicating the node number of the network
        :param dof: degree of freedom
        """
        # Check the input
        assert isinstance(number, int), "Node number must be an integer \n"
        assert isinstance(dof, Dof), "Dof must be a Dof instance \n"

        # Set the attribute
        self.number = number
        self.dof = dof

    def __str__(self):
        return "Node:: ID: {id:5d} \t {dof}".format(id=self.number, dof=self.dof)

    def get_dof(self):
        """ get_dof returns the degree of freedom/s of the node """
        return self.dof

    def set_dof(self, dof):
        """ set_dof sets the degree of freedom/s of the node """
        self.dof = dof


class NodesCollection(UniqueCollections):
    """
    NodesCollection class inherits from collections and helps store
    nodes in the Component class. It serves as a container class.
    """

    def __init__(self, nodes):
        """
        class initializer
        :param nodes: node object, iterable of nodes, or empty list
        """
        # Call superclass method
        super(NodesCollection, self).__init__(nodes, Node, lambda x: x.number == element.number)

    def add_element(self, element, *args):
        """ Override method from parent class"""
        super(NodesCollection, self).add_element(element, Node, lambda x: x.number == element.number)

    def sort(self):
        self.elements_list.sort(key=lambda x: x.number, reverse=False)

    def return_element(self, criteria):
        return next((elem for elem in self.elements_list if elem.number == criteria), None)
