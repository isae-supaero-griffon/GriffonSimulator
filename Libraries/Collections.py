# Collections file contains the declaration of the Collection class
# used to reduce the amount of code declared when handling collections of different objects
# Author: Jose Felix Zapata Usandivaras
# 6/11/2018

# --------------------- IMPORT MODULES ----------------------

from collections import Iterable        # Import Iterable from the Collection Module

# --------------------- CLASSES DEFINITIONS -----------------


class Collections:
    """
    Collections class is an abstract class used to summarize the
    collections of different types of objects.
    Attributes:
        0. elements_list: list of the elements the Collection contains
    """

    def __init__(self, elements, *args):
        """
        class initializer
        :param elements: list of elements, element, or empty list
        :param args: Class of the element acceptable by the collection
        """
        # Set attributes
        self.elements_list = []  # Initialize the list
        if elements:
            if isinstance(elements, Iterable):
                for element in elements:
                    self.add_element(element, type(args[0]))
            else:
                self.add_element(elements, type(args[0]))

    def __str__(self):
        """ return a string with a str representation of the desired object """
        return "\n".join((" --> Element: {num} \n {value}".format(num=i, value=element)
                          for i, element in zip(range(0, len(self.elements_list)), self.elements_list)))

    def add_element(self, element, *args):
        """
        add a new element to the list of satellites
        :param element: element instance
        :param args: Class of the element acceptable by the collection
        :return: nothing
        """
        assert isinstance(element, args[0]), "Only {0} instances can be added to " \
                                             "SatelliteCollection".format(type(element).__name__)
        # Add the new element to the list
        self.elements_list.append(element)

    def is_empty(self):
        """
        is_empty to check if the collection is empty or not
        :return: boolean
        """
        return len(self.elements_list) == 0