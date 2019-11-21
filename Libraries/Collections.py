# Collections file contains the declaration of the Collection class
# used to reduce the amount of code declared when handling collections of different objects
# Author: Jose Felix Zapata Usandivaras
# 6/11/2018

# --------------------- IMPORT MODULES ----------------------

from collections import Iterable        # Import Iterable from the Collection Module
from abc import ABC, abstractmethod     # Import abc

# --------------------- CLASSES DEFINITIONS -----------------


class Collections(ABC):
    """
    Collections class is an abstract class used to summarize the
    collections of different types of objects.
    Attributes:
        0. elements_list: list of the elements the Collection contains
    """

    def __str__(self):
        """ return a string with a str representation of the desired object """
        return "Collection:: \n" + "\n".join(("{num} \t {value}".format(num=i, value=element)
                          for i, element in zip(range(0, len(self.elements_list)), self.elements_list)))

    @abstractmethod
    def add_element(self, element, *args):
        """
        add a new element to the list of satellites
        :param element: element instance
        :param args: Class of the element acceptable by the collection
        :return: nothing
        """
        pass

    def is_empty(self):
        """
        is_empty to check if the collection is empty or not
        :return: boolean
        """
        return len(self.elements_list) == 0

    def contains(self, filter):
        """
        contains searches from the elements_list to see if any object of matching criteria has
        already been instantiated
        :param filter: lambda function
        :return: boolean
        """
        for element in self.elements_list:
            if filter(element):
                return True
        return False

    @abstractmethod
    def sort(self):
        """ sort method is in charge of sorting the elements of the collection based on a given property
        of the instances. It is an abstract method of the class Collections """
        pass


class UniqueCollections(Collections):
    """ UniqueCollections is a variant of collections such that it accepts only new
    elements as defined by the contain criteria """

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
                    self.add_element(element, type(args[0]), args[1])
            else:
                self.add_element(elements, type(args[0]), args[1])

    def add_element(self, element, *args):
        assert isinstance(element, args[0]), "Only {0} instances can be added to " \
                                             "the Collection, tried {1} instead".format(args[0].__name__,
                                                                                        type(element).__name__)

        # Check if object is not present in list, then append, if negative, do nothing
        if not self.contains(filter=args[1]):
            self.elements_list.append(element)

    @abstractmethod
    def return_element(self, criteria):
        pass