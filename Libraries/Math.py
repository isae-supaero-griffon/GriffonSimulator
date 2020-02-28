# Math.py library helps introduce mathematical functionality required.
# @author: Jose Felix Zapata Usandivaras
# Date: 24/01/2020
# ISAE - SUPAERO Space Section, Project Griffon, Combustion Team.

# -------------------------- IMPORT MODULES -------------------------------

import numpy as np


# ------------------------- FUNCTION DEFINITIONS ---------------------------

def lagrange(x, i, xm):
    """Evaluates  the i-th  Lagrange  polynomial  at xbased on grid  data xm"""
    n = len(xm) - 1
    y = 1.
    for j in range(n + 1):
        if i != j:
            y *= (x - xm[j]) / (xm[i] - xm[j])
    # Return the output
    return y

# -------------------------- CLASS DEFINITIONS -----------------------------


class LagrangeCoefficient:
    """ Represent a Lagrange Coefficient """

    def __init__(self, i, xm):
        """
        class initializer
        :param i: number of polynomial
        :param xm: vector of x values
        """
        # Set the attributes
        self.i = i
        self.xm = xm

    def value(self, x):
        """
        value calculates the values of the lagrange coefficient at point x
        :param x: point of interest
        :return: value
        """
        return lagrange(x, self.i, self.xm)


class RangedPolynomial:
    """ Represent the mathematical interpretation of a polynomial with range on x"""

    def __init__(self, coefs, xmax):
        """
        class initializer
        :param coefs: coefficients of the polynomial
        :param xmax: limit values
        """
        # Set the attributes
        self.pol = np.poly1d(coefs)
        self.x_max = xmax

    def value(self, x):
        """
        value calculates the values of the Ranged Polynomial if outside boundary returns 0
        :param x: point of interest
        :return: value
        """
        if x <= self.x_max:
            return self.pol(x)
        else:
            return 0
