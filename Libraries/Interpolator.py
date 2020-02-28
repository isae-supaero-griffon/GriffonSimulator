# Interpolator.py implements the Interpolator class which is consumed in multiple parts of the
# code that work with tabulated data.
# @author: Jose Felix Zapata Usandivaras
# Date: 21/11/2019
# ISAE - SUPAERO, Space Section Project "Griffon", Combustion Team

# ------------------------- IMPORT MODULES --------------------------- #

import numpy as np                                      # Import numpy
from scipy.interpolate import interp1d                  # Import interp1d

# ------------------------- CLASS DEFINITIONS ------------------------ #


class Interpolator:
    """
    Interpolator class helps in the interpolation of the cross section
    area for the single port.

        Attrs:
            0. x_cor: coordinates of x
            1. y_cor: coordinates of y (cross section area of port)
            2. kind: interpolator kind
    """

    def __init__(self, x_cor=0, y_cor=0, kind='linear'):
        """
        class initializer
        :param x_cor: initial values for the x_cor
        :param y_cor: initial values for the y_cor
        """

        # Set the attributes
        self.x_cor = x_cor
        self.y_cor = y_cor
        self.kind = kind
        self.interp_fun = interp1d(self.x_cor, self.y_cor, kind=self.kind)
        self.extrap_fun = interp1d(self.x_cor, self.y_cor, kind=self.kind, bounds_error=False, fill_value='extrapolate')

    def set_coordinates(self, x_cor, y_cor):
        self.x_cor, self.y_cor = x_cor, y_cor
        self.update_interpolator()

    def update_interpolator(self):
        self.interp_fun = interp1d(self.x_cor, self.y_cor, kind=self.kind)
        self.extrap_fun = interp1d(self.x_cor, self.y_cor, kind=self.kind, bounds_error=False, fill_value='extrapolate')

    def interpolate(self, x):
        """
        perform bounded interpolation for the value of x
        :param x: coordinate x
        :return: interpolated value
        """
        # return value
        f_x = interp1d(self.x_cor, self.y_cor, bounds_error=True)
        return f_x(x)

    def extrapolate(self, x):
        """
        perform unbounded interpolation
        :param x: coordinate x
        :return: interpolated value
        """
        # return value
        return self.extrap_fun(x)
