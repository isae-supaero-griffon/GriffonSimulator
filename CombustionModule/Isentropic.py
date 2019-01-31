# Function library for isentropic flow calculations
# Author: Maxime Sicat
# 28/12/2018

# ------------------------- IMPORT MODULES ----------------------

import math as m

# ------------------------ FUNCTION DEFINITIONS -----------------


def mach_via_pressure(gamma, p_total, p_static):
    """ Return the isentropic mach for a given set of total & static pressures
    as well as the flow's isentropic parameter gamma. Pressures can be expressed
    in any unit as long as they are both the same unit.
    """
    return m.sqrt( 2 / (gamma-1) * ( (p_total/p_static)**(1-1/gamma) - 1 ) )


def expansion_via_mach(gamma, mach):
    """ Return the expansion ratio := nozzle section / critical section
    given the flow's mach and isentropic parameter gamma.
    """

    e_square = ( ( 2 / (gamma+1) * (1 + (gamma-1) / 2 * (mach**2) ) ) ** ( (gamma+1) / (gamma-1) ) ) / (mach**2)
    return m.sqrt(e_square)


def exit_mach_via_expansion(gamma, expansion, supersonic):
    """Return the exit mach of a nozzle given its expansion ratio and the flow's isentropic parameter gamma.
    Function uses Newton's method of finding a function's zero to solve the expansion ratio equation for M.

    :param supersonic: boolean defining the exit flow as either subsonic (=false) or supersonic (=true)
    since both solutions exist for a given nozzle.
    """

    # Special case where flow is sonic

    if expansion == 1:
        return 1

    # Newton's method

    else:
        # Function used in Newton's method is defined as f, along with its derivative f_prime

        def f(M):
            k = (gamma - 1) / (gamma + 1)
            return (expansion * M)**(2*k)/k - M**2 - 2/(gamma-1)

        def f_prime(M):
            k = (gamma - 1) / (gamma + 1)
            return 2/M * ((expansion*M)**(2*k)) - 2*M

        # Scheme is initialized depending on supersonic/subsonic nature of the exit flow

        if supersonic:
            M = 3.0
        else:
            M = 0.3

        # Iterate until M reaches a set convergence criteria

        while abs(f(M)) > 0.00001:
            M = M - (f(M) / f_prime(M))

        return M


def pressure_ratio(gamma, mach):
    """ Return the total pressure by static pressure ratio for a given isentropic flow.
    """

    return  (1 + (gamma-1)/2 * mach**2)**(gamma / (gamma-1) )


def temperature_ratio(gamma, mach):
    """ Return the total temperature by static temperature ratio for a given isentropic flow.
    """

    return 1 + (gamma-1) / 2 * (mach**2)


def density_ratio(gamma, mach):
    """ Return the total density by static density ratio for a given isentropic flow.
    """

    return (1 + (gamma-1)/2 * mach**2)**( 1 / (gamma-1) )