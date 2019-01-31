# Hybrid rocket burn simulation function
# Author: Maxime Sicat
# 28/12/2018

# ------------------------- IMPORT MODULES ----------------------

import Geometries as G
import Nozzle as N
import Isentropic as iso
import math
import matplotlib.pyplot as plt
from Interpolator import *
from JsonInterpreter import JsonInterpreter

# ------------------------ FUNCTION DEFINITIONS -----------------

# Toutes les grandeurs sont exprimees en USI sauf précision autre

# Les simulations ne sont valables à priori qu'aux CNTP

def Simulation_debit_constant_fuel_sliver(ox_flow, geom, nozzle, safety_thickness, dt, show_curves):
    """
    Simulate the hybrid rocket burn. Curves can be shown afterwards.
    Every value is expressed in ISU unless specified otherwise. Simulation is only valid for normal ground
    atmosphere conditions (1 bar, 293 K).
    :param ox_flow: value of input oxidizer flow
    :param geom: grain geometry (implements the Geometry class)
    :param nozzle: Nozzle class object, defines the burn nozzle
    :param safety_thickness: minimum allowed thickness for fuel slivers, which conditions burn termination
    :param dt: time increment
    :param show_curves: boolean defining whether ot not simulation curves are plotted at the end
    :return: TBD
    """

    # Initialize the thermodynamic constants with the json interpreter

    Interpreter = JsonInterpreter("Thermodynamic Data Onera 40 bar H2O2 87_5.json")

    combustion_table = Interpreter.return_combustion_table()

    g0 = combustion_table['g0']
    R = combustion_table["R"]
    pression_atmo = combustion_table["Pa"]
    a = combustion_table["a"]
    n = combustion_table["n"]
    m = combustion_table["m"]
    rho_fuel = combustion_table["rho_fuel"]
    initial_chamber_pressure_bar = combustion_table["P_chamber_bar"]
    initial_chamber_pressure = initial_chamber_pressure_bar * (10 ** 5)

    interpolator = Interpolator("THEORETICAL ROCKET PERFORMANCE ASSUMING EQUILIBRIUM", "Thermodynamic Data Onera 40 bar H2O2 87_5.json")

    # Initialize useful variables. Not strictly necessary in Python, but here for clarity

    regression_rate = 0
    fuel_flow = 0
    total_mass_flow = 0
    OF = 0
    T_chambre = 3000
    gamma = 1.4
    masse_molaire_gaz = 20*10**(-3)
    c_star = 1500
    pression_chambre = initial_chamber_pressure

    # Initialize lists for data over time

    time = []
    thrust_curve = []
    isp_curve = []
    pressure_curve = []
    temperature_curve = []
    radius_curve = []
    OF_curve = []

    k = 0 # Iterator to keep track of loop progress


    # Main simulation loop

    while(geom.min_bloc_thickness() > safety_thickness):

        # Regression rate and mass flow calculations

        regression_rate = geom.Compute_regression_rate(ox_flow, a, n, m)
        print("")
        fuel_flow = geom.Compute_fuel_rate(regression_rate, rho_fuel)
        total_mass_flow = ox_flow + fuel_flow
        OF = ox_flow / fuel_flow

        # Data extraction using the interpolator and OF ratio
        # Tc, gamma, molar mass, cstar

        output = interpolator.interpolate_data(o_f_desired_value=OF, variables=["t", "m", "gammas", "cstar"], mole_fractions="")
        T_chambre = [output["variables"][i]["CHAMBER"] for i in range(0, len(output["variables"])) if output["variables"][i]["name"] == "t"][0]
        gamma  = [output["variables"][i]["CHAMBER"] for i in range(0, len(output["variables"])) if output["variables"][i]["name"] == "gammas"][0]
        masse_molaire_gaz =[output["variables"][i]["CHAMBER"] for i in range(0, len(output["variables"])) if output["variables"][i]["name"] == "m"][0] / 1000 # On divise par 1000 pour obtenir des USI
        CEA_c_star = [output["variables"][i]["THROAT"] for i in range(0, len(output["variables"])) if output["variables"][i]["name"] == "cstar"][0]
        Calculated_c_star = pression_chambre*nozzle.get_throat_area() / total_mass_flow
        print(str(k*dt) + " s")

        r = R / masse_molaire_gaz # Specific gaz constant

        # Calculate chamber conditions and motor performance
        # Chamber pressure is not recalculated as it is assumed constant

        Mach_sortie = iso.exit_mach_via_expansion(gamma, nozzle.get_expansion_ratio(), True)
        pression_sortie = pression_chambre / iso.pressure_ratio(gamma,Mach_sortie)
        T_sortie = T_chambre / iso.temperature_ratio(gamma, Mach_sortie)
        v_exit = math.sqrt(gamma * r * T_sortie) * Mach_sortie

        F = nozzle.get_nozzle_effeciency() * (total_mass_flow*v_exit + (pression_sortie-pression_atmo)*nozzle.get_exit_area())

        Isp = F / total_mass_flow / g0

        # Curve updates

        thrust_curve.append(F)
        isp_curve.append(Isp)
        pressure_curve.append(pression_chambre)
        temperature_curve.append(T_chambre)
        radius_curve.append(geom.get_port_radius())
        OF_curve.append(OF)
        time.append(k* dt)

        # Update the geometry and nozzle

        geom.regress(regression_rate, dt)
        nozzle.erode(dt)

        k += 1

    # Plot curves if required

    if(show_curves) :

        fig, axs = plt.subplots(4, 1, squeeze=True)
        axs[0].plot(time, thrust_curve)
        axs[0].set_title('')
        axs[0].set_xlabel('time (s)')
        axs[0].set_ylabel('Thrust (N)')
        fig.suptitle('Simulation results', fontsize=12)

        axs[1].plot(time, radius_curve, '--')
        axs[1].set_xlabel('time (s)')
        axs[1].set_title('')
        axs[1].set_ylabel('Port Radius (m)')

        axs[2].plot(time, isp_curve, 'r')
        axs[2].set_xlabel('time (s)')
        axs[2].set_title('')
        axs[2].set_ylabel('Isp (s)')

        axs[3].plot(time, OF_curve, 'g--')
        axs[3].set_xlabel('time (s)')
        axs[3].set_title('')
        axs[3].set_ylabel('O/F ratio')

        plt.show()

    avg_Isp = 0
    avg_OF = 0
    total_impulse = 0

    # Several averages I chose to show.
    # Not necessary in the full program as it can be processed later on.

    for k in range(len(time)):

        avg_Isp += isp_curve[k]
        avg_OF += OF_curve[k]
        total_impulse += dt * thrust_curve[k]

    avg_Isp = avg_Isp / len(isp_curve)
    avg_OF = avg_OF / len(OF_curve)

    print("Average Isp : " + str(avg_Isp) + " s")
    print("Average OF ratio : " + str(avg_OF))
    print("Total impulse : " + str(total_impulse) + " N.s")

