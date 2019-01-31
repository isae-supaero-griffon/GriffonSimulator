# -*- coding: utf-8 -*-
"""
Rocket parabolic flight calculator, from Thrust Curve - based on Eric Brunner code.
Author: Iñigo Reizabal Arregui
Date: 14 Dec 2018
"""

# ------------------- IMPORT MODULES ----------------------

import math
import numpy as np
import matplotlib.pyplot as plt

# ----------------- FUNCTION DEFINITIONS ------------------


def rho(h):
    """
    Density of the atmosphere as a function of the altitude
    for the I.S.A. (International standard atmosphere)
    """  
    g0=9.81 # m/s^2
    R=8.314 # N.m/(mol.K) Ideal gas constant
    M=0.02896 # kg/mol, Molar mass of air
    rho_0=1.225 # kg/m3
    T_0=298 # K
    h_0=0 # m
    rho = rho_0 * math.exp(-g0*M*(h-h_0)/(R*T_0)) 
    return rho


def masscalc(Vectorthrust,Isp,g0,m0):
    """
    Calculation of the mass in the point where the thrust curve is defined
    for a constant Isp
    Neccesary variables are
        -Thrust vector, defined as 2 line n column matrix where 
         first line is the thrust in a discrete time and second line is 
         the time
        -m0 is the mass of the rocket at t=0
        -g0=gravity (9,81m/s^2)
        -Isp is the specific impulse in s
    It returns a vector with 3 lines and n colums where 
        -First line is the mass at each time
        -Second line is the time
        -Third line is the dm(kg/s) at each time
    
    """

    # Initialize the values
    dm = 0
    pos = 0
    Vectormass = np.zeros((3, T.shape[1]))
    Vectormass[0, pos] = m0
    Vectormass[1, pos] = Vectorthrust[1, pos]

    # Run the loop
    for time in Vectorthrust[1, 1:]:
        pos = pos+1
        Vectormass[0, pos] = Vectormass[0, pos-1] - dm
        Vectormass[1, pos] = Vectorthrust[1, pos]
        Vectormass[2, pos-1] = dm
        b = Vectorthrust[0, pos]
        a = Vectorthrust[0, pos-1]
        t1 = Vectorthrust[1, pos]
        t0 = Vectorthrust[1, pos-1]
        dm = (b*(t1-t0)/(Isp*g0))
        
    Vectormass[0, pos] = Vectormass[0, pos-1]
    Vectormass[1, pos] = Vectorthrust[1, pos]
    Vectormass[2, pos] = dm

    # Return output
    return Vectormass


def interpol(Vector, line1, line2, pos, x):
    """
    Linear interpolation for a vector where:
        -line1 are the points where "y" is defined
        -line2 are the points where "x" is defined
        -pos is the position the n value in the vector
        -x is the variable that has to be interpolated to obtain y
        y(x)=y(x1)+(y(x2)-y(x1))/(x2-x1)*(x-x1)
    """
    m = (Vector[line1, pos+1] - Vector[line1, pos]) / (Vector[line2, pos+1] - Vector[line2, pos])
    y = Vector[line1, pos] + m*(x - Vector[line2, pos])

    # Return output
    return y
    

def Newt(Vectorthrust, Vectormass, S, Cd, g0, pos, t, rho, h, v, dv):
    """
    Calculation of the acceleration at discrete time=t
        -Thrust vector, defined as 2 line n column matrix where 
         first line is the thrust(N) in a discrete time and second line is 
         the time(s)
        -Mass vector, defined as 3 line n column matrix where 
         first line is the mass(kg) at each time, second line is the time(s)
         and third line is the dm(kg/s) at each time
        -S surface of the rocket for drag calculation (m^2)
        -Cd is the drag coeficient (constant)
        -g0 is gravity (9,81m/s^2)
        -pos is the position of the time in the vector
        -t(s) is the time for which we calculate the acc.
        -rho is the density as a function of h
        -h is the altitude
        -v is the velocity
        -dv is the diferential of the velocity of an intermedium step
    Hypotesis:
        -for each step, rho(n+1)=rho(n)
        -Cd is constant
        
    """

    T = interpol(Vectorthrust, 0, 1, pos, t)

    # Set T to 0 if it yields a value less than 0
    if T < 0:
        T = 0

    # Perform newton's method
    D = -0.5*rho(h)*S*Cd*(v + dv)*(v + dv)
    G = -interpol(Vectormass, 0, 1, pos, t)*g0
    V = -(v + dv)*interpol(Vectormass, 2, 1, pos, t)
    Acc = (T + D + G + V) / interpol(Vectormass, 0, 1, pos, t)
    return Acc
    

def RK4(Vectorthrust, Vectormass, S, Cd, g0, rho):
    """
    -Thrust vector, defined as 2 line n column matrix where 
         first line is the thrust(N) in a discrete time and second line is 
         the time(s)
        -Mass vector, defined as 3 line n column matrix where 
         first line is the mass(kg) at each time, second line is the time(s)
         and third line is the dm(kg/s) at each time
        -S surface of the rocket for drag calculation (m^2)
        -Cd is the drag coeficient (constant)
        -g0 is gravity (9,81m/s^2)
        -rho is the density as a function of h
    It returns a vector with 4 lines and n colums where 
        -First line is the velocity(m/s) at each time (line 0)
        -Second line is the time (line 1)
        -Third line is the acceleration(m/s^2 at each time (line 2)
        -Forth line is the altitude(m) at each time (line 3)
    """

    # Initialize the loop
    V = np.zeros((4, T.shape[1]))
    pos = 0
    v = 0.0
    h = 0.0
    V[0, pos] = v
    V[1, :] = Vectorthrust[1, :]

    # Run the loop within a try-catch block
    for time in Vectorthrust[1, :]:
        try:
            v1 = v
            # Perform integration in time
            deltat = Vectorthrust[1, pos + 1] - time
            dv1 = deltat*Newt(Vectorthrust, Vectormass, S, Cd, g0, pos, time, rho, h, v, 0)
            dv2 = deltat*Newt(Vectorthrust, Vectormass, S, Cd, g0, pos, time + deltat/2, rho, h, v, dv1/2)
            dv3 = deltat*Newt(Vectorthrust, Vectormass, S, Cd, g0, pos, time + deltat/2, rho, h, v, dv2/2)
            dv4 = deltat*Newt(Vectorthrust, Vectormass, S, Cd, g0, pos, time + deltat/2, rho, h, v, dv3)
            v = v + (dv1 + 2*dv2 + 2*dv3 + dv4)/6
            h = h+v*deltat
            pos = pos + 1
            V[0, pos] = v                           # store velocity
            V[2, pos] = (v - v1)/deltat                  # store acceleration
            V[3, pos] = h                              # store altitude

            # Break the cycle when the rocket returns to ground
            if h < 0:
                break

        except (ArithmeticError, ValueError, ZeroDivisionError):
            print("Error encountered while running RK4 loop. \n")
    return V

if __name__ == '__main__':

    # ----------------- Declaration of variables ------------------
    g0 = 9.81 #m/s^2
    Isp = 245 #s
    m0 = 28 #kg
    S = math.pi*(0.16*0.16)/4
    Cd = 0.75

    # ----------------- Thrust vector ------------------
    # This part creates a thrust vector to be able to check the programme
    # For the final programme, this part should be removed and the trust vector
    # has to be defined as a matrix of 2 lines and n columns T(2,n)
    deltat = 0.01             # delta-time in seconds
    simulation_time = 45      # simulation-time, adjust simulation time in order to show the entire chart
    size_T = int(simulation_time/deltat)
    T = np.zeros((2, size_T))
    n = 0
    t = 0
    burn_time = 3.5         # Burn time in seconds
    constant_thrust = 3000  # Thrust value in newtons
    for i in T[1, :]:
        T[1, n] = t
        if t < burn_time:
            T[0, n] = constant_thrust
        t = t + deltat
        n = n + 1

    # ----------------- Main ------------------
    m=masscalc(T, Isp, g0, m0)
    X=RK4(T, m, S, Cd, g0, rho)                      # Velocity, Time, Acceleration & Altitude vector
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(20, 15), facecolor='w')
    ax1.plot(X[1, :], X[3, :], label='altitude', color='green')
    ax2.plot(X[1, :], X[0, :], label='velocity', color='red')
    ax3.plot(X[1, :], X[2, :], label='acceleration', color='blue')
    ax1.set(ylabel='Altitude [m]', title='Altitude Simulation')
    ax2.set(ylabel='Velocity [m/sec]')
    ax3.set(ylabel='Acceleration [m/sec^2]', xlabel='Time [secs]')
    ax1.grid(b=True, axis='both')
    ax2.grid(b=True, axis='both')
    ax3.grid(b=True, axis='both')
    plt.show()
