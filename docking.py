#some code
import matplotlib.pyplot as plt
import numpy as np
import math
from constants import *
from Rocket import *


# Equations

def accel_grav(d):
    """
    Returns: current acceleration due to gravity at a distance r from mars
    """
    g = (GRAVITATIONAL_CONSTANT * m_mars) / (d+r)**2
    return g


# Parameters
r_init = 4*(10**6) # m distance away from planet during landing/docking stage
v0 = -12.1*10**3 # initial velocity (m/s) ### Will have to change
m_mars = 6.39*(10**23)
r = 3.3895*(10**6)
G = 6.67408*(10**-11)
F_thrust1 = 1.46*10**5 # kg*m/s^2 Reverse thrust of rocket
F_thrust2 = 2.73*10**5
F_thrust3 = 1.5*10**5
v_ex = -11*10**(2) # m/s velocity of ejecting fuel
mass_fuel = 88517


#creating Rocket with current position, velocity, and current fuel left
rocket = Rocket(y=r_init, v=v0, a=-accel_grav(r_init+r),
thrust=F_thrust1, fuel = mass_fuel)


# Time for Simulation
dt = 1 # day
tmax = 500 # days
timestep = int(tmax/dt)

#initializers
m,y,v,a,time= np.zeros(timestep), np.zeros(timestep), \
np.zeros(timestep), np.zeros(timestep),np.zeros(timestep)

m[0] = rocket.getMass() + rocket.getFuel()
y[0] = rocket.getY()
v[0] = rocket.getVelocity()
a[0] = rocket.getAcceleration()
time[0] =  0

#Running simulation
ifinal = -1
for i in range(1, timestep):
    dm = rocket.getThrust()/v_ex
    if m[i-1] >= rocket.getMass():
        m[i]=(m[i-1] + dm * dt)
        rocket.setFuel(rocket.getFuel() + dm)
    else:
        m[i]=(m[i-1])
        rocket.setThrust(0)
    wind = rocket.wind(y[i-1])
    f_g = accel_grav(y[i-1])
    thrust = rocket.getThrust()
    a[i]=(a[i-1] * dt - f_g + (thrust+wind) / m[i-1])
    v[i]=(v[i-1] + a[i-1] * dt)
    y[i]=(y[i-1] + v[i-1] * dt)
    time[i]=(time[i-1] + dt)

    if a[i] / 9.8 < -0.75:
        rocket.setThrust(F_thrust2)
    if v[i] / 1000 > -2:
        rocket.setThrust(F_thrust3)

    if y[i-1] <=0:
        break

rocket.visualize(y, v, a, i, dt)
