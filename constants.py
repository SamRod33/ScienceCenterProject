"""
To use:
    add `from constants import *` to header
"""

# Student Constants
ALLOWED_ERROR = .2

# Physical Constants
GRAVITATIONAL_CONSTANT = 6.67 * 10**-11
MASS_EARTH = 5.9724 * 10**24  # kg
MASS_MARS = 0.64171 * 10**24  # kg
RADIUS_EARTH = 6.378 * 10**6  # m
RADIUS_MARS = 3.3895 * 10**6 # m
AU = 146230000000  # m
ESCAPE_VELOCITY = 7800  # m/s

# Rocket Constants
# Source: https://www.spacex.com/sites/spacex/files/making_life_multiplanetary-2017.pdf
MASS_FUEL = 9.98 * 10**5  # kg
MASS_ROCKET = 77.111 * 10**3  # kg
ORTH_SURFACE_AREA = 3.1415 * 4.5**2  # m^2
DRAG_COEFF = 0.25  # Rough estimate from research

BOOSTER_THRUST = 5.27 * 10**7  # Newtowns
ENGINE_AREA = 20.750  # m^2
METHANE_DENSITY = 0.226  # kg /m^3 @ 580 Celcius (ignition temp of CH4)
MACH_1 = 332  # m/s 