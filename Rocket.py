# Documentation of Rocket class

from constants import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

class Rocket:

    def __init__(self, x=0, y=0, v=0, a=0, thrust = 1, fuel=MASS_FUEL, mass=MASS_ROCKET):
        """
        Returns: None

        Constructs the Rocket object based on:
            -X coordinate
            -Y coordinate
            -Velocity of Rocket
            -Acceleration of Rocket
            -Thrust exterted by Rocket
            -Fuel remaining
            -Mass of Rocket
        """
        self.x = x
        self.y = y
        self.v = v
        self.a = a
        self.thrust =thrust
        self.fuel = fuel
        self.mass = mass

    def setX(self, x):
        """
        Returns: None

        Sets new X coordinate to rocket
        """
        assert type(x) == float or type(x) == int
        self.x = x

    def setY(self, y):
        """
        Returns: None

        Sets new Y coordinate to rocket
        """
        assert type(y) == float or type(y) == int
        self.y = y

    def setVelocity(self, v):
        """
        Returns: None

        Sets new velocity to Rocket
        """
        assert type(v) == float or type(v) == int
        self.v = v

    def setAcceleration(self, a):
        """
        Returns: None

        Sets new acceleration to Rocket
        """
        assert type(a) == float or type(a) == int
        self.a = a

    def setThrust(self, thrust):
        """
        Returns: None

        Sets new thrust to Rocket
        """
        assert type(thrust) == float or type(thrust) == int
        self.thrust = thrust

    def setFuel(self, fuel):
        """
        Returns: None

        Sets new Fuel mass to Rocket
        """
        assert type(fuel) == float or type(fuel) == int
        self.fuel = fuel

    def getX(self):
        """
        Returns: X-coordinate
        """
        return self.x

    def getY(self):
        """
        Returns: Y-coordinate
        """
        return self.y

    def getVelocity(self):
        """
        Returns: Velocity of Rocket
        """
        return self.v

    def getAcceleration(self):
        """
        Returns: Acceleration of Rocket
        """
        return self.a

    def getThrust(self):
        """
        Returns: Thrust of rocket
        """
        return self.thrust

    def getFuel(self):
        """
        Returns: Fuel of Rocket
        """
        return self.fuel

    def getMass(self):
        """
        Returns: Mass of Rocket (w/out fuel)
        """
        return self.mass

    def gravity_force(self, altitude):
        """
        Calculates the force of gravity currently acting on the rocket
        """
        return GRAVITATIONAL_CONSTANT * (self.mass + self.getFuel()) * MASS_EARTH / (altitude)**2

    def drag(self, altitude):
        """
        Returns the force of drag curently acting on the rocket

        Parameters
        ----------
        rho : float
            density of fluid rocket travelling through
        """
        # TODO: model rho with equation based on current altitude
        return 0.5 * DRAG_COEFF * ORTH_SURFACE_AREA * self.rho(altitude) * self.getVelocity()**2

    def wind(self, altitude):
        """
        Returns the force of wind acting on the rocket
        """
        wind_velocity = 1
        return 0.5 * self.rho(altitude) * wind_velocity**2 * ORTH_SURFACE_AREA

    def rho(self, altitude):
        """
        Calculates the current air density in kg/m^3

        Source: https://www.grc.nasa.gov/WWW/K-12/airplane/atmosmet.html
        """
        if altitude > 25000:
            T = -131.21 + 0.003 * altitude
            p = 2.488 * ((T + 273.1) / 216.6)**-11.388
        elif altitude < 11000:
            T = 15.04 - 0.00649 * altitude
            p = 101.29 * ((T + 273.1) / 288.08)**5.256
        else:
            T = -56.64
            p = 22.65 * np.exp(1.73 - 0.000157*altitude)
        return p / (0.2869 * (T + 273.1))

    def calculate_g(self, s):
        """
        Generate a list of values of G

        Parameters
        ----------
        s : list
            the acceleration of rocket at each time step
        """
        g_force = np.zeros(len(s))
        for i in range(len(s)):
            g_force[i] = s[i] / 9.8

        return g_force

    def visualize(self, s, v, a, nt, dt):
        """
        Parameters
        ----------
        s : list
            the position of the rocket
        v : list
            the velocity of the rocket
        a : list
            the acceleration of the rocket
        nt : int
            the number of time steps
        dt : int
            the time between each time step
        """
        # Set up subplots
        f, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(14,9))
        f.tight_layout(pad=3.0)
        f.set_facecolor('w')

        # Set up common axis and title names
        f.text(0.5, 0.02, 'Time During Docking [sec.]', ha='center', fontsize=14)
        f.text(0.5, 0.96, 'Rocket Docking Model', ha='center', fontsize=16)

        # Set up individual axis titles
        axs[0].set_ylabel('Altitude [km]')
        axs[0].set_title('Position vs. Time')

        axs[1].set_ylabel('Velocity [km/s]')
        axs[1].set_title('Velocity vs. Time')

        axs[2].set_ylabel('G Force')
        axs[2].set_title('Acceleration vs. Time')


        # Remove trailing 0s if rocket reached terminal velocity by tmax
        if s[-1] <= 0:
            s, v, a = map(lambda x: np.trim_zeros(x, 'b'), [s, v, a])

        # Calculate common x values and g forces throughout launch
        x = np.arange(0, int(len(s)*dt), dt)
        gs = self.calculate_g(a)

        # Combine data from launch into only a few variables
        axes = [axs[0], axs[1], axs[2]]
        data = [s, v, gs]
        scale_factors = [1000, 1000, 1]
        colors = ['b', 'g', 'r']

        # Set the x and y axis values for each subplot
        for ax, d, sf in zip(axes, data, scale_factors):
            if axes == axs[2]:
                break
            ax.set_ylim(0, 1.3 * np.amax(d) / sf)
            ax.set_xlim(0, np.amax(x))
        axs[1].set_ylim(-15 , 0)
        axs[2].set_ylim(-5, 10)
        axs[2].set_xlim(0, np.amax(x))

        # Pseudo-animate the launch
        for i in range(0, len(x)):
            for ax, d, sf, color in zip(axes, data, [1000, 1000, 1], colors):
                ax.plot(x[:i], d[:i] / sf, color)
                
        plt.pause(0.001)
        plt.show()
