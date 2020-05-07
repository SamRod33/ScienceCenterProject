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

    def launch(self, tmax, dt, altitude):
        """
        Launches the rocket and outputs graph of simulated trajectory

        Parameters
        ----------
        tmax : int
            the maximum number of time steps the run the simulation for
        dt : int
            the time between each time step
        """

        x = [self.x]
        y = [self.y]
        v = [self.v]
        a = [self.a]
        fuel = [self.fuel]
        mass = [self.mass]

        launch_fuel = self.fuel  # calculate what percent of total fuel used for launch
        consume_fuel = 100

        for nt in np.arange(dt, tmax, dt):
            if self.getFuel() >= 0:
                self.setFuel(self.getFuel()-consume_fuel)
            if self.getFuel() - consume_fuel <= 0:
                self.setThrust(0)
                self.setFuel(0)

            force = self.getThrust() + self.wind(y[nt-1]) - self.gravity_force(y[nt-1]) + self.drag(y[nt-1])
            a.append(a[nt-1]+(force / (self.mass + self.getFuel())))
            v.append(a[nt-1]*dt + v[nt-1])
            y.append(0.5*a[nt-1]*dt**2 + v[nt-1]*dt + y[nt-1])
            if y[nt] < 0:
                break

        self.visualize(y, v, a, nt, dt)


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
        R = 287.05  # Specific gas constant for dry RADIUS_EARTH
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
            the altitude of rocket at each time step
        """
        return s

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
        # f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize=(8,8))
        # f.tight_layout(pad=3.0)
        # x = np.arange(0,nt+dt,dt)
        #
        # x_data = []
        # y_data = []
        # v_data = []
        # a_data = []
        #
        #
        # pos, = ax1.plot(x, s, 'b')
        # ax1.set_ylabel('Altitude')
        # ax1.set_title('Position vs. Time of Rocket Launch')
        #
        # vel, = ax2.plot(x, v, 'g')
        # ax2.set_ylabel('Velocity')
        # ax2.set_title('Velocity vs. Time of Rocket Launch')
        #
        # accel, = ax3.plot(x, a, 'r')
        # ax3.set_ylabel('Acceleration')
        # ax3.set_xlabel('Time (sec.)')
        # ax3.set_title('Acceleration vs. Time of Rocket Launch')
        #
        # def animate(i):
        #     x_data.append(x[i])
        #     y_data.append(s[i])
        #     v_data.append(v[i])
        #     a_data.append(a[i])
        #
        #     pos.set_xdata(x_data)
        #     pos.set_ydata(y_data)
        #
        #     vel.set_xdata(x_data)
        #     vel.set_ydata(v_data)
        #
        #     accel.set_xdata(x_data)
        #     accel.set_ydata(a_data)
        #
        #
        # animation = FuncAnimation(f, func=animate, frames=len(s), interval=10)
        # Set up subplots
        f, axs = plt.subplots(nrows=3, ncols=2, sharex=True, figsize=(14,9))
        f.tight_layout(pad=3.0)
        f.set_facecolor('w')

        # Set up common axis and title names
        f.text(0.5, 0.02, 'Time After Launch [sec.]', ha='center', fontsize=14)
        f.text(0.5, 0.96, 'Modeling a Rocket Launch', ha='center', fontsize=16)

        # Set up individual axis titles
        axs[0,0].set_ylabel('Altitude [km]')
        axs[0,0].set_title('Position vs. Time')

        axs[1,0].set_ylabel('Velocity [km/s]')
        axs[1,0].set_title('Velocity vs. Time')

        axs[2,0].set_ylabel('G Force')
        axs[2,0].set_title('Acceleration vs. Time')


        # Remove trailing 0s if rocket reached terminal velocity by tmax
        if s[-1] <= 0:
            s, v, a = map(lambda x: np.trim_zeros(x, 'b'), [s, v, a])

        # Calculate common x values and g forces throughout launch
        x = np.arange(0, int(len(s)*dt), dt)
        gs = self.calculate_g(a)

        # Combine data from launch into only a few variables
        axes = [axs[0,0], axs[1,0], axs[2,0]]
        data = [s, v, a]
        scale_factors = [1000, 1000, 1]
        colors = ['b', 'g', 'r']

        # Set the x and y axis values for each subplot
        for ax, d, sf in zip(axes, data, scale_factors):
            ax.set_ylim(0, 1.3 * np.amax(d) / sf)
            ax.set_xlim(0, np.amax(x))

        # Pseudo-animate the launch
        for i in range(0, len(x)):
            for ax, d, sf, color in zip(axes, data, [1000, 1000, 1  ], colors):
                ax.plot(x[:i], d[:i] / sf, color)
            plt.pause(0.001)
        plt.show()
