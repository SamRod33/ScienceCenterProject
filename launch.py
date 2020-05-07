from constants import *
import matplotlib.pyplot as plt
import numpy as np
import time


class Rocket:

    crashed = False
    in_orbit = True

    def __init__(self):
        self.mass_rocket = MASS_ROCKET
        self.mass_fuel = MASS_FUEL
        self.pitch_angle = np.pi / 2


    def student_input(self):

        # Instructions
        print("Please provide all values with at least 3 significant figures.")
        # Ask about value of G
        while True:
            G = raw_input("\nWhat is the value of the gravitational constant?\n   _____ x 10^-11 m^3 k^-1 s^-2 >>  ")
            try:
                if abs(float(G) - GRAVITATIONAL_CONSTANT / 10**-11) > ALLOWED_ERROR:
                    print("Wrong value provided, physical properties malfunctioning. Abort.")
                else:
                    break
            except:
                print("Please provide an actual number")

        # Earth
        while True:
            mass_earth = raw_input("\nWhat is the mass of the Earth?\n    _____ x 10^23 kg >>  ")
            try:
                if abs(float(mass_earth) - MASS_EARTH/ 10**23) > ALLOWED_ERROR:
                    print("Wrong value provided. Earth's orbit unstable. Urgent fix required.")
                else:
                    break
            except:
                print("Please provide an actual number")

        while True:
            radius_earth = raw_input("\nWhat is the radius of the Earth?\n    _____ x 10^6 m >>  ")
            try:
                if abs(float(radius_earth) - RADIUS_EARTH / 10**6) > ALLOWED_ERROR:
                    print("Wrong value provided, Earth's orbit unstable. Urgent fix required.")
                else:
                    break
            except:
                print("Please provide an actual number")

        # Mars
        while True:
            mass_mars = raw_input("\nWhat is the mass of Mars\n    _____ x 10^23 kg >>  ")
            try:
                if abs(float(mass_mars) - MASS_MARS / 10**23) > ALLOWED_ERROR:
                    print("Wrong value provided, Mar's orbit unstable. Urgent fix required.")
                else:
                    break
            except:
                print("Please provide an actual number")

        while True:
            radius_mars = raw_input("\nWhat is the radius of Mars\n    _____ x 10^6 m >>  ")
            try:
                if abs(float(radius_mars) - RADIUS_MARS / 10**6) > ALLOWED_ERROR:
                    print("Wrong value provided, Mar's orbit unstable. Urgent fix required.")
                else:
                    break
            except:
                print("Please provide an actual number")

        # Begin Launch
        while True:
            input = raw_input("\nAre you ready to launch?\n   Input Y for yes and N for no. >>  ")
            if input == 'Y':
                break
            else:
                print("Launch sequence aborted.")
                return

        print("Initiating launch sequence:")
        for i in range(10,0,-1):
            print(i)
            time.sleep(1)

        print("Blast Off!")


    def launch(self, tmax, dt):
        """
        Launches the rocket and outputs graph of simulated trajectory

        Parameters
        ----------
        tmax : int
            the maximum number of time steps the run the simulation for
        dt : int
            the time between each time step
        """
        # self.student_input()

        nt = int(tmax / dt)
        s, v, a = np.zeros(nt), np.zeros(nt), np.zeros(nt)  # Position, velocity, acceleration
        fg, fd, ft = np.zeros(nt), np.zeros(nt), np.zeros(nt)  # Force gravity, drag, thrust

        launch_fuel = self.mass_fuel  # calculate what percent of total fuel used for launch

        for i in range(1, nt):

            # calculate all relevant forces
            fg[i] = self.force_gravity(altitude=s[i-1])
            fd[i] = self.force_drag(altitude=s[i-1], velocity=v[i-1])
            ft[i] = self.force_thrust(time=i*dt, force_gravity=fg[i], force_drag=fd[i], altitude=s[i-1])
            force = ft[i] - fg[i] - fd[i]

            # consumer fuel
            self.consume_fuel(velocity=v[i-1], net_thrust=ft[i-1]-fd[i-1], thrust=ft[i-1], dt=dt)
            if self.mass_fuel < 0:
                crashed = True

            # Begin pitch maneuver between in time interval [160:400] seconds
            self.tilt_maneuver(i*dt, dt)

            # Calculate updated acceleration, velocity
            a[i] = force / (self.mass_rocket + self.mass_fuel)
            v[i] = a[i-1] * dt + v[i-1]

            # Calculate altitude based on pitch angle and current escape velocity
            if v[i-1] > ESCAPE_VELOCITY:
                self.in_orbit = False
            if self.in_orbit:
                s[i] = np.sin(self.pitch_angle) * (0.5 * a[i-1] * dt**2 + v[i-1] * dt) + s[i-1]
            else:
                s[i] = s[i-1]
                break

        return s, v, a, i, dt, fg, fd, ft


    def consume_fuel(self, velocity, net_thrust, thrust, dt):
        """
        Consumes fuel for a single time step

        Parameters
        ----------
        velocity : double
            the current velocity of the rocket in m/s
        force : double
            the current net force on the rocket in N
        dt : double
            the time step for this simulation in seconds
        """
        # https://www.grc.nasa.gov/WWW/K-12/airplane/sfc.html
        # Must take into account flow choking
        # Does not take into account changing MACH speed as temperature decreases
        if net_thrust == 0:
            return
        mass_flow_rate = METHANE_DENSITY * (velocity if velocity < MACH_1 else MACH_1) * ENGINE_AREA
        self.mass_fuel -= mass_flow_rate * thrust / net_thrust  * dt


    def force_gravity(self, altitude):
        """
        Calculates the force of gravity currently acting on the rocket

        Parameters
        ----------
        altitude:
            the current height of the rocket
        """
        return GRAVITATIONAL_CONSTANT * (self.mass_rocket + self.mass_fuel) * MASS_EARTH / (RADIUS_EARTH + altitude)**2


    def force_drag(self, altitude, velocity):
        """
        Returns the force of drag curently acting on the rocket

        Parameters
        ----------
        altitude:
            the current altitude of the rocket
        velocity:
            the current velocity of the rocket in m/s
        """
        return 0.5 * DRAG_COEFF * ORTH_SURFACE_AREA * self.rho(altitude) * velocity**2


    def force_thrust(self, time, force_gravity, force_drag, altitude):
        """
        Returns the current thrust of the rocket

        Source
        ------
        WIKI:
            https://upload.wikimedia.org/wikipedia/commons/2/2c/Apollo17_Ascent_Trajectory.pdf
        """
        if time < 160:
            G_force = 0.017 * time + 1
        elif time < 560:
            G_force = .0018 * time + 0.6
        elif time < 760:
            G_force = 0.001 * time + 0.5
        else:
            G_force = 0.0
        G = GRAVITATIONAL_CONSTANT * MASS_EARTH / (RADIUS_EARTH + altitude)**2
        thrust = (MASS_ROCKET + self.mass_fuel) * G_force * G + force_drag + force_gravity
        return thrust if thrust < BOOSTER_THRUST else BOOSTER_THRUST


    def tilt_maneuver(self, time, dt):
        """
        Changes pitch angle of rocket to orbit earth until reached escape velocity

        Parameters
        ----------
        time :
            the time in seconds after liftoff
        """
        if time > 160 and time < 400:
            self.pitch_angle -= np.pi / (480 / dt)


    def rho(self, altitude):
        """
        Calculates the current air density in kg/m^3

        Parameters
        ----------
        altitude:
            the altitude to find desired air density

        Source
        ------
        NASA:
            https://www.grc.nasa.gov/WWW/K-12/airplane/atmosmet.html
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
            p = 22.65 * np.exp(1.73 - 0.000157 * altitude)
        return p / (0.2869 * (T + 273.1))


    def calculate_g(self, s):
        """
        Generate a list of values of G

        Parameters
        ----------
        s : list
            the altitude of rocket at each time step
        """
        return [GRAVITATIONAL_CONSTANT * MASS_EARTH / (RADIUS_EARTH + altitude)**2 for altitude in s]


# Launch Rocket
rocket = Rocket()
s, v, a, i, dt, fg, fd, ft = rocket.launch(800, 1)

# -----------------
#    VISUALIZE
# -----------------

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

axs[0,1].set_ylabel('Force Gravity [kN]')
axs[0,1].set_title('Force Gravity vs. Time')

axs[1,1].set_ylabel('Force Drag [kN]')
axs[1,1].set_title('Force Drag vs. Time')

axs[2,1].set_ylabel('Force Thrust [kN]')
axs[2,1].set_title('Force Thrust vs. Time')

# Remove trailing 0s if rocket reached terminal velocity by tmax
if not rocket.in_orbit:
    s, v, a, fg, fd, ft = map(lambda x: np.trim_zeros(x, 'b'), [s, v, a, fg, fd, ft])

# Calculate common x values and g forces throughout launch
x = np.arange(0, int(len(s)*dt), dt)
gs = rocket.calculate_g(s)

# Combine data from launch into only a few variables
axes = [axs[0,0], axs[1,0], axs[2,0], axs[0,1], axs[1,1], axs[2,1]]
data = [s, v, a, fg, fd, ft]
scale_factors = [1000, 1000, np.amax(gs), 1000, 1000, 1000]
colors = ['b', 'g', 'r', 'c', 'm', 'k']

# Set the x and y axis values for each subplot
for ax, d, sf in zip(axes, data, scale_factors):
    ax.set_ylim(0, 1.3 * np.amax(d) / sf)
    ax.set_xlim(0, np.amax(x))

# Pseudo-animate the launch
for i in range(0, len(x)):
    for ax, d, sf, color in zip(axes, data, [1000, 1000, gs[:i], 1000, 1000, 1000], colors):
        ax.plot(x[:i], d[:i] / sf, color)
    plt.pause(0.001)

plt.show()
