"""
Science Center Project: Inflight Portion

Authors: Callen Reid, Ashley Fan, Cher Xu, Efrain Munoz
Date: 04/15/20
"""

import turtle
import launch
import math
from constants import *
from turtle import *

SCALE = 225.0 / AU

class planet(Turtle):
    """
    Class to represent a planet

    Attribute xloc: the x location of the planet
    Invariant: xloc is a float

    Attribute yloc: the y location of the planet
    Invariant: yloc is a float

    Attribute name: the name of the planet
    Invariant: name is a non-empty string
    """
    xloc = yloc = 0.0
    name=''
    mass = 0.0

    def getX(self):
        return self.xloc

    def getY(self):
        return self.yloc

class spaceship(planet):
    """
    Class to represent a spaceship
    """

    # Attribute xloc: the x location of the spaceship
    # Invariant: xloc is a float >= 0.0 & <= 800.0
    #
    # Attribute yloc: the y location of the spaceship
    # Invariant: yloc is a float >= 0.0 & <= 800.0
    #
    # Attribute xVel: the magnitude of the x-component of velocity of the spaceship
    # Invariant: xVel is a float >= 0.0
    #
    # Attribute yVel: the magnitude of the y-component of velocity of the spaceship
    # Invariant: yVel is a float >= 0.0
    #
    # Attribute fuel: the remaining fuel in the spaceship
    # Invariant: fuel is a float >= 0.0

    def getX(self):
        """
        Returns the x location of the spaceship object
        """
        return self.xloc

    def setX(self, x):
        """
        Sets the x-location of the spaceship object

        Parameter x: The x-location of the spaceship object
        Precondition: x is a float >= 0.0 & <= 800.0
        """
        assert isinstance(x, float), 'Invalid type for x, x must be a float'

        self.xloc = x

    def getY(self):
        """
        Returns the y location of the spaceship object
        """
        return self.yloc

    def setY(self, y):
        """
        Sets the y location of the spaceship object

        Parameter y: The y-location of the spaceship object
        Precondition: y is a float >= 0.0 & <= 800.0
        """
        assert isinstance(y, float), 'Invalid type for y, y must be a float'
        # assert (y >= 0.0 and y <= 800.0), 'Invalid value for y, y must be between 0 and 800 inclusive'
        self.yloc = y

    def getXVel(self):
        """
        Returns the x-velocity of the spaceship object
        """
        return self.xVel

    def setXVel(self, vx):
        """
        Sets the x-velocity of the spaceship object

        Parameter vx: The x-velocity of the spaceship object
        Precondtion: vx is a float >=0.0
        """
        assert isinstance(vx, float), 'Invalid type for vx, vx must be a float'
        # assert (vx >= 0.0), 'Invalid value for vx, vx must be greater than 0.0'
        self.xVel = vx

    def getYVel(self):
        """
        Returns the y-velocity of the spaceship object
        """
        return self.yVel

    def setYVel(self, vy):
        """
        Sets the y-velocity of the spaceship object

        Parameter vy: The y-velocity of the spaceship object
        Precondtion: vy is a float >=0.0
        """
        assert isinstance(vy, float), 'Invalid type for vy, vy must be a float'
        # assert (vy >= 0.0), 'Invalid value for vy, vy must be greater than 0.0'
        self.yVel = vy

    def getFuel(self):
        """
        Returns the amount of remaining fuel in the spaceship object
        """
        return self.fuel

    def setFuel(self, f):
        """
        Sets the amount of fuel in the spaeship object

        Parameter f: The amount of fuel in the spaceship object
        Precondition: f is a float >=0.0
        """
        assert isinstance(f, float), 'Invalid type for f, f must be a float'
        # assert (vy >= 0.0), 'Invalid value for f, f must be greater than 0.0'
        self.fuel = f

    def __init__(self, xCoord, yCoord, alt, ang, vel, fuel = float(MASS_FUEL)):
        """
        Intiales a space spaceship


        Parameter alt: The initial altitude of the rocket w.r.t Earth
        Precondition: alt is a float > 0.0

        Parameter vel: The magnitude of the inital velocity of the rocket
        Precondition: vel is a float > 0.0

        Parameter ang: The angle of launch with respect to the positive x-axis
        Precondition: ang is a float in the range [0.0..360.0]

        Parameter fuel: The amount of remaining fuel in the rocket
        Preconditoin: Fuel is a float >= 0.0
        """
        assert isinstance(alt,float), 'Invalid type for alt, alt must be a float'
        assert alt > 0.0, 'Invalid value for alt, alt must be greater than 0.0'

        assert isinstance(ang,float), 'Invalid type for alt, alt must be a float'
        assert (ang >= 0.0 and ang <= 360.0), 'Invalid value for alt, alt must be in the range [0.0..360.0]'

        assert isinstance(fuel,float), 'Invalid type for fuel, fuel must be a float'
        assert fuel > 0.0, 'Invalid value for fuel, fuel must be greater than 0.0'

        self.setX(xCoord + math.cos(ang)*alt)
        self.setY(yCoord + math.sin(ang)*alt)
        self.setFuel(fuel)
        self.xVel = vel*math.cos(ang)
        self.yVel = vel*math.cos(ang)
        self.mass = MASS_ROCKET+MASS_FUEL
        self.alt = 0.0
        self.ang = 0.0
        planet.__init__(self)

    def thrust(self):
        self.vx=AU * -0.02 / 86400

    def attraction(self, other):
        rx = other.xloc - self.xloc
        ry = other.yloc-self.yloc
        r = math.sqrt((rx**2+ry**2))

        f = GRAVITATIONAL_CONSTANT*(self.mass*other.mass)/r**2

        theta = math.atan2(ry,rx)
        fx = math.cos(theta)*f
        fy = math.sin(theta)*f

        return fx,fy

def loop(system, saturnV, mars):
    timestep = 1*24*3600
    elapsedRuns = 0
    while enroute(saturnV, mars) and not offScreen(saturnV) and elapsedRuns < 500:      #stop the simulation if the spaceship reaches Mars, if the spaceship goes too far off screen, or if the simulation takes too long
        saturnV.goto(saturnV.xloc*SCALE, saturnV.yloc*SCALE)
        saturnV.pendown()
        total_fx = total_fy = 0.0
        for body in system:
            body.goto(body.xloc*SCALE, body.yloc*SCALE)
            fx, fy = saturnV.attraction(body)
            total_fx += fx
            total_fy += fy

            saturnV.setXVel(saturnV.getXVel()+((total_fx*timestep)/saturnV.mass))
            saturnV.setYVel(saturnV.getYVel()+((total_fy*timestep)/saturnV.mass))

            saturnV.xloc += saturnV.getXVel()*timestep
            saturnV.yloc += saturnV.getYVel()*timestep
            saturnV.goto(saturnV.xloc*SCALE, saturnV.yloc*SCALE)

            elapsedRuns += 1


def enroute(saturnV, mars):
    spaceshipCoord = (saturnV.getX(),saturnV.getY())
    print(spaceshipCoord)
    endingCoord = (mars.getX(), mars.getY())
    rx = spaceshipCoord[0] - endingCoord[0]
    ry = spaceshipCoord[1] - endingCoord[1]
    r = math.sqrt(rx**2+ry**2)*SCALE

    if (r < 50):
        return False
    return True
def offScreen(saturnV):
    spaceshipCoord = (saturnV.getX(),saturnV.getY())
    if abs(spaceshipCoord[0]>4*10**11 or abs(spaceshipCoord[1])>4*10**11):
        return True
    else:
        return False

def main():
    askingUser = True
    while(askingUser):
        spaceshipAngleInput = input(
            "At what angle (in degrees) do you want the rocket to leave Earth? (Directly to the right is 0, and the angles increase counterclockwise from there): ")
        try:
            spaceshipAngle = int(spaceshipAngleInput)
        except:
            print("I couldn't understand your input. Please be sure to input a number.")
        else:
            askingUser = False
    askingUser = True
    while(askingUser):
        spaceshipVelocityInput = input(
            "At what velocity do you want the rocket to leave Earth? (Hint: a velocity of at least 4000 m/s is rquired to escape Earth) ")
        try:
            spaceshipVelocity = int(spaceshipVelocityInput)
        except:
            print("I couldn't understand your input. Please be sure to input a number.")
        else:
            if (spaceshipVelocity<4000):
                print("The velocity wasn't fast enough, please try inputting a faster velocity.")
            else:
                askingUser = False

    turtle.setup(800, 800)          # Set the window size to 800 by 800 pixels
    turtle.bgcolor("black")         # Set up the window with a black background

    sun = planet()
    sun.name = 'Sun'
    sun.mass = 1.98892 * 10**30
    sun.penup()
    sun.color('yellow')
    sun.shape('circle')
    sun.diameter = 1.3914 * 10**6
    sun.shapesize(5.0,5.0,1)

    mercury = planet()
    mercury.name = 'Mercury'
    mercury.mass = 3.302 * 10**23
    mercury.penup()
    mercury.color('gray')
    mercury.shape('circle')
    mercury.shapesize(0.5,0.5,1)
    mercury.diameter = 2440*2
    mercury.yloc = (1 * AU) *  -0.4608453269808703
    mercury.xloc = (1 * AU) *  -0.06333487572394930
    mercury.vy = AU * -0.002399853089908365  / 86400
    mercury.vx = AU * 0.02222816779156590 / 86400

    venus = planet()
    venus.name = 'Venus'
    venus.mass = 48.685 * 10**23
    venus.penup()
    venus.color('Khaki')
    venus.shape('circle')
    venus.shapesize(1.5,1.5,1)
    venus.diameter = 6051.893*2
    venus.yloc = (0.7262658 * AU) *  0.0525483
    venus.xloc = (0.7262658 * AU) *  0.7232002
    venus.vy = AU * 0.0200813  / 86400

    earth = planet()
    earth.name = 'Earth'
    earth.mass = MASS_EARTH
    earth.penup()
    earth.color('green')
    earth.shape('circle')
    earth.shapesize(1.5,1.5,1)
    earth.diameter = 12742
    earth.yloc = (1 * AU) *   0.96756
    earth.xloc = (1 * AU) *  -0.17522

    mars = planet()
    mars.name = 'Mars'
    mars.mass = MASS_MARS
    mars.penup()
    mars.color('red')
    mars.shape('circle')
    mars.shapesize(1.25,1.25,1)
    mars.diameter = 3389.92*2
    mars.yloc = (1 * AU) *  -0.857574644771996
    mars.xloc = (1 * AU) *  -1.320107604952232

    saturnV = spaceship(earth.xloc, earth.yloc, 200.0, float(spaceshipAngle), float(spaceshipVelocity))      #(x location, y location, altitude, angle, velocity)
    saturnV.name = "Saturn V"
    saturnV.penup()
    saturnV.shape('classic')
    saturnV.color('pink')
    saturnV.shapesize(0.5,0.5,1)
    saturnV.yloc = (1.01 * AU) *   0.96756
    saturnV.xloc = (1.01 * AU) *  -0.17522

    loop([sun, earth, mars, venus, mercury],saturnV,mars)


if __name__ == '__main__':          # The code starts here
    main()                          # Goes to the function called main (line 82)
