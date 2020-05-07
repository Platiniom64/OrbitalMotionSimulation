"""
This file contains the class Planet that can be used for animation using matplotlib.
"""

import numpy as np
from numpy.linalg import norm
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


class Planet(object):
    """ this object was specifically made for eht MarsPhobosAnimation clas but is very general
    it contains all the information of a planet and has some functions such as finding the distance between
    planets and getting their force between them. """

    G = 6.674 * (10 ** -11) # the units are Newtons kg^-2 m^2 so the value from the file must correspond to these

    def __init__(self, name, mass, initPos, initVel, radius, colour):

        self.name = name
        self.mass = mass
        self.pos = initPos
        self.vel = initVel
        self.rad = radius
        self.colour = colour

        self.previousAcc = 0
        self.currentAcc = 0 # we intitialise it to zero as the acceleration is zero at time 0

        self.trail = 0
        self.positionsx = []
        self.positionsx.append(self.pos[0])
        self.positionsx.append(self.pos[0])
        self.positionsy = []
        self.positionsy.append(self.pos[1])
        self.positionsy.append(self.pos[1])
        
        # we also reate a patch for that planet
        self.patch = patches.Circle(self.pos, self.rad, color = self.colour, animated = True)

        self.distanceToMars = 0 # for the satelite
        self.orbitalPeriod = 0 # for the planets

    @classmethod
    def calculatePotentialEnergy(Planet, listPlanets):
        """ this class method calclates the total potential energy of all the planets
        in a system. """

        partialSumPotEnergy = 0
        
        for planet1 in listPlanets:
            for planet2 in listPlanets:

                if planet1 is not planet2:
                    partialSumPotEnergy += planet1.potentialEnergy(planet2)

        return partialSumPotEnergy / 2
    
    def potentialEnergy(self, planet2):
        """ this method calculates the potential energy between two planets """
        
        potentialEnergy = - (Planet.G * self.mass * planet2.mass)/(norm(self.distance(planet2)))
        return potentialEnergy

    @classmethod
    def calculateTotalKineticEnergy(Planet, listPlanets):
        """ this class method calculates the total kinetic energy of a system. """

        sumEnergy = 0
        for planet in listPlanets:
            sumEnergy += planet.kineticEngergy()
        return sumEnergy

    def kineticEngergy(self):
        """ this method returns the kinetic energy of the planet. """
        
        kineticEnergy = (1/2) * self.mass * (norm(self.vel)**2)
        
        return kineticEnergy

        

    def distance(self, Planet2):
        """ this method returns the distance from the planet self to the planet2.
        it returns a verctor so calling the function on planet2 to self would have an opposite direction. """

        distanceToPlatet = Planet2.pos - self.pos
        
        return distanceToPlatet

    def gravitationalForce(self, Planet2):
        """ this method returns the gravitational force between the planet self and the planet2.
        it returns a verctor so calling the function on planet2 to self would have an opposite direction. """
        
        # we need the distance between the two planets
        distance = self.distance(Planet2)
        # then we calculate the force
        force = Planet.G * (self.mass * Planet2.mass)/(norm(distance)**3) * distance

        return force


    def calculateACC(self, listPlanets):
        """ this method calculates the acceleration of the planet due to the other planets. Hence
        the method receives the list of planets as argument """

        sumFORCE = np.array([0,0])

        for planet in listPlanets:
            # we ignore the case that the planet we called the method on get to calculate the force between itself
            if self != planet:
                sumFORCE = sumFORCE + self.gravitationalForce(planet)

        # renews the value of the acceleration of the planet
        self.newAcc = sumFORCE / self.mass

    def calculatePOS(self, timeStep):
        """ this method updates the position of the planet based on its speed """
        self.pos = self.pos + (self.vel * timeStep) + (1/6)*(4*self.currentAcc-self.previousAcc)*(timeStep**2)

        
        self.positionsx.append(self.pos[0])
        self.positionsy.append(self.pos[1])

    def calculateVEL(self, timeStep):
        """ This method updates the velocity of the planet """

        self.vel = self.vel + (1/6)*(2*self.newAcc + 5*self.currentAcc - self.previousAcc) * timeStep

    def updateACC(self):
        """ this method simply changes the value of the current and previous acceleration.
        This is required for the Beeman scheme to be applied properly. """
        
        self.previousAcc = self.currentAcc
        self.currentAcc = self.newAcc



    def checkDistMars(self, i, timeStep, listPlanets):
        """ this method checks if the satellite is close to Mars, if so return true. """

        for planet in listPlanets:
                if (planet.name == "Mars"):
                    if (self.distance(planet)[0] < 999999999 and self.distance(planet)[1] < 999999999):
                        self.distanceToMars = i * timeStep # seconds
                        self.distanceToMars = self.distanceToMars / 31556736 # in Earth years
                        return True
                    else:
                        return False

    def checkOrbitalPeriod(self, i, timeStep, listPlanets):
        """ this method checks if the planet has gone around the sun, if so return true """

        if (self.pos[0] > 0 and self.positionsy[-2] < 0 and self.pos[1] > 0):
            self.orbitalPeriod = i * timeStep # seconds
            self.orbitalPeriod = self.orbitalPeriod / 31556736 # in Earth years
            return True
        else:
            return False