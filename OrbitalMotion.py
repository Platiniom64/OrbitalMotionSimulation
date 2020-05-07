"""
this program is a simulation of the inner planets of our solar system (namely the sun, Mercury,
Venus, Earth and Mars). The planets are objects of the class Planet which enables this class (solarSystemAnimation)
to animate them. The information of the planets can be found in the file PropertiesPlanets.
"""

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

from PlanetClass import Planet


class SolarSystemAnimation(object):
    """ this class creates an animation using any planets found in the file PropertiesPlanets.
    It then creates a simulation using these planets. There are several plots to display different things
    such as the total energy of the system, or their orbit. """

    def __init__(self):

        # these two lists will contain the planets and their patches respectively
        self.listPlanets = []
        self.listPatches = []

        # this gets the information from the file
        self.getInfoFile()

        # this is some parameters for the simulation, try to keep the timeStep small for more accuracy
        self.interval = 0
        self.nframes = 9999999999
        self.timeStep = 70000
     
        # this list contains the total energy for each iteration
        self.listTotalEnergy = []
    
    def getInfoFile(self):
        """ this method gets all the information in the file about the planets """
        
        self.f = open("PropertiesPlanets.txt", "r")

        # first we ignore the first paragraph where we explain what the file is for.
        line = self.f.readline()
        while line is not "\n":
            line = self.f.readline()

        # this gets the information about the planets
        self.getInformationPlanets()

        self.f.close()

    def getInformationPlanets(self):
        """ This method gets all the information about the planets from the file. It also
        adds a satellite with custrom parameters (but not from the file). """

        line = self.f.readline()
        while line is not "\n":
            try:
                planetName = line[:line.index("\n")]
                planetMass = self.getNumber()
                planetInitPosx = self.getNumber()
                planetInitPosy = self.getNumber()
                planetInitPos = np.array([planetInitPosx, planetInitPosy])
                planetInitVelx = self.getNumber()
                planetInitVely = self.getNumber()
                planetInitVel = np.array([planetInitVelx, planetInitVely])
                planetRadius = self.getNumber()
                planetColour = self.getString()
            except:
                print("There was a problem while getting information from the file.")
                quit()


            planet = Planet(planetName, planetMass, planetInitPos, planetInitVel, planetRadius, planetColour)
            self.listPlanets.append(planet)
            self.listPatches.append(planet.patch)

            line = self.f.readline()

        # we also include a satellite that we implement here
        satName = "Satellite"
        satMass = 500000
        satInitPos = np.array([1.5e+11+100, 0])
        satInitVelx = 11500
        satInitVely = -800
        satInitVel = np.array([satInitVelx, satInitVely+29748.485675745])
        satRadius = 928e+6
        satColour = "#000000"

        Satellite = Planet(satName, satMass, satInitPos, satInitVel, satRadius, satColour)
        self.listPlanets.append(Satellite)
        self.listPatches.append(Satellite.patch)
        
    def getNumber(self):
        """ This is a helper method that reads a line from the file and removes
        everythin before the : as that is simply the name of the variable in the file.
        it returns the value as a float. """

        line = self.f.readline()
        # here we convert the line into a float where we remove the characters begore the colon
        number = float(line[line.index(':')+1:])
        
        return number

    def getString(self):
        """ This is a helper method that reads a line from the file and removes
        everythin before the : as that is simply the name of the variable in the file.
        it returns the value as a string."""

        line = self.f.readline()
        string = line[line.index(':')+1:line.index('\n')]

        return string



    def calculateTotalEnergy(self):
        """ this method calculates the total energy of the system by adding
        the kinetic and potential energies together. """

        self.potEnergy = Planet.calculatePotentialEnergy(self.listPlanets)
        self.sumKineticEnergy = Planet.calculateTotalKineticEnergy(self.listPlanets)

        self.totalEnergy = self.sumKineticEnergy + self.potEnergy

    def updateTotalEnergyPlot(self):
        """ this method updates the plot of the total energy. """

        # y values of the graph
        self.calculateTotalEnergy()
        self.listTotalEnergy.append(self.totalEnergy)

        # updates the plot
        self.totalEnergyPlot.clear()
        self.totalEnergyPlot.title.set_text("Total Energy of system over time")
        self.totalEnergyPlot.plot(self.listTotalEnergy)

    def printTotalEnergyToFile(self, i):
        """ this method prints the total energy of the system every nth iteration. """

        # reduce the frequency at which the info is written to the file
        n = 50
        if (i % n == 0):
            self.f = open("TotalEnergy.txt", "a+")

            self.f.write("Total energy of system: " + str(self.totalEnergy) + "\n")

            self.f.close()

    def checkDistanceSatMars(self, i):
        """ This method checks whether the satelite is close to Mars. if so, it prints the time of the
        journey in the legend of the plot "traces orbit". """

        for planet in self.listPlanets:
            if planet.name == "Satellite":
                if (planet.checkDistMars(i, self.timeStep, self.listPlanets)):
                    for i in range(len(self.textLegendOrbit)):
                        if self.textLegendOrbit[i] == planet.name:
                            self.textLegendOrbit[i] = planet.name + " time to Mars: " + str(round(planet.distanceToMars, 7))
                            self.tracesOrbitPlot.legend(self.tracesOrbitPlot.lines[1:], self.textLegendOrbit, loc='lower left', bbox_to_anchor=(0.0, -0.6))
        
    def checkOrbitPlanets(self, i):
        """ This method checks if the planet has gone around then sun. If so it displays
        the time it took for that planet to go around the sun in the legend of the plot "traces orbit". """
        
        for planet in self.listPlanets:
            if (planet.name != "Sun"):
                if (planet.checkOrbitalPeriod(i, self.timeStep, self.listPlanets)):
                        for i in range(len(self.textLegendOrbit)):
                            if self.textLegendOrbit[i] == planet.name:
                                self.textLegendOrbit[i] = planet.name + " orbit: " + str(round(planet.orbitalPeriod, 7))
                                self.tracesOrbitPlot.legend(self.tracesOrbitPlot.lines[1:], self.textLegendOrbit, loc='lower left', bbox_to_anchor=(0.0, -0.6))
    
    def updateDisplays(self, i):
        """ this method updates the figures on the animation and everything related to
        the animation. """

        self.updateTotalEnergyPlot()

        self.printTotalEnergyToFile(i)

        # this plots the trace of the orbit for each planet
        self.tracesOrbitPlot.lines = []
        for planet in self.listPlanets:
            planet.trail = self.tracesOrbitPlot.plot(planet.positionsx, planet.positionsy, color=planet.colour, linewidth=0.5)
        
        self.checkDistanceSatMars(i)

        self.checkOrbitPlanets(i)   
            
        self.fig.canvas.draw()


    def stepForwardSimulation(self):
        """ this method will make a step forward in the animation by appliying one step in 
        the Beeman scheme. """
        
        # we first move the planets to their positions
        for i in range(len(self.listPatches)):
            self.listPatches[i].center = self.listPlanets[i].pos


        # * then we calculate the postitions of te planets for the next iteration
        for planet in self.listPlanets:
            planet.calculatePOS(self.timeStep)

        # * then we calculate the new acceleration of the planets according to their new position
        for planet in self.listPlanets:
            planet.calculateACC(self.listPlanets)

        # * then we can calculate the velocity of each planet for the next iteration
        for planet in self.listPlanets:
            planet.calculateVEL(self.timeStep)
        
        # * finally we update the values of acceleration for each planet
        for planet in self.listPlanets:
            planet.updateACC()



    def animate(self, i):
        """ This function is the one that is executed every time a frame is redrawn
        so it calls the method to move the planets and the method that update the plots. """

        self.updateDisplays(i)

        self.stepForwardSimulation()

        return self.listPatches

    def initAnim(self):
        """ This method is executed before the animation start, so it adds the patches
        of the planets to the axes. """

        for patch in self.listPatches:
            self.simulationPlot.add_patch(patch)

        return self.listPatches

    def run(self):
        """ This method launches the animation, it is called outside the class to start the
        animation. """

        # we first create the plot and axes
        self.fig = plt.figure(figsize=(7, 7))

        self.simulationPlot = plt.subplot(2,2,1)
        self.tracesOrbitPlot = plt.subplot(2,2,2)
        self.totalEnergyPlot = plt.subplot(2,2,3)
        
        # set up the axes for simulationPlot
        self.simulationPlot.axis('scaled')
        self.simulationPlot.title.set_text("Simulation of planets")
        maxOrbitalR = 0
        for planet in self.listPlanets:
            if planet.pos[0] > maxOrbitalR:
                maxOrbitalR = planet.pos[0]
        scaleUp = 1.1 * maxOrbitalR
        self.simulationPlot.set_xlim(-scaleUp, scaleUp)
        self.simulationPlot.set_ylim(-scaleUp, scaleUp)

        # set up the axes for tracesOrbitPlot
        self.tracesOrbitPlot.axis('scaled')
        self.tracesOrbitPlot.title.set_text("Traces of planets in orbit")
        self.tracesOrbitPlot.set_xlim(-scaleUp, scaleUp)
        self.tracesOrbitPlot.set_ylim(-scaleUp, scaleUp)
        for planet in self.listPlanets:
            planet.trail = self.tracesOrbitPlot.plot(planet.positionsx, planet.positionsy, color=planet.colour, linewidth=0.5, label='orbit of ' + planet.name)
        self.textLegendOrbit = []
        for planet in self.listPlanets:
            if planet.name != "Sun":
                self.textLegendOrbit.append(planet.name)
        self.tracesOrbitPlot.legend(self.tracesOrbitPlot.lines[1:], self.textLegendOrbit, loc='lower left', bbox_to_anchor=(0.0, -0.6))
        
        # create the animator
        FuncAnimation(self.fig, self.animate, init_func = self.initAnim, frames = self.nframes, repeat = False, interval = self.interval, blit = True)

        # show the plot
        plt.show()
        

def main():
    anim = SolarSystemAnimation()
    anim.run()

main() 