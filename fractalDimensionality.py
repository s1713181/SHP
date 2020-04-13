"""
This function finds the fractal dimensionality and anastomosis of a cluster generated
by the DLACluster.py script, and produces relevant plots.

Only to be used for 1 seed particle in a cluster.

INPUTS:

NONE

OUTPUTS:

NONE
"""

from DLAcluster import DLAcluster 
import numpy
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import random

def main():

    # Average mass values
    massAv = []
    # Average radius values
    radiusAv = []
    # Average anastomosis values
    anastAv = []
    # Average anastomosis per unit mass values
    anastPerAv = []
    # Error on anastomosis per unit mass
    anastPerError = []
    # Average fractal dimensionality values
    dimensionalityAv = []
    # Error on fractal dimensionality
    dimensionalityError = []
    # Errors in average radii
    radiusError = []
    # Errors in anastomosis values
    anastError = []
    # List of layering increments
    layerList = numpy.arange(15, 105, 5)
    # Density average
    densityAv = []
    # Error on density
    densityError = []

    for j in range(3, 21):
        print(j)
        # Radius values
        radiusArray = []
        # Mass values
        massArray = []
        # Anastomosis values
        anastArray = []
        # Anastomosis per unit mass values
        anastPerArray = []
        # Calculate fractal dimensionalities
        dimList = []
        # Density List
        densityList = []

        for i in range (10):
            print("trial " + str(i + 1) + " out of 10")
            massValue, radiusValue, clusterArea, matrix, islands = DLAcluster(200, False, 1000, j*5, 1, 1, 1) 
            massArray.append(massValue)
            radiusArray.append(radiusValue)
            anastPerArray.append((islands - 1)/massValue)
            anastArray.append((islands - 1))
            dimList.append(numpy.log(massValue)/numpy.log(radiusValue))
            densityList.append((massValue)/(clusterArea))

        logRadius = numpy.log(radiusArray)
        logMass = numpy.log(massArray)

        massAv.append(numpy.mean(massArray))

        anastAv.append(numpy.mean(anastArray))

        anastPerAv.append(numpy.mean(anastPerArray))

        anastPerError.append(bootstrap(anastPerArray, 100))

        anastError.append(bootstrap(anastArray, 100))

        radiusError.append(bootstrap(radiusArray, 100))

        dimensionalityAv.append(numpy.mean(dimList))

        dimensionalityError.append(bootstrap(dimList, 100))

        densityAv.append(numpy.mean(densityList))

        densityError.append(bootstrap(densityList, 100))

        print("Bootstrap completed")


    """
    # Plot results
    fig=plt.subplot()
    plt.scatter(logRadiusAv,logMassAv, color='tomato', edgecolors='tomato', s=30)
    plt.plot(logRadiusAv, fitLogFunc(logRadius),color='dodgerblue', lw=3)
    plt.title("Log-log plot, mass vs radius",fontsize=20)
    plt.xlabel("Log radius",fontsize=15)
    plt.ylabel("Log mass",fontsize=15)
    plt.grid(True)
    fig.text(2.6,4.3,'fractal dimensionality:'+num)
    fig.spines["top"].set_visible(False)  
    fig.spines["right"].set_visible(False)  
    plt.savefig('logRadiusMass.png')
    plt.show()
    """
    """

    # Plot of anastomosis against layering increment
    plt.errorbar(layerList, anastAv, yerr=anastError)
    plt.title("Anastomosis against layering increment for KPZ Layering")
    plt.xlabel("Particles added between KPZ layers")
    plt.ylabel("Anastomosis in cluster")
    plt.show()

    with open('anastKPZData.dat', 'w+', ) as dataFile:
        for i in range (len(anastAv)):
            dataFile.write('%lf, %lf, %lf\n' % (layerList[i], anastAv[i], anastError[i]))


    plt.plot(layerList, dimensionalityAv)
    plt.title("Fractal dimensionality against layer increment")
    plt.xlabel("Particles added between layers")
    plt.ylabel("Fractal dimensionality of resulting cluster")
    plt.show()
    """

    with open('anastSNDData.dat', 'w+', ) as dataFile:
        for i in range (len(dimensionalityAv)):
            dataFile.write('%lf, %lf, %lf\n' % (layerList[i], anastAv[i], anastError[i]))

    with open('dimensionalitySNDData.dat', 'w+', ) as dataFile:
        for i in range (len(dimensionalityAv)):
            dataFile.write('%lf, %lf, %lf\n' % (layerList[i], dimensionalityAv[i], dimensionalityError[i]))

    with open('anastPerSNDData.dat', 'w+', ) as dataFile:
        for i in range (len(dimensionalityAv)):
            dataFile.write('%lf, %lf, %lf\n' % (layerList[i], anastPerAv[i], anastPerError[i]))

    with open('densitySNDData.dat', 'w+', ) as dataFile:
        for i in range (len(dimensionalityAv)):
            dataFile.write('%lf, %lf, %lf\n' % (layerList[i], densityAv[i], densityError[i]))


    plt.errorbar(layerList, dimensionalityAv, yerr=dimensionalityError, linestyle='None')
    plt.title("Fractal dimensionality against rate of layering (SND)")
    plt.xlabel("Particles added to cluster between layering")
    plt.ylabel("Fractal dimensionality of resulting cluster")
    plt.show()

    plt.errorbar(layerList, densityAv, yerr=densityError, linestyle='None')
    plt.title("Mass Density against rate of layering (SND)")
    plt.xlabel("Particles added to cluster between layering")
    plt.ylabel("Mass Density of resulting cluster")
    plt.show()

    plt.errorbar(layerList, anastPerAv, yerr=anastPerError, linestyle='None')
    plt.title("Anastomosis per unit mass against rate of layering (SND)")
    plt.xlabel("Particles added to cluster between layering")
    plt.ylabel("Islands in cluster per unit mass")
    plt.show()

    plt.errorbar(layerList, anastAv, yerr=anastError, linestyle='None')
    plt.title("Anastomosis against rate of layering (SND)")
    plt.xlabel("Particles added to cluster between layering")
    plt.ylabel("Islands in cluster")
    plt.show()


def bootstrap(list, bootstrapTrials):
    # Initialise list to contain c values for each list of e values
    meanList = []
    # User specifies number of trials
    for i in range (bootstrapTrials):
        # Initialise list for n selected energy values in each trial
        bootStrapSublist = []
        for j in range (len(list)):
            # Randomly select element in list of energy values
            k = random.randint(0, len(list) - 1)
            # Append random value from energy list n times
            bootStrapSublist.append(list[k])
        # Calculate value of c from list of n random energy values
        mean = sum(bootStrapSublist)/len(bootStrapSublist)
        meanList.append(mean)

    # Calculate overall bootstrap error on value of c
    bootStrapMean = sum(meanList)/len(meanList)
    meanSquared = bootStrapMean**2
    squaredBootList = [p**2 for p in meanList]
    squaredMean = sum(squaredBootList)/len(squaredBootList)

    bootStrapError = math.sqrt(squaredMean - meanSquared)

    return (bootStrapError)

main()



                





