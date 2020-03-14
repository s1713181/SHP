"""
This function finds the fractal dimensionality of a cluster generated
by the DLACluster.py script.

Only to be used for 1 seed particle in a cluster.

INPUTS:

NONE

OUTPUTS:

NONE
"""

from DLAcluster import DLAcluster 
import numpy
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Average mass values
massAv = []
# Average radius values
radiusAv = []
# Average anastomosis values
anastAv = []
# Average fractal dimensionality values
dimensionalityAv = []
# Errors in average radii
radiusError = []
# Errors in anastomosis values
anastError = []
# List of layering increments
layerList = numpy.arange(5, 45, 2)

for j in range(1,11):
    # Radius values
    radiusArray = []
    # Mass values
    massArray = []
    # Anastomosis values
    anastArray = []

    for i in range (5):
        print("trial " + str(i + 1) + " out of 5")
        massValue, radiusValue, matrix, islands = DLAcluster(200, False, 1000, j*5, 0, 1, 1) 
        massArray.append(massValue)
        radiusArray.append(radiusValue)
        anastArray.append(islands/massValue)

    logRadius = numpy.log(radiusArray)
    logMass = numpy.log(massArray)

    massAv.append(numpy.mean(massArray))
    anastAv.append(numpy.mean(anastArray) - 1)
    anastError.append(numpy.std(anastArray)/numpy.size(anastArray))
    radiusError.append(numpy.std(radiusArray)/numpy.size(radiusArray))

    """
    # Fit a log function using numpy polyfit
    fitLog = numpy.polyfit(logRadiusAv, logMassAv,1)
    fitLogFunc = numpy.poly1d(fitLog)

    # Print out results
    print("Parameters for the log fit: slope = ",fitLog[0],"shift: ",fitLog[1])
    print("Parameters from the log fit: form is e^",fitLog[1],"*r^",fitLog[0])
    num=str(numpy.round(fitLog[0],3))
    dimensionalityAv.append(num)

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

# Plot of anastomosis against layering increment
plt.errorbar(layerList, anastAv, yerr=anastError)
plt.title("Anastomosis per unit mass against layering increment")
plt.xlabel("Particles added between layers")
plt.ylabel("Anastomosis per unit mass")
plt.show()

"""
plt.plot(layerList, dimensionalityAv)
plt.title("Fractal dimensionality against layer increment")
plt.xlabel("Particles added between layers")
plt.ylabel("Fractal dimensionality of resulting cluster")
plt.show()
"""


# Plot of fractal dimensionality against layering increment

# Plot of fractal dimensionality against temperature

# Plot of anastomosis against temperature

# Plot of fractal dimensionality against crystalographic bias

# Plot of cluster-covering probability against fractal dimensionality 
                





