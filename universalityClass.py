"""
universalityClass.py - verifies the universality class of layers created via
ballistic deposition.

INPUTS:

KPZMatrix: matrix used in creation of layers (matrix)
blockNumber: number of blocks used in one ballistic deposition layer (int)
squareSize: dimensions of simulation area (int)
existingLayers: Number of layers in simulation (int)
diffProb: Probability of diffusion during ballistic deposition (float)
tempProb: Probability of surface normal deposition in layering (float)

OUTPUTS:

KPZMatrix: matrix used in creation of layers (matrix)
existingLayers: Number of layers in simulation (int)
"""

import math
import matplotlib.pyplot as plt
from matplotlib import colors
import random
import numpy as np

def addLayer(blockNumber, squareSize, diffProb):
    # Matrix to hold values of <h>
    hValsMatrix = []
    # Matrix to hold value of <h**2>
    hSquaredValsMatrix = []
    # Time List for plotting
    timeList = np.arange(0, blockNumber, 1)
    # Average data over trials
    for i in range (100):
        print(i)
        # Values for sum(h) at each time
        hVals = []
        # Values for sum(h**2) at each time
        hSquaredVals = []
        # Initialise empty matrix for deposition
        KPZMatrix = np.zeros((10000, squareSize))
        # blockNumber of blocks in each trial
        squareHeight = 100
        for j in range(blockNumber):
            # Boolean for block finished deposition           
            finished = False
            # Particle starts in random position at top of matrix
            location = [random.randint(0, squareSize - 1), squareHeight]
            # Block starts at top of matrix
            KPZMatrix[location[1]][location[0]] = 4                
            # While not finished - particle continues to fall
            while not finished:
                # Define sticky neighbours of particle
                if location[0] == 0:
                    neighbourUp = KPZMatrix[(location[1] - 1)][location[0]]
                    neighbourRight = KPZMatrix[location[1]][location[0] + 1]
                    # Particle has neighbour or is at bottom of surface
                    if neighbourUp == 4 or neighbourRight == 4 or location[1] == 0:
                        # Stop ballistic deposition process
                        finished = True
                    else:
                    # Continue falling
                        KPZMatrix[location[1] - 1][location[0]] = 4
                        KPZMatrix[location[1]][location[0]] = 0
                        location[1] -= 1
                    
                elif location[0] == squareSize - 1:
                    neighbourUp = KPZMatrix[(location[1] - 1)][location[0]]           
                    neighbourLeft = KPZMatrix[location[1]][location[0] - 1]
                    # Particle has neighbour or is at bottom of surface
                    if neighbourUp == 4 or neighbourLeft == 4 or location[1] == 0:
                        # Stop ballistic deposition process
                        finished = True
                    else:
                    # Continue falling
                        KPZMatrix[location[1] - 1][location[0]] = 4
                        KPZMatrix[location[1]][location[0]] = 0
                        location[1] -= 1
                    
                elif location[0] > 0 and location[0] < squareSize - 1:
                    neighbourUp = KPZMatrix[(location[1] - 1)][location[0]]           
                    neighbourLeft = KPZMatrix[location[1]][location[0] - 1]
                    neighbourRight = KPZMatrix[location[1]][location[0] + 1]
                    # Particle has neighbour or is at bottom of surface
                    if neighbourUp == 4 or neighbourLeft == 4 or neighbourRight == 4 or location[1] == 0:
                        # Stop ballistic deposition process
                        finished = True
                    else:
                        # Continue falling
                        KPZMatrix[location[1] - 1][location[0]] = 4
                        KPZMatrix[location[1]][location[0]] = 0
                        location[1] -= 1

            
            # List of heights of each column for time 
            heightList = []
            # Calculate height in each column 
            for l in range(0, squareSize):
                addedFour = False
                for k in range (squareHeight, 0, -1):
                    if KPZMatrix[k][l] == 4:
                        heightList.append(k)
                        addedFour = True
                        break
                    
                if addedFour == False:
                    heightList.append(0)

            # List of heights squared
            heightSquaredList = [q**2 for q in heightList]
            hVals.append(sum(heightList))
            hSquaredVals.append(sum(heightSquaredList))
            # Update height of release
            squareHeight = (max(heightList) + 100)
        # Append data lists to matrices for averaging
        hValsMatrix.append(hVals)
        hSquaredValsMatrix.append(hSquaredVals)

    wAvVals = []
    hAvVals = []
    hAvSquaredVals = []

    for i in (np.array(hValsMatrix)).T:
        hAvVals.append((np.mean(i))/squareSize)

    for i in (np.array(hSquaredValsMatrix)).T:
        hAvSquaredVals.append((np.mean(i))/squareSize)

    """
    # Lists for average
    wAvVals = []
    hAvVals = []
    hAvSquaredVals = []
    # Transpose matrices and averages
    for i in [list(e) for e in zip(*hValsMatrix)]:
        hAvVals.append(sum(i)/((squareSize)*len(i)))
    for i in [list(e) for e in zip(*hSquaredValsMatrix)]:
        hAvSquaredVals.append(sum(i)/((squareSize)*len(i)))
    # Width of KPZ surface
    """
    
    for i in range(len(hAvSquaredVals)):
        wAvVals.append(math.sqrt((hAvSquaredVals[i] - (hAvVals[i])**2)))
    # Scatter plot of width against time
    plt.scatter(timeList, wAvVals)
    plt.title("(W(t)) against time for ballistic deposition")
    plt.xlabel("time / t")
    plt.ylabel("W(t)")
    plt.show()

    with open('universalityData.dat', 'w+', ) as dataFile:
        for i in range(len(wAvVals)):
            dataFile.write('%lf, %lf\n' % (timeList[i], wAvVals[i]))

addLayer(80000, 100, 0)
        
            
        
                        
                            
        
