"""
addLayer.py - adds a KPZ layer of silica solution to the simulation area, and
replaces the previously added layer of silica solution with a layer of
silica solid. After the addition of one KPZ layer, layering is a stochastic
process governed by temperature (tempProb), in which layers can either be KPZ
layers from ballistic depostion, or a surface normal deposition style 'smoothing'
layer to reduce roughness. 

INPUTS:

KPZMatrix: matrix used in creation of layers (matrix)
matrix: matrix representing simulation area (matrix)
blockNumber: number of blocks used in one ballistic deposition layer (int)
squareSize: dimensions of simulation area (int)
existingLayers: Number of layers in simulation (int)
diffProb: Probability of diffusion during ballistic deposition (float)
tempProb: Probability of surface normal deposition in layering (float)
depMod: Modulating factor in surface-normal deposition (float)
clusterMod: Modulating factor in on-cluster deposition (float)

OUTPUTS:

KPZMatrix: matrix used in creation of layers (matrix)
existingLayers: Number of layers in simulation (int)
"""

### Silica Solution == 4
### Solid Silica Layer == 3
### Solid Silica Layer == 2
### Metal Oxide Particle == 1
### Silica gel == 0

import math
import matplotlib.pyplot as plt
from matplotlib import colors
import random
import numpy as np

# Add solid silica layer to simulation (ballistic deposition or surface normal deposition) 
def addLayer(KPZMatrix, matrix, blockNumber, squareSize, existingLayers, tempProb, depMod, clusterMod):

    
    # If layers exist, remove additional area added for clustering
    if existingLayers > 0:
        # Remove additional area for clustering
        for i in range(squareSize):
            for j in range(squareSize - 1, -1, -1):
                # Highest up site of gel in column
                if KPZMatrix[j][i] == 4:
                    marker = j
                    # Remove 10 blocks of silica gel 
                    for k in range(marker, marker - 10, -1):
                        if KPZMatrix[k][i] == 4:
                            # Site is silica solution
                            KPZMatrix[k][i] = 0
                    break
    
    
    
    # Solution from previous layer solidifies (silica gel becomes solid silica)
    for i in range(squareSize):
        for j in range(squareSize):
            # Alternate layer values for layer pattern in visualisation
            if KPZMatrix[i][j] == 4 and existingLayers%2 != 0:
                KPZMatrix[i][j] = 2
            elif KPZMatrix[i][j] == 4 and existingLayers%2 == 0:
                KPZMatrix[i][j] = 3


    # First layer must be ballistic deposition
    if existingLayers == 0:
        for i in range(blockNumber):
            # Ballistic deposition process            
            finished = False
            # Particle starts in random position at top of matrix
            location = [random.randint(0, squareSize - 1), squareSize - 1]
            KPZMatrix[location[1]][location[0]] = 4                
            # While not finished - particle continues to fall
            # Particle at LHS of simulation area
            if location[0] == 0:
                while not finished:
                    # Look for neighbours
                    neighbourUp = KPZMatrix[(location[1] - 1)%squareSize][location[0]] 
                    neighbourRight = KPZMatrix[location[1]][(location[0] + 1)]
                    if ((location[1]) == 0) or (neighbourUp or neighbourRight) != (0):
                        finished = True
                    # If no neighbours - continue falling
                    else:
                        KPZMatrix[location[1] - 1][location[0]] = 4
                        KPZMatrix[location[1]][location[0]] = 0
                        location[1] -= 1
            # Particle in RHS of simulation area
            elif location[0] == squareSize - 1:
                while not finished:
                    # Look for neighbours
                    neighbourUp = KPZMatrix[(location[1] - 1)%squareSize][location[0]] 
                    neighbourLeft = KPZMatrix[location[1]][(location[0] - 1)]
                    if ((location[1]) == 0) or (neighbourUp or neighbourLeft) != (0):
                        finished = True
                    # If no neighbours - continue falling
                    else:
                        KPZMatrix[location[1] - 1][location[0]] = 4
                        KPZMatrix[location[1]][location[0]] = 0
                        location[1] -= 1
            # Particle not at RHS or LHS of simulation area
            else:
                while not finished:
                    # Look for neighbours
                    neighbourUp = KPZMatrix[(location[1] - 1)%squareSize][location[0]] 
                    neighbourRight = KPZMatrix[location[1]][(location[0] + 1)]
                    neighbourLeft = KPZMatrix[location[1]][(location[0] - 1)]
                    if ((location[1]) == 0) or (neighbourUp or neighbourRight or neighbourLeft) != (0):
                        finished = True
                    # If no neighbours - continue falling
                    else:
                        KPZMatrix[location[1] - 1][location[0]] = 4
                        KPZMatrix[location[1]][location[0]] = 0
                        location[1] -= 1


    elif existingLayers > 0:
        # Decide between ballistic deposition and surface normal deposition
        marker = random.uniform(0, 1)
        # Ballistic deposition
        if tempProb < marker:
            for i in range(blockNumber):
                # Ballistic deposition process            
                finished = False
                # Particle starts in random position at top of matrix
                location = [random.randint(0, squareSize - 1), squareSize - 1]
                KPZMatrix[location[1]][location[0]] = 4                
                # While not finished - particle continues to fall
                # Particle at LHS of simulation area
                if location[0] == 0:
                    while not finished:
                        # Look for neighbours
                        neighbourUp = KPZMatrix[(location[1] - 1)%squareSize][location[0]] 
                        neighbourRight = KPZMatrix[location[1]][(location[0] + 1)]
                        if ((location[1]) == 0) or (neighbourUp or neighbourRight) != (0):
                            finished = True
                        # If no neighbours - continue falling
                        else:
                            KPZMatrix[location[1] - 1][location[0]] = 4
                            KPZMatrix[location[1]][location[0]] = 0
                            location[1] -= 1
                # Particle in RHS of simulation area
                elif location[0] == squareSize - 1:
                    while not finished:
                        # Look for neighbours
                        neighbourUp = KPZMatrix[(location[1] - 1)%squareSize][location[0]] 
                        neighbourLeft = KPZMatrix[location[1]][(location[0] - 1)]
                        if ((location[1]) == 0) or (neighbourUp or neighbourLeft) != (0):
                            finished = True
                        # If no neighbours - continue falling
                        else:
                            KPZMatrix[location[1] - 1][location[0]] = 4
                            KPZMatrix[location[1]][location[0]] = 0
                            location[1] -= 1
                # Particle not at RHS or LHS of simulation area
                else:
                    while not finished:
                        # Look for neighbours
                        neighbourUp = KPZMatrix[(location[1] - 1)%squareSize][location[0]] 
                        neighbourRight = KPZMatrix[location[1]][(location[0] + 1)]
                        neighbourLeft = KPZMatrix[location[1]][(location[0] - 1)]
                        if ((location[1]) == 0) or (neighbourUp or neighbourRight or neighbourLeft) != (0):
                            finished = True
                        # If no neighbours - continue falling
                        else:
                            KPZMatrix[location[1] - 1][location[0]] = 4
                            KPZMatrix[location[1]][location[0]] = 0
                            location[1] -= 1

        # Surface Normal Deposition (with PBCs)
        elif tempProb >= marker:
            # List for parallel update of layer deposition
            addList = []
            # List for parallel update of cluster deposition
            clusterList = []
            for i in range(1, squareSize - 1):
                for j in range(0, squareSize):
                    if existingLayers%2 == 0:
                        if (KPZMatrix[i][j] == 3):
                            # Stochastic cellular automaton
                            if KPZMatrix[(i + 1)%squareSize][j] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.9717/depMod:
                                    addList.append([(i + 1)%squareSize,j])

                            if KPZMatrix[(i - 1)%squareSize][j] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.9717/depMod: 
                                    addList.append([(i - 1)%squareSize,j]) 

                            if KPZMatrix[i][(j + 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.9717/depMod:
                                    addList.append([i,(j + 1)%squareSize])

                            if KPZMatrix[i][(j  - 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.9717/depMod:
                                    addList.append([i,(j - 1)%squareSize])

                            if KPZMatrix[(i + 1)%squareSize][(j + 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.54544/depMod:
                                    addList.append([(i + 1)%squareSize,(j + 1)%squareSize])

                            if KPZMatrix[(i - 1)%squareSize][(j + 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.54544/depMod:
                                    addList.append([(i - 1)%squareSize,(j + 1)%squareSize])

                            if KPZMatrix[(i + 1)%squareSize][(j - 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.54544/depMod:
                                    addList.append([(i + 1)%squareSize,(j - 1)%squareSize])

                            if KPZMatrix[(i - 1)%squareSize][(j - 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.54544/depMod:
                                    addList.append([(i - 1)%squareSize,(j - 1)%squareSize])

                        
                        
                        elif (matrix[i][j] == 1) and (3 in nearestNeighbours(i, j, KPZMatrix, squareSize) or 6 in nearestNeighbours(i, j, KPZMatrix, squareSize)):
                            if KPZMatrix[(i + 1)%squareSize][j] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.9717/clusterMod:
                                    clusterList.append([(i + 1)%squareSize,j])

                            if KPZMatrix[(i - 1)%squareSize][j] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.9717/clusterMod: 
                                    clusterList.append([(i - 1)%squareSize,j])

                            if KPZMatrix[i][(j + 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.9717/clusterMod:
                                    clusterList.append([i,(j + 1)%squareSize])

                            if KPZMatrix[i][(j  - 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.9717/clusterMod:
                                    clusterList.append([i,(j - 1)%squareSize])

                            if KPZMatrix[(i + 1)%squareSize][(j + 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.54544/clusterMod:
                                    clusterList.append([(i + 1)%squareSize,(j + 1)%squareSize])

                            if KPZMatrix[(i - 1)%squareSize][(j + 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.54544/clusterMod:
                                    clusterList.append([(i - 1)%squareSize,(j + 1)%squareSize])

                            if KPZMatrix[(i + 1)%squareSize][(j - 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.54544/clusterMod:
                                    clusterList.append([(i + 1)%squareSize,(j - 1)%squareSize])

                            if KPZMatrix[(i - 1)%squareSize][(j - 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.54544/clusterMod:
                                    clusterList.append([(i - 1)%squareSize,(j - 1)%squareSize])
                    
                        
                                    
                    elif existingLayers%2 != 0:
                        if (KPZMatrix[i][j] == 2):
                            # Stochastic cellular automaton
                            if KPZMatrix[(i + 1)%squareSize][j] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.9717/depMod:
                                    addList.append([(i + 1)%squareSize,j])

                            if KPZMatrix[(i - 1)%squareSize][j] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.9717/depMod: 
                                    addList.append([(i - 1)%squareSize,j])

                            if KPZMatrix[i][(j + 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.9717/depMod:
                                    addList.append([i,(j + 1)%squareSize])

                            if KPZMatrix[i][(j  - 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.9717/depMod: 
                                    addList.append([i,(j - 1)%squareSize])

                            if KPZMatrix[(i + 1)%squareSize][(j + 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.54544/depMod:
                                    addList.append([(i + 1)%squareSize,(j + 1)%squareSize])

                            if KPZMatrix[(i - 1)%squareSize][(j + 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.54544/depMod:
                                    addList.append([(i - 1)%squareSize,(j + 1)%squareSize])

                            if KPZMatrix[(i + 1)%squareSize][(j - 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.54544/depMod:
                                    addList.append([(i + 1)%squareSize,(j - 1)%squareSize])

                            if KPZMatrix[(i - 1)%squareSize][(j - 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.54544/depMod:
                                    addList.append([(i - 1)%squareSize,(j - 1)%squareSize])

                        
                    
                        
                        elif (matrix[i][j] == 1) and (2 in nearestNeighbours(i, j, KPZMatrix, squareSize) or 6 in nearestNeighbours(i, j, KPZMatrix, squareSize)):
                            if KPZMatrix[(i + 1)%squareSize][j] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.9717/clusterMod:
                                    clusterList.append([(i + 1)%squareSize,j])

                            if KPZMatrix[(i - 1)%squareSize][j] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.9717/clusterMod: 
                                    clusterList.append([(i - 1)%squareSize,j])

                            if KPZMatrix[i][(j + 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.9717/clusterMod:
                                    clusterList.append([i,(j + 1)%squareSize])

                            if KPZMatrix[i][(j  - 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.9717/clusterMod:
                                    clusterList.append([i,(j - 1)%squareSize])

                            if KPZMatrix[(i + 1)%squareSize][(j + 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.54544/clusterMod:
                                    clusterList.append([(i + 1)%squareSize,(j + 1)%squareSize])

                            if KPZMatrix[(i - 1)%squareSize][(j + 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.54544/clusterMod:
                                    clusterList.append([(i - 1)%squareSize,(j + 1)%squareSize])

                            if KPZMatrix[(i + 1)%squareSize][(j - 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.54544/clusterMod:
                                    clusterList.append([(i + 1)%squareSize,(j - 1)%squareSize])

                            if KPZMatrix[(i - 1)%squareSize][(j - 1)%squareSize] == 0:
                                marker = random.uniform(0, 1)
                                if marker < 0.54544/clusterMod:
                                    clusterList.append([(i - 1)%squareSize,(j - 1)%squareSize])

            # Parallel update
            for i in addList:
                KPZMatrix[i[0]][i[1]] = 4
            for i in clusterList:
                KPZMatrix[i[0]][i[1]] = 6
                                                      
                                
    # Add to layer variable     
    existingLayers += 1

    
    # Fill in gaps in KPZ layer   
    for i in range(squareSize):
        for j in range(squareSize - 1, -1, -1):
            if KPZMatrix[j][i] == 4:
                markerMax = j
                break
        for j in range(squareSize - 1, -1, -1):
            if KPZMatrix[j][i] == 2 or KPZMatrix[j][i] == 3:
                markerMin = j
                break

        if existingLayers > 1:
            try:
                for k in range(markerMax, markerMin, -1):
                    if KPZMatrix[k][i] == 0:
                        KPZMatrix[k][i] = 4
            except NameError:
                pass
        else:
            try:
                for k in range(markerMax, 0, -1):
                    if KPZMatrix[k][i] == 0:
                        KPZMatrix[k][i] = 4
            except NameError:
                pass

    

    
    # Additional area for clustering
    for i in range(squareSize):
        for j in range(squareSize - 1, -1, -1):
            # Reached layers
            if KPZMatrix[j][i] == 2 or KPZMatrix[j][i] == 3 or KPZMatrix[j][i] == 4:
                marker = j
                # Fill 10 blocks above layer surface with silica gel
                if marker < squareSize - 15:
                    try:
                        for k in range(marker + 1, marker+11):
                            if KPZMatrix[k][i] == 0:
                                KPZMatrix[k][i] = 4
                    except NameError:
                        for k in range(0, 10):
                            if KPZMatrix[k][i] == 0:
                                KPZMatrix[k][i] = 4
                    break
                else:
                    pass

    


    return (KPZMatrix, existingLayers)

# Obtain nearest neighbours (Moore neighbourhood)
def nearestNeighbours(i, j, KPZMatrix, squareSize):
    nearestNeighbours = []

    nearestNeighbours.append(KPZMatrix[(i + 1)%squareSize][j])
    

    nearestNeighbours.append(KPZMatrix[(i - 1)%squareSize][j]) 

                            
    nearestNeighbours.append(KPZMatrix[i][(j + 1)%squareSize]) 

                            
    nearestNeighbours.append(KPZMatrix[i][(j - 1)%squareSize]) 

                          
    nearestNeighbours.append(KPZMatrix[(i + 1)%squareSize][(j + 1)%squareSize]) 


    nearestNeighbours.append(KPZMatrix[(i - 1)%squareSize][(j + 1)%squareSize])

                            
    nearestNeighbours.append(KPZMatrix[(i + 1)%squareSize][(j - 1)%squareSize])

                          
    nearestNeighbours.append(KPZMatrix[(i - 1)%squareSize][(j - 1)%squareSize])

    return(nearestNeighbours)

    
