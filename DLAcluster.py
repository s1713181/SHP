"""
DLACluster.py - this is the main function for the DLA cluster model.

INPUTS: DLACluster(squareSize, needGif, blockNumber, layerStep, tempProb, seedNum, alignProb, depMod, clusterMod)

squareSize: dimensions of simulation square (int)
needGif: determines if GIF is produced (boolean)
blockNumber: number of 'blocks' in each sediment layer (int)
layerStep: number of particles added to cluster in between layering (int)
tempProb: probability of surface normal deposition at each layer (float)
seedNum: number of seed particles in the simulation (int)
alignProb: probability of deposition for non-aligned particles (float)
depMod: moderating factor in surface-normal deposition (float)
clusterMod: moderating factor in on-cluster deposition (float)

OUTPUTS:

mass: number of particles added to metal oxide cluster (int)
clusterRadius: distance from seed origin to furthest particle in cluster (float)
matrix: final matrix representation of the simulation (array)

Aligments: 1 == north
           2 == east
           3 == south
           4 == west
"""

import random
import math
import numpy
import matplotlib.pyplot as plt
import os
from matplotlib import colors
from checkAround import checkAround
from randomAtSurface import randomAtSurface
from addLayer import addLayer
import countIslands

# Main simulation script (DLA-CA Process)
def DLAcluster(squareSize, needGif, blockNumber, layerStep, tempProb, seedNum, alignProb, depMod, clusterMod):

    # Check if folder "images" exists, and if not - create it
    if not os.path.isdir("images"):
        os.mkdir("images")
    # Create GIF of formation process
    if needGif:
        # Import imagio if intend to save gif
        import imageio

    if seedNum > 1:
        # Evenly spaced x coordinates between (squareSize/4, 3*squareSize/4)
        seedX = numpy.linspace(squareSize/4, 3*squareSize/4, num=seedNum, dtype=int)
    else:
        # Seed at center of x-axis
        seedX = [int(squareSize / 2)]
    # y coordinate of a seed (Row 1 of matrix)
    seedY = numpy.ones(seedNum, dtype=int)
    # Location of seed particle
    seedLocation = [seedX, seedY]
    # Matrix represents square
    matrix = numpy.zeros((squareSize, squareSize))
    # KPZ Matrix for surface construction
    KPZMatrix = numpy.zeros((squareSize, squareSize))
    # Matrix of crystallographic alignments
    alignMatrix = numpy.zeros((squareSize, squareSize))

    # Set initial matrix conditions
    for row in range (0, squareSize):
        for col in range (0, squareSize):
            if row == 0:
                # Bottom row is solid silica
                matrix[row][col] = 3
    # Place seed particles on 1st row of matrix
    for i in range (len(seedX)):
        y = int(seedY[i])
        x = int(seedX[i])
        matrix[y][x] = 1
        # Assign random crystallographic orientation to seed
        alignMatrix[y][x] = 3

    cmap = colors.ListedColormap(['navy','white','orange','green','black'], N=5)

    # Initialize the random walker counter
    randomWalkersCount = 0

    # Cluster NOT touching top of square
    completeCluster = False
    # Surface NOT touching top of square
    completeSurface = False
    # Boolean for particle added in last iteration (for layering)
    newAddedCount = True

    # Number of random walkers added to the cluster
    addedCount = 0 
    # Initialize variable for number of added layers
    existingLayers = 0

    # Initialize array for the used interval for graphing
    usedInterval = []

    # Add initial KPZ Layer to the simulation (Cavity edge)
    KPZMatrix, existingLayers = addLayer(KPZMatrix, matrix, blockNumber, squareSize, existingLayers, tempProb, depMod, clusterMod)
    for i in range (1, squareSize - 2):
        for j in range (1, squareSize - 1):
            if matrix[i][j] != 1 and i != 0:
                matrix[i][j] = KPZMatrix[i][j]
    islands = countIslands.countIslands(matrix.tolist(), squareSize)
    print("Layer added, number of islands = ", str(islands))

    # Simulation stops when cluster or surfaces touch top of the square
    while not completeCluster and not completeSurface:

        # Release a walker
        randomWalkersCount += 1
        random.seed()

        # Generate an initial position for a walker on the surface of the square
        location = randomAtSurface(squareSize)

        # Initialize variables for finding friend / leaving square / cluster reaching surface / layers reaching surface / walker in solid / cluster enclosed
        foundFriend = False # Not near other particle
        nearEdge = False # Not near the edge of the square
        topSquare = False # Cluster not at edge of the square
        surfaceAtEdge = False # Surface not at edge of square
        inSolid = False # Particle not in solid silica
        clusterNotEnclosed = True # Cluster is not completely enclosed in silica

        # Set an individual walker out, stop if found a neighbouring particle, give up if it reached the edge of the simulation area
        while (not foundFriend) and (not nearEdge) and (not inSolid) and (not surfaceAtEdge) and (clusterNotEnclosed):
            # Run the walking function
            locationNew, foundFriend, nearEdge, orientation = checkAround(location, squareSize, matrix)
                                    
            # Add layer to simulation if 'layerStep' particles newly added to simulation 
            if addedCount%layerStep == 0 and addedCount != 0 and newAddedCount == True:
                KPZMatrix, existingLayers = addLayer(KPZMatrix, matrix, blockNumber, squareSize, existingLayers, tempProb, depMod, clusterMod)
                for i in range (1, squareSize - 2):
                    for j in range (1, squareSize - 1):
                        if matrix[i][j] != 1 and i != 0:
                            matrix[i][j] = KPZMatrix[i][j]
                # Count matrix 'islands' for measure of anastomosis
                islands = countIslands.countIslands(matrix.tolist(), squareSize)
                print("Layer added, number of islands = ", str(islands))
                
                # Check if surface near top of square
                for i in range(squareSize - 15, squareSize - 1):
                    for j in range (0, squareSize - 1):
                        if matrix[i][j] == (2 or 3 or 4):
                            surfaceAtEdge = True

                # Check if cluster is enclosed
                indexList = []
                for i in range(0, squareSize - 1):
                    for j in range(squareSize - 1, 0, -1):
                        if matrix[j][i] != 0 and matrix[j][i] != 4:
                            indexList.append(matrix[j][i])
                            break

                # If cluster not enclosed by layer - continue
                if 1 not in indexList:
                    clusterNotEnclosed = False
                              
        
            # Add to the cluster if neighbouring a particle in the cluster
            if foundFriend:
                # Randomly choose crystalographic aligment of walker
                cell = random.choice([1, 2, 3, 4])
                # Number to determine deposition of non-aligned particles
                marker = random.uniform(0,1)
                # Current location not near top of square, replace with 1 and stop
                if matrix[location[1]][location[0]] == 4 and location[1] < (squareSize - 5) :
                    if orientation == 'down':
                        # Walker aligned with nearest neighbour
                        if alignMatrix[location[1] + 1][location[0]] == cell and cell == 1:
                            matrix[location[1]][location[0]] = 1
                            # Update matrix of alignments
                            alignMatrix[location[1]][location[0]] = cell
                            addedCount += 1
                            newAddedCount = True
                        # If not aligned - deposition is stochastic
                        else:
                            if marker < alignProb:
                                matrix[location[1]][location[0]] = 1
                                # Update matrix of alignments
                                alignMatrix[location[1]][location[0]] = cell
                                addedCount += 1
                                newAddedCount = True
                    elif orientation == 'up':
                        # Walker aligned with nearest neighbour
                        if alignMatrix[location[1] - 1][location[0]] == cell and cell == 3:
                            matrix[location[1]][location[0]] = 1
                            # Update matrix of alignments
                            alignMatrix[location[1]][location[0]] = cell
                            addedCount += 1
                            newAddedCount = True
                        # If not aligned - deposition is stochastic
                        else:
                            if marker < alignProb:
                                matrix[location[1]][location[0]] = 1
                                alignMatrix[location[1]][location[0]] = cell
                                addedCount += 1
                                newAddedCount = True
                    elif orientation == 'left':
                        # Walker aligned with nearest neighbour
                        if alignMatrix[location[1]][location[0] - 1] == cell and cell == 2:
                            matrix[location[1]][location[0]] = 1
                            # Update matrix of alignments
                            alignMatrix[location[1]][location[0]] = cell
                            addedCount += 1
                            newAddedCount = True
                        # If not aligned - deposition is stochastic
                        else:
                            if marker < alignProb:
                                matrix[location[1]][location[0]] = 1
                                # Update matrix of alignments
                                alignMatrix[location[1]][location[0]] = cell
                                addedCount += 1
                                newAddedCount = True
                    elif orientation == 'right':
                        # Walker aligned with nearest neighbour
                        if alignMatrix[location[1]][location[0] + 1] == cell and cell == 4:
                            matrix[location[1]][location[0]] = 1
                            # Update matrix of alignments
                            alignMatrix[location[1]][location[0]] = cell
                            addedCount += 1
                            newAddedCount = True
                        # If not aligned - deposition is stochastic
                        else:
                            if marker < alignProb:
                                matrix[location[1]][location[0]] = 1
                                # Update matrix of alignments
                                alignMatrix[location[1]][location[0]] = cell
                                addedCount += 1
                                newAddedCount = True

                    # Finish conditions for non-constant radius
                    """
                    # Cluster touches top of square   
                    for i in range(squareSize - 10, squareSize):
                        for j in range (0, squareSize):
                            if matrix[i][j] == 1:
                                completeCluster = True

                    # Cluster touches RHS of square
                    for i in range(0, squareSize):
                        for j in range (squareSize - 5, squareSize):
                            if matrix[i][j] == 1:
                                completeCluster = True

                    # Cluster touches LHS of square
                    for i in range(0, squareSize):
                        for j in range (0, 5):
                            if matrix[i][j] == 1:
                                completeCluster = True
                    """
                    
                    for i in range(squareSize):
                        for j in range(squareSize):
                            if matrix[i][j] == 1 and (i**2 + (j-squareSize/2)**2)**(1/2) > (squareSize/2 - 5):
                                completeCluster = True
                    
                    
            # Otherwise, save the location
            else:
                # Forbid walking through solid silica
                if matrix[locationNew[1]][locationNew[0]] == (2 or 3 or 6):
                    inSolid = True
                else:
                    location = locationNew
                    newAddedCount = False
        
        # Print update 
        intervalSavePic=range(2,4000000000,500)
        if randomWalkersCount in intervalSavePic:
            """
            print("Walkers added: ", randomWalkersCount, " Walkers added to cluster: ", addedCount)
            """
        if needGif:
            if randomWalkersCount in intervalSavePic:
                """
                print("Saved picture")
                """
                usedInterval.append(randomWalkersCount) #append to the used count
                label=str(randomWalkersCount)
                plt.matshow(matrix, interpolation='nearest',cmap=cmap)#plt.cm.Blues) #ocean, Paired
                plt.axis('off')
                plt.savefig("images/cluster{}.png".format(label), dpi=200)
                plt.close()
       
        # Prevent infinite simulation loop
        if randomWalkersCount == 20000000:
            print("CAUTION: had to break the cycle, taking too many iterations")
            completeCluster = True

        # Stop when cluster reaches edge of square
        if completeCluster == True:
            print("Walkers added to cluster: ", addedCount)
            print("Cluster reached edge of square")

        # Stop when surface reaches top of square or cluster is enclosed by layers
        if (surfaceAtEdge == True) or (clusterNotEnclosed == False):
            print("Walkers added to cluster: ", addedCount)
            print("Finished, surface reached edge or encloses the cluster entirely ")
            completeSurface = True

    # Constant radius
    clusterRadius = (squareSize/2 - 5)

    # Not constant radius
    """
    
    clusterRadii = []
    for i in range(0, squareSize - 1):
        for j in range (0, squareSize - 1):
            # Update until found furthest distance from seed particle
            if matrix[i][j] == 1:
                clusterRadii.append((i**2 + (j - squareSize/2)**2)**(1/2))
    
    clusterRadius = max(clusterRadii)
    """

    # Calculate cluster area (max(x)-min(x))*(max(y)-min(y))
    minxList = []
    for i in range(squareSize):
        for j in range(squareSize):
            if matrix[i][j] == 1:
                minxList.append(j)
                break
    minx = min(minxList)

    minyList = []
    for i in range(squareSize):
        for j in range(squareSize):
            if matrix[j][i] == 1:
                minyList.append(j)
                break
    miny = min(minyList)

    maxxList = []
    for i in range(squareSize):
        for j in range(squareSize - 1, -1, -1):
            if matrix[i][j] == 1:
                maxxList.append(j)
                break
    maxx = max(maxxList)

    maxyList = []
    for i in range(squareSize):
        for j in range(squareSize - 1, -1, -1):
            if matrix[j][i] == 1:
                maxyList.append(j)
                break
    maxy = max(maxyList)

    clusterArea = (maxx-minx)*(maxy-miny)
    
        
    # Generate final image of cluster and GIF of simulation
    plt.matshow(matrix, interpolation='nearest',cmap=cmap)
    plt.axis('off')
    plt.savefig("images/cluster.png", dpi=200)
    plt.close()

    if needGif:
        with imageio.get_writer('images/movie.gif', mode='I') as writer:
            for i in usedInterval:
                filename = "images/cluster" + str(i) + ".png"
                image = imageio.imread(filename)
                writer.append_data(image)
                os.remove(filename)
            image = imageio.imread("images/cluster.png")
            writer.append_data(image)

    # Return walkers in cluster / cluster radius / final simulation matrix
    return (addedCount, clusterRadius, clusterArea, matrix, islands)

