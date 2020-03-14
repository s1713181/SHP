"""
randomAtSurface.py - generates a random position at the top
of the simulation area for a metal oxide particle to begin a
random walk in DLACluster.py.

INPUTS: randomAtSurface(squareSize)

squareSize: dimensions of simulation square (int)

OUTPUTS:

location: starting position of random walker (list)

"""

import random

def randomAtSurface(squareSize):
    
    marker = random.uniform(0, 1)
    # Add 5 to parameters to avoid edges of square
    
    if marker < float(4/6):
    
        # Release at top of simulation area
        x = random.randint(5, squareSize - 6) # x coordinate
        y = (squareSize - 5) # y coordinate 
        location = [x, y] 
    
    elif marker > float(5/6):
        # Release at left of simulation area
        x = 5
        y = random.randint(5, squareSize - 6)       
        location = [x, y] 

    else:
        # Release at right of simulation area
        x = squareSize - 6
        y = random.randint(5, squareSize - 6)       
        location = [x, y]
        
    # Position of released particle 
    return (location)

    
