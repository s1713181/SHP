"""
countIslands.py - count the number of silica islands in simulation area for
measure of anastomosis, contains function to implement bfs algorithm
iteratively. 

INPUTS: countIslands(matrix, squareSize)

matrix: matrix representing simulation area (matrix)
squareSize: dimensions of simulation square (int)


OUTPUTS:

islands: number of distinct silica islands in the simulation (int)

"""

import math
import numpy
import random


# Count islands in matrix (to quantify anastomosis)
def countIslands(matrix, squareSize):
    # Make matrix binary (0 = silica [any state] , 1 = metal-oxide)
    for i in range(squareSize):
        for j in range(squareSize):
            if matrix[i][j] != 1:
                matrix[i][j] = 0
    if not matrix:
        return 0

    # Dictionary of visited sites
    visited = {}

    # Initiate variable for islands
    islands = 0

    row = squareSize
    col = squareSize
    for i in range (squareSize):
        for j in range (squareSize):
            # Unvisited silica site - add one to island variable
            # Fill Island with bfs algorithm
            if matrix[i][j] == (0) and (i,j) not in visited:
                bfs(matrix, i, j, row, col, visited)
                islands +=1

    # Number of silica 'islands' in simulation 
    return islands

# Iterative breadth first search algorithm for use in countIslands
def bfs(matrix, i, j, row, col, visited):

    # Sites to visit 
    queue = [(i,j)]

    # Conduct bfs algorithm
    while queue:
        # Sites to visit next
        next = []
        for u in queue:
            # Site must not have already been visited
            if u not in visited:
                visited[u] = 1
                # Adjacent unvisited silica sites to visit next
                if u[1] + 1 < col and matrix[u[0]][u[1]+1] == (0):
                    next.append((u[0], u[1]+1))
                if u[1] - 1 > 0 and matrix[u[0]][u[1]-1] == (0):
                    next.append((u[0], u[1]-1))
                if u[0] + 1 < col and matrix[u[0] + 1][u[1]] == (0):
                    next.append((u[0] + 1, u[1]))
                if u[0] + 1 > 0 and matrix[u[0] - 1][u[1]] == (0):
                    next.append((u[0] - 1, u[1]))
        # Update sites to visit
        queue = next
        
        
        
        
        


