"""
checkAround.py - given a  particle location, if particle lies near edge
of simulation area. If not, checks for neighbouring particles. If no
neighbouring particles, particle randomly walks.

INPUTS: checkAround(location, squareSize, matrix)

location: location of a random walker (list)
squareSize: dimensions of simulation square (int)
matrix: matrix representation of simulation (matrix)


OUTPUTS:

location: location of a random walker (list)
foundFriend: random walker touching cluser (bool)
nearEdge: random walker near edge of square (bool)
orientation: direction of nearest neighbour particle (str)
"""

import random

def checkAround(location, squareSize, matrix):

        # Found another particle
	foundFriend = False 
        # Near the edge of the square
	nearEdge = False
	# Default direction of nearest neighbour is down
	orientation = 'down'
	
        # Check if a walker is near the edge
	if ((location[1] + 1) > squareSize - 5) or ((location[1] - 1) < 1) or \
	((location[0] + 1) > squareSize - 1) or ((location[0] - 1) < 1):
            nearEdge = True

        # If not near the edge, check if the walker is near a neighbor 
	if not nearEdge:
		neighborDown = matrix[location[1] + 1, location[0]]
		if neighborDown == 1:
                        # Found friend below particle
			foundFriend = True
			orientation = 'down'

		neighborUp = matrix[location[1] - 1, location[0]]
		if neighborUp == 1:
                        # Found friend above particle
			foundFriend = True
			orientation = 'up'

		neighborRight = matrix[location[1], location[0]+1]
		if neighborRight == 1:
                        # Found friend to right of particle
			foundFriend = True
			orientation = 'right'

		neighborLeft = matrix[location[1], location[0]-1]
		if neighborLeft == 1:
                        # Found friend to left of particle
			foundFriend = True
			orientation = 'left'

        # After checking location, if no nearest neighbour (Von Neumann): randomly walk one step
	if not foundFriend and not nearEdge:
		decide = random.random()
		if decide < 0.25:
			location = [location[0] - 1, location[1]]
		elif decide < 0.5:
			location = [location[0] + 1, location[1]]
		elif decide < 0.75:
			location = [location[0], location[1] + 1]
		else:
			location = [location[0], location[1] - 1]

	return (location, foundFriend, nearEdge, orientation)
