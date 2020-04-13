"""
runner.py - runs the code from DLACluster.py to simulate the
formation of metal oxide clusters in banded silica within
agates.

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
"""

# Import main DLACluster script
from DLAcluster import DLAcluster 

# Import mass, radius of cluster and matrix representing simulation
mass, clusterRadius, clusterArea, matrix, islands = DLAcluster(200, True, 1000, 50, 0, 1, 1, 1, 9)

