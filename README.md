Code for simulating Agate formation in Senior Honours Project                        
                                                                                                                                                                       
Main simulation of agate genesis occurs in DLACluster.py module

Main simulation (DLACluster.py) uses: CheckAround.py
                                      randomAtSurface.py
                                      addLayer.py
                                      countIslands.py

Main simulation is executed in runner.py module:

INPUTS:

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
clusterArea: area of cluster (int)
matrix: final matrix representation of the simulation (array)

Simulation produces a snapshot of the BD-DLA process after 5000 simulation iterations. 
                                   
######################################################################################

universalityClass.py is used to verify the universality class of surfaces produced in 
ballistic deposition. 

INPUTS: 

blockNumber: number of 'blocks' in each sediment layer (int)
squareSize: dimensions of lattice in which ballistic deposition occurs
trialNumber: number of trials to average results over






