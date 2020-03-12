#------------------------------------------------------------------------------
# An agent-based model of the flow-coupled migration of vascular endothelial cells
# during developmental angiogenesis
#
# Lowell Taylor Edgar
# University of Edinburgh
# 2019

import random as rand

from Input import *
from Vessel import *

import Geom
import Network
import Flow
import Output

#------------------------------------------------------------------------------


# Set the random seed
#rand.seed(123456789)
rand.seed(rseed)

# Initialize storage lists
nodes = []                                                                      # List of nodes (initially empty)
vessels = []                                                                    # List of vessel segments (intially empty)
node_degree = []                                                                # List of nodal degree (initially empty)
BCs = []                                                                        # List of boundary conditions (intially empty)
bifur_nodes = []                                                                # List of bifurcation nodes (intially empty)
bifur = []                                                                      # List of bifurcation state (initially empty)

# Create the network
Geom.create_network(nodes, vessels, node_degree, bifur_nodes, bifur, BCs)

# Solve for flow
Flow.solve_for_flow(nodes, vessels, node_degree, BCs, bifur_nodes, bifur)

#Output.print_vessels(vessels)

# Initial polarity realignment
Network.realign_polarity(vessels)

#Output.print_cells(vessels)

# Print out the time step
print("Time step 0...".format(0, Nt))

# Initialize the ouput file
out = Output.OutFile(out_filename, nodes, vessels, node_degree, BCs, bifur_nodes, bifur)

# Step through time
for t in range(Nt):
    # Print out the time step
    print("Time step {0} out of {1}...".format(t+1, Nt))
    
    # If not the first time step, realign polarity
    if (t != 0):
        Network.realign_polarity(vessels)

    # Determine cell migration
    Network.cell_migration(vessels, bifur_nodes, bifur)
        
    # Solve for flow in the new network configuration
    Flow.solve_for_flow(nodes, vessels, node_degree, BCs, bifur_nodes, bifur)
    
    # Save the vessel state for output
    out.save_vessels(vessels, bifur, t+1)


# Write the output files  
out.write_output_file()    

# Terminate the simulation
print("Done!")