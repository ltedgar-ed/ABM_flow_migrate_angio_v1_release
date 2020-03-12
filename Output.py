#------------------------------------------------------------------------------
# The Output class handles all printing, output, writing to files, etc
# *** NOTE: The Output class requires a folder inside the working directory,
# 'Simulation Data,' which it uses to store the output files
# Lowell Taylor Edgar
# University of Edinburgh
# 2019

import scipy.io as sio
import numpy as np
import matlab.engine

from Input import *

import os
#------------------------------------------------------------------------------
# The out_file class stores all the information for writing the output file
class OutFile:
    # OutFile Constructor
    def __init__(self, filename, nodes, vessels, node_degree, BCs, bifur_nodes, bifur):
        # Name of the ouput file
        self.filename = filename
        
        # Container of nodal positions
        self.nodes = nodes.copy()
        
        # Number of vessels
        v = len(vessels)
        
        # Container of vessel connectivity
        self.vess_conn = np.zeros((v, 2))

        for vess in vessels:
            self.vess_conn[vess.ID-1][0] = vess.n0
            self.vess_conn[vess.ID-1][1] = vess.n1
        
        # Container of nodal pressures
        self.nodal_pressures = np.zeros((len(nodes), Nt+1))
        
        # Container of vessel diameter
        self.vess_diameter = np.zeros((v, Nt+1))
        
        # Container of vessel flow
        self.vess_flow = np.zeros((v, Nt+1))
        
        # Container of vessel wall shear stress
        self.vess_WSS = np.zeros((v, Nt+1))
        
        # Container of vessel cell number
        self.vess_num_cells = np.zeros((v, Nt+1))
        
        # Container for the boundary conditions
        self.BCs = BCs
        
        # Containter for the bifurcation nodes
        self.bifur_nodes = bifur_nodes
        
        # Container for bifurcation status
        self.bifur_status = np.zeros((len(bifur), Nt+1))
        
        # Initialize the containers
        for vess in vessels:
            self.nodal_pressures[vess.n0][0] = vess.P0
            self.nodal_pressures[vess.n1][0] = vess.P1
                
            self.vess_diameter[vess.ID-1][0] = vess.D
            self.vess_flow[vess.ID-1][0] = vess.Q
            self.vess_WSS[vess.ID-1][0] = vess.tau
            self.vess_num_cells[vess.ID-1][0] = vess.num_cells
                
        for i in range(len(bifur)):
            self.bifur_status[i][0] = bifur[i]

    
    # Save the vessel state to the containers
    def save_vessels(self, vessels, bifur, time_step):
        for vess in vessels:
            self.nodal_pressures[vess.n0][time_step] = vess.P0
            self.nodal_pressures[vess.n1][time_step] = vess.P1
                
            self.vess_diameter[vess.ID-1][time_step] = vess.D
            self.vess_flow[vess.ID-1][time_step] = vess.Q
            self.vess_WSS[vess.ID-1][time_step] = vess.tau
            self.vess_num_cells[vess.ID-1][time_step] = vess.num_cells
            

        for i in range(len(bifur)):
            self.bifur_status[i][time_step] = bifur[i]
            
        
    # Write the final output file once the simulation has ended    
    def write_output_file(self):        
        # Save the output of the simulation as a MATLAB .mat file in the Simulation Data directory
        curr_path = os.getcwd()
        save_path = curr_path + "\Simulation Data"
        os.chdir(save_path)        
        sio.savemat(self.filename, {'nodes' : self.nodes, 'vess_conn' : self.vess_conn, 'nodal_pressures' : self.nodal_pressures, 'vess_diameter' : self.vess_diameter, 'vess_flow' : self.vess_flow, 'vess_WSS' : self.vess_WSS, 'vess_num_cells' : self.vess_num_cells, 'PBCs' : self.BCs, 'bifur_nodes' : self.bifur_nodes, 'bifur' : self.bifur_status})
        os.chdir(curr_path)
        
        # Open the MATLAB plotting script to visualize the results
        if (plot_animation == True):    
            print("Plotting animation...")
            
            mateng = matlab.engine.start_matlab()
            
            if (sim_type == "A branch"):
                mateng.plot_ABM_output_A_branch(self.filename, nargout=0)
            elif (sim_type == "Y branch"):
                mateng.plot_ABM_output_Y_branch(self.filename, nargout=0)
            elif (sim_type == "ideal cap bed"):
                mateng.plot_ABM_output_ideal_cap(self.filename, nargout=0)
            elif (sim_type == "ideal cap bed2"):
                mateng.plot_ABM_output_ideal_cap2(self.filename, nargout=0)
            elif (sim_type == "ideal sprout front"):
                mateng.plot_ABM_output_ideal_sprout_front(self.filename, nargout=0)
            
            mateng.quit()
        
        
#------------------------------------------------------------------------------
# Print the nodal positions to the console
def print_nodes(nodes):
    # Print the nodal array     
    for node in nodes:
        print("({0}, {1})".format(round(node[0] *1e6, 1), round(node[1] *1e6, 1)))


#------------------------------------------------------------------------------
# Print the vessels status to the console
def print_vessels(vessels):
    # Print the vessels to the console
    print()
    
    for vess in vessels:
        print(vess)
        

#------------------------------------------------------------------------------
# Print info on vessel neighbours to the console
def print_vess_neighbours(vessels):
    # Print info on neighbours to console for verification            
    for vess in vessels:
        print("Vessel {0} has upstream neighbors {1} and downstream neighbors {2}".format(vess.ID, vess.neigh0, vess.neigh1))


#------------------------------------------------------------------------------  
# Print information on the cells to the console
def print_cells(vessels):
    # Print the cells to the console
    print()
    for vess in vessels:
        for cell in vess.cells:
            print(cell)

#------------------------------------------------------------------------------
# Print the number of cells in each vessel to the consoles
def print_num_cells(vessels):
    # Print the number of cells in each vessel to the console
    print()
    for vess in vessels:
        print("Vessel {0} has {1} cells with cell number set as {2}".format(vess.ID, len(vess.cells), vess.num_cells))

    
    
    
    
    
    