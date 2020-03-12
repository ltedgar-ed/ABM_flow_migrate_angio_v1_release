#------------------------------------------------------------------------------
# The Input class contains all the input parameters, can be imported into any module
# that requires these parameters using import *
#
# Lowell Taylor Edgar
# University of Edinburgh
# 2019

import sys
from random_seed_numbers import *

#------------------------------------------------------------------------------
# Input parameters

# Uncomment if running from bash script in the command line
branch_alpha = float(sys.argv[1])
run = int(sys.argv[2])

# Uncomment if running as stand alone in Spyder
#run = 500
#branch_alpha = 0.45

# Random seed generated from rand.randint(1,1000000000)
rseed = rseeds[run-1]

# Type of simulation to run ("A branch", "ideal cap bed")
sim_type = "A branch"
A_branch_type = 0
A_branch_dirchlet = False
#sim_type = "Y branch"
#Y_branch_type = 0
#Y_branch_dirchlet = False
#sim_type = "ideal cap bed"
#sim_type = "ideal cap bed2"
#sim_type = "ideal sprout front"


# Number of time steps
Nt = 72                                                                       

# Dynamic viscosity of blood (Pa-s)
mu = 3.5 * 1e-3                              


# Width of each cell (m)                                             
cell_size = 5 * 1e-6                                                           


# Turn vessel smoothing on or off
yes_smoothing = True


# Create the animation of the simulation
plot_animation = False                                                         # Turn on the MATLAB plotting at the end of the simuation


# Inlet flow condition
inlet_flow_on = False


# Polarization re - alignment weights: w1 + w2 + w3 = 1
w2 = 1.0                                                                        # Polarization weight - flow component
w3 = 0.0                                                                        # Polarization weight - random walk component
w1 = 1 - w2 - w3                                                                # Polarization weight - persistence component


# Turn on shear weighting for polarization realignment
yes_shear_polar = False

if (yes_shear_polar == True):
    w1 = 0.5
    tau_max = 1.0


# -------------------------
# Bifurcation rule parameters
branch_rule = 7                                                                
#branch_alpha = 1.0
yes_new_branch = False

bp_u = 1;
bp_v = 0;


# -------------------------
# Ideal sprouting front parameters

if (sim_type == "ideal sprout front"):
    # Length of each vessel segment/length of each cell (m)
    vess_length = 6.25 * 1e-6   
    
    # Pressure boundary conditions [Partin, Partout, Pveinout, Pveinin] (Pa)
    PBC = [7546, 0]
    
    # Number of honeycombs long and high
    network_num_hc_long = 5
    network_num_hc_high = 5
    
    # Number of cells in the artery, vein, and capillaries
    num_cell_art = 10
    num_cell_vein = 20
    num_cell_cap = 5
    
    # Name of the MATLAB file containing the network geometry
    in_cap_bed_filename = "ideal_sprout_front.mat"
    
    # Apply cell number Dirichlet conditions at the cell inlets of artery and vein
    yes_dirichlet_BCs = False
    
    # Apply cell number Dirichlet conditions along the whole artery and vein
    yes_art_vein_dirichlet = True
    
    # Apply special bifurcation options along the vein (1 or 2)
    vein_bf_option = 1
    
    
# -------------------------
# Ideal capillary bed parameters

if (sim_type == "ideal cap bed"):
    # Length of each vessel segment/length of each cell (m)
    vess_length = 5 * 1e-6   
    
    # Width of each cell (m)                                             
    cell_size = 5 * 1e-6   
    
    # Pressure boundary conditions [Partin, Partout, Pveinout, Pveinin] (Pa)
    PBC = [55*133.322, 50*133.322, 0*133.322, 5*133.322]
    
    # Number of honeycombs long and high
    network_num_hc_long = 7
    network_num_hc_high = 7
    
    # Number of cells in the artery, vein, and capillaries
    num_cell_art = 10
    num_cell_vein = 20
    num_cell_cap = 5
    
    # Name of the MATLAB file containing the network geometry
    in_cap_bed_filename = "ideal_cap_bed.mat"
    
    # Apply cell number Dirichlet conditions at the cell inlets of artery and vein
    yes_dirichlet_BCs = False
    
    # Apply cell number Dirichlet conditions along the whole artery and vein
    yes_art_vein_dirichlet = True
    
    # Apply special bifurcation options along the vein (1 or 2)
    vein_bf_option = 1


# -------------------------
# Ideal capillary bed parameters

if (sim_type == "ideal cap bed2"):
    # Length of each vessel segment/length of each cell (m)
    vess_length = 5 * 1e-6   
    
    # Width of each cell (m)                                             
    cell_size = 5 * 1e-6   
    
    # Pressure boundary conditions [Part, Pvein] (Pa)
    Part = 55*133.322
    Pvein = 0*133.322
    
    # Number of honeycombs long and high
    network_num_hc_long = 7
    network_num_hc_high = 7
    
    # Number of cells in the and capillaries
    num_cell_cap = 5
    
    # Name of the MATLAB file containing the network geometry
    in_cap_bed_filename = "ideal_cap_bed2.mat"
    
    # Apply cell number Dirichlet conditions at the cell inlets of artery and vein
    yes_dirichlet_BCs = False
    
    # Apply cell number Dirichlet conditions along the whole artery and vein
    yes_art_vein_dirichlet = True
    
    # Apply special bifurcation options along the vein (1 or 2)
    vein_bf_option = 1
    
    
# -------------------------
# A branch model parameters

if (sim_type == "A branch"):    
    # Length of each vessel segment/length of each cell (m)
    vess_length = 2*cell_size
    
    # Inlet pressure (Pa)
    Pin = 200       
    
    # Outlet pressure (Pa)
    Pout = 0                             
    
    # Initial number of cells in vessel
    num_cell = 8                                                                    
    
    if (A_branch_type == 0):
        # Number of nodes
        Nn = 40          
        # Number of vessel segments                                                               
        Nseg = 40    
    elif (A_branch_type == 1):
        Nn = 32
        Nseg = 32
    elif (A_branch_type == 2):
        Nn = 24
        Nseg = 24                                                                   

# -------------------------
# Y branch model parameters

if (sim_type == "Y branch"):    
    # Length of each vessel segment/length of each cell (m)
    vess_length = 2*cell_size
    
    # Inlet pressure left (Pa)
    Pleft = 100       
    
    # Inlet pressure right (Pa)
    Pright = 100
    
    # Outlet pressure (Pa)
    Pout = 0                             
    
    # Initial number of cells in vessel
    num_cell = 8                                                                    
    
    if (Y_branch_type == 0):
        # Number of nodes
        Nn = 31     
        # Number of vessel segments                                                               
        Nseg = 30    
    elif (Y_branch_type == 1):
        # Number of nodes
        Nn = 41     
        # Number of vessel segments                                                               
        Nseg = 40    
    elif (Y_branch_type == 2):
        # Number of nodes
        Nn = 16     
        # Number of vessel segments                                                               
        Nseg = 15  
        

# -------------------------
# Output file name

# Output file for A branch simulation
if (sim_type == "A branch"):
    out_filename = "ABM_output_A_branch" + str(A_branch_type) + ("_Nt_" + str(Nt)) + ("_Ncell_" + str(num_cell)) +  ("_bifrule_" + str(branch_rule))
    
    if (branch_rule == 7):
        out_filename += "_alpha_" + "{0:.3f}".format(branch_alpha)
        

# Output file for Y branch simulation
if (sim_type == "Y branch"):
    out_filename = "ABM_output_Y_branch" + ("_Nt_" + str(Nt)) + ("_Ncell_" + str(num_cell)) +  ("_bifrule_" + str(branch_rule))
    
    if (branch_rule == 7):
        out_filename += "_alpha_" + "{0:.3f}".format(branch_alpha)


# Output file for ideal capillary bed simulation
if ((sim_type == "ideal cap bed") or (sim_type == "ideal cap bed2")):
    out_filename = "ABM_output_ideal_cap_bed" + ("_" + str(network_num_hc_long) + "by" + str(network_num_hc_high)) + ("_Ncell_cap_" + str(num_cell_cap)) + ("_Nt_" + str(Nt)) +  ("_bifrule_" + str(branch_rule))
    
    if (yes_dirichlet_BCs == True and yes_art_vein_dirichlet != True):
        out_filename += "_inlet_dirichlet"
    
    if (branch_rule == 7):
        out_filename += "_alpha_" + "{0:.3f}".format(branch_alpha)
               
    out_filename += "_vein_bf" + str(vein_bf_option)
    

# Output file for ideal capillary bed simulation
if (sim_type == "ideal sprout front"):
    out_filename = "ABM_output_ideal_sprout_front" + ("_" + str(network_num_hc_long) + "by" + str(network_num_hc_high)) + ("_Ncell_cap_" + str(num_cell_cap)) + ("_Nt_" + str(Nt)) +  ("_bifrule_" + str(branch_rule))
    
    if (yes_dirichlet_BCs == True and yes_art_vein_dirichlet != True):
        out_filename += "_inlet_dirichlet"
    
    if (branch_rule == 7):
        out_filename += "_alpha_" + str(branch_alpha)
        
    out_filename += "_vein_bf" + str(vein_bf_option)

   
# Modify output file name is old branching rule is used
#if (yes_new_branch == False):
#    out_filename += "_old_branch_rule"

    
# Add run number to the output file name
out_filename += ("_run" + str(run)) 