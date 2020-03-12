#------------------------------------------------------------------------------
# The Geom class contains the methods needed to construct the vascular network
# as a collection of nodes, vessel segments, and boundary conditions
#
# Lowell Taylor Edgar
# University of Edinburgh
# 2019

import scipy.io as sio

from math import cos
from math import sin
from math import radians

from Vessel import *
from Tensor3D import *
from Input import *
from ECell import *

import matlab.engine

#------------------------------------------------------------------------------
# Create the vessel network
def create_network(nodes, vessels, node_degree, bifur_nodes, bifur, BCs):

    if (sim_type == "A branch"):
        if (A_branch_type == 0):
            create_old_A_branch_network(nodes, vessels, node_degree, bifur_nodes, bifur, BCs)
        if (A_branch_type == 1):
            create_small_A_branch_network1(nodes, vessels, node_degree, bifur_nodes, bifur, BCs)
        if (A_branch_type == 2):
            create_small_A_branch_network2(nodes, vessels, node_degree, bifur_nodes, bifur, BCs)
        if (A_branch_type == 3):
            create_longer_distal_A_branch_network(nodes, vessels, node_degree, bifur_nodes, bifur, BCs)
           
    elif (sim_type == "Y branch"):
        if (Y_branch_type == 0):
            create_Y_branch_network(nodes, vessels, node_degree, bifur_nodes, bifur, BCs)
        if (Y_branch_type == 1):
            create_Y_branch_network2(nodes, vessels, node_degree, bifur_nodes, bifur, BCs)
        if (Y_branch_type == 2):
            create_small_Y_branch_network(nodes, vessels, node_degree, bifur_nodes, bifur, BCs)
            
    elif (sim_type == "ideal cap bed"):
        load_ideal_cap_bed(nodes, vessels, node_degree, bifur_nodes, bifur, BCs)
    
    elif (sim_type == "ideal cap bed2"):
        load_ideal_cap_bed2(nodes, vessels, node_degree, bifur_nodes, bifur, BCs)
        
    elif (sim_type == "ideal sprout front"):
        load_ideal_sprout_front(nodes, vessels, node_degree, bifur_nodes, bifur, BCs)
    
  
#------------------------------------------------------------------------------
# Create the vessel network (load the ideal capillary bed from the MATLAB file)
def load_ideal_cap_bed(nodes, vessels, node_degree, bifur_nodes, bifur, BCs):
    # Generate the ideal capillary bed using the MATLAB script
    mateng = matlab.engine.start_matlab()
    mateng.generate_ideal_cap_bed(network_num_hc_long, network_num_hc_high, num_cell_art, num_cell_vein, num_cell_cap, nargout=0)
    mateng.quit()
    
    # Load the network geometry from the MATLAB file
    ideal_cap_bed = sio.loadmat(in_cap_bed_filename)
    vess_seg = ideal_cap_bed['vess_seg']
    
    filtered_vess_seg = []                                                      # Container of filtered segments
    rp = 8                                                                      # Rounding precision
    
    # Filter out any repeated vessel segments from the raw data
    for seg in vess_seg:
         
        yes_append = True
        
        # Optimize this if possible, must be a better way to do this
        for fvess_seg in filtered_vess_seg:
            
            if (round(seg[0], rp) == round(fvess_seg[0], rp) and round(seg[1], rp) == round(fvess_seg[1], rp) and round(seg[2], rp) == round(fvess_seg[2], rp) and round(seg[3], rp) == round(fvess_seg[3], rp)):
                yes_append = False
                
        if (yes_append == True):
            filtered_vess_seg.append(seg)   
    
    # Create the sub-vessel segments for the ABM based on the input segment size    
    for seg in filtered_vess_seg:
        x0 = float(seg[0] *1e-6)
        y0 = float(seg[1] *1e-6)
        x1 = float(seg[2] *1e-6)
        y1 = float(seg[3] *1e-6)
        
        seg_vect = Vect(x1 - x0, y1 - y0, 0.)
        seg_length = seg_vect.length()
        seg_vect.unit()
        
        num_sub_vess = int(round(seg_length/vess_length))
        seg_vess_length = seg_length/num_sub_vess
        
        for i in range(num_sub_vess):
            node1 = (round(x0 + i*(seg_vess_length*seg_vect.x), rp), round(y0 + i*(seg_vess_length*seg_vect.y), rp))
            node2 = (round(x0 + (i+1)*(seg_vess_length*seg_vect.x), rp), round(y0 + (i+1)*(seg_vess_length*seg_vect.y), rp))
            
            if (node1 not in nodes):
                nodes.append(node1)
                
            if (node2 not in nodes):
                nodes.append(node2)
                
            index1 = nodes.index(node1)
            index2 = nodes.index(node2)
            
            vess = Vessel(len(vessels)+1, index1, index2, int(seg[4]), mu)
            
            if (vess.num_cells == num_cell_art):
                vess.type = "artery"
            elif (vess.num_cells == num_cell_vein):
                vess.type = "vein"
                
            vess.calc_length(node1, node2)
            vess.update_diameter()
            vess.update_conductance()
            
            vessels.append(vess)
    
    # Update the cell number identifiers within all vessel segments
    cell_count = 1
    
    for vess in vessels:        
        for cell in vess.cells:
            cell.ID = cell_count
            cell_count += 1
            
    # Give vessels information about their neighbours
    BC_nodes = []
    
    find_vessel_neighbours(nodes, vessels, BC_nodes)
    
    # Calculate the degree of each node
    for node in nodes:
        node_degree.append(0)
    
    for vess in vessels:
        node_degree[vess.n0] += 1
        node_degree[vess.n1] += 1
    
    # Store the nodes that are bifurcation points
    for i in range(len(node_degree)):
        if (node_degree[i] == 3):
            bifur_nodes.append(i)
    
    # Initialize the bifurcation status array
    for i in range(len(bifur_nodes)):
        bifur.append(0)
        
    # Output network info
#    Output.print_nodes(nodes)
#    Output.print_vessels(vessels)
    
    # Prescribe boundary conditions
    assert (len(BC_nodes) == len(PBC)), "Number of boundary nodes and number of prescribed pressures is not the same"
    
    for i in range(len(BC_nodes)):
        BCs.append((BC_nodes[i], PBC[i]))
        
    # Apply Dirichlet cell conditions at the boundary nodes
    x_min, y_min = min(nodes, key=lambda item:item)
    x_max, y_max = max(nodes, key=lambda item:item)
    
    node_top_left = nodes.index((x_min, y_max))
    node_bottom_right = nodes.index((x_max, y_min))
    
    if (yes_dirichlet_BCs):
        for vess in vessels:
            if (vess.n0 == node_top_left or vess.n1 == node_top_left):
                vess.dirichlet = vess.num_cells
            if (vess.n0 == node_bottom_right or vess.n1 == node_bottom_right):
                vess.dirichlet = vess.num_cells
    
    # Apply Dirichlet cell conditions along the whole artery and bein
    if (yes_art_vein_dirichlet): 
        for vess in vessels:
            if (vess.type == "artery" or vess.type == "vein"):
                vess.dirichlet = vess.num_cells


#------------------------------------------------------------------------------
# Create the vessel network (load the ideal capillary bed from the MATLAB file)
def load_ideal_cap_bed2(nodes, vessels, node_degree, bifur_nodes, bifur, BCs):
    # Generate the ideal capillary bed using the MATLAB script
#    mateng = matlab.engine.start_matlab()
#    mateng.generate_ideal_cap_bed2(network_num_hc_long, network_num_hc_high, num_cell_cap, nargout=0)
#    mateng.quit()
    
    # Load the network geometry from the MATLAB file
    ideal_cap_bed = sio.loadmat(in_cap_bed_filename)
    vess_seg = ideal_cap_bed['vess_seg']
    
    filtered_vess_seg = []                                                      # Container of filtered segments
    rp = 8                                                                      # Rounding precision
    
    # Filter out any repeated vessel segments from the raw data
    for seg in vess_seg:
         
        yes_append = True
        
        # Optimize this if possible, must be a better way to do this
        for fvess_seg in filtered_vess_seg:
            
            if (round(seg[0], rp) == round(fvess_seg[0], rp) and round(seg[1], rp) == round(fvess_seg[1], rp) and round(seg[2], rp) == round(fvess_seg[2], rp) and round(seg[3], rp) == round(fvess_seg[3], rp)):
                yes_append = False
                
        if (yes_append == True):
            filtered_vess_seg.append(seg)   
    
    # Create the sub-vessel segments for the ABM based on the input segment size    
    for seg in filtered_vess_seg:
        x0 = float(seg[0] *1e-6)
        y0 = float(seg[1] *1e-6)
        x1 = float(seg[2] *1e-6)
        y1 = float(seg[3] *1e-6)
        
        seg_vect = Vect(x1 - x0, y1 - y0, 0.)
        seg_length = seg_vect.length()
        seg_vect.unit()
        
        num_sub_vess = int(round(seg_length/vess_length))
        seg_vess_length = seg_length/num_sub_vess
        
        for i in range(num_sub_vess):
            node1 = (round(x0 + i*(seg_vess_length*seg_vect.x), rp), round(y0 + i*(seg_vess_length*seg_vect.y), rp))
            node2 = (round(x0 + (i+1)*(seg_vess_length*seg_vect.x), rp), round(y0 + (i+1)*(seg_vess_length*seg_vect.y), rp))
            
            if (node1 not in nodes):
                nodes.append(node1)
                
            if (node2 not in nodes):
                nodes.append(node2)
                
            index1 = nodes.index(node1)
            index2 = nodes.index(node2)
            
            vess = Vessel(len(vessels)+1, index1, index2, int(seg[4]), mu)
                            
            vess.calc_length(node1, node2)
            vess.update_diameter()
            vess.update_conductance()
            
            vessels.append(vess)
    
    # Update the cell number identifiers within all vessel segments
    cell_count = 1
    
    for vess in vessels:        
        for cell in vess.cells:
            cell.ID = cell_count
            cell_count += 1
            
    # Give vessels information about their neighbours
    BC_nodes = []
    
    find_vessel_neighbours(nodes, vessels, BC_nodes)
    
    # Calculate the degree of each node
    for node in nodes:
        node_degree.append(0)
    
    for vess in vessels:
        node_degree[vess.n0] += 1
        node_degree[vess.n1] += 1
    
    # Store the nodes that are bifurcation points
    for i in range(len(node_degree)):
        if (node_degree[i] == 3):
            bifur_nodes.append(i)
    
    # Initialize the bifurcation status array
    for i in range(len(bifur_nodes)):
        bifur.append(0)
        
    # Output network info
#    Output.print_nodes(nodes)
#    Output.print_vessels(vessels)
    
    # Prescribe boundary conditions
    #assert (len(BC_nodes) == len(PBC)), "Number of boundary nodes and number of prescribed pressures is not the same"
    
    x_min, y_min = min(nodes, key=lambda item:item)
    x_max, y_max = max(nodes, key=lambda item:item)   
    
    for i in range(len(nodes)):
        if (nodes[i][0] == x_min):
            BCs.append((i, Part))
        elif (nodes[i][0] == x_max):
            BCs.append((i, Pvein))
        
    # Apply Dirichlet cell conditions at the boundary nodes   
    node_top_left = nodes.index((x_min, y_max))
    node_bottom_right = nodes.index((x_max, y_min))
    
#    if (yes_dirichlet_BCs):
#        for vess in vessels:
#            if (vess.n0 == node_top_left or vess.n1 == node_top_left):
#                vess.dirichlet = vess.num_cells
#            if (vess.n0 == node_bottom_right or vess.n1 == node_bottom_right):
#                vess.dirichlet = vess.num_cells
#    
#    # Apply Dirichlet cell conditions along the whole artery and bein
#    if (yes_art_vein_dirichlet): 
#        for vess in vessels:
#            if (vess.type == "artery" or vess.type == "vein"):
#                vess.dirichlet = vess.num_cells


#------------------------------------------------------------------------------
# Create the vessel network (load the ideal sprout front from the MATLAB file)
def load_ideal_sprout_front(nodes, vessels, node_degree, bifur_nodes, bifur, BCs):
    # Generate the ideal capillary bed using the MATLAB script
    mateng = matlab.engine.start_matlab()
    mateng.generate_ideal_sprout_front(network_num_hc_long, network_num_hc_high, num_cell_art, num_cell_vein, num_cell_cap, nargout=0)
    mateng.quit()
    
    # Load the network geometry from the MATLAB file
    ideal_cap_bed = sio.loadmat(in_cap_bed_filename)
    vess_seg = ideal_cap_bed['vess_seg']
    
    filtered_vess_seg = []                                                      # Container of filtered segments
    rp = 8                                                                      # Rounding precision
    
    # Filter out any repeated vessel segments from the raw data
    for seg in vess_seg:
         
        yes_append = True
        
        # Optimize this if possible, must be a better way to do this
        for fvess_seg in filtered_vess_seg:
            
            if (round(seg[0], rp) == round(fvess_seg[0], rp) and round(seg[1], rp) == round(fvess_seg[1], rp) and round(seg[2], rp) == round(fvess_seg[2], rp) and round(seg[3], rp) == round(fvess_seg[3], rp)):
                yes_append = False
                
        if (yes_append == True):
            filtered_vess_seg.append(seg)   
    
    # Create the sub-vessel segments for the ABM based on the input segment size    
    for seg in filtered_vess_seg:
        x0 = float(seg[0] *1e-6)
        y0 = float(seg[1] *1e-6)
        x1 = float(seg[2] *1e-6)
        y1 = float(seg[3] *1e-6)
        
        seg_vect = Vect(x1 - x0, y1 - y0, 0.)
        seg_length = seg_vect.length()
        seg_vect.unit()
        
        num_sub_vess = int(round(seg_length/vess_length))
        seg_vess_length = seg_length/num_sub_vess
        
        for i in range(num_sub_vess):
            node1 = (round(x0 + i*(seg_vess_length*seg_vect.x), rp), round(y0 + i*(seg_vess_length*seg_vect.y), rp))
            node2 = (round(x0 + (i+1)*(seg_vess_length*seg_vect.x), rp), round(y0 + (i+1)*(seg_vess_length*seg_vect.y), rp))
            
            if (node1 not in nodes):
                nodes.append(node1)
                
            if (node2 not in nodes):
                nodes.append(node2)
                
            index1 = nodes.index(node1)
            index2 = nodes.index(node2)
            
            vess = Vessel(len(vessels)+1, index1, index2, int(seg[4]), mu)
            
            if (vess.num_cells == num_cell_art):
                vess.type = "artery"
            elif (vess.num_cells == num_cell_vein):
                vess.type = "vein"
                
            vess.calc_length(node1, node2)
            vess.update_diameter()
            vess.update_conductance()
            
            vessels.append(vess)
    
    # Update the cell number identifiers within all vessel segments
    cell_count = 1
    
    for vess in vessels:        
        for cell in vess.cells:
            cell.ID = cell_count
            cell_count += 1
            
    # Give vessels information about their neighbours
    BC_nodes = []
    
    find_vessel_neighbours(nodes, vessels, BC_nodes)
    
    # Calculate the degree of each node
    for node in nodes:
        node_degree.append(0)
    
    for vess in vessels:
        node_degree[vess.n0] += 1
        node_degree[vess.n1] += 1
    
    # Store the nodes that are bifurcation points
    for i in range(len(node_degree)):
        if (node_degree[i] == 3):
            bifur_nodes.append(i)
    
    # Initialize the bifurcation status array
    for i in range(len(bifur_nodes)):
        bifur.append(0)
        
    # Output network info
#    Output.print_nodes(nodes)
#    Output.print_vessels(vessels)
    
    # Prescribe boundary conditions
    assert (len(BC_nodes) == len(PBC)), "Number of boundary nodes and number of prescribed pressures is not the same"
    
    for i in range(len(BC_nodes)):
        BCs.append((BC_nodes[i], PBC[i]))
        
    # Apply Dirichlet cell conditions at the boundary nodes
    x_min, y_min = min(nodes, key=lambda item:item)
    x_max, y_max = max(nodes, key=lambda item:item)
    
    node_top_left = nodes.index((x_min, y_max))
    node_bottom_right = nodes.index((x_max, y_min))
    
    if (yes_dirichlet_BCs):
        for vess in vessels:
            if (vess.n0 in BC_nodes):
                vess.dirichlet = vess.num_cells
    
    # Apply Dirichlet cell conditions along the whole artery and bein
    if (yes_art_vein_dirichlet): 
        for vess in vessels:
            if (vess.type == "artery" or vess.type == "vein"):
                vess.dirichlet = vess.num_cells
                
    if (vein_bf_option == 1):
        for vess in vessels:
            if (vess.type == "capillary"):
                if (len(vess.neigh1) == 1):
                    if (vessels[vess.neigh1[0]-1].type == "vein"):
                        if (vessels[vess.neigh1[0]-1].neigh1[0] == vess.ID):
                            vess.dirichlet = vess.num_cells        
                    
                
#------------------------------------------------------------------------------
# Create the vessel network (the original A branch configuration)
def create_old_A_branch_network(nodes, vessels, node_degree, bifur_nodes, bifur, BCs):
    # List of Nodes
    nodes : list[tuple]
    # List of Vessels
    vessels: list[Vessel]
    # List of boundary conditions (known pressures at inlets and outlets)
    BCs :  list[tuple]    
    
    # Create the first vessel section (lower left)
    nodes.append((0, 0))

    for i in range(5):
        node = (0, vess_length*(i+1))
        nodes.append(node)
        
        vess = Vessel(i+1, i, i+1, num_cell, mu)
        vessels.append(vess)
    
    # Create the second vessel section (proximal horizontal branch)
    for i in range(6, 16):
        node = (vess_length*(i-5), 50.e-6)
        nodes.append(node)
        
        vess = Vessel(i, i-1, i, num_cell, mu)
        vessels.append(vess)
    
    # Create the third vessel section (lower right)
    for i in range(16, 21):
        node = (100.e-6, 50.e-6 - vess_length*(i-15))
        nodes.append(node)
        
        vess = Vessel(i, i-1, i, num_cell, mu)
        vessels.append(vess)
    
    # Create the fourth vessel section (upper left)
    for i in range(21, 26):
        node = (0, 50.e-6 + vess_length*(i-20))
        nodes.append(node)
        
        if (i == 21):
            vess = Vessel(i, 5, i, num_cell, mu)
        else:
            vess = Vessel(i, i-1, i, num_cell, mu)
            
        vessels.append(vess)
    
    # Create the fifth vessel segment (distal horizontal branch)
    for i in range(26, 36):
        node = (vess_length*(i-25), 100.e-6)
        nodes.append(node)
        
        vess = Vessel(i, i-1, i, num_cell, mu)
        vessels.append(vess)
    
    # Create the sixth vessel segment (lower right)    
    for i in range(36, 41):
        node = (100.e-6, 100.e-6 - vess_length*(i-35))
        if (i != 40):
            nodes.append(node)
        
        if (i == 40):
            vess = Vessel(i, i-1, 15, num_cell, mu)
        else:
            vess = Vessel(i, i-1, i, num_cell, mu)
    
        vessels.append(vess)
    
    # Print the nodal array     
    #Output.print_nodes(nodes)
    
    # Give vessels information about their neighbours
    BC_nodes = []
    
    
    find_vessel_neighbours(nodes, vessels, BC_nodes)
    
    # Assign the Dirchlet cell number condition to the outlet vessel
    if (A_branch_dirchlet == True):
        vessels[19].dirichlet = vessels[19].num_cells
    
    # Calculate the degree of each node
    for node in nodes:
        node_degree.append(0)
    
    for vess in vessels:
        node_degree[vess.n0] += 1
        node_degree[vess.n1] += 1
    
    # Store the nodes that are bifurcation points
    for i in range(len(node_degree)):
        if (node_degree[i] == 3):
            bifur_nodes.append(i)
    
    # Initialize the bifurcation state array
    for i in range(len(bifur_nodes)):
        bifur.append(0)
        
    # Update the cell number identifiers within all vessel segments
    cell_count = 1
    
    for vess in vessels:
        vess.calc_length(nodes[vess.n0], nodes[vess.n1])
        vess.update_conductance()
        
        for cell in vess.cells:
            cell.ID = cell_count
            cell_count += 1
    
    
    # Determine the vessels with Dirchlet conditions
    # Apply the boundary conditions
    BCs.append((0, Pin))
    BCs.append((20, Pout))
    

#------------------------------------------------------------------------------
# Create the vessel network (the original A branch configuration)
def create_longer_distal_A_branch_network(nodes, vessels, node_degree, bifur_nodes, bifur, BCs):
    # List of Nodes
    nodes : list[tuple]
    # List of Vessels
    vessels: list[Vessel]
    # List of boundary conditions (known pressures at inlets and outlets)
    BCs :  list[tuple]    
    
    # Create the first vessel section (lower left)
    nodes.append((0, 0))

    for i in range(5):
        node = (0, vess_length*(i+1))
        nodes.append(node)
        
        vess = Vessel(i+1, i, i+1, num_cell, mu)
        vessels.append(vess)
    
    # Create the second vessel section (proximal horizontal branch)
    for i in range(6, 16):
        node = (vess_length*(i-5), 50.e-6)
        nodes.append(node)
        
        vess = Vessel(i, i-1, i, num_cell, mu)
        vessels.append(vess)
    
    # Create the third vessel section (lower right)
    for i in range(16, 21):
        node = (100.e-6, 50.e-6 - vess_length*(i-15))
        nodes.append(node)
        
        vess = Vessel(i, i-1, i, num_cell, mu)
        vessels.append(vess)
    
    # Create the fourth vessel section (upper left)
    for i in range(21, 66):
        node = (0, 50.e-6 + vess_length*(i-20))
        nodes.append(node)
        
        if (i == 21):
            vess = Vessel(i, 5, i, num_cell, mu)
        else:
            vess = Vessel(i, i-1, i, num_cell, mu)
            
        vessels.append(vess)
    
    # Create the fifth vessel segment (distal horizontal branch)
    for i in range(66, 76):
        node = (vess_length*(i-65), 500.e-6)
        nodes.append(node)
        
        vess = Vessel(i, i-1, i, num_cell, mu)
        vessels.append(vess)
    
    # Create the sixth vessel segment (lower right)    
    for i in range(76, 121):
        node = (100.e-6, 500.e-6 - vess_length*(i-75))
        if (i != 120):
            nodes.append(node)
        
        if (i == 120):
            vess = Vessel(i, i-1, 15, num_cell, mu)
        else:
            vess = Vessel(i, i-1, i, num_cell, mu)
    
        vessels.append(vess)
    
    # Print the nodal array     
    #Output.print_nodes(nodes)
    
    # Give vessels information about their neighbours
    BC_nodes = []
    
    
    find_vessel_neighbours(nodes, vessels, BC_nodes)
    
    # Assign the Dirchlet cell number condition to the outlet vessel
    if (A_branch_dirchlet == True):
        vessels[19].dirichlet = vessels[19].num_cells
    
    # Calculate the degree of each node
    for node in nodes:
        node_degree.append(0)
    
    for vess in vessels:
        node_degree[vess.n0] += 1
        node_degree[vess.n1] += 1
    
    # Store the nodes that are bifurcation points
    for i in range(len(node_degree)):
        if (node_degree[i] == 3):
            bifur_nodes.append(i)
    
    # Initialize the bifurcation state array
    for i in range(len(bifur_nodes)):
        bifur.append(0)
        
    # Update the cell number identifiers within all vessel segments
    cell_count = 1
    
    for vess in vessels:
        vess.calc_length(nodes[vess.n0], nodes[vess.n1])
        vess.update_conductance()
        
        for cell in vess.cells:
            cell.ID = cell_count
            cell_count += 1
    
    
    # Determine the vessels with Dirchlet conditions
    # Apply the boundary conditions
    BCs.append((0, Pin))
    BCs.append((20, Pout))
    

#------------------------------------------------------------------------------
# Create the vessel network (the smaller A branch configuration)
def create_small_A_branch_network1(nodes, vessels, node_degree, bifur_nodes, bifur, BCs):
    # List of Nodes
    nodes : list[tuple]
    # List of Vessels
    vessels: list[Vessel]
    # List of boundary conditions (known pressures at inlets and outlets)
    BCs :  list[tuple]    
    
    # Create the first vessel section (lower left)
    nodes.append((0, 0))

    for i in range(4):
        node = (0, vess_length*(i+1))
        nodes.append(node)
        
        vess = Vessel(i+1, i, i+1, num_cell, mu)
        vessels.append(vess)
    
    # Create the second vessel section (proximal horizontal branch)
    for i in range(5, 13):
        node = (vess_length*(i-4), 40.e-6)
        nodes.append(node)
        
        vess = Vessel(i, i-1, i, num_cell, mu)
        vessels.append(vess)
    
    # Create the third vessel section (lower right)
    for i in range(13, 17):
        node = (80.e-6, 40.e-6 - vess_length*(i-12))
        nodes.append(node)
        
        vess = Vessel(i, i-1, i, num_cell, mu)
        vessels.append(vess)
    
    # Create the fourth vessel section (upper left)
    for i in range(17, 21):
        node = (0, 40.e-6 + vess_length*(i-16))
        nodes.append(node)
        
        if (i == 17):
            vess = Vessel(i, 4, i, num_cell, mu)
        else:
            vess = Vessel(i, i-1, i, num_cell, mu)
            
        vessels.append(vess)
    
    # Create the fifth vessel segment (distal horizontal branch)
    for i in range(21, 29):
        node = (vess_length*(i-20), 80.e-6)
        nodes.append(node)
        
        vess = Vessel(i, i-1, i, num_cell, mu)
        vessels.append(vess)
    
    # Create the sixth vessel segment (lower right)    
    for i in range(29, 33):
        node = (80.e-6, 80.e-6 - vess_length*(i-28))
        if (i != 32):
            nodes.append(node)
        
        if (i == 32):
            vess = Vessel(i, i-1, 12, num_cell, mu)
        else:
            vess = Vessel(i, i-1, i, num_cell, mu)
    
        vessels.append(vess)
    
    # Print the nodal array     
    #Output.print_nodes(nodes)
    
    # Give vessels information about their neighbours
    BC_nodes = []
    
    
    find_vessel_neighbours(nodes, vessels, BC_nodes)
    
    # Calculate the degree of each node
    for node in nodes:
        node_degree.append(0)
    
    for vess in vessels:
        node_degree[vess.n0] += 1
        node_degree[vess.n1] += 1
    
    # Store the nodes that are bifurcation points
    for i in range(len(node_degree)):
        if (node_degree[i] == 3):
            bifur_nodes.append(i)
    
    # Initialize the bifurcation state array
    for i in range(len(bifur_nodes)):
        bifur.append(0)
        
    # Update the cell number identifiers within all vessel segments
    cell_count = 1
    
    for vess in vessels:
        vess.calc_length(nodes[vess.n0], nodes[vess.n1])
        vess.update_conductance()
        
        for cell in vess.cells:
            cell.ID = cell_count
            cell_count += 1
    
    
    # Determine the vessels with Dirchlet conditions
    # Apply the boundary conditions
    BCs.append((0, Pin))
    BCs.append((16, Pout))
    
    
#------------------------------------------------------------------------------
# Create the vessel network (the smaller A branch configuration)
def create_small_A_branch_network2(nodes, vessels, node_degree, bifur_nodes, bifur, BCs):
    # List of Nodes
    nodes : list[tuple]
    # List of Vessels
    vessels: list[Vessel]
    # List of boundary conditions (known pressures at inlets and outlets)
    BCs :  list[tuple]    
    
    # Create the first vessel section (lower left)
    nodes.append((0, 0))

    for i in range(3):
        node = (0, vess_length*(i+1))
        nodes.append(node)
        
        vess = Vessel(i+1, i, i+1, num_cell, mu)
        vessels.append(vess)
    
    # Create the second vessel section (proximal horizontal branch)
    for i in range(4, 10):
        node = (vess_length*(i-3), 30.e-6)
        nodes.append(node)
        
        vess = Vessel(i, i-1, i, num_cell, mu)
        vessels.append(vess)
    
    # Create the third vessel section (lower right)
    for i in range(10, 13):
        node = (60.e-6, 30.e-6 - vess_length*(i-9))
        nodes.append(node)
        
        vess = Vessel(i, i-1, i, num_cell, mu)
        vessels.append(vess)
    
    # Create the fourth vessel section (upper left)
    for i in range(13, 16):
        node = (0, 30.e-6 + vess_length*(i-12))
        nodes.append(node)
        
        if (i == 13):
            vess = Vessel(i, 3, i, num_cell, mu)
        else:
            vess = Vessel(i, i-1, i, num_cell, mu)
            
        vessels.append(vess)
    
    # Create the fifth vessel segment (distal horizontal branch)
    for i in range(16, 22):
        node = (vess_length*(i-15), 60.e-6)
        nodes.append(node)
        
        vess = Vessel(i, i-1, i, num_cell, mu)
        vessels.append(vess)
    
    # Create the sixth vessel segment (lower right)    
    for i in range(22, 25):
        node = (60.e-6, 60.e-6 - vess_length*(i-21))
        if (i != 24):
            nodes.append(node)
        
        if (i == 24):
            vess = Vessel(i, i-1, 9, num_cell, mu)
        else:
            vess = Vessel(i, i-1, i, num_cell, mu)
    
        vessels.append(vess)
    
    # Print the nodal array     
    #Output.print_nodes(nodes)
    
    # Give vessels information about their neighbours
    BC_nodes = []
    
    
    find_vessel_neighbours(nodes, vessels, BC_nodes)
    
    # Calculate the degree of each node
    for node in nodes:
        node_degree.append(0)
    
    for vess in vessels:
        node_degree[vess.n0] += 1
        node_degree[vess.n1] += 1
    
    # Store the nodes that are bifurcation points
    for i in range(len(node_degree)):
        if (node_degree[i] == 3):
            bifur_nodes.append(i)
    
    # Initialize the bifurcation state array
    for i in range(len(bifur_nodes)):
        bifur.append(0)
        
    # Update the cell number identifiers within all vessel segments
    cell_count = 1
    
    for vess in vessels:
        vess.calc_length(nodes[vess.n0], nodes[vess.n1])
        vess.update_conductance()
        
        for cell in vess.cells:
            cell.ID = cell_count
            cell_count += 1
    
    
    # Determine the vessels with Dirchlet conditions
    # Apply the boundary conditions
    BCs.append((0, Pin))
    BCs.append((12, Pout))
    
 
#------------------------------------------------------------------------------
# Create the vessel network (the original A branch configuration)
def create_Y_branch_network(nodes, vessels, node_degree, bifur_nodes, bifur, BCs):
    # List of Nodes
    nodes : list[tuple]
    # List of Vessels
    vessels: list[Vessel]
    # List of boundary conditions (known pressures at inlets and outlets)
    BCs :  list[tuple]    
    
    # Create the first vessel section (lower left)
    nodes.append((0, 0))

    v2 = Vect(cos(radians(45)), sin(radians(45)), 0.)
    v3 = Vect(-cos(radians(45)), sin(radians(45)), 0.)
    
    for i in range(10):
        node = (0, vess_length*(i+1))
        nodes.append(node)
        
        vess = Vessel(i+1, i, i+1, num_cell, mu)
        vessels.append(vess)
        
    # Create the left branch
    for i in range(11, 21):
        node = (vess_length*(i-10)*v2.x, (10*vess_length) + vess_length*(i-10)*v2.y)
        nodes.append(node)
        
        vess = Vessel(i, i-1, i, num_cell, mu)
        vessels.append(vess)
        
    # Create the right branch
    for i in range(21, 31):
        node = (vess_length*(i-20)*v3.x, (10*vess_length) + vess_length*(i-20)*v3.y)
        nodes.append(node)
        
        if (i == 21):
            vess = Vessel(i, 10, i, num_cell, mu)
        else:
            vess = Vessel(i, i-1, i, num_cell, mu)
        
        vessels.append(vess)

    # Print the nodal array     
    #Output.print_nodes(nodes)
    
    # Give vessels information about their neighbours
    BC_nodes = []
      
    find_vessel_neighbours(nodes, vessels, BC_nodes)
    
    # Assign the Dirchlet cell number condition to the outlet vessel
    if (Y_branch_dirchlet == True):
        vessels[0].dirichlet = vessels[0].num_cells
    
    # Calculate the degree of each node
    for node in nodes:
        node_degree.append(0)
    
    for vess in vessels:
        node_degree[vess.n0] += 1
        node_degree[vess.n1] += 1
    
    # Store the nodes that are bifurcation points
    for i in range(len(node_degree)):
        if (node_degree[i] == 3):
            bifur_nodes.append(i)
    
    # Initialize the bifurcation state array
    for i in range(len(bifur_nodes)):
        bifur.append(0)
        
    # Update the cell number identifiers within all vessel segments
    cell_count = 1
    
    for vess in vessels:
        vess.calc_length(nodes[vess.n0], nodes[vess.n1])
        vess.update_conductance()
        
        for cell in vess.cells:
            cell.ID = cell_count
            cell_count += 1
    
    # Determine the vessels with Dirchlet conditions
    # Apply the boundary conditions
    BCs.append((0, Pout))
    BCs.append((20, Pright))
    BCs.append((30, Pleft))


#------------------------------------------------------------------------------
# Create the vessel network (the original A branch configuration)
def create_Y_branch_network2(nodes, vessels, node_degree, bifur_nodes, bifur, BCs):
    # List of Nodes
    nodes : list[tuple]
    # List of Vessels
    vessels: list[Vessel]
    # List of boundary conditions (known pressures at inlets and outlets)
    BCs :  list[tuple]    
    
    # Create the first vessel section (lower left)
    nodes.append((0, 0))

    v2 = Vect(cos(radians(45)), sin(radians(45)), 0.)
    v3 = Vect(-cos(radians(45)), sin(radians(45)), 0.)
    
    for i in range(10):
        node = (0, vess_length*(i+1))
        nodes.append(node)
        
        vess = Vessel(i+1, i, i+1, num_cell, mu)
        vessels.append(vess)
        
    # Create the left branch
    for i in range(11, 31):
        node = (vess_length*(i-10)*v2.x, 100.e-6 + vess_length*(i-10)*v2.y)
        nodes.append(node)
        
        vess = Vessel(i, i-1, i, num_cell, mu)
        vessels.append(vess)
        
    # Create the right branch
    for i in range(31, 41):
        node = (vess_length*(i-30)*v3.x, 100.e-6 + vess_length*(i-30)*v3.y)
        nodes.append(node)
        
        if (i == 31):
            vess = Vessel(i, 10, i, num_cell, mu)
        else:
            vess = Vessel(i, i-1, i, num_cell, mu)
        
        vessels.append(vess)

    # Print the nodal array     
    #Output.print_nodes(nodes)
    
    # Give vessels information about their neighbours
    BC_nodes = []
      
    find_vessel_neighbours(nodes, vessels, BC_nodes)
    
    # Calculate the degree of each node
    for node in nodes:
        node_degree.append(0)
    
    for vess in vessels:
        node_degree[vess.n0] += 1
        node_degree[vess.n1] += 1
    
    # Store the nodes that are bifurcation points
    for i in range(len(node_degree)):
        if (node_degree[i] == 3):
            bifur_nodes.append(i)
    
    # Initialize the bifurcation state array
    for i in range(len(bifur_nodes)):
        bifur.append(0)
        
    # Update the cell number identifiers within all vessel segments
    cell_count = 1
    
    for vess in vessels:
        vess.calc_length(nodes[vess.n0], nodes[vess.n1])
        vess.update_conductance()
        
        for cell in vess.cells:
            cell.ID = cell_count
            cell_count += 1
    
    # Determine the vessels with Dirchlet conditions
    # Apply the boundary conditions
    BCs.append((0, Pout))
    BCs.append((30, Pright))
    BCs.append((40, Pleft))
    
    
#------------------------------------------------------------------------------
# Create the vessel network (the original A branch configuration)
def create_small_Y_branch_network(nodes, vessels, node_degree, bifur_nodes, bifur, BCs):
    # List of Nodes
    nodes : list[tuple]
    # List of Vessels
    vessels: list[Vessel]
    # List of boundary conditions (known pressures at inlets and outlets)
    BCs :  list[tuple]    
    
    # Create the first vessel section (lower left)
    nodes.append((0, 0))

    v2 = Vect(cos(radians(45)), sin(radians(45)), 0.)
    v3 = Vect(-cos(radians(45)), sin(radians(45)), 0.)
    
    for i in range(5):
        node = (0, vess_length*(i+1))
        nodes.append(node)
        
        vess = Vessel(i+1, i, i+1, num_cell, mu)
        vessels.append(vess)
        
    # Create the left branch
    for i in range(6, 11):
        node = (vess_length*(i-5)*v2.x, 50.e-6 + vess_length*(i-5)*v2.y)
        nodes.append(node)
        
        vess = Vessel(i, i-1, i, num_cell, mu)
        vessels.append(vess)
        
    # Create the right branch
    for i in range(11, 16):
        node = (vess_length*(i-10)*v3.x, 50.e-6 + vess_length*(i-10)*v3.y)
        nodes.append(node)
        
        if (i == 11):
            vess = Vessel(i, 5, i, num_cell, mu)
        else:
            vess = Vessel(i, i-1, i, num_cell, mu)
        
        vessels.append(vess)

    # Print the nodal array     
    #Output.print_nodes(nodes)
    
    # Give vessels information about their neighbours
    BC_nodes = []
      
    find_vessel_neighbours(nodes, vessels, BC_nodes)
    
    # Calculate the degree of each node
    for node in nodes:
        node_degree.append(0)
    
    for vess in vessels:
        node_degree[vess.n0] += 1
        node_degree[vess.n1] += 1
    
    # Store the nodes that are bifurcation points
    for i in range(len(node_degree)):
        if (node_degree[i] == 3):
            bifur_nodes.append(i)
    
    # Initialize the bifurcation state array
    for i in range(len(bifur_nodes)):
        bifur.append(0)
        
    # Update the cell number identifiers within all vessel segments
    cell_count = 1
    
    for vess in vessels:
        vess.calc_length(nodes[vess.n0], nodes[vess.n1])
        vess.update_conductance()
        
        for cell in vess.cells:
            cell.ID = cell_count
            cell_count += 1
    
    # Determine the vessels with Dirchlet conditions
    # Apply the boundary conditions
    BCs.append((0, Pout))
    BCs.append((10, Pright))
    BCs.append((15, Pleft))
    

#------------------------------------------------------------------------------
# Give vessels information on their upstream and downstream neighbours
def find_vessel_neighbours(nodes, vessels, BC_nodes):

    for vess in vessels:
        for vess2 in vessels:
            if (vess.ID != vess2.ID):
                if ((vess.n0 == vess2.n0) or (vess.n0 == vess2.n1)):
                    vess.neigh0.append(vess2.ID)
                    
                if ((vess.n1 == vess2.n0) or (vess.n1 == vess2.n1)):
                    vess.neigh1.append(vess2.ID)
    
    # Apply periodic conditions at the boundary
    free_n0 = []
    free_n1 = []
    
    for vess in vessels:
        if len(vess.neigh0) == 0:
            free_n0.append(vess.ID)
            
        if len(vess.neigh1) == 0:
            free_n1.append(vess.ID)
    
    if ((sim_type != "ideal sprout front") and (sim_type != "Y branch")):
        assert (len(free_n0) == len(free_n1)), "Periodic Boundary Error: Number of free ends is not equal"
        
        for i in range(len(free_n0)):
            vessels[free_n0[i]-1].neigh0.append(free_n1[i])
            vessels[free_n1[i]-1].neigh1.append(free_n0[i])
    elif (sim_type == "ideal sprout front"):
        vessels[free_n0[0]-1].neigh0.append(free_n0[1])
        vessels[free_n0[1]-1].neigh0.append(free_n0[0])
    elif (sim_type == "Y branch"):
        vessels[free_n0[0]-1].neigh0.append(free_n1[0])
        vessels[free_n0[0]-1].neigh0.append(free_n1[1])
        vessels[free_n1[0]-1].neigh1.append(free_n0[0])
        vessels[free_n1[1]-1].neigh1.append(free_n0[0])
        
    for n0 in free_n0:
        BC_nodes.append(vessels[n0-1].n0)
        
    for n1 in free_n1:
        BC_nodes.append(vessels[n1-1].n1)
    
    BC_nodes.sort()
    
    # Print info on neighbours to console for verification            
    #Output.print_vess_neighbours(vessels)