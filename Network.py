#------------------------------------------------------------------------------
# The Network class contains the methods that require knowledge of the vascular network as a whole
#
# Lowell Taylor Edgar
# University of Edinburgh
# 2019

from math import pi
from math import acos
from math import cos
from math import sin

from numpy import sign
import random as rand

from Vessel import *
from Tensor3D import *
from Input import *
from ECell import *


#------------------------------------------------------------------------------    
# Realign polarity vectors in cells based on weights
def realign_polarity(vessels):
    
    for vess in vessels:
        for cell in vess.cells:
            
            # Original polarity vector
            pol_old = cell.polarity.copy()
            
            # Flow vector
            flow_vect = -float(sign(vess.Q))*vess.unit.copy()
            
            # Random walk component
            rand_vect = Vect(rand.uniform(-1, 1), rand.uniform(-1, 1), 0.)
            rand_vect.unit()
            
#            if (flow_vect.length() == 0.):
#                if (vessels[vess.n1-1].num_cells == 0):
#                    flow_vect = -vess.unit.copy()
#                    
#                if (vessels[vess.n0-1].num_cells == 0):
#                    flow_vect = vess.unit.copy()
#
#                if (vessels[vess.n1-1].num_cells == 0 and vessels[vess.n0-1].num_cells == 0):
#                    flow_vect = -vess.unit.copy()
                    
            # Phi1 - Realignment angle due to persistence
            phi1 = 0.
            
            # Phi2 - Realignment angle due to flow
            phi2_dot = pol_old * flow_vect
            
            if (phi2_dot > 1.):
                phi2_dot = 1.
            elif (phi2_dot < -1.):
                phi2_dot = -1.
            
            phi2 = acos(phi2_dot)
            
            if (pol_old.cross(flow_vect).z < 0):
                phi2 = -phi2            
            
            # Phi3 - Realignment angle due to random walk component
            phi3_dot = pol_old * rand_vect
            
            if (phi3_dot > 1.):
                phi3_dot = 1.
            elif (phi3_dot < -1.):
                phi3_dot = -1.
            
            phi3 = acos(phi3_dot)
            
            if (pol_old.cross(rand_vect).z < 0):
                phi3 = -phi3
            
#            # If shear weighting of polarization is enabled
#            if (yes_shear_polar):
#                if (vess.tau <= tau_max):
#                    eta = vess.tau/tau_max
#                else:
#                    eta = 1.0
#                    
#                w2 = (1 - w1)*eta
#                w3 = (1 - w1)*(1 - eta)
                
            
            # Calculate new realignment angle from weighted average
            theta = w1*phi1 + w2*phi2 + w3*phi3
            
            # Calculate the new polarity vector by rotating the old vector by theta
            Q = Tensor2O([cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 0])
            
            pol_new = Q*pol_old
            pol_new.unit()
            
            cell.polarity = pol_new.copy()
            

#------------------------------------------------------------------------------            
# Perform cell migration
def cell_migration(vessels, bifur_nodes, bifur):
    
    # Determine if a cell is migrating out of a vessel segment by checking the 
    # angle between polarity and the vessel's unit vector
    for vess in vessels:
        for cell in vess.cells:
            dot_prod = vess.unit * cell.polarity
            
            if (dot_prod > 1.):
                dot_prod = 1.
            elif (dot_prod < -1.):
                dot_prod = -1.
                
            pol_angle = acos(dot_prod)
            
            cell.migrate = 0
            
            if (0. <= pol_angle <= pi/3.):
                cell.migrate = -1                                               # Migrate along the vessel's unit vector
            elif ((2*pi)/3. <= pol_angle <= pi):
                cell.migrate = 1                                                # Migrate against the vessel's unit vector
    
    
    # Diffusion-based smoothing scheme
    if (yes_smoothing == True):
        vessel_smoothing(vessels)
    
    
    # Move the cells that are migrating to their new vessel        
    for vess in vessels:
        cells_that_migrated = []
        
        for i in range(len(vess.cells)):
            
            # If cell is migrating upstream, move to upstream neighbour
            if (vess.cells[i].migrate == 1):
                migrating_cell = vess.cells[i]
                
                if (len(vess.neigh0) == 1):
                    migrating_cell.vessID = vess.neigh0[0]
                    vessels[vess.neigh0[0]-1].cells.append(migrating_cell)
                
                # If more than one neighbour...
                elif (len(vess.neigh0) == 2):
                    if (yes_new_branch == False):
                        handle_bifurcation(migrating_cell, vess, vess.n0, vessels[vess.neigh0[0]-1], vessels[vess.neigh0[1]-1])
                    else:
                        new_handle_bifurcation(migrating_cell, vess, vess.n0, vessels[vess.neigh0[0]-1], vessels[vess.neigh0[1]-1], bifur_nodes, bifur)
                    
                if (migrating_cell.migrate != 0):
                    cells_that_migrated.append(i)
                    migrating_cell.migrate = 0
            
            
            # If cell is migrating downstream, move to downstream neighbour
            if (vess.cells[i].migrate == -1):
                migrating_cell = vess.cells[i]
                
                if (len(vess.neigh1) == 1):
                    migrating_cell.vessID = vess.neigh1[0]
                    vessels[vess.neigh1[0]-1].cells.append(migrating_cell)
                               
                # If more than one neighbour...
                elif (len(vess.neigh1) == 2):
                    if (yes_new_branch == False):
                        handle_bifurcation(migrating_cell, vess, vess.n1, vessels[vess.neigh1[0]-1], vessels[vess.neigh1[1]-1])
                    else:
                        new_handle_bifurcation(migrating_cell, vess, vess.n1, vessels[vess.neigh1[0]-1], vessels[vess.neigh1[1]-1], bifur_nodes, bifur)
                
                if (migrating_cell.migrate != 0):
                    cells_that_migrated.append(i)
                    migrating_cell.migrate = 0
                    
                    
        # Remove the cells that left the vessel (or rather, keep only the cells that didn't leave)
        new_cells = []
        
        for i in range(len(vess.cells)):
            if i not in cells_that_migrated:
                new_cells.append(vess.cells[i])
                
        vess.cells = new_cells

    # Update cell number in vessels after migration
    for vess in vessels:
        vess.num_cells = len(vess.cells)
        
    # Apply Dirichlet boundary conditions
    max_cell_num = 0
    
    # Determine the max cell ID number for adding new cells
    for vess in vessels:
        for cell in vess.cells:
            if (cell.ID > max_cell_num):
                max_cell_num = cell.ID
    
    for vess in vessels:
        if (vess.dirichlet != 0 and vess.num_cells != vess.dirichlet):
            # If cell number is less than dirichlet condition, add cells
            if (vess.num_cells < vess.dirichlet):
                while (vess.num_cells < vess.dirichlet):
                    max_cell_num += 1
                    vess.cells.append(ECell(max_cell_num, vess.ID, Vect(0, 0, 0)))
                    vess.num_cells = len(vess.cells)
                    pass
            # If cell number is more than dirichlet condition, remove cells randomly
            elif (vess.num_cells > vess.dirichlet):
                while (vess.num_cells > vess.dirichlet):
                    vess.cells.pop(rand.randrange(len(vess.cells)))
                    vess.num_cells = len(vess.cells)
                    pass
    
    # Update vessels after migration
    for vess in vessels:
        vess.num_cells = len(vess.cells)
        vess.update_diameter()
        vess.update_conductance()


#------------------------------------------------------------------------------
# Determine cell behavior at a bifurcation
def vessel_smoothing(vessels):        
    for vess in vessels:
                
        # If flow in the vessel is positive, check downstream neighbour
        if (vess.Q > 0.):
            # If current vessel has only one downstream neighbour
            if (len(vess.neigh1) == 1):
                if (vessels[vess.neigh1[0]-1].num_cells < vess.num_cells) and (vessels[vess.neigh1[0]-1].num_cells != 0):
                #if (vessels[vess.neigh1[0]-1].num_cells < vess.num_cells):
                    
                    smoothed = False
                    
                    while (smoothed == False):
                        rand_cell = rand.randint(0, vess.num_cells-1)
                    
                        if (vess.cells[rand_cell].migrate != 0):
                            vess.cells[rand_cell].migrate = 0
                            smoothed = True                  
                            
#            # If current vessel has more than one downstream neighbour
#            if (len(vess.neigh1) == 2):
#                if ((vessels[vess.neigh1[0]-1].num_cells < vess.num_cells) and (vessels[vess.neigh1[0]-1].num_cells != 0)) or ((vessels[vess.neigh1[1]-1].num_cells < vess.num_cells) and (vessels[vess.neigh1[1]-1].num_cells != 0)):
#                    
#                    smoothed = False
#                    
#                    while (smoothed == False):
#                        rand_cell = rand.randint(0, vess.num_cells-1)
#                    
#                        if (vess.cells[rand_cell].migrate != 0):
#                            vess.cells[rand_cell].migrate = 0
#                            smoothed = True 
                            
        # If flow in the vessel is negative, check upstream neighbour                   
        if (vess.Q < 0.):
            # If current vessel has only one upstream neighbour
            if (len(vess.neigh0) == 1):
                if (vessels[vess.neigh0[0]-1].num_cells < vess.num_cells) and (vessels[vess.neigh0[0]-1].num_cells != 0):
                #if (vessels[vess.neigh0[0]-1].num_cells < vess.num_cells):
                     
                    smoothed = False
                    
                    while (smoothed == False):
                        rand_cell = rand.randint(0, vess.num_cells-1)
                    
                        if (vess.cells[rand_cell].migrate != 0):
                            vess.cells[rand_cell].migrate = 0
                            smoothed = True                  
                            
#            # If current vessel has more than one upstream neighbour
#            if (len(vess.neigh0) == 2):
#                if ((vessels[vess.neigh0[0]-1].num_cells < vess.num_cells) and (vessels[vess.neigh0[0]-1].num_cells != 0)) or ((vessels[vess.neigh0[1]-1].num_cells < vess.num_cells) and (vessels[vess.neigh0[1]-1].num_cells != 0)):
#                    
#                    smoothed = False
#                    
#                    while (smoothed == False):
#                        rand_cell = rand.randint(0, vess.num_cells-1)
#                    
#                        if (vess.cells[rand_cell].migrate != 0):
#                            vess.cells[rand_cell].migrate = 0
#                            smoothed = True      
        
        
#------------------------------------------------------------------------------
# Handle a bifurcation point while migrating
def handle_bifurcation(migrating_cell, parent_vess, bif_node, branch1, branch2):
    # Determine if flow within branches are incoming or outgoing of the bifurcation
    if (branch1.n1 == bif_node):
        if (branch1.Q > 0):
            branch1_status = "incoming"
        else:
            branch1_status = "outgoing"
            
    if (branch1.n0 == bif_node):
        if (branch1.Q > 0):
            branch1_status = "outgoing"
        else:
            branch1_status = "incoming" 
            
    if (branch2.n1 == bif_node):
        if (branch2.Q > 0):
            branch2_status = "incoming"
        else:
            branch2_status = "outgoing"
            
    if (branch2.n0 == bif_node):
        if (branch2.Q > 0):
            branch2_status = "outgoing"
        else:
            branch2_status = "incoming" 
    
    # If flow is only incoming in one branch, always choose that branch (ie, against flow)
    if (branch1_status == "incoming" and branch2_status == "outgoing"):
        migrating_cell.vessID = branch1.ID
        branch1.cells.append(migrating_cell)
        return
    
    if (branch2_status == "incoming" and branch1_status == "outgoing"):
        migrating_cell.vessID = branch2.ID
        branch2.cells.append(migrating_cell)
        return
    
    # If a flow-converging, cell-diverging bifurcation, randomly pick one
    if (branch1_status == "incoming" and branch2_status == "incoming"):        
        #bifurcation_rule(migrating_cell, parent_vess, branch1, branch2)
        new_bifurcation_rule(migrating_cell, parent_vess, branch1, branch2)
        return


#------------------------------------------------------------------------------
# Determine cell behavior using the specified bifurcation rule
def bifurcation_rule(migrating_cell, parent_vess, branch1, branch2):
    
    # Bifurcation Rule 1 - Choose branch with the highest flow
    if (branch_rule == 1):
        if (abs(branch1.Q) > abs(branch2.Q)):
            migrating_cell.vessID = branch1.ID
            branch1.cells.append(migrating_cell)
            return
        else:
            migrating_cell.vessID = branch2.ID
            branch2.cells.append(migrating_cell)
            return
    
    
    # Bifurcation Rule 3 - Choose one branch at random, equal probability        
    if (branch_rule == 3):
        if (rand.uniform(0, 1) > 0.5):
            migrating_cell.vessID = branch1.ID
            branch1.cells.append(migrating_cell)
            return
        else:
            migrating_cell.vessID = branch2.ID
            branch2.cells.append(migrating_cell)
            return
    
    
    # Bifurcation Rule 4 - Assign probability based on number of cells        
    if (branch_rule == 4):
        if (branch1.num_cells < parent_vess.num_cells) and (branch2.num_cells < parent_vess.num_cells):
            prob1 = (branch1.num_cells/parent_vess.num_cells)
        else:
            prob1 = (branch1.num_cells/(branch1.num_cells + branch2.num_cells))
        
        if (rand.uniform(0, 1) < prob1):
            migrating_cell.vessID = branch1.ID
            branch1.cells.append(migrating_cell)
            return
        else:
            migrating_cell.vessID = branch2.ID
            branch2.cells.append(migrating_cell)
            return
        
        
    # Bifurcation Rule 6 - Weighted average of flow and cell number
    if (branch_rule == 6):
        Q_ratio = branch2.Q/(branch1.Q + branch2.Q)
        cell_ratio = branch2.num_cells/(branch1.num_cells + branch2.num_cells)
        
        prob2 = branch_alpha*Q_ratio + (1 - branch_alpha)*cell_ratio
        
        if (rand.uniform(0, 1) > prob2):
            migrating_cell.vessID = branch1.ID
            branch1.cells.append(migrating_cell)
            return
        else:
            migrating_cell.vessID = branch2.ID
            branch2.cells.append(migrating_cell)
            return
        
        
    # Bifurcation Rule 7 - Weighted average of WSS and cell number
    if (branch_rule == 7):
        if ((branch1.tau + branch2.tau) != 0. and (branch1.num_cells + branch2.num_cells) != 0.):
        
            tau_ratio = branch2.tau/(branch1.tau + branch2.tau)
            cell_ratio = branch2.num_cells/(branch1.num_cells + branch2.num_cells)
            
            tau1 = branch1.tau/(branch1.tau + branch2.tau)
            cell1 = branch1.num_cells/(branch1.num_cells + branch2.num_cells)
            
            tau2 = branch2.tau/(branch1.tau + branch2.tau)
            cell2 = branch2.num_cells/(branch1.num_cells + branch2.num_cells)
            
            #prob1 = branch_alpha*tau1 + (1 - branch_alpha)*cell1
            #prob2 = branch_alpha*tau2 + (1 - branch_alpha)*cell2
            
            prob1 = branch_alpha*tau1 + (1 - branch_alpha)*cell2
            prob2 = branch_alpha*tau2 + (1 - branch_alpha)*cell1
            
            assert (round(prob1, 8) == round((1 - prob2), 8)), "p1 = {0} does not match {1}".format(prob1, (1 - prob2))
            
            if (rand.uniform(0, 1) > prob2):
                migrating_cell.vessID = branch1.ID
                branch1.cells.append(migrating_cell)
                return
            else:
                migrating_cell.vessID = branch2.ID
                branch2.cells.append(migrating_cell)
                return   
    
    
#------------------------------------------------------------------------------
# Handle a bifurcation point while migrating (new)
def new_handle_bifurcation(migrating_cell, parent_vess, bif_node, branch1, branch2, bifur_nodes, bifur):
    
    # Set the minimum value of flow, less than this and cells won't consider entering that branch
    flow_min = (1e-4)/(3.6e12)
       
    # Determine if flow within branches are incoming or outgoing of the bifurcation
    if (branch1.n0 == bif_node):
        if (abs(branch1.Q) < flow_min):
            branch1_status = "zero"
        else:
            if (branch1.Q > 0):
                branch1_status = "outgoing"
            else:
                branch1_status = "incoming" 
    
    
    if (branch1.n1 == bif_node):
        if (abs(branch1.Q) < flow_min):
            branch1_status = "zero"
        else:
            if (branch1.Q > 0):
                branch1_status = "incoming"
            else:
                branch1_status = "outgoing"
            

    if (branch2.n0 == bif_node):
        if (abs(branch2.Q) < flow_min):
            branch2_status = "zero"
        else:
            if (branch2.Q > 0):
                branch2_status = "outgoing"
            else:
                branch2_status = "incoming" 
    
    
    if (branch2.n1 == bif_node):
        if (abs(branch2.Q) < flow_min):
            branch2_status = "zero"
        else:
            if (branch2.Q > 0):
                branch2_status = "incoming"
            else:
                branch2_status = "outgoing"
            

    # If flow is zero in one of the branches, chose the other one
    if (branch1_status != "zero" and branch2_status == "zero"):
        migrating_cell.vessID = branch1.ID
        branch1.cells.append(migrating_cell)
        return
    
    if (branch1_status == "zero" and branch2_status != "zero"):
        migrating_cell.vessID = branch2.ID
        branch2.cells.append(migrating_cell)
        return
    
    # If flow is only incoming in one branch, always choose that branch (ie, against flow)
    if (branch1_status == "incoming" and branch2_status == "outgoing"):
        migrating_cell.vessID = branch1.ID
        branch1.cells.append(migrating_cell)
        return
    
    if (branch2_status == "incoming" and branch1_status == "outgoing"):
        migrating_cell.vessID = branch2.ID
        branch2.cells.append(migrating_cell)
        return
    
    # If a flow-converging, cell-diverging bifurcation, enact bifurcation rules
    if (branch1_status == "incoming" and branch2_status == "incoming"):
        new_bifurcation_rule(migrating_cell, parent_vess, branch1, branch2)
                    
                    
#------------------------------------------------------------------------------
# Determine cell behavior using the specified bifurcation rule (new)
def new_bifurcation_rule(migrating_cell, parent_vess, branch1, branch2):
    
    # Handle the case of capillaries splitting off from the vein
    if (parent_vess.type == "vein" and branch1.type == "vein") or (parent_vess.type == "vein" and branch2.type == "vein"):
            
            # Vein bifurcation option 1 - Probabilty of splitting off vein given by initial capillary and vein size
            if (vein_bf_option == 1):
                if (branch1.type == "capillary"):
                    prob1 = num_cell_cap/num_cell_vein
                    
                    if (rand.uniform(0,1) <= prob1):
                        migrating_cell.vessID = branch1.ID
                        branch1.cells.append(migrating_cell)
                        return
                    
                if (branch2.type == "capillary"):
                    prob2 = num_cell_cap/num_cell_vein
                    
                    if (rand.uniform(0,1) <= prob2):
                        migrating_cell.vessID = branch2.ID
                        branch2.cells.append(migrating_cell)
                        return
            
            
            # Vein bifurcation option 2 - Implement form of capillary bifurcation rules at vein
            if (vein_bf_option == 2):
                
                # Bifurcation Rule 3 - Choose a branch at random, equal chance
                if (branch_rule == 3):
                    if (branch1.type == "capillary"):
                        if (rand.uniform(0, 1) <= 0.5):
                            migrating_cell.vessID = branch1.ID
                            branch1.cells.append(migrating_cell)
                            return
                        
                    if (branch2.type == "capillary"):
                        if (rand.uniform(0, 1) <= 0.5):
                            migrating_cell.vessID = branch2.ID
                            branch2.cells.append(migrating_cell)
                            return
                
                
                # Bifurcation Rule 7 - Weighted average of WSS and cell number
                if (branch_rule == 7):
                    if (branch1.type == "capillary"):
                            tau_ratio = branch1.tau/parent_vess.tau
                            cell_ratio = branch1.num_cells/parent_vess.num_cells
                            
                            prob1 = branch_alpha*tau_ratio + (1 - branch_alpha)*cell_ratio
                            
                            if (rand.uniform(0,1) <= prob1):
                                migrating_cell.vessID = branch1.ID
                                branch1.cells.append(migrating_cell)
                                return
                    
                    if (branch2.type == "capillary"):
                            tau_ratio = branch2.tau/parent_vess.tau
                            cell_ratio = branch2.num_cells/parent_vess.num_cells
                            
                            prob2 = branch_alpha*tau_ratio + (1 - branch_alpha)*cell_ratio
                            
                            if (rand.uniform(0,1) <= prob2):
                                migrating_cell.vessID = branch2.ID
                                branch2.cells.append(migrating_cell)
                                return                    
        
        
    # If all vessels in question are capillaries, enact bifurcation rule      
    if (branch1.type == "capillary" and branch2.type == "capillary"):
        
        # Bifurcation Rule 1 - Choose branch with the highest flow
        if (branch_rule == 1):
            if (abs(branch1.Q) > abs(branch2.Q)):
                migrating_cell.vessID = branch1.ID
                branch1.cells.append(migrating_cell)
                return
            else:
                migrating_cell.vessID = branch2.ID
                branch2.cells.append(migrating_cell)
                return
        
        
        # Bifurcation Rule 2 - Choose the branch with most similar polarity
        if (branch_rule == 2):
            if (migrating_cell.polarity*branch1.unit < migrating_cell.polarity*branch2.unit):
                migrating_cell.vessID = branch1.ID
                branch1.cells.append(migrating_cell)
                return
            else:
                migrating_cell.vessID = branch2.ID
                branch2.cells.append(migrating_cell)
                return
            
        
        # Bifurcation Rule 3 - Choose one branch at random, equal probability        
        if (branch_rule == 3):
            if (rand.uniform(0, 1) > 0.5):
                migrating_cell.vessID = branch1.ID
                branch1.cells.append(migrating_cell)
                return
            else:
                migrating_cell.vessID = branch2.ID
                branch2.cells.append(migrating_cell)
                return
        
        
        # Bifurcation Rule 4 - Choose one branch at random, equal probability        
        if (branch_rule == 4):
            if (rand.uniform(0, 1) > 0.3):
                migrating_cell.vessID = branch1.ID
                branch1.cells.append(migrating_cell)
                return
            else:
                migrating_cell.vessID = branch2.ID
                branch2.cells.append(migrating_cell)
                return
            
        
        # Bifurcation Rule 5 - Assign probability based on number of cells        
        if (branch_rule == 5):
            if (branch1.num_cells < parent_vess.num_cells) and (branch2.num_cells < parent_vess.num_cells):
                prob1 = (branch1.num_cells/parent_vess.num_cells)
            else:
                prob1 = (branch1.num_cells/(branch1.num_cells + branch2.num_cells))
            
            if (rand.uniform(0, 1) < prob1):
                migrating_cell.vessID = branch1.ID
                branch1.cells.append(migrating_cell)
                return
            else:
                migrating_cell.vessID = branch2.ID
                branch2.cells.append(migrating_cell)
                return
            
            
        # Bifurcation Rule 6 - Weighted average of flow and cell number
        if (branch_rule == 6):
            Q_ratio = branch2.Q/(branch1.Q + branch2.Q)
            cell_ratio = branch2.num_cells/(branch1.num_cells + branch2.num_cells)
            
            prob2 = branch_alpha*Q_ratio + (1 - branch_alpha)*cell_ratio
            
            if (rand.uniform(0, 1) > prob2):
                migrating_cell.vessID = branch1.ID
                branch1.cells.append(migrating_cell)
                return
            else:
                migrating_cell.vessID = branch2.ID
                branch2.cells.append(migrating_cell)
                return
            
            
        # Bifurcation Rule 7 - Weighted average of WSS and cell number
        if (branch_rule == 7):
            if ((branch1.tau + branch2.tau) != 0. and (branch1.num_cells + branch2.num_cells) != 0.):
            
                tau_ratio = branch2.tau/(branch1.tau + branch2.tau)
                cell_ratio = branch2.num_cells/(branch1.num_cells + branch2.num_cells)
                
                tau1 = branch1.tau/(branch1.tau + branch2.tau)
                cell1 = branch1.num_cells/(branch1.num_cells + branch2.num_cells)
                
                tau2 = branch2.tau/(branch1.tau + branch2.tau)
                cell2 = branch2.num_cells/(branch1.num_cells + branch2.num_cells)
                
                prob1 = branch_alpha*tau1 + (1 - branch_alpha)*cell1
                prob2 = branch_alpha*tau2 + (1 - branch_alpha)*cell2
                
                if (branch1.first_bprob == True):
                    branch1.bprob = prob1
                    branch1.first_bprob = False
                else:
                    prob1 = (bp_u*prob1 + bp_v*branch1.bprob)/(bp_u + bp_v)
                    branch1.bprob = prob1
                    
                if (branch2.first_bprob == True):
                    branch2.bprob = prob2
                    branch2.first_bprob = False
                else:
                    prob2 = (bp_u*prob2 + bp_v*branch2.bprob)/(bp_u + bp_v)
                    branch2.bprob = prob2
                    
                #assert (round(prob1, 8) == round((1 - prob2), 8)), "p1 = {0} does not match {1}".format(prob1, (1 - prob2))
                
                if (rand.uniform(0, 1) > prob2):
                    migrating_cell.vessID = branch1.ID
                    branch1.cells.append(migrating_cell)
                    return
                else:
                    migrating_cell.vessID = branch2.ID
                    branch2.cells.append(migrating_cell)
                    return
                
                    
