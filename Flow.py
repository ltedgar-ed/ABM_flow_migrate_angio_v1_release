#------------------------------------------------------------------------------
# The Flow class contains the methods needed to solve for flow within a generalized
# vascular network consisting of Vessels
#
# Lowell Taylor Edgar
# University of Edinburgh
# 2019

from Vessel import *
import numpy

#------------------------------------------------------------------------------
def solve_for_flow(nodes, vessels, node_degree, BCs, bifur_nodes, bifur):
    # List of Nodes
    nodes : list[tuple]
    # List of Vessels
    vessels : list[Vessel]
    # List of boundary conditions (known pressures at inlets and outlets)
    BCs :  list[tuple]
    
    n = len(nodes)                                                              # Number of nodes       
    v = len(vessels)                                                            # Number of Vessels
    nBCs = len(BCs)                                                             # Number of boundary conditions (known pressures)
    
    p = numpy.zeros(n)                                                          # Nodal pressure array
    A = numpy.zeros((n,n))                                                      # Conductance matrix
    b = numpy.zeros(n)                                                          # Solution array
    

    # Construct the conductance matrix and solution array
    for i in range(n):
          for vess in vessels:
    
              if (vess.n1 == i):
                   if (vess.n0 in BCs):
                       b[i] -= vess.G*p[vess.n0]
                       A[i][i] -= vess.G
                   else:
                       A[i][i] -= vess.G
                       A[i][vess.n0] += vess.G
                     
              if (vess.n0 == i):
                   if (vess.n1 in BCs):
                       b[i] -= vess.G*p[vess.n1]
                       A[i][i] -= vess.G
                   else:
                       A[i][i] -= vess.G
                       A[i][vess.n1] += vess.G
    
    if (inlet_flow_on == False):
        # Apply the pressure boundary conditions
        for BC in BCs:
            A[BC[0]][BC[0]] = 1
            b[BC[0]] = BC[1]
            p[BC[0]] = BC[1]    
    
    if ((sim_type == "A branch") and (inlet_flow_on == True)):
        # A branch inlet flow BC    
        b[0] = -1.1057668948525408e-13
        A[BCs[1][0]][BCs[1][0]] = 1
        b[BCs[1][0]] = BCs[1][1]
        p[BCs[1][0]] = BCs[1][1]
    
    if ((sim_type == "Y branch") and (inlet_flow_on == True)):
        # Y branch inlet flow BCs
        qinR = -6.143149415847469e-14
        qinL = qinR
        #qinL = qinR/10
        #qinL = qinR/20
        b[20] = qinR
        b[30] = qinL
        A[BCs[0][0]][BCs[0][0]] = 1
        b[BCs[0][0]] = BCs[0][1]
        p[BCs[0][0]] = BCs[0][1]
    
    # Solve the system of equations for unknown pressures
    p = numpy.linalg.solve(A, b)
    #A_pinv = numpy.linalg.pinv(A)
    
    #p = numpy.dot(A_pinv, b)
    
    # Calculate flow within each Vessel
    for vess in vessels:
        vess.P0 = p[vess.n0]
        vess.P1 = p[vess.n1]
        
        vess.calc_flow()
    
    
    # Determine if bifurcations are flow-converging or not
         
    # Find all vessels involved with each bifurcation
    for i in range(len(bifur_nodes)):
        bifur[i] = 0
        bif_node = bifur_nodes[i]
        vess_n0 = []
        vess_n1 = []
        
        for vess in vessels:
            if (vess.n0 == bif_node):
                vess_n0.append(vess.ID)
                
            if (vess.n1 == bif_node):
                vess_n1.append(vess.ID)
        
        # Determine if the bifurcation if flow-converging or flow-diverging
        if (len(vess_n0) == 2 and len(vess_n1) == 1):
            v01 = vess_n0[0]-1
            v02 = vess_n0[1]-1
            v11 = vess_n1[0]-1
            
            if (vessels[v01].Q > 0 and vessels[v02].Q > 0 and vessels[v11].Q > 0):
                bifur[i] = -1
                #print("Flow at node {} is diverging".format(bif_node))
                
            if ((vessels[v01].Q > 0 and vessels[v02].Q < 0) or (vessels[v01].Q < 0 and vessels[v02].Q > 0)) and vessels[v11].Q > 0:
                bifur[i] = 1
                #print("Flow at node {} is converging".format(bif_node))
            
            if (vessels[v01].Q < 0 and vessels[v02].Q < 0 and vessels[v11].Q < 0):
                bifur[i] = 1
                #print("Flow at node {} is converging".format(bif_node))
            
            if ((vessels[v01].Q > 0 and vessels[v02].Q < 0) or (vessels[v01].Q < 0 and vessels[v02].Q > 0)) and vessels[v11].Q < 0:
                bifur[i] = -1
                #print("Flow at node {} is diverging".format(bif_node))
                
        
        if (len(vess_n0) == 1 and len(vess_n1) == 2):
            v01 = vess_n0[0]-1
            v11 = vess_n1[0]-1
            v12 = vess_n1[1]-1
            
            if (vessels[v01].Q > 0 and vessels[v11].Q > 0 and vessels[v12].Q > 0):
                bifur[i] = 1
                #print("Flow at node {} is converging".format(bif_node))
                
            if (vessels[v01].Q < 0 and vessels[v11].Q < 0 and vessels[v12].Q < 0):
                bifur[i] = -1
                #print("Flow at node {} is diverging".format(bif_node))
                
            if vessels[v01].Q > 0 and ((vessels[v11].Q < 0 and vessels[v12].Q > 0) or (vessels[v11].Q > 0 and vessels[v12].Q < 0)):
                bifur[i] = -1
                #print("Flow at node {} is diverging".format(bif_node))
                
            if vessels[v01].Q < 0 and ((vessels[v11].Q < 0 and vessels[v12].Q > 0) or (vessels[v11].Q > 0 and vessels[v12].Q < 0)):
                bifur[i] = 1
                #print("Flow at node {} is converging".format(bif_node))