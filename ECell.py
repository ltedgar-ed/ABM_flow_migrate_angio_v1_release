#------------------------------------------------------------------------------
# The ECell class represents a migratory endothelial cell within the blood vessel wall
#
# Lowell Taylor Edgar
# University of Edinburgh
# 2019

import random as rand
from Tensor3D import *

#------------------------------------------------------------------------------
class ECell:
    # ECell constructor
    def __init__(self, ID =0, vessID =0, polvect =Vect(0., 0., 0.,)):
        
        # Cell number identifier
        ID : int 
        self.ID = ID
        
        # Identifier for which Vessel the cell resides in
        vessID : int
        self.vessID = vessID
        
        # Cell polarity vector
        polarity : Vect
        polvect = Vect(rand.uniform(-1, 1), rand.uniform(-1, 1), 0.)
        polvect.unit()
        self.polarity = polvect
        
        # Migration idicator (+1 for upstream/against vessel unit vector, -1 for downstream/along vessel unit vector, 0 for not migrating)
        migrate : int
        self.migrate = 0
        
        
    def __str__(self):
        return "Cell {0} is in Vessel {1}; polarity = {2}; migrate = {3}".format(self.ID, self.vessID, round(self.polarity, 2), self.migrate)

    
 #------------------------------------------------------------------------------       
# main program
if __name__ == "__main__":
    pass
