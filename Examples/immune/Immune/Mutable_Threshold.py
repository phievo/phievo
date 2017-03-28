
# importing all relevant modules.
import phievo.Networks.classes_eds2 as classes_eds2
import phievo.Networks.mutation as mutation
import copy
import phievo.Networks.deriv2 as deriv2
import random

# overwritting the Phospho tag from classes_eds2 to keep track of the degree of phosphorylation.
# Tags_Species['Phospho'] = ['n_Phosho']
# does not seem to work.

# parameters needed to perform evolution.
mutation.dictionary_ranges['Mutable_Threshold.concentration']=mutation.C

# definition of class phosphorylation: sub-class of interaction class (see classes_eds2).
# contains the parameter quantifying the interaction as well as the input and output of the reaction. 
# these characteristics define the interaction completely. they are "free standing" initially in that 
# they are not appended to the graph per se. see new_Phosphorylation below to see how to actually 
# add the interaction to the network.

class Mutable_Threshold(classes_eds2.Interaction):
    """ This "interaction" only contains the threshold for activation. """

    def __init__(self,thresh=0):
        classes_eds2.Node.__init__(self)
        self.thresh = thresh
        self.label='Mutable_Threshold'
        self.input=[]
        self.output=[]
        self.removable = False # by default the mutable threshold node cannot be removed.


    def outputs_to_delete(self,net):
        """ Returns the species to delete when deleting Mutable_Threshold."""
        return None



####### Attributes attached to Network for  Phosphorylations/Dephosphorylations############



# important function that adds the interaction to the actual graph. 

def new_Mutable_Threshold(self,thresh):

        Mut_Thresh=Mutable_Threshold(thresh) #creates the interaction (see function above)  
        self.add_Node(Mut_Thresh) # adding the various nodes and edges corresponding to the interaction.            
        return Mut_Thresh


setattr(classes_eds2.Network,'new_Mutable_Threshold',new_Mutable_Threshold)
