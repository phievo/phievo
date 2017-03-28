
# importing all relevant modules.
import phievo.Networks.classes_eds2 as classes_eds2
import phievo.Networks.mutation as mutation
import copy        
from phievo.initialization_code import *
import phievo.Networks.deriv2 as deriv2

# parameters needed to perform evolution.
mutation.dictionary_ranges['Initial_Concentration.concentration']=mutation.C


class Initial_Concentration(classes_eds2.Interaction):
    """ Fixes concentration of kinase or phosphatase of a given type, conc : concentration of the species of interest. """

    # This is almost a trivial class, which essentially only holds a concentration.
    def __init__(self,c=0):
        classes_eds2.Node.__init__(self)
        self.conc=c
        self.label='Initial_Concentration'
        self.input=['Species']
        self.output=['Species']

    def __str__(self):
        return "{0.id} Initial concentration: Conc. = {0.conc:.2f}".format(self)

    def outputs_to_delete(self,net):
        """ Returns the species to delete when deleting a Initial_Concentration"""
        #return net.Init_Conc_to_remove(self)
        return net.graph.predecessors(self)   
        
#def Init_Conc_to_remove(self,interaction):
#        """Returns the species whose initial concentration is specified by Initial_Concentration."""
#        listOut=self.graph.predecessors(interaction)        
#        return listOut

# important function that adds the interaction to the actual graph. 
def new_Initial_Concentration(self,species,conc):

        Init_Conc=Initial_Concentration(conc) #creates the interaction (see function above, an instance of the Node classe)
        if Init_Conc.check_grammar([species],[species]):
            self.add_Node(Init_Conc)
            self.graph.add_edge(species,Init_Conc)
            self.graph.add_edge(Init_Conc,species)

setattr(classes_eds2.Network,'new_Initial_Concentration',new_Initial_Concentration)


def new_Species(self,list,bool_conc = True,manual_conc = 0):
        """Create new Species instance with argument parameters of form [ [], [],.. ], and add to graph """
        S=classes_eds2.Species(list)
        self.add_Node(S)
        if ((S.isinstance('Phosphatase') or S.isinstance('Kinase')) and not S.isinstance('pMHC') ):
            if bool_conc:
                conc = mutation.sample_dictionary_ranges('Initial_Concentration.concentration',self.Random)
            else:
                conc = manual_conc
            self.new_Initial_Concentration(S,conc)
        return S
        
setattr(classes_eds2.Network,'new_Species',new_Species)


