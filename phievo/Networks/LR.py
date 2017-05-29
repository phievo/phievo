"""
Definition of Ligand-Receptor Interaction
Creation: unknown
Last edition: 2016-10-26
"""
from phievo import __silent__,__verbose__
if __verbose__:
    print("Execute LR (Interaction Template)")

from . import classes_eds2
from . import deriv2
from . import mutation
import copy

#default range
mutation.dictionary_ranges['LR.association'] = 0.0 / mutation.T
mutation.dictionary_ranges['LR.threshold'] = 0.0 * mutation.C

#######################################
### LigandReceptor Class definition ###
#######################################

class LR(classes_eds2.Interaction):
    """Ligand-Receptor interaction between two species

Note that LR Interaction work on an enzymatic frame, (there is a
    threshold but no dissociation rate). For as assoc./dissoc. frame,
    use PPI instead.

    Attributes:
        association (float): the association rate
        threshold (float): the Michaelis Menten constant
        label (str): 'LR Interaction' by default
        input (list): list of input types: ['Ligand','Receptor']
        output (list): list of output types: ['Species']
    """
    def __init__(self,association=0,threshold=0):
        """Constructor of a new Ligand-Receptor

        Args:
            association (float): the association rate
            threshold (float): the Michaelis Menten constant
        """
        classes_eds2.Node.__init__(self)
        self.association=association
        self.threshold=threshold
        self.label='LR Interaction'
        self.input=['Ligand','Receptor']
        self.output=['Species']

    def __str__(self):
        return "{0.id} LR interaction: assoc. = {0.association:.2f}, thres. = {0.threshold:.2f}".format(self)

    def outputs_to_delete(self,net):
        """Returns the Ligand-Receptor complex to delete when removing the LR"""
        return net.graph.successors(self)

########## Attributes attached to Network for  LR interaction ##########

def number_LR(self):
        """Return the number of possible LR interactions"""
        nL=self.number_nodes('Ligand')
        nR=self.number_nodes('Receptor')
        n_LR=self.number_nodes('LR')
        return nL*nR-n_LR

def new_LR(self,ligand,receptor,association,threshold,types):
        """Create a new LR, its associated complex and add then to the network.

        Args:
            ligand (Species): ligand species
            receptor (Species): receptor species
            association (float): -
            threshold (float): -
            types (list): the types of the complex species

        Return:
            list: of the form [LR interaction,complex created]
            or None if an error occured
        """
        complex=classes_eds2.Species(types)
        lr_inter=LR(association,threshold)
        if lr_inter.check_grammar([ligand,receptor], [complex]):
            self.add_Node(complex)
            self.add_Node(lr_inter)
            self.graph.add_edge(ligand,lr_inter)
            self.graph.add_edge(receptor,lr_inter)
            self.graph.add_edge(lr_inter,complex)
            return [lr_inter,complex]
        else:
            print("Error in grammar, new_LR")
            return None

# Add the corepromoter functions to the Network class
setattr(classes_eds2.Network,'number_LR',number_LR)
setattr(classes_eds2.Network,'new_LR',new_LR)

########## Attributes attached to Mutable_Network for  LR interaction ##########

def new_random_LR(self, ligand, receptor):
    """Creates a LR with random parameters between the Species

    Args:
        ligand (Species): ligand species
        receptor (Species): receptor species

    Return:
        list: of the form [lr interaction,complex created]
    """
    types=[ ['Degradable', mutation.sample_dictionary_ranges('Species.degradation',self.Random)],
            ['TF',int(self.Random.random()*2)],['Complexable'],['Phosphorylable'] ]
    association = mutation.sample_dictionary_ranges('LR.association',self.Random)
    threshold = mutation.sample_dictionary_ranges('LR.threshold',self.Random)
    return self.new_LR(ligand,receptor,association,threshold,types)

def random_LR(self):
    """Create new random LR among all those possible


    Return:
        list: of the form [lr interaction,complex created]
        or None if an error occured
    """
    if 'Ligand' in self.list_types and 'Receptor' in self.list_types:
        #Evaluate all possible LR interactions
        possible_LR = [(lig,rec) for lig in self.list_types['Ligand']
                                 for rec in self.list_types['Receptor']
                                 if not self.check_existing_binary([lig,rec],'LR Interaction')]
        n_pLR=len(possible_LR)
        if not (n_pLR==self.number_LR()):
            print("Potential Bug : Inconsistency in Computation of number of LR")
        #Draw a random couple and create the interaction
        if (n_pLR==0):
            print("In random_LR : No other posible LRs")
            return None
        else:
            [L,R] = possible_LR[int(self.Random.random()*n_pLR)]
            return self.new_random_LR(L, R)
    else:
        print("Error in random_LR (try to create a LR from non exsiting pieces)")
        return None

# Add the corepromoter functions to the Mutable_Network class
setattr(mutation.Mutable_Network,'random_LR',random_LR)
setattr(mutation.Mutable_Network,'new_random_LR',new_random_LR)

########## Integration C Tools ##########

def compute_LR(net):
    """gives the string corresponding to LR for integration

    Return:
        str: a single string for all LR in the network
    """
    func="\n/**************LR interactions*****************/\n"
    func+="void LRinC(double s[],double ds[],double ligands[]){\n"
    func = func + "    double increment=0;\n    double rate=0;\n"
    if ('LR' in net.list_types):
        for index in net.list_types['LR']:
            C=net.graph.successors(index)[0]#finds the product of LR interaction
            [P1,P2]=net.graph.predecessors(index) #find the components
            L,R = (P1,P2) if P1.isinstance('Ligand') else (P2,P1) #determine the ligand and the receptor
            arate="%f*"%index.association+"ligand"+L.id+"*"+R.id+"/("+R.id+"+"+str(index.threshold)+")"
            func=func+deriv2.compute_leap([R.id],[C.id],arate)
    func=func+"}\n \n"

    return func

#update deriv2
deriv2.interactions_deriv_inC["LR"] = compute_LR
