## Create a methylation interaction
from phievo import __silent__,__verbose__
if __verbose__:
    print("Execute Methyl (Interaction Template)")

from phievo.Networks import mutation
from phievo.Networks import deriv2
from phievo.Networks import classes_eds2
import copy


## Define the function that assigns new parameters to a new methylable species 
mutation.species_types["Methylable"] = lambda random_generator:[
    ["Methylable"],
    ['Diffusible',mutation.sample_dictionary_ranges('Species.diffusion',random_generator)]
]
classes_eds2.Species.Tags_Species["Methylable"] = []


## Define the default dictionary_range
mutation.dictionary_ranges['Methyl.methyl'] = 0.0/(mutation.C*mutation.T)
mutation.dictionary_ranges['Methyl.demethyl'] = 0.0/mutation.T

class Methyl(classes_eds2.Interaction):
    """
    Methylation interaction

    Args:
        Methyl.methyl(float): binding rate of a methyl group
        Methyl.demethyl(float): unbinding rate of a methyl group
        label(str): Methylation
        input(list): Type of the inputs
        output(list): Type of the outputs
    """
    def __init__(self,methyl=0,demethyl=0):
        classes_eds2.Node.__init__(self)
        self.methyl=methyl
        self.demethyl=demethyl
        self.label='Methylation'
        self.input=['Methylable']
        self.output=['Species']

    def __str__(self):
        """
        Used by the print function to display the interaction.
        """
        return "{0.id} Methylation: methyl. = {0.methyl:.2f}, demethyl = {0.demethyl:.2f}".format(self)

    def outputs_to_delete(self,net):
        """
        Returns the methylated form of the species to delete when the reaction is deleted.
        """
        return net.graph.successors(self)
    
#################################################
#### Functions to add to the Mutable_Network ####
#################################################

def number_Methyl(self):
    """
    Returns the number of possible methylation in the current network.
    Note: this function is optional, it is used to check the consistency of
    the random_Methyl function.
    """
    n = self.number_nodes('Methylable')
    n_Methyl = self.number_nodes('Methyl')
    return n-n_Methyl

def new_Methyl(self,S,methyl,demethyl,parameters):
    """
    Creates a new :class:`Networks.Methyl.Methyl` and the species methylated for in the the network.

    Args:
        S: species to methylate
        methyl(float): binding rate of a methyl group
        demethyl(float): unbinding rate of a methyl group
        parameters(list): Parameters of the methylated species
    Returns:
        [methyl_inter,S_methyl]: returns a Methyl interaction and a methylated species.
    """

    S_methyl = classes_eds2.Species(parameters)
    meth_inter = Methyl(methyl,demethyl)
    assert meth_inter.check_grammar([S],[S_methyl]),"Error in grammar, new Methylation"
    self.add_Node(S_methyl)
    self.add_Node(meth_inter)
    self.graph.add_edge(S,meth_inter)
    self.graph.add_edge(meth_inter,S_methyl)
    return [meth_inter,S_methyl]


def new_random_Methyl(self,S):
    """
    Creates a methylation with random parameters.
        
    Args:
        S: Species to methylate
    Returns:
        [methyl_inter,S_methyl]:returns a Methyl interaction and a methylated species.
    """
    parameters = {}
    if S.isinstance("TF"):
        parameters['TF'] = self.Random.random()*2
    for tt in S.types:
        if tt not in ["TF","Methylable","Input","Output"]:
            parameters[tt] = [mutation.sample_dictionary_ranges('Species.{}'.format(attr),self.Random) for attr in S.Tags_Species[tt]]

    # Transform to fit phievo list structure
    parameters = [[kk]+val if val else [kk] for kk,val in parameters.items()]
    methyl = mutation.sample_dictionary_ranges('Methyl.methyl',self.Random)
    demethyl = mutation.sample_dictionary_ranges('Methyl.demethyl',self.Random)
    return self.new_Methyl(S,methyl,demethyl,parameters)
    
def random_Methyl(self):
    """
    Evaluates the species that can be phosphorilated, picks one an create a random
    methylation. The random mutation is made using :func:`new_random_Methyl <phievo.Networks.classes_eds2.new_random_Methyl>`.

    Returns:
        [methyl_inter,S_methyl]: returns a Methyl interaction and a methylated species.
    """
    try:
        list_methylable=self.dict_types["Methylable"]
    except KeyError:
        print("\tThe network contain no Methylacble species.")
        raise
    list_possible_methylable = []
    for S in list_methylable:
        if not self.check_existing_binary([S],"Methyl"):
            list_possible_methylable.append(S)
    n_possible = len(list_possible_methylable)
    assert n_possible==self.number_Methyl(),"The number of possible new methylation does not match its theoretical value."
    if n_possible==0:
        if __verbose__:
            print("No more possible methylation.")
        return None
    else:
        S = list_possible_methylable[int(self.Random.random()*n_possible)]
        return self.new_random_Methyl(S)
        
def Methyl_deriv_inC(net):
    """
    Function called to generate the string corresponding to in a methylation in C.
    """
    func_str = "\n/************** Methylations *****************/\n"
    methylations = net.dict_types.get("Methyl",[])
    for methyl_inter in methylations:
        S = net.graph.predecessors(methyl_inter)[0]
        S_meth = net.graph.successors(methyl_inter)[0]
        f_rate = "{M.methyl}*{S.id}".format(M=methyl_inter,S=S)
        b_rate = "{M.demethyl}*{S_m.id}".format(M=methyl_inter,S_m=S_meth)

        func_str += deriv2.compute_leap([S.id],[S_meth.id],f_rate)
        func_str += deriv2.compute_leap([S_meth.id],[S.id],b_rate)
    return func_str

## Add the current the new functions to the network.
setattr(classes_eds2.Network,"number_Methyl",number_Methyl)
setattr(classes_eds2.Network,"new_Methyl",new_Methyl)
setattr(classes_eds2.Network,"new_random_Methyl",new_random_Methyl)
setattr(classes_eds2.Network,"random_Methyl",random_Methyl)
deriv2.interactions_deriv_inC["Methyl"] = Methyl_deriv_inC
