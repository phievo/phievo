"""
Definition of Protein-Protein-Interaction
Creation: unknown
Last edition: 2016-10-26
"""
from phievo import __silent__,__verbose__
if __verbose__:
    print("Execute PPI (Interaction Template)")

from . import classes_eds2
from . import mutation
from . import deriv2
import copy

#default range
mutation.dictionary_ranges['PPI.association'] = 0.0/(mutation.C*mutation.T)
mutation.dictionary_ranges['PPI.disassociation'] = 0.0/mutation.T

############################
### PPI Class definition ###
############################

class PPI(classes_eds2.Interaction):
    """Protein-protein interaction between two species

    Args:
        association (float): the association rate
        disassociation (float):  the dissociation rate fo the complex
        label (str): 'PP Interaction' by default
        input (list): list of input types: ['Complexable','Complexable']
        output (list): list of output types: ['Species']
    """
    def __init__(self,association=0,disassociation=0):
        classes_eds2.Node.__init__(self)
        self.association=association
        self.disassociation=disassociation
        self.label='PP Interaction'
        self.input=['Complexable','Complexable']
        self.output=['Species']

    def __str__(self):
        return "{0.id} PPI: assoc. = {0.association:.2f}, dissoc. = {0.disassociation:.2f}".format(self)

    def outputs_to_delete(self,net):
        """Return the complex to delete when removing the LR"""
        return net.graph.list_successors(self)
    
    def check_grammar(self,input_list,output_list):
        """checks the grammar for the interactions (custom for PPI)

        Args:
            input_list (list): nodes to be checked
            output_list (list): nodes to be checked

        Return:
            Boolean for the consistency of up and downstream grammar
        """
        if len(input_list) == 1: #auto-phosphorylation
            input_check = classes_eds2.check_consistency(['Complexable'],input_list)
            output_check = classes_eds2.check_consistency(self.output,output_list)
            return input_check and output_check
        else:
            return classes_eds2.Interaction.check_grammar(self,input_list,output_list)

########## Attributes attached to Network for PPI ##########

def number_PPI(self):
    """Return the number of possible PPI in network"""
    n=self.number_nodes('Complexable')
    n_PPI=self.number_nodes('PPI')
    return n*(n+1)//2 - n_PPI #n(n-1)/2+n - self interactions

def new_PPI(self,P1,P2,assoc,dissoc,types):
    """Create a new :class:`Networks.PPI.PPI`, its associated complex and add then to the network.

    Args:
        P1 (:class:`Species <phievo.Networks.classes_eds2.Species>`): First Protein
        P2 (:class:`Species <phievo.Networks.classes_eds2.Species>`): Second Protein
        assoc (float): the association rate
        dissoc (float): the dissociation rate of the complex
        types (list): the types of the complex species
    Returns:
        list of the form [`ppi`,`complex created`] with:
            - `ppi`: :class:`PPI <phievo.Networks.PPI.PPI>`
            - `complex created`: :class:`Species <phievo.Networks.classes_eds2.Species>`
    """
    complex = classes_eds2.Species(types)
    ppi=PPI(assoc,dissoc)
    if ppi.check_grammar([P1,P2], [complex]):
        self.add_Node(complex)
        self.add_Node(ppi)
        self.graph.add_edge(P1,ppi)
        self.graph.add_edge(P2,ppi)
        self.graph.add_edge(ppi,complex)
        return [ppi,complex]
    else:
        print("Error in grammar, new_Complex")
        return None

def duplicate_PPI(self,species,D_species,interaction,module,D_module):
    """function to duplicate a PPI interaction

    Args:
        species (:class:`Species <phievo.Networks.classes_eds2.Species>`): the original species
        D_species (:class:`Species <phievo.Networks.classes_eds2.Species>`): the new species
        interaction (:class:`PPI <phievo.Networks.PPI.PPI>`): the interaction you want to duplicate
        module (:class:`TModule <phievo.Networks.classes_eds2.TModule>`): the original module
        D_module (:class:`TModule <phievo.Networks.classes_eds2.TModule>`): the new module

    """
    # Copy the interaction and add it to the network
    D_interaction=copy.deepcopy(interaction)
    D_interaction.mutable=1
    D_interaction.removable=True
    self.add_Node(D_interaction)

    # Copy the complex and add it to the network
    Complex=self.graph.list_successors(interaction)[0]
    D_Complex=copy.deepcopy(Complex)
    D_Complex.mutable=1
    D_Complex.removable=True
    self.add_Node(D_Complex)

    self.graph.add_edge(D_interaction,D_Complex)#copies the PPI
    PPI_components=self.graph.list_predecessors(interaction)
    PPI_components.remove(species) # removes the duplicated component
    component=PPI_components[0] #other PPI component
    self.graph.add_edge(component,D_interaction)
    self.graph.add_edge(D_species,D_interaction)
    self.duplicate_downstream_interactions(Complex,D_Complex,module,D_module)
    if (component==species): #if species complexes with itself we must create a new dimer
        D_Complex_2=copy.deepcopy(Complex)#copies the complex
        D_Complex_2.mutable=1
        D_Complex_2.removable=True
        D_interaction_2=copy.deepcopy(interaction)
        D_interaction_2.mutable=1
        D_interaction_2.removable=True
        self.add_Node(D_interaction_2)
        self.graph.add_edge(D_interaction_2,D_Complex_2)#copies the PPI
        self.graph.add_edge(D_species,D_interaction_2)
        self.graph.add_edge(D_species,D_interaction_2)
        self.duplicate_downstream_interactions(Complex,D_Complex_2,module,D_module)

# Add the corepromoter functions to the Network class
setattr(classes_eds2.Network,'number_PPI',number_PPI)
setattr(classes_eds2.Network,'new_PPI',new_PPI)
setattr(classes_eds2.Network,'duplicate_PPI',duplicate_PPI)

########## Attributes attached to Mutable_Network for PPI ##########

def new_random_PPI(self, P1, P2):
    """Creates a PPI with random parameters between the Species

    Args:
        P1 (:class:`Species <phievo.Networks.classes_eds2.Species>`): First  protein
        P2 (:class:`Species <phievo.Networks.classes_eds2.Species>`): Second  protein

    Returns:
        list of the form [`ppi`,`complex created`] with:
            - `ppi`: :class:`PPI <phievo.Networks.PPi.PPI>`
            - `complex created`: :class:`Species <phievo.Networks.classes_eds2.Species>`
    """
    types=[['Degradable', mutation.sample_dictionary_ranges('Species.degradation',self.Random)],
          ['Complex'],['Phosphorylable']]
    #Then we check the types of P1 and P2 and ensure that the complex has them too
    if (P1.isinstance('TF') or P2.isinstance('TF')):
        activity=int(self.Random.random()*2)
        types.append(['TF',activity])
    if (P1.isinstance('Receptor') and P2.isinstance('Receptor')):
        types.append(['Receptor'])
    if (P1.isinstance('Kinase') or P2.isinstance('Kinase')):
        types.append(['Kinase'])
    assoc = mutation.sample_dictionary_ranges('PPI.association',self.Random)
    dissoc = mutation.sample_dictionary_ranges('PPI.disassociation',self.Random)
    return self.new_PPI(P1,P2,assoc,dissoc,types)

def random_PPI(self):
    """Create new random PPI among all those possible

    Returns:
        list of the form [`ppi`,`complex created`] with:
            - `ppi`: :class:`PPI <phievo.Networks.PPI.PPI>`
            - `complex created`: :class:`Species <phievo.Networks.classes_eds2.Species>`
    """
    if 'Complexable' in self.dict_types:
        list_complexable=self.dict_types['Complexable']
        n=len(list_complexable)
        list_possible_PPI = []
        for ip1 in range(n):
            P1=list_complexable[ip1]
            if not self.check_existing_binary([P1],'PPI') and not self.check_existing_binary([P1,P1],'PPI'):
                list_possible_PPI.append([P1,P1])
            for ip2 in range(ip1):
                    P2=list_complexable[ip2]
                    if not self.check_existing_binary([P1,P2],'PPI'):
                        list_possible_PPI.append([P1,P2])
        n_pPPI=len(list_possible_PPI)
        if not (n_pPPI==self.number_PPI()):
            print("Potential Bug : Inconsistency in Computation of number of PPI")
        if (n_pPPI==0):
            print("In random_PPI : No other posible PPIs")
            return None
        else:
            [P1,P2]=list_possible_PPI[int(self.Random.random()*n_pPPI)]
            [PPI,C]=self.new_random_PPI(P1,P2)
            return [PPI,C]
    else:
        print("Error in random_PPI (try to create a PPI from non existing pieces)")
        return None

# Add the corepromoter functions to the Mutable_Network class
setattr(mutation.Mutable_Network,'random_PPI',random_PPI)
setattr(mutation.Mutable_Network,'new_random_PPI',new_random_PPI)

########## Integration C Tools ##########

def PPI_deriv_inC(net):
    """gives the string corresponding to :class:`Networks.PPI.PPI` for integration

    Return:
        str a single string for all :class:`Networks.PPI.PPI` in the network
    """
    func="\n/**************Protein protein interactions*****************/\n"
    if ('PPI' in net.dict_types):
        for index in net.dict_types['PPI']:
            C=net.graph.list_successors(index)[0]#finds the complex
            #print net.graph.list_predecessors(index)
            list_Pi=net.graph.list_predecessors(index) #find the components
            P1=list_Pi[0]
            #we need a special case in the new version of networkx for self complexation
            if (len(list_Pi)==1):
                P2=P1
            else:
                P2=list_Pi[1]
            arate="%f * %s * %s"%(index.association,P1.id,P2.id)
            drate="%f * %s"%(index.disassociation,C.id)
            func=func+deriv2.compute_leap([P1.id,P2.id],[C.id],arate)
            func=func+deriv2.compute_leap([C.id],[P1.id,P2.id],drate)
    return func

#update deriv2
deriv2.interactions_deriv_inC["PPI"] = PPI_deriv_inC
