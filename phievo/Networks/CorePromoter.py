"""
Definition of CorePromoter Interaction
Creation: unknown
Last edition: 2016-10-25
"""
from phievo import __silent__,__verbose__
if __verbose__:
    print("Execute CorePromoter (Interaction Template)")

from . import classes_eds2
from . import mutation
import copy

#default range
mutation.dictionary_ranges['CorePromoter.delay'] = 0   # convert to int(T/dt) in run_evolution.py

#####################################
### CorePromoter Class definition ###
#####################################

class CorePromoter(classes_eds2.Interaction):
    """Core promoter for transcription of one species

    The CorePromoter serve as a link between the TModule and the Species
    to preserve the bipartite nature of the network.

    Attributes:
        delay (int):
        label (str): 'transcription' by default
        input (list): list of input types: ['TModule']
        output (list): list of output types: ['Species']
    """
    def __init__(self, delay=0):
        """Constructor of a new CorePromoter

        Args:
            delay (int): -
        """
        classes_eds2.Node.__init__(self)
        self.delay=int(delay)
        self.label='transcription'
        self.input=['TModule']
        self.output=['Species']

    def __str__(self):
        return "{0.id} CorePromoter: delay = {0.delay}".format(self)

    def outputs_to_delete(self,net):
        """indicate the Nodes to remove when deleting the CorePromoter

        Args:
            net (Networks): The network to which the CP belongs

        Return:
            list: All the predec. and succ. of self in net
        """
        toDelete = net.graph.predecessors(self)
        successorNode = net.graph.successors(self)[0]
        ## Search corepromoters in successorNode predecessors
        types = [type(inter) is type(self) for inter in net.graph.predecessors(successorNode)]
        if sum(types)==1:
            ## Delete successorNode if self is its last corepromoter
            toDelete.append(successorNode)
        return toDelete

    def string_param(self):
        """Self description of the Interaction"""
        return "Delay=%i"%self.delay

########## Attributes attached to Network for CorePromoter ##########

def add_CorePromoter2Species(self,inter,output):
    """Add a CorePromoter Interaction and its output to the network

    Args:
        inter (:class:`Networks.CorePromoter.CorePromoter`): the CorePromoter to be added
        output (:class:`Networks.classes_eds2.Species`): the CorePromoter output

    Return:
        None: in place modification
    """
    if inter.isinstance('CorePromoter'):
        self.add_Node(inter)
        self.add_Node(output)
        self.graph.add_edge(inter,output)
    else:
        print("Error in grammar in add_CorePromoter2Species")

def add_TModule2CorePromoter(self,module,inter):
    """Add a CorePromoter Interaction and its TModule to the network

    Args:
        module (:class:`Networks.classes_eds2.TModule`): the CorePromoter module
        inter (:class:`Networks.CorePromoter.CorePromoter`): the CorePromoter to be added

    Return:
        None: in place modification
    """
    if module.isinstance('TModule'):
        self.add_Node(inter)
        self.add_Node(module)
        self.graph.add_edge(module,inter)
    else:
        print("Error in grammar in add_TModule2CorePromoter")

def new_gene(self, rate, delay, parameters,basal=0.):
    """Create a complete new gene (TModule, CorePromoter and Species)

    Args:
        rate (float): the rate production of the TModule
        delay (int): the delay of the CorePromoter
        parameters (list): the species parameter (see Network.new_Species)
        basal (float): the basal production of the TModule (default to 0.)

    Return:
        list: of the form [:class:`Networks.classes_eds2.TModule`, :class:`Networks.CorePromoter.CorePromoter`, :class:`Networks.classes_eds2.Species`]
        or None if an error occured
    """
    species = self.new_Species(parameters)
    module = classes_eds2.TModule(rate,basal)
    prom = CorePromoter(delay)
    if prom.check_grammar([module],[species]):
        self.add_CorePromoter2Species(prom,species)
        self.add_TModule2CorePromoter(module,prom)
        return [module, prom, species]
    else:
        print("Error in grammar, new_gene")
        return None

def duplicate_gene(self,species):
    """Duplicate a gene, i.e. a triplet Tmodule/CorePromoter/Species

    Args:
        species (:class:`Networks.classes_eds2.Species`): Species to duplicate

    Return:
        list: of the form [`new_TModule`, `new_CorePromoter`, `new_Species`, `old_TModule`]
            - `new_TModule`: :class:`Networks.classes_eds2.TModule`
            - `new_CorePromoter`: :class:`Networks.CorePromoter.CorePromoter`
            - `new_Species`: :class:`Networks.classes_eds2.Species`
            - `old_TModule`: :class:`Networks.classes_eds2.TModule`

        or None if an error occured
    """
    if species.isinstance('Species'):
        #species copy
        D_species=copy.deepcopy(species)
        D_species.mutable=1
        D_species.removable=True
        if self.remove_output_when_duplicate:
            D_species.clean_type('Output') #Remove Output Types
        self.add_Node(D_species)

        #CorePromoter copy
        listIn=self.graph.predecessors(species)
        listIn.sort(key=classes_eds2.compare_node)#to be deterministic
        for interaction in listIn:
            if interaction.isinstance('CorePromoter'):
                promoter=interaction
        D_promoter=copy.deepcopy(promoter)
        D_promoter.mutable=1
        D_promoter.removable=True
        self.add_Node(D_promoter)

        #TModule copy (the only pred. of the CorePromoter)
        module=self.graph.predecessors(promoter)[0]
        D_module=copy.deepcopy(module)
        D_module.mutable=1
        D_module.removable=True
        self.add_Node(D_module)

        #Link everybody in self.graph
        self.graph.add_edge(D_module,D_promoter)
        self.graph.add_edge(D_promoter,D_species)
        return [D_module,D_promoter,D_species,module]
    else:
        print("Error : tried to duplicate something that is not a species !!")
        return None

def new_enhancer(self, species, rate, delay, parameters,basal=0.):
    """Create a complete new gene (TModule, CorePromoter and Species)

    Args:
        species(Species):
        rate (float): the rate production of the TModule
        delay (int): the delay of the CorePromoter
        parameters (list): the species parameter (see Network.new_Species)
        basal (float): the basal production of the TModule (default to 0.)

    Return:
        list: of the form [:class:`Networks.classes_eds2.TModule`, :class:`Networks.CorePromoter.CorePromoter`]
        or None if an error occured
    """
    module = classes_eds2.TModule(rate,basal)
    prom = CorePromoter(delay)
    if prom.check_grammar([module],[species]):
        self.add_CorePromoter2Species(prom,species)
        self.add_TModule2CorePromoter(module,prom)
        return [module, prom]
    else:
        print("Error in grammar, new_gene")
        return None



# Add the corepromoter functions to the Network class
setattr(classes_eds2.Network,'add_CorePromoter2Species',add_CorePromoter2Species)
setattr(classes_eds2.Network,'add_TModule2CorePromoter',add_TModule2CorePromoter)
setattr(classes_eds2.Network,'new_gene',new_gene)
setattr(classes_eds2.Network,'new_enhancer',new_enhancer)
setattr(classes_eds2.Network,'duplicate_gene',duplicate_gene)

########## Attributes attached to Mutable_Network for CorePromoter ##########

def random_gene(self,Type='Species'):
    """Create a new random gene with a species of type Type

    Args:
        Type (list): following the traditional template ['type', param]

    Return:
        list: of the form [`tmodule`, `core_promoter`, `species`] with:
            - `tModule`: :class:`Networks.classes_eds2.TModule`
            - `core_promoter`: :class:`Networks.CorePromoter.CorePromoter`
            - `species`: :class:`Networks.classes_eds2.Species`
    """
    rate = mutation.sample_dictionary_ranges('TModule.rate',self.Random)
    basal = mutation.sample_dictionary_ranges('TModule.basal',self.Random)
    delay=int(mutation.sample_dictionary_ranges('CorePromoter.delay',self.Random))
    parameters=mutation.random_parameters(Type,self.Random)
    return self.new_gene(rate, delay, parameters,basal)

def random_enhancer(self,Type='TModule'):
    """Create a new random enhancer. It includes a TModule and a CorePromoter.

    Args:
        Type (list): following the traditional template ['type', param]

    Return:
        list: of the form [`tmodule`, `core_promoter`] with:
            - `tModule`: :class:`Networks.classes_eds2.TModule`
            - `core_promoter`: :class:`Networks.CorePromoter.CorePromoter`
    """
    rate = mutation.sample_dictionary_ranges('TModule.rate',self.Random)
    basal = mutation.sample_dictionary_ranges('TModule.basal',self.Random)
    delay=int(mutation.sample_dictionary_ranges('CorePromoter.delay',self.Random))
    parameters=mutation.random_parameters(Type,self.Random)
    species = list(self.list_types['Species'])
    species = species[int(self.Random.random()*len(species))]

    return self.new_enhancer(species,rate, delay, parameters,basal)


# Add the corepromoter functions to the Mutable_Network class
setattr(mutation.Mutable_Network,'random_gene',random_gene)
setattr(mutation.Mutable_Network,'random_enhancer',random_enhancer)
