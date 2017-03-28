"""
Define the main class used to describe the evolved networks
The class ierarchy is the following:
Node
....Species
....TModule
....Interaction
........Corepromoter
........PPI
........<other interactions>
Network

Types: Should be just the class, but for Species we have mutliple types (eg TF, Complex, Kinase, Phosphatase, Input, Output), several of which can apply at once, so the class defn in python not general enough.

IO: species of type = 'Input' has a defined time course supplied within the integrator.c.  Type = 'Output' are the species whos time course is used by the fitness function, they are numbered and created as genes ie with TModule and CorePromoter nodes attached to them.

Grammar: the rules as for what can interact with what, depends on the type of interaction and the types of inputs and outputs.  All the data for checking grammar given in the class defn of interaction.  The grammar also enters the function, Network.remove_Node().

Class Network, then defines a bipartite graphs with adjacent nodes either 'physical-objects, segments of the genome (eg Species or TModule) or interactions. Network class then has lots of methods to add and remove nodes and edges, check the grammar rules, and output the network either as C-code or as dot diagram .


Caps:  Classes begin as caps, abbreviations eg TF in CAPS.  Functions within classes lc, '_' to separate names, retain Caps for embedded class names.

Arguments to functions are in order implied by directed graph, eg
check_grammar( nodes_in, node_tested, nodes_out)
add_interaction( upstream_species, interaction, downstream species)
"""
from phievo import __silent__,__verbose__
if __verbose__:
    print("Execute classes_eds2.py")
from phievo.initialization_code import display_error
import networkx as NX    # keep name spaces distinct
import numpy as np
import string,copy,sys
import pickle
import os

### Global parameters
debugging_node_removal = False
list_unremovable=['Input'] #list of attributes that are unremovable : a species with these types can not be removed durong the evolution process

#############################
### Node class definition ###
#############################

class Node(object):
    """Superclass for all nodes object
    """
    id = 'None' # inherited by all derived classes and instances

    def __init__(self):
        self.label='Generic Node'
        self.id= 'None'
        self.order=0 #order of the node in the network
        self.removable=True#tag to know if we can remove a node,1  by default means the node can be removed
        self.mutable=1 # 1 by default means the parameters can be changed, used in mutation.rand_modify()

    def __str__(self):
        return str(self.id)+' : Generic Node Object'

    def print_node(self):
        """print a full description of the current node

        Return:
            None: everything is diplayed through stdout
        """
        str_class = str(self.__class__).split("'")[-2]
        label = self.__dict__.get('label',None)
        print('class=', str_class, '\t', 'label=', label) #main information
        other = [str(k)+'='+str(v) for k,v in self.__dict__.items() if k != 'label']
        print(*other,sep='\t') #other information

    def isinstance(self,name):
        """check the type of a node

        Customed the builtin isinstance(derived_class, base_class)
        for the general case

        Args:
            name: the type to be tested

        Return:
            bool: return True if self is of type name
        """
        return name==self.__class__.__name__

    def list_types(self):
        """Return the list of types associated to a node"""
        return [self.__class__.__name__]

    def int_id(self):
        """ extract the integer identifer computed in Network.write_id()

        Return:
            (int): the identifier of the Node
            None: when valid int not found eg if write_id not called
        """
        # extract the int from the 's[1]' format
        index = self.id.split(']')[0]
        index = index.split('[')[-1]
        try:
            return int(index)
        except Exception:
            display_error('Node.int_id() failed to extract integer index from id={}'.format(self.id))
            return None

    def outputs_to_delete(self,net):
        """Indicates a list of objects to delete when removing the node from the network

        Needs to be tuned specifically by all derived classes.

        Args:
            net (Network): the network self belongs to
        Return:
            (list): default empty list
        """
        return []

    def isremovable(self,net,list_nodes_loop,verbose=debugging_node_removal):
        """Check if a Node can be removed from the network

        Args:
            net (Network): the network self belongs to
            list_nodes_loop (list): to handle non tree like network
            debug (bool): Flag to activate a prolix version
        Return:
            (bool): if self is removable or not
        """

        if self in list_nodes_loop:
            # We face a cycle within the network.
            # Break out of this by saying that we can remove the whole thing.
            return True

        list_nodes_loop.append(self)
        Bool=self.removable
        for type in list_unremovable:
            Bool=Bool*(1-self.isinstance(type))

        if verbose:
            print("\nAs of now the node:")
            for property, value in vars(self).items():
                print(property, ": ", value)
            print("Can be removed? %d" %Bool)

        list=self.outputs_to_delete(net)

        if verbose:
            if not list:
                print("No outputs to remove, should be ok.")
        if Bool:
            if verbose:
                if list:
                    print("\nNeed to check if its outputs can be deleted though.")
            for Species in list:
                if verbose:
                    print("Outputs to delete (from isremovable interaction):")
                    for property, value in vars(Species).items():
                        print(property, ": ", value)
                node_removal_bool = Species.isremovable(net,list_nodes_loop)
                Bool=Bool*node_removal_bool

                if verbose:
                    print("\nSo, is it possible to delete the output:")
                    for property, value in vars(Species).items():
                        print(property, ": ", value)
                    print("\nof:")
                    for property, value in vars(self).items():
                        print(property, ": ", value)
                    print("\n??? Answer: %d" %node_removal_bool)

        if Bool:
            list_nodes_loop.remove(self)
        return Bool

    def string_param(self):
        """Returns a function with parameters for the nodes

        Mainly here to be customized in subclasses

        Return:
            string: default '.'
        """
        return " "

"""****************************************************************************
Definitions of physical objects on chromosome, the Species are to be time stepped,
For species input list of [Type, parameters] eg
[ [Degradable, degradation], [TF, activity], [Complex,], ... see Tags_Species dict
****************************************************************************"""

class Species(Node):
    """Class for any type of species, or list of species of various types
       Input list of lists eg [ [Degradation, degradation], [TF, activity], [Complex,],
       [Kinase],..
       [Output, n_put], [Input, n_put] ]  where n_put is an integer enumerating IO
       The first tag of ['Species'] is assumed and should not be input
    """
    #Dictionnary of possible tags with list of corresponding attributes
    Tags_Species=dict(Species = [],
                      Degradable = ['degradation'], #degration constant of species
                      TF = ['activity'], #activity of TF; 1 is activator, 0 repressor
                      Kinase = [],
                      Phosphatase = [],
                      Output = ['n_put'], #index as an output
                      Input = ['n_put'], #index as an input
                      Complexable = [],
                      Complex = [],
                      Ligand = [], # ligand diffuses only if diffusible tag is on
                      Receptor = [],
                      Phospho = ['n_phospho'], #Phosphorylated species #For immune case
                      Phosphorylable = [], #only species with phosphorylated tags can be phosphorylated
                      Diffusible = ['diffusion'],
                      pMHC = [],
                      # tags specific to the IL2 model
                      Common = ['common'], # 1 if the species is a common good, 0 otherwise.
                      Linear_Producer=[])
    label='Generic Species'

    def __init__(self,listtypes=[]):
        """Initialize a new Species node

        Include a check to avoid that a ligand be complexable or receptor

        Args:
            listtypes (list): following the template ['type', param.]
        """
        Node.__init__(self)
        self.types=['Species']

        #Forbidden mutliple types
        #for simplicity, we do not want that a Ligand is Complexable or is a Receptor
        if 'Ligand' in listtypes and 'Complexable' in listtypes:
            print("Warning : a ligand with a Complexable Tag was created. I remove the Complexable tag")
            listtypes.remove('Complexable')
        if 'Ligand' in listtypes and 'Receptor' in listtypes:
            print("Warning : a ligand with a Receptor Tag was created. I remove the Receptor tag")
            listtypes.remove('Complexable')

        #Now we add the various tags
        for index in listtypes:
            if index == 'Species': continue # present by default
            if index[0] == 'Input': self.removable=False
            if index[0] in self.Tags_Species:
                self.types.append(index[0])
                for i in range(0,len(self.Tags_Species[index[0]])):
                    try:
                        #updates the attributes corresponding to the types
                        setattr(self,self.Tags_Species[index[0]][i],index[i+1])
                    except Exception:
                        display_error('Error in Species definition tag={0} : no attributes defined'.format(index[0]))
            else:
                print("Error in Species definition : no Tags "+index[0])

    def __str__(self):
        def cutter(arg):
            if isinstance(arg,float):
                return '{0:.2f}'.format(arg)
            elif isinstance(arg,int):
                return '{0}'.format(arg)
            elif isinstance(arg,str):
                return '{0:.5}.'.format(arg)
            else:
                return arg

        res = ''
        for key in self.types:
            if key == 'Species':continue
            res += key
            if self.Tags_Species[key]:
                res+=' ('+', '.join([cutter(name)+"="+cutter(getattr(self,name)) for name in self.Tags_Species[key]])+'), '
            else:
                res += ', '
        return str(self.id)+" Species: "+res[:-2]

    def isinstance(self,name):
        """check the type of a node

        Customed the builtin isinstance(derived_class, base_class)
        for the Species class

        Args:
            name: the type to be tested

        Return:
            bool: return True if self is of type name
        """
        return name in self.types

    def list_types(self):
        """Return the list of types associated to a node (custom for Species)"""
        return self.types

    def def_label(self):
        """Function to write labels for graphical representation

        Return:
            None: this function update the label attribute
        """
        self.label=""
        for k in self.types:
            self.label += ", "+k
            if k in self.Tags_Species:
                for item in self.Tags_Species[k]:
                    self.label += ", "+str(getattr(self,item))
            else:
                print("Error in label definition : no Tags "+k)

    def clean_type(self,Type):
        """Removes a type and corresponding attributes from a species"""
        if Type in self.types:
            self.types.remove(Type) #remove the type
            for item in self.Tags_Species[Type]: #remove the corresponding attributes
                delattr(self,item)

    def change_type(self,Type,parameters):
        """ Change the parameters of a type

        Args:
            Type (string): name of the type to modify
            parameters (list): list of the new parameters as defined in the Tag_Species
        """
        if len(parameters) == len(self.Tags_Species[Type]):
            for i,param in enumerate(self.Tags_Species[Type]):
                setattr(self,param,parameters[i])
        else:
            raise ValueError("parameters must be a list of the same length as Tag_Species. Make sure you gave the correct number of parameters.")


    def add_type(self,Type):
        """add Type and its corresponding parameters

        Several layer of check are done before the core function
        to insure that Type is correctly added
        Also used to add the output/input tag to species. e.g.:species.add_type(['Output',n_put])

        Args:
            Type (list): must be provided in a list of the form
            ['Tag',parameter1,parameter2] as defined in Tags_Species
        Return:
            1: if everything is done properly
            None: if an error occur during the process
        """
        if not isinstance(Type,list): #catch the case where Type is not a list
            print("Error in Species.add_type : "+str(Type)+" is not a list")
            return None

        if Type[0] in self.types: #catch the case where Type is already present
            print("Warning in Species.add_type : Species already of Type "+Type[0]+"; doing nothing")
            return None

        if not Type[0] in self.Tags_Species: #catch the case where Type is unknown
            print("Error in Species.add_type : no Type with name= "+Type[0])
            print("Allowed Types are:", list(self.Tags_Species.keys()))
            return None

        self.types.append(Type[0])
        for i,item in enumerate(self.Tags_Species[Type[0]]):
            try: #updates the attributes corresponding to the types
                setattr(self,item,str(Type[i+1]))
            except Exception:
                display_error('Error in Species.add_type tag={0} require attributes {1} input as [Type, a1,..]'.format(Type[0],self.Tags_Species[Type[0]]))
        return 1

class TModule(Node):
    """Definition of the Transcription Factor (TF)

    A TModule regulate the production of a Species, it generally binds
    upstream to a CorePromoter (direct production) or a TFHill
    (regulation) and downstream to another TFHill which point to the
    product Species.

    Parameters:
        rate (float): the production rate to be regulated
        basal (float): the basal production rate
    """
    def __init__(self,rate=0,basal=0):
        Node.__init__(self)
        self.rate = rate
        self.basal= basal #basal rate from version 1.3
        self.label='TModule'

    def __str__(self):
        return "{0.id} TModule: rate = {0.rate:.2f}, basal = {0.basal:.2f}".format(self)

    def string_param(self):
    	return "Rate=%f"%self.rate

class Interaction(Node):
    """Interaction class derived from Node, defines interaction between Species or TModule"""
    def __str__(self):
        return "{} : {}".format(str(self.id),self.label)

    def check_grammar(self,input_list,output_list):
        """checks the grammar for the interactions

        Args:
            input_list (list): nodes to be checked
            output_list (list): nodes to be checked

        Return:
            bool: the consistency of up and downstream grammar
        """
        input_check = check_consistency(self.input,input_list)
        output_check = check_consistency(self.output,output_list)
        return input_check and output_check

################################
### Node function definition ###
################################

def compare_node(x):
    """Used to order nodes in arbitrary but deterministic order when needed"""
    return x.order

def check_consistency(list_types,list_nodes):
    """Check the consistency between a list of types and a list of nodes

    Typically used when constructing an interaction to check the
    biochemical grammar. For each type, it recursively checks if there
    is a corresponding node in list_nodes.

    Args:
        list_types: the desired type of each node
        list_nodes: the list of nodes

    Return:
        bool: Indicating if the consistency is OK
    """
    if not (len(list_types)==len(list_nodes)): return False
    if (len(list_types)<=1):
        try:
            return list_nodes[0].isinstance(list_types[0])
        except Exception:
            display_error("Problem in check_consistency at the end of recursion")
            return False
    else:
        for current_type in list_types:
            for k,current_node in enumerate(list_nodes):
                if current_node.isinstance(current_type):
                    #copy of the list of types without the first occurence of current_type
                    list_types_short = copy.deepcopy(list_types)
                    list_types_short.remove(current_type)
                    #copy of the list of nodes without the node to remove
                    list_nodes_short = [list_nodes[l] for l in range(len(list_nodes)) if l!=k]
                    # consistency is checked if the rest of the lists are consistent
                    if check_consistency(list_types_short,list_nodes_short): return True
        return False #if ends the loop, that means that there is no type consistency

################################
### Network class definition ###
################################

class Network(object):
    """ Complete description of a network of interactions.

    It is represented as a bipartite graph between the biochemical species and
    the interactions. The very description is stored in the `graph` attribute.

    Note that each interaction import add new methods to the Network class.

    Attributes:
        graph (networkx.MultiDiGraph): the network properly speaking
        order_node (int): index to keep track of the order of the nodes
        list_types (dict): a dictionary indicating the Nodes of a given type (types are the keys)
        hash_topology (int): to index the topologies (see __hash_net_topology__)
        title (str): for graphing network and to hold misc info
        Cseed (int): random seed for the integration in C
        remove_output_when_duplicate (bool): if you want to remove Output tag when duplicating genes
        activator_required (bool): if an activator is required to get any gene product
        fixed_activity_for_TF (bool): if a TF either an activator or repressor (if False, they can do both)

    Main functions:
        add_* methods just add objects to the graph
        new_* create and add objects (usually by calling add_* method)
    """
    def __init__(self):
        """The constructor of the Network, default settings
        See Network for complete doc
        """
        self.versionnx=NX.__version__
        if (int(self.versionnx[0])<1):
            self.graph = NX.XDiGraph(selfloops=True,multiedges=True)
        else:
            self.graph = NX.MultiDiGraph(selfloops=True,multiedges=True)
        self.order_node=0
        self.list_types = dict(Output = [], Input = []) # to filled later with __build_list_types__()
        self.hash_topology = 0
        self.title = ""
        self.Cseed=0
        self.remove_output_when_duplicate=False
        self.activator_required=False
        self.fixed_activity_for_TF= True

    def __str__(self):
        """Brief summary of the network items"""
        Nodes = sorted(self.graph.nodes(),key=lambda X:X.id)
        Species = '\n'.join([str(node) for node in Nodes if 'Species' in node.label])
        Interactions = '\n'.join([str(node) for node in Nodes if not 'Species' in node.label])
        return '### Species ###\n'+Species+'\n\n### Interactions ###\n'+Interactions+'\n'

    def __repr__(self):
        """Copy the __str__ method"""
        return self.__str__()

    def add_Node(self, node):
        """add_node to graph unless already present

        Args:
            node (:class:`Networks.classes_eds2.Node`): The node to be added

        Return:
            bool: indicate if the node has effectively been added
        """
        if self.graph.has_node(node):
            return False
        else:
            node.order=self.order_node
            self.order_node+=1
            self.graph.add_node(node)
            self.__build_list_types__()
            return True

    def new_Species(self,types):
        """Create a new Species instance and add it to the network

        Args:
            types (list): the list types of the Species (see Species.__init__)

        Return:
            Species: the species which have been created
        """
        S=Species(types)
        self.add_Node(S)
        return S

###### Various generic utilities for interactions ######
    def number_nodes(self,Type):
        """count the number of Nodes of type Type

        Args:
            Type (str): the type you are looking for
        Return:
            int: the number of Nodes of types Type in list_types
        """
        return len(self.list_types[Type]) if Type in self.list_types else 0

    def catal_data(self,interaction):
        """Find the reactants, catalysors, products for a catalytic interaction

        Args:
            interaction: the Interaction you are interested in

        Return:
            list: of the form [catalyst,reactants,products]
        """
        listIn=self.graph.predecessors(interaction)
        if (len(listIn)==1): #special case of autocatalysis
            listIn.append(listIn[0])
        listOut=self.graph.successors(interaction)
        for x in listIn:
            if x in listOut:
                catalyst=x #what if their is two catalysts ???
        if 'catalyst' in locals():
            listIn.remove(catalyst)
            listOut.remove(catalyst)
            return [catalyst,listIn,listOut]
        else:
            return [[],listIn,listOut]

    def check_existing_binary(self,list,Type):
        """Check if a specific binary interaction of type Type already exists

        typically used for PPI or LR

        Args:
            list (list): the reactants (Nodes) you are looking for
            Type (str): the type of Interaction you are looking for
        Return:
            bool
        """
        list.sort(key=compare_node)
        if Type in self.list_types: #goes through the list of interactions
            for inter in self.list_types[Type]:
                inputs=self.graph.predecessors(inter)#check the inputs
                inputs.sort(key=compare_node)
                if (list==inputs):#if inputs are the same, the interaction already exists
                    return True
        return False

    def check_existing_link(self,list,Type):
        """Check if a specific interaction of type Type already exists between the elements of list

        Args:
            list (list): the reactant/product couple (Nodes) you are looking for
            Type (str): the type of Interaction you are looking for
        Return:
            bool
        """
        list.sort(key=compare_node)
        if Type in self.list_types:#goes through the list of interactions
            for inter in self.list_types[Type]:
                inputs=self.graph.predecessors(inter)#check the inputs
                inputs.append(self.graph.successors(inter)[0])#adds the first successor
                inputs.sort(key=compare_node)
                if (list==inputs):#if inputs are the same, the interaction already exists
                    return True
        return False

    def verify_IO_numbers(self):
        """Redetermine all the input/output index

        label_them run through the list and give the correct index to all
        the items

        Return:
            None: in place modification
        """
        def label_them(liszt):
            for index,species in enumerate(liszt):
                species.n_put = index

        self.__build_list_types__()
        label_them(self.list_types['Output'])
        label_them(self.list_types['Input'])

###### Duplication tools ######
    def duplicate_downstream_interactions(self,species,D_species,module,D_module):
        """Called in case of gene duplication to copy the downstream interactions

        Args:
            species (Species): the 'mother' species
            D_species (Species): the 'daughter' species
            module (TModule): the 'father' module
            D_module (TModule): the 'son' module

        Return:
            None: in place modification
        """
        listOut = sorted(self.graph.successors(species),key=compare_node) #careful, for self PPI, counted twice
        already_seen_PPI = [] #to keep a list of the PPI already considered
        for interaction in listOut:
            if interaction.isinstance('TFHill'):
                self.duplicate_TFHill(D_species,interaction,module,D_module)
            if interaction.isinstance('PPI') and not interaction in already_seen_PPI:
                self.duplicate_PPI(species,D_species,interaction,module,D_module)
                already_seen_PPI.append(interaction)

    def duplicate_species_and_interactions(self,species):
        """Called to duplicate a species with its interactions

        Right now only duplicates downstream TFHills and PPI and
        upstream TFHills.
        The input&output tags are removed from duplicate gene
        (see self.remove_output_when_dulicate)

        Args:
            species (Species): the mother species

        Return:
            A list [D_module,D_promoter,D_species]
            D_module (TModule): the duplicate TModule
            D_promoter (CorePromoter): the duplicate CorePromoter
            D_species (Species): the duplicate Species
        """
        print("Duplicate")
        #one first starts to duplicate the gene
        [D_module,D_promoter,D_species,module] = self.duplicate_gene(species)
        if D_species.isinstance('Output') and self.remove_output_when_duplicate:
            D_species.clean_type('Output') #clean output tag
        if D_species.isinstance('Input'):
            D_species.clean_type('Input') #clean input tag

        ###duplicate the DOWNSTREAM interactions#####
        self.duplicate_downstream_interactions(species,D_species,module,D_module)

        ###duplicate the UPSTREAM TFHills#####
        listIn=self.graph.predecessors(module) #look at the predecessors of the module before the duplication
        listIn.sort(key=compare_node)#to be deterministic
        for interaction in listIn:
            if interaction.isinstance('TFHill'):
                D_interaction=copy.deepcopy(interaction)
                D_interaction.mutable=1
                D_interaction.removable=True
                self.add_Node(D_interaction)
                #One looks for the TF upstream of this TFHill
                predecessor=self.graph.predecessors(interaction)[0]
                self.graph.add_edge(predecessor,D_interaction)
                self.graph.add_edge(D_interaction,D_module)
        #self.write_id()
        return [D_module,D_promoter,D_species]

###### Indexation tools ######
    def __build_list_types__(self):
        """Update the list_type dictionary of the network

        Note that it include all types (Node, Species and all Species type),
        one object can thus appear in several lists.

        Return:
            None: in place modification
        """
        self.list_types=dict(Node=self.graph.nodes())
        for index in self.graph.nodes():
            names = index.list_types()
            for name in names:
                self.list_types.setdefault(name,[]).append(index)
            if isinstance(index,Interaction):
                self.list_types.setdefault('Interaction',[]).append(index)
        for key in self.list_types:
            self.list_types[key].sort(key=compare_node) #sort the lists

    def __write_id__(self):
        """Write the ids for the network

        Update the id of all Nodes with a form n[int] for nodes
        and s[int] for species.

        Return:
            None: in place modification
        """
        for index,node in enumerate(self.list_types['Node']):
            node.id = "n[%i]"%index
        for index,species in enumerate(self.list_types['Species']):
            species.id = "s[%i]"%index
            species.def_label()
            species.label += " Node #%i"%species.order

    def __hash_net_topology__(self):
        """Update the hash key reflecting the network topology,

        Run over nodes ordered by net.order and build hash hey from Node.label
        and Species.types (it ignore all numerical parameters)
        PRIVATE METHOD, called by __build_list_types__() which must be current

        Return:
            int: a number comprise between 0 and sys.maxint
        """
        hh = 1
        for index in self.list_types['Node']:
            if isinstance(index, Species):
                kk = ''.join(sorted(index.types))
                hh *= hash(kk) % sys.maxsize
            else:
                hh *= hash(index.label) % sys.maxsize
        self.hash_topology = hh

    def write_id(self):
        """Update all indexations of the network

        Return:
            int: a number comprise between 0 and sys.maxint
        """
        self.__build_list_types__()
        self.__hash_net_topology__()
        self.__write_id__()

###### Removal tools ######
    def check_Node(self,node,list_nodes_loop):
        """Check if a Node can be removed from the network

        Delegate to Node.isremovable
        check if node is not an input/output or a node uniquely
        and directly upstream of a nonremovable species
        (eg part of output gene)

        Args:
            node (:class:`Networks.classes_eds2.Node`): the node to be checked
            list_node_loops (list): to handle non tree like network
        Return
            bool: indicate if node can be safely removed
        """
        return node.isremovable(self,list_nodes_loop)

    def version_remove_Node(self,node):
        """Remove a node form self.graph
        Used to handle the diff. version of networkx

        Args:
            node (:class:`Networks.classes_eds2.Node`): the node to be deleted

        Return:
            None: in place modification
        """
        if (int(self.versionnx[0])<1):
            self.graph.delete_node(node)
        else:
            self.graph.remove_node(node)

    def remove_Node(self,Node,verbose=debugging_node_removal):
        """remove node from the network graph

        In case of interactions, also remove any phys objects (eg species,
        TModule) that are no more defined in absence of this interaction
        In the course of evolution, only interactions should be explicitly
        removed, then the other nodes are managed with the help of clean_nodes

        Args:
            node (:class:`Networks.classes_eds2.Node`): The node to be removed
            verbose (bool): Flag to activate the prolix mode

        Return:
            bool: indicating the completion of the process
        """
        list_nodes_loop = []
        if (self.graph.has_node(Node)) and (self.check_Node(Node,list_nodes_loop)):
            if verbose:
                print("\n\n\n")
                print("===================================")
                print("===== We are in remove_Node =======")
                print("===================================")
                print("\nWe are removing:")
                for property, value in vars(Node).items():
                    print(property, ": ", value)
                print("\nAs well as all its output. I.e.:")

            for phys_obj in Node.outputs_to_delete(self): #remove the products of the interaction
                if verbose:
                    for property, value in vars(phys_obj).items():
                        print(property, ": ", value)
                self.version_remove_Node(phys_obj)
            self.version_remove_Node(Node) #remove the node
            if verbose:
                print("===================================")
                print("===== We leave remove_Node ========")
                print("===================================")
                print("\n\n\n")
            return True
        return False

    def clean_Nodes(self,verbose=debugging_node_removal):
        """remove nodes from the network until all nodes pass the check_grammar test

        Args:
            verbose (bool): Flag to activate the prolix mode

        Return:
            bool: indicating the completion of the process

        Delete any node with incorrect grammar until all remaining nodes pass test
        Currently implemented to check grammar on interaction nodes only, thus need
        remove_Node function that kills species and other phys objects that are not defined
        in absence of interaction
        """
        self.__build_list_types__()
        modification,nloop = True,0
        while modification and nloop<1000:
            modification=False
            nloop+=1
            for inter in self.list_types.get('Interaction',[]):
                if self.graph.has_node(inter):
                    listOut=self.graph.successors(inter)
                    listIn=self.graph.predecessors(inter)
                    if not inter.check_grammar(listIn, listOut):
                        if verbose:
                            print("We have some inconsistent grammar with node:")
                            for property, value in vars(inter).items():
                                print(property, ": ", value)
                        for l in inter.outputs_to_delete(self):
                            self.remove_Node(l)
                        self.remove_Node(inter)
                        modification=True

        if (nloop>900):
            print("Potential infinite loop")
            P=draw_Network(self.graph)
            name="Infiniteloop.jpg"
            P.write_jpg(name)

    def draw(self,file=None,edgeLegend=False,extended=False):
        """Draw the network in a matplotlib framework

        Delegate to :class:`Networks.lovelyGraph.pretty_graph`

        Args:
            file (str): save the picture in file,
                        or print it on screen if file is None

        Returns:
            None

        Examples:
            my_Network.draw('my_lovely_network.pdf')
        """
        from io import StringIO
        from matplotlib import pyplot as plt
        from matplotlib import image as mpimg
        from phievo.Networks.lovelyGraph import pretty_graph
        self.write_id()
        graph = pretty_graph(self,extended=extended)
        graph.draw(file,edgeLegend=edgeLegend)

    def run_dynamics(self,project_dir,inits,trial=1,erase_buffer=True):
        """ This function will call the C functions written for the evolution in order to generate a new dynamics.

        Args:
            file : None

        Returns:
            data (dict): One of the key returns the times steps and the others are the index of the corresponding run. A key corresponding to a run returns a other dictionnary where each key stands for a cell. The data stored for a run cell is the time time course of the species.
        """
        if project_dir[-1] != os.sep:
            project_dir = project_dir + os.sep
        workplace_dir = project_dir + "Workplace/"
        working_dir = os.getcwd() + os.sep

        cfile = inits.cfile
        prmt = inits.prmt
        prmt["ntries"] =trial
        from importlib import import_module
        deriv2 = import_module("phievo.Networks.deriv2")

        here, this_filename = os.path.split(__file__)
        here += "/"
        # Define default directory for cfile then overwrite with information from inits
        deriv2.cfile['header'] = here +  '../CCodes/integrator_header.h'
        deriv2.cfile['utilities'] = here + '../CCodes/utilities.c'
        deriv2.cfile['geometry'] = here + '../CCodes/linear_geometry.c'
        deriv2.cfile['integrator'] = here + '../CCodes/euler_integrator.c'
        deriv2.cfile['main'] = here + '../CCodes/main_general.c'

        for k, v in cfile.items():
            path = working_dir + v
            if os.path.isfile(path):
                deriv2.cfile[k] = path
            else:
                raise FileNotFoundError("ERROR to find the C code:\n{} doesn't match a file.".format(path))
                sys.exit()

        deriv2.workplace_dir = workplace_dir
        if ('langevin_noise' in prmt):
            if (prmt['langevin_noise'] > 0):
                deriv2.noise_flag = 1
        deriv2.compile_and_integrate(self,prmt,1000,True)

        data = {"time":np.arange(0,prmt["dt"]*(prmt["nstep"]),prmt["dt"])}
        N_species = len(self.list_types['Species'])
        N_cell = prmt["ncelltot"]
        for i in range(trial):
            temp = np.genfromtxt('Buffer%d'%i, delimiter='\t')[::,1:]
            data[i] = {cell:temp[::,cell:cell+N_species] for cell in range(N_cell)}
            if erase_buffer:
                os.remove("Buffer%d"%i)
            else:
                os.rename("Buffer{0}".format(i),"{to},Buffer{index}".format(to=project_dir,index=i))

        return data

    def store_to_pickle(self,filename):
        """Save the whole network in a pickle object named filename

        Args:
            filename (str): the directory where the object is saved

        Returns:
            None: in place saving
        """
        with open(filename,'wb') as my_file:
            pickle.dump(self,my_file)
        print("Network save at: {}".format(filename))

    def network_to_string(self):
        """Return a string defining the network net as a python file

        Return:
            str: the network description
        """
        term="import random\ng=random.Random(42)\nL=Mutable_Network(g)\nL.Cseed=%i\n"%self.Cseed #initialization of the network
        name=""
        listarrows=""
        listnodes=self.graph.nodes()
        listnodes.sort(key=compare_node)#to be deterministic
        for index in listnodes:
            str_class = str(index.__class__).split('.')[-1]
            for k,v in index.__dict__.items():
                if k == 'id':
                    name=str(v)
                    name=name.replace('[','')
                    name=name.replace(']','')
                    term=term+name+"="+str_class+"()\n"#initializes the nodes
            for k,v in index.__dict__.items() :
                if (k == 'label') or (k=='id'):
                    continue
                term=term+name+"."+k+'='+str(v)+"\n"#give attributes to the nodes
            term=term+"L.add_Node("+name+")\n"#adds nodes
            listsuccessors=self.graph.successors(index)
            listsuccessors.sort(key=compare_node)#to be deterministic
            for successor in listsuccessors:
                outname=str(successor.id).replace('[','')
                outname=outname.replace(']','')
                listarrows=listarrows+"L.graph.add_edge("+name+","+outname+")\n"#adds links
        if hasattr(self,"activator_required"):
            term=term+"L.activator_required=%i\n"%self.activator_required
        if hasattr(self,"noise_level"):
            term=term+"L.noise_level=%i\n"%self.noise_level
        if hasattr(self,"fixed_activity_for_TF"):
            term=term+"L.fixed_activity_for_TF=%i\n"%self.fixed_activity_for_TF
        return term+listarrows

####################################
### Printing function definition ###
####################################

# functions to be deprecated, using shelve object instead
def retrieve_from_pickle(filename,verbose=True):
    """Retrieve a whole network from a pickle object named filename

    Args:
        filename (str): the directory where the object is saved

    Returns:
        Network: the object having been stored
    """
    with open(filename,'rb') as my_file:
        net = pickle.load(my_file)
    if verbose:
        print("Network retrieve from: {}".format(filename))
    return net

def print_Network(net):
    """Return a string defining the network net as a python file

    delegate to Network.__repr__

    Args:
        net (Network): the network you want to write

    Return:
        (str): the network description
    """
    return net.__repr__()

def str2Network(net_str):
    """Return the network object described by net_str

    Args:
        net_str (str): A description of a network as produced by print_Network

    Return:
        Network: the corresponding network
    """
    exec_dict = {}
    exec_str = 'from phievo.Networks.classes_eds2 import *\nfrom phievo.Networks.mutation import *\nfrom phievo.Networks.interaction import *\nfrom phievo.Networks.mutation import Mutable_Network\n' + ''.join( net_str )
    exec(exec_str, exec_dict) # the network will be L from best net file
    net = exec_dict['L']   # L is the name of network used in best net file
    net.write_id()
    return net

if __name__ == "__main__":
    # traditionnal canned test
    X = Network()
    X.store_to_pickle('dummy.ntk')
    Y = retrieve_from_pickle('dummy.ntk')
