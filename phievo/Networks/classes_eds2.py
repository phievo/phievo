"""
Defines the main class used to describe the evolved networks
The class hierarchy is the following:

:class:`Network <phievo.Networks.classes_eds2.Network>`
 - :class:`Node <phievo.Networks.classes_eds2.Node>`:
     - :class:`Species <phievo.Networks.classes_eds2.Species>`
     - :class:`TModule <phievo.Networks.classes_eds2.TModule>`
     - :class:`Interaction <phievo.Networks.classes_eds2.Interaction>`:
        - :class:`CorePromoter <phievo.Networks.CorePromoter.CorePromoter>`
        - :class:`TFHill <phievo.Networks.TFHill.TFHill>`
        - :class:`PPI <phievo.Networks.PPI.PPI>`
        - :class:`Phosphorylation <phievo.Networks.Phosphorylation.Phosphorylation>`
        - other interactions

**Types:** Should be just the class, but for Species we have mutliple types (eg TF, Complex, Kinase, Phosphatase, Input, Output), several of which can apply at once, so the class defn in python not general enough.

**IO:** species of type = 'Input' has a defined time course supplied within the integrator.c.  Type = 'Output' are the species whos time course is used by the fitness function, they are numbered and created as genes ie with TModule and CorePromoter nodes attached to them.

**Grammar:** the rules as for what can interact with what, depends on the type of interaction and the types of inputs and outputs.  All the data for checking grammar given in the class defn of interaction.  The grammar also enters the function, Network.remove_Node().

**Class Network:** defines a bipartite graphs with adjacent nodes either 'physical-objects, segments of the genome (eg Species or TModule) or interactions. Network class then has lots of methods to add and remove nodes and edges, check the grammar rules, and output the network either as C-code or as dot diagram .

**Caps:**  Classes begin as caps, abbreviations eg TF in CAPS.  Functions within classes lc, '_' to separate names, retain Caps for embedded class names.

Arguments to functions are in order implied by directed graph, eg
check_grammar( nodes_in, node_tested, nodes_out)
add_interaction( upstream_species, interaction, downstream species)


----------------------------
"""
from phievo import __silent__,__verbose__
if __verbose__:
    print("Execute classes_eds2.py")
from phievo.initialization_code import display_error
from importlib import import_module
import phievo.networkx as nx
import numpy as np
import string,copy,sys
import pickle
import os

### Global parameters
list_unremovable=['Input'] #list of attributes that are unremovable : a species with these types can not be removed durong the evolution process

#############################
### Node class definition ###
#############################

class Node(object):
    """
    Superclass for all nodes object
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
             returns True if self is of type name
        """
        return name==self.__class__.__name__

    def list_types(self):
        """Return the list of types associated to a node"""
        return [self.__class__.__name__]

    def int_id(self):
        """ extract the integer identifer computed in Network.write_id()

        Return:
            int - the identifier of the Node
            None when valid int not found eg if write_id not called
        """
        # extract the int from the 's[1]' format
        index = self.id.split(']')[0]
        index = index.split('[')[-1]
        try:
            return int(index)
        except ValueError:
            if self.isinstance('Node'):
                return None
            else:
                display_error('Node.int_id() failed to extract integer index from id={}'.format(self.id))
                return None

    def outputs_to_delete(self,net):
        """Indicates a list of objects to delete when removing the node from the network

        Needs to be tuned specifically by all derived classes.

        Args:
            net (Network): the network self belongs to
        Return:
            list - default empty list
        """
        return []

    def isremovable(self,net,list_nodes_loop,verbose=False):
        """Check if a Node can be removed from the network

        Args:
            net (Network): the network self belongs to
            list_nodes_loop (list): to handle non tree like network
            debug (bool): Flag to activate a prolix version
        Return:
            Boolean removable or not
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
            string , default '.'
        """
        return " "

"""****************************************************************************
Definitions of physical objects on chromosome, the Species are to be time stepped,
For species input list of [Type, parameters] eg
[ [Degradable, degradation], [TF, activity], [Complex,], ... see Tags_Species dict
****************************************************************************"""

class Species(Node):
    """
    Class for any type of species, or list of species of various types
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
                      Phosphorylable = [], #only species with phosphorylated tags can be phosphorylated
                      Diffusible = ['diffusion'],
                      pMHC = [],# tags specific to the IL2 model
                      Linear_Producer=[])
    default_tags = parameters=['Degradable','Phosphorylable',"Diffusible"]
    label='Generic Species'

    def __init__(self,listtypes=[]):
        """Initialize a new Species node

        Include a check to avoid that a ligand be complexable or receptor

        Args:
            listtypes (list): following the template ['type', param.]
        """
        Node.__init__(self)
        self.types=['Species']
        for index in listtypes: #add the various tags
            if index == 'Species': continue # present by default
            if index[0] == 'Input': self.removable=False
            assert index[0] in self.Tags_Species,"Error in Species definition : no Tags "+index[0]
            self.types.append(index[0])
            for i in range(0,len(self.Tags_Species[index[0]])):
                setattr(self,self.Tags_Species[index[0]][i],index[i+1])
                # try:
                #     #updates the attributes corresponding to the types
                #     setattr(self,self.Tags_Species[index[0]][i],index[i+1])
                # except Exception:
                #     display_error('Error in Species definition tag={0} : no attributes defined'.format(index[0]))

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

        Return: return True if self is of type name
        """
        return name in self.types

    def list_types(self):
        """Return the list of types associated to a node (custom for Species)"""
        return self.types

    def def_label(self):
        """Function to write labels for graphical representation
        """
        self.label=""
        for k in self.types:
            self.label += ", "+k
            if k in self.Tags_Species:
                for item in self.Tags_Species[k]:
                    self.label += ", "+str(getattr(self,item))
            else:
                print("Error in label definition : no Tags "+k)
                return False
        return True

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
            1 if everything is done properly
            None if an error occur during the process
        """
        if not isinstance(Type,list):
            raise TypeError("Error in Species.add_type : "+str(Type)+" is not a list")

        # if Type[0] in self.types:
        #     raise Warning("in Species.add_type : Species already of Type {} doing nothing".format(Type[0]))

        if not Type[0] in self.Tags_Species:
            raise ValueError("Error in Species.add_type : no Type with name= {}".format(Type[0]))

        self.types.append(Type[0])
        
        for i,item in enumerate(self.Tags_Species[Type[0]]):
            
            try: #updates the attributes corresponding to the types
                setattr(self,item,Type[i+1])
            except Exception:
                display_error('Error in Species.add_type tag={0} require attributes {1} input as [Type, a1,..]'.format(Type[0],self.Tags_Species[Type[0]]))
        return True

class TModule(Node):
    """A TModule regulate the production of a Species, it generally binds
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
            Boolean for the consistency of up and downstream grammar
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

def check_consistency(lTypes,lNodes):
    """Check the consistency between a list of types and a list of nodes

    Typically used when constructing an interaction to check the
    biochemical grammar. For each type, it recursively checks if there
    is a corresponding node in list_nodes.

    Args:
        lTypes: the desired type of each node
        lNodes: the list of nodes

    Return:
        Boolean indicating if the consistency is OK
    """
    if len(lTypes) != len(lNodes): return False
    # cut return M without line i and column j
    cut = lambda M,i,j: [[M[x][y] for x in range(len(M)) if x != i]
                                  for y in range(len(M[0])) if y != j]
    # test if there exist a permutation p such that M[x][p(x)] is always True
    test = lambda M: (M[0][0] if len(M) == 1 else
                     sum(M[0][i] and test(cut(M,0,i)) for i in range(len(M))))
    return test([[nod.isinstance(typ) for nod in lNodes] for typ in lTypes])

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
        dict_types (dict): a dictionary indicating the Nodes of a given type (types are the keys)
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
        self.graph = nx.MultiDiGraph(selfloops=True,multiedges=True)
        self.nodes = self.graph.nodes #proxy
        self.order_node=0
        self.dict_types = dict(Output = [], Input = []) # to filled later with __build_dict_types__()
        self.hash_topology = 0
        self.title = ""
        self.Cseed=0
        self.remove_output_when_duplicate=True
        self.activator_required=False
        self.fixed_activity_for_TF= True

    def __str__(self):
        """Brief summary of the network items"""
        Nodes = sorted(self.graph.nodes(),key=lambda X:X.id)
        Species = '\n'.join([str(node) for node in Nodes if 'Species' in node.label])
        Interactions = '\n'.join([str(node) for node in Nodes if not 'Species' in node.label])
        return '### Species ###\n'+Species+'\n\n### Interactions ###\n'+Interactions+'\n'

    def __repr__(self):
        """Short object representation"""
        return '<Network object at {0}>'.format(hex(id(self)))

    def add_Node(self, node):
        """add_node to graph unless already present

        Args:
            node (:class:`Node <phievo.Networks.classes_eds2.Node>`): The node to be added

        Return:
            boolean indicating if the node has effectively been added
        """
        if self.graph.has_node(node):
            return False
        else:
            node.order=self.order_node
            self.order_node+=1
            self.graph.add_node(node)
            self.__build_dict_types__()
            return True

    def new_Species(self,types):
        """Create a new Species instance and add it to the network

        Args:
            types (list): the list types of the Species (see Species.__init__)

        Return:
            The :class:`Species <phievo.Networks.classes_eds2.Species>` which have been created
        """
        S=Species(types)
        self.add_Node(S)
        return S

###### Various generic utilities for interactions ######
    def number_nodes(self,Type):
        """count the number of Nodes of type Type

        Args:
            Type (str): the type you are looking for
        Return: The number of Nodes of types Type in dict_types
        """
        return len(self.dict_types[Type]) if Type in self.dict_types else 0

    def catal_data(self,interaction):
        """Find the reactants, catalysors, products for a catalytic interaction

        Args:
            interaction: the Interaction you are interested in

        Return:
            list of the form [catalyst,reactants,products]
        """
        listIn=self.graph.list_predecessors(interaction)
        listOut=self.graph.list_successors(interaction)
        listCata = [spc for spc in listIn if spc in listOut]
        #special case of autocatalysis
        if len(listIn)==1 and listIn[0] in listOut:
            listIn.append(listIn[0])
        for spc in listCata:
            listIn.remove(spc)
            listOut.remove(spc)
        return [listCata,listIn,listOut]

    def check_existing_binary(self,list,Type):
        """Check if a specific binary interaction of type Type already exists

        typically used for PPI

        Args:
            list (list): the reactants (Nodes) you are looking for
            Type (str): the type of Interaction you are looking for
        Return: bool
        """
        list.sort(key=compare_node)
        for inter in self.dict_types.get(Type,[]):
            inputs=self.graph.list_predecessors(inter)
            inputs.sort(key=compare_node)
            if list==inputs: return True
        return False

    def check_existing_link(self,list,Type):
        """Check if a specific interaction of type Type already exists between the elements of list

        Args:
            list (list): the reactant/product couple (Nodes) you are looking for
            Type (str): the type of Interaction you are looking for
        Return: bool
        """
        list.sort(key=compare_node)
        for inter in self.dict_types.get(Type,[]):
            inputs=self.graph.list_predecessors(inter)
            inputs.append(self.graph.list_successors(inter)[0])#add first successor
            inputs.sort(key=compare_node)
            if list==inputs: return True
        return False

    def verify_IO_numbers(self):
        """Redetermine all the input/output index

        label_them run through the list and give the correct index to all
        the items
        """
        def label_them(liszt):
            for index,species in enumerate(liszt):
                species.n_put = index

        self.__build_dict_types__()
        label_them(self.dict_types['Output'])
        label_them(self.dict_types['Input'])

###### Duplication tools ######
    def duplicate_downstream_interactions(self,species,D_species,module,D_module):
        """Called in case of gene duplication to copy the downstream interactions

        Args:
            species (:class:`Species <phievo.Networks.classes_eds2.Species>`): the 'mother' species
            D_species (:class:`Species <phievo.Networks.classes_eds2.Species>`): the 'daughter' species
            module (:class:`TModule <phievo.Networks.classes_eds2.TModule>`): the 'father' module
            D_module (:class:`TModule <phievo.Networks.classes_eds2.TModule>`): the 'son' module

        """
        listOut = sorted(self.graph.list_successors(species),key=compare_node) #careful, for self PPI, counted twice
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
            species (:class:`Species <phievo.Networks.classes_eds2.Species>`): the mother species

        Return:
            A list [D_module,D_promoter,D_species]
            D_module (:class:`TModule <phievo.Networks.classes_eds2.TModule>`): the duplicate TModule
            D_promoter (:class:`CorePromoter <phievo.Networks.CorePromoter.CorePromoter>`): the duplicate CorePromoter
            D_species (:class:`Species <phievo.Networks.classes_eds2.Species>`): the duplicate Species
        """
        #one first starts to duplicate the gene
        [D_module,D_promoter,D_species,module] = self.duplicate_gene(species)
        ###duplicate the DOWNSTREAM interactions#####
        self.duplicate_downstream_interactions(species,D_species,module,D_module)

        ###duplicate the UPSTREAM TFHills#####
        listIn=self.graph.list_predecessors(module) #look at the list_predecessors of the module before the duplication
        listIn.sort(key=compare_node)#to be deterministic
        for interaction in listIn:
            predecessor=self.graph.list_predecessors(interaction)[0]
            if interaction.isinstance('TFHill') and not self.check_existing_link([predecessor,D_module],'TFHill'):
                D_interaction=copy.deepcopy(interaction)
                D_interaction.mutable=1
                D_interaction.removable=True
                self.add_Node(D_interaction)
                #One looks for the TF upstream of this TFHill
                self.graph.add_edge(predecessor,D_interaction)
                self.graph.add_edge(D_interaction,D_module)
        #self.write_id()
        return [D_module,D_promoter,D_species]

###### Indexation tools ######
    def __build_dict_types__(self):
        """Update the dict_types dictionary of the network

        Note that it include all types (Node, Species and all Species type),
        one object can thus appear in several lists.

        """
        self.dict_types=dict(Node=self.graph.list_nodes())
        for node in self.graph.list_nodes():
            
            names = node.list_types() 
            if isinstance(node,Interaction): names.append('Interaction')
            for name in names:
                self.dict_types.setdefault(name,[]).append(node)
        for key in self.dict_types:
            self.dict_types[key].sort(key=compare_node)

    def __write_id__(self):
        """Write the ids for the network

        Update the id of all Nodes with a form n[int] for nodes
        and s[int] for species.
        """
        for index,node in enumerate(self.dict_types['Node']):
            node.id = "n[%i]"%index
        for index,species in enumerate(self.dict_types['Species']):
            species.id = "s[%i]"%index
            species.def_label()
            species.label += " Node #%i"%species.order

    def __hash_net_topology__(self):
        """Update the hash key reflecting the network topology,

        Run over nodes ordered by net.order and build hash hey from Node.label
        and Species.types (it ignore all numerical parameters)
        PRIVATE METHOD, called by __build_dict_types__() which must be current

        Return: a number comprise between 0 and sys.maxint
        """
        hh = 1
        for index in self.dict_types['Node']:
            if isinstance(index, Species):
                kk = ''.join(sorted(index.types))
                hh *= hash(kk) % sys.maxsize
            else:
                hh *= hash(index.label) % sys.maxsize
        self.hash_topology = hh

    def write_id(self):
        """Update all indexations of the network

        Return: a number comprise between 0 and sys.maxint
        """
        self.__build_dict_types__()
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
            node (:class:`Node <phievo.Networks.classes_eds2.Node>`): the node to be checked
            list_node_loops (list): to handle non tree like network
        Return: Boolean indicates if node can be safely removed
        """
        return node.isremovable(self,list_nodes_loop)

    def remove_Node(self,Node):
        """remove node from the network graph

        In case of interactions, also remove any phys objects (eg species,
        TModule) that are no more defined in absence of this interaction
        In the course of evolution, only interactions should be explicitly
        removed, then the other nodes are managed with the help of clean_nodes

        Args:
            node (:class:`Node <phievo.Networks.classes_eds2.Node>`): The node to be removed
            verbose (bool): Flag to activate the prolix mode

        Return: Boolean indicating the completion of the process
        """
        list_nodes_loop = []
        if Node in self.nodes() and self.check_Node(Node,list_nodes_loop):
            for childrens in Node.outputs_to_delete(self):
                self.graph.remove_node(childrens)
            self.graph.remove_node(Node)
            return True
        return False

    def clean_Nodes(self,verbose=False):
        """remove nodes from the network until all nodes pass the check_grammar test

        Args:
            verbose (bool): Flag to activate the prolix mode

        Return: Boolean indicating the completion of the process

        Delete any node with incorrect grammar until all remaining nodes pass test
        Currently implemented to check grammar on interaction nodes only, thus need
        remove_Node function that kills species and other phys objects that are not defined
        in absence of interaction
        """
        self.__build_dict_types__()
        modification,nloop = True,0
        while modification:
            modification=False
            nloop+=1
            if nloop > 1000: raise RuntimeError('Maximum recursion reach in clean_Nodes')
            for inter in self.dict_types.get('Interaction',[]):
                if self.graph.has_node(inter):
                    listOut=self.graph.list_successors(inter)
                    listIn=self.graph.list_predecessors(inter)
                    if not inter.check_grammar(listIn, listOut):
                        self.remove_Node(inter)
                        modification=True

    def delete_clean(self,id,target = 'interaction',verbose=False):
        """
        Remove a node according to its id and clean the network
        Warning: This operation renames all the nodes (and changes the id)

        Args:
            id: integer id of the node
            target: string either interaction or species, the type of the node to delete
        """

        #for node in self.graph.list_nodes():
        #    if node.int_id() == id:
        #        bRemove = self.remove_Node(node)
        #        if not bRemove:
        #            if verbose: print('Error while removing the node')
        #            return False
        #        bClean = self.clean_Nodes(verbose)
        #        if not bClean and verbose:
        #            if verbose: print('Error while cleaning the network')
        #        return bClean
        node=self.get_node(id,target)
        try :
            if (target=='species'):
                node.clean_type('Output')
                try:
                    node=self.graph.list_predecessors(node)[0]
                except IndexError:
                    node = node
            self.check_Node(node,[])
        except:
            print("Node=",node,"id=",id)
        if self.check_Node(node,[]):
            bRemove = self.remove_Node(node)
            if not bRemove:
                if verbose: print('Error while removing the node')
                return False
            bClean = self.clean_Nodes(verbose)
            if not bClean and verbose:
                    if verbose: print('Error while cleaning the network')
            return bClean
        if verbose: print('Node {0} not found!'.format(id))
        return False

    def get_node(self,id,target = 'interaction'):
        """
        Return the node correspoding to the specified id and target

        Args:
            id: integer id of the node
            target: string either interaction or species, the type of the node to search
        """
        if target == 'interaction':
            id_target = "n[{}]".format(id)
        elif target == 'species':
            id_target = "s[{}]".format(id)
        else:
            return None

        for node in self.graph.nodes():
            if node.id == id_target:
                return node

        return None

### Other tools ###
    def draw(self,file=None,edgeLegend=False,extended=False,display=True,return_graph=False):
        """Draw the network in a matplotlib framework

        Delegate to :func:`pretty_graph <phievo.Networks.lovelyGraph.pretty_graph>`

        Args:
            file (str): save the picture in file,
                        or print it on screen if file is None
            edgeLegend (bool): Label the graph edges
            extended (bool): Display inner modules (ex: TModules)
            display (bool): Display the figure
            return_graph(bool): Returns a graph object rather than a figure

        Examples:
            my_Network.draw('my_lovely_network.pdf')
        """
        from io import StringIO
        from matplotlib import pyplot as plt
        from matplotlib import image as mpimg
        from phievo import Networks
        if not hasattr(Networks,"pretty_graph"):
            import phievo.Networks.lovelyGraph as pretty_graph
            #setattr(Networks,"pretty_graph",lovelyGraph)
        else:
            pretty_graph = import_module(Networks.pretty_graph)
        self.write_id()
        graph = pretty_graph.pretty_graph(self,extended=extended)
        if return_graph:
            return graph
        fig = graph.draw(file,edgeLegend=edgeLegend)
        return fig

    def store_to_pickle(self,filename):
        """Save the whole network in a pickle object named filename

        Args:
            filename (str): the directory where the object is saved
        """
        with open(filename,'wb') as my_file:
            pickle.dump(self,my_file)
        print("Network save at: {}".format(filename))
