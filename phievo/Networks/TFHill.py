"""
Definition of TFHill interaction
TFHill are mainly a convenient link between TModule and their regulating
species. It is used to conserve the bipartite nature of the network.

--------
"""
from phievo import __silent__,__verbose__
if __verbose__:
    print("Execute TFHill (Interaction Template)")

from phievo.initialization_code import * #is it necessary ?
from . import classes_eds2
from . import mutation
from . import deriv2
import copy

#default range
mutation.dictionary_ranges['TFHill.hill'] = 0.0
mutation.dictionary_ranges['TFHill.threshold'] = 0.0 * mutation.C

###############################
### TFHill Class definition ###
###############################

class TFHill(classes_eds2.Interaction):
    """Implement the link between TModule and the TF

    Args:
        hill (float): the hill coefficient of the reaction
        threshold (float): the Michaelis-Menten constant
        activity (int): flag for activation (1) or repression (0)
        label (str): 'transcription' by default
        input (list): list of input types: ['TModule']
        output (list): list of output types: ['Species']
    """
    def __init__(self,hill=0,threshold=0,activity=0):
        classes_eds2.Node.__init__(self)
        self.hill=hill
        self.threshold=threshold
        self.activity=activity
        self.label = 'TF Hill Fn'
        self.input=['TF']
        self.output=['TModule']

    def __str__(self):
        act = 'activator' if self.activity else 'repressor'
        return "{0.id} TFHill {1}: H_coeff = {0.hill:.2f}, thres. = {0.threshold:.2f}".format(self,act)

    def string_param(self):
        """Self description of the Interaction"""
        return "Hill coeff=%f, Threshold=%f"%(self.hill,self.threshold)

########## Attributes attached to Network for TFHill ##########

# The transcription related interaction nodes are added along with edge connecting
# it to species, either up or down stream of it.

def add_TFHill(self, tf, inter, module):
    """Add a TF, a TModule and a :class:`TFHill <phievo.Networks.TFHill.TFHill>` interaction to the network

    Args:
        tf (:class:`Species <phievo.Networks.classes_eds2.Species>`): with the 'TF' tag
        inter (:class:`TFHill <phievo.Networks.TFHill.TFHill>`): will link tf and module
        module (:class:`TModule <phievo.Networks.classes_eds2.TModule>`): TModule to link the TFHill to

    """
    if inter.isinstance('TFHill') and inter.check_grammar([tf], [module]):
        self.add_Node(tf)
        self.add_Node(inter)
        self.add_Node(module)
        self.graph.add_edge(inter, module)
        self.graph.add_edge(tf, inter)
    else:
        print("Error in grammar add_TFHill")

def number_TFHill(self):
    """Return the number of possible TFHill"""
    n_TFHill=self.number_nodes('TFHill')
    n_TF=self.number_nodes('TF')
    n_TM=self.number_nodes('TModule')
    return n_TF*n_TM - n_TFHill

def propagate_activity_TFHill(self):
    """Ensure that TFHill activity correspond to the one of their predecessor - done for compatibility with older versions

    """
    self.write_id()
    if self.fixed_activity_for_TF:
        for tfh in self.dict_types['TFHill']:
            tf=self.graph.list_predecessors(tfh)
            tfh.activity=tf[0].activity

def new_TFHill(self, tf, hill, threshold, module, activity=0):
    """Create a new TFHill with given parameters and link it to the network.

    Args:
        tf (:class:`Species <phievo.Networks.classes_eds2.Species>`): the upstrem Species
        hill (float): the hill coefficient of the reaction
        threshold (float): the Michaelis-Menten constant
        module (:class:`TModule <phievo.Networks.classes_eds2.TModule>`): the downstream TModule
        activity (int): if fixed_activity_for_TF is True, always use the activity of tf

    Return:
        :class:`TFHill <phievo.Networks.TFHill.TFHill>`: return the new interaction or None if an error occured
    """
    if self.fixed_activity_for_TF:
        activity_tfh=tf.activity
    else:
        activity_tfh=activity
    r = TFHill(hill, threshold,activity_tfh)
    if r.check_grammar([tf], [module]):
        self.add_TFHill(tf, r, module)
    else:
        print("Error in grammar, new_TFHill")
        r = None
    return r

def duplicate_TFHill(self,D_species,interaction,module,D_module):
    """duplicate a TFHill interaction

    Args:
        D_species (:class:`Species <phievo.Networks.classes_eds2.Species>`): the new species
        interaction (:class:`TFHill <phievo.Networks.TFHill.TFHill>`): the interaction you want to duplicate
        module (:class:`TModule <phievo.Networks.classes_eds2.TModule>`): the original module
        D_module (:class:`TModule <phievo.Networks.classes_eds2.TModule>`): the new module

    """
    #copy the TFHill
    D_interaction=copy.deepcopy(interaction)
    D_interaction.mutable=1
    D_interaction.removable=True
    self.add_Node(D_interaction)
    #handle the links with the reminder of the network
    self.graph.add_edge(D_species,D_interaction)
    successors=self.graph.list_successors(interaction)
    self.graph.add_edge(D_interaction,successors[0])#TFHill only one successor so this is OK
    if (successors[0]==module):#specific case for auto regulatory feedback loop, one also needs to plug the TFHill back on the duplicated Tmodule
    # POTENTIAL BUG : check conflict with classes_eds2.duplicate_species_and_interactions.
        D_interaction_2=copy.deepcopy(interaction)
        D_interaction_2.mutable=1
        D_interaction_2.removable=True
        self.add_Node(D_interaction_2)
        self.graph.add_edge(D_species,D_interaction_2)
        self.graph.add_edge(D_interaction_2,D_module)

# Add the corepromoter functions to the Network class
setattr(classes_eds2.Network,'add_TFHill',add_TFHill)
setattr(classes_eds2.Network,'number_TFHill',number_TFHill)
setattr(classes_eds2.Network,'propagate_activity_TFHill',propagate_activity_TFHill)
setattr(classes_eds2.Network,'new_TFHill',new_TFHill)
setattr(classes_eds2.Network,'duplicate_TFHill',duplicate_TFHill)

########## Attributes attached to Mutable_Network for TFHill ##########

def new_random_TFHill(self, tf, module):
    """Creates a TFHill between tf and module with random parameters

    Args:
        tf (:class:`Species <phievo.Networks.classes_eds2.Species>`): must have the 'TF' tag
        module (:class:`TModule <phievo.Networks.classes_eds2>`): TModule associated to the TFHill

    Return:
        :class:`TFHill <phievo.Networks.TFHill.TFHill>`: return the new interaction or None if an error occured
    """
    hill = mutation.sample_dictionary_ranges('TFHill.hill',self.Random)
    threshold = mutation.sample_dictionary_ranges('TFHill.threshold',self.Random)
    activity=int(2*self.Random.random())
    return self.new_TFHill(tf, hill, threshold, module,activity)

def random_TFHill(self):
    """Creates a new :class:`TFHill <phievo.Networks.TFHill.TFHill>` among all possibles

    Return:
        :class:`TFHill <phievo.Networks.TFHill.TFHill>`: return the new interaction or None if an error occured
    """
    if 'TF' in self.dict_types and 'TModule' in self.dict_types:
        #Evaluate all possible TFHill
        possible_TFHill=[(tf,module) for module in self.dict_types['TModule']
                                     for tf in self.dict_types['TF']
                                     if not self.check_existing_link([tf,module],'TFHill')]
        n_pTFH=len(possible_TFHill)
        if not (n_pTFH==self.number_TFHill()):
            print("Potential Bug : Inconsistency in Computation of number of TFHill")
            print(n_pTFH)
            print(self.number_nodes('TFHill'))
            print(self.number_nodes('TF'))
            print(self.number_nodes('TModule'))
            for pair in possible_TFHill:
                print(pair[0].int_id(),pair[1].int_id())
            for tfh in self.dict_types['TFHill']:
                successors=self.graph.list_successors(tfh)
                predecessors=self.graph.list_predecessors(tfh)
                print(predecessors,successors)
                print(predecessors[0].int_id(),successors[0].int_id(),tfh.hill,tfh.threshold,tfh.activity)
            
        if (n_pTFH==0):
            print("In random_TFHill : No other posible TFHill")
            return None
        else:
            [tf,module]=possible_TFHill[int(self.Random.random()*n_pTFH)]
            return self.new_random_TFHill(tf, module)
    else:
        print("Error in random_TFHill (try to create TFHill from non existing pieces)")
        return None

# Add the corepromoter functions to the Mutable_Network class
setattr(mutation.Mutable_Network,'random_TFHill',random_TFHill)
setattr(mutation.Mutable_Network,'new_random_TFHill',new_random_TFHill)

########## Integration C Tools ##########

def compute_transcription(net,module):
    """Determine the transcription rate of a given module

    Used for integration in transcription_deriv_inC

    Args:
        module (:class:`TModule <phievo.Networks.classes_eds2.TModule>`): TModule to compute .

    Return:
        string the algebraic transcription rate of module
    """
    listactivator=[]
    listrepressor=[]
    if isinstance(module,classes_eds2.TModule):
        for index in net.graph.in_edges(module):
            reg=index[0] #detect the corresponding regulations
            current_activity=reg.activity
            if (current_activity==0):
                listrepressor.append("HillR(history[%i][memory][ncell],%f,%f)"%(net.graph.list_predecessors(reg)[0].int_id(),reg.threshold,reg.hill))
            else:
                listactivator.append("HillA(history[%i][memory][ncell],%f,%f)"%(net.graph.list_predecessors(reg)[0].int_id(),reg.threshold,reg.hill))
        l=len(listactivator)
        term = ""
        if(l==0):
            term="%f"%module.rate
            if hasattr(net,"activator_required"): #tests if we want to turn one genes by default from version 1.4.2
                if (net.activator_required==1):
                    term="0.00"

        if (l==1):
            term="%f*"%module.rate+listactivator[0]
        if (l>1):
            term="%f*"%module.rate
            for index in range(l-1):
                term=term+"MAX("
            term=term+listactivator[0]
            for index in range(l-1):
                term=term+","+listactivator[index+1]+")"

        if hasattr(module, "basal"): #tests on the existenc of a basal rate from version 1.3
            term="MAX("+term+",%f)"%module.basal

        if (len(listrepressor)>0):
            term=term+"*"+"*".join(listrepressor)
        return term
    else:
        print("Error in ComputeTranscription")

def transcription_deriv_inC(net):
    """gives the string corresponding to transcription for integration

    Return: A single string for all transcriptions in the network
    """
    func="\n/**************Transcription rates*****************/\n"
    func=func+" \t int k,memory=-1;\n"
    net.write_id()
    if ('TModule' in net.dict_types):
        for index in net.dict_types['TModule']:
            if isinstance(index,classes_eds2.TModule):
                trans=net.graph.list_successors(index)    #find the CorePromoter
                output=net.graph.list_successors(trans[0])    #find the transcribed protein
                func=func+"\t memory=step-%i;\n"%trans[0].delay #trans[0].delay must be an integer
                func=func+"\t if(memory>=0){\n"
                func=func+deriv2.compute_leap([],[output[0].id],compute_transcription(net,index))
                func=func+"\t}\n"
    return func

#update deriv2
deriv2.compute_transcription=compute_transcription
deriv2.interactions_deriv_inC["TFHill"] = transcription_deriv_inC
