"""
Definition of degradation catalysis
Coder: M.Hemery
Creation: 2016-08-16
Last edition: 2016-10-25
"""
from phievo import __silent__,__verbose__
if __verbose__:
    print("Execute Degradation (Interaction Template)")

from . import classes_eds2
from . import mutation
import copy
from . import deriv2

#default range
mutation.dictionary_ranges['Degradation.rate'] = 0.0/mutation.T

####################################
### Degradation Class definition ###
####################################

class Degradation(classes_eds2.Interaction):
    """Catalyse the degradation of a given species

    Attributes:
        rate (float): the degradation constant
        label (str): 'Degradation' by default
        input (list): list of input types: ['Species','Degradable']
        output (list): list of output types: ['Species']
    """
    def __init__(self,rate=0.):
        """Constructor of a new Degradation

        Args:
            rate (float): the degradation constant
        """
        classes_eds2.Node.__init__(self)
        self.rate=rate
        self.label='Degradation'
        self.input=['Species','Degradable']
        self.output=['Species']

    def __str__(self):
        return "{0.id} Degradation: rate = {0.rate:.2f}".format(self)


    def outputs_to_delete(self,net):
        """indicate the Nodes to remove when deleting the Degradation

        Args:
            net (Networks): The network to which the interaction belongs

        Return:
            list: here an empty list
        """
        return []

    def check_grammar(self,input_list,output_list):
        """checks the grammar for the interactions (custom for degradation)

        Args:
            input_list (list): nodes to be checked
            output_list (list): nodes to be checked

        Return:
            bool: the consistency of up and downstream grammar
        """
        if len(input_list) != 1: return False
        if len(output_list) != 1: return False
        return output_list[0].isinstance('Degradable')

########## Attributes attached to Network for Degradation ##########

def new_Degradation(self,Input1, Input2, rate):
    """Create a new Degradation and add it to the network

    Args:
        Input1 (Species): the 'enzyme'
        Input2 (Species): the species degraded (have to be Degradable)
        rate (float): the degradation rate

    Return:
        list: of the form [Degradation]
        or None if an error occured
    """
    p=Degradation(rate)
    if p.check_grammar([Input1,Input2], [Input1]):
        self.add_Node(p)
        self.graph.add_edge(Input1,p)
        self.graph.add_edge(p,Input2)
        return [p]
    else:
        print("Error in grammar in creation of new_Degradation ")
    return None

def check_existing_Degradation(self,i1,i2):
    """Check if a Degradation exists between species i1 and i2

    Args:
        i1 (Species): the 'enzyme'
        i2 (Species): the species degraded

    Return:
        bool: if i1 is known to degrade i2
    """
    if 'Degradation' in self.list_types: #goes through the list of interactions
        for reaction in self.list_types['Degradation']:
            test = self.graph.predecessors(reaction)+self.graph.successors(reaction)
            if test == [i1,i2]:
                return True
    return False

def list_possible_Degradation(self):
    """Return the list of all possible new degradations"""
    input1_list = self.list_types['Species']
    input2_list = self.list_types['Degradable']
    list_possible = []
    for i1 in input1_list:
        for i2 in input2_list:
            if i1 != i2 and not self.check_existing_Degradation(i1,i2):
                list_possible.append([i1,i2])
    return list_possible

def number_Degradation(self):
    """Computes the number of possible Degradations"""
    return len(self.list_possible_Degradation())

# Add the corepromoter functions to the Network class
setattr(classes_eds2.Network,'new_Degradation',new_Degradation)
setattr(classes_eds2.Network,'check_existing_Degradation',check_existing_Degradation)
setattr(classes_eds2.Network,'list_possible_Degradation',list_possible_Degradation)
setattr(classes_eds2.Network,'number_Degradation',number_Degradation)

########## Attributes attached to Mutable_Network for Degradation ##########

def new_random_Degradation(self, Input1, Input2):
    """Creates a Degradation with random parameters between the Species

    Args:
        Input1 (Species): the 'enzyme'
        Input2 (Species): the species degraded (have to be Degradable)

    Return:
        list: of the form [Degradation]
    """
    rate = mutation.sample_dictionary_ranges('Degradation.rate',self.Random)
    [Deg] = self.new_Degradation(Input1,Input2,rate)
    return [Deg]

def random_Degradation(self):
    """Create new random Degradation among all possible ones

    Args:
        -

    Return:
        list: of the form [Degradation]
        or None if an error occured
    """
    if 'Degradable' in self.list_types:
        #First generate a list of all possible couples of Inputs for non existing Degradation
        list_possible = self.list_possible_Degradation()

        if list_possible:
            [Input1,Input2]=self.Random.choice(list_possible)
            [Deg]=self.new_random_Degradation(Input1,Input2)
            return [Deg]
        else:
            print("In random_Degradation : no (other) possible random_Degradation, Error")
            return None
    else:
        print("Error in random_Dummy_Interaction (try to create a Dummy_Interaction from non existing pieces)")
        return None

# Add the corepromoter functions to the Mutable_Network class
setattr(mutation.Mutable_Network,'random_Degradation',random_Degradation)
setattr(mutation.Mutable_Network,'new_random_Degradation',new_random_Degradation)

########## Integration C Tools ##########

def Degradation_deriv_inC(net):
    """gives the string corresponding to degradations for integration

    Return:
        str: a single string for all degradation in the network
    """
    if ('Degradation' in net.list_types):
        func="\n/**************Degradation interactions*****************/\n"
        for reaction in net.list_types['Degradation']:
            Input1 = net.graph.predecessors(reaction)[0]
            Input2 = net.graph.successors(reaction)[0]
            #defines interaction rate
            rate = str(reaction.rate)+' * '+Input1.id+' * '+Input2.id
            #writes down the equation in C
            func += deriv2.compute_leap([Input2.id],[],rate)
        return func
    else:
        return '' #Empty string if no degradation in net

#update deriv2
deriv2.Degradation_deriv_inC = Degradation_deriv_inC
