"""
Definition of Phosphorylation interaction

/!\ WARNING: IF USING THIS CLASS PUT config.multiple_phospho to 0, otherwise you might have bugs (for now)
TODO: in New Phosphorylation, test on n_phospho; if it is 1 (or higher than something)
then remove Phosphorylable. Also update n_phospho accordingly when phosphorylated

Creation: unknown
Last edition: 2016-10-26
"""
from phievo import __silent__,__verbose__
if __verbose__:
    print("Execute Phosphorylation (Interaction Template)")

from . import classes_eds2
from . import mutation
from . import deriv2
import copy

#default range
mutation.dictionary_ranges['Phosphorylation.rate'] = 0.0/mutation.T
mutation.dictionary_ranges['Phosphorylation.hill'] = 0.0
mutation.dictionary_ranges['Phosphorylation.threshold'] = 0.0 * mutation.C
mutation.dictionary_ranges['Phosphorylation.dephosphorylation'] = 0.0 / mutation.T

########################################
### Phosphorylation Class definition ###
########################################

class Phosphorylation(classes_eds2.Interaction):
    """Phosphorylation interaction

    Attributes:
        rate (float): the phosphorylation rate
        threshold (float): the Michaelis-Menten constant
        hill (float): the hill coefficient of the reaction
        dephosphorylation (float): the dephosphorylation rate
        label (str): 'Phosphorylation' by default
        input (list): list of input types: ['Kinase','Phosphorylable']
        output (list): list of output types: ['Kinase','Phospho']
    """
    def __init__(self,rate=0,threshold=1,hill_coeff=1,dephospho_rate=1):
        classes_eds2.Node.__init__(self)
        self.rate=rate
        self.threshold=threshold
        self.hill=hill_coeff
        self.dephosphorylation=dephospho_rate
        self.label='Phosphorylation'
        self.input=['Kinase','Phosphorylable']
        self.output=['Kinase','Phospho']

    def __str__(self):
        return "{0.id} Phosphorylation: rate = {0.rate:.2f}, thres. = {0.threshold:.2f}, H_coeff. = {0.hill:.2f}, dissoc. = {0.dephosphorylation:.2f}".format(self)

    def outputs_to_delete(self,net):
        """Return the phosphorylated species to delete when deleting a Phosphorylation"""
        listIn = net.graph.predecessors(self)
        listOut = net.graph.successors(self)
        return [out for out in listOut if out not in listIn] #to avoid the kinase

########## Attributes attached to Network for Phosphorylations/Dephosphorylations ##########

def check_existing_Phosphorylation(self,signature):
    """check if a particular phosphorylation exists in the network

    Args:
        signature (list): The signature of the phospho in the form [Kinase,Input]

    Return:
        bool: if this phosphorylation exist
    """
    if 'Phosphorylation' in self.list_types:#goes through the list of interactions
        for inter in self.list_types['Phosphorylation']:
            [catalyst,listIn,listOut] = self.catal_data(inter)
            if (catalyst==signature[0]) and (listIn[0]==signature[1]):
                return True
    return False

def number_Phosphorylation(self):
    """Return the number of possible Phosphorylations"""
    nK=self.number_nodes('Kinase')
    nS=self.number_nodes('Phosphorylable')
    n_P=self.number_nodes('Phosphorylation')#number of existing Phosphorylations
    return nK*nS - n_P

def new_Phosphorylation(self,kinase,species,rate,threshold,hill,dephospho):
    """Create a new Phosphorylation, its associated product and add them to the network.

    Args:
        kinase (Species): -
        species (Species): -
        rate (float): the association rate
        threshold (float): the Michaelis-Menten constant
        hill (float): the hill coefficient of the reaction
        dephospho (float): the dephosphorylation rate of the product

    Return:
        list: of the form [Phosphorylation,phosphorylated_Species]
        or None if an error occured
    """
    phospho=Phosphorylation(rate,threshold,hill,dephospho)#creates the interaction
    species_P=copy.deepcopy(species)# the phosphorylated species has the same properties
    species_P.clean_type('Input')#Remove Input, output types from the copied species
    species_P.clean_type('Output')#Remove Input, output types from the copied species
    if species.isinstance('Phospho'):#allow at most two phosphorylations : if species is already phosphorylated, the product is no longer phosphorylable
        species_P.clean_type('Phosphorylable')#Phosphorylable species can not be more phosphorylated
    species_P.clean_type('Phospho')
    species_P.types.append('Phospho')#  phosphorylates
    species_P.n_phospho=1
    if phospho.check_grammar([kinase,species],[kinase,species_P]):
        self.add_Node(phospho)
        self.graph.add_edge(kinase,phospho)
        self.graph.add_edge(species,phospho)
        self.add_Node(species_P)
        self.graph.add_edge(phospho,species_P)
        self.graph.add_edge(phospho,kinase)
        return [species_P,phospho]
    else:
        print("Error in grammar : new Phosphorylation")
        return None

# Add the corepromoter functions to the Network class
setattr(classes_eds2.Network,'check_existing_Phosphorylation',check_existing_Phosphorylation)
setattr(classes_eds2.Network,'number_Phosphorylation',number_Phosphorylation)
setattr(classes_eds2.Network,'new_Phosphorylation',new_Phosphorylation)

########## Attributes attached to Mutable_Network for Phosphorylations/Dephosphorylations ##########

def new_random_Phosphorylation(self, kinase, species):
    """Creates a Phosphorylation of species by kinase with random parameters

    Args:
        kinase (Species): the kinase came first
        species (Species): -

    Return:
        list: of the form [Phosphorylation,phosphorylated_Species]
        or None if an error occured
    """
    r = mutation.sample_dictionary_ranges('Phosphorylation.rate',self.Random)
    t = mutation.sample_dictionary_ranges('Phosphorylation.threshold',self.Random)
    h = mutation.sample_dictionary_ranges('Phosphorylation.hill',self.Random)
    d = mutation.sample_dictionary_ranges('Phosphorylation.dephosphorylation',self.Random)
    return self.new_Phosphorylation(kinase,species,r,t,h,d)

def random_Phosphorylation(self):
    """Creates a new Phosphorylation among all possibles


    Return:
        list: of the form [Phosphorylation,phosphorylated_Species]
        or None if an error occured
    """
    if 'Kinase' in self.list_types and 'Phosphorylable' in self.list_types:
        #List all possible phosphorylations
        possible_Phospho=[(kinase,species) for kinase in self.list_types['Kinase']
                                           for species in self.list_types['Phosphorylable']
                                           if not self.check_existing_Phosphorylation([kinase,species])]
        n_pP=len(possible_Phospho)
        if not (n_pP==self.number_Phosphorylation()):
            print("Potential Bug : Inconsistency in Computation of number of Phosphorylations")
            print(n_pP,self.number_Phosphorylation())
            print(possible_Phospho)
            print(self.list_types['Kinase'])
            print(self.list_types['Phosphorylable'])
            p=self.list_types['Phosphorylable'][0]
            p.print_node()
            print(self.list_types['Phosphorylation'])
        if (n_pP==0):
            print("In random_Phosphorylation : No other posible Phosphorylationss")
            return None
        else:
            [K,S]=possible_Phospho[int(self.Random.random()*n_pP)]
            return self.new_random_Phosphorylation(K,S)
    else:
        print("Error in random_Phosphorylation (try to create Phosphorylation from non existing pieces)")
        return None

# Add the corepromoter functions to the Mutable_Network class
setattr(mutation.Mutable_Network,'random_Phosphorylation',random_Phosphorylation)
setattr(mutation.Mutable_Network,'new_random_Phosphorylation',new_random_Phosphorylation)

########## Integration C Tools ##########

def Phospho_deriv_inC(net):
    """gives the string corresponding to Phosphorylation for integration

    Return:
        str: a single string for all Phosphorylations in the network
    """
    func="\n/**************Phosphorylation*****************/\n float total;\n"
    if ('Phosphorylation' in net.list_types):
        dict_kinase={}#dictionnary to keep track of multiple phosphorylations by same kinase; each entry is a list; 1st term to put in the denominator, 2nd term is the list of all equations  where the kinase plays a role
        for node in net.list_types['Kinase']:
            dict_kinase[node]=["1",""]
        for reaction in net.list_types['Phosphorylation']:
            [kinase,species,species_P]=net.catal_data(reaction)
            species=species[0]
            species_P=species_P[0]
            term="POW(%s/%f,%f)"%(species.id , reaction.threshold , reaction.hill) #computes the numerator corresponding to this specific phophorylation
            dict_kinase[kinase][0]+="+"+term#adds to denominator

            prate="%f*%s*(%s/total)"%(reaction.rate , kinase.id , term) #writes the rate
            dephosphorate="%f*%s"%(reaction.dephosphorylation , species_P.id)
            dict_kinase[kinase][1]=dict_kinase[kinase][1]+"\t \t/*Phosphorylation*/\n"
            dict_kinase[kinase][1]=dict_kinase[kinase][1]+deriv2.compute_leap([species.id],[species_P.id],prate)
            dict_kinase[kinase][1]=dict_kinase[kinase][1]+"\t \t /*Dehosphorylation*/\n"
            dict_kinase[kinase][1]=dict_kinase[kinase][1]+deriv2.compute_leap([species_P.id],[species.id],dephosphorate)
            dict_kinase[kinase][1]=dict_kinase[kinase][1]+"\n"

        for kinase in net.list_types['Kinase']:
            if not (dict_kinase[kinase][0]=="1"):
                func=func+"\ntotal="+dict_kinase[kinase][0]+";\n"+dict_kinase[kinase][1]#writes the rates for each kinase
    return func

#update deriv2
deriv2.interactions_deriv_inC["Phospho"] = Phospho_deriv_inC
