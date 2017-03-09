""" This program contains or calls all the functions used for mutations which are defined to act on a subclass of network
called Mutable_Network defined here
For mutation/removal, effective mutation rate will be the reference mutation rate times the number of instances of the considered Type.
Attributes:
    C,L,T (float):
    dictionary_mutation (dict): referenced mutation rates and associated command as key
    dictionary_ranges (dict): referenced parameters that can change and their ranges
    list_create (list): list of Nodes subject to creation
    list_mutate (list): list of Nodes subject to mutation
    list_remove (list): list of Nodes subject to removal
    list_types_output (list): list of the possible types for the output
"""
print('Execute mutation.py')

#from phievo.initialization_code import *
import random
import copy
import collections
from math import log,exp
from . import classes_eds2
from .deriv2 import compile_and_integrate
from phievo.initialization_code import multiple_phospho
"""
    The dictionary_mutation[] has ; mutation will exec the key.
    
  The rates referring to 'mutate' change parameters.   See initialization file for samples
    The  are in dictionary_ranges, whose keys are of form Node.attribute, see rand_modify().
    The lists below are names of Nodes that are subject to indicated operations, needed in various routines below
"""

##########################
### DEFAULT PARAMETERS ###
##########################

T=1.0 #typical time scale
C=1.0 #typical concentration
L=1.0 #typical size for diffusion

dictionary_ranges = {}
dictionary_ranges['Species.degradation']=1.0/T
dictionary_ranges['Species.diffusion']=L   # for ligands diffusion
dictionary_ranges['TModule.rate']=C/T
dictionary_ranges['TModule.basal']=0.0

dictionary_mutation={}

list_types_output=['Species']
list_mutate=[]
list_remove=[]
list_create=[]

#########################
### Routine functions ###
#########################

def build_lists(mutation_dict):
    """Construct the index of Species types subject to various operation
    
    Args:
        mutation_dict (dict): the dictionary listing the various operation (typically inits.dictionary_mutations)
    
    Return:
        None: modify list_create - list_mutate - list_remove as global variable
    """
    global list_mutate, list_remove, list_create
    # when a term is found, we used split on ', the type should be the second term
    for index in mutation_dict:
        if 'mutate_Node' in index:
            list_mutate.append(index.rsplit('\'')[1])
        if 'remove_Interaction' in index:
            list_remove.append(index.rsplit('\'')[1])
        if 'random_Interaction' in index:
            list_create.append(index.rsplit('\'')[1])

def sample_dictionary_ranges(key,random_generator):
    """Draw a random value for a parameter of type key
    
    Look on dictionary_range, if the attribute to key is:
    a real number, a list or tuple of two reals
    defining min-max of range, and sample accordingly.
    
    Args:
        key (str): the type of parameter you want
        random_generator: a random number generator (.random() called here)

    Return:
        float: a random value
        or int if key is CorePromoter.delay
        or None if an error occured
    """
    if key not in dictionary_ranges:
        print('in sample_dictionary_ranges(), invalid key=', key,'for dictionary_ranges')
        return None
    
    interval = dictionary_ranges[key]
    rrand = random_generator.random()
    if isinstance(interval,collections.Sequence): #when interval is a list or tuple
        dice = interval[0] + rrand*(interval[1] - interval[0])
    else: #when it is a number
        dice = rrand * interval
    
    if key == 'CorePromoter.delay':
        return int(dice)
    else:
        return dice

def random_parameters(Type,random_generator):
    """Create a set of new random parameters for a Species instance of type Type
    
    This used only for initialization and adds attributes to various types.
    Some of which may not be mutable later
    
    Args:
        Type (str): a species type
        random_generator: a random number generator (.random() called here)
    
    Return:
        list: a list of random parameters that can create a new Species
        or None if an error occured
    """
    if not Type in classes_eds2.Species.Tags_Species: #look directly at the class attribute 
        print("Try to create a  not allowed random species of type "+Type)
        return None

    list=[['Degradable', sample_dictionary_ranges('Species.degradation',random_generator) ],['Phosphorylable']]
    
    if multiple_phospho: list.append(['Phospho',0])

    # We start with Ligand because it is the trickiest one
    if Type=='Ligand':
        list.append(['Ligand'])
        external=int(2*random_generator.random())
        if (external==1):
            parameter=sample_dictionary_ranges('Species.diffusion',random_generator)
            list.append(['Diffusible',parameter])
    elif Type=='TF':
        parameter=sample_dictionary_ranges('Species.diffusion',random_generator)
        list.append(['Diffusible',parameter])
        list.append(['TF',int(2*random_generator.random())])
        list.append(['Complexable'])
    elif Type=='Kinase' or Type=='Receptor':
        list.append([Type])
        list.append(['Complexable'])
    return list

def rand_modify(self,random_generator):
    """modify every parameters of the node self
    
    This subroutine is then export to the :class:`Networks.classes_eds2.Node` class and used as a method
    Called the sample_dictionary_ranges subroutine when needed
    
    Args:
        self (:class:`Networks.classes_eds2.Node`): the node you want to modify
        random_generator: a random number generator (.random() called here)

    Return:
        None: in place modification
    """
    name=self.__class__.__name__ #take the name of the class
    if not self.mutable: return None #do nothing
        
    for k in self.__dict__: #take all the attributes of a class
        entry=name+"."+k  #builds the key to check
        if entry in dictionary_ranges:
            if 'relative_variation' in dictionary_ranges:
                # in that case, we suppose that each parameter is basically max_value*exp(-Energy) where max_value is the maximum authorized
                # value of the parameter; we then modify randomly the energy, a bias might be included to increase probability of "killing" interactions
                previous = getattr(self,k)

                # takes the range of variation
                modif=dictionary_ranges['relative_variation']
                # takes the upper boundary of the authorized interval
                interval=dictionary_ranges[entry]
                if isinstance(interval,float) or isinstance(interval,int):
                    interval=[0,dictionary_ranges[entry]]
                min_param=interval[0]
                if (min_param==0):
                    min_param=0.01 #minimum value to have a maximum energy, otherwise the range of energies is infinite
                E_I=log(interval[1]/min_param) #range of possible energies
                if (previous>0):
                   current_E=log(interval[1]/previous) #current energy
                else:
                    current_E=E_I
                if 'bias' in dictionary_ranges:
                    bias=dictionary_ranges['bias']
                else:
                    bias=0
                if (k=='threshold'):
                    bias=-bias#one wants to bias thresholds to increase them
                if (k=='hill') or (k=='diffusion') or (k=='delay'):
                    bias=0#no bias for these parameters

                #computes randomly the new energy; the bigger the energy, the smaller the parameter
                new_energy=current_E+(2*(random_generator.random()-0.5)+bias)*E_I*modif

                next=interval[1]*exp(-new_energy)
                #print entry,previous,interval[1],next
                if isinstance(interval[1],int):
                    next=int(next)
                if (interval[1]==0):
                    next=0

                setattr(self,k,next)

                if (next>interval[1]): setattr(self,k,interval[1])
                if (next<interval[0]): setattr(self,k,interval[0])
            else:
                setattr(self,k,sample_dictionary_ranges(entry,random_generator))
        else:
            pass

# update classes_eds2.Node
setattr(classes_eds2.Node,'rand_modify',rand_modify)

########################################
### Mutable_Network Class Definition ###
########################################

class Mutable_Network(classes_eds2.Network):
    """Expand the Network class with all functions related to mutation
    
    the random_Type() routines are the only ones called by evolution to sample all possible
    links on graph that can give rise to given interaction Type and then choses one.
    The assignment of random interaction parametes and types of output, packaged in
    separate routines new_random_Type(), that can be used independently to generate
    specific topologies with random parameters
    Attributes:
        fitness (float): the fitness of the Network, None by default (worst than everyother number)
        dlt_fitness (float): the change of fitness at the last generation
        data_evolution (list): keep various information such as fitness variance, average…
        data_next_mutation (list): field to keep the data on the next mutation
        Random (Random): defines the local random generator number
    
    Main functions:
    """
    def __init__(self,generator = random.Random()):
        """Constructor of the Mutable_Network class
        
        Args:
            generator (Random): the local random generator number
        """
        classes_eds2.Network.__init__(self)
        self.fitness = None
        self.dlt_fitness = 0
        self.data_evolution = []
        self.data_next_mutation = [0,""]
        self.Random = generator

    def compute_Cseed(self):
        """Return a random integer to determine the integrator seed"""
        return int(self.Random.random()*100000)

########## Tools to add and remove Nodes ##########

    def random_Species(self, Type='Species'):
        """Create a new random species instance of a given type
        
        Args:
            Type (str): the desired type of the new Species instance 
        
        Return:
            Species: note that it is automatically added to the network
        """
        return self.new_Species(random_parameters(Type,self.Random))

    def random_Interaction(self,Interaction_Type):
        """ create a new (and unique) interaction
        
        Args:
            Interaction_Type (str): the type of interaction you want
        
        Return:
            None
        """
        try:
            getattr(self,'random_'+Interaction_Type)()
        except Exception:
            display_error("Error when creating randomize Interaction "+Interaction_Type)

    def remove_Interaction(self,Type):
        """Randomly removes a Node of a given Type
        
        Args:
            Type (str): the type you want to remove (e.g. 'Interaction', :class:`Networks.TFHill.TFHill`, ...)
        
        Return:
            bool: indicate if something is effectively removed
        """
        self.write_id()
        if Type in self.list_types:
            condemn = self.Random.choice(self.list_types[Type])
            if not (isinstance(condemn,classes_eds2.Interaction)):
                print("Try to remove something else than an interaction :"+Type)
                return False
            Bool=self.remove_Node(condemn)
            self.clean_Nodes()
            return Bool
        else:
            print("Nothing to remove")
            return False

########## Tools to manage the Output Tags ##########

    def random_remove_output(self):
        """Removes at random an output tag on some species
        
        Outputs are always index 0,1,2...; not possible to have 0,1,3 for instance
        """
        if 'Output' in self.list_types:
            to_remove=self.Random.choice(self.list_types['Output'])
            to_remove.clean_type('Output') #cleans previous output
            self.list_types['Output'].remove(to_remove)
            #recomputes the index of outputs
            for index,species in enumerate(self.list_types['Output']):
                species.n_put = index

    def random_add_output(self):
        """Randomly adds an output tag to a random species"""
        noutput = len(self.list_types.get('Output',[]))
        Type = self.Random.choice(list_types_output)
        
        test_output = lambda S: not S.isinstance('Input') and not S.isinstance('Output')
        list_possible_outputs = [species for species in self.list_types.get(Type,[]) if test_output(species)]

        if list_possible_outputs:
            species = self.Random.choice(list_possible_outputs)
            species.add_type(['Output',noutput])
            return True
        else:
            return False

    def random_change_output(self):
        """Function that changes one output by adding then removing a TAG output"""
        if self.random_add_output():
            self.random_remove_output()
            self.write_id()

########## Duplication Tools ##########

    def random_duplicate(self):
        """Routine to duplicate gene and its interactions
        
        Currently the classes_eds2.duplicate_* only implemented for TF & PPI interactions
        If duplicating an output gene, add a new output tag to duplicated species, irrespective of other dictionary_mutation['*output*'] values in initialization.
        
        Return:
            bool: indicate if a duplication has been finally done
        """
        possible_duplicate = [self.graph.successors(interaction)[0] for interaction in self.list_types['CorePromoter']]
        possible_duplicate.sort(key = classes_eds2.compare_node) #to be deterministic
        
        if possible_duplicate:
            species= self.Random.choice(possible_duplicate)
            [D_module,D_promoter,D_species] = self.duplicate_species_and_interactions(species)
            if species.isinstance('Output'):
                D_species.add_type(['Output',len(self.list_types['Output'])])
            self.write_id()
            return True
        else:
            return False

########## Global Mutation Tools ##########

    def mutate_Node(self,Type):
        """randomly selects then mutates a Node of a given Type
        
        Args:
            Type (str): the Type to mutate (e.g. ``Species``: :class:`Networks.classes_eds2.Species`, ``TFHill``: :class:`Networks.TFHill.TFHill`, ``Node``: :class:`Networks.classes_eds2.Node`...)
        
        Return:
            bool: if something is mutated
        """
        self.write_id()
        possible_mutated = self.list_types.get(Type,[])
        if possible_mutated:
            self.Random.choice(possible_mutated).rand_modify(self.Random)
            return True
        else:
            print("Nothing to mutate "+Type)
            return False

    def build_mutations(self):
        """builds a dictionary with relative mutation rates for a specific network
        
        This method is based on dictionary_mutation
        
        
        Returns:
            dict: with the rate of each events for the network
        """
        self.write_id()
        dictionary=copy.deepcopy(dictionary_mutation) #takes the predefined dictionary
        
        def update_dict(self,dictionary,key,name):
            """Subroutine for build_mutations
            Multiply dictionary[key] with the length of self.list_types[name]
            if non zero or delete the key otherwise
            """
            if (name in self.list_types):
                dictionary[key] *= len(self.list_types[name]) #update the corresponding mutation rate
            else: #delete corresponding key if does not exist
                del dictionary[key]
        
        for name in list_remove:
            key='remove_Interaction(\''+name+'\')'
            update_dict(self,dictionary,key,name)

        for name in list_mutate:
            key='mutate_Node(\''+name+'\')'
            update_dict(self,dictionary,key,name)
        
        update_dict(self,dictionary,'random_remove_output()','Output')
        update_dict(self,dictionary,'random_change_output()','Output')
        
        for name in list_create:
            key ='random_Interaction(\''+name+'\')' #find the name of the command
            l = getattr(self,'number_'+name)()
            dictionary[key] *= l #compute the corresponding mutation rate
            if (l==0): del dictionary[key]
        
        key='random_duplicate()'
        if key in dictionary:
            if 'CorePromoter' in self.list_types:
                #l=len(self.list_types['CorePromoter'])
                dictionary[key]*=1 #compute the corresponding mutation rate
            else:
                del dictionary[key]
        
        return dictionary

    def compute_next_mutation(self):
        """determine the time and type of next mutation for the gillespie algo.
        
        Return:
            float: time to next mutation
            
        given a network, computes the time of the next mutation and the command to execute to perform the mutation
        for the gillispie algorithm
        """
        dictionary=self.build_mutations()
        list_commands=list(dictionary.keys())
        list_commands.sort() #sort to keep something deterministic
        a0=0
        for keys in list_commands:
            #print keys,dictionary[keys]
            a0=a0+dictionary[keys] #compute sum of the rates
        #print a0
        r=self.Random.random()
        tau=-1.0/a0*log(r) #computes the next time of mutation

        mutation=self.Random.random()*a0 #quantity to compute the next mutation
        a1=0
        index_keys=0
        while mutation>a1:
            a1=a1+dictionary[list_commands[index_keys]]
            index_keys+=1# we go out of the loop when mutation<=a1, it means that we have gone too far of one key
        return [tau,list_commands[index_keys-1]]

################## mutation/integration tools for one network #####################

    def mutate_and_integrate(self,prmt,nnetwork,tgeneration,mutation=True):
        """ function to mutate, integrate and update the fitness
        
        Note that compile_and_integrate is defined in Networks/deriv2.py
        
        Args:
            prmt (dict): 
            nnetwork (int): an id for the C-file
            tgeneration (float): the time before the next gen.
            mutation (bool): if False, no mutation will be made
        
        Returns:
            (list): [n_mutations,nnetwork,self,result] where:
                - n_mutations (int): the numbre of mutation performed
                - nnetwork (int): same as args
                - self (Network): the Mutable_Network object itself
                - result (list): output of treatment_fitness (see compile_and_integrate)
        """
        n_mutations,age = 0,0
        if mutation:
            while True:
                tau,next_mutation = self.compute_next_mutation()
                age += tau
                if age > tgeneration: break #exit the loop when enough time has passed
                exec("self."+next_mutation)
                n_mutations+=1
            age -= tgeneration
            self.data_next_mutation[0:2] = [age,next_mutation]  #keeps track of the time and type of the next mutation
        
        self.Cseed = self.compute_Cseed()
        result = compile_and_integrate(self,prmt,nnetwork,0,self.Cseed)
        return [n_mutations,nnetwork,self,result]




