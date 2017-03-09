
# importing all relevant modules.
import classes_eds2 
import mutation 
import copy        
import random
import config
exec('import '+config.name_deriv2+' as deriv2')

# parameters needed to perform evolution.
mutation.dictionary_ranges['Simple_Phosphorylation.rate_p']=1.0/mutation.T/mutation.C
mutation.dictionary_ranges['Simple_Phosphorylation.rate_d']=1.0/mutation.T/mutation.C

# definition of class phosphorylation: sub-class of interaction class (see classes_eds2).
# contains the parameter quantifying the interaction as well as the input and output of the reaction. 
# these characteristics define the interaction completely. they are "free standing" initially in that 
# they are not appended to the graph per se. see new_Phosphorylation below to see how to actually 
# add the interaction to the network.

class Simple_Phosphorylation(classes_eds2.Interaction):
    """ Phosphorylation, init : rate of phosphorylation. """

    def __init__(self,p=0,d=0):
        classes_eds2.Node.__init__(self)
        self.rate_p=p
        self.rate_d=d
        self.label='Simple_Phosphorylation'
        self.input=['Kinase','Phosphatase','Phosphorylable']		# these are required to check if the required species within the network are present to start with.
        self.output=['Kinase','Phosphatase','Phospho']


    def outputs_to_delete(self,net):
        """ Returns the phosphorylated species to delete when deleting a Phosphorylation"""
        return net.find_Simple_phospho_to_remove(self)



####### Attributes attached to Network for  Phosphorylations/Dephosphorylations############

# determines if two species (kinase and protein) are related by phosphorylation interaction (first 
# species should be the kinase).
def check_existing_Simple_Phosphorylation(self,list):
        """ special function to check if a simple phosphorylation exists : order in the list should be Kinase first, second species and third phosphatase."""
        if 'Simple_Phosphorylation' in self.list_types:	#goes through the list of interactions
            for inter in self.list_types['Simple_Phosphorylation']:
                [catalyst,listIn,listOut]=self.catal_data(inter)
                if (catalyst[0]==list[0]) and (listIn[0]==list[1]) and (catalyst[1]==list[3]):
                    return True
                if (catalyst[1]==list[0]) and (listIn[0]==list[1]) and (catalyst[0]==list[3]):
                    return True
        return False

def number_Simple_Phosphorylation(self):
        """ Computes the number of possible Simple_Phosphorylations"""
        if 'Simple_Phosphorylation' in self.list_types:
            return len(self.list_types['Simple_Phosphorylation'])
        else:
            return 0
        
        
def find_Simple_phospho_to_remove(self,interaction):
        """Finds products to remove when deleting a Simple_Phosphorylation interaction."""

        listIn=self.graph.predecessors(interaction)
        listOut=self.graph.successors(interaction)
        Bool=False
        for x in listIn:
            if x in listOut:
                catalyst=x
                Bool=True
                #there is this Bool business to intercept the case where the catalyst has been already removed by clean_Nodes, we then just remove the output list
        if Bool:
            listOut.remove(catalyst)
            return listOut

        else:
            return []
        

def check_(self,species,Type):
    list_successors = self.graph.successors(species)
    for i in range(len(list_successors)):
        if list_successors[i].isinstance(Type):
            return True
    return False

# important function that adds the interaction to the actual graph. 
def new_Simple_Phosphorylation(self,kinase,phosphatase,species,p,d):
        # we allow for any number of phosphorylations and keep track of the order of phosphorylation through attribute n_phospho of tag 'Phospho'.

        phospho=Simple_Phosphorylation(p,d) #creates the interaction (see function above)
        species_P=copy.deepcopy(species) # the phosphorylated species has the same properties (need to use deep copy and not just copy: we do not want the pointer but to create an altogether completely new species)
        species_P.clean_type('Input')#Remove Input, output types from the copied species
        species_P.clean_type('Output')#Remove Input, output types from the copied species
        species.add_type('Phospho')
          
        if phospho.check_grammar([kinase,phosphatase,species],[kinase,phosphatase,species_P]): 		# verifying that the "grammar" is correct.
            self.add_Node(phospho)				# adding the various nodes and edges corresponding to the interaction.
            self.graph.add_edge(kinase,phospho)
            self.graph.add_edge(species,phospho)
            self.add_Node(species_P)
            self.graph.add_edge(phospho,species_P)
            self.graph.add_edge(phospho,kinase)
            self.graph.add_edge(phosphatase,phospho)
            self.graph.add_edge(phospho,phosphatase)

            return [species_P,phospho]
    
        else:
            print("Error in grammar : new Simple_Phosphorylation")
            return None
             


setattr(classes_eds2.Network,'check_existing_Simple_Phosphorylation',check_existing_Simple_Phosphorylation)
setattr(classes_eds2.Network,'number_Simple_Phosphorylation',number_Simple_Phosphorylation)
setattr(classes_eds2.Network,'find_Simple_phospho_to_remove',find_Simple_phospho_to_remove)
setattr(classes_eds2.Network,'new_Simple_Phosphorylation',new_Simple_Phosphorylation)



####### Attributes attached to Mutable_Network for  Phosphorylations/Dephosphorylations############


def new_random_Simple_Phosphorylation(self, K, P, S):
        """ Create new random Simple_Phosphorylation interaction and return interaction and product. Note that the inputs must be given (as opposed to random_Simple_Phosphorylation)."""
        r = mutation.sample_dictionary_ranges('Simple_Phosphorylation.rate_p',self.Random)
        d = mutation.sample_dictionary_ranges('Simple_Phosphorylation.rate_d',self.Random)
        [SP,Phospho]=self.new_Simple_Phosphorylation(K,P,S,r,d)
        return [SP,Phospho]



def random_Simple_Phosphorylation(self):
        """Create new random Simple_Phosphorylations (includes directly dephosphorylation by a phosphatase)."""
        
        if 'Kinase' in self.list_types and 'Phosphorylable' in self.list_types and 'Phosphatase' in self.list_types:
                    
            list_possible_Phospho=[]
            nK=len(self.list_types['Kinase'])
            nS=len(self.list_types['Phosphorylable'])
            nP=len(self.list_types['Phosphatase'])
            
            for iS in range(nS):
                S=self.list_types['Phosphorylable'][iS]
                Out_S=self.graph.successors(S)
                bool_already_phospho = False
                
                for iOut_S in range(len(Out_S)):
                    if not Out_S[i].isinstance('Simple_Phosphorylation'):
                        [catalyst,species,species_P]=net.catal_data(Out_S[i])
                        if species==S:
                            bool_already_phospho = True
                if not bool_already_phospho:
                    for iK in range(nK):
                        K=self.list_types['Kinase'][iK]
                        for iP in range(nP):
                            P=self.list_types['Phosphatase'][iP]
                            list_possible_Phospho.append([K,P,S])
                            
            n_pP=len(list_possible_Phospho)
            if (n_pP==0):
                print("In random_Simple_Phosphorylation : No other posible Simple_Phosphorylationss")
                return None
            else:
                [K,P,S]=list_possible_Phospho[int(self.Random.random()*n_pP)]
                [SP,Phospho]=self.new_random_Simple_Phosphorylation(K,P,S)
                return [SP,Phospho]        
        else:
            print("Error in random_Phosphorylation (try to create Simple_Phosphorylation from non exsiting pieces)")
            return None   

    
setattr(mutation.Mutable_Network,'random_Simple_Phosphorylation',random_Simple_Phosphorylation)
setattr(mutation.Mutable_Network,'new_random_Simple_Phosphorylation',new_random_Simple_Phosphorylation)


#############Integration C Tools #################



def Simple_Phosphorylation_deriv_inC(net):
    func="\n/**** Simple_Phosphorylation and Dephosphorylation ****/\n" 
 
    list_phospho = []
    if ('Simple_Phosphorylation' in net.list_types):
        for reaction in net.list_types['Simple_Phosphorylation']:
		    [catalyst,species,species_P]=net.catal_data(reaction)
		    species = species[0]
		    species_P = species_P[0]
		    if catalyst[0].isinstance('Kinase'):
		        kinase = catalyst[0]
		        phosphatase = catalyst[1]
		    else: 
		        phosphatase = catalyst[0]
		        kinase = catalyst[1]
		        
		    if (kinase.common==0) and (phosphatase.common==0) and (species.common==1):
		        rate_p = "N_cell_*%f*"%reaction.rate_p+kinase.id+"*"+species.id
		        rate_d = "N_cell_*%f*"%reaction.rate_d+phosphatase.id+"*"+species_P.id
		        func = func+deriv2.compute_leap([species.id],[species_P.id],rate_p)
		        func = func+deriv2.compute_leap([species_P.id],[species.id],rate_d)
		    elif (kinase.common==0) and (phosphatase.common==1) and (species.common==1):
		        rate_p = "N_cell_*%f*"%reaction.rate_p+kinase.id+"*"+species.id
		        rate_d = "%f*"%reaction.rate_d+phosphatase.id+"*"+species_P.id
		        func = func+deriv2.compute_leap([species.id],[species_P.id],rate_p)
		        func = func+deriv2.compute_leap([species_P.id],[species.id],rate_d)
		    elif (kinase.common==1) and (phosphatase.common==0) and (species.common==1):
		        rate_p = "%f*"%reaction.rate_p+kinase.id+"*"+species.id
		        rate_d = "%f*"%reaction.rate_d+phosphatase.id+"*"+species_P.id
		        func = func+deriv2.compute_leap([species.id],[species_P.id],rate_p)
		        func = func+deriv2.compute_leap([species_P.id],[species.id],rate_d)	
		    else:
		        rate_p = "%f*"%reaction.rate_p+kinase.id+"*"+species.id
		        rate_d = "%f*"%reaction.rate_d+phosphatase.id+"*"+species_P.id
		        func = func+deriv2.compute_leap([species.id],[species_P.id],rate_p)
		        func = func+deriv2.compute_leap([species_P.id],[species.id],rate_d)	

    return func

deriv2.Simple_Phosphorylation_deriv_inC=Simple_Phosphorylation_deriv_inC