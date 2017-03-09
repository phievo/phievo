
# importing all relevant modules.
import classes_eds2 
import mutation 
import copy  
import config      
exec('import '+config.name_deriv2+' as deriv2')
import random


# parameters needed to perform evolution.
mutation.dictionary_ranges['Linear_Production.rate']=1.0/mutation.T


class Linear_Production(classes_eds2.Interaction):
    """ Phosphorylation, init : rate of phosphorylation. """

    def __init__(self,p=0,d=0):
        classes_eds2.Node.__init__(self)
        self.rate=r
        self.label='Linear_Production'
        self.input=['Linear_Producer']		# these are required to check if the required species within the network are present to start with.
        self.output=['Species','Linear_Producer']


    def outputs_to_delete(self,net):
        """ Returns the phosphorylated species to delete when deleting a Phosphorylation"""
        return net.find_Linearly_Produced_to_remove(self)



####### Attributes attached to Network for  Phosphorylations/Dephosphorylations############

# determines if two species (kinase and protein) are related by phosphorylation interaction (first 
# species should be the kinase).
def check_existing_Linear_Production(self,list):
        """ special function to check if a Linear_Production exists : order in the list should be Linear_Producer first, second species."""
        if 'Linear_Production' in self.list_types:	#goes through the list of interactions
            for inter in self.list_types['Linear_Production']:
                [catalyst,listIn,listOut]=self.catal_data(inter)
                if (catalyst==list[0]) and (listOut[0]==list[1]):
                    return True
        return False

def number_Linear_Production(self):
        """ Computes the number of possible Simple_Phosphorylations"""
        if 'Linear_Production' in self.list_types:
            return len(self.list_types['Linear_Production'])
        else:
            return 0
        
        
def find_Linearly_Produced_to_remove(self,interaction):
        """Finds products to remove when deleting a Linear_Production interaction."""

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
def new_Linear_Production(self,linear_producer,species,r):
        # we allow for any number of phosphorylations and keep track of the order of phosphorylation through attribute n_phospho of tag 'Phospho'.

        lin_prod=Linear_Production(r) #creates the interaction (see function above)
          
        if phospho.check_grammar([linear_producer],[linear_producer,species]): 		# verifying that the "grammar" is correct.
            self.add_Node(lin_prod)				# adding the various nodes and edges corresponding to the interaction.
            self.graph.add_edge(linear_producer,lin_prod)
            self.graph.add_edge(lin_prod,linear_producer)
            self.graph.add_edge(lin_prod,species)

            return lin_prod
    
        else:
            print("Error in grammar : new_Linear_Production")
            return None
             


setattr(classes_eds2.Network,'check_existing_Linear_Production',check_existing_Linear_Production)
setattr(classes_eds2.Network,'number_Linear_Production',number_Linear_Production)
setattr(classes_eds2.Network,'find_Linearly_Produced_to_remove',find_Linearly_Produced_to_remove)
setattr(classes_eds2.Network,'new_Linear_Production',new_Linear_Production)



####### Attributes attached to Mutable_Network for  Phosphorylations/Dephosphorylations############


def new_random_Linear_Production(self, L, S):
        """ Create new random Simple_Phosphorylation interaction and return interaction and product. Note that the inputs must be given (as opposed to random_Simple_Phosphorylation)."""
        r = mutation.sample_dictionary_ranges('Linear_Production.rate',self.Random)
        Lin_Prod=self.new_Linear_Production(L,S,r)
        return Lin_Prod



def random_Linear_Production(self):
        """Create new random Simple_Phosphorylations from list of possible kinase substrates."""
        
        if 'Linear_Producer' in self.list_types and 'Species' in self.list_types:
                    
            list_possible_Lin_Prod=[]
            nL=len(self.list_types['Linear_Producer'])
            nS=len(self.list_types['Species'])
            
            for iL in range(nL):
                L=self.list_types['Linear_Producer'][iL]
                                
                for iS in range(nS):
                    S=self.list_types['Species']
                    if not (S==L):
                        list_possible_Lin_Prod.append([L,S])
                
            n_LP=len(list_possible_Lin_Prod)
            if (n_LP==0):
                print("In random_Linear_Production : No other posible Linear_Production")
                return None
            else:
                [L,S]=list_possible_Lin_Prod[int(self.Random.random()*n_LP)]
                Lin_Prod=self.new_random_Linear_Production(L,S)
                return Lin_Prod        
        else:
            print("Error in random_Linear_Production (try to create Linear_Production from non exsiting pieces)")
            return None   

    
setattr(mutation.Mutable_Network,'random_Linear_Production',random_Linear_Production)
setattr(mutation.Mutable_Network,'new_random_Linear_Production',new_random_Linear_Production)


#############Integration C Tools #################



def Linear_Production_deriv_inC(net):
    func="\n/**** Linear_Production ****/\n" 
 
    list_Linear_Production = []
    if ('Linear_Production' in net.list_types):
        for reaction in net.list_types['Linear_Production']:
            lin_producer = net.graph.predecessor(reaction)
            succ = net.graph.successor(reaction)
            if succ[0] == lin_producer:
                species = succ[1]
            else:
                species = succ[0]
            
            if (species.common==1) and (lin_producer.common==0):
                rate = "N_cell_*%f*"%reaction.rate+lin_producer.id
                func = func+deriv2_pMHC.compute_leap([species.id],[],rate)
            else:
                rate = "%f*"%reaction.rate+lin_producer.id
                func = func+deriv2_pMHC.compute_leap([species.id],[],rate)	
    return func

deriv2.Linear_Production_deriv_inC=Linear_Production_deriv_inC