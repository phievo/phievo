# importing all relevant modules.
import phievo.Networks.classes_eds2 as classes_eds2
import phievo.Networks.mutation as mutation
import copy        
#import config
from phievo.initialization_code import *
exec('import '+name_deriv2+' as deriv2')
import random

# overwritting the Phospho tag from classes_eds2 to keep track of the degree of phosphorylation.
# Tags_Species['Phospho'] = ['n_Phosho']
# does not seem to work.

# parameters needed to perform evolution.
mutation.dictionary_ranges['Simple_Phosphorylation.rate']=1.0/mutation.T

# definition of class phosphorylation: sub-class of interaction class (see classes_eds2).
# contains the parameter quantifying the interaction as well as the input and output of the reaction. 
# these characteristics define the interaction completely. they are "free standing" initially in that 
# they are not appended to the graph per se. see new_Phosphorylation below to see how to actually 
# add the interaction to the network.

class Simple_Phosphorylation(classes_eds2.Interaction):
    """ Phosphorylation, init : rate of phosphorylation. """

    def __init__(self,r=0,d=0):
        classes_eds2.Node.__init__(self)
        self.rate=r
        self.label='Simple_Phosphorylation'
        self.input=['Kinase','Phosphorylable']        # these are required to check if the required species within the network are present to start with.
        self.output=['Kinase','Phospho']

    def __str__(self):
        return "{0.id} Phosphorylation: rate = {0.rate:.2f}".format(self)

    def outputs_to_delete(self,net):
        """ Returns the phosphorylated species to delete when deleting a Phosphorylation"""
        return net.find_Simple_phospho_to_remove(self)



####### Attributes attached to Network for  Phosphorylations/Dephosphorylations############

# determines if two species (kinase and protein) are related by phosphorylation interaction (first 
# species should be the kinase).
def check_existing_Simple_Phosphorylation(self,list):
        """ special function to check if a simple phosphorylation exists : order in the list should be Kinase first"""
        if 'Simple_Phosphorylation' in self.list_types:    #goes through the list of interactions
            for inter in self.list_types['Simple_Phosphorylation']:
                [catalyst,listIn,listOut]=self.catal_data(inter)
                if (catalyst==list[0]) and (listIn[0]==list[1]):
                    return True
        return False

def number_Simple_Phosphorylation(self):
        """ Computes the number of possible Simple_Phosphorylations"""
        counter = 0
        nK=len(self.list_types['Kinase'])
        nS=len(self.list_types['Phosphorylable'])
        for iK in range(nK):
            for iS in range(nS):  
                K=self.list_types['Kinase'][iK]
                S=self.list_types['Phosphorylable'][iS]
                Out_S=self.graph.successors(S)
                bool_kinase = True
                bool_K_is_Sp = False
                for i in range(len(Out_S)):
                    if Out_S[i].isinstance('Simple_Phosphorylation'):
                        [catalyst,In_phospho,Out_phospho]=self.catal_data(Out_S[i])
                        Sp = Out_phospho[0]
                        if K == catalyst:
                            bool_kinase = False
                            break
                        if Sp == K:
                            bool_K_is_Sp = True
                            break           
                if not bool_K_is_Sp:  # We need to make sure that the phosphorylation between S and S_p is not catalyzed by Sp!
                    if not K == S: # we do not want phosphorylation of a species by itself.
                        if bool_kinase: # verifying that the species is not already phosphorylated by the same kinase.
                            if not (K.isinstance('pMHC') and S.isinstance('pMHC')): # we do not want a kinase from the cascade to phosphorylate another member in the cascade.
                                counter += 1
        return counter
        
        
        
def find_Simple_phospho_to_remove(self,interaction):
        """Finds products to remove when deleting a Simple_Phosphorylation interaction."""
        
        [kinase,species,species_P]=self.catal_data(interaction)
        already_phosphorylated = False
        
        if species:  # the species could already be deleted.
            species = species[0]
            Out_S=self.graph.successors(species)
            for i in range(len(Out_S)):
                if Out_S[i].isinstance('Simple_Phosphorylation'):
                    [kinase2,species2,species2_P]=self.catal_data(Out_S[i])
                    if species2:
                        species2 = species2[0]
                        if (Out_S[i] != interaction) and (species2 == species):
                            already_phosphorylated = True
                            break
        elif species_P:
           species_P = species_P[0]
           To_Sp = self.graph.predecessors(species_P)
           for j in range(len(To_Sp)):
               if To_Sp[j].isinstance('Simple_Phosphorylation'):
                   [kinase2,species2,species2_P]=self.catal_data(To_Sp[j])
                   if species2_P:
                       species2_P = species2_P[0]
                       if (To_Sp[j] != interaction) and (species_P == species2_P):
                           already_phosphorylated = True
                           break
                    
    
        if not already_phosphorylated:
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
def new_Simple_Phosphorylation(self,kinase,species,rate):
        # we allow for any number of phosphorylations and keep track of the order of phosphorylation through attribute n_phospho of tag 'Phospho'.

        # we do not want a species tAttio phosphorylate itself.
        if kinase == species:
            print("Error in new Simple_Phosphorylation: trying to phosphorylate a species with itself as kinase.")
        
        Out_S=self.graph.successors(species)
        already_phosphorylated = False
        for i in range(len(Out_S)):
            if Out_S[i].isinstance('Simple_Phosphorylation'):
                [catalyst,In_phospho,Out_phospho]=self.catal_data(Out_S[i])
                if In_phospho[0] == species:
                    already_phosphorylated = True
                    species_P = Out_phospho[0]
                    Old_phospho = Out_S[i]
                    break
        
        if not already_phosphorylated:
            phospho=Simple_Phosphorylation(rate) #creates the interaction (see function above)
            species_P=copy.deepcopy(species) # the phosphorylated species has the same properties (need to use deep copy and not just copy: we do not want the pointer but to create an altogether completely new species)
            species_P.clean_type('Input')#Remove Input, output types from the copied species
            species_P.clean_type('Output')#Remove Input, output types from the copied species

            # phosphorylation (keeping track of degree of phosphorylation).
            species_P.n_phospho = species.n_phospho + 1    
            if phospho.check_grammar([kinase,species],[kinase,species_P]):         # verifying that the "grammar" is correct.
                self.add_Node(phospho)                # adding the various nodes and edges corresponding to the interaction.
                self.graph.add_edge(kinase,phospho)
                self.graph.add_edge(species,phospho)
                self.add_Node(species_P)
                self.graph.add_edge(phospho,species_P)
                self.graph.add_edge(phospho,kinase)

                # keeping track of the "pMHC cascade".
                if species.isinstance('pMHC'):
                    #if species.n_phospho == 0: # this is if we take C_0 not a kinase.
                    #    species_P.add_type(['Kinase'])
                    # connecting the newly phosphorylated pMHC components to the ligand/receptor
                    L = self.list_types['Ligand'] #Assuming the correct ligand and receptor are the first in the list (there should be only one anyway...)
                    R = self.list_types['Receptor']
                    self.new_KPR_Unbinding(L[0],R[0],species_P)
                
                
                # we want to transfer all catalytic activity of the species to the phosphorylated species at first.
                # To do so, we track all interactions for which the species is a catalyst and add the corresponding interactions for S_p.
                for i in range(len(Out_S)):
                    if species.isinstance('Kinase'):
                        if Out_S[i].isinstance('Simple_Phosphorylation'):
                            [catalyst,In_phospho,Out_phospho]=self.catal_data(Out_S[i])
                            if catalyst == species:
                                self.new_Simple_Phosphorylation(species_P,In_phospho[0],Out_S[i].rate)
                    if species.isinstance('Phosphatase'):
                        if Out_S[i].isinstance('Simple_Dephosphorylation'):
                            [catalyst,In_phospho,Out_phospho]=self.catal_data(Out_S[i])
                            if catalyst == species:
                                self.new_Simple_Dephosphorylation(species_P,In_phospho[0],Out_phospho[0],Out_S[i].rate)
                                
                return [species_P,phospho]
    
            else:
                print("Error in grammar : new Simple_Phosphorylation")
                return None
             
           
        elif already_phosphorylated:
            phospho=Simple_Phosphorylation(rate) #creates the interaction (see function above)
            if phospho.check_grammar([kinase,species],[kinase,species_P]):         # verifying that the "grammar" is correct.
                self.add_Node(phospho)                # adding the various nodes and edges corresponding to the interaction.
                self.graph.add_edge(kinase,phospho)
                self.graph.add_edge(species,phospho)
                self.graph.add_edge(phospho,species_P)
                self.graph.add_edge(phospho,kinase)
                return [species_P,phospho]
            else:
                print("Error in grammar : new Simple_Phosphorylation")
                return None


setattr(classes_eds2.Network,'check_existing_Simple_Phosphorylation',check_existing_Simple_Phosphorylation)
setattr(classes_eds2.Network,'number_Simple_Phosphorylation',number_Simple_Phosphorylation)
setattr(classes_eds2.Network,'find_Simple_phospho_to_remove',find_Simple_phospho_to_remove)
setattr(classes_eds2.Network,'new_Simple_Phosphorylation',new_Simple_Phosphorylation)



####### Attributes attached to Mutable_Network for  Phosphorylations/Dephosphorylations############


def new_random_Simple_Phosphorylation(self, K, S):
        """ Create new random Simple_Phosphorylation interaction and return interaction and product. Note that the inputs must be given (as opposed to random_Simple_Phosphorylation)."""
        r = mutation.sample_dictionary_ranges('Simple_Phosphorylation.rate',self.Random)
        [SP,Phospho]=self.new_Simple_Phosphorylation(K,S,r)
        return [SP,Phospho]



def random_Simple_Phosphorylation(self):
        """Create new random Simple_Phosphorylations from list of possible kinase substrates."""
        
        if 'Kinase' in self.list_types and 'Phosphorylable' in self.list_types and 'Phosphatase' in self.list_types:
                    
            list_possible_Phospho=[]
            nK=len(self.list_types['Kinase'])
            nS=len(self.list_types['Phosphorylable'])
            for iK in range(nK):
                for iS in range(nS):
                    
                    K=self.list_types['Kinase'][iK]
                    S=self.list_types['Phosphorylable'][iS]
                    
                    Out_S=self.graph.successors(S)
                    bool_kinase = True
                    bool_K_is_Sp = False
                    already_phosphorylated = False
                    for i in range(len(Out_S)):
                        if Out_S[i].isinstance('Simple_Phosphorylation'):
                            [catalyst,In_phospho,Out_phospho]=self.catal_data(Out_S[i])
                            Sp = Out_phospho[0]
                            if S == In_phospho[0]:
                                already_phosphorylated = True
                            if K == catalyst:
                                bool_kinase = False
                                break
                            if Sp == K:
                                bool_K_is_Sp = True
                                break                                           
                
                    if not bool_K_is_Sp:  # We need to make sure that the phosphorylation between S and S_p is not catalyzed by Sp!
                        if not K==S: # we do not want phosphorylation of a species by itself.
                            if bool_kinase: # verifying that the species is not already phosphorylated by the same kinase.
                                if not (K.isinstance('pMHC') and S.isinstance('pMHC')): # we do not want a kinase from the cascade to phosphorylate another member in the cascade.
                                    if not already_phosphorylated:
                                        nP=len(self.list_types['Phosphatase'])
                                        list_possible_phosphatase = []
                                        for m in range(nP):
                                            P = self.list_types['Phosphatase'][m]
                                            if not P==S:
                                                list_possible_phosphatase.append(P)
                                        n_Phosphatase = len(list_possible_phosphatase)
                                        chosen_P = list_possible_phosphatase[int(self.Random.random()*n_Phosphatase)]
                                        list_possible_Phospho.append([K,chosen_P,S])                                        
                                    else:
                                        list_possible_Phospho.append([K,None,S])
                                        #print "We just appended kinase to the list:"
                                        #for property, value in vars(K).iteritems():
                                        #    print property, ": ", value
                                        #print "We just appended species to the list:"
                                        #for property, value in vars(S).iteritems():
                                        #    print property, ": ", value
                                        #print "======================================="
            n_pP=len(list_possible_Phospho)
            if (n_pP==0):
                print("In random_Simple_Phosphorylation : No other posible Simple_Phosphorylationss")
                return None
            else:
                [K,P,S]=list_possible_Phospho[int(self.Random.random()*n_pP)]
                
                if not (P == None):
                    [SP,Phospho]=self.new_random_Simple_Phosphorylation(K,S)
                    dephospho = self.new_random_Simple_Dephosphorylation(P,SP,S)
                    return [SP,Phospho] #[ ,[dephospho]] 
                elif P==None:
                    [SP,Phospho]=self.new_random_Simple_Phosphorylation(K,S)
                    return [SP,Phospho]        
        else:
            print("Error in random_Phosphorylation (try to create Simple_Phosphorylation from non exsiting pieces)")
            return None   


def isphosphorylated_by_K(self,kinase,S):
    Out_S=self.graph.successors(S)
    for i in range(len(Out_S)):
        if Out_S[i].isinstance('Simple_Phosphorylation'):
            [catalyst,In_phospho,Out_phospho]=self.catal_data(Out_S[i])
            if kinase == catalyst:
                return True
    return False



def random_shift_Simple_Phosphorylation(self):
    """Shifts the kinase responsible for a given phosphorylation."""
    
    number_simple_phosphorylation = 0
    if 'Simple_Phosphorylation' in self.list_types:
        number_simple_phosphorylation = len(self.list_types['Simple_Phosphorylation'])
    
    if number_simple_phosphorylation == 0:
        print("No possible Simple_Phosphorylation to shift")
        return None
    
    else:
        simple_phospho = self.list_types['Simple_Phosphorylation'][random.randint(0,number_simple_phosphorylation-1)]
        [kinase,unphosphorylated,phosphorylated]=self.catal_data(simple_phospho)
        unphosphorylated = unphosphorylated[0]
        phosphorylated = phosphorylated[0]
        list_possible_phospho = []
        
        if unphosphorylated.isinstance('pMHC'):
            # we want to find all kinases which are not pMHC molecules.
            for species in self.list_types['Species']:
                if species.isinstance('Kinase') and (not species.isinstance('pMHC')) and (species != phosphorylated) and (species != unphosphorylated):
                    if not self.isphosphorylated_by_K(species,unphosphorylated):
                        list_possible_phospho.append(species)

            if len(list_possible_phospho) == 0:
                print("No possible Simple_Phosphorylation to shift.")
                return None
        
            new_kinase = list_possible_phospho[random.randint(0,len(list_possible_phospho)-1)]
            self.graph.add_edge(new_kinase,simple_phospho)
            self.graph.add_edge(simple_phospho,new_kinase)
            self.graph.delete_edge(simple_phospho,kinase)
            self.graph.delete_edge(kinase,simple_phospho)

            return [simple_phospho,new_kinase]
        else:
            # any kinase will do here.
            # we want to find all kinases (not the species itself or the phosphorylated one) which are not pMHC molecules.
            for species in self.list_types['Species']:
                if species.isinstance('Kinase') and (species != phosphorylated) and (species != unphosphorylated):
                    if not self.isphosphorylated_by_K(species,unphosphorylated):
                        list_possible_phospho.append(species)

            if len(list_possible_phospho) == 0:
                print("No possible Simple_Phosphorylation to shift.")
                return None
            
            new_kinase = list_possible_phospho[random.randint(0,len(list_possible_phospho)-1)]
            self.graph.add_edge(new_kinase,simple_phospho)
            self.graph.add_edge(simple_phospho,new_kinase)
            self.graph.delete_edge(simple_phospho,kinase)
            self.graph.delete_edge(kinase,simple_phospho)
    
            return [simple_phospho,new_kinase]
    
    
    
setattr(mutation.Mutable_Network,'random_Simple_Phosphorylation',random_Simple_Phosphorylation)
setattr(mutation.Mutable_Network,'new_random_Simple_Phosphorylation',new_random_Simple_Phosphorylation)
setattr(mutation.Mutable_Network,'random_shift_Simple_Phosphorylation',random_shift_Simple_Phosphorylation)
setattr(mutation.Mutable_Network,'isphosphorylated_by_K',isphosphorylated_by_K)
#############Integration C Tools #################



def SimplePhospho_deriv_inC(net):
    func="\n/**** Simple_Phosphorylation (without spontaneous dephosphorylation) ****/\n" 
 
    list_phospho = []
    if ('Simple_Phosphorylation' in net.list_types):
        for reaction in net.list_types['Simple_Phosphorylation']:
            [kinase,species,species_P]=net.catal_data(reaction)
            species = species[0]
            species_P = species_P[0]
         
            # agonist term
            rate = "%f*"%reaction.rate+kinase.id+"*"+species.id
            func = func+deriv2.compute_leap([species.id],[species_P.id],rate)     

            # self term
            if species.isinstance('pMHC'):
                n = species.n_phospho
                number_species = len(net.list_types['Species'])
                rate = "%f*"%reaction.rate+kinase.id+"*s[%d]"%(number_species+n+1)
                func = func+deriv2.compute_leap(["s[%d]"%(number_species+n+1)],["s[%d]"%(number_species+n+2)],rate)     
            elif kinase.isinstance('pMHC'):
                n = kinase.n_phospho
                number_species = len(net.list_types['Species'])
                rate = "%f*"%reaction.rate+"s[%d]"%(number_species+n+1)+"*"+species.id
                func = func+deriv2.compute_leap([species.id],[species_P.id],rate)
    return func

deriv2.SimplePhospho_deriv_inC=SimplePhospho_deriv_inC


def compute_gillespie_Simple_Phosphorylation(net,n_reactions):
    """ Create a function computing the Simple_Phosphorylations."""
    proba="\n\t/*****************Simple_Phosphorylations*****************/\n"
    action="\n\t/*****************Simple_Phosphorylations*****************/\n"
    
    if ('Simple_Phosphorylation' in net.list_types):
        for index in net.list_types['Simple_Phosphorylation']:  # looping over all Simple_Phosphorylation reactions
            [K,S,S_p] = net.catal_data(index) # the kinase, initial species and phosphorylated species.
            S = S[0]
            S_p = S_p[0]
            
            # Agonist.
            # the phosphorylation probability (enzymatic)
            proba=proba + "\t\t p[" + str(n_reactions) + "]=%f*floor(s["%index.rate + "%d][ncell])*floor(s["%S.int_id() +"%d][ncell]);\n"%K.int_id()
            # the effect of a phosphorylation (S --> S-1  &  S_p --> S_p + 1 )
            action=action + "\tif (index_action==" + str(n_reactions) + ") {\n\t\ts[%d][ncell] -= INCREMENT;\n"%S.int_id()
            action=action + "\t\ts[%d][ncell] += INCREMENT;\n\t}\n"%S_p.int_id()
                
            n_reactions+=1
            
            # Self.
            if S.isinstance('pMHC'):
                n = S.n_phospho
                number_species = len(net.list_types['Species'])
                
                # the phosphorylation probability (enzymatic)
                proba=proba + "\t\t p[" + str(n_reactions) + "]=%f*floor(s["%index.rate + "%d][ncell])*floor(s["%(number_species+n+1) +"%d][ncell]);\n"%K.int_id()            
                # the effect of a phosphorylation (S --> S-1  &  S_p --> S_p + 1 )
                action=action + "\tif (index_action==" + str(n_reactions) + ") {\n\t\ts[%d][ncell] -= INCREMENT;\n"%(number_species+n+1)
                action=action + "\t\ts[%d][ncell] += INCREMENT;\n\t}\n"%(number_species+n+2)
                
                n_reactions+=1
                
                
            elif K.isinstance('pMHC'):
                n = K.n_phospho
                number_species = len(net.list_types['Species'])
                
                # the phosphorylation probability (enzymatic)
                proba=proba + "\t\t p[" + str(n_reactions) + "]=%f*floor(s["%index.rate + "%d][ncell])*floor(s["%S.int_id() +"%d][ncell]);\n"%(number_species+n+1)           
                # the effect of a phosphorylation (S --> S-1  &  S_p --> S_p + 1 )
                action=action + "\tif (index_action==" + str(n_reactions) + ") {\n\t\ts[%d][ncell] -= INCREMENT;\n"%S.int_id()
                action=action + "\t\ts[%d][ncell] += INCREMENT;\n\t}\n"%S_p.int_id()
                
                n_reactions+=1
            
            
    return [proba,action,n_reactions]
    
#gillespie_pMHC.compute_gillespie_Simple_Phosphorylation=compute_gillespie_Simple_Phosphorylation
