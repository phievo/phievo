# importing all relevant modules.
import phievo.Networks.classes_eds2 as classes_eds2
import phievo.Networks.mutation as mutation
import copy
from phievo.initialization_code import *
import phievo.Networks.deriv2 as deriv2


# parameters needed to perform evolution.
mutation.dictionary_ranges['Simple_Dephosphorylation.rate']=1.0/mutation.T

# definition of class phosphorylation: sub-class of interaction class (see classes_eds2).
# contains the parameter quantifying the interaction as well as the input and output of the reaction. 
# these characteristics define the interaction completely. they are "free standing" initially in that 
# they are not appended to the graph per se. see new_Phosphorylation below to see how to actually 
# add the interaction to the network.

class Simple_Dephosphorylation(classes_eds2.Interaction):
    """ Phosphorylation, init : rate of phosphorylation. """

    def __init__(self,r=0):
        classes_eds2.Node.__init__(self)
        self.rate=r
        self.label='Simple_Dephosphorylation'
        self.input=['Phosphatase','Phospho']        # these are required to check if the required species within the network are present to start with.
        self.output=['Phosphatase','Phosphorylable']

    def __str__(self):
        return "{0.id} Dephosphorylation: rate = {0.rate:.2f}".format(self)

    def outputs_to_delete(self,net):
        """ Returns the phosphorylated species to delete when deleting a Phosphorylation"""
        return net.find_phosphorylated_to_remove(self)



####### Attributes attached to Network for  Phosphorylations/Dephosphorylations############

# determines if two species (kinase and protein) are related by phosphorylation interaction (first 
# species should be the kinase).
def check_existing_Simple_Dephosphorylation(self,list):
        """ special function to check if a dephosphorylation exists : order in the list should be phosphatase first"""
        if 'Simple_Dephosphorylation' in self.dict_types:    #goes through the list of interactions
            for inter in self.dict_types['Simple_Dephosphorylation']:
                [catalyst,listIn,listOut]=self.catal_data(inter)
                if (catalyst==list[0]) and (listIn[0]==list[1]):
                    return True
        return False

def is_phospho_successor(self,Sp,S):
    listOut = self.graph.successors(S)
    for i in range(len(listOut)):
        if listOut[i].isinstance('Simple_Phosphorylation'):
            [catalyst,In_phospho,Out_phospho]=self.catal_data(listOut[i])
            if catalyst == S: # Avoids the situation where S is the kinase.
                continue
            if Sp == Out_phospho[0]:
                return True
    return False

def number_Simple_Dephosphorylation(self):
        """ Computes the number of possible Simple_Dephosphorylations"""
        counter = 0
        nP=self.number_nodes('Phosphatase')
        nSp=self.number_nodes('Phospho')
        nS=self.number_nodes('Phosphorylable')
        for iP in range(nP):
            for iSp in range(nSp):
                for iS in range(nS):
                    P=self.dict_types['Phosphatase'][iP]
                    Sp=self.dict_types['Phospho'][iSp]
                    S=self.dict_types['Phosphorylable'][iS]
                    if self.is_phospho_successor(Sp,S): # Need the two species to be connected by a phosphorylation.
                        if P != S: # We do not want the phosphatase to be the species itself.
                            if P != Sp: # We do not want the phosphatase to be the phosphorylated species.
                                if not(Sp.isinstance('pMHC') != S.isinstance('pMHC')): # Not really necessary from our construction of Simple_Phosphorylation, but a double safety...
                                    if not self.check_existing_Simple_Dephosphorylation([P,Sp,S]):
                                        counter += 1
        return counter

def find_phosphorylated_to_remove(self,interaction):
        """Finds products to remove when deleting a Simple_Dephosphorylation interaction. The output to remove is S_p when the dephosphorylation between S_p and S is unique. """
        
        [phosphatase,species_P,species]=self.catal_data(interaction)
        already_dephosphorylated = False
        
        if species_P:
            species_P = species_P[0]
            Out_Sp=self.graph.successors(species_P)
            for i in range(len(Out_Sp)):
                if Out_Sp[i].isinstance('Simple_Dephosphorylation'):
                    [phosphatase2,species2_P,species2]=self.catal_data(Out_Sp[i])
                    if species2_P:
                        species2_P = species2_P[0]
                        if (Out_Sp[i] != interaction) and (species2_P == species_P):
                            already_dephosphorylated = True
                            break
        elif species:
            species = species[0]
            To_S=self.graph.predecessors(species);
            for j in range(len(To_S)):
                if To_S[j].isinstance('Simple_Dephosphorylation'):
                    [phosphatase2,species2_P,species2]=self.catal_data(To_S[j])
                    if species2:
                        species2 = species2[0]
                        if (To_S[j] != interaction) and (species2 == species):
                            already_dephosphorylated = True
                            break
                    
                    
        
        if not already_dephosphorylated:
            listIn=self.graph.predecessors(interaction)
            listOut=self.graph.successors(interaction)
            Bool=False
            for x in listIn:
                if x in listOut:
                    catalyst=x
                    Bool=True
                    #there is this Bool business to intercept the case where the catalyst has been already removed by clean_Nodes, we then just remove the output list
            if Bool:
                listIn.remove(catalyst)
            return listIn

        else:
            return []
            
        #listIn=self.graph.predecessors(interaction)
        #listOut=self.graph.successors(interaction)
        #Bool=False
        #for x in listIn:
        #    if x in listOut:
        #        catalyst=x
        #        Bool=True
        #        #there is this Bool business to intercept the case where the catalyst has been already removed by clean_Nodes, we then just remove the output list
        #if Bool:
        #    listOut.remove(catalyst)
        #return listOut


# important function that adds the interaction to the actual graph. 
def new_Simple_Dephosphorylation(self,phosphatase,species1,species2,rate):

        dephospho=Simple_Dephosphorylation(rate) #creates the interaction (see function above, an instance of the Node classe)
        assert abs(species1.n_phospho - species2.n_phospho) ==1, "grammar error"
        
        if(species1.n_phospho - species2.n_phospho == 1):#verifies that the degree of phosphorylation of the two species differ by one.            
            if dephospho.check_grammar([phosphatase,species1],[phosphatase,species2]):                
                # adding the interaction node
                self.add_Node(dephospho)
                # what goes in the interaction node
                self.graph.add_edge(phosphatase,dephospho)
                self.graph.add_edge(species1,dephospho)
                # what comes out of the interaction node
                self.graph.add_edge(dephospho,species2)
                self.graph.add_edge(dephospho,phosphatase)
                return dephospho
                
        elif(species2.n_phospho - species1.n_phospho == 1):
           if dephospho.check_grammar([phosphatase,species2],[phosphatase,species1]):
                # adding the interaction node
                self.add_Node(dephospho)
                # what goes in the interaction node
                self.graph.add_edge(phosphatase,dephospho)
                self.graph.add_edge(species2,dephospho)
                # what comes out of the interaction node
                self.graph.add_edge(dephospho,species1)
                self.graph.add_edge(dephospho,phosphatase)
                return dephospho
           
        else:
            print("Error in grammar : new Simple_Dephosphorylation")
            return None


setattr(classes_eds2.Network,'check_existing_Simple_Dephosphorylation',check_existing_Simple_Dephosphorylation)
setattr(classes_eds2.Network,'number_Simple_Dephosphorylation',number_Simple_Dephosphorylation)
setattr(classes_eds2.Network,'find_phosphorylated_to_remove',find_phosphorylated_to_remove)
setattr(classes_eds2.Network,'new_Simple_Dephosphorylation',new_Simple_Dephosphorylation)



####### Attributes attached to Mutable_Network for  Phosphorylations/Dephosphorylations############


def new_random_Simple_Dephosphorylation(self, K, Sp,S):
        """ create new random Simple_Dephosphorylation interaction and return interaction."""
        
        r = mutation.sample_dictionary_ranges('Simple_Dephosphorylation.rate',self.Random)
        D=self.new_Simple_Dephosphorylation(K,Sp,S,r)
        return D
    
    
def random_Simple_Dephosphorylation(self):
        """Create new random  Phosphorylations from list of possible kinase substrates"""
        if 'Phosphatase' in self.dict_types and 'Phosphorylable' in self.dict_types and 'Phospho' in self.dict_types:
            list_possible_Simple_Dephosphorylation=[]
            nP=self.number_nodes('Phosphatase')
            nSp=self.number_nodes('Phospho')
            nS=self.number_nodes('Phosphorylable')
            
            for iP in range(nP):
                for iSp in range(nSp):
                    for iS in range(nS):
                        P=self.dict_types['Phosphatase'][iP]
                        Sp=self.dict_types['Phospho'][iSp]
                        S=self.dict_types['Phosphorylable'][iS]
                    
                        if self.is_phospho_successor(Sp,S): # Need the two species to be connected by a phosphorylation.
                            if P != S: # We do not want the phosphatase to be the species itself.
                                if P != Sp: # We do not want the phosphatase to be the phosphorylated species.
                                    if not(Sp.isinstance('pMHC') != S.isinstance('pMHC')): # Not really necessary from our construction of Simple_Phosphorylation, but a double safety...
                                        if not self.check_existing_Simple_Dephosphorylation([P,Sp,S]):
                                            list_possible_Simple_Dephosphorylation.append([P,Sp,S])
                            
            n_pP=len(list_possible_Simple_Dephosphorylation)
            #if not (n_pP==self.number_Phosphorylation()):
            #    print "Potential Bug : Inconsistency in Computation of number of Phosphorylations"
            #    print n_pP,self.number_Phosphorylation()
            if (n_pP==0):
                print("In random_Simple_Dephosphorylation : No other posible Simple_Dephosphorylation")
                return None
            else:
                [P,Sp,S]=list_possible_Simple_Dephosphorylation[int(self.Random.random()*n_pP)]
                D=self.new_random_Simple_Dephosphorylation(P,Sp,S)
                return D    
        else:
            print("Error in random_Phosphorylation (try to create Phosphorylation from non exsiting pieces)")
            return None   




setattr(mutation.Mutable_Network,'random_Simple_Dephosphorylation',random_Simple_Dephosphorylation)
setattr(mutation.Mutable_Network,'is_phospho_successor',is_phospho_successor)
setattr(mutation.Mutable_Network,'new_random_Simple_Dephosphorylation',new_random_Simple_Dephosphorylation)


#############Integration C Tools #################


## Deterministic stuff.
def SimpleDephospho_deriv_inC(net):
 func="\n/***********Simple_Dephosphorylation**************/\n"
 if ('Simple_Dephosphorylation' in net.dict_types):
     for reaction in net.dict_types['Simple_Dephosphorylation']:
         [phosphatase,species_P,species]=net.catal_data(reaction)
         species = species[0]
         species_P = species_P[0]
         phosphatase = phosphatase[0]
         # agonist term
         rate = "%f*"%reaction.rate+phosphatase.id+"*"+species_P.id
         func = func+deriv2.compute_leap([species_P.id],[species.id],rate)
         # self term
         if species.isinstance('pMHC'):
             n = species.n_phospho
             number_species = len(net.dict_types['Species'])
             rate = "%f*"%reaction.rate+phosphatase.id+"*s[%d]"%(number_species+n+2)
             func = func+deriv2.compute_leap(["s[%d]"%(number_species+n+2)],["s[%d]"%(number_species+n+1)],rate)
         elif phosphatase.isinstance('pMHC'):
             n = phosphatase.n_phospho
             rate = "%f*"%reaction.rate+"s[%d]"%(number_species+n+1)+"*"+species_P.id
             func = func+deriv2.compute_leap([species_P.id],[species.id],rate)
         
         
 return func

deriv2.SimpleDephospho_deriv_inC=SimpleDephospho_deriv_inC


## Gillespie stuff.
def compute_gillespie_Simple_Dephosphorylation(net,n_reactions):
    """ Create a function computing the Simple_Dephosphorylations."""
    proba="\n\t/*****************Simple_Dephosphorylations*****************/\n"
    action="\n\t/*****************Simple_Dephosphorylations*****************/\n"
    
    if ('Simple_Dephosphorylation' in net.dict_types):
        for index in net.dict_types['Simple_Dephosphorylation']:  # looping over all Simple_Dephosphorylation reactions
            [P,S_p,S] = net.catal_data(index) # the phosphatase, initial phosphorylated species and species.
            S = S[0]
            S_p = S_p[0]
            # Agonist.
            # the dephosphorylation probability (enzymatic)
            proba=proba + "\t \t p[" + str(n_reactions) + "]=%f*floor(s["%index.rate + "%d][ncell])*floor(s["%S_p.int_id() + "%d][ncell]);\n"%P.int_id()
            # the effect of a dephosphorylation (S_p --> S_p-1  &  S --> S + 1 )
            action=action + "\tif (index_action==" + str(n_reactions) + ") {\n\t\ts[%d][ncell] -= INCREMENT;\n"%S_p.int_id()
            action=action + "\t\ts[%d][ncell] += INCREMENT;\n\t}\n"%S.int_id()
                
            n_reactions+=1
            
            # Self.
            if S.isinstance('pMHC'):
                n = S.n_phospho
                number_species = len(net.dict_types['Species'])
                # the dephosphorylation probability (enzymatic)
                proba=proba + "\t \t p[" + str(n_reactions) + "]=%f*floor(s["%index.rate + "%d][ncell])*floor(s["%(number_species+n+2) + "%d][ncell]);\n"%P.int_id()
                # the effect of a dephosphorylation (S_p --> S_p-1  &  S --> S + 1 )
                action=action + "\tif (index_action==" + str(n_reactions) + ") {\n\t\ts[%d][ncell] -= INCREMENT;\n"%(number_species+n+2)
                action=action + "\t\ts[%d][ncell] += INCREMENT;\n\t}\n"%(number_species+n+1)
                
                n_reactions+=1
                
            elif P.isinstance('pMHC'):
                n = P.n_phospho
                number_species = len(net.dict_types['Species'])
                # the dephosphorylation probability (enzymatic)
                proba=proba + "\t \t p[" + str(n_reactions) + "]=%f*floor(s["%index.rate + "%d][ncell])*floor(s["%S_p.int_id() + "%d][ncell]);\n"%(number_species+n+1)
                # the effect of a dephosphorylation (S_p --> S_p-1  &  S --> S + 1 )
                action=action + "\tif (index_action==" + str(n_reactions) + ") {\n\t\ts[%d][ncell] -= INCREMENT;\n"%S_p.int_id()
                action=action + "\t\ts[%d][ncell] += INCREMENT;\n\t}\n"%S.int_id()
                
                n_reactions+=1
                     
    return [proba,action,n_reactions]
    
#gillespie_pMHC.compute_gillespie_Simple_Dephosphorylation = compute_gillespie_Simple_Dephosphorylation   
 

