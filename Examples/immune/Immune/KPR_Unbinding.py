import phievo.Networks.classes_eds2 as classes_eds2
import phievo.Networks.mutation as mutation
import copy
import phievo.Networks.deriv2 as deriv2

mutation.dictionary_ranges['KPR_Unbinding.association']=1.0/mutation.T
mutation.dictionary_ranges['KPR_Unbinding.dissociation']=mutation.C


class KPR_Unbinding(classes_eds2.Interaction):
    """ Ligand-Receptor interaction, init :disassociation"""
    def __init__(self): #,d=0
        classes_eds2.Node.__init__(self)
        #self.dissociation=d
        self.label='KPR_Unbinding'
        self.input=['Species']
        self.output=['Ligand','Receptor']

    def __str__(self):
        return "{0.id} Unbinding (KPR)".format(self)

     #def outputs_to_delete(self,net):
     #     """ Returns Ligand-Receptor complex as physical object to delete for a KPR_Unbinding"""
     #     return net.graph.successors(self)

####### Attributes attached to Network for  KPR_Unbinding interaction ######################################


def number_KPR_Unbinding(self):
        """ Computes the number of possible KPR_Unbinding interactions"""
        npMHC=self.number_nodes('pMHC')
        n_KPR_Unbinding=self.number_nodes('KPR_Unbinding')
        return npMHC-n_KPR_Unbinding


def new_KPR_Unbinding(self,L,R,C):
        """Create a new KPR_Unbinding-interaction and the complex it creates.
        args: species1, species2, complex (disassoc)"""
        
        p=KPR_Unbinding()
        if p.check_grammar([C], [L,R]):
            self.graph.add_edge(C,p)
            self.graph.add_edge(p,L)
            self.graph.add_edge(p,R)
            return p
        else:
            print("Error in grammar, new_KPR_Unbinding")
            return None


setattr(classes_eds2.Network,'number_KPR_Unbinding',number_KPR_Unbinding)
setattr(classes_eds2.Network,'new_KPR_Unbinding',new_KPR_Unbinding)


####### Attributes attached to Mutable_Network for  KPR_Unbinding interaction ######################################



def random_KPR_Unbinding(self):
        """Create new random KPR_Unbinding from list of those possible """
        if 'Ligand' in self.list_types and 'Receptor' in self.list_types:
            list_possible_KPR_Unbinding=[]
            nL=len(self.list_types['Ligand'])
            nR=len(self.list_types['Receptor'])
            for iL in range(nL):
                for iR in range(nR):
                    L=self.list_types['Ligand'][iL]
                    R=self.list_types['Receptor'][iR]
                    if not self.check_existing_binary([L,R],'KPR_Unbinding'):
                        list_possible_KPR_Unbinding.append([L,R])
            n_pKPR_Unbinding=len(list_possible_KPR_Unbinding)
            if not (n_pKPR_Unbinding==self.number_KPR_Unbinding()):
                print("Potential Bug : Inconsistency in Computation of number of KPR_Unbinding")
            if (n_pKPR_Unbinding==0):
                print("In random_KPR_Unbinding : No other posible KPR_Unbindings")
                return None
            else:
                [L,R]=list_possible_KPR_Unbinding[int(self.Random.random()*n_pKPR_Unbinding)]
                [KPR_Unbinding,C]=self.new_random_KPR_Unbinding(L, R)
                return [KPR_Unbinding,C]     
        else:
            print("Error in random_KPR_Unbinding (try to create a KPR_Unbinding from non exsiting pieces)")
            return None        

def new_random_KPR_Unbinding(self, L, R):
        """create random KPR_Unbinding, and return interaction Node and output complex """
        
        list=[['Degradable', mutation.sample_dictionary_ranges('Species.degradation',self.Random) ],['TF',int(self.Random.random()*2)],['Complexable'],['Phosphorylable']]
        a = mutation.sample_dictionary_ranges('KPR_Unbinding.association',self.Random)
        conc = mutation.sample_dictionary_ranges('KPR_Unbinding.dissociation',self.Random)
        [KPR_Unbinding,C]=self.new_KPR_Unbinding(L,R,a,conc,list)
        return [KPR_Unbinding,C]

    


setattr(mutation.Mutable_Network,'random_KPR_Unbinding',random_KPR_Unbinding)
setattr(mutation.Mutable_Network,'new_random_KPR_Unbinding',new_random_KPR_Unbinding)




############Integration C Tools ##############


def compute_KPR_Unbinding(net):
    """ Create a function that computes the derivative for Receptor and product of KPR_Unbinding interaction for a given cell"""
    func="\n/***********KPR_Unbinding**************/\n"

    if ('KPR_Unbinding' in net.list_types):
        for index in net.list_types['KPR_Unbinding']:
            C=net.graph.predecessors(index)[0] #finds the complex
            [P1,P2]=net.graph.successors(index) #finds the receptor and ligands
            if P1.isinstance('Ligand'):
                L=P1
                R=P2
            else:
                R=P1
                L=P2
            # agonist term
            drate="1.0/tau*"+C.id
            func=func+deriv2.compute_leap([C.id],[L.id,R.id],drate)
            # self term
            n = C.n_phospho
            number_species = len(net.list_types['Species'])
            drate="1.0/TAU_SELF*s[%d]"%(number_species+n+1)
            func=func+deriv2.compute_leap(["s[%d]"%(number_species+n+1)],["s[%d]"%number_species,R.id],drate)
    return func


deriv2.compute_KPR_Unbinding=compute_KPR_Unbinding


                
def compute_gillespie_KPR_Unbinding(net,n_reactions):
    """ Create a function computing the KPR_Unbinding."""
    proba="\n\t/*****************KPR_Unbinding*****************/\n"
    action="\n\t/*****************KPR_Unbinding*****************/\n"
    
    if ('KPR_Unbinding' in net.list_types):
        for index in net.list_types['KPR_Unbinding']:  # looping over all Simple_Dephosphorylation reactions
            C = net.graph.predecessors(index)[0]  # a complex
            [P1,P2] = net.graph.successors(index)  # the receptor and ligand
            if (P1.isinstance('Ligand')):
                L = P1
                R = P2
            else:
                L = P2
                R = P1
            # Agonist.    
            # the unbinding probability
            proba=proba + "\t \t p[" + str(n_reactions) + "]=1.0/tau*floor(s[%d][ncell]);\n"%C.int_id()
            # the effect of unbinding (R,L --> R+1,L+1  &  C --> C - 1 )
            action=action + "\tif (index_action==" + str(n_reactions) + ") {\n\t\ts[%d][ncell] += INCREMENT;\n"%R.int_id()
            action=action + "\t\ts[%d][ncell] += INCREMENT;\n"%L.int_id()
            action=action + "\t\ts[%d][ncell] -= INCREMENT;\n\t}\n"%C.int_id()
                
            n_reactions+=1
            
            # Self.
            # the unbinding probability
            n = C.n_phospho
            number_species = len(net.list_types['Species'])
            proba=proba + "\t \t p[" + str(n_reactions) + "]=1.0/TAU_SELF*floor(s[%d][ncell]);\n"%(number_species+n+1)
            # the effect of unbinding (R,L --> R+1,L+1  &  C --> C - 1 )
            action=action + "\tif (index_action==" + str(n_reactions) + ") {\n\t\ts[%d][ncell] += INCREMENT;\n"%R.int_id()
            action=action + "\t\ts[%d][ncell] += INCREMENT;\n"%number_species
            action=action + "\t\ts[%d][ncell] -= INCREMENT;\n\t}\n"%(number_species+n+1)
                
            n_reactions+=1
    return [proba,action,n_reactions]

#gillespie_pMHC.compute_gillespie_KPR_Unbinding = compute_gillespie_KPR_Unbinding
