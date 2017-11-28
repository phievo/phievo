import phievo.Networks.classes_eds2 as classes_eds2
import phievo.Networks.mutation as mutation
import copy
import phievo.Networks.deriv2 as deriv2

mutation.dictionary_ranges['KPR_Binding.association']=1.0/mutation.T
mutation.dictionary_ranges['KPR_Binding.dissociation']=mutation.C



# Ligand tag
mutation.species_types["Ligand"] = lambda random_generator:[
    ["Ligand"],
]
classes_eds2.Species.Tags_Species["Ligand"] = []

# Receptor tag
mutation.species_types["Receptor"] = lambda random_generator:[
    ["Receptor"],
]
classes_eds2.Species.Tags_Species["Receptor"] = []

class KPR_Binding(classes_eds2.Interaction):
    """ Ligand-Receptor interaction, init : rates of association"""
    def __init__(self,a=0):
        classes_eds2.Node.__init__(self)
        self.association=a
        self.label='KPR_Binding'
        self.input=['Ligand','Receptor']
        self.output=['Species']

    def __str__(self):
        return "{0.id} Binding (KPR): assoc. = {0.association:.2f}".format(self)

    def outputs_to_delete(self,net):
        """ Returns Ligand-Receptor complex as physical object to delete for a KPR_Binding"""
        return net.graph.successors(self) 


####### Attributes attached to Network for  KPR_Binding interaction ######################################


def number_KPR_Binding(self):
        """ Computes the number of possible KPR_Binding interactions"""
        nL=self.number_nodes('Ligand')
        nR=self.number_nodes('Receptor')
        n_KPR_Binding=self.number_nodes('KPR_Binding')
        return nL*nR-n_KPR_Binding
     

def new_KPR_Binding(self,L,R,a,list):
        """Create a new KPR_Binding-interaction and the complex it creates.
        args: species1, species2 (objects), assoc, list of characteristic of the complex"""
        
        c=classes_eds2.Species(list)
        p=KPR_Binding(a)
        if p.check_grammar([L,R], [c]):
            self.add_Node(c)
            self.add_Node(p)
            self.graph.add_edge(L,p)
            self.graph.add_edge(R,p)
            self.graph.add_edge(p,c)
            return [p,c]
        else:
            print("Error in grammar, new_KPR_Binding")
            return None




setattr(classes_eds2.Network,'number_KPR_Binding',number_KPR_Binding)
setattr(classes_eds2.Network,'new_KPR_Binding',new_KPR_Binding)


####### Attributes attached to Mutable_Network for  KPR_Binding interaction ######################################



def random_KPR_Binding(self):
        """Create new random KPR_Binding from list of those possible """
        if 'Ligand' in self.dict_types and 'Receptor' in self.dict_types:
            list_possible_KPR_Binding=[]
            nL=len(self.dict_types['Ligand'])
            nR=len(self.dict_types['Receptor'])
            for iL in range(nL):
                for iR in range(nR):
                    L=self.dict_types['Ligand'][iL]
                    R=self.dict_types['Receptor'][iR]
                    if not self.check_existing_binary([L,R],'KPR_Binding'):
                        list_possible_KPR_Binding.append([L,R])
            n_pKPR_Binding=len(list_possible_KPR_Binding)
            if not (n_pKPR_Binding==self.number_KPR_Binding()):
                print("Potential Bug : Inconsistency in Computation of number of KPR_Binding")
            if (n_pKPR_Binding==0):
                print("In random_KPR_Binding : No other posible KPR_Bindings")
                return None
            else:
                [L,R]=list_possible_KPR_Binding[int(self.Random.random()*n_pKPR_Binding)]
                [KPR_Binding,C]=self.new_random_KPR_Binding(L, R)
                return [KPR_Binding,C]     
        else:
            print("Error in random_KPR_Binding (try to create a KPR_Binding from non exsiting pieces)")
            return None        

def new_random_KPR_Binding(self, L, R):
        """create random KPR_Binding, and return interaction Node and output complex """
        
        list=[['Degradable', mutation.sample_dictionary_ranges('Species.degradation',self.Random) ],['TF',int(self.Random.random()*2)],['Complexable'],['Phosphorylable']]
        a = mutation.sample_dictionary_ranges('KPR_Binding.association',self.Random)
        conc = mutation.sample_dictionary_ranges('KPR_Binding.dissociation',self.Random)
        [KPR_Binding,C]=self.new_KPR_Binding(L,R,a,conc,list)
        return [KPR_Binding,C]

    
setattr(mutation.Mutable_Network,'random_KPR_Binding',random_KPR_Binding)
setattr(mutation.Mutable_Network,'new_random_KPR_Binding',new_random_KPR_Binding)




############Integration C Tools ##############




def compute_KPR_Binding(net):
    """ Create a function that computes the derivative for Receptor and product of KPR_Binding interaction for a given cell"""
    func="\n/***********KPR_Binding**************/\n"
    if ('KPR_Binding' in net.dict_types):
        for index in net.dict_types['KPR_Binding']:
            C=net.graph.successors(index)[0]#finds the product of KPR_Binding interaction
            [P1,P2]=net.graph.predecessors(index) #find the components
            if P1.isinstance('Ligand'):
                L=P1
                R=P2
            else:
                R=P1
                L=P2
            # agonist term
            arate="%.10f*"%index.association+L.id+"*"+R.id
            func=func+deriv2.compute_leap([R.id,L.id],[C.id],arate)
            # self term
            number_species = len(net.dict_types['Species'])
            arate="%.10f*"%index.association+"s[%d]"%number_species+"*"+R.id;
            func = func+deriv2.compute_leap([R.id,"s[%d]"%number_species],["s[%d]"%(number_species+1)],arate);
    return func
deriv2.compute_KPR_Binding=compute_KPR_Binding


                            
def compute_gillespie_KPR_Binding(net,n_reactions):
    """ Create a function computing the KPR_Binding."""
    proba="\n\t/*****************KPR_Binding*****************/\n"
    action="\n\t/*****************KPR_Binding*****************/\n"
    if ('KPR_Binding' in net.dict_types):
        for index in net.dict_types['KPR_Binding']:  # looping over all Simple_Dephosphorylation reactions
            C = net.graph.successors(index)[0]  # the complex
            [P1,P2] = net.graph.predecessors(index)  # the receptor and ligand
            if (P1.isinstance('Ligand')):
                L = P1
                R = P2
            else:
                L = P2
                R = P1
            # Agonist.
            # the binding probability
            proba=proba + "\t \t p[" + str(n_reactions) + "]=%f*floor(s["%index.association + "%d][ncell])*floor(s["%R.int_id()  + "%d][ncell]);\n"%L.int_id()
            # the effect of binding (R,L --> R-1,L-1  &  C --> C + 1 )
            action=action + "\tif (index_action==" + str(n_reactions) + ") {\n\t\ts[%d][ncell] -= INCREMENT;\n"%R.int_id()
            action=action + "\t\ts[%d][ncell] -= INCREMENT;\n"%L.int_id()
            action=action + "\t\ts[%d][ncell] += INCREMENT;\n\t}\n"%C.int_id()
            n_reactions+=1 
            
            # Self.
            number_species = len(net.dict_types['Species'])
            # the binding probability
            proba=proba + "\t \t p[" + str(n_reactions) + "]=%.10f*floor(s["%index.association + "%d][ncell])*floor(s["%R.int_id()  + "%d][ncell]);\n"%number_species
            # the effect of binding (R,L --> R-1,L-1  &  C --> C + 1 )
            action=action + "\tif (index_action==" + str(n_reactions) + ") {\n\t\ts[%d][ncell] -= INCREMENT;\n"%R.int_id()
            action=action + "\t\ts[%d][ncell] -= INCREMENT;\n"%number_species
            action=action + "\t\ts[%d][ncell] += INCREMENT;\n\t}\n"%(number_species+1)
            n_reactions+=1 
    return [proba,action,n_reactions]
#gillespie_pMHC.compute_gillespie_KPR_Binding=compute_gillespie_KPR_Binding


