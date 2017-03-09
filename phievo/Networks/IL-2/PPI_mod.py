import classes_eds2
import mutation
import copy
import config
exec('import '+config.name_deriv2+' as deriv2')

mutation.dictionary_ranges['PPI_mod.association']=1.0/(mutation.C*mutation.T)
mutation.dictionary_ranges['PPI_mod.disassociation']=1.0/mutation.T


class PPI_mod(classes_eds2.Interaction):
    """ Protein-protein interaction, init : rates of association, disassociation"""
    
    def __init__(self,a=0,d=0):
        classes_eds2.Node.__init__(self)
        self.association=a
        self.disassociation=d
        self.label='PP Interaction'
        self.input=['Complexable','Complexable']
        self.output=['Species']
    

    def outputs_to_delete(self,net):
        """ Returns complex as physical object to delete for a PPI_mod"""
        return net.graph.successors(self)




####### Attributes attached to Network for PPI_mod############


def number_PPI_mod(self):
    """return number of possible PPI_mods in network"""
    n=self.number_nodes('Complexable')
    n_PPI_mod=self.number_nodes('PPI_mod')
    return n*(n+1)/2-n_PPI_mod #n(n-1)/2+n - self interactions


def new_PPI_mod(self,P1,P2,a,d,list):
        """Create a new PP-interaction and the complex it creates + associated edges
           args: species1, species2 (objects), assoc, dis-assoc (rates, PPI_mod),
           list of characteristic of the complex"""
        c=classes_eds2.Species(list)
        p=PPI_mod(a,d)
        if p.check_grammar([P1,P2], [c]):
            self.add_Node(c)
            self.add_Node(p)
            self.graph.add_edge(P1,p)
            self.graph.add_edge(P2,p)
            self.graph.add_edge(p,c)
            return [p,c]
        else:
            print("Error in grammar, new_Complex")
        return None

def duplicate_PPI_mod(self,species,D_species,interaction,module,D_module):
        """function to duplicate a PPI_mod interaction, D_Species is the duplicated component of the complex, species is the original protein,  D_module is the Tmodule that is duplicated in this gene duplication, the original being module"""
        D_interaction=copy.deepcopy(interaction)
        D_interaction.mutable=1
        D_interaction.removable=True
        self.add_Node(D_interaction)
        Complex=self.graph.successors(interaction)[0]
        D_Complex=copy.deepcopy(Complex)#copies the complex
        self.add_Node(D_Complex)
        D_Complex.mutable=1
        D_Complex.removable=True
        #D_Complex.print_node()
        self.graph.add_edge(D_interaction,D_Complex)#copies the PPI_mod
        PPI_mod_components=self.graph.predecessors(interaction)
        PPI_mod_components.remove(species) # removes the duplicated component
        component=PPI_mod_components[0] #other PPI_mod component
        self.graph.add_edge(component,D_interaction)
        self.graph.add_edge(D_species,D_interaction)
        self.duplicate_downstream_interactions(Complex,D_Complex,module,D_module)
        if (component==species): #if species complexes with itself we must create a new dimer
            D_Complex_2=copy.deepcopy(Complex)#copies the complex
            D_Complex_2.mutable=1
            D_Complex_2.removable=True
            D_interaction_2=copy.deepcopy(interaction)
            D_interaction_2.mutable=1
            D_interaction_2.removable=True
            self.add_Node(D_interaction_2)
            self.graph.add_edge(D_interaction_2,D_Complex_2)#copies the PPI_mod
            self.graph.add_edge(D_species,D_interaction_2)
            self.graph.add_edge(D_species,D_interaction_2)
            self.duplicate_downstream_interactions(Complex,D_Complex_2,module,D_module)

setattr(classes_eds2.Network,'number_PPI_mod',number_PPI_mod)
setattr(classes_eds2.Network,'new_PPI_mod',new_PPI_mod)
setattr(classes_eds2.Network,'duplicate_PPI_mod',duplicate_PPI_mod)


####### Attributes attached to Mutable_Network for PPI_mod############



def random_PPI_mod(self):
        """Create new random PPI_mod from list of those possible """
        if 'Complexable' in self.list_types:
            list_complexable=self.list_types['Complexable']
            n=len(list_complexable)
            list_possible_PPI_mod=[]
            for ip1 in range(n):
                for ip2 in range(ip1+1):
                    P1=list_complexable[ip1]
                    P2=list_complexable[ip2]
                    if not self.check_existing_binary([P1,P2],'PPI_mod'):
                        list_possible_PPI_mod.append([P1,P2])
            n_pPPI_mod=len(list_possible_PPI_mod)
            if not (n_pPPI_mod==self.number_PPI_mod()):
                print("Potential Bug : Inconsistency in Computation of number of PPI_mod")
            if (n_pPPI_mod==0):
                print("In random_PPI_mod : No other posible PPI_mods")
                return None
            else:
                [P1,P2]=list_possible_PPI_mod[int(self.Random.random()*n_pPPI_mod)]
                [PPI_mod,C]=self.new_random_PPI_mod(P1,P2)
                return [PPI_mod,C]   
        else:
            print("Error in random_PPI_mod (try to create a PPI_mod from non existing pieces)")
            return None


def new_random_PPI_mod(self, P1, P2):
        """ creates new PPI_mod with random parameters and output types"""
        
        list=[['Degradable', mutation.sample_dictionary_ranges('Species.degradation',self.Random) ],['Complex'],['Phosphorylable'],['Complexable']]
        if (P1.isinstance('TF') or P2.isinstance('TF')):
            activity=int(self.Random.random()*2)
            list.append(['TF',activity])
        if (P1.isinstance('Kinase') or P2.isinstance('Kinase')) and ((not P1.isinstance('Phosphatase')) or (not P2.isinstance('Phosphatase'))):
            list.append(['Kinase'])
        if (P1.isinstance('Phosphatase') or P2.isinstance('Phosphatase')) and ((not P1.isinstance('Kinase')) or (not P2.isinstance('Kinase'))):
            list.append(['Phosphatase'])
        if (P1.isinstance('Linear_Producer') or P2.isinstance('Linear_Producer')):
            list.append(['Linear_Producer'])
        if ((P1.common == 1) and (P2.common == 1)):
            list.append(['Common',1])
        else:
            list.append(['Common',0])
        a = mutation.sample_dictionary_ranges('PPI_mod.association',self.Random)
        d = mutation.sample_dictionary_ranges('PPI_mod.disassociation',self.Random)
        [PPI_mod,C]=self.new_PPI_mod(P1,P2,a,d,list)
        return [PPI_mod,C]  


setattr(mutation.Mutable_Network,'random_PPI_mod',random_PPI_mod)
setattr(mutation.Mutable_Network,'new_random_PPI_mod',new_random_PPI_mod)


###########Integration C Tools ##################



def PPI_mod_deriv_inC(net):
 func="\n/**************Complexation interactions*****************/\n"
 if ('PPI_mod' in net.list_types):
  for index in net.list_types['PPI_mod']:
     C=net.graph.successors(index)[0]#finds the complex
     [P1,P2]=net.graph.predecessors(index) #find the components
     
     # Different cases depending on which species is common
     if not (P1.common or P2.common) :	# in the case both the complexable species are of the same type.
         arate="%.10f*"%index.association+P1.id+"*"+P2.id
         drate="%.10f*"%index.disassociation+C.id
         func=func+deriv2.compute_leap([P1.id,P2.id],[C.id],arate)
         func=func+deriv2.compute_leap([C.id],[P1.id,P2.id],drate)
     else: 	# this is the interesting case where one of the species is common and the other specific.
         if P1.common:
             arate_P1="%.10f*"%index.association+P1.id+"*"+P2.id+"*N_cell_"
             drate_P1="%.10f*"%index.disassociation+C.id+"*N_cell_"
             func=func+deriv2.compute_leap([P1.id],[],arate_P1)
             func=func+deriv2.compute_leap([],[P1.id],drate_P1)
             
             arate="%.10f*"%index.association+P1.id+"*"+P2.id
             drate="%.10f*"%index.disassociation+C.id
             func=func+deriv2.compute_leap([P2.id],[C.id],arate)
             func=func+deriv2.compute_leap([C.id],[P2.id],drate)
         else:
             arate_P2="%.10f*"%index.association+P1.id+"*"+P2.id+"*N_cell_"
             drate_P2="%.10f*"%index.disassociation+C.id+"*N_cell_"
             func=func+deriv2.compute_leap([P2.id],[],arate_P2)
             func=func+deriv2.compute_leap([],[P2.id],drate_P2)
             
             arate="%.10f*"%index.association+P1.id+"*"+P2.id
             drate="%.10f*"%index.disassociation+C.id
             func=func+deriv2.compute_leap([P1.id],[C.id],arate)
             func=func+deriv2.compute_leap([C.id],[P1.id],drate)
 return func

deriv2.PPI_mod_deriv_inC=PPI_mod_deriv_inC

