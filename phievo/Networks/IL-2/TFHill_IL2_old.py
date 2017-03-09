import classes_eds2
import mutation


import copy

import config
exec('import '+config.name_deriv2+' as deriv2')  


print("Importing TFHill_IL2")    

mutation.dictionary_ranges['TFHill.hill']=5.0
mutation.dictionary_ranges['TFHill.threshold']=mutation.C


class TFHill(classes_eds2.Interaction):
    """Couple TF to promoter, init: hill coef, threshold """

    def __init__(self,hill=0,threshold=0,activity=0):
        classes_eds2.Node.__init__(self)
        self.hill=hill
        self.threshold=threshold
        self.activity=activity
        self.label = 'TF Hill Fn'
        self.input=['TF']
        self.output=['TModule']

        


####### Attributes attached to Network for TFHill ######################################

# The transcription related interaction nodes are added along with edge connecting
# it to species, either up or down stream of it.


        
def add_TFHill(self, tf, inter, module):  # name logically TFHill2TModule ??
        """Add a TF to TModule interaction and two unique edges to graph"""
        
        if inter.isinstance('TFHill') * inter.check_grammar([tf], [module]):
            self.add_Node(tf)
            self.add_Node(inter)
            self.add_Node(module)
            self.graph.add_edge(inter, module)
            self.graph.add_edge(tf, inter)
        else:
            print("Error in grammar add_TFHill")

                  



def number_TFHill(self):
        """return the number of possible TFHIll, requires current list_types
        """
        n_TFHill=self.number_nodes('TFHill')
        n_TF=self.number_nodes('TF')
        n_TM=self.number_nodes('TModule')           
        return n_TF*n_TM-n_TFHill
       
def propagate_activity_TFHill(self):
        """copies the activity of the predecessors of the TFHill into the current TFHill - done for compatibility with older versions
        """
        self.build_list_types()
        if (self.fixed_activity_for_TF==1):
            for tfh in self.list_types['TFHill']:
                tf=self.graph.predecessors(tfh)
                tfh.activity=tf[0].activity
            
        
def new_TFHill(self, tf, hill, threshold, module,activity=0):
        """Create a new TFHill with parameters(hill, thresh) and add 2 links
        args: TF(object), hill, thresh (float) TModule(object)"""
        if (self.fixed_activity_for_TF==1):
            activity_tfh=tf.activity
        else:
            activity_tfh=activity
        
        r = TFHill(hill, threshold,activity_tfh)
        if r.check_grammar([tf], [module]):
            self.add_TFHill(tf, r, module)
        else:
            print("Error in grammar, new_TFHill")
            r = None

        return r
               

def duplicate_TFHill(self,D_species,interaction,module,D_module):
        #function to duplicate a TFHill interaction, D_Species is the duplicated input of this TFHill, D_module is the Tmodule that is duplicated in this gene duplication, the original being module
        D_interaction=copy.deepcopy(interaction)
        D_interaction.mutable=1
        D_interaction.removable=True
        self.add_Node(D_interaction)
        self.graph.add_edge(D_species,D_interaction)
        successors=self.graph.successors(interaction)
        self.graph.add_edge(D_interaction,successors[0])#TFHill only one successor so this is OK
        if (successors[0]==module):#specific case for auto regulatory feedback loop, one also needs to plug the TFHill back on the duplicated Tmodule
            D_interaction_2=copy.deepcopy(interaction)
            D_interaction_2.mutable=1
            D_interaction_2.removable=True
            self.add_Node(D_interaction_2)
            self.graph.add_edge(D_species,D_interaction_2)
            self.graph.add_edge(D_interaction_2,D_module)



setattr(classes_eds2.Network,'add_TFHill',add_TFHill)
setattr(classes_eds2.Network,'number_TFHill',number_TFHill)
setattr(classes_eds2.Network,'propagate_activity_TFHill',propagate_activity_TFHill)
setattr(classes_eds2.Network,'new_TFHill',new_TFHill)
setattr(classes_eds2.Network,'duplicate_TFHill',duplicate_TFHill)




####### Attributes attached to Mutable_Network for TFHill ######################################


def random_TFHill(self):
        """Create new random TFHill instance from list of those possible """
        if 'TF' in self.list_types and 'TModule' in self.list_types:
            list_possible_TFHill=[]
            nTF=len(self.list_types['TF'])
            nTM=len(self.list_types['TModule'])
            for iTF in range(nTF):
                for iTM in range(nTM):
                    tf=self.list_types['TF'][iTF]
                    module=self.list_types['TModule'][iTM]
                    if not self.check_existing_link([tf,module],'TFHill'):
                        list_possible_TFHill.append([tf,module])#list of possible interactions
            n_pTFH=len(list_possible_TFHill)
            if not (n_pTFH==self.number_TFHill()):
                print("Potential Bug : Inconsistency in Computation of number of TFHill")
            
            if (n_pTFH==0):
                print("In random_TFHill : No other posible TFHill")
                return None
            else:
                [tf,module]=list_possible_TFHill[int(self.Random.random()*n_pTFH)]
                return self.new_random_TFHill(tf, module)
            return None

def new_random_TFHill(self, tf, module):
       hill = mutation.sample_dictionary_ranges('TFHill.hill',self.Random)
       threshold = mutation.sample_dictionary_ranges('TFHill.threshold',self.Random)
       activity=int(2*self.Random.random())
       return self.new_TFHill(tf, hill, threshold, module,activity) 

setattr(mutation.Mutable_Network,'random_TFHill',random_TFHill)
setattr(mutation.Mutable_Network,'new_random_TFHill',new_random_TFHill)


############## Integration C Tools ################################



def compute_transcription(net,module):
        """Input a TModule object in network and return a string with an algebraic expression
        for the transciption rate.  """
        
        listactivator=[]
        listrepressor=[]
        if isinstance(module,classes_eds2.TModule):
            for index in net.graph.in_edges(module):
                reg=index[0] #detect the corresponding regulations
		current_activity=reg.activity
                if (current_activity==0):
                    listrepressor.append("HillR(s[%i],%f,%f)"%(net.graph.predecessors(reg)[0].int_id(),reg.threshold,reg.hill))
                else:

                    listactivator.append("HillA(s[%i],%f,%f)"%(net.graph.predecessors(reg)[0].int_id(),reg.threshold,reg.hill))
            l=len(listactivator)
            term = ""
            if(l==0):
                term="%f"%module.rate
		if hasattr(net,"activator_required"): #tests if we want to turn one genes by default from version 1.4.2
			if (net.activator_required==1):
				term="0.00"
		
            if (l==1):
                term="%f*"%module.rate+listactivator[0]
            if (l>1):
                term="%f*"%module.rate
                for index in range(l-1):
                    term=term+"MAX("
                term=term+listactivator[0]
                for index in range(l-1):
                    term=term+","+listactivator[index+1]+")"
	    
	    if hasattr(module, "basal"): #tests on the existenc of a basal rate from version 1.3
		    term="MAX("+term+",%f)"%module.basal

	    if (len(listrepressor)>0):
                term=term+"*"+"*".join(listrepressor)

	    
            return term
        else:
            print("Error in ComputeTranscription")



def transcription_deriv_inC(net):
 """computing the string encoding the function that compute the transcription rate from a network
 """
 func="\n/**************Transcription rates*****************/\n"
 net.write_id()
 if ('TModule' in net.list_types):
  for index in net.list_types['TModule']:
   if isinstance(index,classes_eds2.TModule):
    trans=net.graph.successors(index)    #find the CorePromoter
    output=net.graph.successors(trans[0])    #find the transcribed protein
    
    if (output[0].common == 0) :
        func=func+deriv2.compute_leap([],[output[0].id],compute_transcription(net,index))
    else:
        rate = "N_cell_*"+compute_transcription(net,index)
        func=func+deriv2.compute_leap([],[output[0].id],rate)
 return func

deriv2.compute_transcription=compute_transcription
deriv2.transcription_deriv_inC=transcription_deriv_inC

# We disabled any delay.
#func=func+"\t memory=dt-%i;\n"%trans[0].delay #trans[0].delay must be an integer
#func=func+"\t if(memory>=0){\n"