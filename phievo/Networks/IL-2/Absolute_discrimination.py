import classes_eds2
import mutation
import config
exec('import '+config.name_deriv2+' as deriv2')


import copy

print("Importing Absolute discrimination")        



class Absolute_discrimination(classes_eds2.Interaction):
    """Couple TF to promoter, init: hill coef, threshold """

    def __init__(self,LR,tf, N, phitau):
        classes_eds2.Node.__init__(self)
        #self.hill=hill
        #self.threshold=threshold
        #self.activity=activity
        self.N=N
        self.phitau=phitau
        self.label = 'Absolute_discrimination'
        self.input=['Species']
        self.output=['TF']
    
    def string_param(self):
    	return "Absolute_discrimination"

        


####### Attributes attached to Network for Absolute_discrimination ######################################

            
        
def new_Absolute(self,LR,tf, N,phitau):
        """Create a new Absolute discrimination black_box between LR and TF"""
 
        
        r = Absolute_discrimination(LR,tf,N,phitau)
        self.add_Node(r)
        self.graph.add_edge(LR,r)
        self.graph.add_edge(r,tf)        
        return r
               





setattr(classes_eds2.Network,'new_Absolute',new_Absolute)




def absolute_discrimination_deriv_inC(net):
 """computing the string encoding the function that compute the transcription rate from a network
 """
 func="\n/**************Absolute Discrimination*****************/\n"
 net.write_id()
 if ('Absolute_discrimination' in net.list_types):
  for absolute in net.list_types['Absolute_discrimination']:
   if isinstance(absolute,Absolute_discrimination):
    output=net.graph.successors(absolute)[0]    #find the Output
    index_LR=net.graph.predecessors(absolute)[0].int_id()   #find the LR complex
    AS="POW(%f,%f)*s[%i]"%(absolute.phitau,absolute.N,index_LR)
    print(AS)
    func=func+deriv2.compute_leap([],[output.id],AS)
    func=func+deriv2.compute_leap([output.id],[],output.id)
 return func


deriv2.absolute_discrimination_deriv_inC=absolute_discrimination_deriv_inC
