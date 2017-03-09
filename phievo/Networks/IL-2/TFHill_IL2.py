import classes_eds2
import mutation


import copy
from TFHill import *
import config
exec('import '+config.name_deriv2+' as deriv2')  


print("Importing TFHill_IL2")    

mutation.dictionary_ranges['TFHill.hill']=5.0
mutation.dictionary_ranges['TFHill.threshold']=mutation.C



#This class rewrites C integration tools of TFHill to make them compativle with RK integrator 
# and accounts for release of "common" proteins in the environment.

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

