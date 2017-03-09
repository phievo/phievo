
print("importing interaction_IL2")
from CorePromoter import *
from PPI_mod import *
#from Phospho_Dephospho import *
#from Linear_Production import *
from Initial_Concentration import *
from TFHill_IL2 import *
from Absolute_discrimination import *
import config
exec('import '+config.name_deriv2+' as deriv2')

'''Note the ligand dynamics is hard-encoded in the C integrator for now adaptive_RK4_integrator.c'''



def compute_deriv_inC(net):
 func="void derivC(double s[],double ds[], double N_cell_, double dt){\n int index;"
 func=func+"\t for (index=0;index<SIZE;index++) ds[index]=0;//initialization\n"
 func=func+"\t double increment=0;\n"
 func=func+"\t double rate=0;\n"
 func=func+deriv2.degrad_deriv_inC(net)#add degradation rates
 func=func+deriv2.transcription_deriv_inC(net)#add transcription rates
 func=func+deriv2.PPI_mod_deriv_inC(net)#add PPIs
 func=func+deriv2.absolute_discrimination_deriv_inC(net)
 #func=func+deriv2.Linear_Production_deriv_inC(net)#adds linear production
 #func=func+deriv2.Simple_Phosphorylation_deriv_inC(net)# adds linear phosphorylation
 func=func+"\n\n}"
 return func


deriv2.compute_deriv_inC=compute_deriv_inC
