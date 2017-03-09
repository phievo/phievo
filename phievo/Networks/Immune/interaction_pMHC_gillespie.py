from CorePromoter import *
from TFHill import *
from PPI import *
from LR import *
from Phosphorylation import *

# interactions specific to modelling of the TC-pMHC chemical pathway.
from Simple_Phosphorylation import *
from Simple_Dephosphorylation import *
from KPR_Binding import *
from KPR_Unbinding import *
from Initial_Concentration import *
from Mutable_Threshold import *
import gillespie_pMHC



# overwritting the compute_gillespie_inC function to be found in gillespie_pMHC. 
# This function essentially writes the system of differential equation
# (appropriate for numerical integration) of the corresponding network.

# This makes the code more modular.

# added the variable tau at the end so that we can vary the 
def compute_gillespie_inC(net,prmt):

 prog_proba="void compute_probabilities(double s[][NCELLTOT], double tau, int ncell){\n"
 prog_action="void update_state(double s[][NCELLTOT], int index_action, int ncell){\n"
 n_reactions=0
 
 [proba,action,n_reactions]=gillespie_pMHC.compute_gillespie_Simple_Phosphorylation(net,n_reactions)  # add Simple_Phosphorylations
 prog_proba=prog_proba+proba
 prog_action=prog_action+action

 [proba,action,n_reactions]=gillespie_pMHC.compute_gillespie_Simple_Dephosphorylation(net,n_reactions)  # add Simple_Dephosphorylations
 prog_proba=prog_proba+proba
 prog_action=prog_action+action

 [proba,action,n_reactions]=gillespie_pMHC.compute_gillespie_KPR_Binding(net,n_reactions)  # add KPR_Binding
 prog_proba=prog_proba+proba
 prog_action=prog_action+action

 [proba,action,n_reactions]=gillespie_pMHC.compute_gillespie_KPR_Unbinding(net,n_reactions)  # add KPR_Unbinding
 prog_proba=prog_proba+proba
 prog_action=prog_action+action

 prog_proba=prog_proba+"\n}\n\n"
 prog_action=prog_action+"\n}\n\n"
 
 return [prog_proba,prog_action]

# overwriting the function from gillespie_pMHC.
gillespie_pMHC.compute_gillespie_inC=compute_gillespie_inC
