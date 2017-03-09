

from phievo.Networks.CorePromoter import *
from phievo.Networks.TFHill import *
from phievo.Networks.PPI import *
from phievo.Networks.LR import *
from phievo.Networks.Phosphorylation import *
# interactions specific to modelling of the TC-pMHC chemical pathway.
from phievo.Networks.Immune.Simple_Phosphorylation import *
from phievo.Networks.Immune.Simple_Dephosphorylation import *
from phievo.Networks.Immune.KPR_Binding import *
from phievo.Networks.Immune.KPR_Unbinding import *
from phievo.Networks.Immune.Initial_Concentration import *
from phievo.Networks.Immune.Mutable_Threshold import *
from phievo.initialization_code import *
exec('import '+name_deriv2+' as deriv2')

# overwritting the compute_deriv_inC function to be found in deriv2.py. 
# This function essentially writes the system of differential equation
# (appropriate for numerical integration) of the corresponding network.

# added the variable tau at the end so that we can vary the 
def compute_deriv_inC(net):
    func="void derivC(double s[], double ds[], double tau, double dt){\n int index;"
    func=func+"\t for (index=0;index<SIZE;index++) ds[index]=0;//initialization\n"
    func=func+"\t double increment=0;\n"
    func=func+"\t double rate=0;\n"
    #func=func+deriv2_pMHC.degrad_deriv_inC(net)#add degradation rates # Nothing is degradable here.
    func=func+deriv2.SimplePhospho_deriv_inC(net)
    func=func+deriv2.SimpleDephospho_deriv_inC(net)
    func=func+deriv2.compute_KPR_Binding(net)
    func=func+deriv2.compute_KPR_Unbinding(net)
    func=func+"}\n\n"
    return func

deriv2.compute_deriv_inC=compute_deriv_inC
