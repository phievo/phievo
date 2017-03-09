"""
This file import all the different interaction and update the deriv2 module
to incorporate the equation in the integration file.

If you want to add a new interaction, you have to add:
    an import at the beginning
    an add() call in the write_deriv_inC function
"""
print("Execute interaction.py")

#import the different interaction module
from .CorePromoter import *
from .TFHill import *
from .PPI import *
from .LR import *
from .Phosphorylation import *
from .Degradation import *
from . import deriv2

def write_deriv_inC(net,programm_file):
    """Write the integration equation in the C-file
    
    This function is an update from the one in deriv2
    
    Args:
        net (Network): the network under study
        programm_file (TextIOWrapper): the built_integrator file
    Return:
        None: directly write the string in the C-file
    """
    start="void derivC(double s[],double history[][NSTEP][NCELLTOT],int step, double ds[],double memories[],int ncell){\n int index;"

    ## Updates species IDs
    net.write_id()
    add = programm_file.write #create a bound method for readibility
    add(start)
    add("\t for (index=0;index<SIZE;index++) ds[index]=0;//initialization\n")
    add("\t double increment=0;\n")
    add("\t double rate=0;\n")    
    add(deriv2.degrad_deriv_inC(net))#add degradation rates
    add(deriv2.transcription_deriv_inC(net))#add transcription rates
    add(deriv2.PPI_deriv_inC(net))#add PPIs
    add(deriv2.Phospho_deriv_inC(net))#add phosphorylation
    add(deriv2.Degradation_deriv_inC(net))#add degradation
    add("}\n\n")
    add(deriv2.compute_LR(net))#add LR Interaction (why the hell is it different from all the other ???)

#updates deriv2
deriv2.write_deriv_inC=write_deriv_inC
