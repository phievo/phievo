"""
This file import all the different interaction and update the deriv2 module
to incorporate the equation in the integration file.

If you want to add a new interaction, you have to add:
    an import at the beginning
    an add() call in the write_deriv_inC function
"""
from phievo import __silent__,__verbose__
if __verbose__:
    print("Execute interaction.py")

#import the different interaction module
from .CorePromoter import *
from .TFHill import *
from .PPI import *
from .Phosphorylation import *
from .Degradation import *
from . import deriv2

def write_deriv_inC(net,programm_file):
    """Write the integration equation in the C-file

    This function is an update from the one in deriv2

    Args:
        net (:class:`Mutable_Network <phievo.Networks.mutation.Mutable_Network>`): the network under study
        programm_file (TextIOWrapper): the built_integrator file
    """
    start="void derivC(double s[],double history[][NSTEP][NCELLTOT],int step, double ds[],double memories[],int ncell){\n int index;"

    ## Updates species IDs
    net.write_id()
    add = programm_file.write #create a bound method for readibility
    add(start)
    add("\t for (index=0;index<SIZE;index++) ds[index]=0;//initialization\n")
    add("\t double increment=0;\n")
    add("\t double rate=0;\n")
    for ind,deriv_inC in deriv2.interactions_deriv_inC.items():
        add(deriv_inC(net))
    add("}\n\n")

#updates deriv2
deriv2.write_deriv_inC=write_deriv_inC
