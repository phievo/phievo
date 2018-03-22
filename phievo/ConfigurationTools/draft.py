from phievo.Networks import mutation
from phievo.Networks import initialization
from phievo.initialization_code import ccode_dir
import os

PRINT = False
pfile = {"deriv2" : "phievo.Networks.deriv2",
             "interaction" : "phievo.Networks.interaction",
             "pretty_graph": "phievo.Networks.lovelyGraph",
             "plotdata" : "phievo.Networks.plotdata"}

cfile = {
    "header" : os.path.join(ccode_dir,'integrator_header.h'),
    "utilities" : os.path.join(ccode_dir,'utilities.c'),
    "geometry" : os.path.join(ccode_dir,'linear_geometry.c'),
    "integrator" : os.path.join(ccode_dir,'euler_integrator.c'),
    "main" : os.path.join(ccode_dir,'main_general.c')
}

prmt = initialization.prmt
restart = prmt.pop("restart")

configurations = {
    "dictionary_ranges":mutation.dictionary_ranges,
    "dictionary_mutation":initialization.dictionary_mutation,
    "cfile":cfile,
    "pfile":pfile,
    "prmt":prmt,
    "restart":restart
}


if PRINT:
    print("Kinetic parameters:")
    for key,val in mutation.dictionary_ranges.items():
        print("\t",key,":",val)

    print("\nMutation parameters:")
    for key,val in initialization.dictionary_mutation.items():
        print("\t",key,":",val)    

    print("\npfile")
    for key,val in pfile.items():
        print("\t",key,":",val)

    print("\ncfile")
    for key,val in cfile.items():
        print("\t",key,":",val)

    print("\nprmt")
    for key,val in prmt.items():
        print("\t",key,":",val)

    print("\nrestart")
    for key,val in restart.items():
        print("\t",key,":",val)
