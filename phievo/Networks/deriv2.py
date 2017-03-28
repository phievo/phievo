"""Here are the tools to convert a Network object to a C-file that will
be compiled and run.
The C-file goes into workplace_dir/built_integrator*.c along with executable
The C-file is assembled with several pieces:

    - header, utilities, geometry, integrator and main: see initialization_code.init_deriv2
    - for each interaction: see interaction.interaction_deriv_inC (bottom of file)
    - see also Networks.interaction.py and the cfile dictionary

All these pieces are assembled by compute_program(), and then compiled with
compile_and_integrate().

The c-code files passed only once in form of dictionary cfile.  The numerical parameters
need to find dimensions of arrays, integration steps, input as argments to functions

Attributes:
    workplace_dir (str): the directory where build_integrator*.c will go
    Ccompiler (str): 'gcc' by default
    cfile (dict): where the generic c-code are found (can be reset to fit problem)
    noise_flag (bool): flag to know if we integrate or not with noise

TODO:  it would be nice to include in header.h declaration of all C functions used
so that they can then be loaded in any order, currently order constrained by declare
before use.
"""
from phievo import __silent__,__verbose__
if __verbose__:
    print("Execute deriv2")

from phievo.initialization_code import display_error
from phievo.Networks.classes_eds2 import *
from math import sqrt
import numpy
import os, sys, select, random
import subprocess

# Parameters
workplace_dir = './Workplace/'
Ccompiler = 'gcc'
cfile = {}  # see initialization_code.init_deriv2 for the whole definition
noise_flag = False

########## Routine Functions ##########

def compute_leap(list_input_id, list_output_id, rate):
    """Routine to compute strings for derivative in C associated to an interaction

    if noise_flag, adds a Langevin noise term which scaled with concentration

    Args:
        list_input_id (list): contains id of the input, i.e. the depleted species
        list_output_id (list): contains id of the created species
        rate (str): the rate, should be positive

    Return:
        str: a C-formatted string
    """
    func = "\t \t rate=" + rate + ";\n"
    func += "\t \t increment=compute_noisy_increment(rate);\n" if noise_flag else "\t \t increment=rate;\n"
    func += ''.join("\t \t d" + id + "-=increment;\n" for id in list_input_id)
    func += ''.join("\t \t d" + id + "+=increment;\n" for id in list_output_id)

    return func

def track_variable(net, name):
    """Return a list of the indices of the species with type name

    This is way of keeping track of fixed IO variables.
    Use this function only if the output or input are fixed in the
    algorithm, otherwise, use track_changing_variable

    Args:
        net (Networks): -
        name (str): a Species tag, usually 'Input' or 'Output'

    Return:
        list: the id species list ordered by growing n_put
    """
    if name not in net.list_types:
        return []
    track = {s.n_put:s.int_id() for s in net.list_types[name]}
    # verify that the n_put attributes on the IO variables are numbered consecutively from 0
    track_list = []
    for ii in range(len(net.list_types[name])):
        try:
            track_list.append(track[ii])
        except Exception:
            display_error("mapping IO variable number of id in deriv2.track_variable() failed for IO="+str(name))
            for s in net.list_types[name]:
                s.print_node()
            return
    return track_list

def track_changing_variable(net, name):
    """Return a list of the indices of the species with type name

    Use this function when Output or Input may be added
    (we do not care about their order)

    Args:
        net (Networks): -
        name (str): a Species tag, usually 'Input' or 'Output'

    Return:
        list: the id species list ordered by growing n_put
    """
    return [s.int_id() for s in net.list_types.get(name,[])]

########## Writing Functions ##########
# Here are the functions which explicitely construct the C-file

def degrad_deriv_inC(net):
    """gives the string corresponding to the degradation integration

    Return:
        str: a single string for all degradations in the network
    """
    if 'Degradable' in net.list_types:
        func = "\n/**************degradation rates*****************/\n"
        for species in net.list_types['Degradable']:
            rate = '{0}*{1}'.format(species.degradation,species.id)
            func += compute_leap([species.id], [], rate)
        return func
    else:
        return "\n"

def write_deriv_inC(net,programm_file):
    """Write the integration equations in the C-file

    This function is a default and should be updated in Networks/interaction.py

    Args:
        net (Network): the network under study
        programm_file (TextIOWrapper): the built_integrator file
    Return:
        None: directly write the string in the C-file
    """
    start="void derivC(double s[],double history[][NSTEP][NCELLTOT],int step, double ds[],double memories[],int ncell){\n int index;"
    add = programm_file.write #create a bound method for readibility
    add(start)
    add("\t for (index=0;index<SIZE;index++) ds[index]=0;//initialization\n")
    add("\t double increment=0;\n")
    add("\t double rate=0;\n")
    add(deriv2.degrad_deriv_inC(net))#add degradation rates
    add("}\n\n")

def all_params2C(net, prmt, print_buf, Cseed=0):
    """ Collect all the numerical constants and format them to C like

    neelocalneig,diff,index_ligand,ded

    Args:
        net (Network): -
        prmt (dict): dictionary from initialization file
        print_buf (bool): control printing of time history by C codes
        Cseed (int): seed for the integrator random number generator

    Return:
        str: a C formated string of parameters
    """
    hdr = [] # collect lines of output as list then join, speed issue

    # various sizes/lengths mostly from prmt dict
    hdr.append("#define SIZE %i" % len(net.list_types['Species']))
    hdr.append("#define NSTEP %i" % prmt['nstep'])
    hdr.append("#define NCELLTOT %i" % prmt['ncelltot'])
    hdr.append("#define NNEIGHBOR %i" % prmt['nneighbor'])
    hdr.append("#define NTRIES %i" % prmt['ntries'])
    if 'langevin_noise' in prmt:
        hdr.append("#define  CONCENTRATION_SCALE %f" % prmt['langevin_noise'])
    else:
        hdr.append("#define  CONCENTRATION_SCALE 1.0")
    # optional generic parameters for specific C subroutines as dict or list.
    # define NFREE_PRMT is flag in Ccode that free_prmt as list is being used
    if 'free_prmt' in prmt:
        if isinstance(prmt['free_prmt'], dict):
            hdr.append("#define NFREE_PRMT 0")
            for key, value in list(prmt['free_prmt'].items()):
                hdr.append("#define %s %s" % ( key.upper(), value))
        elif isinstance(prmt['free_prmt'], list):
            hdr.append("#define NFREE_PRMT %i" % len(prmt['free_prmt']))
            str_pp = ', '.join([str(nn) for nn in prmt['free_prmt']])
            hdr.append("static double free_prmt[] = {%s};" % str_pp)
        else:
            print('unexpected type of prmt[free_prmt] in deriv2.py.  Skipping', prmt['free_prmt'])
    else:
        hdr.append("#define NFREE_PRMT 0")  # need define dummy length variable incase test in C
        hdr.append("static double free_prmt[] = {}; ")  # dummy defintion for C

    # IO stuff,
    trackin = track_variable(net, 'Input')
    trackout = track_changing_variable(net, 'Output')  # this routine does not use n_put attribute of Species.
    # Ligand + Diffusible species
    tracklig = track_changing_variable(net, 'Ligand')
    trackdiff = track_changing_variable(net, 'Diffusible')
    #print print_Network(net)
    hdr.append("#define	NOUTPUT %i" % len(trackout))
    hdr.append("#define	NINPUT %i" % len(trackin))
    hdr.append("#define	NLIGAND %i" % len(tracklig))
    hdr.append("#define	NDIFFUSIBLE %i" % len(trackdiff))
    # misc other numbers
    if (Cseed == 0):
        hdr.append("#define SEED %i" % (int(random.random() * 1000000)))  #seed for the C rand
    else:
        hdr.append("#define SEED %i" % Cseed)

    hdr.append("#define PRINT_BUF %i" % print_buf)
    hdr.append("#define DT %f" % prmt['dt'])

    # the mapping of input/output indices
    str_in = ', '.join([str(nn) for nn in trackin])
    hdr.append("static int trackin[] = {%s};" % str_in)
    str_out = ', '.join([str(nn) for nn in trackout])
    hdr.append("static int trackout[] = {%s};" % str_out)
    str_lig = ', '.join([str(nn) for nn in tracklig])
    hdr.append("static int tracklig[] = {%s};" % str_lig)

    str_diff = ', '.join([str(nn) for nn in trackdiff])
    hdr.append("static int trackdiff[] = {%s};" % str_diff)
    str_diff_constant = ', '.join([str(net.list_types['Species'][nn].diffusion) for nn in trackdiff])
    hdr.append(
        "static double diff_constant[] = {%s};" % str_diff_constant)  #table containing diffusion constants of ligands
    list_ext = []
    if 'Ligand' in net.list_types:
        for nn in net.list_types['Ligand']:
            if nn.isinstance('Diffusible'):
                list_ext.append(str(1))
            else:
                list_ext.append(str(0))

    str_ext = ', '.join(list_ext)
    hdr.append("static int externallig[] = {%s};" % str_ext)  #table containing external tags for ligands
    hdr.append("\n\n")

    return '\n'.join(hdr)  # note added the \n here between all elements of hdr

def write_program(programm_file,net, prmt, print_buf, Cseed=0):
    """Write the built_integrator of the network in the C file

    Collect python encoded C and the stored files selected via cfile
    dictionary and write them in the correct order.

    Args:
        programm_file (TextIOWrapper): the built_integrator file
        net (Network): -
        prmt (dict): passed to all_params2C
        print_buf (bool): passed to all_params2C
        Cseed (int): passed to all_params2C

    Return:
        str: the C programm as a python string
    """
    # these have to be loaded in this order due to implicit type def's
    required_files2 = ['fitness', 'geometry', 'init_history', 'input', 'integrator', 'main']
    programm_file.write(all_params2C(net, prmt, print_buf, Cseed))
    programm_file.write(open(cfile['header']).read())
    programm_file.write(open(cfile['utilities']).read())
    programm_file.write('/***** end of header, begining of python computed functions ***/\n\n')
    write_deriv_inC(net,programm_file) #define in Networks/interaction.py
    programm_file.write('/***** end of python computed functions, beginning problem specific fns ***/\n\n')
    for file_name in required_files2:
        if file_name in cfile and cfile[file_name].endswith('.c'):  # omit files = ' ' etc
            programm_file.write(open(cfile[file_name]).read())

########## Program Functions ##########

def compile_and_integrate(network, prmt, nnetwork, print_buf=False, Cseed=0):
    """Compile and integrate a network

    Run the code in 'workplace_dir defined top of this file.
    Wait for process completion before launching another integration
    See https://www.python.org/dev/peps/pep-0324/ for interface to run C code

    Args:
        network (Network): -
        prmt (dict): dictionary from initialization file
        nnetwork (int): an id to separate the different C-file
        print_buf (bool): control printing of time history by C codes to a file
        Cseed (int): seed for the integrator random number generator

    Return:
        list: corresponding to the different line of the output of treatment_fitness
        (see your fitness.c file) or None if an error occured
    """
    network.write_id()
    cfile_directory = workplace_dir+'built_integrator'+str(nnetwork)

    # check for outputs
    if 'Output' not in network.list_types:
        print("No Output for network %i" % nnetwork)
        return None
    # Write the program in a c file
    with open(cfile_directory+'.c','w') as cfile:
        write_program(cfile,network, prmt, print_buf, Cseed)

    # Compile the program
    # cmd contains the command in the same order as they would be on a full bash commans
    # ex: cmd = ["gcc", "-o", "run",  "test.c"] for "gcc -o run test.c"
    cmd = [Ccompiler , cfile_directory+".c" , "-lm" , "-o" , cfile_directory]
    out = subprocess.Popen(cmd,stderr=subprocess.PIPE).communicate()


    if out[1]:
        print('bug in Ccompile for', cfile_directory, 'err=', out[1], 'BYE')
        sys.exit(1)

    # Execute the programm
    out = subprocess.Popen(cfile_directory, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

    if out[1] or len(out[0]) < 1:  # some floating exceptions do not get to stderr, but loose stdout
        print('bug during run (or no stdout) for', cfile_directory, out[1], 'BYE')
        sys.exit(1)
    else:
        out_str = out[0].strip().decode()
        out_list = out_str.split('\n')
        for arg in out_list:
            if not arg: return None
        return out_list
