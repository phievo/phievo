"""Convert a network defined in python to set of C expressions that can be compiled and run
   The C-code goes into workplace_dir/built_integrator*.c along with executable
   The C-code is assembled from .... (source)
	all numerical parameters	all_params2C()
		Note constants encoded as compiler directives, so that they can dimensionalize
		stored arrays input generically in the header.h file
	header.h	fixed file read from /ccode_dir
	utilities.c	fixed file read from /ccode_dir
	all translated network fns	compute_deriv_inC() and others
	appl specific C		read from /ccode_dir, must have one program from each of
		several categories (defined in compute_program()).  The actual instances
		defined by dictionary cfile in initialization.  Skip file type by type = '' or ' '
		or omitting from cfile{} dictionary.

   All these pieces assembled by compute_program(), and then compiled etc by
   compile_and_integrate().

   The c-code files passed only once in form of dictionary cfile.  The numerical parameters
   need to find dimensions of arrays, integration steps, input as argments to functions

   TODO:  it would be nice to include in header.h declaration of all C functions used
   so that they can then be loaded in any order, currently order constrained by declare
   before use.

"""
print("importing deriv2_IL2")
from classes_eds2 import *

from math import sqrt

import numpy

import os,sys,select,random

import subprocess

# default location of directory where build_integrator*.c will go
workplace_dir = './Workplace/'

# default compiler
# Ccompiler = '/usr/local/intel/cc/10.0.023/bin/icc'
Ccompiler = 'cc'

# cfile defines where the generic c-code is found, and names of all files.  These files can be
# reset to problem dependent ones.

cfile = {}  # need define before fillin initialization


noise_flag=0 #flag to know if we integrate or not with noise
"""  information input via run_evolution.py and initialization file
ccode_dir = './CCodes/'
cfile['core'] = ccode_dir + 'integrator.core.c'
etc
"""


def compute_leap(list_input_id,list_output_id,rate):
        """Routine to compute strings for derivative in C associated to an interaction. list_input_id contains id of the input, i.e. the depleted species, list_output_id contains id of the created species. Rate is a string containing the rate, should be positive. If flag is 1, one computes the derivative with the Tau_leap algorithm, i.e. one adds the Langevin noise term scaled with concentration """
	func=""
	func=func+"\t \t rate="+rate+";\n"
	if (noise_flag==1):
		func=func+"\t \t increment=compute_noisy_increment(rate);\n"
	else:
		func=func+"\t \t increment=rate;\n"


	for id in list_input_id:
		func=func+"\t \t d"+id+"-=increment;\n"
	for id in list_output_id:
		func=func+"\t \t d"+id+"+=increment;\n"

	return func



def degrad_deriv_inC(net):
 """start computing the string encoding the function that compute the degradation rateq from a network
 """
 func="\n/**************degradation rates*****************/\n"
 if 'Degradable' in net.list_types:
     for index in net.list_types['Degradable']:
         rate="%f*" %index.degradation+index.id
         func=func+compute_leap([index.id],[],rate)
         #func=func + "\t d" + index.id + "=-%f*" %index.degradation+index.id + ";\n"
 return func




def track_variable(net,name):
    """For name = Input or Output, return a list of the indices of the species with this type, ordered by
    the attribute values of n_put = 0,1,2...These indices communicated to integrator and fitness function
    (for Output).  This is way of keeping track of fixed IO variables. Assumes attribute n_put

    Use this function only if the output or input are fixed in the algorithm, otherwise, use track_changing_variable
    """
    if name not in net.list_types:
        return []

    track = {}
    for s in net.list_types[name]:
        track[s.n_put] = s.int_id()

    # verify that the n_put attributes on the IO variables are numbered consecutively from 0
    track_list = []
    for ii in range( len(net.list_types[name]) ):
        try:
	    track_list.append( track[ii] )
	except:
	    print("mapping IO variable number of id in deriv2.track_variable() failed for IO=", name)
	    for s in net.list_types[name]:
	        s.print_node()
	    return
    return track_list


def track_changing_variable(net,name):
    """For name = Input or Output, return a list of the indices of the species with this type. These indices communicated to integrator and fitness function
    (for Output)

    Use this function when Output or input may be added ( we do not care about their order)
    """
    if name not in net.list_types:
        return []

    track_list = []
    try:
	    for s in net.list_types[name]:
		    track_list.append(s.int_id())
    except:
	    print("mapping IO variable number of id in deriv2.track_variable() failed for IO=", name)
	    for s in net.list_types[name]:
		    s.print_node()
	    return
    return track_list


def all_params2C(net, prmt, print_buf,Cseed=0):
    """ collect all the numerical constants neelocalneig,diff,index_ligand,ded in C code, for a network, net,
        prmt dictionary from initialziation, and boolean print_buf to control printing
	of time history by C codes
    """
    hdr = []   # collect lines of output as list then join, speed issue

    # various sizes/lengths mostly from prmt dict
    size=len(net.list_types['Species'])
    hdr.append("#define SIZE %i" %size )
    #hdr.append("#define NSTEP %i" %prmt['nstep'] )
    #hdr.append("#define NCELLTOT %i" %prmt['ncelltot'] )
    hdr.append("#define NNEIGHBOR %i" %prmt['nneighbor'] )
    hdr.append("#define NTRIES %i" %prmt['ntries'] )
    if 'langevin_noise' in prmt:
	    hdr.append("#define  CONCENTRATION_SCALE %f" %prmt['langevin_noise'] )
    else:
	    hdr.append("#define  CONCENTRATION_SCALE 1.0" )
    # optional generic parameters for specific C subroutines as dict or list.
    # define NFREE_PRMT is flag in Ccode that free_prmt as list is being used
    if 'free_prmt' in prmt:
        if isinstance(prmt['free_prmt'], dict):
	    hdr.append("#define NFREE_PRMT 0" )
            for key,value in list(prmt['free_prmt'].items()):
	        hdr.append("#define %s %s" %( key.upper(), value) )
	elif isinstance(prmt['free_prmt'], list):
	    hdr.append("#define NFREE_PRMT %i" %len(prmt['free_prmt']) )
	    str_pp = ', '.join( [str(nn) for nn in prmt['free_prmt'] ] )
	    hdr.append("static double free_prmt[] = {%s};" %str_pp )
	else:
	    print('unexpected type of prmt[free_prmt] in deriv2.py.  Skipping', prmt['free_prmt'])
    else:
        hdr.append("#define NFREE_PRMT 0" )  # need define dummy length variable incase test in C
	hdr.append("static double free_prmt[] = {}; " )  # dummy defintion for C

    # IO stuff,
    trackin  = track_variable(net, 'Input')
    trackout = track_changing_variable(net, 'Output')  # this routine does not use n_put attribute of Species.
    #Ligand + Diffusible species
    tracklig = track_changing_variable(net, 'Ligand' )
    trackdiff = track_changing_variable(net, 'Diffusible' )
    hdr.append("#define	NOUTPUT %i" %len(trackout) )
    hdr.append("#define	NINPUT %i" %len(trackin) )
    hdr.append("#define	NLIGAND %i" %len(tracklig) )
    hdr.append("#define	NDIFFUSIBLE %i" %len(trackdiff) )
    # misc other numbers
    if (Cseed==0):
        hdr.append("#define SEED %i" %(int(random.random()*1000000)) ) #seed for the C rand
    else:
        hdr.append("#define SEED %i" %Cseed )

    if 'print_buf' in prmt:
        hdr.append("#define PRINT_BUF %i" %prmt['print_buf'] )
    else:
        hdr.append("#define PRINT_BUF %i" %print_buf )

    hdr.append("#define DT %f" %prmt['dt'] )

    # the mapping of input/output indices
    str_in = ', '.join( [str(nn) for nn in trackin] )
    hdr.append("static int trackin[] = {%s};" %str_in )
    str_out = ', '.join( [str(nn) for nn in trackout] )
    hdr.append("static int trackout[] = {%s};" %str_out )
    str_lig = ', '.join( [str(nn) for nn in tracklig] )
    hdr.append("static int tracklig[] = {%s};" %str_lig )

    str_diff = ', '.join( [str(nn) for nn in trackdiff] )
    hdr.append("static int trackdiff[] = {%s};" %str_diff )
    str_diff_constant=', '.join([str(net.list_types['Species'][nn].diffusion) for nn in trackdiff])
    hdr.append("static double diff_constant[] = {%s};" %str_diff_constant )#table containing diffusion constants of ligands
    list_ext=[]
    if 'Ligand' in net.list_types:
	    for nn in net.list_types['Ligand']:
		    if nn.isinstance('Diffusible'):
			    list_ext.append(str(1))
		    else:
			    list_ext.append(str(0))

    str_ext=', '.join(list_ext)
    hdr.append("static int externallig[] = {%s};" %str_ext )#table containing external tags for ligands
    hdr.append("\n\n")

    # Added parameters in the specific present problem
    hdr.append("#define N_PARTITION_OUTPUT %i" %prmt['partition_output'])
    hdr.append("#define NFUNCTIONS %i" %prmt['nfunctions']) 	# number of functions (could be >1 in pareto optimization)

    # parameters specific to the model
    hdr.append("#define TAU_AG %f" %prmt['tau_Ag'])
    hdr.append("#define INV_LIGAND_AFFINITY %f" %prmt['inv_ligand_affinity'])
    hdr.append("#define RECEPTOR %f" %prmt['receptor'])
    #hdr.append("#define ON_RATE_INPUT %f" %prmt['pMHC_on_rate'])
    #hdr.append("#define TAU_INPUT %f" %prmt['pMHC_binding_time'])


    # specifically related to the adaptive integrator.
    hdr.append("#define MAX_TIME %f" %(prmt['max_time']+prmt['equilibration_lag']))
    hdr.append("#define TOL %.15f" %prmt['tolerance'])
    hdr.append("#define BELOW %f" %prmt['scalebelow'])
    hdr.append("#define ABOVE %f" %prmt['scaleabove'])
    hdr.append("#define TINY %.15f" %prmt['tiny'])
    hdr.append("#define EQUILIBRATION_LAG %f" %prmt['equilibration_lag'])

    # Related to concentrations and mechanics of the simulation.
    hdr.append("#define N_N_CELL_ %i" %len(prmt['Number_cells']) ) # Number of different dissociation times.
    hdr.append("#define N_DOSES %i" %len(prmt['Dose_Input']) ) # Number of different ligand concentrations.
    hdr.append("#define NCELLTOT N_N_CELL_*N_DOSES") # Total number of cells given by the number of dissociation times multiplied by the number of ligand concentrations.
    hdr.append("static double DOSES[N_DOSES];")
    hdr.append("static double NUMBER_CELLS[N_N_CELL_];")

    # required for printing the full history (in verifications).
    if 'total_number_steps' in prmt:
        hdr.append("#define NSTEPTOT %i" %prmt['total_number_steps'])
    if 'print_loss' in prmt:
        hdr.append("#define PRINT_LOSS %i" %prmt['print_loss'])

    S_list = net.list_types['Species']
    counter = 0
    for i in range(len(S_list)):
        if (S_list[i].isinstance('Kinase') or S_list[i].isinstance('Phosphatase') or S_list[i].isinstance('Linear_Producer')) and not S_list[i].isinstance('Input'):
            out = net.graph.successors(S_list[i])
            for j in range(len(out)):
                if out[j].isinstance('Initial_Concentration'):
                    counter += 1
    counter +=1
    hdr.append("#define IC_LENGTH %i" %counter )
    hdr.append("static int IC_LIST[IC_LENGTH];")
    hdr.append("static double IC_CONC[SIZE];")

    return '\n'.join(hdr)    # note added the \n here between all elements of hdr


def compute_program(net, prmt, print_buf,Cseed=0):
    """Collect python encoded C and the stored files seleced via cfile dictionary and
       combine in the correct order to create the built_integrator as one string
    """

    # these have to be loaded in this order due to implicit type def's
    required_files2 = ['fitness', 'geometry', 'init_history', 'input', 'integrator', 'main' ]
    pgm = []
    pgm.append( all_params2C(net, prmt, print_buf,Cseed) )
    pgm.append( open(cfile['header']).read() )
    pgm.append( open(cfile['utilities']).read() )
    pgm.append( '/***** end of header, begiining of python computed functions ***/\n\n' )

    # The function below declares the various lists required to make the integration (initial conditions, number of cells present, etc.).
    pgm.append('void declare_list(){\n')
    for j in range(len(prmt['Dose_Input'])):
        pgm.append("DOSES[%i] = %f;\n" %(j, prmt['Dose_Input'][j]) )
    for j in range(len(prmt['Number_cells'])):
        pgm.append("NUMBER_CELLS[%i] = %f;\n" %(j, prmt['Number_cells'][j] ))

    # Specifically: the initial concentrations
    S_list = net.list_types['Species']
    counter = 0
    for i in range(len(S_list)):
        #if (S_list[i].isinstance('Kinase') or S_list[i].isinstance('Phosphatase') or S_list[i].isinstance('Linear_Producer')):
        out = net.graph.successors(S_list[i])
        for j in range(len(out)):
            if out[j].isinstance('Initial_Concentration'):
                pgm.append("IC_LIST[%i] = %i;\n" %(counter, S_list[i].int_id()) )
                pgm.append("IC_CONC[%i] = %f;\n" %(S_list[i].int_id(), out[j].conc) )
                counter += 1
    pgm.append("}\n\n")

    pgm.append( compute_deriv_inC(net) ) #+ compute_LR(net)
    #pgm.append( input_deriv_inC(net) )
    pgm.append( '/***** end of python computed functions, beginning problem specific fns ***/\n\n')
    for file in required_files2:
        if( file in cfile and cfile[file].endswith('.c') ):  # omit files = ' ' etc
	        exec("string = open(cfile[ '%s' ]).read()" %file)
		pgm.append( string )
		# err_3tuple = sys.exc_info() this works only with try: except: and
		# in except clause this command will grap error and traceback and can test and exit.
		# Just try: except will loose where the error occured

    return ' '.join( pgm )

def compile_and_integrate(net, prmt, nnetwork, bool,Cseed=0):
    """Compile and integrate a network=net, with prmt parameter dictionary,
    bool (print history to file or not), and give Ccode an ID of nnetwork (str or int)
    Run the code in 'workplace_dir defined top of this file.
    Wait for process completion before launching another integration
    See http://pydoc.org/2.4.1/subprocess.htm for interface to run C code
    """
    net.write_id()

    # check for outputs, if none return huge fitness
    if 'Output' in net.list_types:
        # if .c files missing this subroutine will fail and we die
        program=compute_program(net,prmt,bool,Cseed)
	integrator_ = workplace_dir + 'built_integrator' + str(nnetwork)
	open(integrator_ + '.c','w').write(program)
	string = Ccompiler + ' ' + integrator_ + '.c -lm -o ' +  integrator_

       	# See http://pydoc.org/2.4.1/subprocess.htm for interface to run C code
	out = subprocess.Popen(string, shell=True, stderr=subprocess.PIPE).communicate()
	# print 'first out', out
	if out[1]:
            print('bug in Ccompile for',integrator_, 'err=', out[1], 'BYE')
	    sys.exit(1)

	string = 'chmod +x ' + integrator_ + '\n' + integrator_ + '\n'
	out = subprocess.Popen(string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	# print 'second out', out
	if out[1] or len(out[0])<1:    # some floating exceptions do not get to stderr, but loose stdout
            print('bug during run (or no stdout) for', integrator_, out[1], 'BYE')
	    sys.exit(1)
	else:
            out_list = out[0].strip().split('\n')
	    if out_list:
                return out_list
	    else:
                return [str(1000000000)]

    else:
        print("No Output for network %i"%nnetwork)
	return [str(900000),str(900000),str(900000)]
