"""
This module contains the default for the creation of a mutable network.
All of the present parameter can be changed in the initialization file present in
the project directory.
This module also contains default init_network and fitness_treatment functions.
"""

# ranges over which parmeters in classes_eds2 classes can vary. Defined in mutation.py
# NB use 1.0 etc to keep params real not int.

T=1.0 #typical time scale
C=1.0 #typical concentration
L=1.0 #typical size for diffusion

#dictionary key is the Class.attribute whose range is give an a real number or as a LIST
# indicating the interval over which it has to vary

dictionary_ranges={}
dictionary_ranges['Species.degradation']=1.0/T
dictionary_ranges['Species.diffusion']=L   # for ligands diffusion
dictionary_ranges['TModule.rate']=C/T
dictionary_ranges['TModule.basal']=0.0
dictionary_ranges['CorePromoter.delay']=0   # convert to int(T/dt) in run_evolution.py
dictionary_ranges['TFHill.hill']=5.0
dictionary_ranges['TFHill.threshold']=C
dictionary_ranges['PPI.association']=1.0/(C*T)
dictionary_ranges['PPI.disassociation']=1.0/T
dictionary_ranges['Phosphorylation.rate']=1.0/T
dictionary_ranges['Phosphorylation.hill']=5.0
dictionary_ranges['Phosphorylation.threshold']=C
dictionary_ranges['Phosphorylation.dephosphorylation']=1.0/T


# when key='relative_variation' is present, use this fraction of parameter range listed
# above, to vary parameter for each mutation.  Comment out line to pick parameter
# randomly from range defined above.  See mutation.py for details
#dictionary_ranges['relative_variation']=0.1

# Bias is a parameter in [0,1] that is in effect when 'relative_variation' is set.
# A postive value means rates tend to decrease more than increase and
# binding affinities are less strong.  bias=1 implies one sided change
#dictionary_ranges['bias']=0.1

#################################################################################
# names of c-code files needed by deriv2.py (. for extension only)
# normally need 'header', 'utilities' (loaded before python derived routines and
# ['fitness', 'geometry', 'init_history', 'input', 'integrator', 'main' ] placed afterwards
# skip by setting cfile[] = ' ' or ''

cfile = {}
cfile['header'] = 'integrator_header.h'
cfile['utilities'] = 'utilities.c'
cfile['geometry'] = 'linear_geometry.c'
cfile['integrator'] = 'euler_integrator.c'
cfile['main'] = 'main_general.c'

#################################################################################
# mutation rates
dictionary_mutation={}

# Rates for nodes to add
dictionary_mutation['random_gene()']=0.01
dictionary_mutation['random_gene(\'TF\')']=0.01
dictionary_mutation['random_gene(\'Kinase\')']=0.01
dictionary_mutation['random_gene(\'Ligand\')']=0.00
dictionary_mutation['random_gene(\'Receptor\')']=0.00

dictionary_mutation['random_Interaction(\'TFHill\')']=0.002
dictionary_mutation['random_Interaction(\'PPI\')']=0.002
dictionary_mutation['random_Interaction(\'Phosphorylation\')']=0.002


# Rates for nodes to remove
dictionary_mutation['remove_Interaction(\'TFHill\')']=0.005
dictionary_mutation['remove_Interaction(\'PPI\')']=0.005
dictionary_mutation['remove_Interaction(\'CorePromoter\')']=0.01
dictionary_mutation['remove_Interaction(\'Phosphorylation\')']=0.005

# Rates to change parameters for a node
dictionary_mutation['mutate_Node(\'Species\')']=0.1
dictionary_mutation['mutate_Node(\'TFHill\')']=0.1
dictionary_mutation['mutate_Node(\'CorePromoter\')']=0.1
dictionary_mutation['mutate_Node(\'TModule\')']=0.1
dictionary_mutation['mutate_Node(\'PPI\')']=0.1
dictionary_mutation['mutate_Node(\'Phosphorylation\')']=0.1

#rates to change output tags.  See list_types_output array below
dictionary_mutation['random_add_output()']=0.0
dictionary_mutation['random_remove_output()']=0.0
dictionary_mutation['random_change_output()']=0.1


# duplicates gene and some interactions, see mutation.random_duplicate()
# dictionary_mutation['random_duplicate()']=0.01

#############################################################################
# parameters in various modules, created as one dict, so that can then be passed as argument

# Needed in deriv2.py to create C program from network
prmt = {}
prmt['nstep'] =30000         #number of time steps
prmt['ncelltot']=1           #number of cells in an organism
prmt['nneighbor'] = 3  # must be >0, whatever geometry requires, even for ncelltot=1
prmt['ntries'] = 10    # number of initial conditions tried in C programs
prmt['dt'] = 0.05      # time step
prmt['noutput'] = 1    # to define initial network, recomputed by python from network
prmt['ninput'] = 1     # to define init network, does not change during evolution, computed in python from network
# if >0, turns on langevin integration and value becomes CONCENTRATION_SCALE in CCode. should be ~100 ie # molecules
prmt['langevin_noise'] = 0

# Generic parameters, transmitted to C as list or dictionary.
# dict encoded in C as: #define KEY value (KEY converted to CAPS)
# list encoded in C as: static double free_prmt[] = {}
# The dict is more readable, but no way to test in C if given key supplied.  So
# always include in C: #define NFREE_PRMT int. if(NFREE_PRMT) then use free_prmt[n]
# prmt['free_prmt'] = { 'step_egf_off':20 }  # beware of int vs double type
# prmt['free_prmt'] = [1,2]

# Needed in evolution_gill to define evol algorithm and control output
prmt['npopulation'] =30     # number of species in population
prmt['ngeneration'] =301    # number of generations of evolution (+1 to get stats last step)
prmt['tgeneration']=0.01    # initial generation time (for gillespie), reset during evolution
prmt['freq_plot'] = None     # plot time course every generation%freq_plot = 0
prmt['freq_stat'] = 5       # print stats every freq_stat generations
prmt['frac_mutate'] = 0.5   # fraction of networks to mutate, keep at 1/2
prmt['redo'] = 1            # rerun the networks that do not change to compute fitness for different IC

# used in run_evolution,
prmt['nseed'] = 20		# number of times entire evol procedure repeated, see main program.
prmt['firstseed'] = 0	# first seed, defines random number seed, so can exactly repeat

# Any species with these tags can not be removed.  MUST include 'Input'.  If fixed number of
# outputs add 'Output' to list. Overwrites default list in classes_eds2.py
list_unremovable=['Input','Output']

# Controls the types of outputs, when the output tags can move or be created.  MUST
# be at least 'Species'.  NB 'TF' alone is ok since subset of 'Species'.  Overwrites default
# list in mutation.py
list_types_output=['TF']

# used in evolution_gillespie,
# =0 serial processing, only one C job running at a time
# =1 threaded, multiple C jobs started on one machine to use multi-core capabilities
# =2 multiple computers, cluster, using mpirun See HowTo in /Doc
prmt['multipro_level'] = 1

# turns on the pareto module in run_evolution.py
prmt['pareto']=0 # set to 1 to run pareto optimization
prmt['npareto_functions']= 1 # number of pareto functions
prmt['rshare']= 0 #radius to use for fitness sharing
				  # set to zero to turn off

# Control the network and species(time) plotting from evolution_gillespie.py (off for cluster)
prmt['plot']=0

# To restart from a Restart_* file created in Seed* directories:
#    Move the Restart file into a new model directory (see run_evolution -m option)
#    Adjust the entries below as desired and run_evolution module.
#    Note option to either exactly reproduce prior evolution simulation or just initialize with some
# prior population and then do indepedent evolution

# The parameters for restart related functions, grouped into subdictionary, used in evolution_gillespie.py
# The flag to restart from file is file name, None skips restart file and uses new instances of init functions
# If restart the numbering of generations continuous with that of restart file.
# To redo a saved population with different random number, use ['same_seed'] = False
# Note the 'freq' parameter below, if ~= npopulation then Restart file size ~1/2 Bests file
prmt['restart'] = {}
prmt['restart']['file'] =  None     # None skips restart, otherwise name of file in model directory (-m option run_evolution.py)
prmt['restart']['kgeneration'] = 2  # restart from After this generation number (see loop in Population.evolution)
prmt['restart']['freq'] = 500       # save population every freq generations
prmt['restart']['same_seed'] = True # get seed of random() from restart file to reproduce prior data.

# necessary imports to define following functions
import random
from . import mutation

# generic init network function that just creates desired number of input/output species
def init_network():
   """
   Generate a initial network to start the algorithm with.

   Return:
        initial :class:`Mutable_Network <phievo.Networks.mutation.Mutable_Network>`
   """
   seed=int(random.random()*100000)
   g=random.Random(seed)
   L=mutation.Mutable_Network(g)

   # see the add_type function in classes_eds2.Species for options
   for k in range(0, prmt['ninput']):
      spec = L.random_Species()
      spec.add_yype(['Input', k])

   for k in range(0, prmt['noutput']):
      [tm, prom, spec] = L.random_gene()
      spec.add_type(['Output', k])

   # checks consecutive numbering for IO species.
   L.verify_IO_numbers()

   # examples of network wide flags to control properties of genes and interactions. See Doc/ or Network class defns
   L.activator_required=1
   L.fixed_activity_for_TF=0
   L.write_id()

   return L


def fitness_treatment(population):
    """
    Function to change the fitness of the networks.
    It can add a small random number to the computed fitness for example.

    Args:
        population (:class:`Population <phievo.Population_Types.evol_gillespie.Population>`): Population to modify.
    """


    for nnetwork in range(population.npopulation):
            population.genus[nnetwork].fitness=(float(population.genus[nnetwork].data_evolution[1])+0.01)/(float(population.genus[nnetwork].data_evolution[2])+0.001)
