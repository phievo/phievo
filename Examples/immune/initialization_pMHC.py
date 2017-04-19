""" initialization script, import and copy data to like named variables in other modules.
"""
import numpy as np
# ranges over which parmeters in classes_eds2 classes can vary. Defined in mutation.py
# NB use 1.0 etc to keep params real not int.
T=1.0 #typical time scale
C=1.0 #typical concentration
L=1.0 #typical size for diffusion
R=10000 # typical number of receptors
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
dictionary_ranges['Simple_Phosphorylation.rate']=1.0/(C*T)
dictionary_ranges['Simple_Phosphorylation.spontaneous_dephospho']=1.0/T
dictionary_ranges['Simple_Dephosphorylation.rate']=1.0/(C*T)
dictionary_ranges['Initial_Concentration.concentration'] = C
dictionary_ranges['Mutable_Threshold.thresh'] = 1.0
dictionary_ranges['KPR_Binding.association'] = 2.0/(R*T)
dictionary_ranges['LR.association']=1.0/T
dictionary_ranges['LR.threshold']=C

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
cfile['header'] = 'integrator_header_improved.h'
cfile['utilities'] = 'utilities_pMHC.c'
cfile['fitness'] = 'fitness_pMHC_mutinfo2_ago.c'
cfile['geometry'] = 'linear_geometry.c'
cfile['init_history'] = 'init_history_pMHC.c'
cfile['integrator'] = 'integrator_pMHC_improved.c'
cfile['main'] = 'main_pMHC.c'
cfile['input'] =  'input.c'
# this one was commented : 'immune/input_adaptation1.c'

pfile = {}
pfile["deriv2"] = "immune.Immune.deriv2_pMHC_modifier"
pfile["interaction"] = "immune.Immune.interaction_pMHC"
pfile["pretty_graph"] = "immune.Immune.pretty_graph_2_pMHC"
pfile["plotdata"] = "immune.Immune.plotdata_pMHC"

#################################################################################

# mutation rates
dictionary_mutation={}

# Rates for nodes to add
dictionary_mutation['random_gene()']=0.00
dictionary_mutation['random_gene(\'TF\')']=0.00
dictionary_mutation['random_gene(\'Kinase\')']=0.00
dictionary_mutation['random_gene(\'Ligand\')']=0.00
dictionary_mutation['random_gene(\'Receptor\')']=0.00

dictionary_mutation['random_Interaction(\'TFHill\')']=0.000
dictionary_mutation['random_Interaction(\'PPI\')']=0.000
dictionary_mutation['random_Interaction(\'Phosphorylation\')']=0.000
dictionary_mutation['random_Interaction(\'Simple_Phosphorylation\')']=0.05
dictionary_mutation['random_Interaction(\'Simple_Dephosphorylation\')']=0.05
dictionary_mutation['random_Interaction(\'LR\')']=0.0

# Rates for nodes to remove
dictionary_mutation['remove_Interaction(\'TFHill\')']=0.000
dictionary_mutation['remove_Interaction(\'PPI\')']=0.000
dictionary_mutation['remove_Interaction(\'CorePromoter\')']=0.00
dictionary_mutation['remove_Interaction(\'Phosphorylation\')']=0.00
dictionary_mutation['remove_Interaction(\'LR\')']=0.0
dictionary_mutation['remove_Interaction(\'Simple_Phosphorylation\')']=0.5
dictionary_mutation['remove_Interaction(\'Simple_Dephosphorylation\')']=0.5
dictionary_mutation['remove_Interaction(\'Initial_Concentration\')']=0.00
dictionary_mutation['remove_Interaction(\'KPR_Binding\')']=0.00
dictionary_mutation['remove_Interaction(\'KPR_Unbinding\')']=0.00

# Rates to change parameters for a node
dictionary_mutation['mutate_Node(\'Species\')']=0.5
dictionary_mutation['mutate_Node(\'TFHill\')']=0.0
dictionary_mutation['mutate_Node(\'CorePromoter\')']=0.0
dictionary_mutation['mutate_Node(\'TModule\')']=0.0
dictionary_mutation['mutate_Node(\'PPI\')']=0.0
dictionary_mutation['mutate_Node(\'LR\')']=0.0
dictionary_mutation['mutate_Node(\'Phosphorylation\')']=0.0
dictionary_mutation['mutate_Node(\'Simple_Phosphorylation\')']=0.5
dictionary_mutation['mutate_Node(\'Simple_Dephosphorylation\')']=0.5
dictionary_mutation['mutate_Node(\'Initial_Concentration\')']=0.5
dictionary_mutation['mutate_Node(\'Mutable_Threshold\')']=0.0
dictionary_mutation['mutate_Node(\'KPR_Binding\')']=0.5

#rates to change output tags.  See list_types_output array below
dictionary_mutation['random_add_output()']=0.0
dictionary_mutation['random_remove_output()']=0.0
dictionary_mutation['random_change_output()']=0.5

# Rates for random creation of a species
dictionary_mutation['random_Species(\'Kinase\')'] = 0.05
dictionary_mutation['random_Species(\'Phosphatase\')'] = 0.05

# duplicates gene and some interactions, see mutation.random_duplicate()
# dictionary_mutation['random_duplicate()']=0.01
# Rate for shifting kinase on Simple_Phosphorylations.
dictionary_mutation['random_shift_Simple_Phosphorylation()']=0.0

#############################################################################

# parameters in various modules, created as one dict, so that can then be passed as argument

# Needed in deriv2.py to create C program from network
prmt = {}
prmt['nstep'] =2        #number of time steps
prmt['nneighbor'] = 3  # must be >0, whatever geometry requires, even for ncelltot=1
prmt['ntries'] = 1    # number of initial conditions tried in C programs
prmt['dt'] = 0.0001      # Initial time step
prmt['noutput'] = 1    # to define initial network, recomputed by python from network
prmt['ninput'] = 1     # to define init network, does not change during evolution, computed in python from network

# if >0, turns on langevin integration and value becomes CONCENTRATION_SCALE in CCode. should be ~100 ie # molecules
prmt['langevin_noise'] = 0

##########################################
# Parameters added for the integration.
prmt['max_time'] = 150       # Maximum time for integration
prmt['tolerance'] = 1E-9     # Desired tolerance for integration.
prmt['tiny'] = 1E-15
prmt['scalebelow'] = 0.3     # Smallest factor by which the time step can decrease in integration.
prmt['scaleabove'] = 3       # Largest factor by which the time step can increase in integration.

##########################################
# Parameters added which are specific to the model.
prmt['tau_off'] = [3,10]              # Dissociation times considered.
prmt['Ligands'] = [int(n) for n in np.logspace(0,3,20)] # Initial ligand concentrations considered.
prmt['ntau'] = len(prmt['tau_off'])	      # Number of different dissociation times.
prmt['nLigands'] = len(prmt['Ligands'])   # Number of different ligand concentrations.
prmt['ncelltot']=prmt['ntau']*prmt['nLigands']        # Number of cells in an organism (in the present case, each cell represents a different ligand concentration and dissociation time).
prmt['Receptor_Conc'] = 30000             # Initial receptor concentration (kappa can be found in the initialized network below).

#########################################
import math
prmt['n_partition_output'] = math.ceil(len(prmt['Ligands'])/2) # Related to mutual information fitness (number of bins).

#########################################
# Related to self-resistance.
prmt['Ligand_self'] = [0]			# number of self ligands.
prmt['tau_self'] = 0.3				# binding time of antagonists.

############################################

# Generic parameters, transmitted to C as list or dictionary.
# dict encoded in C as: #define KEY value (KEY converted to CAPS)
# list encoded in C as: static double free_prmt[] = {}
# The dict is more readable, but no way to test in C if given key supplied.  So
# always include in C: #define NFREE_PRMT int. if(NFREE_PRMT) then use free_prmt[n]
# prmt['free_prmt'] = { 'step_egf_off':20 }  # beware of int vs double type
# prmt['free_prmt'] = [1,2]

# Needed in evolution_gill to define evol algorithm and control output
prmt['npopulation'] =50  # number of species in population
prmt['ngeneration'] =500   # number of generations of evolution (+1 to get stats last step)
prmt['tgeneration']=0.01    # initial generation time (for gillespie), reset during evolution
prmt['freq_plot'] = 10000     # plot time course every generation%freq_plot = 0
prmt['freq_stat'] = 5       # print stats every freq_stat generations
prmt['frac_mutate'] = 0.5   # fraction of networks to mutate, keep at 1/2
prmt['redo'] = 0            # rerun the networks that do not change to compute fitness for different IC

# used in run_evolution,
prmt['nseed'] = 5		# number of times entire evol procedure repeated, see main program.
prmt['firstseed'] = 0	# first seed, defines random number seed, so can exactly repeat
prmt['agefitness'] = 0

# Any species with these tags can not be removed.  MUST include 'Input'.  If fixed number of

# outputs add 'Output' to list. Overwrites default list in classes_eds2.py
list_unremovable=['Input','Receptor','Ligand','Mutable_Threshold']

# Controls the types of outputs, when the output tags can move or be created.  MUST
# be at least 'Species'.  NB 'TF' alone is ok since subset of 'Species'.  Overwrites default
# list in mutation.py
list_types_output=['Kinase','Phosphatase']

# used in evolution_gillespie,
# =0 serial processing, only one C job running at a time
# =1 threaded, multiple C jobs started on one machine to use multi-core capabilities
# =2 multiple computers, cluster, using mpirun See HowTo in /Doc
prmt['multipro_level'] = 1

# Parameters for pareto evolution module. Pareto optimization allows one to use multiple
# fitness functions. When running pareto optimization, network X is only considered better
# than network Y if the all fitness functions evaluated for X are less than or equal to those for
# Y (with at least one that is less than, not equal). phievo.Networks of pareto rank 1 are those for which there is
# no other network which is considered better
#To increase the diversity of networks, fitness sharing is implemented. If two networks are closer than ['rshare']
# in the space of fitness functions, their pareto ranks are increased (but are still consider better than those of pareto rank 2)
#this gives a selective advantage to rarer networks in population which is pareto rank 1

prmt['pareto']=0 # a flag, set to 1 to run pareto
prmt['npareto_functions']= 2 # number of functions to use for pareto optimization.
prmt['rshare']= 15 #radius to use for fitness sharing, set to zero to turn off

# Control the network and species(time) plotting from evolution_gillespie.py (off for cluster)
prmt['plot']=0

## prmt['restart']['kgeneration'] tells phievo to save a complete generation
## at a given frequency in a restart_file. The algorithm can relaunched at
## backuped generation by turning  prmt['restart']['activated'] = True and
## setting the generation and the seed. If the two latters are not set, the
## algorithm restart from the highest generation  in the highest seed.
prmt['restart'] = {}
prmt['restart']['activated'] = False #indicate if you want to restart or not
prmt['restart']['freq'] = 50  # save population every freq generations
#prmt['restart']['seed'] =  0 # the directory of the population you want to restart from
#prmt['restart']['kgeneration'] = 50  # restart from After this generation number
#prmt['restart']['same_seed'] = True # Backup the random generator to the same seed

# necessary imports to define following functions
import random
from phievo.Networks import mutation

# generic init network function that just creates desired number of input/output species

def init_network():

    seed=int(random.random()*510)
    g=random.Random(seed)
    L=mutation.Mutable_Network(g)

    conc = 1000.0
    T = 1.0

    # the ligand (the pMHC on the antigen presenting cell)
    parameters=[['Ligand']]
    parameters.append(['Input',0])
    Lig=L.new_Species(parameters)

    # the receptor (on the T cell)
    parameters=[['Receptor']]
    R = L.new_Species(parameters)

    # a kinase.
    parameters=[['Kinase']]
    parameters.append(['Phosphorylable'])
    parameters.append(['Phospho',0])
    K1 = L.new_Species(parameters)

    # a kinase.
    parameters=[['Kinase']]
    parameters.append(['Phosphorylable'])
    parameters.append(['Phospho',0])
    K2 = L.new_Species(parameters)

    # a kinase.
    parameters=[['Kinase']]
    parameters.append(['Phosphorylable'])
    parameters.append(['Phospho',0])
    K3 = L.new_Species(parameters)

    # a phosphatase.
    parameters=[['Phosphatase']]
    parameters.append(['Phosphorylable'])
    parameters.append(['Phospho',0])
    P1 = L.new_Species(parameters)

    # a phosphatase.
    parameters=[['Phosphatase']]
    parameters.append(['Phosphorylable'])
    parameters.append(['Phospho',0])
    P2 = L.new_Species(parameters)

    # a phosphatase.
    parameters=[['Phosphatase']]
    parameters.append(['Phosphorylable'])
    parameters.append(['Phospho',0])
    P3 = L.new_Species(parameters)

    # the complex generated: it is a kinase that is phosphorylable and with the label pMHC to indicate it is in the cascade.
    parameters=[['Kinase']]
    parameters.append(['Phosphorylable'])
    parameters.append(['Phospho',0])
    parameters.append(['pMHC'])

    # the first binding of the receptor and ligand generating the complex.
    kappa = 1E-4
    [Binding, C] = L.new_KPR_Binding(Lig,R,kappa,parameters)
    unbinbing = L.new_KPR_Unbinding(Lig,R,C)

    # the output tag placed on the unphosphorylated complex arising from binding of receptor and ligand.
    C.add_type(['Output',0])
    L.write_id()

    # checks consecutive numbering for IO species.
    L.verify_IO_numbers()
    return L


def fitness_treatment(population):
    """Function to change the fitness of the networks"""
    import random
    for ind in population:
        if ind.data_evolution:
            ind.fitness = float(ind.data_evolution[0])+0.01*random.random()
        else:
            ind.fitness = None
