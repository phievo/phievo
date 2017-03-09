""" initialization script, import and copy data to like named variables in other modules.
"""

# ranges over which parmeters in classes_eds2 classes can vary. Defined in mutation.py
# NB use 1.0 etc to keep params real not int.  

T=1.0 #typical time scale
C=1.0 #typical concentration
L=1.0 #typical size for diffusion

# dictionary key is the Class.attribute whose range is given as a real number or as a LIST
# indicating the interval over which it has to vary

dictionary_ranges={}
dictionary_ranges['Species.degradation']=[0.1,1.0/T]
dictionary_ranges['Species.diffusion']=L   # for ligands diffusion 
dictionary_ranges['TModule.rate']=C/T
dictionary_ranges['TModule.basal']=0.0*C/T #basal rate
dictionary_ranges['CorePromoter.delay']=0.0   # convert to int(T/dt) in run_evolution.py
dictionary_ranges['TFHill.hill']= [0.5, 5.0]
dictionary_ranges['TFHill.threshold']=[0.0,C]
dictionary_ranges['PPI.association']=1.0/(C*T)
dictionary_ranges['PPI.disassociation']=1.0/T
dictionary_ranges['Phosphorylation.rate']=1.0/T
dictionary_ranges['Phosphorylation.hill']=5.0
dictionary_ranges['Phosphorylation.threshold']=[0.01,C]
dictionary_ranges['Phosphorylation.dephosphorylation']=1.0/T
dictionary_ranges['LR.association']=1.0/T
dictionary_ranges['LR.concentration']=C

# when key='relative_variation' is present, use this fraction of parameter range listed
# above, to vary parameter for each mutation.  Comment out line to pick parameter
# randomly from range defined above.
#dictionary_ranges['relative_variation']=0.1

#################################################################################
# names of c-code files needed by deriv2.py (. for extension only)
# normally need 'header', 'utilities' (loaded before python derived routines and
# ['fitness', 'geometry', 'init_history', 'input', 'integrator', 'main' ] placed afterwards
# skip by setting cfile[] = ' ' or ''
 
cfile = {}
cfile['header'] = 'integrator_header.h'
cfile['utilities'] = 'utilities.c'
cfile['fitness'] = 'Examples/SlidingHox/fitness_hox3_information.c'
cfile['geometry'] = 'linear_geometry.c'
cfile['init_history'] = 'Examples/SlidingHox/init_history_0.c'
cfile['input'] =  'Examples/SlidingHox/input_Hox_euler.c'
cfile['integrator'] = 'euler_integrator.c'
cfile['main'] = 'main_general.c'

#################################################################################
# mutation rates
dictionary_mutation={}

# Rates for speices nodes to add
dictionary_mutation['random_gene()']=0.0
dictionary_mutation['random_gene(\'TF\')']=0.03
dictionary_mutation['random_gene(\'Kinase\')']=0.0
dictionary_mutation['random_gene(\'Ligand\')']=0.0
dictionary_mutation['random_gene(\'Receptor\')']=0.0

# Rates for interaction nodes to add
dictionary_mutation['random_Interaction(\'TFHill\')']=0.05
dictionary_mutation['random_Interaction(\'PPI\')']=0.0
dictionary_mutation['random_Interaction(\'Phosphorylation\')']=0.0
dictionary_mutation['random_Interaction(\'LR\')']=0.0

# Rates for nodes to remove. NB 'CorePromoter' rate kills all Species nodes
dictionary_mutation['remove_Interaction(\'TFHill\')']=0.4
dictionary_mutation['remove_Interaction(\'PPI\')']=0.0
dictionary_mutation['remove_Interaction(\'CorePromoter\')']=0.1
dictionary_mutation['remove_Interaction(\'Phosphorylation\')']=0.0
dictionary_mutation['remove_Interaction(\'LR\')']=0.0

# Rates to change parameters for a node
dictionary_mutation['mutate_Node(\'Species\')']=0.5
dictionary_mutation['mutate_Node(\'TFHill\')']=0.5
dictionary_mutation['mutate_Node(\'CorePromoter\')']=0.5
dictionary_mutation['mutate_Node(\'TModule\')']=0.5
dictionary_mutation['mutate_Node(\'PPI\')']=0.0
dictionary_mutation['mutate_Node(\'LR\')']=0.0
dictionary_mutation['mutate_Node(\'Phosphorylation\')']=0.0

#rates to change output tags.  See list_types_output array below  
dictionary_mutation['random_add_output()']=0.6
dictionary_mutation['random_remove_output()']=0.1
dictionary_mutation['random_change_output()']=0.01



dictionary_mutation['random_duplicate()']=0.01
#############################################################################
# parameters in various modules, created as one dict, so that can then be passed as argument

# Needed in deriv2.py to create C program from network
prmt = {}
prmt['nstep'] =3000         #number of time steps
prmt['ncelltot']=20            #number of cells in an organism
prmt['nneighbor'] = 3 # must be >0, whatever geometry requires, even for ncelltot=1
prmt['ntries'] = 1    # number of initial conditions tried in C programs
prmt['dt'] = 0.05     # time step

# Generic parameters, transmitted to C as list or dictionary.
# dict encoded in C as: #define KEY value (KEY converted to CAPS)
# list encoded in C as: static double free_prmt[] = {}
# The dict is more readable, but no way to test in C if given key supplied.  So
# always include in C: #define NFREE_PRMT int. if(NFREE_PRMT) then use free_prmt[n]
# prmt['free_prmt'] = { 'step_egf_off':20 }  # beware of int vs double type
# prmt['free_prmt'] = [1,2] 

# Needed in evolution_gill to define evol algorithm and create initial network
prmt['npopulation'] =10
prmt['ngeneration'] =25001
prmt['tgeneration']=0.01       #initial generation time (for gillespie), reset during evolution
prmt['noutput']=5    # to define initial network
prmt['ninput']=1
prmt['freq_plot'] = 100  #plot time course every generation%freq_plot = 0
prmt['freq_stat'] = 5     # print stats every freq_stat generations
prmt['frac_mutate'] = 0.5 #fraction of networks to mutate
prmt['redo'] = 1   # rerun the networks that do not change to compute fitness for different IC

# used in run_evolution, 
prmt['nseed'] = 20   # number of times entire evol procedure repeated, see main program.  
prmt['firstseed'] = 0  #first seed

# used in evolution_gillespie.  =1 if using parallel processing(threading or pypar) =0 for serial processing
prmt['multipro_level']=0
prmt['pareto']=0
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
prmt['restart']['file'] =  None  # None skips restart, otherwise name of file in model directory (-m option run_evolution.py)
prmt['restart']['kgeneration'] = 2  # restart from After this generation number (see loop in Population.evolution) 
prmt['restart']['freq'] = 50  # save population every freq generations
prmt['restart']['same_seed'] = True  # get seed of random() from restart file to reproduce prior data.

list_unremovable=['Input']

list_types_output=['TF'] #only TF can be outputs, if remove('Species')


# necessary imports to define following functions
import random                                              
import mutation

# Two optional functions to add input species or output genes, with IO index starting from 0.
# Will overwrite default variants in evol_gillespie if supplied below



def init_network():
   seed=int(random.random()*1000000)
   
   from base import L
   import copy
   net=copy.deepcopy(L)
   net.Random=random.Random(seed)
   return net

def fitness_treatment(population): 
    # Function to change the fitness of the networks
   
   
    for nnetwork in range(population.npopulation):
       population.genus[nnetwork].fitness-=0.001*random.random()
          
       
