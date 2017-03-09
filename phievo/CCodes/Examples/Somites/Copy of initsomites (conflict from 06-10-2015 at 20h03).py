""" initialization script, import and copy data to like named variables in other modules.
"""

# ranges over which parmeters in classes_eds2 classes can vary. Defined in mutation.py
# NB use 1.0 etc to keep params real not int.  

T=1.0 #typical time scale
C=1.0 #typical concentration
L=1.0 #typical size for diffusion

#dictionary key is the Class.attribute whose range is give an a real number or as a LIST
# indicating the interval over which it has to vary



dictionary_ranges={}
dictionary_ranges['Species.degradation']=[0.1,1.0/T]
dictionary_ranges['Species.diffusion']=L   # for ligands diffusion 
dictionary_ranges['TModule.rate']=C/T
dictionary_ranges['TModule.basal']=0.0*C/T #basal rate
dictionary_ranges['CorePromoter.delay']=10.0   # convert to int(T/dt) in run_evolution.py
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




#################################################################################
# names of c-code files needed by deriv2.py (. for extension only)
# normally need 'header', 'utilities' (loaded before python derived routines and
# ['fitness', 'geometry', 'init_history', 'input', 'integrator', 'main' ] placed afterwards
# skip by setting cfile[] = ' ' or ''

cfile = {}
cfile['header'] = 'integrator_header.h'
cfile['utilities'] = 'utilities.c'
cfile['fitness'] = 'Examples/Somites/fitness_somites.c'
cfile['geometry'] = 'linear_geometry.c'
cfile['init_history'] = 'Examples/Somites/init_history_0.c'
cfile['input'] =  'Examples/Somites/input_somites.c'
cfile['integrator'] = 'euler_integrator.c'
cfile['main'] = 'main_general.c'

#################################################################################
# mutation rates
dictionary_mutation={}

# Rates for nodes to add
dictionary_mutation['random_gene()']=0.0
dictionary_mutation['random_gene(\'TF\')']=0.1
dictionary_mutation['random_gene(\'Kinase\')']=0.0
dictionary_mutation['random_gene(\'Ligand\')']=0.0
dictionary_mutation['random_gene(\'Receptor\')']=0.0


dictionary_mutation['random_Interaction(\'TFHill\')']=0.1
dictionary_mutation['random_Interaction(\'PPI\')']=0.0
dictionary_mutation['random_Interaction(\'Phosphorylation\')']=0.0
dictionary_mutation['random_Interaction(\'LR\')']=0.0

# Rates for nodes to remove
dictionary_mutation['remove_Interaction(\'TFHill\')']=0.1
dictionary_mutation['remove_Interaction(\'PPI\')']=0.0
dictionary_mutation['remove_Interaction(\'CorePromoter\')']=0.1
dictionary_mutation['remove_Interaction(\'Phosphorylation\')']=0.0
dictionary_mutation['remove_Interaction(\'LR\')']=0.0

# Rates to change parameters for a node
dictionary_mutation['mutate_Node(\'Species\')']=0.1
dictionary_mutation['mutate_Node(\'TFHill\')']=0.1
dictionary_mutation['mutate_Node(\'CorePromoter\')']=0.1
dictionary_mutation['mutate_Node(\'TModule\')']=0.1
dictionary_mutation['mutate_Node(\'PPI\')']=0.0
dictionary_mutation['mutate_Node(\'LR\')']=0.0
dictionary_mutation['mutate_Node(\'Phosphorylation\')']=0.0
#rates to change outputs
dictionary_mutation['random_add_output()']=0.0
dictionary_mutation['random_remove_output()']=0.0
dictionary_mutation['random_change_output()']=0.1


#############################################################################
# parameters in various modules, created as one dict, so that can then be passed as argument

# Needed in deriv2.py to create C program from network
prmt = {}
prmt['nstep'] =12000         #number of time steps
prmt['ncelltot']=40            #number of cells in an organism
prmt['nneighbor'] = 3 # must be >0, whatever geometry requires, even for ncelltot=1
prmt['ntries'] = 1    # number of initial conditions tried in C programs
prmt['dt'] = 0.01      # time step

# Needed in evolution_gill to define evol algorithm and create initial network
prmt['npopulation'] =50
prmt['ngeneration'] =501
prmt['tgeneration']=1.0       #initial generation time (for gillespie)
prmt['noutput']=3     # to define initial network
prmt['ninput']=1
prmt['freq_plot'] = 40    #plot time course every generation%freq_plot = 0
prmt['freq_stat'] = 5     # print stats every freq_stat generations
prmt['frac_mutate'] = 0.5 #fraction of networks to mutate
prmt['redo'] = 0
# number of times entire evol procedure repeated, see main program.
prmt['nseed'] = 20     



prmt['multipro_level']=2

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



# reset the range for delay to integer units of dt, and optionally freeze the output variables

list_unremovable=['Input','Output']
list_types_output=['TF']



#dictionary_ranges['CorePromoter.delay'] = int( 0.5 + dictionary_ranges['CorePromoter.delay']/prmt['dt'] )

# list_IO_words = ['Input', 'Output'] this does not work, must modify existing object

# necessary imports to define following functions
import random                                              
import mutation

# Two optional functions to add input species or output genes, with IO index starting from 0.
# Will overwrite default variants in evol_gillespie if supplied below




def init_network():
   seed=int(random.random()*100000)
   g=random.Random(seed)
   L=mutation.Mutable_Network(g)
   parameters=[['Degradable', mutation.sample_dictionary_ranges('Species.degradation',random) ]]
   parameters.append(['TF',1])
   parameters.append(['Input',0])
   TF=L.new_Species(parameters)
   for k in xrange(2):
       [tm, prom, o1] = L.random_gene('TF')
       L.add_output2species(k, o1)
   L.activator_required=1
   L.fixed_activity_for_TF=0
   L.write_id()
   return L

def fitness_treatment(population): 
    # Function to slightly change the fitness of the networks
   
   
    for nnetwork in xrange(population.npopulation):
       population.genus[nnetwork].fitness-=0.001*random.random()
          
       

