""" initialization script, import and copy data to like named variables in other modules.
"""
# ranges over which parmeters in classes_eds2 classes can vary. Defined in mutation.py
# NB use 1.0 etc to keep params real not int.
T=1.0 #typical time scale
C=1.0 #typical concentration
L=1.0 #typical size for diffusion

###########################################
## SET THE PARAMETER RANGES OF VARIATION ## 
###########################################

## Set the allowed range of variations for the interaction parameters
## If one number is give, it sets the upper bound, if it's an array of
## two it sets the lower and upper bounds.
dictionary_ranges={}
dictionary_ranges['Species.degradation'] = 1.0/T
##dictionary_ranges['CorePromoter.delay']= 0   
dictionary_ranges['TModule.rate'] = C/T
dictionary_ranges['TModule.basal'] = 0.0
dictionary_ranges['TFHill.hill'] = [1.0,5.0]
dictionary_ranges['TFHill.threshold'] = C

##################################################
## LINK TO THE C FILES USED FOR THE INTEGRATION ##
##################################################
cfile = {}
cfile['fitness'] = 'fitness.c'
cfile['init_history'] = 'init_history.c'
cfile['input'] =  'input.c'

###############################
## SET THE RATES OF MUTATION ##
###############################

## Rate of the possible mtation rates
## Note that the rates a normalized to have an everage of 1 mutation per generation per network
## in the Gillespie
##    Adding genes
dictionary_mutation={}
dictionary_mutation['random_gene()']=0.01
dictionary_mutation['random_gene(\'TF\')']=0.01

##    Adding interactions:
dictionary_mutation['random_Interaction(\'TFHill\')']=0.002

#     Removing nodes
dictionary_mutation['remove_Interaction(\'TFHill\')']=0.005

#     Mutate the parameters
dictionary_mutation['mutate_Node(\'Species\')']=0.1
dictionary_mutation['mutate_Node(\'TFHill\')']=0.1
dictionary_mutation['mutate_Node(\'TModule\')']=0.1

#     Change the ouput genes
dictionary_mutation['random_add_output()']=0.0
dictionary_mutation['random_remove_output()']=0.0
dictionary_mutation['random_change_output()']=0.1

##########################
## SET MODEL PARAMETERS ##
##########################

## Integration parameters
prmt = {}
prmt['nstep'] =30000 # number of time steps
prmt['dt'] = 0.05 # time step
prmt['ncelltot']=1 # number of cells in an organism
prmt['ntries'] = 10 # number of initial conditions tried in C programsprmt['npopulation'] =50

## Other settings
prmt['ngeneration'] =101
prmt['npopulation'] = 50
prmt['nneighbor'] = 3 # must be >0, whatever geometry requires, even for ncelltot=1
prmt['tgeneration']=1.0       #initial generation time (for gillespie), reset during evolution
prmt['noutput']=1    # to define initial network
prmt['ninput']=1
prmt['freq_stat'] = 5     # print stats every freq_stat generations
prmt['frac_mutate'] = 0.5 #fraction of networks to mutate
prmt['redo'] = 1   # rerun the networks that do not change to compute fitness for different IC

## Software parameters:
prmt['nseed'] = 5  # number of times entire evol procedure repeated, see main program.
prmt['firstseed'] = 0  #first seed
prmt['multipro_level']=1
prmt['pareto']=0
prmt['npareto_functions'] = 2
prmt['rshare']= 0.0
prmt['plot']=0
prmt['langevin_noise']=0


prmt['restart'] = {}
prmt['restart']['activated'] = False ## If false False, start the integration from scratch.
prmt['restart']['freq'] = 50  # save the complete population of networks every freq generations
## Customize restart point:
#prmt['restart']['seed'] =  0 # the directory of the population you want to restart from
#prmt['restart']['kgeneration'] = 50  # restart from After this generation number
#prmt['restart']['same_seed'] = True # Backup the random generator to the same seed

list_unremovable=['Input','Output']
list_types_output=['TF']

###############################
## DEFINE AN INITIAL NETWORK ##
###############################
import random
from phievo.Networks import mutation

def init_network():
   seed=int(random.random()*100000)
   g=random.Random(seed)
   net=mutation.Mutable_Network(g)
   ## Input
   parameters=[]
   parameters.append(['TF',1])
   parameters.append(['Input',0])
   TF=net.new_Species(parameters)

   ## Add one output
   N_output = 1
   for k in range(2):
       [tm, prom, o1] = net.random_gene('TF')
       o1.add_type(['Output',k])
   # for k in range(N_output):
   #     [tm, prom, o1] = net.random_gene(['TF',1])
   #     o1.add_type(['Output',k])
   net.write_id()
   return net

def fitness_treatment(population):
    """Function to change the fitness of the networks"""
    ## Does nothing by default
    pass
    
