"""
Define the Class Population with her principal method, evolution, which
evolve a set of networks. All initialization done from an initialization.py
file. All the modules are initialized through run_evolution.py.

The initial networks to evolve, can be built from just the input/output genes,
a predefined newtork, or restarted from any saved population from a previous
run. (See initialization file for details)

The time between generations is variable, and about the same for all species,
we sample the mutation rates with a gillespie like algorithm, hence the name

The evolution method will write the following files in the namefolder given
as argument to Population.__init__ stdout basic info each generation:
* Bests = for generation, the network with best fitness in text form to edit or process with stat_best_net.py
* Restart* = binary dbm type file with data to restart evolution at selected generation numbers
* graphic files with time course and best network diagram at selected generations
"""
print('Execute evolution_gillespie.py')

import phievo.Networks.classes_eds2 as classes_eds2
import phievo.Networks.plotdata as plotdata
import phievo.Networks.mutation as mutation
import gc # Garbage collector
from math import log,sqrt
import copy,os,random,sys,glob
import shelve
import time, pickle, dbm  # for restart's
import phievo.Populations_Types.population_stat as pop_stat

#########################
### Global Parameters ###
#########################

# read various parameters for ODE generation and evolution from input file into dictionary called prmt
prmt = {}

# dictionary of line labels and commands to define statistics extracted from networks see Class NetworkStat
stat_dict = {}
stat_dict['Node'] = lambda net: net.list_types['Node']
stat_dict['Species'] = lambda net: net.list_types['Species']

count_interactions = ['TFHill', 'PPI', 'LR', 'Phosphorylation']
count_keys = list(mutation.dictionary_mutation.keys())
for i in count_interactions:
    rates0 = [mutation.dictionary_mutation[kk]==0 for kk in count_keys if kk.find(i) >= 0]
    if len(rates0) > 0:
        continue
    stat_dict[i] = 'list_types[{0}]'.format(i)

#######################
### Dummy Functions ###
#######################

def fitness_treatment(population):
    """default function for fitness treatment
    
    If necessary, should be implemented in the init*.py file
    """
    pass

def init_network(mutation):
    """Default function to create network
    
    It must be overwritten with function from the init*.py file
    otherwise stop the programm
    """
    raise NotImplementedError("must supply init_network() in initialization.py")

#########################
### General Functions ###
#########################
def restart(directory, generation, verbose = True):
    """Allow the user to restart an old run

        Args:
            directory (str): the directory of the restart file
            generation (int): the generation number
        
        Returns:
            rprmt (dict): the parameters of the run
            genus (list): the list of individuals of the population
        """
    restart_file = os.path.join(directory,"Restart_file")
    with shelve.open(restart_file) as restart_data:
        rprmt, genus = restart_data[str(generation)]
        if verbose:
            print('successfully restarted from file= ', dir, 'generation= ', generation)
            print('header=', rprmt['header'])
        return rprmt, genus

###################################
### Class Population Definition ###
###################################

class Population(object):
    """Define a population as a list of networks called Population. genus and a principal method evolution; object means it is a newstyle class ! See the web for distinction between new and olds style class, important for inheritance"""

    """Define a population as a list of networks called Population.genus
    
    Attributes:
        best_fitness (float): keep trace of the best fitness in the population
        genus (list): the list of individuals of the population
        same_seed (bool): indicate if the file is a restart or not
        tgeneration (float): starting hop time for the gillespie algorithm
        npopulation (int): size of te population
        bests_file (str): directory to save the data of evolution
    
    Main methods:
        evolution: launch the evolutionary algorithm
        pop_mutate_and_integrate: update the whole population
    """
    def __init__(self, namefolder):
        """ Copy a few parameters from prmt-dictionary that logically belong to the population, then
        setup various files and from initialization file decide how to initialize networks in genus.
        """
        self.best_fitness = sys.maxsize  # tag to determine if fitness increasing each generation
        self.best_fitness_counter = 0   # number of generations that best fitness has not changed.
        self.same_seed = False  # False -> starting from new data, reset to True below if replicating restart file
        self.generation0 = 0  # number of first generation, if replicating restart file set>0
        self.tgeneration = prmt['tgeneration']
        self.npopulation = prmt['npopulation']
        self.namefolder = namefolder   # directory where all data going
        self.n_mutations = 0 #number of mutations per generation

        #file to hold best network each generation
        self.data_file = namefolder + os.sep + 'data'
        self.bests_file = namefolder + os.sep + 'Bests_{}.net'

        # unique file name to save restart data,
        self.restart_file = namefolder + os.sep + 'Restart_file'

        # reset crucial parameters for loop over generations in evolution.
        # Set self.same_seed = True if want to exactly
        # recreate evolution that lead to restart data
        if prmt['restart']['activated']:
            rprmt, self.genus = restart(prmt['restart']['dir'], prmt['restart']['kgeneration'] )
            self.tgeneration = rprmt['tgeneration']
            if prmt['restart']['same_seed']:
                random.setstate( rprmt['state'] )
                self.same_seed = True
                self.generation0 = 1 + prmt['restart']['kgeneration']
            return None

        # no restart, generate randomized list of networks from init file or routines supplied here.
        self.genus=[]
        for i in range(self.npopulation):
            L = init_network()
            L.write_id()
            self.genus.append(L)

    def __getitem__(self,index):
        """Allow the population to be indexed as a list"""
        return self.genus[index]

    def __len__(self):
        """Overload the population size"""
        return self.npopulation

    def save_bests_file(self, data_str):
        """Write the data_str at the bottom of self.bests_file
        
        Deprecated function, use storing instead
        """
        ff = open(self.bests_file, 'a')
        ff.write(data_str)
        ff.close()

    def storing(self,t_gen,net):
        """Store the work and various data for later analysis
        
        Network object are stored in individual pickle file in Seed{}/data
        Data are stored in a shelve called the Seed{}/Bests_{}.net
        
        Args:
            t_gen: the key (normally the generation number)
            net (Network): the object to be saved
        
        Return:
            None
        """
        def add(data,key,value):
            data[key] = data.get(key,[])+[value]
        
        with shelve.open(self.data_file) as data:
            add(data,'generation',t_gen)
            add(data,'fitness',net.fitness)
            add(data,'n_interactions',len(net.list_types['Interaction']))
            add(data,'n_species',len(net.list_types['Species']))
            
        with open(self.bests_file.format(t_gen),'wb') as freezer:
            pickle.dump(net,freezer)

    def save_restart_file(self, kgeneration, header, tgeneration):
        """Save a dbm file, keyed by the generation number (a string!) and with value a
        [parameter dictionary, genus].  Might be more transparent to write out Poulation instance and
        forget header, and be sure to update tgeneration
        """
        rprmt = dict(header = header,
                     state = random.getstate(),
                     tgeneration = tgeneration)
        with shelve.open(self.restart_file) as restart_data:
            restart_data[str(kgeneration)] = (rprmt, self.genus)

        print('restart file saved after generation=', kgeneration, 'next tgeneration=', tgeneration)

    def pop_sort(self):
        """Sort the population with respect to fitness"""
        self.genus.sort(key=lambda X: X.fitness if X.fitness is not None else 9999)

    def update_fitness(self,nnetwork,integration_result):
        """Update (in place) the fitness and the dlt_fitness
        
        Args:
            nnetwork (int): the index of the network in the population
            integration_result (list): the output of compile_and_integrate
        
        Returns:
            None: in place modification
        """
        if integration_result:
            current_fitness = float(integration_result[0])
            self.genus[nnetwork].data_evolution = [data for data in integration_result] #stores various data on evolution
        else: #catches the None fitness
            current_fitness = None
            self.genus[nnetwork].data_evolution = None
            
        #update dlt_fitness and handle the None fitness option
        if current_fitness and self.genus[nnetwork].fitness:
            self.genus[nnetwork].dlt_fitness = current_fitness - self.genus[nnetwork].fitness
        elif current_fitness:
            self.genus[nnetwork].dlt_fitness = 9999
        elif self.genus[nnetwork].fitness:
            self.genus[nnetwork].dlt_fitness = -9999
        else:
            self.genus[nnetwork].dlt_fitness = 0
        self.genus[nnetwork].fitness = current_fitness

    def genus_mutate_and_integrate(self,prmt,nnetwork,mutation=True):
        """mutate, and update the fitness of one individual
        
        Args:
            prmt (dict): the inits parameters for integration
            nnetwork (int): the index of the network in the population
            mutation (bool): a flag to activate mutation
        
        Returns:
            int: the number of mutation
            int: the index of the network in the population
            Network: The resulting network after mutation
        """
        [n_mutations,nnetwork,mutated_net,result]=self.genus[nnetwork].mutate_and_integrate(prmt,nnetwork,self.tgeneration,mutation)
        self.update_fitness(nnetwork,result)
        self.n_mutations+=n_mutations
        return [n_mutations,nnetwork,mutated_net]

    def pop_mutate_and_integrate(self,initial,first_mutated,last_mutated,prmt,net_stat):
        """ Recompute the fitness for half the population and mutate/compute the fitness for the rest. Save all the data in net_stat
        
        Args:
            initial (int): index of the first individual in population
            first_mutated (int): index of the first mutated individual in population
            last_mutated (int): index of the last mutated individual in population
            prmt (dict): the inits parameters for integration
            net_stat (NetworkStat): to store the population data
        
        Returns:
            None: in place modification
        """
        self.n_mutations=0
        for nnetwork in range(initial,first_mutated):
            self.genus_mutate_and_integrate(prmt,nnetwork,mutation=False)
        for nnetwork in range(first_mutated,last_mutated):
            self.genus_mutate_and_integrate(prmt,nnetwork,mutation=True)
        for nnetwork in range(initial,last_mutated):
            net_stat.add_net(self.genus[nnetwork])
        return None

    def evolution(self):
        """Main method to evolve population
        
        Args:
            -
        
        Return:
            None
        """
        first_mutated = int( self.npopulation * (1-prmt['frac_mutate']) )
        net_stat = pop_stat.NetworkStat(stat_dict)
        gen_stat = pop_stat.GenusStat()

        #initialize attributs for each network needed in loop over generations
        if self.same_seed:
            print('Best fitness prior to mutations=', self.genus[0].fitness)
        else:
            self.pop_mutate_and_integrate(0,self.npopulation,self.npopulation-1,prmt,net_stat)
            for nnetwork in range( self.npopulation ):
                self.genus[nnetwork].data_next_mutation=self.genus[nnetwork].compute_next_mutation()
            self.pop_sort()
            print('Best/worst fitness prior to mutation=', self.genus[0].fitness, self.genus[-1].fitness)

        # MAIN EVOLUTIONARY LOOP
        for t_gen in range(self.generation0, prmt['ngeneration'] + self.generation0):
            net_stat = pop_stat.NetworkStat(stat_dict)
            gen_stat = pop_stat.GenusStat()

            # mutate a fraction of networks in population (those least fit)
            if (prmt['redo']==1): # in this case, recompute the fitness of non-mutated ind.
                self.pop_mutate_and_integrate(0,first_mutated,self.npopulation,prmt,net_stat)
            else: #only mutation
                self.pop_mutate_and_integrate(first_mutated,first_mutated,self.npopulation,prmt,net_stat)
            print("Total number of mutations in the population :%i"%self.n_mutations)
            
            # Adjust the tgeneration time to have roughly one mutation per individual in pop
            if (self.n_mutations>0):
                self.tgeneration=self.tgeneration*self.npopulation*prmt['frac_mutate']/self.n_mutations
            else:
                self.tgeneration=2*self.tgeneration

            fitness_treatment(self)
            self.pop_sort()
            gen_stat.process_sorted_genus( self )

            # print info after mutation step so built_integrator*.c consistent with Bests file
            header = "\nAfter generation {0:d} Best fitness={1}\ndata=[]\n".format(t_gen,self.genus[0].fitness)
            for data in self.genus[0].data_evolution:
                if data and len(data) > 0:
                    header=header+"data.append("+data+")\n"
            print(header)
            print("New generation time: %f"%self.tgeneration)
            sys.stdout.flush()
            self.storing(t_gen,self.genus[0])
            
            # Handling of different options
            if prmt['pareto'] == 1 and t_gen % prmt['freq_plot'] == 0:
                self.pop_print_pareto(self.namefolder+'/pareto'+str(t_gen),self.namefolder+'/rank1_nets'+str(t_gen))
            
            # Selection step, replace less fit networks by the fitter ones.
            for nnetwork in range( self.npopulation//2 ):
                self.genus[-1-nnetwork]=copy.deepcopy(self.genus[nnetwork]) # duplicates best half

            for individual in self.genus:
                new_seed = int(random.random()*100000) #generates new seed  to be sure not to overlap
                individual.Random=random.Random(new_seed) #reinitializes the random generator of every network

            # print statistics for this generation.  Need wrap all these variables into generic stat.
            if(t_gen%prmt['freq_stat'] == 0):
                net_stat.output()
                print("Total number of mutations: %i"%self.n_mutations)
                gen_stat.output()

            #plot the best network and its time course
            if (t_gen%prmt['freq_plot'] == 0) and (prmt['plot']==1):
                print("Plotting")
                plotdata.net_test_plot(self.genus[0], prmt, self.namefolder, t_gen)

            # save an exact copy of genus and relevant parameters for continuing loop
            if( t_gen%prmt['restart']['freq'] == 0):
                
                self.save_restart_file( t_gen, header, self.tgeneration )

            sys.stdout.flush()
