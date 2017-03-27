import pickle
import shelve
import glob,os,sys
import re
from phievo.AnalysisTools import main_functions as MF
import matplotlib.pyplot as plt

import numpy as np
import importlib.util
from phievo.Networks import mutation,classes_eds2,plotdata
from phievo import  initialization_code

class Simulation:
    """
    The simulation class is a container in which the informations about a simulation are unpacked. This is used for easy access to a simulation results.
    """
    def  __init__(self, path,verbose=False):
        self.root = path
        if self.root[-1] != os.sep:
            self.root += os.sep

        ## Upload Run parameters
        model_files = os.listdir(self.root)
        (model_dir , self.inits , init_file) =tuple(initialization_code.check_model_dir(self.root))

        searchSeed  = re.compile("\d+$") ## integer ## at the end of the string "project_root/Seed##"
        seeds = [int(searchSeed.search(seed).group(0)) for seed in glob.glob(self.root+"Seed*")]
        seeds.sort()
        if self.inits.prmt["pareto"]:
            self.type = "pareto"
            nbFunctions = self.inits.prmt["npareto_functions"]
            self.seeds = {seed:Seed_Pareto(self.root+"Seed%d"%seed,nbFunctions=nbFunctions) for seed in seeds}
        else:
            self.type = "default"
            self.seeds = {seed:Seed(self.root+"Seed%d"%seed) for seed in seeds}

        self.buffer_data = None



    def show_fitness(self,seed,smoothen=0,**kwargs):
        """Plot the fitness as a function of time

        Args:
            seed (int): the seed-number of the run

        Returns:
            list: fitness as a function of time
        """
        fig = self.seeds[seed].show_fitness(smoothen,kwargs)
        return fig

    def custom_plot(self,seed,X,Y):
        """Plot the Y as a function of X. X and Y can be chosen in the list ["fitness","generation","n_interactions","n_species"]

        Args:
            seed (int): number of the seed to look at
            X (str): x-axis observable
            Y (str): y-axis observable
        """
        x_val = self.seeds[seed].custom_plot(X,Y)

    def get_best_net(self,seed,generation):
        """ The functions returns the best network of the selected generation

        Args:
            seed (int): number of the seed to look at
            generation (int): number of the generation

        Returns:
            Networks : the best network for the selected generation
        """

        return self.seeds[seed].get_best_net(generation)

    def run_dynamics(self,seed=None,generation=None,trial=1,net=None,erase_buffer=False):
        """
        Run Dynamics for the selected network. The function either needs the network as an argument or the seed and generation information to select it. If a network is provided, seed and generation are ignored.

        Args:
            seed (int): seed number
            generation (int): generation number
            trial (int): Number of independent simulation to run
            net (Networks): network to simulate
        Returns: data (dict) : dictionnary conatining the time steps
            at the "time" key and the corresponding time series for
            the indexes 0...nb_trials.

        """
        if net is None:
            net = self.seeds[seed].get_best_net(generation)
        self.buffer_data = net.run_dynamics(self.root,self.inits,trial=trial,erase_buffer=erase_buffer)
        self.buffer_data["net"] = net
        return self.buffer_data


    def clear_buffer(self):
        """
        Clears the variable self.buffer_data.
        """
        self.buffer_data = None



    def Plot_Data(self,trial_index,cell=0):
        net = self.buffer_data["net"]
        nstep = self.inits.prmt['nstep']
        size = len(net.list_types['Species'])
        try:
            plotdata.Plot_Data(self.root+"Buffer%d"%trial_index,cell, size, nstep)
        except FileNotFoundError:
            print("Make sure you have run the function run_dynamics with the correct number of trials.")
            raise

    def Plot_Profile(self,trial_index,time=0):
        net = self.buffer_data["net"]
        nstep = self.inits.prmt['nstep']
        size = len(net.list_types['Species'])
        ncelltot = self.inits.prmt['ncelltot']

        try:
            plotdata.Plot_Profile(self.root+"Buffer%d"%trial_index, ncelltot,size,time)
        except FileNotFoundError:
            print("Make sure you have run the function run_dynamics with the correct number of trials.")
            raise


class Seed:
    """
    This is a container to load the information about a Simulation seed. It contains mainly the indexes of the generations and some extra utilities to analyse them.
    """


    def  __init__(self, path):
        self.root = path
        if self.root[-1] != os.sep:
            self.root += os.sep
        self.name = re.search("[^/.]*?(?=/?$)",path).group(0)

        data = shelve.open(self.root+"data")
        indexes = data["generation"]
        interac = data['n_interactions']
        species = data['n_species']
        fitness = data['fitness']

        self.generations = {
            i:{
                "n_interactions" : interac[i],
                "n_species" : species[i],
                "fitness" : fitness[i]
                }
            for i in indexes}
        self.indexes = indexes



    def show_fitness(self,smoothen=0,**kwargs):
        """Plot the fitness as a function of time

        Args:
            seed (int): the seed-number of the run

        Returns:
            list: fitness as a function of time
        """
        gen = self.get_observable("generation")
        fit = self.get_observable("fitness")
        if smoothen:
            fit = MF.smoothing(fit,smoothen)
        fig = plt.figure()
        ax = fig.gca()
        ax.plot(gen,fit,lw=2,color='#B41111',label='fitness')
        ax.set_xlabel('Generation')
        ax.set_ylabel('Fitness')
        fig.show()
        return fig

    def custom_plot(self,X,Y):
        """Plot the Y as a function of X. X and Y can be chosen in the list ["fitness","generation","n_interactions","n_species"]

        Args:
            seed (int): number of the seed to look at
            X (str): x-axis observable
            Y (str): y-axis observable
        """
        print("Hey")
        x_val = self.get_observable(X)
        y_val = self.get_observable(Y)
        fig = plt.figure()
        ax = fig.gca()
        ax.plot(x_val,y_val,lw=2,color='#B41111',label='fitness')
        ax.set_xlabel(X)
        ax.set_ylabel(Y)
        fig.show()
        return fig

    def get_observable(self,observable):
        """
        Generates the list of all the values taken by the observable during the run. The observable is chosen in ["fitness","generation","n_interactions","n_species"].

        Args:
            str : name of the observable

        Returns:
            list : list of the fitnesses for the best run for every generation
        """

        if observable == "generation":
            return list(self.indexes)
        else:
            return [gen[observable] for i,gen in self.generations.items()]

    def get_best_net(self,generation):
        """ The functions returns the best network of the selected generation

        Args:
            seed (int): number of the seed to look at

        Returns:
            Networks : the best network for the selected generation
        """

        return classes_eds2.retrieve_from_pickle(self.root+"Bests_%d.net"%generation,verbose=False)

    def compute_best_fitness(self,generation):
        pass


class Seed_Pareto(Seed):
    def __init__(self,path,nbFunctions):
        super(Seed_Pareto, self).__init__( path)
        restart_path = self.root + "Restart_file"
        self.nbFunctions = nbFunctions
        with shelve.open(restart_path) as data:
            self.restart_generations = sorted([int(xx) for xx in data.dict.keys()])
        print(self.restart_generations)

    def show_fitness(self,smoothen=0,index=None):
        """Plot the fitness as a function of time

        Args:
            seed (int): the seed-number of the run
            index(array): index of of the fitness to plot. If None, all the fitnesses are ploted
        Returns:
            list: fitness as a function of time
        """
        gen = self.get_observable("generation")

        fit = np.array(self.get_observable("fitness"))

        if not index :
            index = range(fit.shape[1])



        fig = plt.figure()
        ax = fig.gca()

        for i in index:
            if smoothen:
                fit = MF.smoothing(fit,smoothen)

            ax.plot(gen,fit[:,i],lw=2,label='fitness%d'%i)
        ax.set_xlabel('Generation')
        ax.set_ylabel('Fitness')
        ax.legend(loc='upper right')
        fig.show()
        return fig


    def pareto_scatter(self,generation):
        """Display one generation as pareto fronts

        Run on 2D and 3D and use the Restart_file to retrieve whole population

        Args:
           generation (int): must be a whole generation saved (from restart file)
        Returns:
           dict: rank -> [[fitness1],[fitness2],…]
        """
        restart_path = self.root + "Restart_file"
        with shelve.open(restart_path) as data:
            dummy,pop_list = data[str(generation)]
        fitness_dico = {}

        fitnesses = self.get_observable("fitness")
        for ind in pop_list:
            try:
                fitness_dico[ind.prank].append([ind.fitness[i] for i in range(self.nbFunctions)])
            except KeyError:
                fitness_dico[ind.prank] = [[ind.fitness[i] for i in range(self.nbFunctions)]]

        if self.nbFunctions == 2:
            pareto_plane(fitness_dico,fitnesses)
        elif self.nbFunctions == 3:
            pareto_space(fitness_dico,fitnesses)
        else:
            print('Error, too many fitnesses to plot them all')
        return fitness_dico

## Functions
def pareto_plane(fitness_dico,fitnesses):
    """2d plotting subroutine of pareto_scatter"""
    for rank,points in fitness_dico.items():
        F1,F2 = list(zip(*points))
        #plt.plot(F1,F2,'k')
        plt.plot(F1,F2,'d',label='rank {}'.format(rank))
    for func,index in zip([plt.xlabel,plt.ylabel],fitnesses):
        func("fitness_{}".format(index))
    plt.legend(loc=0)
    plt.show()

def pareto_space(fitness_dico,fitnesses):
    """3d plotting subroutine of pareto_scatter"""
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for rank,points in fitness_dico.items():
        F1,F2,F3 = list(zip(*points))
        ax.plot(F1,F2,F3,label='rank {}'.format(rank))

    for func,index in zip([plt.xlabel,plt.ylabel,plt.zlabel],fitnesses):
        func("fitness_{}".format(index))
    plt.legend(loc=0)
    plt.show()
