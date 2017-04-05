import pickle
import shelve
import glob,os,sys
import re
from phievo.AnalysisTools import main_functions as MF
import matplotlib.pyplot as plt
from matplotlib import pylab,colors
import matplotlib.patches as mpatches
from  matplotlib.lines import Line2D
import numpy as np
import importlib.util
from phievo.Networks import mutation,classes_eds2
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
        self.deriv2 = initialization_code.init_deriv2(self.inits,self.root,self.inits.prmt)
        self.plotdata = initialization_code.import_module(self.inits.pfile['plotdata'])

        searchSeed  = re.compile("\d+$") ## integer ## at the end of the string "project_root/Seed##"
        seeds = [int(searchSeed.search(seed).group(0)) for seed in glob.glob(os.path.join(self.root,"Seed*"))]
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
        self.inits.prmt["ntries"] = trial
        prmt = dict(self.inits.prmt)
        N_cell = prmt["ncelltot"]
        N_species = len(net.list_types['Species'])
        self.buffer_data = {"time":np.arange(0,prmt["dt"]*(prmt["nstep"]),prmt["dt"])}
        prmt["ntries"] = trial
        self.deriv2.compile_and_integrate(net,prmt,1000,True)
        for i in range(trial):
            temp = np.genfromtxt('Buffer%d'%i, delimiter='\t')[::,1:]
            self.buffer_data[i] = {cell:temp[::,cell:cell+N_species] for cell in range(N_cell)}
            if erase_buffer:
                os.remove("Buffer%d"%i)
            else:
                os.rename("Buffer{0}".format(i),os.path.join(self.root,"Buffer{0}".format(i)))
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
            self.plotdata.Plot_Data(self.root+"Buffer%d"%trial_index,cell, size, nstep)
        except FileNotFoundError:
            print("Make sure you have run the function run_dynamics with the correct number of trials.")
            raise

    def Plot_Profile(self,trial_index,time=0):
        net = self.buffer_data["net"]
        nstep = self.inits.prmt['nstep']
        size = len(net.list_types['Species'])
        ncelltot = self.inits.prmt['ncelltot']

        try:
            self.plotdata.Plot_Profile(self.root+"Buffer%d"%trial_index, ncelltot,size,time)
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

        self.observables = {
            "generation":lambda:list(self.indexes),
            "n_interactions":lambda:[gen["n_interactions"] for i,gen in self.generations.items()],
            "n_species":lambda:[gen["n_species"] for i,gen in self.generations.items()],
            "fitness":lambda:[gen["fitness"] for i,gen in self.generations.items()],
        }
        self.default_observable = "fitness"

    def show_fitness(self,smoothen=0,**kwargs):
        """Plot the fitness as a function of time
        """
        self.custom_plot("generation","fitness")



    def custom_plot(self,X,Y):
        """Plot the Y as a function of X. X and Y can be chosen in the keys of
            self.observables.

        Args:
            seed (int): number of the seed to look at
            X (str): x-axis observable
            Y (list): list (or string) of y-axis observable
        """
        x_val = self.observables[X]()
        if isinstance(Y,str):
            Y = [Y]
        Y_val = {y:self.observables[y]() for y in Y}

        NUM_COLORS = len(Y)
        cm = pylab.get_cmap('gist_rainbow')
        color_l= {Y[i]:colors.rgb2hex(cm(1.*i/NUM_COLORS)) for i in range(NUM_COLORS)}
        fig = plt.figure()
        ax = fig.gca()
        for label,y_val in Y_val.items():
            if y_val is not None:
                ax.plot(x_val,y_val,lw=2,color=color_l[label],label=label)
        ax.set_xlabel(X)
        ax.set_ylabel(Y[0])
        ax.legend()
        fig.show()
        return fig

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

        self.observables.pop("fitness")

        for ff in range(self.nbFunctions):
            self.observables["fitness{0}".format(ff)] = lambda ff=ff:[gen["fitness"][ff] for i,gen in self.generations.items()]
        self.default_observable = "fitness1"

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

    def plot_pareto_fronts(self,generations):
        """
            Plots the pareto fronts for a selected list of generations.

            Args:
                generations: list of generations indexes
        """

        ## Load fitness data for the selected generations and format them to be
        ## understandable by plot_multiGen_front2D
        data = load_generation_data(generations,self.root+"Restart_file")
        generation_fitness = {}
        for gen in generations:
            fitness_dico = {}
            for ind in data[gen]:
                try:
                    fitness_dico[ind.prank].append([ind.fitness[i] for i in range(self.nbFunctions)])
                except KeyError:
                    fitness_dico[ind.prank] = [[ind.fitness[i] for i in range(self.nbFunctions)]]
            generation_fitness[gen] = fitness_dico
        ## Obvious: plot
        fig = plot_multiGen_front2D(generation_fitness)
        return fig

## Functions
def load_generation_data(generations,restart_file):
    """
        Searches in the restart file the the informations that has been backed up
        up about the individuals at  a given generations.

        Args:
            generations (list): index of the generations to load_generation_data
            restart_file: path of the restart_file
        Returns:
            dictionary where each key contains the informations about one generation.
    """
    gen_data = {}

    with shelve.open(restart_file) as data:
        restart_generations = sorted([int(xx) for xx in data.dict.keys()])
        for gen in generations:
            if gen not in restart_generations:
                limit_print = 20
                err_str = ""
                err_str += "Generation {0} is not saved in the  restart file.\n".format(gen)
                err_str += "Please choose among the following generations:\n"
                if len(restart_generations)<limit_print:
                    err_str+=", ".join([str(x) for x in restart_generations[:limit_print]])+"\n"
                else:
                    err_str+=", ".join([str(x) for x in restart_generations[:limit_print]])+", ...\n"
                raise AssertionError(err_str)
            dummy,gen_data[gen] = data[str(gen)]
    return gen_data

def plot_multiGen_front2D(generation_fitness):
    """
        Uses the fitness data for multiple generations to represent the pareto fronts
        of those multiple generations.

        Args:
            generation_fitness: nested dictionnaries:
                                level0 keys: generation
                                level1 keys: rank of the fitness (1,2,...)
                                index : index of the fitness doublet (they might be
                                        multiple fitnesses with identical rank).

    """
    NUM_COLORS = len(generation_fitness)
    shapes = ["o","s","^"]
    cm = pylab.get_cmap('gist_rainbow')
    color_l= [colors.rgb2hex(cm(1.*i/NUM_COLORS)) for i in range(NUM_COLORS)]
    legend_patches = []
    #plt.legend(handles=[red_patch])
    i = 0
    fig = plt.figure()
    ax = fig.gca()
    for gen in sorted(generation_fitness.keys()):
        gen_dico = generation_fitness[gen]
        legend_patches.append(mpatches.Patch(color=color_l[i], label='Generation {0}'.format(gen)))
        color = color_l[i]
        i +=1
        for rank,points in gen_dico.items():
            F1,F2 = list(zip(*points))
            shape = shapes[rank-1] if rank<3 else shapes[-1]
            ax.scatter(F1,F2,c=color,edgecolor=color,s=50,marker=shape)
    ax.set_xlabel('Fitness 1')
    ax.set_ylabel('Fitness 2')
    ax.legend(handles=legend_patches)
    legend_patches.append(Line2D([0], [0], linestyle="none", marker=shapes[0], markersize=10,markerfacecolor="black",label="Rank 1"))
    legend_patches.append(Line2D([0], [0], linestyle="none", marker=shapes[1], markersize=10,markerfacecolor="black",label="Rank 2"))
    legend_patches.append(Line2D([0], [0], linestyle="none", marker=shapes[2], markersize=10,markerfacecolor="black",label="Rank≥3"))
    ax.legend(handles=legend_patches)
    plt.show()
    return fig




def pareto_plane(fitness_dico,fitnesses):
    """2d plotting subroutine of pareto_scatter"""
    ax = plt.gca()#(111, projection='polar')
    for rank,points in fitness_dico.items():
        F1,F2 = list(zip(*points))
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
