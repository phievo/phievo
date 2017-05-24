import pickle
import shelve
import glob,os,sys
import re
from phievo.AnalysisTools import main_functions as MF
from phievo.AnalysisTools  import palette
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
    def  __init__(self, path,mode="default"):
        """
        Creates a Simulation object. When th mode is not default (ex:test), the seeds are not loaded.
        usefull for prerun test.
        Args:
            path: directory of the project
            mode: Allows different mode to load the project.
        """
        self.root = path
        if self.root[-1] != os.sep:
            self.root += os.sep
        ## Upload Run parameters
        model_files = os.listdir(self.root)
        (model_dir , self.inits , init_file) =tuple(initialization_code.check_model_dir(self.root))
        self.inits.prmt["workplace_dir"] = os.path.join(self.inits.model_dir,"Workplace")
        self.deriv2 = initialization_code.init_networks(self.inits)
        self.plotdata = initialization_code.import_module(self.inits.pfile['plotdata'])
        if mode in ["default"]:
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

        try:
            palette.update_default_colormap(self.inits.prmt["palette"]["colormap"])
        except KeyError:
            pass
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

    def get_backup_net(self,seed,generation,index):
        """
            Get network from the backup file(or restart). In opposition to the best_net file
            the restart file is note stored at every generation but it contains a full
            population. This funciton allows to grab any individual of the population when
            the generation is stored

            Args:
                seed : index of the seed
                generation : index of the generation (must be a stored generation)
                index : index of the network within its generation
            Return:
                the selected network object
        """
        return self.seeds[seed].get_backup_net(generation,index)

    def stored_generation_indexes(self,seed):
        """
        Return the list of the stored generation indexes

        Args:
            seed (int): Index of Seed, you want the stored generation for.
        Return:
            list of the stored generation indexes
        """
        return self.seeds[seed].stored_generation_indexes()

    def run_dynamics(self,net=None,trial=1,erase_buffer=False,return_treatment_fitness=False):
        """
        Run Dynamics for the selected network. The function either needs the network as an argument or the seed and generation information to select it. If a network is provided, seed and generation are ignored.

        Args:
            net (Networks): network to simulate
            trial (int): Number of independent simulation to run
        Returns: data (dict) : dictionnary containing the time steps
            at the "time" key, the network at "net" and the corresponding
            time series for index of the trial.
             - net : Network
             - time : time list
             - outputs: list of output indexes
             - inputs: list of input indexes
             - 0 : data for trial 0
                - 0 : array for cell 0:
                       g0 g1 g2 g3 ..
                    t0  .
                    t1     .
                    t2        .
                    .
                    .
        """
        if net is None:
            net = self.seeds[seed].get_best_net(generation)
        self.inits.prmt["ntries"] = trial
        prmt = dict(self.inits.prmt)
        N_cell = prmt["ncelltot"]
        N_species = len(net.list_types['Species'])
        self.buffer_data = {"time":np.arange(0,prmt["dt"]*(prmt["nstep"]),prmt["dt"])}
        prmt["ntries"] = trial
        treatment_fitness = self.deriv2.compile_and_integrate(net,prmt,1000,True)
        col_select = np.arange(N_species)
        for i in range(trial):
            temp = np.genfromtxt('Buffer%d'%i, delimiter='\t')[:,1:]
            self.buffer_data[i] = {cell:temp[:,col_select + cell*N_species] for cell in range(N_cell)}
            if erase_buffer:
                os.remove("Buffer%d"%i)
            else:
                os.rename("Buffer{0}".format(i),os.path.join(self.root,"Buffer{0}".format(i)))

        self.buffer_data["net"] = net
        get_species = re.compile("s\[(\d+)\]")
        self.buffer_data["outputs"] = [int(get_species.search(species.id).group(1)) for species in net.list_types["Output"]]
        self.buffer_data["inputs"] = [int(get_species.search(species.id).group(1)) for species in net.list_types["Input"]]

        if return_treatment_fitness:
            return treatment_fitness
        return self.buffer_data


    def clear_buffer(self):
        """
        Clears the variable self.buffer_data.
        """
        self.buffer_data = None



    def Plot_TimeCourse(self,trial_index,cell=0,list_species=[]):
        """
        Searches in the data last stored in the Simulation buffer for the time course
        corresponding to the trial_index and the cell and plot the gene time series

        Args:
            trial_index: index of the trial you. Refere to run_dynamics to know how
            many trials there are.
            cell: Index of the cell to plot
        Return:
            figure
        """
        net = self.buffer_data["net"]
        nstep = self.inits.prmt['nstep']
        size = len(net.list_types['Species'])

        try:
            self.plotdata.Plot_Data(self.root+"Buffer%d"%trial_index,cell, size, nstep,list_species=list_species)
        except FileNotFoundError:
            print("Make sure you have run the function run_dynamics with the correct number of trials.")
            raise

    def Plot_Profile(self,trial_index,time=0):
        """
        Searches in the data last stored in the Simulation buffer for the time course
        corresponding to the trial_index and plot the gene profile along the cells at
        the selected time point.

        Args:
            trial_index: index of the trial you. Refere to run_dynamics to know how
            many trials there are.
            time: Index of the time to select
        Return:
            figure
        """
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

        self.restart_path = self.root + "Restart_file"
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

        with shelve.open(self.restart_path) as data:
            self.pop_size = len(data["0"][1])
            self.restart_generations = sorted([int(xx) for xx in data.dict.keys()])

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

        return MF.read_network(self.root+"Bests_%d.net"%generation,verbose=False)

    def get_backup_net(self,generation,index):
        """
        Get network from the backup file(or restart). In opposition to the best_net file
        the restart file is note stored at every generation but it contains a full
        population. This funciton allows to grab any individual of the population when
        the generation is stored

        Args:
            generation : index of the generation (must be a stored generation)
            index : index of the network within its generation
        Return:
            the selected network object
        """
        with shelve.open(self.restart_path) as data:
            dummy,nets = data[str(generation)]
        return(nets[index])

    def stored_generation_indexes(self):
        """
        Return the list of the stored generation indexes

        Return:
            list of the stored generation indexes
        """
        return list(self.restart_generations)

    def compute_best_fitness(self,generation):
        pass


class Seed_Pareto(Seed):
    def __init__(self,path,nbFunctions):
        super(Seed_Pareto, self).__init__( path)
        self.nbFunctions = nbFunctions

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

        with shelve.open(self.restart_path) as data:
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

    def plot_pareto_fronts(self,generations,with_indexes = False):
        """
            Plots the pareto fronts for a selected list of generations.
            Args:
                generations: list of generations indexes
        """

        ## Load fitness data for the selected generations and format them to be
        ## understandable by plot_multiGen_front2D
        data = MF.load_generation_data(generations,self.root+"Restart_file")

        generation_fitness = {}
        generation_indexes = {}
        for gen in generations:

            fitness_dico = {}
            index_dico = {}
            for i,ind in enumerate(data[gen]):
                try:
                    index_dico[ind.prank].append(i)
                    fitness_dico[ind.prank].append(list(ind.fitness))
                except KeyError:
                    index_dico[ind.prank] = [i]
                    fitness_dico[ind.prank] = [list(ind.fitness)]
            generation_fitness[gen] = fitness_dico
            generation_indexes[gen] = index_dico
        ## Obvious: plot
        fig = MF.plot_multiGen_front2D(generation_fitness,generation_indexes if with_indexes else None)
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
