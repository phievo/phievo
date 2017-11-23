import pickle
import shelve
import glob,os,sys
import re
from phievo.AnalysisTools import main_functions as MF
from phievo.AnalysisTools  import palette
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from  matplotlib.lines import Line2D
import numpy as np
import importlib.util
from phievo.Networks import mutation,classes_eds2
from phievo import  initialization_code
import phievo

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
        setattr(phievo.Networks.mutation,"dictionary_ranges",self.inits.dictionary_ranges)

        self.deriv2 = initialization_code.init_networks(self.inits)
        self.plotdata = initialization_code.import_module(self.inits.pfile['plotdata'])
        if mode in ["default"]:
            searchSeed  = re.compile("Seed(.+)$") ## integer ## at the end of the string "project_root/Seed##"
            seeds = [searchSeed.search(seed).group(1)  for seed in glob.glob(os.path.join(self.root,"Seed*"))]

            seeds.sort()
            seeds = [int(ss) for ss in seeds if ss.isdigit()] + [ss for ss in seeds if not ss.isdigit()]          
            if self.inits.prmt["pareto"]:
                self.type = "pareto"
                nbFunctions = self.inits.prmt["npareto_functions"]
                self.seeds = {seed:Seed_Pareto(os.path.join(self.root,"Seed{}".format(seed)),nbFunctions=nbFunctions) for seed in seeds}
            else:
                self.type = "default"
                self.seeds = {seed:Seed(os.path.join(self.root,"Seed{}".format(seed))) for seed in seeds}

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
            matplotlib figure
        """
        fig = self.seeds[seed].show_fitness(smoothen,**kwargs)
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
            The best network for the selected generation
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
                The selected network object
        """
        return self.seeds[seed].get_backup_net(generation,index)

    def get_backup_pop(self,seed,generation):
        """
        Cf get_backup_net. Get the complete population of networks for a generation that
        was backuped.

        Args:
            seed : index of the seed
            generation : index of the generation (must be a stored generation)
        Return:
            List of the networks present in the population at the selected generation
        """
        return self.seeds[seed].get_backup_pop(generation)

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

        Returns: 
            data (dict) dictionnary containing the time steps
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
        N_species = len(net.dict_types['Species'])
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
        self.buffer_data["outputs"] = [int(get_species.search(species.id).group(1)) for species in net.dict_types["Output"]]
        self.buffer_data["inputs"] = [int(get_species.search(species.id).group(1)) for species in net.dict_types["Input"]]

        if return_treatment_fitness:
            return treatment_fitness
        return self.buffer_data


    def clear_buffer(self):
        """
        Clears the variable self.buffer_data.
        """
        self.buffer_data = None

    def PlotData(self,data,xlabel,ylabel,select_genes=[],no_popup=False,legend=True,lw=1,ax=None):
        """
        Function in charge of the call to matplotlib for both Plot_TimeCourse and Plot_Profile.
        """
        fig = plt.figure()
        if not ax:ax = fig.gca()
        Ngene = data.shape[1]
        colors = palette.color_generate(Ngene)
        for gene in range(Ngene):
            if select_genes!=[] and gene not in select_genes:
                continue
            ls = "--"
            label = "Species {}"
            if gene in self.buffer_data["outputs"]:
                ls = "-"
                label = "Output {}"
            if gene in self.buffer_data["inputs"]:
                ls = "-"
                label = "Input {}"
            ax.plot(data[:,gene],ls=ls,label=label.format(gene),lw=lw,color=colors[gene])
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if legend:ax.legend()
        return fig
    
    def Plot_TimeCourse(self,trial_index,cell=0,select_genes=[],no_popup=False,legend=True,lw=1,ax=None):
        """
        Searches in the data last stored in the Simulation buffer for the time course
        corresponding to the trial_index and the cell and plot the gene time series

        Args:
            trial_index: index of the trial you. Refere to run_dynamics to know how
            many trials there are.
            cell: Index of the cell to plot
            select_genes: list of gene indexes to plot
            no_popup: False by default. Option used to forbid matplotlib popup windows
                     Useful when saving figures to a file.
        Return:
            figure
        """
        data = self.buffer_data[trial_index]#[cell]
        try:
            data = data[cell]
        except KeyError:
            assert cell<0
            cell = sorted(list(data.keys()))[-1]
            data = data[cell]
        fig = self.PlotData(data,"Time","Concentration",select_genes=select_genes,no_popup=no_popup,legend=legend,lw=lw,ax=ax)        
        return fig

    def Plot_Profile(self,trial_index,time=0,select_genes=[],no_popup=False,legend=True,lw=1,ax=None):
        """
        Searches in the data last stored in the Simulation buffer for the time course
        corresponding to the trial_index and plot the gene profile along the cells at
        the selected time point.

        Args:
            trial_index: index of the trial you. Refere to run_dynamics to know how
            many trials there are.
            time: Index of the time to select
            select_genes: list of gene indexes to plot
            no_popup: False by default. Option used to forbid matplotlib popup windows
                     Useful when saving figures to a file.
        Return:
            figure
        """
        data = []
        for key,dynamics in sorted(self.buffer_data[trial_index].items()):
            data.append(dynamics[time,:])
        data = np.array(data)
        fig = self.PlotData(data,"Cell index","Concentration",select_genes=select_genes,no_popup=no_popup,legend=legend,lw=lw,ax=ax)
        return fig
        
    def load_Profile_data(self,trial_index,time):
        """
        Loads the data from the simulation and generate ready to plot data.
        Args:
            trial_index: index of the trial you. Refere to run_dynamics to know how
            many trials there are.
            time: Index of the time to select
        """
        data = []
        for key,dynamics in sorted(self.buffer_data[trial_index].items()):
            data.append(dynamics[time,:])
        return np.array(data)
        
    def get_genealogy(self,seed):
        return Genealogy(self.seeds[seed])

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
        color_l= {Y[i]:col for i,col in zip(range(NUM_COLORS),palette.color_generate(NUM_COLORS))}
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

    def get_backup_pop(self,generation):
        """
        Cf get_backup_net. Get the complete population of networks for a generation that
        was backuped.

        Args:
            generation : index of the generation (must be a stored generation)
        Return:
            List of the networks present in the population at the selected generation
        """
        with shelve.open(self.restart_path) as data:
            dummy,nets = data[str(generation)]
        return nets

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
            Matplolib figure
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

    def pareto_generate_fit_dict(self,generations,max_rank=1):
        """
        Load fitness data for the selected generations and format them to be
        understandable by plot_pareto_fronts
        """
        data = MF.load_generation_data(generations,self.root+"Restart_file")        
        fitnesses = {gen:
                     [
                         [net.fitness
                          for ind,net in enumerate(networks) if net.prank==rank+1]
                      for rank in range(min(max_rank,3))]                      
                     for gen,networks in data.items()}
        net_info  = {gen:
                     [
                         [dict(gen=gen,net=ind,rank=net.prank,F1=net.fitness[0],F2=net.fitness[1],ind=net.identifier)
                          for ind,net in enumerate(networks) if net.prank==rank+1]
                      for rank in range(min(max_rank,3))]                      
                     for gen,networks in data.items()}
        return net_info,fitnesses
    
    def plot_pareto_fronts(self,generations,max_rank=1,with_indexes=False,legend=False,xlim=[],ylim=[],colors=[],gradient=[],xlabel="F_1",ylabel="F_2",s=50,no_popup=False):
        """
        Plot every the network of the selected generations in the (F_1,F_2) fitness space.

        Args:
            generations (list): list of the selected generations
            max_rank (int): In given population plot only the network of rank <=max_rank
            with_indexes(bool): NotImplemented 
            legend(bool): NotImplemented
            xlim (list): [xmax,xmin]
            ylim (list): [ymax,ymin]
            colors (list): List of html colors, one for each generation
            gradient (list): List of colors to include in the gradient
            xlabel(str): Label of the xaxis
            ylabel(str): Label of the yaxis
            s (float): marker size
            no_popup(bool): prevents the popup of the plot windows 
        
        Returns:
            matplotlib figure
        """
        
        net_info,fitnesses = self.pareto_generate_fit_dict(generations,max_rank)
        if not colors and not gradient:colors = {gen:col for gen,col in zip(fitnesses.keys(),palette.color_generate(len(fitnesses)))}
        if gradient:
            colors = {gen:col for gen,col in zip(generations,palette.generate_gradient(generations,gradient))}
            
        shapes = ["o","s","^"]
        fig = plt.figure()
        ax = fig.gca()
        for gen,ranks in sorted(fitnesses.items(), key=lambda x:x[1]):
            for ind,rank in enumerate(ranks):
                if not rank: continue
                F1,F2 = zip(*rank)
                ax.scatter(F1,F2,s=s,color=colors[gen],marker=shapes[ind])                
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if xlim:ax.set_xlim(xlim)
        if ylim:ax.set_ylim(ylim)
        if not no_popup: plt.show()
        return fig
        
        

class Genealogy:
    def __init__(self,seed):
        self.generations = seed.stored_generation_indexes()
        assert self.generations == list(range(min(self.generations),max(self.generations)+1)),"\n\tA Genealogy object cannot be created from a seed where not all the generations were stored.\n\tMake sure prmt['restart']['freq'] = 1 in the init file."
        self.root = seed.root
        self.restart_path = seed.restart_path
        self.seed = seed
        self.networks = {}

    def sort_networks(self,verbose=False,write_pickle=True):
        """
        Order the networks, by the label_ind, in a dictionary.
        The dictonary contains the most useful information but takes last space.
        The information dictionaries is easier to handle than the actual networks.

        Args:
            verbose: print information during sorting
            write_pickle: backup the sorting information in a pickle file
        Return:
            dictionary. A key is associated to each network
        """
        networks = {}
        with shelve.open(self.restart_path) as data:            
            for gen in self.generations:
                dummy,population = data[str(gen)]
                
                for i,net in enumerate(population):
                    net_id = net.identifier
                    try:
                        networks[net_id]
                    except KeyError:                        
                        networks[net_id] = dict(
                            ind = net_id,
                            gen = gen,
                            par = net.parent,
                            fit = net.fitness,
                            pos = i,
                            las = net.last_mutation if net.last_mutation else [] 
                        )
                if verbose and (gen%100==0):
                    print("Generation\t{}/{} done.".format(gen,self.generations[-1]))
                    
        if write_pickle:
            with open(os.path.join(self.root,"networks_info.pkl"),"wb") as pkl_networks:
                pickle.dump(networks,pkl_networks)
        self.networks = networks
        return networks

    def load_sort_networks(self):
        """Loads an existing network classification"""
        with open(os.path.join(self.root,"networks_info.pkl"),"rb") as pkl_networks:
            networks = pickle.load(pkl_networks)
        self.networks = networks
        return networks

    def search_ancestors(self,network):
        if type(network) is int:
            network = self.networks[network]
        ancestors = [network]        
        while True:    
            try:
                network = self.networks[network["par"]]
                ancestors = [network] + ancestors
            except KeyError:
                break
        return ancestors
    
    
    def plot_front_genealogy(self,generations,extra_networks_info=[],filename=""):
        """
        Uses the seed plot_pareto_fronts function to display the pareto fronts.
        In addition, the function allows to plots extra networks in the fitness plan
        
        Args:
            generations: list of generation indexes
            extra_networks_indexes: list of extra network informatino dictionaries.

        """
        from phievo.AnalysisTools import plotly_graph
        fig = self.seed.plot_pareto_fronts(generations,with_indexes=True)
        if extra_networks_info:
            fit0 = [net_inf["fit"][0] for net_inf in extra_networks_info]
            fit1 = [net_inf["fit"][1] for net_inf in extra_networks_info]

            trace = plotly_graph.go.Scatter(x=fit0,y=fit1,mode = 'markers',name="Extra networks",
                                            marker= dict(size=9,color= "black",symbol="square"),
                                            hoverinfo="text",
                                            legendgroup = "Extra networks",
                                            text=["net #{}\nmutation:{}".format(net_inf["ind"]," - ".join(net_inf["las"] if net_inf["las"] else [])) for net_inf in extra_networks_info]
            )
            fig.data.append(trace)
            if filename:
                plotly_graph.py.plot(fig,filename=filename)
            else:
                plotly_graph.py.plot(fig)

    def plot_mutation_fitness_deviation(self,only_one_mutation=True,networks=None,ploted_ratio=1):
        """
        Plot the deviation of fitness in the fitness space caused by a generation's mutation.
        
        Arg:
           only_one_mutation (bool): If True, plot only the networks that undergone only a single mutation durign a generation.

        """
        from phievo.AnalysisTools import plotly_graph
        dict_data = {}
        if networks:
            if type(networks)==list:
                networks = {net["ind"]:net  for net in networks}
        else:
            networks = self.networks
            
        for net_ind,net_inf in networks.items():
            if net_inf["las"] is None or (only_one_mutation and len(net_inf["las"])!=1):
                continue
            label = "-".join(net_inf["las"])
            try:
                parent = networks[net_inf["par"]]
            except KeyError:
                ## Parent not provided in the dict
                continue
            diff = [net_inf["fit"][0]-parent["fit"][0],net_inf["fit"][1]-parent["fit"][1]]
            data =  {"diff":diff,"label":"Net #{}\nparent #{}\nmutation: {}\nfitness change:{}".format(net_inf["ind"],parent["ind"],label,str(diff)) }
            dict_data[label] = dict_data.get(label,[]) + [data]

        plot_list = []
        colors = palette.color_generate(len(dict_data))

        for mut_ind,mut_name in enumerate(dict_data.keys()):
            mutation = dict_data[mut_name]
            x_val = [mut["diff"][0] for mut in mutation]
            y_val = [mut["diff"][1] for mut in mutation]
            hover_info = [mut["label"] for mut in mutation]
            L = len(x_val)
            if ploted_ratio!=1:
                selected_indexes = np.random.choice(range(L),int(L*ploted_ratio),replace=False)
                x_val = np.array(x_val)[selected_indexes]
                y_val = np.array(y_val)[selected_indexes]
                hover_info = np.array(hover_info)[selected_indexes]
            trace = plotly_graph.go.Scatter(x = x_val,y=y_val,mode = 'markers',name=mutation,
                                            marker= dict(size=14,color=colors[mut_ind]),
                                            hoverinfo="text",
                                            text=hover_info,
            )
            plot_list.append(trace)
            
        plotly_graph.py.plot(plot_list)

    def get_network_from_identifier(self,net_ind):
        try:
            network = self.networks[net_ind]
        except KeyError:
            raise KeyError("Index {} corresponds to no stored network.".format(net_ind))
        net = self.seed.get_backup_net(network["gen"],network["pos"])
        assert net.identifier==net_ind,"\tBug in get_network_from_identifier.\n\tThe obtained network has indec {} instead of {}.".format(net.identifier,net_ind)
        return net
    
    def plot_lineage_fitness(self,line,formula="{}",highlighted_mutations = []):
        from phievo.AnalysisTools import plotly_graph
        generations = [net_inf["gen"] for net_inf in line]
        fitness = [eval(formula.format(net_inf["fit"])) for net_inf in line]
        generate_hoverinfo = lambda net:"net #{}\nmutations:{}".format(net["ind"],"-".join(net["las"] if net["las"] else [] ))
        hover_info = [generate_hoverinfo(net_inf)  for net_inf in line]
        colors = palette.color_generate(len(highlighted_mutations)+1)
        trace = plotly_graph.go.Scatter(
            x = generations,
            y = fitness,
            text = hover_info,
            mode = "markers+lines",
            marker = dict(color=colors[0]),
        )
        data_to_plot = [trace]
        if highlighted_mutations:
            dict_highlighted = {mut:plotly_graph.go.Scatter(
                x = [],
                y = [],
                text = [],
                mode = "markers",
                marker = dict(color=colors[i+1]),
            ) for i,mut in enumerate(highlighted_mutations)}

            for net_inf in line:
                try:
                    ind = "-".join(net_inf["las"]) if net_inf["las"] else ""
                    trace = dict_highlighted[ind]
                except KeyError:
                    continue
                trace["x"].append(net_inf["gen"])
                trace["y"].append(eval(formula.format(net_inf["fit"])))
                trace["text"].append(generate_hoverinfo(net_inf))
            data_to_plot += list(dict_highlighted.values())
        plotly_graph.py.plot(data_to_plot)
        
    def plot_compare_multiple_networks(self,sim,indexes,cell=0):
        """
        Print a svg figure of the cell profile,time series and the network layout in
        the seed folder.
        """
        
        fig_format = "svg"
        for filename in glob.glob(os.path.join(self.root,"*svg")):
            os.remove(filename)
        fig_name = lambda xx,ind: os.path.join(self.root,"{}{}.{}".format(xx,ind,fig_format))
        for i,net_ind in enumerate(indexes):
            net = self.get_network_from_identifier(net_ind)
            net.draw(fig_name("net",net.identifier))
            res = sim.run_dynamics(net)
            fig_profile=sim.Plot_Profile(0,time=2999,no_popup=True)
            fig_profile.savefig(fig_name("profile",net.identifier))
            fig_time_course=sim.Plot_TimeCourse(0,cell=cell,no_popup=True)
            fig_time_course.savefig(fig_name("timecourse",net.identifier))        
            del fig_profile,fig_time_course

    def compare_ss_wrt_parent(self,sim,child,parent):
        from phievo.AnalysisTools import plotly_graph
        child =  self.get_network_from_identifier(child)
        parent = self.get_network_from_identifier(parent)
        res = sim.run_dynamics(child)
        child_profile = sim.Plot_Profile(0,time=2999,no_popup=True)
        res = sim.run_dynamics(parent)
        parent_profile = sim.Plot_Profile(0,time=2999,no_popup=True)
        for dat_par,dat_chi in zip(child_profile["data"],parent_profile["data"]):
            chi_y = np.array(dat_chi["y"])
            par_y = np.array(dat_par["y"])
            relative_y =  np.array(chi_y - par_y)/np.where(par_y==np.array(0),1,par_y)
            dat_chi["y"] = relative_y
        plotly_graph.py.plot(parent_profile)
        
    def scatter_pareto_accross_generations(self,generation,front_to_plot,xrange,yrange,step=1):
        from phievo.AnalysisTools import plotly_graph
        plotly_graph.py.init_notebook_mode(connected=True)
        generation_identifiers = [net.identifier for net in self.seed.get_backup_pop(generation)]
        generation_info = [self.networks[net_ind] for net_ind in generation_identifiers]
        
        data_to_plot = []
        steps = []
        full_scatter =dict(
            name = generation,
            mode = "markers",
            x = [],
            y = [],
            marker = dict(color="black",size=2),
        )
        for gen in front_to_plot:
            generation_identifiers = [net.identifier for net in self.seed.get_backup_pop(gen)]
            networks = [self.networks[net_ind] for net_ind in generation_identifiers]
            fit0 = [net["fit"][0] for net in networks]
            fit1 = [net["fit"][1] for net in networks]
            full_scatter["x"]+=fit0
            full_scatter["y"]+=fit1
            
        while True:
            for net_pos in range(len(generation_info)):
                net_inf = generation_info[net_pos]
                
                while net_inf["gen"] > generation:
                    generation_info[net_pos] = self.networks[net_inf["par"]]
                    assert generation_info[net_pos]["ind"] == net_inf["par"]                    
                    assert generation_info[net_pos], "\n\tscatter_pareto_accross_generations found no parent for network #{}(generation {}).".format(net_inf["ind"],generation)
                    net_inf = self.networks[net_inf["par"]]
                    
            fit0 = [net_inf["fit"][0] for net_inf in generation_info]
            fit1 = [net_inf["fit"][1] for net_inf in generation_info]
            networks_info = ["#{} parent:#{} fitness:{}".format(net_inf["ind"],net_inf["par"],net_inf["fit"].__str__()) for net_inf in generation_info]
            generation_dict =dict(
                name = generation,
                mode = "markers",
                x = fit0,
                y = fit1,
                text=networks_info
            )

            slider_step = {'args': [
                [generation],
                {'frame': {'duration': 300, 'redraw': False},
                 'mode': 'immediate',
                 'transition': {'duration': 300}}
            ],
                           'label': generation,
                           'method': 'animate'}
            data_to_plot.append(generation_dict)
            steps.append(slider_step)
            
            
            generation -= step
            if generation < 0:
                break
        figure = {
            'data': [full_scatter,data_to_plot[0]],
            'layout': {
                'xaxis' : {'range': xrange, 'title': 'Fitness 1'},
                'yaxis' : {'range': yrange, 'title': 'Fitness 2'},
                'hovermode':'closest',
                'updatemenus':[
                    {
                        'buttons': [
                            {
                                'args': [None, {'frame': {'duration': 500, 'redraw': False},
                                                'fromcurrent': True, 'transition': {'duration': 300, 'easing': 'quadratic-in-out'}}],
                                'label': 'Play',
                                'method': 'animate'
                            },
                            {
                                'args': [[None], {'frame': {'duration': 0, 'redraw': False}, 'mode': 'immediate',
                                                  'transition': {'duration': 0}}],
                                'label': 'Pause',
                                'method': 'animate'
                            }
                        ],
                        'direction': 'left',
                        'pad': {'r': 10, 't': 87},
                        'showactive': False,
                        'type': 'buttons',
                        'x': 0.1,
                        'xanchor': 'right',
                        'y': 0,
                        'yanchor': 'top'
                    } 
                ],
                'sliders':[{
                    'active': 0,
                    'yanchor': 'top',
                    'xanchor': 'left',
                    'currentvalue': {
                        'font': {'size': 20},
                        'prefix': 'Generation:',
                        'visible': True,
                        'xanchor': 'right'
                    },
                    'transition': {'duration': 300, 'easing': 'cubic-in-out'},
                    'pad': {'b': 10, 't': 50},
                    'len': 0.9,
                    'x': 0.1,
                    'y': 0,
                    'steps': steps[::-1],
                }]
            },
            'frames': [{"data":[full_scatter,dat],"name":dat["name"]}
                for dat in data_to_plot[::-1]
            ],
        }
        
        
        plotly_graph.py.plot(figure, filename='pareto_accross_generations.html')
        return figure
        
