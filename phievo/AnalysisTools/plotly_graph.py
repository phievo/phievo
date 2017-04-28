import phievo.AnalysisTools
from phievo.AnalysisTools import Simulation
from phievo.AnalysisTools import Seed
from phievo.AnalysisTools import main_functions as MF
from phievo.AnalysisTools import palette
import numpy as np
import plotly.offline as py
py.init_notebook_mode()
import plotly.graph_objs as go
import plotly.tools as tls
from matplotlib import pylab,colors
from string import Template


linewidth = 5
plotfunc = py.plot

def run_in_nb():
    global plotfunc
    plotfunc = py.iplot

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
    data = []
    for label,y_val in Y_val.items():
        trace = go.Scatter(
            name=label,
            x = x_val,
            y = y_val,
            line = dict(width=linewidth)
        )
        data.append(trace)
        #clear_output()
    layout = go.Layout(
    title="Custom plot",
        xaxis=dict(title=X),
        yaxis=dict(title=Y[0]),
        )
    fig = go.Figure(data=data, layout=layout)
    plotfunc(fig)
    return fig
setattr(Seed,"custom_plot",custom_plot)

def plot_multiGen_front2D(generation_fitness,generation_indexes=None):
    """
        Uses the fitness data for multiple generations to represent the pareto fronts
        of those multiple generations.

        Args:
            generation_fitness: nested dictionnaries:
                                level0 keys: generation
                                level1 keys: rank of the fitness (1,2,etc.)
                                index : index of the fitness doublet (they might be
                                        multiple fitnesses with identical rank).
            generation_indexes: Same dictionnary structure as generation_fitness.
                                Contains the index of each network in its population

    """
    NUM_COLORS = len(generation_fitness)
    # https://plot.ly/python/reference/
    shapes = ["circle","square","triangle-up"]
    color_l = palette.color_generate(NUM_COLORS)
    legend_patches = []

    i = 0
    data = []
    for gen in sorted(generation_fitness.keys()):
        gen_dico = generation_fitness[gen]
        #legend_patches.append(mpatches.Patch(color=color_l[i], label='Generation {0}'.format(gen)))
        color = color_l[i]
        i +=1
        for rank,points in gen_dico.items():
            F1,F2 = list(zip(*points))
            shape = shapes[rank-1] if rank<3 else shapes[-1]
            trace =  go.Scatter(x = F1,y = F2,mode = 'markers',name="G{0}-rank {1}".format(gen,rank),
                marker= dict(size=14,color= color,symbol=shape),
                hoverinfo="text",
                legendgroup = "group{}".format(gen)
                )
            #ax.scatter(F1,F2,c=color,edgecolor=color,s=50,marker=shape)
            if generation_indexes:
                str_label = Template("Generation: $gen\nnetwork: $net\nrank: $rank\nfitness: ($F1 , $F2)")
                ind_list = generation_indexes[gen][rank]
                label_dict = { l:dict(gen=gen,net=ind_list[l],rank=rank,F1=F1[l],F2=F2[l]) for l in range(len(ind_list))}

                trace["text"]= [str_label.substitute(label_dict[l]) for l in range(len(ind_list))]
            if generation_indexes and False:
                ind_list = generation_indexes[gen][rank]
                for l in range(len(ind_list)):
                    ax.text(F1[l],F2[l],'%d' % ind_list[l],ha='center', va='bottom')
            data.append(trace)
    layout = go.Layout(
        title="Pareto fronts",
        xaxis=dict(title="Fitness 1"),
        yaxis=dict(title="Fitness 2"),
        autosize=True,
        hovermode='closest',
        )
    fig = go.Figure(data=data, layout=layout)
    plotfunc(fig)
    return fig
#import pdb;pdb.set_trace()
setattr(MF,"plot_multiGen_front2D",plot_multiGen_front2D)


def Plot_Profile(self,trial_index,time=0):
    """
    Searches in the data last stored in the Simulation buffer for the cell profile
    corresponding to the time point "time" and plot the profile.

    Args:
        trial_index: index of the trial you. Refere to run_dynamics to know how
        many trials there are.
        time: Index of the time point to plot
    Return:
        figure
    """

    data = []
    nstep = self.inits.prmt['nstep']
    nSpecies = self.buffer_data[trial_index][0].shape[1]
    nCell = len(self.buffer_data[trial_index])
    cells = list(range(nCell))
    self.buffer_data[trial_index][0]
    profiles = np.array([self.buffer_data[trial_index][cc][time,:] for cc in range(nCell)]).T

    species_labels = ["Species {}".format(i) for i in range(nSpecies)]
    for ind in self.buffer_data["outputs"]:
        species_labels[ind] = species_labels[ind].replace("Species","Output")
    for ind in self.buffer_data["inputs"]:
        species_labels[ind] = species_labels[ind].replace("Species","Inputs")
    dt = self.inits.prmt["dt"]
    colors = palette.color_generate(nSpecies)
    for i in range(nSpecies):
        trace = go.Scatter(
            x=cells,
            y=profiles[i,:],
            name=species_labels[i],
            line = dict(color=colors[i],width=linewidth)
        )
        data.append(trace)
    layout = go.Layout(
        title="Cell profile",
        xaxis=dict(title="Cell, time={}".format(time)),
        yaxis=dict(title="Concentration"),
        autosize=True,
        hovermode='closest',
        )
    fig = go.Figure(data=data, layout=layout)
    plotfunc(fig)
    return fig
setattr(Simulation,"Plot_Profile",Plot_Profile)

def Plot_TimeCourse(self,trial_index,cell=0):
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

    data = []
    nstep = self.inits.prmt['nstep']
    time_course = self.buffer_data[trial_index][cell]
    nSpecies = time_course.shape[1]
    species_labels = ["Species {}".format(i) for i in range(nSpecies)]
    for ind in self.buffer_data["outputs"]:
        species_labels[ind] = species_labels[ind].replace("Species","Output")
    for ind in self.buffer_data["inputs"]:
        species_labels[ind] = species_labels[ind].replace("Species","Inputs")
    dt = self.inits.prmt["dt"]
    nstep = self.inits.prmt["nstep"]
    time_axis = np.arange(0,nstep*dt,dt)
    colors = palette.color_generate(nSpecies)
    for i in range(nSpecies):
        trace = go.Scatter(
            x=time_axis,
            y=time_course[:,i],
            name=species_labels[i],
            line = dict(color=colors[i],width=linewidth)
        )
        data.append(trace)
    layout = go.Layout(
        title="Time course",
        xaxis=dict(title="Time, cell={}".format(cell)),
        yaxis=dict(title="Concentration"),
        autosize=True,
        hovermode='closest',
        )
    fig = go.Figure(data=data, layout=layout)
    plotfunc(fig)
    return fig
setattr(Simulation,"Plot_TimeCourse",Plot_TimeCourse)
