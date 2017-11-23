import phievo.AnalysisTools
from phievo.AnalysisTools import Simulation
from phievo.AnalysisTools import Seed,Seed_Pareto
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
            plotly figure
        """    
    net_info,fitnesses = self.pareto_generate_fit_dict(generations,max_rank)
    if not colors and not gradient:colors = {gen:col for gen,col in zip(fitnesses.keys(),palette.color_generate(len(fitnesses)))}
    if gradient:
        colors = {gen:col for gen,col in zip(generations,palette.generate_gradient(generations,gradient))}
            
    shapes = ["circle","square","triangle-up"]
    data = []
    for gen,ranks in sorted(fitnesses.items(), key=lambda x:x[1]):
        for ind,rank in enumerate(ranks):
            if not rank: continue
            F1,F2 = zip(*rank)
            trace =  go.Scatter(x = F1,y = F2,mode = 'markers',name="G{0}-rank {1}".format(gen,ind),
                                marker= dict(size=14,color= colors[gen],symbol=shapes[ind]),
                                hoverinfo="text",
                                legendgroup = "group{}".format(gen)
                                
            )
            infos = net_info[gen][ind]
            label_template = Template("Generation: $gen<br>net: $net<br>rank: $rank<br>fitness: ($F1 , $F2)")
            trace["text"] = [label_template.substitute(inf) for inf in infos]
            
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

setattr(Seed_Pareto,"plot_pareto_fronts",plot_pareto_fronts)


def Plot_Profile(self,trial_index,time=0,no_popup=False):
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
    if not no_popup:
        plotfunc(fig)
    return fig
setattr(Simulation,"Plot_Profile",Plot_Profile)

def Plot_TimeCourse(self,trial_index,cell=0,no_popup=False):
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
    
    if not no_popup:
        plotfunc(fig)
    return fig
setattr(Simulation,"Plot_TimeCourse",Plot_TimeCourse)
