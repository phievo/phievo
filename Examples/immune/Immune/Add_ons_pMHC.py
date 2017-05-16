from phievo.AnalysisTools.Notebook import Notebook,CellModule
from phievo.AnalysisTools import Simulation
from IPython.display import display,HTML,clear_output, Image
from  ipywidgets import widgets
import numpy as np
import glob,os,sys
import re
import matplotlib.pyplot as plt

from importlib import import_module
pretty_graph = import_module('immune.Immune.pretty_graph_2_pMHC')


def Plot_pMHC(self):
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
    pMHC_size = len(net.list_types['pMHC'])
    ntau=len(self.inits.prmt['tau_off'])
    #print(size)
    #print(pMHC_size)
    total_size=size+pMHC_size+1
    ncelltot=total_size*ntau
    #ncelltot = self.inits.prmt['ncelltot']
    list_ligands=self.inits.prmt['Ligands']
    outputs=self.buffer_data["outputs"]
    trial_index=0
    try:
        self.plotdata.Plot_pMHC(self.root+"Buffer%d"%trial_index, ncelltot,total_size, ntau,list_ligands,list_output=outputs)
    except FileNotFoundError:
        print("Make sure you have run the function run_dynamics with the correct number of trials.")
        raise
        



def run_dynamics_pMHC(self,net=None,trial=1,erase_buffer=False,return_treatment_fitness=False):
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
    #if return_treatment_fitness:
    #    return treatment_fitness
    #col_select = np.arange(N_species)
    #for i in range(trial):

    #    temp = np.genfromtxt('Buffer%d'%i, delimiter='\t')[:,1:]
    #    self.buffer_data[i] = {cell:temp[:,col_select + cell*N_species] for cell in range(N_cell)}
    #    if erase_buffer:
    #        os.remove("Buffer%d"%i)
    #    else:
    i=0
    os.rename("Buffer{0}".format(i),os.path.join(self.root,"Buffer{0}".format(i)))
    self.buffer_data["net"] = net
    get_species = re.compile("s\[(\d+)\]")

    self.buffer_data["outputs"] = [int(get_species.search(species.id).group(1)) for species in net.list_types["Output"]]
    self.buffer_data["inputs"] = [int(get_species.search(species.id).group(1)) for species in net.list_types["Input"]]

    return self.buffer_data
    
    
setattr(Simulation,'Plot_pMHC',Plot_pMHC)
setattr(Simulation,'run_dynamics_pMHC',run_dynamics_pMHC)


class Run_Dynamics_pMHC(CellModule):
    def __init__(self,Notebook):
        super(Run_Dynamics_pMHC, self).__init__(Notebook)
        self.widget_nputs = widgets.IntText(value=1,description='N :',disabled=False)
        self.button_launchRun = widgets.Button(description="Run dynamics",disabled=True)
        self.notebook.dependencies_dict["generation"].append(self)
        self.notebook.dependencies_dict["dynamics"] = []
        self.notebook.extra_variables = {"ntries":None}
    def launch_dynamics(self,button):
        self.notebook.sim.run_dynamics_pMHC(net=self.notebook.net,erase_buffer=False,trial=self.widget_nputs.value)
        self.notebook.extra_variables["ntries"] = self.widget_nputs.value
        for cell in self.notebook.dependencies_dict["dynamics"]:
            cell.update()
    def update(self):
        if self.notebook.generation is None:
            self.button_launchRun.disabled = True
            self.notebook.sim.buffer_data = None
            self.notebook.extra_variables["ntries"] = None
        else:
            self.notebook.extra_variables["ntries"] = None
            self.notebook.sim.buffer_data = None
            self.button_launchRun.disabled = False
        for cell in self.notebook.dependencies_dict["dynamics"]:
            cell.update()

    def display(self):
        self.button_launchRun.on_click(self.launch_dynamics)
        display(widgets.HBox([self.widget_nputs,self.button_launchRun]))

class Plot_pMHC(CellModule):
    def __init__(self,Notebook):
        super(Plot_pMHC, self).__init__(Notebook)
        self.notebook.dependencies_dict["dynamics"].append(self)
        self.notebook.dependencies_dict["generation"].append(self)
        #self.widget_selectInput = widgets.IntSlider(value = 0,min=0,max=0,description = 'Input:',disabled=True)
        #self.widget_selectCell = widgets.IntSlider(value = 0,min=0,max=0,description = 'Cell:',disabled=True)
        self.button_plotdynamics = widgets.Button(description="Plot Response Curve",disabled=True)

    def plot_pMHC(self,button):
        plt.close()
        clear_output()
        self.notebook.sim.Plot_pMHC()

    def update(self):
       self.button_plotdynamics.disabled = False

    def display(self):
        self.button_plotdynamics.on_click(self.plot_pMHC)
        display(widgets.HBox([self.button_plotdynamics]))
        
class Plot_Layout_Immune(CellModule):
    def __init__(self,Notebook):
        super(Plot_Layout_Immune, self).__init__(Notebook)
        self.notebook.dependencies_dict["seed"].append(self)
        self.notebook.dependencies_dict["generation"].append(self)
        self.button_plotLayout = widgets.Button(description="Plot network layout",disabled=True)

    def plot_layout_immune(self,button):
        plt.close()
        clear_output()
        network=self.notebook.net
        #print(self.notebook.sim.inits)
        graph=pretty_graph.pretty_graph(network)
        #print(graph)
        graph.write_png('current_graph.png')
        display(Image(filename='current_graph.png'))

    def update(self):
        if self.notebook.generation is None:
            self.button_plotLayout.disabled = True
        else:
            self.button_plotLayout.disabled = False

    def display(self):
        self.button_plotLayout.on_click(self.plot_layout_immune)
        display(self.button_plotLayout)

