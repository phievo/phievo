import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
from  ipywidgets import widgets
from ipywidgets import interact, interactive, fixed
from IPython.display import display,HTML,clear_output
import os
from phievo.AnalysisTools import Simulation
## Load plotly_graph in order to use plotly in the notebook
import phievo.AnalysisTools.plotly_graph as plotly_graph
plotly_graph.run_in_nb()

found_str = "<p style=\"color:#31B404;font-size: 30px;\">✔</p>"
notfound_str = "<p style=\"color:#DF3A01;font-size: 30px;\">✘</p>"

class Notebook(object):
    """
        Wrapper that contains both the the widgets and  the simulation results.
        This way it is easy to update the state of the widgets when you load a
        new simulation
    """
    def __init__(self):
        self.sim = None
        rcParams['figure.figsize'] = (9.0, 8.0)
        self.project = None
        self.dependencies_dict = {
            "project" : [],
            "seed" : [],
            "generation" : []
        } ## List of cell objects to update when a key change changes. New dependencies
        ## may be added with new functions.
        self.seed = None
        self.generation = None
        self.net = None
        self.type = None ## Ex of type: "pareto"
        self.extra_variables = {} ## Allows to add new variables with a new cell
        ## object
        self.select_project = Select_Project(self)
        self.select_seed = Select_Seed(self)
        self.plot_pareto_fronts = Plot_Pareto_Fronts(self)
        self.plot_evolution_observable = Plot_Evolution_Observable(self)
        self.select_generation = Select_Generation(self)
        self.plot_layout = Plot_Layout(self)
        self.run_dynamics = Run_Dynamics(self)
        self.plot_dynamics = Plot_Dynamics(self)
        self.plot_cell_profile = Plot_Cell_Profile(self)

class CellModule():
    """
    Template class from which a module should inheritate.
    """
    def __init__(self,Notebook):
        self.notebook = Notebook
    def display(self):
        raise NotImplemented("Please define a display function for your module.")
    def update(self):
        raise NotImplemented("Please define a update function for your module.")

class Select_Project(CellModule):
        def __init__(self,Notebook):
            super(Select_Project, self).__init__(Notebook)
            self.widget_select_project = widgets.Text(value='',placeholder='Name of project directory',description='Directory:',disabled=False)
            self.widget_loadDir = widgets.Button(description="Load Run",disabled=True)
            self.foundFile_widget = widgets.HTML("")

        def inspect_run(self,path):
            """
            Test if the dir name exists

            Args:
                path (str): path of the directory
            """
            self.update()
            if os.path.isdir(path) and not self.notebook.project:
                self.foundFile_widget.value = found_str
                self.widget_loadDir.disabled=False
            else:
                self.foundFile_widget.value = notfound_str
                self.widget_loadDir.disabled=True

        def load_project(self,button):
            """
            Load the project directory in a somulation object.

            Args:
                directory : Path of the project
            """
            self.notebook.sim = Simulation(self.widget_select_project.value)
            self.notebook.project = self.widget_select_project.value
            self.widget_loadDir.button_style = "success"
            self.widget_loadDir.disabled=True
            self.notebook.type = self.notebook.sim.type
            print("To load a new project, please restart the kernel")
            print("Kernel > Restart & Run All")
            for cell in self.notebook.dependencies_dict["project"]:
                cell.update()


        def display(self):
            self.widget_loadDir.on_click(self.load_project)
            interactive(self.inspect_run,path=self.widget_select_project);
            main_options = widgets.VBox([widgets.HBox([self.widget_select_project,self.foundFile_widget]),self.widget_loadDir])
            display(main_options)
            #display(self.widget_select_project)

        def update(self):
            """
                Clears what needs to be cleared when the directory is changed.
            """
            self.widget_loadDir.disabled = True
            self.widget_loadDir.button_style = ""

class Select_Seed(CellModule):
    def __init__(self,Notebook):
        super(Select_Seed, self).__init__(Notebook)
        self.notebook = Notebook
        self.widget_select_seed = widgets.Dropdown(options={"None":None},value=None,description='Seed:',disabled=True)
        self.notebook.dependencies_dict["project"].append(self)

    def read_seed(self,seed_name):
        self.notebook.seed =  seed_name
        for cell in self.notebook.dependencies_dict["seed"]:
            cell.update()

    def display(self):
        interactive(self.read_seed,seed_name=self.widget_select_seed)
        display(self.widget_select_seed)

    def update(self):
        if self.notebook.project is None:
            self.notebook.seed = None
            self.widget_select_seed.options = {"None":None}
            self.widget_select_seed.disabled = True
            self.widget_select_seed.value = None
        else:
            self.widget_select_seed.disabled = False
            self.widget_select_seed.options = {"Seed {}".format(i):i for i,seed in self.notebook.sim.seeds.items()}            
            self.widget_select_seed.value = self.widget_select_seed.options[list(self.widget_select_seed.options.keys())[0]]
            self.notebook.seed = self.widget_select_seed.value

class Plot_Evolution_Observable(CellModule):
    def __init__(self,Notebook):
        super(Plot_Evolution_Observable, self).__init__(Notebook)
        self.widget_Xobs = widgets.Dropdown(options=[None],value=None,description='x-axis:',disabled=True)
        self.widget_Yobs = widgets.Dropdown(options=[None],value=None,description='y-axis:',disabled=True)
        self.widget_replot_observable = widgets.Button(description="Plot",disabled=True)
        self.notebook.dependencies_dict["seed"].append(self)
    def replot_observable(self,b):
        plt.close()
        clear_output()
        self.notebook.sim.custom_plot(self.notebook.seed,self.widget_Xobs.value,[self.widget_Yobs.value])


    def display(self):
        self.widget_replot_observable.on_click(self.replot_observable)
        plot_observable_options = widgets.VBox([widgets.HBox([self.widget_Xobs,self.widget_Yobs]),widgets.HBox([self.widget_replot_observable])])
        display(plot_observable_options)

    def update(self):
        if self.notebook.seed is None:
            self.widget_Xobs.disabled = self.widget_Yobs.disabled = self.widget_replot_observable.disabled = True
            self.widget_Xobs.value = self.widget_Yobs.value = None
        else:
            self.widget_Xobs.disabled = self.widget_Yobs.disabled = self.widget_replot_observable.disabled = False
            self.widget_Xobs.options = list(self.notebook.sim.seeds[self.notebook.seed].observables.keys())
            self.widget_Yobs.options = list(self.widget_Xobs.options)
            self.widget_Xobs.value = "generation"
            self.widget_Yobs.value = self.notebook.sim.seeds[self.notebook.seed].default_observable


class Select_Generation(CellModule):
    def __init__(self,Notebook):
        super(Select_Generation, self).__init__(Notebook)
        self.notebook.dependencies_dict["seed"].append(self)
        self.widget_gen = widgets.IntSlider(value = 0,min=0,max=0,description = 'Gen:',disabled=True)
        self.widget_restart_gen = widgets.IntSlider(value = 0,min=0,max=0,description = 'Restart Gen:',disabled=True)
        self.widget_restart_net = widgets.IntSlider(value = 0,min=0,max=0,description = 'Network:',disabled=True)
    def read_best_generation(self,gen_index):
        if not self.widget_gen.disabled:
            self.notebook.generation = gen_index
            self.notebook.net = self.notebook.sim.get_best_net(self.notebook.seed,self.notebook.generation)
            for cell in self.notebook.dependencies_dict["generation"]:
                cell.update()
    def read_restart_generation(self,gen_index,net_index):
        if not self.widget_restart_gen.disabled:
            self.notebook.generation = gen_index
            self.notebook.net = self.notebook.sim.get_backup_net(self.notebook.seed,gen_index,net_index)
            for cell in self.notebook.dependencies_dict["generation"]:
                cell.update()
    def display(self):

        widget1 = widgets.VBox([widgets.HTML("Select Best Network"),interactive(self.read_best_generation,gen_index=self.widget_gen)])
        widget2 = widgets.VBox([widgets.HTML("Select Backup Network"),interactive(self.read_restart_generation,gen_index=self.widget_restart_gen,net_index=self.widget_restart_net)])


        to_display = widgets.HBox([widget1,widget2])
        display(to_display)

    def update(self):
        if self.notebook.seed is None:
            self.widget_gen.value = 0
            self.widget_gen.min = self.widget_gen.max = 0
            self.widget_gen.disabled = True
            self.widget_restart_gen.disabled = True
            self.widget_restart_net.disabled = True
            self.notebook.generation = None
            self.notebook.net = None
        else:
            self.notebook.generation = None
            self.notebook.net = None
            self.widget_gen.value = 0
            self.widget_gen.disabled = False
            self.widget_gen.max = len(self.notebook.sim.seeds[self.notebook.seed].generations)-1
            self.widget_restart_gen.disabled = False
            restart_generations = list(self.notebook.sim.seeds[self.notebook.seed].restart_generations)
            step = restart_generations[1]
            self.widget_restart_gen.max = restart_generations[-1]
            self.widget_restart_gen.step = step
            self.widget_restart_net.disabled = False
            self.widget_restart_net.max = self.notebook.sim.seeds[self.notebook.seed].pop_size -1

class Plot_Layout(CellModule):
    def __init__(self,Notebook):
        super(Plot_Layout, self).__init__(Notebook)
        self.notebook.dependencies_dict["seed"].append(self)
        self.notebook.dependencies_dict["generation"].append(self)
        self.button_plotLayout = widgets.Button(description="Plot network layout",disabled=True)

    def plot_layout(self,button):
        plt.close()
        clear_output()
        self.notebook.net.draw()

    def update(self):
        if self.notebook.generation is None:
            self.button_plotLayout.disabled = True
        else:
            self.button_plotLayout.disabled = False

    def display(self):
        self.button_plotLayout.on_click(self.plot_layout)
        display(self.button_plotLayout)

class Run_Dynamics(CellModule):
    def __init__(self,Notebook):
        super(Run_Dynamics, self).__init__(Notebook)
        self.widget_nputs = widgets.IntText(value=1,description='N :',disabled=False)
        self.button_launchRun = widgets.Button(description="Run dynamics",disabled=True)
        self.notebook.dependencies_dict["generation"].append(self)
        self.notebook.dependencies_dict["dynamics"] = []
        self.notebook.extra_variables = {"ntries":None}
    def launch_dynamics(self,button):
        self.notebook.sim.run_dynamics(net=self.notebook.net,erase_buffer=False,trial=self.widget_nputs.value)
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

class Plot_Dynamics(CellModule):
    def __init__(self,Notebook):
        super(Plot_Dynamics, self).__init__(Notebook)
        self.notebook.dependencies_dict["dynamics"].append(self)
        self.notebook.dependencies_dict["generation"].append(self)
        self.widget_selectInput = widgets.IntSlider(value = 0,min=0,max=0,description = 'Input:',disabled=True)
        self.widget_selectCell = widgets.IntSlider(value = 0,min=0,max=0,description = 'Cell:',disabled=True)
        self.button_plotdynamics = widgets.Button(description="Plot dynamics",disabled=True)

    def plot_dynamics(self,button):
        plt.close()
        clear_output()
        self.notebook.sim.Plot_TimeCourse(self.widget_selectInput.value,cell=self.widget_selectCell.value)

    def update(self):
        if self.notebook.extra_variables.get("ntries",None) is None or self.notebook.generation is None:
            self.widget_selectInput.value=self.widget_selectInput.min=self.widget_selectInput.max = 0
            self.widget_selectCell.value=self.widget_selectCell.min=self.widget_selectCell.max = 0
            self.button_plotdynamics.disabled = True
        else:
            self.widget_selectInput.max = self.notebook.extra_variables["ntries"]-1
            self.widget_selectCell.max = self.notebook.sim.inits.prmt["ncelltot"]-1
            self.widget_selectInput.disabled = self.widget_selectCell.disabled = False
            self.button_plotdynamics.disabled = False

    def display(self):
        self.button_plotdynamics.on_click(self.plot_dynamics)
        display(widgets.HBox([self.widget_selectInput,self.widget_selectCell,self.button_plotdynamics]))

class Plot_Cell_Profile(CellModule):
    def __init__(self,Notebook):
        super(Plot_Cell_Profile, self).__init__(Notebook)
        self.notebook.dependencies_dict["dynamics"].append(self)
        self.notebook.dependencies_dict["generation"].append(self)
        self.widget_selectInput = widgets.IntSlider(value = 0,min=0,max=0,description = 'Input:',disabled=True)
        self.widget_selectTime = widgets.IntSlider(value = 0,min=0,max=0,description = 'Time:',disabled=True)
        self.button_plotdynamics = widgets.Button(description="Plot profile",disabled=True)

    def plot_dynamics(self,button):
        plt.close()
        clear_output()
        self.notebook.sim.Plot_Profile(trial_index=self.widget_selectInput.value,time=self.widget_selectTime.value)

    def update(self):
        if self.notebook.extra_variables.get("ntries",None) is None or self.notebook.generation is None:
            self.widget_selectInput.value=self.widget_selectInput.min=self.widget_selectInput.max = 0
            self.widget_selectTime.value=self.widget_selectTime.min=self.widget_selectTime.max = 0
            self.button_plotdynamics.disabled = True
        else:
            self.widget_selectInput.max = self.notebook.extra_variables["ntries"]-1
            self.widget_selectTime.max = self.notebook.sim.inits.prmt["nstep"]-1
            self.widget_selectInput.disabled = self.widget_selectTime.disabled = False
            self.button_plotdynamics.disabled = False

    def display(self):
        self.button_plotdynamics.on_click(self.plot_dynamics)
        display(widgets.HBox([self.widget_selectInput,self.widget_selectTime,self.button_plotdynamics]))


class Plot_Pareto_Fronts(CellModule):
    def __init__(self,Notebook):
        super(Plot_Pareto_Fronts, self).__init__(Notebook)
        self.notebook.dependencies_dict["seed"].append(self)
        self.widget_selectGenerations = widgets.SelectMultiple(options=[None],value=[None],description='Generation',disabled=True)
        self.widget_selectText = widgets.Text(value="",placeholder="List of generations separated with commas.",disabled=True)
        self.widget_plot = widgets.Button(description="Plot Pareto Fronts",disabled=True)
        self._widget_with_indexes  = widgets.Checkbox(value=False,description='Display indexes',disabled=False)


    def plot_function(self,button):
        plt.close()
        clear_output()
        gen = self.widget_selectGenerations.value
        if self.widget_selectText.value:
            gen = [int(xx) for xx in self.widget_selectText.value.split(",")]
            #gen = [int(xx) for xx in self.widget_selectText.value.split(",")]

        self.notebook.sim.seeds[self.notebook.seed].plot_pareto_fronts(gen,True)#self._widget_with_indexes.value)
    def update(self):
        if self.notebook.seed is None or self.notebook.type!="pareto":
            self.widget_selectGenerations.options = [None]
            self.widget_selectGenerations.value = [None]
            self.widget_selectGenerations.disabled = True
            self.widget_selectText.value = ""
            self.widget_selectText.disabled = True
            self.widget_plot.disabled = True
            self._widget_with_indexes.value = False
        else:
            self._widget_with_indexes.value = False
            self.widget_selectGenerations.disabled = False
            self.widget_plot.disabled = False
            self.widget_selectGenerations.options = self.notebook.sim.seeds[self.notebook.seed].restart_generations
            self.widget_selectGenerations.value = []
            self.widget_selectText.value = ""
            self.widget_selectText.disabled = False

    def display(self):

        #interactive(self.read_selected,generations=self.widget_selectGenerations)
        self.widget_plot.on_click(self.plot_function)
        instructions  = widgets.HTML("<p>Press <i>ctrl</i>, <i>cmd</i>, or <i>shift</i>  for multi-select</p>")
        to_display = widgets.VBox([instructions,widgets.HBox([self.widget_selectGenerations,self.widget_selectText]),self.widget_plot])
        #to_display = widgets.VBox([self.widget_plot])
        display(to_display)
