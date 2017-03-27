import matplotlib.pyplot as plt
import numpy as np
from  ipywidgets import widgets
from ipywidgets import interact, interactive, fixed
from IPython.display import display,HTML,clear_output
import os
from phievo.AnalysisTools import Simulation

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

        self.project = None
        self.depends_on_project = [] ## List of cell objects to update when project changes
        self.seed = None
        self.depends_on_seed = [] ## List of cell objects to update when seed changes
        self.generation = None
        self.depends_on_generation = [] ## List of cell objects to update when generation changes

        self.select_project = Select_Project(self)
        self.select_seed = Select_Seed(self)
        self.plot_evolution_observable = Plot_Evolution_Observable(self)
        self.select_generation = Select_Generation(self)
        self.plot_layout = Plot_Layout(self)


class Select_Project(object):
        def __init__(self,Notebook):
            self.widget_select_project = widgets.Text(value='',placeholder='Name of project directory',description='Directory:',disabled=False)
            self.widget_loadDir = widgets.Button(description="Load Run",disabled=True)
            self.foundFile_widget = widgets.HTML("")
            self.notebook = Notebook
        def inspect_run(self,path):
            """
            Test if the dir name exists

            Args:
                path (str): path of the directory
            """
            self.update()
            if os.path.isdir(path):
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
            for cell in self.notebook.depends_on_project:
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

class Select_Seed:
    def __init__(self,Notebook):
        self.notebook = Notebook
        self.widget_select_seed = widgets.Dropdown(options={"None":None},value=None,description='Seed:',disabled=True)
        self.notebook.depends_on_project.append(self)

    def read_seed(self,seed_name):
        self.notebook.seed =  seed_name
        for cell in self.notebook.depends_on_seed:
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
            seeds = self.widget_select_seed.options = {"Seed%d"%i:i for i,seed in self.notebook.sim.seeds.items()}
            self.widget_select_seed.disabled = False
            self.widget_select_seed.options = {"Seed%d"%i:i for i,seed in self.notebook.sim.seeds.items()}
            self.widget_select_seed.value = self.widget_select_seed.options[list(self.widget_select_seed.options.keys())[0]]
            self.notebook.seed = self.widget_select_seed.value

class Plot_Evolution_Observable:
    def __init__(self,Notebook):
        self.notebook = Notebook
        self.widget_Xobs = widgets.Dropdown(options=[None],value=None,description='x-axis:',disabled=True)
        self.widget_Yobs = widgets.Dropdown(options=[None],value=None,description='y-axis:',disabled=True)
        self.widget_replot_observable = widgets.Button(description="Plot",disabled=True)
        self.notebook.depends_on_seed.append(self)
    def replot_observable(self,b):
        plt.close()
        clear_output()
        self.notebook.sim.custom_plot(self.notebook.seed,self.widget_Xobs.value,self.widget_Yobs.value)


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
            self.widget_Xobs.options = list(self.notebook.sim.seeds[self.notebook.seed].observables)
            self.widget_Yobs.options = list(self.widget_Xobs.options)
            self.widget_Xobs.value = "generation"
            self.widget_Yobs.value = "fitness"


class Select_Generation:
    def __init__(self,Notebook):
        self.notebook = Notebook
        self.notebook.depends_on_seed.append(self)
        self.widget_gen = widgets.IntSlider(value = 0,min=0,max=0,description = 'Gen:',disabled=True)

    def read_generation(self,gen_index):
        self.notebook.generation = gen_index
        for cell in self.notebook.depends_on_generation:
            cell.update()

    def display(self):
        interactive(self.read_generation,gen_index=self.widget_gen)
        display(self.widget_gen)

    def update(self):
        if self.notebook.seed is None:
            self.widget_gen.value = 0
            self.widget_gen.min = self.widget_gen.max = 0
            self.widget_gen.min.disabled = True
        else:
            self.widget_gen.value = 0
            self.widget_gen.disabled = False
            self.widget_gen.max = len(self.notebook.sim.seeds[self.notebook.seed].generations)-1

class Plot_Layout:
    def __init__(self,Notebook):
        self.notebook = Notebook
        self.notebook.depends_on_generation.append(self)
        self.button_plotLayout = widgets.Button(description="Plot network layout",disabled=True)

    def plot_layout(self,button):
        net = self.notebook.sim.get_best_net(self.notebook.seed,self.notebook.generation)
        plt.close()
        clear_output()
        net.draw()

    def update(self):
        if self.notebook.generation is None:
            self.button_plotLayout.disabled = True
        else:
            self.button_plotLayout.disabled = False

    def display(self):
        self.button_plotLayout.on_click(self.plot_layout)
        display(self.button_plotLayout)
