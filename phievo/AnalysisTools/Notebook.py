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
        self.load_project = Load_Project(self)
        self.select_seed = Select_Seed(self)
        self.project = None
        self.seed = None
        self.generation = None

class Load_Project(object):
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
            self.notebook.select_seed.update()
            self.notebook.generation = None
            self.widget_loadDir.button_style = "success"
            self.widget_loadDir.disabled=True

            #seed = {"Seed%d"%i:i for i,seed in self.notebook.seeds.items()}

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

    def read_seed(self,seed_name):
        self.notebook.seed =  seed_name

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
        self.widget_Xobs = widgets.Dropdown(options=choices,value="generation",description='x-axis:',disabled=True)
        self.widget_Yobs = widgets.Dropdown(options=choices,value="fitness",description='y-axis:',disabled=True)
        self.widget_replot_observable = widgets.Button(description="Plot",disabled=True)

    def replot_fitness(self,b):
        plt.close()
        clear_output()
        self.notebook.sim.show_fitness(seed)


    def display(self):
        widget_replot_observable.on_click(self.replot_fitness)
        plot_observable_options = widgets.VBox([widgets.HBox([widget_Xobs,widget_Yobs]),widgets.HBox([widget_replot_observable])])
        display(plot_observable_options)

    def update(self):
        NotImplemented
