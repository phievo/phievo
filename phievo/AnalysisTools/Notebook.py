import matplotlib.pyplot as plt
import numpy as np
from  ipywidgets import widgets
from ipywidgets import interact, interactive, fixed
from IPython.display import display,HTML,clear_output
import os
from phievo.AnalysisTools.Simulation import Simulation

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
        self.seed =None
        self.load_project = Load_Project(self)

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
            ## Search for the Simulation directory
            if os.path.isdir(path):
                # self.foundSeed_widget.value = ""
                # foundDirFlag = True
                self.foundFile_widget.value = found_str
                self.widget_loadDir.disabled=False
            else:
                #self.foundDirFlag = False
                self.foundFile_widget.value = notfound_str
                self.widget_loadDir.disabled=True


        def display(self):
            # widget_loadDir.on_click(onLoadRunClicked)
            interactive(self.inspect_run,path=self.widget_select_project);
            main_options = widgets.VBox([widgets.HBox([self.widget_select_project,self.foundFile_widget]),self.widget_loadDir])
            display(main_options)
            #display(self.widget_select_project)
