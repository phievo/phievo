import matplotlib.pyplot as plt
import numpy as np
from  ipywidgets import widgets
from ipywidgets import interact, interactive, fixed
from IPython.display import display,HTML,clear_output
import os
from phievo.AnalysisTools.Simulation import Simulation


class Notebook:
    """
        Wrapper that contains both the the widgets and  the simulation results.
        This way it is easy to update the state of the widgets when you load a
        new simulation
    """
    def __init__(self):
        self.sim = None
        self.frame

class __frame__:
    NotImplemented
