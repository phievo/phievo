import ipywidgets as w
from IPython.display import display

import phievo.ConfigurationTools.widgetfunctions as wf
import phievo.ConfigurationTools.containers as wc

from phievo.Networks import mutation
from phievo.Networks import initialization
from phievo.initialization_code import ccode_dir
import os

PRINT = False
pfile = {"deriv2" : "phievo.Networks.deriv2",
             "interaction" : "phievo.Networks.interaction",
             "pretty_graph": "phievo.Networks.lovelyGraph",
             "plotdata" : "phievo.Networks.plotdata"}

cfile = {
    "fitness":"",
    "init_history":"",
    "input":"",
    "header" : os.path.join(ccode_dir,'integrator_header.h'),
    "utilities" : os.path.join(ccode_dir,'utilities.c'),
    "geometry" : os.path.join(ccode_dir,'linear_geometry.c'),
    "integrator" : os.path.join(ccode_dir,'euler_integrator.c'),
    "main" : os.path.join(ccode_dir,'main_general.c')
}

prmt = initialization.prmt
restart = prmt.pop("restart")

configurations = {
    "dictionary_ranges":mutation.dictionary_ranges,
    "dictionary_mutation":initialization.dictionary_mutation,
    "cfile":cfile,
    "pfile":pfile,
    "prmt":prmt,
    "restart":restart
}

class App:
    def __init__(self):
        self.tabs = {
            "dictionary_mutation" : wc.w_table("dictionary_mutation",configurations["dictionary_mutation"],"float_range_widget","<p>Contains the parameter related to mutations.</p>"),
            "dictionary_ranges" : wc.w_table("dictionary_ranges",configurations["dictionary_ranges"],"float_range_widget","<p>Contains the kinetic parameter ranges.</p>"),
            "restart" : wf.w_restart()
        }
        self.w_tab = w.Tab()
        self.w_tab.children = [self.tabs[key].get_widget() for key in self.tabs.keys()]
        for i,name in enumerate(list(self.tabs.keys())):
            self.w_tab.set_title(i,name)
    def get_widget(self):
        return self.w_tab
