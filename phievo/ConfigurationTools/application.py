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

def to_dict(dict_name,key,val):
    return "{}[\"{}\"] = {}".format(dict_name,key,val)

class App:
    def __init__(self):
        self.project_name = w.Text(description="Project",value="",placeholder="Enter the project name")
        self.project_exists = w.Label("A project cannot be an empty string or an existing directory.")
        w.Valid(value=False,description=" ")
        def update_validity(path):
            if os.path.isdir(path) or path=="":
                self.project_exists.value = "A project cannot be an empty string or an existing directory."
            else:
                self.project_exists.value = "Valid new project name."
        w.interactive(update_validity,path=self.project_name)
        self.tabs = {
            "dictionary_mutation" : wc.w_table("dictionary_mutation",configurations["dictionary_mutation"],"float_range_widget","<p>Contains the parameter related to mutations.</p>"),
            "dictionary_ranges" : wc.w_table("dictionary_ranges",configurations["dictionary_ranges"],"float_range_widget","<p>Contains the kinetic parameter ranges.</p>"),
            "restart" : wf.w_restart(),
            "prmt":wc.prmt_widget(values=configurations["prmt"]),
            "codes":wc.widget_initialization("codes","<h2>Initializations codes</h2>"),
            "cfile":wc.ccode_widget("cfile",configurations["cfile"],infos="<h2>cfiles</h2>\n<p>You may leave these setting as default. In this case the blank files will be created in the project directory and can be updated before starting a simulation. The files that already have a setting should be modified only by advanced users.</p><p><b>Note:</b> It is important that you update the fitness.c since its default value returns 1 for all the networks.</p>")
        }
        self.w_tab = w.Tab()
        self.w_tab.children = [self.tabs[key].get_widget() for key in self.tabs.keys()]
        for i,name in enumerate(list(self.tabs.keys())):
            self.w_tab.set_title(i,name)

        self.create_button = w.Button(description="Write project",button_style="info")
        self.create_button.on_click(self.write)
    def get_widget(self):
        pname_box = w.HBox([self.project_name,self.project_exists])
        return w.VBox([pname_box,self.w_tab,self.create_button])

    def get_values(self):
        return {key:self.tabs[key].get_values() for key in self.tabs.keys()}

    def write(self):
        data = self.get_values()
        to_write = []

        to_write += ["## Initialization file for project {}.".format(self.project_name.value)]
        to_write += ["\n## Kinetic parameters\ndictionary_ranges = {}"]
        for key,val in data["dictionary_ranges"].items():
            to_write.append(to_dict("dictionary_ranges",key,val))
        to_write += ["\n## Mutation parameters\ndictionary_mutation = {}"]
        for key,val in data["dictionary_mutation"].items():
            to_write.append(to_dict("dictionary_mutation",key,val))
        to_write += ["\n## Cfiles\n cfile = {}"]
        for key,val in data["cfile"].items():
            if val=="":
                val = key+".c"
            to_write.append(to_dict("cfile",key,val))

        tags = dict(list_types_output=data["prmt"].pop("list_types_output"),list_unremovable=data["prmt"].pop("list_unremovable"))
        to_write += ["\n## General simulation parameters\nprmt = {}"]
        for key,val in data["prmt"].items():
            to_write.append(to_dict("prmt",key,val))
            
        to_write += ["\n## General simulation parameters\nprmt[\"restart\"] = {}"]        
        for key,val in data["restart"].items():
            to_write.append("prmt[\"restart\"]"+to_dict("",key,val))

        to_write += ["\n## Outputs and unremovables"]
        for key,val in tags.items():
            to_write.append("{} = {}".format(key,val))
            
        to_write += ["\n## Initialize networks"]
        to_write.append(data["codes"]["init_network"])
        to_write += ["\n## Fitness treatment function"]
        to_write.append(data["codes"]["fitness_treatment"])
        return to_write
        
