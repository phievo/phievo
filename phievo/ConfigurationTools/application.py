import ipywidgets as w
from IPython.display import display

import phievo.ConfigurationTools.widgetfunctions as wf
import phievo.ConfigurationTools.containers as wc

from phievo.Networks import mutation
from phievo.Networks import initialization
from phievo.initialization_code import ccode_dir
import os
from shutil import copyfile


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
    if type(val) is str:
        to_return = "{}[\"{}\"] = \"{}\"".format(dict_name,key,val)
    else:
        to_return = "{}[\"{}\"] = {}".format(dict_name,key,val)
    return to_return

tab_labels = {
    "dictionary_mutation" : "Mutation parameters",
    "dictionary_ranges" : "Kinetic parameters",
    "restart":"Restart",
    "prmt":"General simulation parameters",
    "codes":"Initializations codes",
    "cfile":"C files"
}

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
            "dictionary_mutation" : wc.w_table("dictionary_mutation",configurations["dictionary_mutation"],"float_range_widget","<h2>Mutation parameters</h2>"),
            "dictionary_ranges" : wc.w_table("dictionary_ranges",configurations["dictionary_ranges"],"float_range_widget","<h2>Kinetic parameters</h2>"),
            "restart" : wf.w_restart(infos="<h2>Restart</h2><p>For now restart needs to be defined manually in the initialization file of an existing project.</p><p>This tab can be used to set the frequence at which a complete generation is saved.</p>"),
            "prmt":wc.prmt_widget(values=configurations["prmt"],infos="<h2>General simulation parameters</h2>"),
            "codes":wc.widget_initialization("codes","<h2>Initializations codes</h2>"),
            "cfile":wc.ccode_widget("cfile",configurations["cfile"],infos="<h2>cfiles</h2>\n<p>You may leave these setting as default. In this case the blank files will be created in the project directory and can be updated before starting a simulation. The files that already have a setting should be modified only by advanced users.</p><p><b>Note:</b> It is important that you update the fitness.c since its default value returns 1 for all the networks.</p>")
        }
        self.w_tab = w.Tab()
        self.w_tab.children = [self.tabs[key].get_widget() for key in self.tabs.keys()]
        for i,name in enumerate(list(self.tabs.keys())):
            self.w_tab.set_title(i,tab_labels[name])

        self.create_button = w.Button(description="Write project",button_style="info")
        self.create_button.on_click(self.write)
    def get_widget(self):
        pname_box = w.HBox([self.project_name,self.project_exists])
        return w.VBox([pname_box,self.w_tab,self.create_button])

    def get_values(self):
        return {key:self.tabs[key].get_values() for key in self.tabs.keys()}

    def write(self,button):
        data = self.get_values()
        to_write = []
        proj_dir = self.project_name.value
        os.makedirs(proj_dir)
        to_write += ["## Initialization file for project {}.".format(self.project_name.value)]
        to_write += ["\n## Kinetic parameters\ndictionary_ranges = {}"]
        for key,val in data["dictionary_ranges"].items():
            to_write.append(to_dict("dictionary_ranges",key,val))
        to_write += ["\n## Mutation parameters\ndictionary_mutation = {}"]
        for key,val in data["dictionary_mutation"].items():
            to_write.append(to_dict("dictionary_mutation",key,val))

        dir_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),"default_cfiles")
        print(dir_path)
        to_write += ["\n## Cfiles\ncfile = {}"]
        for key,val in data["cfile"].items():            
            if val=="":
                val = key+".c"
                copyfile(os.path.join(dir_path,val),os.path.join(proj_dir,val))
            elif val==configurations["cfile"][key]:
                continue
            else:
                old_path = val
                val = os.path.split(val)[-1]
                copyfile(old_path,os.path.join(proj_dir,val))
            to_write.append(to_dict("cfile",key,val))
            print("Wrote {}".format(os.path.join(proj_dir,val)))

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
        with open(os.path.join(proj_dir,"initialization.py"),"w") as init_file:
            init_file.write("\n".join(to_write))
            print("Wrote {}".format(os.path.join(proj_dir,"initialization.py")))

        
