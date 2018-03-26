import phievo.ConfigurationTools.widgetfunctions as wf
import ipywidgets as w

class w_table:
    
    def __init__(self,name,parameters,w_type,infos=""):
        self.name = name
        self.infos = infos
        if type(w_type) is str:
            w_type = {key:w_type for key in parameters}
        self.obj_dict = {key:getattr(wf,w_type[key])(key,val) for key,val in parameters.items()}
    def get_widget(self):
        infos = w.HTML(value="{}".format(self.infos))
        return w.VBox([infos]+[self.obj_dict[key].get_widget() for key in self.obj_dict.keys()])
    def get_values(self):
        return {key:self.obj_dict[key].get_values() for key in self.obj_dict.keys()}

    def set_values(self,values):
        for key,val in values.items():
            self.obj_dict[key].set_values(val)

prmt_order = ["nseed","firstseed","ngeneration","ncelltot","npopulation","nneighbor","frac_mutate","ninput","noutput","ntries","dt","nstep","langevin_noise","tgeneration","redo","pareto","npareto_functions","rshare","multipro_level","freq_stat"]

class drange_woidget:
    def __init__(self)
    
prmt_types ={
    "nseed":"int_widget",
    "firstseed":"int_widget",
    "ngeneration":"int_widget",
    "ncelltot":"int_widget",
    "npopulation":"int_widget",
    "nneighbor":"int_widget",
    "frac_mutate":"float_widget",
    "ninput":"int_widget",
    "noutput":"int_widget",
    "ntries":"int_widget",
    "dt":"float_widget",
    "nstep":"int_widget",
    "langevin_noise":"float_widget",
    "tgeneration":"float_widget",
    "redo":"bool_widget",
    "pareto":"bool_widget",
    "npareto_functions":"int_widget",
    "rshare":"float_widget",
    "multipro_level":"bool_widget",
    "freq_stat":"int_widget",
}

prmt_descriptions ={
    "nseed":"Number of seeds",
    "firstseed":"First seed",
    "ngeneration":"Number of generations",
    "ncelltot":"Number of cells",
    "npopulation":"Population size",
    "nneighbor":"Number of neighbors",
    "frac_mutate":"Fraction mutated per gen",
    "ninput":"Number of Inputs",
    "noutput":"Number of Outputs",
    "ntries":"Number of trials",
    "dt":"Time step dt",
    "nstep":"Number of time steps",
    "langevin_noise":"Langevin noise value",
    "tgeneration":"Gillespie generation time",
    "redo":"Recompute networks",
    "pareto":"Pareto simulation",
    "npareto_functions":"Number of pareto functions",
    "rshare":"Pareto penalty radius",
    "multipro_level":"Multiple threads",
    "freq_stat":"Generation printing frequency",
}
tag_choices = ["Species","Degradable","TF","Kinase","Phosphatase","Output","Input","Complexable","Complex","Phosphorylable","Diffusible"]

class prmt_widget:
    def __init__(self,values = None,infos=""):
        self.infos = w.HTML(value=infos)
        self.obj_dict = {key:getattr(wf,prmt_types[key])(prmt_descriptions[key],1) for key in prmt_order}
        if values: 
            self.set_values(values)
        def activate_pareto(action):
            for key in ["npareto_functions","rshare"]:
                self.obj_dict[key].value.disabled = not self.obj_dict["pareto"].get_values()
        w.interactive(activate_pareto,action=self.obj_dict["pareto"].value)
        activate_pareto(1)
        self.obj_dict["list_unremovable"] = wf.tags_widgets("Unremovable species",tag_choices,["Input","Output"])
        tag_choices.remove("Output")
        tag_choices.remove("Input")
        self.obj_dict["list_types_output"] = wf.tags_widgets("Possible outputs",tag_choices,["TF"])
        
    def get_widget(self):
        return w.VBox([self.infos]+[self.obj_dict[key].get_widget() for key in prmt_order+["list_unremovable","list_types_output"]])
    
    def set_values(self,values):
        for key,val in values.items():
            self.obj_dict[key].set_values(val)
            
    def get_values(self):
        return {key:self.obj_dict[key].get_values() for key in self.obj_dict.keys()}

default_initialization_code = "import random\n\
from phievo.Networks import mutation\n\
\n\
def init_network():\n\
   seed = int(random.random()*100000)\n\
   g = random.Random(seed)\n\
   net = mutation.Mutable_Network(g)\n\
   for i_ind in range(prmt[\"ninput\"]):\n\
       parameters=[['Degradable', mutation.sample_dictionary_ranges('Species.degradation',random)],['TF',1],['Input',i_ind]]\n\
       TF=net.new_Species(parameters)\n\
   for o_ind in range(prmt[\"noutput\"]):\n\
       [tm, prom, o1] = net.random_gene('TF')\n\
       o1.add_type(['Output',o_ind])\n\
   net.activator_required=1\n\
   net.fixed_activity_for_TF=0\n\
   net.write_id()\n\
   return net\n\
"
default_fitness_treatment = "def fitness_treatment(population):\n\
    \"\"\"\n\
        This function can be used to artificially transform the\n\
        the True fitness of the networks in the population. \n\
    \"\"\"\n\
    # Uncomment the next two lines to add radomness to the fitness: \n\
\n\
    # for nnetwork in range(population.npopulation):\n\
    #     population.genus[nnetwork].fitness += 0.001*random.random()"

class widget_initialization:
    def __init__(self,name,infos=""):
        self.infos = w.HTML(value=infos)
        self.init_network = wf.code_widget("init_code",default_initialization_code,"Network initialization code")
        self.fitness_treatment = wf.code_widget("fitness_treatment",default_fitness_treatment,"Fitness treatment",height="200px")
    def get_widget(self):
        return w.VBox([self.infos,self.init_network.get_widget(),self.fitness_treatment.get_widget()])
    def get_values(self):
        return dict(init_network=self.init_network.get_values(),fitness_treatment=self.fitness_treatment.get_values())

class ccode_widget:
    def __init__(self,name,values,infos = ""):
        self.name = name
        self.infos = w.HTML(value=infos)
        self.obj_dict = {key:wf.code_path_widget(key,val) for key,val in values.items()}
    def get_widget(self):
        return w.VBox([self.infos]+[self.obj_dict[key].get_widget() for key in self.obj_dict.keys()])
    def get_values(self):
        return {key:self.obj_dict[key].get_values() for key in self.obj_dict.keys()}
    def set_values(self,values):
        for key,val in values.items():
            self.obj_dict[key].set_values(val)
