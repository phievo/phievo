import ipywidgets as w
from IPython.display import display
import random
from phievo.Networks import mutation
import re 
tag_choices = {
    "Species":(True,),
    "Degradable":(True,"rate","Float",0.5,"Species.degradation"),
    "TF":(True,"activity","Int",1),
    "Kinase":(False,),
    "Phosphatase":(False,),
    "Output":(False,"track index","Int"),
    "Input":(False,"track index","Int"),
    "Complexable":(False,),
    "Complex":(False,),
    "Phosphorylable":(False,),
    "Diffusible":(False,"rate","Float",0,"Species.diffusion")
}
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
def tags_to_str(tags):
    array = ["[\"{}\",{}]".format(*tt) if len(tt)==2 else "[\"{}\"]".format(*tt) for tt in tags]
    
    return "["+",".join(array)+"]"
class add_tag_widget:
    def __init__(self,tag,state=False,param=None,ptype=None,defaultvalue=0,draw_random=""):
        self.tag = tag
        self.tag_w = w.HTML(description=" ",value="<p><b>{}</b></p>".format(tag),layout=w.Layout(width="15%", height="40px"))
        self.is_on = w.Checkbox(value=bool(state),description=" ",layout=w.Layout(width="15%", height="40px"))
        self.is_on_d = self.is_on.value
        
        if param:
            self.custom_p = w.Checkbox(description="custom parameters",value=False,disabled=not state)
            self.param = getattr(w,"Bounded{}Text".format(ptype))(value=defaultvalue,description=param,disabled=not state or not self.custom_p.value,layout=w.Layout(width="20%", height="40px"))
            
            if tag=="TF":
                self.param = w.Dropdown(options={"repressor by default":0,"activator by default":1},value=1,disabled=not self.is_on.value)
                self.custom_p.value = True
            
            def activate(b):
                if self.tag!="TF":
                    self.custom_p.disabled = not self.is_on.value
                self.param.disabled = not self.is_on.value or not self.custom_p.value
                
            
            w.interactive(activate,b=self.is_on)
            w.interactive(activate,b=self.custom_p)
            self.param_d = self.param.value
        self.draw_random = 'mutation.sample_dictionary_ranges(\"{}\",net.Random)'.format(draw_random)
    
    def get_widget(self):
        inhbox = [self.is_on,self.tag_w]
        if hasattr(self,"param"):
            inhbox.append(self.param)
            if self.tag!="TF":
                inhbox.append(self.custom_p)
        return w.HBox(inhbox)
    def enable(self):
        if hasattr(self,"param"):
            self.is_on.disabled = False
            self.custom_p.disabled = not self.is_on.value
            self.param.disabled = not self.is_on.value or not self.custom_p.value
    def disable(self):
        if hasattr(self,"param"):
            self.is_on.value = False
            self.is_on.disabled = True
            self.custom_p.disabled = True
            self.param.disabled = True
    
            
    def get_tag(self):
        if self.is_on.value:
            toreturn = [self.tag]
            if hasattr(self,"param"):
                if self.custom_p.value:
                    pval = self.param.value
                else:
                    if self.tag=="Input":
                        pval = "{ind_i}"
                    elif self.tag=="Output":
                        pval = "{ind_o}"
                    else:
                        pval = self.draw_random
                toreturn.append(pval)
        else:
            toreturn = None
        return toreturn
        
    def clear(self):
        self.is_on.value = self.is_on_d
        if hasattr(self,"param"):
            #self.param.value = self.param_d
            if self.tag!="TF":
                self.custom_p.value=False
                self.param.disabled = True

class add_species_widget:
    def __init__(self):
        #self.select_ = 
        self.isGene = w.Checkbox(value=True,description="gene",layout=w.Layout(width="15%", height="40px"))
        self.tags = {key:add_tag_widget(key,*val) for key,val in tag_choices.items()}
        self.tags["Input"].is_on.disabled = True
        self.add_button = w.Button(description="Add",button_style="info")
        self.israndom = w.Checkbox(value=True,description=" ")

        self.output_counter = w.HTML("outputs: ()")
        self.input_counter = w.HTML("inputs: ()")
        self.g_basal= w.BoundedFloatText(description="basal:",value=0,step=0.001,layout=w.Layout(width="15%", height="40px"),disabled=True)
        self.g_rate = w.BoundedFloatText(description="rate:",value=1,step=0.001,layout=w.Layout(width="15%", height="40px"),disabled=True)
        self.g_delay = w.BoundedFloatText(description="delay:",value=0,step=0.001,layout=w.Layout(width="15%", height="40px"),disabled=True)
        self.g_custom = w.Checkbox(description="Custom parameters",value=False,disabled=False)
        
        self.message =w.HTML("")
        def input_not_gene1(val):
            self.isGene.disabled = self.tags["Input"].is_on.value
            
            self.tags["Input"].enable()            
            if self.tags["Input"].is_on.value:
                self.tags["Output"].disable()
            else:
                self.tags["Output"].enable()
        def input_not_gene2(val):
            self.g_basal.disabled = not self.isGene.value or not self.g_custom.value
            self.g_delay.disabled = not self.isGene.value or not self.g_custom.value
            self.g_rate.disabled = not self.isGene.value or not self.g_custom.value            
            self.g_custom.disabled = not self.isGene.value
            if self.isGene.value:
                self.tags["Input"].disable()
            else:
                self.tags["Input"].enable()
        def inp_not_out(val):
            self.tags["Output"].enable()            
            if self.tags["Output"].is_on.value:
                self.tags["Input"].disable()
            else:
                self.tags["Input"].enable()
        w.interactive(input_not_gene1,val=self.tags["Input"].is_on)
        w.interactive(input_not_gene2,val=self.isGene)
        w.interactive(input_not_gene2,val=self.g_custom)
        w.interactive(inp_not_out,val=self.tags["Output"].is_on)
        
    def get_widget(self):
        return w.VBox(
                      [self.input_counter]\
                      +[self.output_counter]\
                      +[w.HTML(value="<b>gene settings:</b><p>If \"gene\" is selected, the program also adds a TModule and a CorePromoter and the species is regulated.</p><p>The parameters that are not customized will be drawn randomly in the space of the possible values.</p><p>An input species cannot be a gene.</p>")]
                      +[w.HBox([self.isGene,self.g_rate,self.g_basal,self.g_delay,self.g_custom])]\
                      +[w.HTML(value="<b>add tags:</b><p>The parameters that are not customized will be drawn randomly in the space of the possible values or will take the next available index (for Inputs and Outputs).</p>")]\
                      +[self.tags[key].get_widget() for key in self.tags.keys()]\
                      +[w.HBox([self.add_button,self.message])])
    def clear(self):
        self.add_button.description = "Add"
        self.add_button.button_style = "info"
        self.isGene.value = True
        self.g_rate.disabled = True
        self.g_basal.disabled = True
        self.g_delay.disabled = True
        self.g_custom.value = False
        self.g_custom.disabled = False
        for key in self.tags.keys():
            self.tags[key].clear()
        self.message.value = ""
    def error(self,message):
        self.add_button.description = "Clear"
        self.add_button.button_style = "danger"
        self.message.value = "Error: "+message
    def get_command(self):
        tags = [tag for tag in [self.tags[key].get_tag() for key in self.tags.keys()] if tag]
        param = "\tparam = {}".format(tags_to_str(tags))
        if self.isGene.value:
            if self.g_custom.value:
                function = "\ttm_d[ind],prom_d[ind],s_d[ind] = net.new_gene({rate},{delay},param,{basal})".format(rate=self.g_rate.value,delay=self.g_delay.value,basal=self.g_basal.value)
            else:
                function = "\ttm_d[ind],prom_d[ind],s_d[ind] = net.new_custom_random_gene(param)"
        else:
            function = "\ts_d[ind] = net.new_Species(param)"
        return [param,function],tags
                

def empty_net():
    import random
    from phievo.Networks import mutation
    seed = int(random.random()*100000)
    g = random.Random(seed)
    net = mutation.Mutable_Network(g)
    return net
def convert_code_selfnet(code):
    if type(code) is list:
        code = ("\n".join(code))
    code = code.replace("\t","")
    for key in ["tm_d","prom_d","s_d","net"]:
        code = code.replace(key,"self.{}".format(key))
    return code


class init_network_widget:
    def __init__(self,infos=""):
        self.infos = w.HTML(value=infos)
        self.add_species_w = add_species_widget()
        self.nb_inputs = 1
        self.nb_outputs = 1
        self.table = w.HTML(value="")
        self.default_network = w.ToggleButton(description="Create a default network",value=False,layout=w.Layout(width="30%", height="40px"))
        self.fixed_activity_for_TF = w.Checkbox(description="Fixed Activity for TF",value=False)
        self.clear(1,clearall=False)
        self.add_inter_w = interaction_widget(self)
        
        self.add_species_w.add_button.on_click(self.add_species)
        self.add_inter_w.add_button.on_click(self.add_inter)
        self.clear_button = w.Button(description="Clear",button_style="danger")
        self.clear_button.on_click(self.clear)
        self.update_counters()
        self.infos = w.HTML("<h2>Create an initial network</h2><p>The network is used to initialize the first population before starting the simulation.</p><p>The initial network can either be created manually or left to default by activating the <code>Create a default network</code> button (In the latter case, the network generated with the tool is not taken into account). A default networks has only the appropriate number of inputs and outputs and no interactions.</p>")
        
    def clear(self,button,clearall=True):
        self.net = empty_net()
        self.code = []
        self.inputs = []
        self.outputs = []
        self.species = []
        self.types = {"Genes":[],"Complexable":[],"Kinase":[],"Phosphorylable":[],"TF":[]}
        self.tm_d = {}
        self.prom_d = {}
        self.s_d = {}
        self.species_list = []
        self.inter_list = []
        self.table.value = ""
        self.default_network.value = False
        self.fixed_activity_for_TF.value = False
        if clearall:
            self.add_species_w.clear()
            self.add_inter_w.clear()
        self.update_counters()
        
    def add_species(self,button):
        ind = len(self.species)
        
        if self.add_species_w.add_button.description=="Clear":
            self.add_species_w.clear()
            return None
        species_tags = []
        code,data = self.add_species_w.get_command()
        #print(code)
        species_tags = re.findall("\"\w+\"",code[0])
        species_tags = [tag.replace("\"","") for tag in species_tags]
        code[0] = code[0].format(ind_o=len(self.outputs),ind_i=len(self.inputs))
        #print(code[0])
        code[1] = code[1].replace("[ind]","[{}]".format(ind))
        #print(code[1])
        result = re.search("\"(Output|Input)\",(\d+)",code[0])
        #print(result)
        if not result:
            if (len(self.inputs)<self.nb_inputs or len(self.outputs)<self.nb_outputs):          
                self.add_species_w.error("All inputs and outputs must be added before adding other species.")
                return None
        else:
            # why do we need the following lines ?
            #trackind = int(result.group(2))
            #if trackind>(self.nb_outputs+self.nb_inputs):
            #    self.add_species_w.error("Track trackindex ({}) cannot be larger than the total number of inputs plus outputs.".format(trackind))
            #    return None
            #if trackind in self.inputs:
            #    self.add_species_w.error("Track trackindex {} is already an input.".format(trackind))
            #    return None
            #if trackind in self.outputs:
            #    self.add_species_w.error("Track trackindex {} is already an output.".format(trackind))
            #    return None
            
            if result.group(1)=="Input":
                if len(self.inputs)==self.nb_inputs:
                    self.add_species_w.error("There is already enough inputs in the network.")
                    return None
                self.inputs.append(ind)
                
            elif result.group(1)=="Output":
                if len(self.outputs)==self.nb_outputs:
                    self.add_species_w.error("There is already enough outputs in the network.")
                    return None
                self.outputs.append(ind)
                
        
        if "gene" in code[1]:
            self.types["Genes"].append(ind)
        for key in self.types.keys():
            if key in code[0]:
                self.types[key].append(ind)
        self.species.append(ind)
        self.species_list.append("<p>S{} ({})</p>".format(ind,",".join(species_tags)))
        self.code += code
        temp_code = convert_code_selfnet(self.code)
        
        exec(temp_code)
        self.add_inter_w.update_widgets()
        self.update_counters()
        self.add_species_w.clear()
        self.write_network()
        
    def write_network(self):
        self.table.value = "<b>Network:</b>\n"+"\n".join(self.species_list) + "\n" + "\n".join(self.inter_list)
    def add_inter(self,button):
        code = self.add_inter_w.get_command()
        self.code += [code]
        exec(convert_code_selfnet(code))
        self.write_network()
    def update_counters(self):
        self.add_species_w.input_counter.value = "inputs: ({}/{})".format(len(self.inputs),self.nb_inputs)
        self.add_species_w.output_counter.value = "outputs: ({}/{})".format(len(self.outputs),self.nb_outputs)
    def get_widget(self):
        self.accordion = w.Accordion(children=[self.add_species_w.get_widget(),self.add_inter_w.get_widget()])
        self.accordion.set_title(0, 'Add Species')
        self.accordion.set_title(1, 'Add Interaction')
        infos_manual = w.HTML("<p>Before building a network manually, we recommend that you read the <a href=\"http://phievo.readthedocs.io/en/latest/create_new_project.html#build-a-network-manually\">documentation</a>.</p>")
        fixed_activity_info = w.HTML("When a TF has a fixed activity, only the default type of the TF is considered (activator or repressor). The type of the TFHill does not matter.")
        return w.VBox([self.infos,self.default_network,infos_manual,self.accordion,fixed_activity_info,self.fixed_activity_for_TF,self.table,self.clear_button])

    def get_values(self):
        if self.default_network.value:
           return default_initialization_code
        else:
            code = "import random\nfrom phievo.Networks import mutation\n\ndef init_network():\n\tseed = int(random.random()*100000)\n\tg = random.Random(seed)\n\tnet = mutation.Mutable_Network(g)\n\ts_d,tm_d,prom_d = {},{},{}\n"
        code= code + "\n\tnet.activator_required=0\n\tnet.fixed_activity_for_TF={}\n".format(self.fixed_activity_for_TF.value)
        return code + "\n".join(self.code)+"\n\treturn net\n"    
    
class interaction_widget:
    def __init__(self,net):
        self.net = net
        self.interactions = {
            "TFHill":TFHill_widget(net,self),
            "Phosphorylation":Phosphorylation_widget(net,self),
            "PPI":PPI_widget(net,self)
        }
        self.add_button = w.Button(description="Add",button_style="info")
    def get_widget(self):
        to_return = []
        to_return += [self.interactions["TFHill"].get_widget()]
        to_return += [w.HTML("<hr>")]
        to_return += [self.interactions["Phosphorylation"].get_widget()]
        to_return += [w.HTML("<hr>")]
        to_return += [self.interactions["PPI"].get_widget()]
        to_return += [self.add_button]
        return w.VBox(to_return)
    def update_widgets(self):
        for key in self.interactions.keys():
            
            self.interactions[key].update_inter_state(False)
    def get_command(self):
        for key in self.interactions.keys():
            if self.interactions[key].is_on.value:
                command = self.interactions[key].get_command()
                self.interactions[key].clear()
                return command
    def clear(self):
        for key in self.interactions.keys():            
            self.interactions[key].clear()

class TFHill_widget:
    def __init__(self,net_widget,interactions=None):
        self.key = "TFHill"
        self.net = net_widget
        self.interactions = interactions
        self.is_on = w.ToggleButton(description="TFHill",value=False,layout=w.Layout(width="15%", height="40px"))
        self.s1 = w.Dropdown(description="TF",options={"S{}".format(key):key for key in self.net.types["TF"]})
        self.s2 = w.Dropdown(description="Gene",options={"S{}".format(key):key for key in self.net.types["Genes"]})
        self.custom_p = w.Checkbox(description="custom parameters",value=False)
        
        self.n = w.FloatText(description="n",value=1,step=0.0001)
        self.h = w.FloatText(description="h",value=1,step=0.0001)
        self.activity = w.Dropdown(description="Type",options={"activation":1,"repression":0})
        self.update_param_state(False)
        self.update_inter_state(False)
        
        w.interactive(self.update_inter_state,button=self.is_on)
        w.interactive(self.update_param_state,button=self.custom_p)
        
    def update_param_state(self,button):
        state = not (self.custom_p.value and self.is_on.value)
        self.n.disabled = state
        self.h.disabled = state
        self.activity.disabled = state
    def update_inter_state(self,button):
        if self.is_on.value:
            for key in self.interactions.interactions.keys():
                if key!=self.key:
                    self.interactions.interactions[key].clear()
        self.s1.disabled = not self.is_on.value
        self.s2.disabled = not self.is_on.value
        self.s1.options = {"S{}".format(key):key for key in self.net.types["TF"]}
        self.s2.options = {"S{}".format(key):key for key in self.net.types["Genes"]}
        self.custom_p.disabled = not self.is_on.value
        
        self.update_param_state(1)
    def clear(self):
        self.is_on.value = False
        self.n.value=1
        self.h.value=1
        self.is_on.value = False
        self.custom_p.value = False
        self.update_param_state(1)
        self.update_inter_state(1)
    def get_widget(self):
        to_return = [w.HBox([self.is_on,self.s1,self.s2,self.activity])]
        to_return += [w.HBox([self.custom_p,self.h,self.n])]
        return w.VBox(to_return)
    def get_command(self):
        if self.custom_p.value:
            n = self.n.value
            h = self.h.value
            command = "\tnet.new_TFHill(s_d[{TF}],{n},{h},tm_d[{TM}],activity={a})".format(TF=self.s1.value,TM=self.s2.value,n=n,h=h,a=self.activity.value)

        else:
            command = "\tnet.new_random_TFHill(s_d[{TF}],tm_d[{TM}])".format(TF=self.s1.value,TM=self.s2.value)
        self.net.inter_list.append("<p>TFHill (S{} -> S{})</p>".format(self.s1.value,self.s2.value))
        return command

class Phosphorylation_widget:
    def __init__(self,net_widget,interactions=None):
        self.key = "Phosphorylation"
        self.net = net_widget
        self.interactions = interactions
        self.is_on = w.ToggleButton(description="Phosphorylation",value=False,layout=w.Layout(width="15%", height="40px"))
        self.s1 = w.Dropdown(description="Phosphorylable",options={"S{}".format(key):key for key in self.net.types["Phosphorylable"]})
        self.s2 = w.Dropdown(description="Kinase",options={"S{}".format(key):key for key in self.net.types["Kinase"]})
        self.custom_p = w.Checkbox(description="custom parameters",value=False,layout=w.Layout(width="23%", height="40px"))
        
        self.kp = w.FloatText(description="k_p",value=1,step=0.0001,layout=w.Layout(width="15%", height="40px"))
        self.n = w.FloatText(description="n",value=1,step=0.0001,layout=w.Layout(width="15%", height="40px"))
        self.h = w.FloatText(description="h",value=1,step=0.0001,layout=w.Layout(width="15%", height="40px"))
        self.kd = w.FloatText(description="k_d",value=1,step=0.0001,layout=w.Layout(width="15%", height="40px"))
        
        self.update_param_state(False)
        self.update_inter_state(False)
        w.interactive(self.update_inter_state,button=self.is_on)
        w.interactive(self.update_param_state,button=self.custom_p)
        
    def update_param_state(self,button):
        state = not (self.custom_p.value and self.is_on.value)
        self.n.disabled = state
        self.h.disabled = state
        self.kp.disabled = state
        self.kd.disabled = state
        
    def update_inter_state(self,button):
        if self.is_on.value:
            for key in self.interactions.interactions.keys():
                if key!=self.key:
                    self.interactions.interactions[key].clear()
        self.s1.disabled = not self.is_on.value
        self.s2.disabled = not self.is_on.value
        self.s1.options = {"S{}".format(key):key for key in self.net.types["Phosphorylable"]}
        self.s2.options = {"S{}".format(key):key for key in self.net.types["Kinase"]}
        self.custom_p.disabled = not self.is_on.value
        self.update_param_state(1)
        
    def clear(self):
        self.is_on.value = False
        self.n.value=1
        self.h.value=1
        self.kp.value=1
        self.kd.value=1
        self.is_on.value = False
        self.custom_p.value = False
        self.update_param_state(1)
        self.update_inter_state(1)
    def get_widget(self):
        to_return = [w.HBox([self.is_on,self.s1,self.s2])]
        to_return += [w.HBox([self.custom_p,self.kp,self.h,self.n,self.kd])]
        return w.VBox(to_return)
    def get_command(self):
        ind = len(self.net.species)
        if self.custom_p.value:
            n = self.n.value
            h = self.h.value
            kp = self.kp.value
            kd = self.kd.value
            
            command = "\ts_d[{ind}],temp = net.new_Phosphorylation(s_d[{s2}],s_d[{s1}],{kp},{h},{n},{kd})".format(s1=self.s1.value,s2=self.s2.value,n=n,h=h,kp=kp,kd=kd,ind=ind)
        else:
            command = "\ts_d[{ind}],temp = net.new_random_Phosphorylation(s_d[{s2}],s_d[{s1}])".format(s1=self.s1.value,s2=self.s2.value,ind=ind)
        self.net.species.append(ind)
        types = eval("self.net.s_d[{s1}].list_types()".format(s1=self.s1.value))

        species_tags = list(set(types))
        species_tags.append("Phospho")
        if "Input" in species_tags:species_tags.remove("Input")
        if "Output" in species_tags:species_tags.remove("Output")
        for key in self.net.types.keys():
            if key in types:
                self.net.types[key].append(ind)
        
        self.net.species.append(ind)
        self.net.species_list.append("<p>S{} ({})</p>".format(ind,",".join(species_tags)))
        self.net.inter_list.append("<p>Phosphorylation (S{} -(S{})-> S{})</p>".format(self.s1.value,self.s2.value,ind))        
        return command
    
class PPI_widget:
    def __init__(self,net_widget,interactions = None):
        self.key = "PPI"
        self.net = net_widget
        self.interactions = interactions
        self.is_on = w.ToggleButton(description="PPI",value=False,layout=w.Layout(width="15%", height="40px"))
        self.s1 = w.Dropdown(description="Species 1",options={"S{}".format(key):key for key in self.net.types["Complexable"]})
        self.s2 = w.Dropdown(description="Species 2",options={"S{}".format(key):key for key in self.net.types["Complexable"]})
        self.custom_p = w.Checkbox(description="custom parameters",value=False,layout=w.Layout(width="23%", height="40px"))
        
        self.kp = w.FloatText(description="k+",value=1,step=0.0001,layout=w.Layout(width="15%", height="40px"))
        self.km = w.FloatText(description="k-",value=1,step=0.0001,layout=w.Layout(width="15%", height="40px"))
        
        self.update_param_state(False)
        self.update_inter_state(False)
        w.interactive(self.update_inter_state,button=self.is_on)
        w.interactive(self.update_param_state,button=self.custom_p)
  
        
    def update_param_state(self,button):
        state = not (self.custom_p.value and self.is_on.value)
        self.kp.disabled = state
        self.km.disabled = state\
        
    def update_inter_state(self,button):
        if self.is_on.value:
            for key in self.interactions.interactions.keys():
                if key!=self.key:
                    self.interactions.interactions[key].clear()
        self.s1.disabled = not self.is_on.value
        self.s2.disabled = not self.is_on.value
        self.s1.options = {"S{}".format(key):key for key in self.net.types["Complexable"]}
        self.s2.options = {"S{}".format(key):key for key in self.net.types["Complexable"]}
        self.custom_p.disabled = not self.is_on.value
        self.update_param_state(1)
    def clear(self):
        self.is_on.value = False
        self.kp.value=1
        self.km.value=1
        self.is_on.value = False
        self.custom_p.value = False
        self.update_param_state(False)
        self.update_inter_state(False)
    def get_widget(self):
        to_return = [w.HBox([self.is_on,self.s1,self.s2])]
        to_return += [w.HBox([self.custom_p,self.kp,self.km])]
        return w.VBox(to_return)
    def get_command(self):
        ind = len(self.net.species)
        if self.custom_p.value:
            kp = self.kp.value
            km = self.km.value            
            command = "\ttemp,s_d[{ind}] = net.new_PPI(s_d[{s2}],s_d[{s1}],{kp},{km})".format(s1=self.s1.value,s2=self.s2.value,kp=kp,km=km,ind=ind)
        else:
            command = "\ttemp,s_d[{ind}] = net.new_random_PPI(s_d[{s2}],s_d[{s1}])".format(s1=self.s1.value,s2=self.s2.value,ind=ind)
        types1 = eval("self.net.s_d[{s1}].list_types()".format(s1=self.s1.value))
        types2 = eval("self.net.s_d[{s1}].list_types()".format(s1=self.s2.value))
        species_tags = list(set(types1 + types2 + ["Complex","Degradable","Phosphorylable"]))
        species_tags.remove("Complexable")
        if "Input" in species_tags:species_tags.remove("Input")
        if "Output" in species_tags:species_tags.remove("Output")
        self.net.species.append(ind)
        self.net.inter_list.append("<p>PPI (S{} + S{} -> S{})</p>".format(self.s1.value,self.s2.value,ind))
        self.net.species_list.append("<p>S{} ({})</p>".format(ind,",".join(species_tags)))
        return command
