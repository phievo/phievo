import ipywidgets as w
from IPython.display import display
import random
from phievo.Networks import mutation

tag_choices = {
    "Species":(True,),
    "Degradable":(True,"rate","Float",0.5,"Species.degradation"),
    "TF":(True,"activity","Int",1),
    "Kinase":(False,),
    "Phosphatase":(False,),
    "Output":(False,"index","Int"),
    "Input":(False,"index","Int"),
    "Complexable":(False,),
    "Complex":(False,),
    "Phosphorylable":(False,),
    "Diffusible":(False,"rate","Float",0,"Species.diffusion")
}
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
        
        self.g_basal= w.BoundedFloatText(description="basal:",value=0,step=0.001,layout=w.Layout(width="15%", height="40px"),disabled=True)
        self.g_rate = w.BoundedFloatText(description="rate:",value=1,step=0.001,layout=w.Layout(width="15%", height="40px"),disabled=True)
        self.g_delay = w.BoundedFloatText(description="delay:",value=0,step=0.001,layout=w.Layout(width="15%", height="40px"),disabled=True)
        self.g_custom = w.Checkbox(description="Custom parameters",value=False,disabled=False)
        
        self.message =w.HTML("")
        def input_not_gene1(val):
            self.isGene.disabled = self.tags["Input"].is_on.value
            #self.tags["Input"].param.disabled = not val
        def input_not_gene2(val):
            self.tags["Input"].is_on.disabled = self.isGene.value 
            self.g_basal.disabled = not self.isGene.value or not self.g_custom.value
            self.g_delay.disabled = not self.isGene.value or not self.g_custom.value
            self.g_rate.disabled = not self.isGene.value or not self.g_custom.value
            self.g_custom.disabled = not self.isGene.value 
        w.interactive(input_not_gene1,val=self.tags["Input"].is_on)
        w.interactive(input_not_gene2,val=self.isGene)
        w.interactive(input_not_gene2,val=self.g_custom)
                        
    def get_widget(self):
        return w.VBox(
                      [w.HTML(value="<h3>gene settings:</h3><p>If \"gene\" is selected, the program also adds a TModule and a CorePromoter and the species is regulated.</p><p>The parameters that are not customized will be drawn randomly in the space of the possible values.</p><p>An input species cannot be a gene.</p>")]
                      +[w.HBox([self.isGene,self.g_rate,self.g_basal,self.g_delay,self.g_custom])]\
                      +[w.HTML(value="<h3>add tags:</h3><p>The parameters that are not customized will be drawn randomly in the space of the possible values or will take the next available index (for Inputs and Outputs).</p>")]\
                      +[self.tags[key].get_widget() for key in self.tags.keys()]\
                      +[w.HBox([self.add_button,self.message])])
    def clear(self,b):
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
            
    def get_command(self):
        tags = [tag for tag in [self.tags[key].get_tag() for key in self.tags.keys()] if tag]
        param = "\tparam = {}".format(tags_to_str(tags))
        if self.isGene.value:
            if self.g_custom.value:
                function = "\ttm_d[ind],prom_d[ind],s_d[ind] = net.new_gene({rate},{delay},param,{basal})".format(rate=self.g_rate.value,delay=self.g_delay.value,basal=self.g_basal.value)
            else:
                function = "\ttm_d[ind],prom_d[ind],s_d[ind] = net.random_gene(param)"
        else:
            function = "\ts_d[ind] = net.new_Species(param)"
        return [param,function],tags
                
