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
        infos = w.HTML(value="<h2>{}</h2>{}".format(self.name,self.infos))
        return w.VBox([infos]+[self.obj_dict[key].get_widget() for key in self.obj_dict.keys()])
    def get_values(self):
        return {key:self.obj_dict[key].get_values() for key in self.obj_dict.keys()}

    def set_values(self,values):
        for key,val in values.items():
            self.obj_dict[key].set_values(val)
