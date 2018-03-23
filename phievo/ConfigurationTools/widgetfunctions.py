import ipywidgets as w
from IPython.display import display

box_width = "20%"
box_width2 = "10%"
description_width="25%"
box_height = "30px"

def range_test_arg(arg):
        bad_arg_msg = "A range widget must be initialize with a list of size two or a positive float."
        if type(arg) is list:
            assert len(arg)==2,bad_arg_msg
            m,M = arg
        else:
            assert arg>=0,bad_arg_msg
            m,M = 0,arg
        return m,M
    
class float_range_widget:
    def __init__(self,name,arg):
        m,M  = range_test_arg(arg)
       
        self.name = name
        name_w = w.HTML(value="<p>{}:</p>".format(name),layout=w.Layout(width=description_width, height=box_height))
        self.m_w = w.BoundedFloatText(value=m,description='min',step=0.001,layout=w.Layout(width=box_width, height=box_height),min=0,max=M)
        self.M_w = w.BoundedFloatText(value=M,description='max',step=0.001,layout=w.Layout(width=box_width, height=box_height),min=0,max=100000)
        self.link = w.jslink((self.m_w,"max"),(self.M_w,"value"))
        
        self.widget = w.HBox([name_w,self.m_w,self.M_w])
    def set_values(self,values):
        m,M = range_test_arg(values)        
        self.M_w.value = M
        self.m_w.max = M
        self.m_w.value = m
        
    def get_widget(self):
        return self.widget
    def get_values(self):
        if self.m_w.value==0:
            return self.M_w.value
        else:
            return [self.m_w.value,self.M_w.value]
    
    
class int_range_widget(float_range_widget):
    def __init__(self,name,arg):
        m,M  = range_test_arg(arg)
        self.name = name
        name_w = w.HTML(value="<p>{}:</p>".format(name),layout=w.Layout(width=description_width, height=box_height))
        self.m_w = w.BoundedIntText(value=m,description='min',step=1,layout=w.Layout(width=box_width, height=box_height),min=0,max=M)
        self.M_w = w.BoundedIntText(value=M,description='max',step=1,layout=w.Layout(width=box_width, height=box_height),min=0,max=100000)
        self.link = w.jslink((self.m_w,"max"),(self.M_w,"value"))        
        self.widget = w.HBox([name_w,self.m_w,self.M_w])
 

        
    
class bool_widget:
    def __init__(self,name,value=False):
        self.value = w.Checkbox(description=" ",value=bool(value),layout=w.Layout(wheight=box_width2))
        self.name = w.HTML(value="<p>{}:</p>".format(name),layout=w.Layout(width=description_width, height=box_height))
    def get_widget(self):
        return w.HBox([self.name,self.value])
    def get_values(self):
        return self.value.value
    def set_values(self,value):
        self.value.value = bool(value)
        
class int_widget(bool_widget):
    def __init__(self,name,value=1):
        super(int_widget,self).__init__(name,value)
        self.value = w.BoundedIntText(value=value,description=' ',step=1,layout=w.Layout(width=box_width, height=box_width2),min=0,max=1000000)
    def set_values(self,value):
        self.value.value = int(value)
        
class float_widget(bool_widget):
    def __init__(self,name,value=1):
        super(float_widget,self).__init__(name,value)
        self.value = w.BoundedFloatText(value=value,description=' ',step=0.001,layout=w.Layout(width=box_width, height=box_width2),min=0,max=1000000)
    def set_values(self,value):
        self.value.value = float(value)
        
class tags_widgets:
    def __init__(self,name,choices,values=""):
        self.name = w.HTML(value="<p>{} (press shift for multiple):</p>".format(name),layout=w.Layout(width=description_width, height=box_height))
        self.tags = w.SelectMultiple(description=" ",value=values,options=choices,layout=w.Layout(wheight=box_width2))
    def get_widget(self):
        return w.HBox([self.name,self.tags])
    def get_values(self):
        return list(self.tags.value)
    def set_values(self,values):
        self.tags.value = values

class w_restart:
    def __init__(self,values=None):
        self.name = "restart"
        self.widgets = {
            "activated":w.Checkbox(description="activated",value=False),
            "freq":w.BoundedIntText(description="Saving frequency",value=50,min=0,max=1000000,disabled=False),
            "kgeneration":w.HBox([w.BoundedIntText(description="Restart generation",value=0,min=0,max=1000000,disabled=True),
                          w.Checkbox(description=" ",value=False,disabled=True)]),
            "seed":w.HBox([w.BoundedIntText(description="Restart Seed",value=0,min=0,max=1000000,disabled=True),
                          w.Checkbox(description=" ",value=False,disabled=True)]),
            "same_seed":w.Checkbox(description="Restart with the same seed",value=True,disabled=True),
        }
        
        def run_activate(value):
            for key in ["kgeneration","seed"]:
                self.widgets[key].children[1].disabled = not value
            self.widgets["same_seed"].disabled = not value
        def activate_seed(value):
            self.widgets["seed"].children[0].disabled = not self.widgets["seed"].children[1].value
        def activate_kgeneration(value):
            self.widgets["kgeneration"].children[0].disabled = not self.widgets["kgeneration"].children[1].value
        w.interactive(run_activate,value=self.widgets["activated"])
        w.interactive(activate_seed,value=self.widgets["seed"].children[1])   
        w.interactive(activate_kgeneration,value=self.widgets["kgeneration"].children[1]) 
        if values:
            self.set_values(values)
    def set_values(self,values):
        self.widgets["freq"].value = values["freq"]
        self.widgets["activated"].value = values["activated"]
        if values["activated"]:
            values.setdefault("seed",None)
            values.setdefault("kgeneration",None)
            if values["seed"]:
                self.widgets["seed"].children[0].value=values["seed"]
                self.widgets["seed"].children[1].value=True
            else:
                self.widgets["seed"].children[1].value=False
            if values["kgeneration"]:
                self.widgets["kgeneration"].children[0].value=values["kgeneration"]
                self.widgets["kgeneration"].children[1].value=True
            else:
                self.widgets["kgeneration"].children[1].value=False
            self.widgets["same_seed"].value = values["same_seed"]
                
    def get_widget(self):
        return w.VBox([self.widgets[key] for key in self.widgets.keys()])
    def get_values(self):
        values = {}
        for key in self.widgets.keys():
            try:
                values[key] = self.widgets[key].value
            except AttributeError:
                values[key] = self.widgets[key].children[0].value if self.widgets[key].children[1].value else None
        return values

class code_widget:
    def __init__(self,name, default_value="",title="",height="300px"):
        self.header= w.HTML(value="<h3>{}:</h3>".format(title),layout=w.Layout(width="40%", height="50px"))
        self.code = w.Textarea(value=default_value,description="",layout=w.Layout(width="70%", height=height))
        self.button_reset = w.Button(description="Reset default")
        self.button_load = w.Button(description="Load",disabled=True)
        self.path = w.Text(description="file",placeholder="Path to a python file",value="")
        def check_file(path):
            if os.path.isfile(path) and path[-3:]==".py":
                self.button_load.disabled = False
                self.button_load.button_style="success"
            else:
                self.button_load.disabled = True
                self.button_load.button_style=""
        def load_file(button):
            with open(self.path.value,"rb") as code_file:
                self.code.value = code_file.read()
        w.interactive(check_file,path=self.path)
        self.button_load.on_click(load_file)
        def reset_code(button):
            self.code.value = default_initialization_code
        self.button_reset.on_click(reset_code)
    def get_widget(self):
        return w.VBox([self.header,self.code,w.HBox([self.path,self.button_load,self.button_reset])])
    def get_values(self):
        return self.code.value

class code_path_widget:
    def __init__(self,name,path=""):
        self.default = path
        self.name = w.HTML(value="<p>{}:</p>".format(name),layout=w.Layout(width="10%", height="50px"))
        self.value = w.Text(value=path,placeholder="Custom file",description=" ",disabled=True,layout=w.Layout(width="50%", height="50px"))
        self.activate = w.Button(description="Unlock customize",button_style="info",layout=w.Layout(height="35px"))
        def unlock(button):
            if self.value.disabled:
                self.value.disabled = False
                self.activate.description = "Reset"
            self.value.value = self.default
        def check_file(path):
            if path == self.default:
                self.activate.button_style = "info"
            elif os.path.isfile(path):
                self.activate.button_style = "success"
            else:
                self.activate.button_style = "danger"
        
        self.activate.on_click(unlock)
        w.interactive(check_file,path=self.value)
    def get_values(self):
        return self.value.value
    def get_widget(self):
        return w.HBox([self.name,self.value,self.activate])

    def set_values(self,value):
            self.value.value = value
