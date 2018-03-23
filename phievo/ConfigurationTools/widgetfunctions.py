import ipywidgets as w
from IPython.display import display



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
        name_w = w.HTML(value="<p>{}:</p>".format(name),layout=w.Layout(width='20%', height='40px'))
        self.m_w = w.BoundedFloatText(value=m,description='min',step=0.001,layout=w.Layout(width='20%', height='30px'),min=0,max=M)
        self.M_w = w.BoundedFloatText(value=M,description='max',step=0.001,layout=w.Layout(width='20%', height='30px'),min=0,max=100000)
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
        name_w = w.HTML(value="<p>{}:</p>".format(name))
        self.m_w = w.BoundedIntText(value=m,description='min',step=1,layout=w.Layout(width='20%', height='30px'),min=0,max=M)
        self.M_w = w.BoundedIntText(value=M,description='max',step=1,layout=w.Layout(width='20%', height='30px'),min=0,max=100000)
        self.link = w.jslink((self.m_w,"max"),(self.M_w,"value"))        
        self.widget = w.HBox([name_w,self.m_w,self.M_w])
 

    
class int_range_widget:
    def __init__(self,name,arg):
        bad_arg_msg = "A range widget must be initialize with a list of size two or a positive int."
        if type(arg) is list:
            assert len(arg)==2,bad_arg_msg
            m,M = arg
        else:
            assert arg>=0,bad_arg_msg
            m,M = 0,arg
        self.name = name
        name_w = w.HTML(value="<p>{}:</p>".format(name))
        self.m_w = w.BoundedIntText(value=m,description='min',step=1,layout=w.Layout(width='20%', height='30px'),min=0,max=M)
        self.M_w = w.BoundedIntText(value=M,description='max',step=1,layout=w.Layout(width='20%', height='30px'),min=0,max=100000)
        self.link = w.jslink((self.m_w,"max"),(self.M_w,"value"))
        
        self.widget = w.HBox([name_w,self.m_w,self.M_w])
    def get_widget(self):
        return self.widget
    def get_values(self):
        if self.m_w.value==0:
            return self.M_w.value
        else:
            return [self.m_w.value,self.M_w.value]

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
