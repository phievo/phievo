import ipywidgets as w
from IPython.display import display

class float_range_widget:
    def __init__(self,name,arg):
        bad_arg_msg = "A range widget must be initialize with a list of size two or a positive float."
        if type(arg) is list:
            assert len(arg)==2,bad_arg_msg
            m,M = arg
        else:
            assert arg>0,bad_arg_msg
            m,M = 0,arg
        self.name = name
        name_w = w.HTML(value="<p>{}:</p>".format(name))
        self.m_w = w.BoundedFloatText(value=m,description='min',step=0.1,layout=w.Layout(width='20%', height='80px'),min=0,max=M)
        self.M_w = w.BoundedFloatText(value=M,description='max',step=0.1,layout=w.Layout(width='20%', height='80px'),min=0,max=100000)
        self.link = w.jslink((self.m_w,"max"),(self.M_w,"value"))
        
        self.widget = w.HBox([name_w,self.m_w,self.M_w])
    def get_widget(self):
        return self.widget
    def get_values(self):
        if self.m_w.value==0:
            return self.M_w.value
        else:
            return [self.m_w.value,self.M_w.value]

class int_range_widget:
    def __init__(self,name,arg):
        bad_arg_msg = "A range widget must be initialize with a list of size two or a positive int."
        if type(arg) is list:
            assert len(arg)==2,bad_arg_msg
            m,M = arg
        else:
            assert arg>0,bad_arg_msg
            m,M = 0,arg
        self.name = name
        name_w = w.HTML(value="<p>{}:</p>".format(name))
        self.m_w = w.BoundedIntText(value=m,description='min',step=1,layout=w.Layout(width='20%', height='80px'),min=0,max=M)
        self.M_w = w.BoundedIntText(value=M,description='max',step=1,layout=w.Layout(width='20%', height='80px'),min=0,max=100000)
        self.link = w.jslink((self.m_w,"max"),(self.M_w,"value"))
        
        self.widget = w.HBox([name_w,self.m_w,self.M_w])
    def get_widget(self):
        return self.widget
    def get_values(self):
        if self.m_w.value==0:
            return self.M_w.value
        else:
            return [self.m_w.value,self.M_w.value]


