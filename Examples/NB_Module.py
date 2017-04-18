from  phievo.AnalysisTools.Notebook import Notebook,CellModule
from ipywidgets import interact, interactive, widgets
from IPython.display import display

class DisplayFitness(CellModule):
    def __init__(self,Notebook):
        super(DisplayFitness, self).__init__(Notebook)
        self.button = widgets.Button(description="Display fitness",disabled=True)
        self.display_area = widgets.HTML(value=None, placeholder='<p></p>',description='Fitness:')
        self.notebook.dependencies_dict["seed"].append(self)
        self.notebook.dependencies_dict["generation"].append(self)
        self.notebook.dependencies_dict["project"].append(self)
    def update_display(self,button):
        """
        Custom function that handles the button click and wrtie the fitness in the HTML widget.
        """
        seed = self.notebook.seed
        gen = self.notebook.generation
        fit = str(self.notebook.sim.seeds[seed].generations[gen]["fitness"])
        self.display_area.value = "<p>{0}</p>".format(fit)
    def update(self):
        """
        Clear the HTML text and when the seed or the generation is updated.
        """
        if self.notebook.sim is None or self.notebook.seed is None or self.notebook.generation is None:
            self.button.disabled=True
        else:
            self.button.disabled=False
        self.display_area.value="<p></p>"
    def display(self):
        """
        Display the button and the display area on one row.
        """
        self.button.on_click(self.update_display)
        display(widgets.HBox([self.button,self.display_area]))
