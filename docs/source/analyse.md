# Results and Analysis Tools

φ-evo has a module dedicated to the analysis of the results. The results are stored in a *Simulation* object that contains a set of method that give a quick access to the most relevant observables of a run. To start analyzing the `evo_dir` project, you need to create a *Simulation* object associated to it.

```python
from phievo.AnalysisTools import Simulation

sim = Simulation("evo_dir")
```

From there it is pretty straight forward to explore the architecture of the results. A simulation contains Seeds which themselves contain Networks. In order not to overload the memory, the Seeds only store a link to the networks. As an example, here is how you would load the best network for generation 350 in the seed number 2:

```python
sim = Simulation("evo_dir")
best_net_2_350 = sim.Seeds[2].get_best_net(350)
# Equivalent to
sim = Simulation("evo_dir")
best_net_2_350 = sim.get_best_net(2,350)
```

## Organization  of the results

If you want to understand why the *Simulation* object is organized the way it is and how to go beyond its possibilities, you need to have an idea of how φ-evo stores the results of a simulation.

By default, for every generation *g* only one *Network* is stored using pickle in a file labelled `Bests_g.net`. When the simulation has only one fitness objective, this network is the one with the best fitness in the population. However when the evolution is run using a multiobjective criterium (like pareto optimisation), the best net is chosen randomly among the network of rank 1.

The former storing method limits the disk space usage. However you might want to store the whole population either for restarting the algorithm from a given generation or to analyze every member of the generation. To add this feature, you can specify a storing period by setting the `prmt['restart']['freq']` parameter in the initialization file before launching the simulation. For example, if you set it to 50, the complete population will be stored every 50 generations in a python *shelve* named `restart_file`.

Other files created:

 - `data` is a quick access shelve file to certain informations stored as lists at the following keys:
    - *generation*: index of the generation
    - *fitness*: fitness of the best network
    - *n_species*: number of species in the best network
    - *n_interactions*: number of interaction in the best network
- `parameters` is a copy of the parameter dictionnaries (defined for the non default in the initialization file) that were used during the simulation.
- `log_#.c` Copy of the input, fitness, history, etc. **c** files used for the simulation.
- `log_init_file.py` Copy of the init file used for the simulation

## Analysis Tools

In this section we will explore the built-in functions that are bound to a *Simulation* object.


### custom_plot

Plots two observables one against each other for a given seed. The available observables are the ones contained in the `data` file ("generation", "fitness", "n_species", "n_interactions").

```python
sim.seeds[1].custom_plot("generation","fitness")
# Similarly you can use the shortcut
sim.custom_plot(1,"generation","fitness")
```

### plot_fitness
There also exists a method to plot the fitness directly:

```python
sim.seeds[1].show_fitness()
# or
sim.show_fitness(1)
```
### get_best_net
Get the best net found in a given generation (the function reads the `Bests_g.net` file and return the Network object)

```python
bestnet_g5_seed3 = sim.seeds[3].get_best_net(5)
# or
bestnet_g5_seed3 = sim.get_best_net(3,5)
```

### get_backup_net

If you want to extract a network from a entirely stored generation, you can use *get_backup_net*. Be careful though, not every population is stored in the `restart_file`. You can use the `stored_generation_indexes` to check which generation has been stored.
```python
net8_g50_seed3 = sim.seeds[3].get_backup_net(50,8)
# Or
net8_g50_seed3 = sim.get_backup_net(3,50,8)
```

### stored_generation_indexes
 The *stored_generation_indexes* is method that returns the list of stored generations.

```python
list_stored = sim.seeds[1].stored_generation_indexes()
# Or
list_stored = sim.stored_generation_indexes(1)
```

### Read a network from the pickle file
The simulation stores the best networks of every generation in the name `Bests_#.net`. This is only a pickle file and can be read manually using the pickle library:

```python3
import pickle

with open("Bests_#.net","rb") as net_file:
    net = pickle.load(net_file)
```
Or using the φ-evo function:
```python3
import phievo

phievo.read_network("Bests_#.net")
```

### Running a network's dynamics

By construction φ-evo does not allow to quickly run the dynamics of a network. Because the dynamics is computed in C (for performance reason), a python Network object does not have a method that directly returns the derivative at a given state of gene quantities. However φ-evo has the method `run_dynamics` to symplify the run of a dynmics for a given network based on the history and inputs defined in *init_history.c* and *input.c* respectively.

```python
net = sim.get_best_net(3,5)
dyn_buffer = sim.run_dynamics(net=net,trial=1)
```

 You can specify the number of trial you want to run (if the dynamics is stochastic for example). The buffer returned by the function is a dictionary where the "time" and "net" keys give you access to  the time vector and the network used for the run respectively. The other keys are the index of the trial for which you want to access the data. Note that the buffer is also stored in the *Simulation.buffer_data*, the latter is erased every time you run a new set of dynamics for *Simulations*.

### Plotting the results of a dynamics

The simulation object allows you to plot the two results you would like to see after running a dynamics:

1) The time course of the genes in a given cell with *Plot_TimeCourse*
2) The evolution of the genes along the system at a given time point with *Plot_Profile*

```python
sim.Plot_TimeCourse(trial_index=1,cell=1)
sim.Plot_Profile(trial_index=1,time=1)
```

### Draw a network's layout
The network object contains a function to draw the layout of its gene interactions:
```python
net = sim.get_best_net(3,5)
net.draw()
```

## Notebook

To facilitate the use of the former functions, φ-evo as a class *Notebook* that is used to run them in a [jupyter notebook](https://jupyter.org).

All the functions described previously can be used directly in a jupyter notebook but the *Notebook* class improves the usability by handling the dependencies between widgets. For instance you want the module in charge of plotting a network's layout to be disabled as long as a Seed and a Network have not been selected.

A Notebook  object serves as a container for all the available modules you can use in the jupyter notebook. A module contains the material to handle a cell: its widgets, some update functions and a display function that displays the widgets in the jupyter notebook. In the end, the user  only needs to run `myNotebook.myModule.display()` to create a jupyter elementary app in a cell. Then the module should be able to handle the expected inputs from the user.

### Creating a custom module

Every module of contained in the *Notebook* inherits from the *CellModule* class. The latter is a minimal template used to constrain the requirements a module must have:

- `__init__(self,Notebook)` : The init function takes the *Notebook* it is contained in as an argument.
- `update(self)` : If the module has dependencies, this function must be defined. When dependency is updated, this function is called.
- `display(self)`: The function must be redefined to display the widgets and to handle the relation between them.

#### \_\_init\_\_
This is the function where you define the different widgets for the module. It is also here that you define the dependencies of the module or create a new ones. The dependencies system allows communication between different *CellModules*.

```python
## Inform the notebook that MyModule depends on the Seed
self.notebook.dependencies_dict["seed"].append(self)
## Creates a dependencies
self.notebook.dependencies_dict["dep_name"] = []
```  
Note that if you create a new dependency, you should make sure that you also handle the updates when the dependency changes:

```python
for cell in self.notebook.dependencies_dict["dep_name"]:
    cell.update()
```

#### update

Every module, particularly those with dependencies, should have an update function. This is the function to call when the dependency is changed. The update function can do whatever you want but mostly its purpose is to unable/disabled the widgets when a dependency is changed or to reset their options.

In Addition to the *self.notebook.dependencies_dict*, a module can access the dictionnary *self.notebook.extra_variables* to pass values between *CellWidgets*.

#### display

The display function is here to contain the interaction and display code you would normally put in a jupyter notebook to handle the communication of the widgets with the functions.

The philosophy of the *CellModule* is to create an elementary app in charge of one action (plotting a curve, setting the seed, etc.). Using a module's display method in a cell gives access to the app at this location.

#### Other functions

The *update* and *dispay* functions are usually not enough to run the *CellModule*. You will need to define custom methods for your module to handle the widget interactions(for instance, what happens when a widget is clicked?).

#### Example: DisplayFitness

Here is a little practical example on how to include a custom *CellModule* that displays the best fitness of the selected generation when the button is clicked.

Create a module file *NB_Module.py*  and import the Notebook module and some widget libraries:

```python
from  phievo.AnalysisTools.Notebook import Notebook,CellModule
from ipywidgets import interact, interactive, widgets
from IPython.display import display
```

Then create the *CellModule* object:

```python
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
```

Save the file and open the notebook to associate the newly created module to a notebook object.

```python
...
from  phievo.AnalysisTools.Notebook import Notebook
import NB_Module

notebook = Notebook()
setattr(notebook,"display_fitness",NB_Module.DisplayFitness(notebook))
```

Now the display_fitness module can be used as any other *CellModule* by creating a new cell and running:

```python
notebook.display_fitness.display()
```

A copy of the [*NB_Module.py*](https://raw.githubusercontent.com/phievo/phievo/master/Examples/NB_Module.py) file is available in the *Examples/* directory.
