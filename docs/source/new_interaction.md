## Create a new interaction
φ-evo allows you to add a custom interaction that is not available in the default list.  
To do so  you need to write an interaction file. 

To make make the explanation clearer, we will explain how to build a new interaction on a real example of a methylation interaction.

The methylation adds a methyl group to a species $S$. The methylated species is denoted with a ${}^{*}$ symbol:
$$ S \leftrightarrow S^{*} $$

We choose the simplest kinetics for this reaction:

$$ \frac{d S^{*}}{d t} = -\frac{d S}{d t} = k_f S - k_b S^{*} $$

Let us start by creating the *Methyl.py* in a project directory.

### Imports

Every interaction depends on the following φ-evo modules:

 - classes_eds2: for the core structure of the intaction
 - mutation: to handle mutation
 - deriv2: to explain how to generate the C code associated to the new mutation
 
```python
# In Methyl.py

from phievo import __silent__,__verbose__
if __verbose__:
    print("Execute Methyl (Interaction Template)")

from phievo.Networks import mutation
from phievo.Networks import deriv2
from phievo.Networks import classes_eds2
import copy
```
### Define a new type of species
Only methylable species can be methylated. For now φ-evo does not know how to create a methylable species and what are its characteristics. There should be a few line telling how to do it:

```python
# In Methyl.py
mutation.species_types["Methylable"] = lambda random_generator:[
    ["Methylable"],
    ['Diffusible',mutation.sample_dictionary_ranges('Species.diffusion',random_generator)]
]
classes_eds2.Species.Tags_Species["Methylable"] = []
```
In the above lines, we tell φ-evo that a Methylable species has two characteristics:

 - Methylable: obviously
 - Diffusable: An extra characteristic is added to show how one would add a characteristics that comes with a parameter. A lambda function allows the program to generate new parameters when a new species is created.
	
**Note:** You can use the `mutation.sample_dictionary_ranges` to sample a random variable whose range has been define in `dictionary_ranges` in the *init* file.
### Set the default ranges for the parameters

```python
# In Methyl.py

### Define the default dictionary_range
mutation.dictionary_ranges['Methyl.methyl'] = 0.0/(mutation.C*mutation.T)
mutation.dictionary_ranges['Methyl.demethyl'] = 0.0/mutation.T

```

### Define the *Methyl* class

Every interaction in φ-evo inherits from the *classes_eds2.Interaction*:

```python
# In Methyl.py
class Methyl(classes_eds2.Interaction):
    """
    Methylation interaction

    Args:
        Methyl.methyl(float): binding rate of a methyl group
        Methyl.demethyl(float): unbinding rate of a methyl group
        label(str): Methylation
        input(list): Type of the inputs
        output(list): Type of the outputs
    """
    def __init__(self,methyl=0,demethyl=0):
        classes_eds2.Node.__init__(self)
        self.methyl=methyl
        self.demethyl=demethyl
        self.label='Methylation'
        self.input=['Methylable']
        self.output=['Species']

    def __str__(self):
        """
        Used by the print function to display the interaction.
        """
        return "{0.id} Methylation: methyl. = {0.methyl:.2f}, demethyl = {0.demethyl:.2f}".format(self)

    def outputs_to_delete(self,net):
        """
        Returns the methylated form of the species to delete when the reaction is deleted.
        """
        return net.graph.successors(self)

```

The interaction's methods are the following:

 - `__init__`: Creates the interaction object
 - `__str__`: Produces the string used by the print function
 - `outputs_to_delete`: Function that tells what are the species that were added to the network when the interaction was built and that need to be deleted when the interaction is removed.
 
### Handling the mutation

The program needs five functions to tell φ-evo how to add the mutation via a mutation

#### number_Methyl

Evaluate the number of possible interactions of type *Methyl* that can be added to the network. This number is used to verify that the actual number of possible mutation found in `random_Methyl` is consistant with our intuition.

```python
# In Methyl.py

def number_Methyl(self):
    """
    Returns the number of possible methylation in the current network.
    Note: this function is optional, it is used to check the consistency of
    the random_Methyl function.
    """
    n = self.number_nodes('Methylable')
    n_Methyl = self.number_nodes('Methyl')
    return n-n_Methyl
```
 
#### new_Methyl

This is the function that adds the *Methyl* interaction to the Network. It creates both a *Methyl* interaction and a *methylated species*.

```python
# In Methyl.py
def new_Methyl(self,S,methyl,demethyl,parameters):
    """
    Creates a new :class:`Networks.Methyl.Methyl` and the species methylated for in the the network.

    Args:
        S: species to methylate
        methyl(float): binding rate of a methyl group
        demethyl(float): unbinding rate of a methyl group
        parameters(list): Parameters of the methylated species
    Returns:
        [methyl_inter,S_methyl]: returns a Methyl interaction and a methylated species.
    """

    S_methyl = classes_eds2.Species(parameters)
    meth_inter = Methyl(methyl,demethyl)
    assert meth_inter.check_grammar([S],[S_methyl]),"Error in grammar, new Methylation"
    self.add_Node(S_methyl)
    self.add_Node(meth_inter)
    self.graph.add_edge(S,meth_inter)
    self.graph.add_edge(meth_inter,S_methyl)
    return [meth_inter,S_methyl]
```

**Note:** Then function needs a list of characteristics for the methylated species created. It is provide via `parameters`.

#### new_random_Methyl

Wrapping of the `new_Methyl` function. It generates randomly the rate of the methylation and the parameters of the methylated species created.


```python
# In Methyl.py

def new_random_Methyl(self,S):
    """
    Creates a methylation with random parameters.
        
    Args:
        S: Species to methylate
    Returns:
        [methyl_inter,S_methyl]:returns a Methyl interaction and a methylated species.
    """
    parameters = {}
    if S.isinstance("TF"):
        parameters['TF'] = self.Random.random()*2
    for tt in S.types:
        if tt not in ["TF","Methylable","Input","Output"]:
            parameters[tt] = [mutation.sample_dictionary_ranges('Species.{}'.format(attr),self.Random) for attr in S.Tags_Species[tt]]

    # Transform to fit phievo list structure
    parameters = [[kk]+val if val else [kk] for kk,val in parameters.items()]
    methyl = mutation.sample_dictionary_ranges('Methyl.methyl',self.Random)
    demethyl = mutation.sample_dictionary_ranges('Methyl.demethyl',self.Random)
    return self.new_Methyl(S,methyl,demethyl,parameters)
    

```

#### random_Methyl

Function called by the φ-evo to add a new Methylation interaction to the network during the evolution. It chooses a methylable species randomly and calls `new_random_Methyl` to add a methylation to this species.

```python
# In Methyl.py

def random_Methyl(self):
    """
    Evaluates the species that can be phosphorilated, picks one an create a random
    methylation. The random mutation is made using :func:`new_random_Methyl <phievo.Networks.classes_eds2.new_random_Methyl>`.

    Returns:
        [methyl_inter,S_methyl]: returns a Methyl interaction and a methylated species.
    """
    try:
        list_methylable=self.dict_types["Methylable"]
    except KeyError:
        print("\tThe network contain no Methylacble species.")
        raise
    list_possible_methylable = []
    for S in list_methylable:
        if not self.check_existing_binary([S],"Methyl"):
            list_possible_methylable.append(S)
    n_possible = len(list_possible_methylable)
    assert n_possible==self.number_Methyl(),"The number of possible new methylation does not match its theoretical value."
    if n_possible==0:
        if __verbose__:
            print("No more possible methylation.")
        return None
    else:
        S = list_possible_methylable[int(self.Random.random()*n_possible)]
        return self.new_random_Methyl(S)
        

```

#### Methyl_deriv_inC

Function that generates the C code string of the interaction kinetics.

```python
# In Methyl.py

def Methyl_deriv_inC(net):
    """
    Function called to generate the string corresponding to in a methylation in C.
    """
    func_str = "\n/************** Methylations *****************/\n"
    methylations = net.dict_types.get("Methyl",[])
    for methyl_inter in methylations:
        S = net.graph.predecessors(methyl_inter)[0]
        S_meth = net.graph.successors(methyl_inter)[0]
        f_rate = "{M.methyl}*{S.id}".format(M=methyl_inter,S=S)
        b_rate = "{M.demethyl}*{S_m.id}".format(M=methyl_inter,S_m=S_meth)

        func_str += deriv2.compute_leap([S.id],[S_meth.id],f_rate)
        func_str += deriv2.compute_leap([S_meth.id],[S.id],b_rate)
    return func_str
```


### Bind the code to φ-evo

The last step is to add all the functions written previously to the default `Mutable_Network`.

```python
# In Methyl.py
setattr(classes_eds2.Network,"number_Methyl",number_Methyl)
setattr(classes_eds2.Network,"new_Methyl",new_Methyl)
setattr(classes_eds2.Network,"new_random_Methyl",new_random_Methyl)
setattr(classes_eds2.Network,"random_Methyl",random_Methyl)
deriv2.interactions_deriv_inC["Methyl"] = Methyl_deriv_inC
```

You can download [Methyl.py](https://github.com/phievo/phievo/raw/master/Examples/Methyl.py) from φ-evo's examples
### Edit the init file to load Methyl

The top of the init file should now be able to load the Methyl module with an import if the two files are in the same directory:

```python
# In initialization.py
import Methyl
```

Now the new mutation settings are made similarly to any of the default interaction:

```python
# In initialization.py

mutation.dictionary_ranges['Methyl.methyl'] = [0.1,1]
mutation.dictionary_ranges['Methyl.demethyl'] = [0.1,0.5]

dictionary_mutation['random_gene(\'Methylable\')'] = 0.1
dictionary_mutation['random_Interaction(\'Methyl\')']=0.1 
dictionary_mutation['remove_Interaction(\'Methyl\')']=0.01
....
```
