## Create a new interaction
 φ-evo allows you to add custom interactions that are not available in the default list.  
 To do so  you need to write an interaction file, here is a step by step explanation
 of how to do so

First create you `custom_interaction.py` file and write down the libraries to import:

```python
import classes_eds2
import mutation
import copy
import deriv2
```

An interaction inherits from the *classes_eds2.Interaction*:

```python
class Dummy_Interaction(classes_eds2.Interaction):
    """ A Dummy_Interaction """

    def __init__(self,a=0,d=0):
        """
        Initialize all the parameters and the allowed structure of
        the interaction (input and output)
        """
        classes_eds2.Node.__init__(self)
        self.parameter1=a
        self.parameter2=d
        self.label='Dummy_Interaction'
        # the input and ouput attributes are used by classes_eds2.check_grammar()
        # Typical values are 'Species', 'Degradable' or 'TModule'
        self.input=[#<'One Tag or TModule per input '>#]
        self.output=[#<'One Tag or TModule per input '>#]

    def outputs_to_delete(self,net):
        """
        Returns the objects to delete when removing the Dummy_Interaction,
        typically some of the outputs of the interaction or nothing
        """
        return [] #by default
        return net.graph.successors(self) # use this one when the interaction
                                          # creates a new species
                                          # (PPI, LR or Phospho typically)
def number_Dummy_Interaction(self):
    """
    Computes the number of possible new Dummy_Interactions,
    used when the algorithm determine the next mutation
    """
    # By default but you may want to optimize
    return len(self.list_possible_Dummy_Interaction())
```
This creates the  dummy interaction but does not tell the network how to use it.
(how to add it, how to delete it, etc.)

The following functions tell a network how to add the newly created interaction

- `new_Dummy_Interaction`: Function in charge of adding the new interaction. It involves
linking the species together and possibly creating new ones.
- `duplicate_Dummy_Interaction`: Reproduces an existing interaction between another set of species.
- `number_Dummy_Interaction`: Function used to get the total number of possible mutation new_Dummy_Interaction
in the network.

```python
def new_Dummy_Interaction(self,Input1, Input2,parameter1,parameter2,output_parameters):
    """
    Create a new Dummy_Interaction and its Outputs if necessary and add them
    to self (Network), should return a list of the interaction and
    the species created
    """
    output = classes_eds2.Species(output_parameters)
    interaction = Dummy_Interaction(parameter1,parameter2)
    if interaction.check_grammar([Input1,Input2], [output]):
        self.add_Node(output)
        self.add_Node(interaction)
        self.graph.add_edge(Input1,interaction)
        self.graph.add_edge(Input2,interaction)
        self.graph.add_edge(interaction,output)
        return [interaction,output]
    else:
        print "Error in grammar in creation of new_Dummy_Interaction "
        return None
## Attach the interaction
setattr(classes_eds2.Network,'new_Dummy_Interaction',new_Dummy_Interaction)

def duplicate_Dummy_Interaction(self,species,new_species,interaction,module,new_module):
    """function to duplicate an existing Dummy_Interaction interaction.
    """
    new_interaction=copy.deepcopy(interaction)
    ## ...
    ## Creates a the clone of the interaction between species between
    ## new_species in the network `self`
    ##
    return None
setattr(classes_eds2.Network,'duplicate_Dummy_Interaction',duplicate_Dummy_Interaction)

def number_Dummy_Interaction(self):
    """
    Computes the number of possible new Dummy_Interactions,
    used when the algorithm determine the next mutation
    """
    return len(self.list_possible_Dummy_Interaction())
setattr(classes_eds2.Network,'number_Dummy_Interaction',number_Dummy_Interaction)
```

In its architecture φ-evo needs a layer that comes over the *classes_eds2.Network*
to manage the mutation events. This class that inherits from *classes_eds2.Network*
is *mutation.Mutable_Network*. It also needs some adds on to know how to use the
new interaction.

```%%python
##### Attributes attached to the Mutable_Network Class #####


def random_Dummy_Interaction(self):
    """Creates a new random Dummy_Interaction among all possible ones by calling the new_random_Dummy_Interaction() method"""
    if self.list_types.has_key('key allowing this type of interaction'):
        list_possible = self.list_possible_Dummy_Interaction()
        if list_possible:
            #create randomly one Dummy_Interaction among those possible
            [Input1,Input2] = self.Random.choice(list_possible)

            a = mutation.sample_dictionary_ranges('Dummy_Interaction.parameter1',self.Random)
            d = mutation.sample_dictionary_ranges('Dummy_Interaction.parameter2',self.Random)
            new_Dummy_Interaction = self.new_Dummy_Interaction(Input1,Input2,a,d,output_parameters)
            return new_Dummy_Interaction
        else:
            print "In random_Dummy_Interaction : no other possible random_dummy_Interaction, Error"
            return None
    else:
        print "Error in random_Dummy_Interaction (try to create a Dummy_Interaction from non existing pieces)"
        return None

setattr(mutation.Mutable_Network,'random_Dummy_Interaction',new_random_Dummy_Interaction)
```

Finally you need to tell the *deriv2* how to write the c code associated to the new interaction.

``` python
def Dummy_Interaction_deriv_inC(net):
    """Return a string of C-code describing the equation for all the Dummy_Interactions in net"""
    func="\n/**************Dummy_Interaction interactions*****************/\n"
    # Loop over all Dummy_Interaction if there is at least one
    if net.list_types.has_key('Dummy_Interaction'):
        for index in net.list_types['Dummy_Interaction']:
            Output = net.graph.successors(index) #finds the Outputs
            [Input1,Input2] = net.graph.predecessors(index) #find the Inputs

            # defines interaction rate, input_id_list will decreases and output_id_list increases at rate rate1
            rate1 = String combining Input1.id, Input2.id, parameter1, parameter2
            func += deriv2.compute_leap([input_id_list],[output_id_list],rate1)

            rate2 = String combining Input1.id, Input2.id, parameter1, parameter2
            func += deriv2.compute_leap([input_id_list],[output_id_list],rate2)
    else:
        func += "No Dummy_Interaction in this network\n"
    return func

##### Update of the deriv2 method #####
deriv2.Dummy_Interaction_deriv_inC = Dummy_Interaction_deriv_inC
```
