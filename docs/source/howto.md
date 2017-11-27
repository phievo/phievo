HowTo
=====

This tutorial lists a series examples on how to perform common tasks with φ-evo.

## Build a network manually


Before even starting a simulation, let us build a network manually in
order to get familiar with the way they are encoded in the program. Most
of the code is written in python [^1], let us call our first file
**HowTo\_manualNetwork.py** is provided in the example directory.

``` python
# Import libraries
from phievo.Networks import mutation
import random

# Create a random generator and a network
seed = 20160225
g = random.Random(seed) # This define a new random number generator
L = mutation.Mutable_Network(g) # Create an empty network
```

We have created a first network, **L**, that can be used as a container for
the species and insteractions. For now **L** is still empty, we
can add a new species as follows

``` python
parameters=[['Degradable',0.5]] ## The species is degradable with a rate 0.5
parameters.append(['Input',0]) ## The species serves as an input referenced by the index 0 in the evolution algorithm.
parameters.append(['Complexable']) ## The species can be involved in a complex
parameters.append(['Kinase']) ## The specise can phosphorylate another species.
parameters.append(['TF',1]) ## 1 for activator 0 for repressor

## Create a species and add it to the network
## S1 is a reference to access quickly to the newly created species latter in the code
S1 = L.new_Species(parameters)
```

All the characteristics we want to associate with the species are listed
followed by their parameters. The list is then sent to the *new\_Species*
function to create the species. This methods is used when adding a external species (such as an input) that is not produced by the network itself.

In most cases a species comes with its transciptional machinery (*Species* + *CorePromoter* + *TModule*). The species and its related componant are added via the *new_gene* function.

Similarly a *PPI*(protein-protein interaction) is added with the complexation reaction and a phosphorylated species is added with the phosphorylation interaction.

Adding these functions to a code would look like this

``` python
parameters=[['Degradable',0.5]]
parameters.append(['TF',1])
parameters.append(['Complexable'])
TM0,prom0,S0 = L.new_gene(0.5,5,parameters)

parameters=[['Degradable',0.5]]
parameters.append(['TF',0])
parameters.append(['Complexable'])
TM1,prom1,S1 = L.new_gene(0.5,5,parameters)

parameters=[['Degradable',0.5]]
parameters.append(['TF',1])
parameters.append(['Phosphorylable'])
TM2,prom2,S2 = L.new_gene(0.5,5,parameters)


parameters=[['Degradable',0.5]]
parameters.append(['TF',0])
TM3,prom3,S3 = L.new_gene(0.5,5,parameters)

## Add complexation between S0 and S1.
parameters.append(['Kinase'])
ppi,S4 = L.new_PPI(S0 , S1 , 2.0 , 1.0 , parameters)

## Add a phosphorylation of S2 by S4
S5,phospho = L.new_Phosphorylation(S4,S2,2.0,0.5,1.0,3)
S5.change_type("TF",[1]) # Note this is already the default value for a phosphorilated species

## Regulate the production of S1 by S3 and S5
tfhill1 = L.new_TFHill( S3, 1, 0.5, TM1,activity=1)
tfhill2 = L.new_TFHill( S5, 1, 0.5, TM1,activity=1)
```

To display the layout of the former network, the program provides draw
function :

``` python
L.draw()
```

## Run a simulation

A φ-evo project is stored in a directory named as the project. 

``` bash
mkdir lac_operon
```

It contains all the configuration files of the project

-   initialization.py (name must start with "init"): Contains the initialization parameters, the path to the C files and optionally an inial network. If the former is not described in the initialyzation file, it will be generated randomly.

-   a fitness **C** file code used to compute the fitness. After an integration, the dynamics is stored in an array `history[SPECIES][TIME][CELL]`. You need to create a custom set function that analyse this array. In the end, the function *treatment\_fitness* should print the fitness of the network. 

- An init history file that contains the code that sets `history[SPECIES][t=0][CELL]` wrapped in a function called *init_history*.

- An init input file creates an *input* function. The input function is called at every time step to modify the `history` if necessary. 

### initialization.py

This file stores the informations about the evolution such as the ranges
of variation for the parameters, the mutation rates, the paths to the C
files, or the algorithm parameters.

The dictionary *dictionary\_ranges* sets the range of values a parameter
can take. If only one value Max is given, then the the range is
\[0,Max\]. To specify the the minimal value for a parameter, you have to
provide an array \[Min,Max\]

``` python
## The hill coefficient of a TFhill can varry between 1 and 5.
dictionary_ranges['TFHill.hill']= [1., 5.0]
## The rate of a TModule can varry between 0 and 2.
dictionary_ranges['TModule.rate']= 2
```

The dictionary *cfile* contains the path of the C files

``` python
cfile['fitness'] = fitness.c
cfile['init_history'] = init_history.c
cfile["inputc] = input.c
```

The dictionary *dictionary\_mutation* contains the rates at which a
mutation in the network appears. Note that the alorithm gathers the
rates provided and normalizes them in order to have an average of one
mutation per new generation during the evolution.

``` python
## Rate of appearance of the new transcription factor
dictionary_mutation['random_gene(\'TF\')']=0.02
```

The *prmt* dictionary contains the parameters related to the functioning
of the program and the algorithm.

``` python
## Number of integration step in the Euler integrator
prmt['nstep'] =3000
## time step during the integration
prmt['dt'] = 0.05
## Setting prmt['restart']['activated'] to False allows to start a fresh simulation
prmt['restart'] = {
  "activated": False,
  "freq": 50 # Generation frequency for saving the complete population
}
## Define the compiler (gcc by default)
prmt["compiler"] = "g++"

prmt['langevin_noise'] = 0 # Intensity of the langevin noise for stochastic simulation
prmt['multipro_level'] = 1 # Use multiprocess if one 1. If 0, singlethread.
## 
```


You may also specify the type of output you want and to prevent deleting species with a specific tag:
```python
list_unremovable=['Input','Output']
list_types_output=['TF']
```

We can choose an intial network to start the simulation with. This is
done through the *init\_network* function. The construction of the
initial network follows the steps presented in [Build a network manually](#build-a-network-manually).


### fitness.c

This file contains a C function *treatment\_fitness* used by the
algorithm to compute the fitnesses during the runs. After the integration, the algotithm reads the fitness(es) preinted by this function. You are free to add more analysis functions and to redefine *treatment_fitness* as long as it prints the network's fitness and has the following prototype:

``` c
void treatment_fitness(double history[NGENE][NSTEP][NCELLTOT], int trackout[])
    ...
    printf("%f",fitness)
```

The `trackout` lists the indexes of the outputs in the networks. You can also decide to use the global list `trackin` which contains the indexes of the ouputs.

### init_history.c

Before every integration, the algorithm reads the array
`history[NGENE][0][NCELLTOT]` to set the initial conditions of the run. You can use the *init_history.c* file to edit the first time step, this way it will be used as a initial condition.

Note that you can be more specific by using the two lists `trackin` and
`trackout` that contain the indexes for the inputs and outputs
respectively.

``` c
void init_history()  {
 int ncell,n_gene;
   for (ncell=0;ncell<NCELLTOT;ncell++){
     for (n_gene=0;n_gene<SIZE;n_gene++){
       history[n_gene][0][ncell]=0;
     }
   }
 }
```


### input.c

Sometime it is necessary to add artificial inputs during an integration. This is done via the *input* function. The *input* function is called at every time step and for every cell before computing the species derivatives. Since the derivatives for the species at time *t* are computed based on the values `history[NGENE][t][NCELLTOT]`, you can use *input* to modify the `history` array.

```c
void inputs(int time,int cell,int trial){
	...
}
```

To get more precise informations, we recommand you to have to look at how
*Examples/lac_operon/* project  is built.

### Launching a run

The program is launched with the *run_evolution.py* script:

``` bash
./run_evolution.py -m lac_operon/
```

The script loads the parameters and launches the run.
 
*run_evolution.py* should be placed in the same project directory as the project directory:
```bash
	|
     --- run_evolution.py
	 --- (Analyse Run.ipynb)
     --- example_project/
	              |
				   --- initialization.py
				   --- fitness.c
				   --- init_history.c
				   --- input.c
			
```

**Note:** *run_evolution.py* is not installed with phievo and must be downloaded manually from [here](https://raw.githubusercontent.com/phievo/phievo/master/run_evolution.py) or by running the command `phievo.download_tools()` in a python shell.

To restart a new run, one must provide the *#* of the run (or seed index). By default,
the run number is 0. To prevent errasing a run by mistake, the code will
not start if you do not provide a new run number in the initialization file. You can also tell the program explicitly to clear the Seeds with the "-c" or "--clear" option.

``` bash
./run_evolution.py -cm lac_operon/
```

## Restart an evolution

Every *k* generations, the algorithm saves a complete generation in a file called *Restart_file* in the Seed's directory. If interrupted, you can use this *Restart_file* to restart from a backup generation. You can set the restart generation in the initialization file:

```python
prmt['restart'] = {
  "activated": True, ## Activate restart
  "seed": 0, ## Index of the restart seed
  "kgeneration": 50, # Generation where to restart the algorithm
  "same_seed": True,
  "freq": 50 # Keep the same saving frequency
}
```
When the seed and the generation is not set or `None`, φ-evo will uses the last backup-ed generation in the seed with highest index.

## Pareto evolution

To start a pareto(multi-objectives) optimization with φ-evo, extra paremeters need to be defined in the initialization file:
```python
prmt['pareto']=True ## Activates pareto evolution
prmt['npareto_functions']=2 ## Number of fitness components
prmt['rshare']=0 ## Radius under which networks are penalysed for being too
                 ## close on the pareto front
```


[^1]: The front interface is coded in **python** (version &gt;3.5). But
    for efficiency reason, the core integration is coded in **C**.
