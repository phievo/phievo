HowTo
=====

This tutorial lists examples of how to perform common tasks with φ-evo.

## Build a network manually


Before even starting a simulation, let us build a network manually in
order to get familiar with the way they are encoded in the program. Most
of the code is written in python [^1], let us call our first file
**HowTo\_manualNetwork.py** (A jupyter notebook) is provided in the example
directory).

``` python
# Import libraries
from Networks import mutation
import random

# Create a random generator and a network
seed = 20160225
g = random.Random(seed) # This define a new random number generator
L = mutation.Mutable_Network(g) # Create an empty network
```

We have created a first network, **L**, that be used as a container for
the genes and the regulatiry elements. For now **L** is still empty, we
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

All characteristics we wnt to associate to the species are listed along
with their parameters. The list is then sent to the *new\_Species*
function to create the species.

This way of adding a species is **not recommanded**. For instance a gene
includes also a CorePromoter and and a TModule, adding a only a new
species makes it difficult to track the relation it has with the rest of
the network. Instead one would prefer using a dedicated functions that
will manage the synchronization between the different elements. Those
functions are

Adding these function to a code would look like this

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

To run a simulation, the first thing is to create a run directory where
to store the configuration files. The run directory is the place in
which the program stores and compiles the C files used by the run. It is
also the place where the the result a generated.

``` bash
mkdir lac_operon
```

The configuration files that must be added to the directory are

-   initialization.py

    Contains the initialyzation parameters, the path to the C files and
    optionally an inial network. If the former is not described in the
    initialyzation file, it will be generated randomly.

-   a fitness **C** file code used to compute the fitness. This file
    must contain a function *treatment\_fitness* that computes the
    network fitness(or calls other function to do so) and communicates
    the fitness to the rest of the program through a standard print. The
    path of this file is given in the *initialyzation.py* file.
-   init\_history.c The code stores th dynamics in a C array called
    *history*. This file contains a function, *init\_history*, in charge
    of setting the initial conditions before every run.

### initialyzation.py

This file stores the informations about the evolution such as the ranges
of variation for the parameters, the mutation rates, the paths to the C
files, or the algorithm parameters.

The dictionary *dictionary\_ranges* sets the range of values a parameter
can take. If only one value Max is given, then the the range is
\[0,Max\]. To specify the the minimal value for a parameter, you must
provide an array \[Min,Max\]

``` python
## The hill coefficient of a TFhill can varry between 1 and 5.
dictionary_ranges['TFHill.hill']= [1., 5.0]
## The rate of a TModule can varry between 0 and 2.
dictionary_ranges['TModule.rate']= 2
```

The dictionary *cfile* contains the path of the C files

``` python
cfile['fitness'] = lac_operon/fitness.c
cfile['init_history'] = lac_operon/init_history.c
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
  "freq": 50 # Backup for restart frequency(in generation)
}
```

We can choose an intial network to start the simulation with. This is
done through the *init\_network* function. The construction of the
initial network follows the steps presented in [Build a network manually](#build-a-network-manually).

### fitness.c

This file contains the function *treatment\_fitness* used by the
algorithm to compute the fitnesses during the runs. The file is written
in C. You are free to define this function as you whish as long as it
has the following prototype:

``` c
void treatment_fitness(double [NGENE][NSTEP][NCELLTOT], int trackout[])
    ...
    printf("%f",fitness)
```

### init\_history.c

Before every integration, the algorithm reads the array
*history\[NGENE\]\[NSTEP\]\[NCELLTOT\]* to init the run. You can use the
*init\_history.c* file to edit the first time step history, this way it
will be used as a initial condition.

Note that you can be more specific by using the two lists *trackin* and
*trackout* that contain the indexes for the input and output
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

To get more precise informations, we recommand you to have to look how
*Examples/lac\_operon/* project is built.

### Launching a run

The program is launched with the *run\_evolution.py*

``` bash
./run_evolution.py -m lac_operon/
```

The script loads the the parameters and launches the run. Along the run,
several files are kept: - The dynamics data are stored in the
*lac\_operon/Buffer\#* - a C file is generated and stored in
*lac\_operon/Workspace/* - The best network for each generation is
stored in *lac\_operon/Seed\#/*

To restart a new run, one must provide the *\#* of the run. By default,
the run number is 0. To prevent errasing a run by mistake, the code will
not start if you do not provide a new run number. You can also tell the
program explicitly to clear the Seeds with the "-c" or "--clear" option.

``` bash
./run_evolution.py -cm lac_operon/
```
## Restart an evolution

Every *k* generations, the algorithm saves a complete generation in file called *Restart_file* in the Seed's directory. If interrupted, you can use this *Restart_file* to restart from a backup generation. You can set the restart generation in the initialization file:

```python
prmt['restart'] = {
  "activated": True, ## Activate restart
  "seed": 0, ## Index of the seed
  "kgeneration": 50, # Generation where to restart the algorithm
  "same_seed": True,
  "freq": 50 # Keep the same saving frequency
}
```
When the seed and the generation are not set or `None`, φ-evo will search for the last backuped generation in the seed with highest index.

## Pareto evolution

To start a pareto optimization with φ-evo, extra paremeters need to be defined in the initialization file:
```python
prmt['pareto']=True ## Activates pareto evolution
prmt['npareto_functions']=2 ## Number of fitness components
prmt['rshare']=0 ## Radius under which networks are penalysed for being too
                 ## close on the pareto front
```

[^1]: The front interface is coded in **python** (version &gt;3.4). But
    for efficiency reason, the core integration is coded in **C**.
