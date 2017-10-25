Presentation
============

This section presents the basics element to understand the structure of the algorithm and
the role of the various python modules.

## An algorithm overview

There is three main blocs in the algorithm that correspond to the three libraries present
at the root of the project.

- *Networks*: gathers all the elements to represent, modify and simulate the evolving
biological networks that are the indidividual level of our population.
- *Population_types*: implements the so-called genetic algorithm and its several variants
through a *Population* class and its subclasses.
- *AnalysisTools*: gathers the tools used after the simulation to analyse, represent and
study the results.

## Network components

*Network* (and its sub-class *Mutable_Network*) represents the individual level of our
evolutionary algorithm. Apart of the methods used to implement the different operations,
the main attribute is *graph*, a networkx.MultiDiGraph object that stores the biochemical
network as a bipartite graph of *Species* (and *TModule*) on one side and *Interaction* on
the other. The organisation of the graph thus relies on the [networkx
package](https://networkx.github.io/).

The subclass called *MutableNetwork* handles the mutations in the *Network*

The *deriv2* module is responsible for reading a *Network*'s  interactions and to generate
a C file that will integrates the differential equations simulating the species and then
computes the fitness of the network that will be used at the genetic algorithm level.

### Species

Species is one of the two major components of a network. A species is a protein that can
have different types (Degradable, Phosphorylable, etc.). Most of the time, those species
will be added automatically by the algorithm when handling the various interactions. For
example when two species are chosen to be part of a new protein-protein interaction, a new
species will be added to simulate the complex thus formed.

However, to manually build the initial network, you may want to add species with some
fixed properties. For this, you need to build a list of lists containing the type name as
a first element and the different parameters (if any) must complete the list in a
pre-defined order (```l_types``` in the example below). For instance a degradable species
comes with its degradation rate. Note that adding a single Species is actually quite rare
as they often came as a whole gene with a CorePromoter and a TModule (see the TModule
picture below).

``` python
l_types = [["Degradation",0.5],["Complexable"],["Output",0]]
mySpecies = my_Network.new_Species(l_types)
mySpecies = my_Network.new_gene(rate, delay, l_types, basal_rate)
```

But see the ```initiation.py``` file of an example to a complete construction of a Network
object.

### Interaction

The Interactions, as suggested by its name, accounts for how species and TModules
interact. Examples of interactions are protein-protein interactions, transcription factor
regulations, etc. See the sections below for various type of preimplemented interactions.

Note also that it is often necessary to implement new interactions tailored for a specific
task. (See ```Examples/immune/``` for an example of such new interactions.)

### TModule

A TModule is the last type of network component. It corresponds to the transcription part of a
gene and is connected to the species it controls via a CorePromoter interaction. In the
φ-evo's framework, a gene is thus represented by three components: *TModule*,
*CorePromoter*, and *Species*.

![](TModule.svg){.align-center width="500px"}

Uphill, the transcription factor species that regulate a gene are connected to the Tmodule
through TFHill interactions.

## Population & Evolution

The evolution algorithm mimics Darwinian selection by simulating a population where the
individuals are in competition to pass their genome to the next generation.

It first generates an initial population the size of which is defined by the user by
cloning an initial Network and from then follow cycles of mutation, fitness computation
and selection. Each cycle thus defines our time step of evolution and will subsequently be
called a generation.

### Elite strategy

By default we choose to use the *elite strategy* because of its robustness and its cheap
computationnal cost. Thus, during the selection step, the worst part of the population is
deleted, while the fittest half of the individuals are directly passed to the next
generation. Then, each of them is copied and this copy is mutated.

Note that this scheme automatically keeps the population size constant. Moreover, it
relies only on the rank of the individuals in the population and not on the quantitative
fitness. This makes it very robust to the possible difficulties and failures of the
fitness implementation.

### Pareto evolution

In the case where the fitness is composed of multiple components, it is not obvious how to
balance the different modules in the global fitness. It may be interesting to have a
multiple objective optimization where all the components of the fitness have the same
importance; only changes improving one component without decreasing the others are
considered as an improvement.

For a fitness splited in $N$ components: $F = \{f_1,f_2,...,f_N\}$. We say that individual
$i$ dominates (strictly) $j$ if and only if the fitnesses $F^i$ and $F^j$ are such that:
$$\forall k\quad f^i_k\geq f^j_k,\quad (\exists k \quad f^i_k>f^j_k)$$

Clearly multiple objective optimisation does not result in one best network in the end but
to a population of highest rank networks called the Pareto front. More information can be
found on [Wikipedia](https://en.wikipedia.org/wiki/Multi-objective_optimization).

From a practical standpoint, the algorithm works similarly to the genetic algorithm with a
modified selection process. As in the genetic algorithm, half of the population is passed
to the next generation and duplicated. Because the only classification criterion is the
network's rank, the cutoff may occur in the middle of a set of equivalent network since
they have the same rank. In such a case the algorithm selects randomly the networks with
the cutoff rank to complete the set of individuals passed to the next generation.

### Results

During the evolution, the results are stored in separate folder for each seed soberly called
_Seed*/_, this folder contains three main type of elements:

- *log_?* — are brute copy of the files used as input for this seed (the correspondance should be obvious).
- *Bests_?.net* — is a pickle of the *Network* object with the best fitness at the corresponding generation, this allows you to trace back the evolution of the individuals in the population
- *data.?* — contains various data about the seed (mean fitness, times, etc.)
- *Restart_file.?* — this shelve object contain a copy of the whole population in case you want
to restart the evolution after the termination of the first run of the program.

## Modelization & Integration

To simulate the dynamics of a species the program first needs to explore the nodes and the
interactions that are connected to it. Then it builds the equations that govern the
dynamic of its concentration. These equations are then written as C code and integrated.

The following sections presents the predefined networks interactions and there
corresponding ordinary differential equations.

### TModule and gene production

There exists two types of TF actions: activition and inhibition. Both types are modelled
using Hill functions but there their effects is included differently to the global
regulation. Only the maximum of all the activations is accounted for whereas the
inhibitions are multiplicative. In some extend activation and repression work respectively
as OR and NOR logical operations.

Next the CorePromoter interaction adds a delay $\tau_P$ to accounts for the protein
synthesis time. Practically, the algorithm considers the state of the system at time
$t-\tau_P$ to estimate the production of $P$ at time $t$.

The following configuration

![](TFHill_interaction.svg){.align-center width="500px"}

leads to the following equation:

$$\frac{d S}{d t} = \left(\max\left\{PR_S \times\max\left\{\frac{A_1^{n_{A1}}}{A_1^{n_{A1}} + h_{A1}^{n_{A1}}}, \frac{A_2^{n_{A2}}}{A_2^{n_{A2}} + h_{A2}^{n_{A2}}}, \ldots \right\},B_S\right \}\times \frac{h_{R1}^{n_{R1}}}{R_1^{n_{R1}} + h_{R1}^{n_{R1}}} \times \ldots \right)_{(t-d_S)}$$

__

In the above equation, the $h$ and $n$ parameters correspond
respectively to the Hill constant and coefficient. The $PR$ is the
production rate of the protein in optimal conditions and $B$ is the
basal rate(in case no activator is present). The overall production is
modulated by the repression.

### Degradation

Every protein $P$ labelled as *degradable* is degraded over time with a
rate $\delta_P$. This

$$\frac{d P}{d t} =  - \delta_P P$$

### Phosphorylation

The phosphorylation is the addition of a phosphate group to a Species by
a kinase. It creates a new phophorylated species. The dynamics of this
mechanism is controlled by a hill function that accounts for the use of
the kinase by all the different species. In the case of of kinase that
catalyses the phosphorilation of two species $S_1$ and $S_2$.

$$\frac{d S_1}{dt} = - \frac{d S_1^{*}}{dt} = k_p^1\frac{K \left(\frac{S_1}{h_1}\right)^{n_1}}{1+\left(\frac{S_1}{h_1}\right)^{n_1} + \left(\frac{S_2}{h_2}\right)^{n_2}} - k_d^1 S_1^{*}$$

$$\frac{d S_2}{dt} = - \frac{d S_2^{*}}{dt} = k_p^2\frac{K \left(\frac{S_2}{h_2}\right)^{n_2}}{1+\left(\frac{S_1}{h_1}\right)^{n_1} + \left(\frac{S_2}{h_2}\right)^{n_2}} - k_d^2 S_2^{*}$$

![](Phospho_interaction.svg){.align-center width="300px"}

Note that by default, there is no mechanism implemented for active dephosphorylation so that they hapen with constant rates $k_d^1$ and $k_d^2$.

### Protein-Protein-Interaction (PPI)

The PPI interaction accounts for the complexation of two single proteins
into one complex.

![](PPI_interaction.svg){.align-center width="300px"}

The rate is obtained from a mass-action dynamics:

$$\frac{d P_1}{dt} = \frac{d P_2}{dt} = - \frac{d C}{dt} = - \text{rate} = - k^{+}P_1P_2 + k^{-} C$$

with $k^{+}$ and $k^{-}$ being respectively the forward and backward
rate constants
