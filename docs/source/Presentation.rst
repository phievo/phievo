Presentation
============

This section presents the basics

Network components
------------------

Species
~~~~~~~

Species is one of the two major components of a network. A species is a
protein that can have different functions. Adding a new type requires a
list containing the type name as a first element. If a type require
parameter, the must complete the list in a pre-defined order. For
instance a degradable species comes with its degradation rate.

Most of the former attributes are handled internally, but nothing
prevent us from adding a type manually:

.. code:: python

    mySpecies.add_type(["Degradation",0.5])
    print(mySpecies.dict_types())

Interaction
~~~~~~~~~~~

The Interactions, as suggested by its name, accounts for how species and
TModules interact. Examples of interactions are protein-protein
interactions, transcription factor regulations, etc.

TModule
~~~~~~~

A TModule is the last type of network component. It models for the
transcription part of a gene and are connected to the species they
control via a CorePromoter interaction. I φ-evo's framework, a gene is
represented by the three components *TModule*, *CorePromoter*, and
*Species*.

The transcription factor species that regulate a gene are connected to
the Tmodule.

.. figure:: TModule.svg
   :alt: 
   :figclass: align-center
   :width: 500px

Network
~~~~~~~

The network class is a container that includes the different Components
(species, TModules, and interactions) presented above. A network is
encoded using a biparpatite where species and TModules are connected to
interactions. The organisation of the *Network*'s graph relies on the
`networkx package <https://networkx.github.io/>`__.

An extra layer called *MutableNetwork* handles the mutations in the
*Network*

The *deriv2* class is responsible for reading a *Network*'s interactions
and to generate a C file with the species differencial equation used for
the integration.

Dynamical components
~~~~~~~~~~~~~~~~~~~~

To simulate the dynamics of a species the program first needs to explore
the nodes and the interactions that are connected and to build the
equations that govern the dynamic of the its concentration. The
equations are exported to c code and integrated.

The following examples presents networks components are converted into
ordinary differential equations.

TModule
^^^^^^^

There exists two types of TF actions: activition and inhibition. Both
types are modelled using Hill functions but there their effects is
included differently to the global regulation. Only the maximum of all
the activation is accounted for whereas the inhibitions are
multiplicative. In some extend activation and repression work
respectively as OR and NAND logic gates.

Next the CorePromoter interaction adds a delay :math:`\tau_P` to
accounts for the protein synthesis time. Practically, the algorithm
considers the state of the system at time :math:`t-\tau_P` to estimate
the production of :math:`P` at time :math:`t`.

The following configuration

.. figure:: TFHill_interaction.svg
   :alt: 
   :figclass: align-center
   :width: 500px

leads to the equation

.. math:: \frac{d S}{d t} = \left(\max\left\{PR_S \times\max\left\{\frac{A_1^{n_{A1}}}{A_1^{n_{A1}} + h_{A1}^{n_{A1}}}, \frac{A_2^{n_{A2}}}{A_2^{n_{A2}} + h_{A2}^{n_{A2}}}, \ldots \right\},B_S\right \}\times \frac{h_{R1}^{n_{R1}}}{R_1^{n_{R1}} + h_{R1}^{n_{R1}}} \times \ldots \right)_{(t-d_S)}

\_\_

In the above equation, the :math:`h` and :math:`n` parameters correspond
respectively to the hill saturation and exponent. The :math:`PR` is the
production rate of the protein in optimal conditions and :math:`B` is
the basal rate(in case no activator is present). The overall production
is modulated by the repression.

Degradation
^^^^^^^^^^^

Every protein :math:`P` labelled as *degradable* is degraded over time
with a rate :math:`\delta_P`. This

.. math:: \frac{d P}{d t} =  - \delta_P P

Phosphorylation
^^^^^^^^^^^^^^^

The phosphorilasion is the addition of a phosphate group to a Species by
a kinase. It creates a new phophorilated species. The dynamics of this
mechanism is controlled by a hill function that accounts for the use of
the kinase by all the different species. In the case of of kinase that
catalyses the phosphorilation of two species :math:`S_1` and
:math:`S_2`.

.. math:: \frac{d S_1}{dt} = - \frac{d S_1^{*}}{dt} = k_p^1\frac{K \left(\frac{S_1}{h_1}\right)^{n_1}}{1+\left(\frac{S_1}{h_1}\right)^{n_1} + \left(\frac{S_2}{h_2}\right)^{n_2}} - k_d^1 S_1^{*}

.. math:: \frac{d S_2}{dt} = - \frac{d S_2^{*}}{dt} = k_p^2\frac{K \left(\frac{S_2}{h_2}\right)^{n_2}}{1+\left(\frac{S_1}{h_1}\right)^{n_1} + \left(\frac{S_2}{h_2}\right)^{n_2}} - k_d^2 S_2^{*}

.. figure:: Phospho_interaction.svg
   :alt: 
   :figclass: align-center
   :width: 300px

Protein-Protein-Interaction (PPI)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The PPI interaction accounts for the complexation of two single proteins
into one complex.

.. figure:: PPI_interaction.svg
   :alt: 
   :figclass: align-center
   :width: 300px

The rate is obtained from a mass-action dynamics:

.. math:: \frac{d P_1}{dt} = \frac{d P_2}{dt} = - \frac{d C}{dt} = - \text{rate} = - k^{+}P_1P_2 + k^{-} C

with :math:`k^{+}` and :math:`k^{-}` being respectively the forward and
backward rate constants

.. raw:: html

   <!-- #### Ligand-Receptor interaction (LR) -->

.. raw:: html

   <!-- This interaction corresponds to the complexation of two species - a -->

.. raw:: html

   <!-- ligand and a receptor - to trigger a response in the system. -->

.. raw:: html

   <!-- ![](LR_interaction.svg){.align-center width="300px"} -->

.. raw:: html

   <!-- The ligand concentration are assumed to be add steady state which allows -->

.. raw:: html

   <!-- to describe the rate using the *Michaelis-Menten-Henri* formalism: -->

.. raw:: html

   <!-- $$\frac{d L}{dt} = \frac{d R}{dt} = - \frac{d C}{dt} = - \text{rate} = - \frac{V\,L\,R}{h + R}$$ -->

.. raw:: html

   <!-- with $V$ and $h$ being respectively the association rate and the -->

.. raw:: html

   <!-- association threshold. -->

.. raw:: html

   <!-- ## Evolution -->

.. raw:: html

   <!-- The evolution algorithm mimics Darwinian selection. It generates an initial population (of constant size size defined by the user) where the individuals are in competition to pass their genome to the next generation. Only the fittest half of the individuals passes to next generation and is allowed do reproduce (by duplication) in order to maintain the population size. -->

Pareto evolution
----------------

In the case where the fitness is composed of multiple components, it is
not obvious how to balance the different modules in the global fitness.
It may be interesting to have a multiple objective optimization where
all the components have the same importance; only changes improving a
component without decreasing the others are kept. The fitness
:math:`F = \{f_1,f_2,...,f_N\}` is of higher rank than
:math:`G = \{g_1,g_2,...,g_N\}` if

.. math:: \forall i\quad f_i\geq g_i

.. math:: \exists k,\quad f_k>g_k

Clearly multiple objective optimisation does not result in one best
network in the end but to a population of highest rank networks called
the Pareto front. More information can be found on
`Wikipedia <https://en.wikipedia.org/wiki/Multi-objective_optimization>`__.

From a practical standpoint, the algorithm works similarly to the
genetic algorithm with a modified selection process. As in the genetic
algorithm, half of the population is passed to the next generation and
duplicated. Because the only classification criterion is the network's
rank, the cutoff may occur in the middle of a set of equivalent network
since they have the same rank. In such a case the algorithm selects
randomly the networks with the cutoff rank to complete the set of
individuals passed to the next generation.
