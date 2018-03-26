Simulation parameters
=====================

This sections presents the different parameters that can be set in the
`initialization file <create_new_project.html#initialization-py>`__.

Kinetic parameters (``dictionary_ranges``)
------------------------------------------

The kinetic parameters are specific to a type of interaction or a type
of species. They are stored in the ``dictionary_ranges`` dictionary. One
can define the range over which they can vary by setting the its range
with a size 2 list (if the minimum is 0, the range can be set with a
float corresponding to the maximum).

In the example of a ``Degradation`` interaction, the range of variation
of the ``rate`` of degration is set with

.. code:: python

    dictionary_ranges["Degradation.rate"] = 0.0

Mutation parameters (``dictionary_mutation``)
---------------------------------------------

The mutation parameters define the rate a which a given mutation occurs.
Note that the evolution rescale the generation time so that a network
undergoes an average of one mutation per generation. The mutation
parameters are defined in the ``dictionary_mutation`` dictionary.

A new mutation function is defined when creating a `new
interaction <new_interaction.html>`__. Each new mutation can have its
own rate defined.

Examples:

-  ``dictionary_mutation["random_gene()"]``: Rate at which
   ``random_gene()`` mutation is executed (with default settings).
-  ``dictionary_mutation["random_gene(type = 'TF')"]``: Rate at which
   ``random_gene()`` mutation is executed with the parameter *type*
   equal to *"TF"* (it creates a species with a tag "TF" corresponding
   to a transcription factor).

General simulation parameters (``prmt``)
----------------------------------------

The general simulation parameters are stored in a dictionary called
``prmt``:

-  Number of seeds (``nseed``): Number of independent evolution to
   simulate.
-  First seed (``firstseed``): Index of the first seed. This index is
   also used to seed the random number generator.
-  Number of generations (``ngeneration``): Number of generation to
   simulate in each independent evolution.
-  Number of cells (``ncelltot``): Number of cells in the organism.
-  Population size (``npopulation``): Number of network in the
   population.
-  Number of neighbors (``nneighbor``): Number of neighbors cell has.
-  Fraction mutated per gen (``frac_mutate``): Fraction of networks in
   the population to mutate at every generation.
-  Number of Inputs (``ninput``): Number of species with an inputs (with
   an input tag) a network should have.
-  Number of Outputs (``noutput``): Number of species with an outputs
   (with an output tag) a network should have.
-  Number of trials (``ntries``): When a fitness depends of a network's
   initial conditions or in the presence of Langevin's noise, it is
   useful to run several independent kinetic integrations. ``ntries``
   determines the number of integrations to run. Note that in the case
   the initialization of each trial should be done with the
   `init\_history <create_new_project.html#init-history-c>`__ function
   and the agregation(e.g. averaging) of the fitnesses coresponding to
   each integration is done with the
   `treatment\_fitness <create_new_project.html#fitness-c>`__ function.
-  Time step dt (``dt``): Size of an integration time step in the Euler
   algorithm.
-  Number of time steps (``nstep``): Number of integration time step in
   the Euler algorithm.
-  Langevin noise value (``langevin_noise``): Level of the langevin
   noise in a stochastic simulation. When 0, the integrations are
   deterministic.
-  Gillespie generation time (``tgeneration``): The computation of the
   next mutation follows a Gillespie algorithm. ``tgeneration`` defines
   the initial time, then the time ``tgeneration`` is updated to have
   roughly one mutation in ``frac_mutate`` of the networks.
-  Recompute networks (``redo``): Should the networks that do not change
   in from a generation to the other be re-integrated ti compute the
   fitness?
-  Pareto simulation (``pareto``): Should we run a Pareto integration?
-  Number of pareto functions (``npareto_functions``): Number of pareto
   functions defined?
-  Pareto penalty radius (``rshare``): This parameter prevents a network
   from being dominated by a networks with fitnesses that fall too close
   to it current position in the fitness space. Increasing ``rshare``
   helps to explore a larger portion of the fitness space. `Warmflash et
   al
   2012 <http://iopscience.iop.org/article/10.1088/1478-3975/9/5/056001/meta>`__.
-  Multiple threads (``multipro_level``): Should the algorithm run in
   parallel?
-  Generation printing frequency (``freq_stat``): During a simulation
   the algorithm regularly prints informations about its current state.
   ``freq_stat`` defines the number of generations between two prints.

Restart parameters (``prmt["restart"]``)
----------------------------------------

To restart a simulation either after it has been stopped or from a
specific seed and generation one can configure the *restart* parameters.
The parameters are hosted in a sub-dictionary or ``prmt``,
``prmt["restart"]``:

-  ``prmt["restart"]["activated"]``: Activate restart
-  ``prmt["restart"]["freq"]``: Frequency at which a complete generation
   is saved.
-  ``prmt["restart"]["kgeneration"]``: Generation at which to restart
   the algorithm
-  ``prmt["restart"]["seed"]``: Seed at which to restart the algorithm
-  ``prmt["restart"]["same_seed"]``: Restart with the same seed

More information is available on the (restart an evolution
section)[create\_new\_project.html#restart-an-evolution].
