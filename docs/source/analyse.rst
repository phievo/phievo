Results and Analysis Tools
==========================

φ-evo has a module dedicated to the analysis of the results. The results
are stored in a *Simulation* object that contains a set of method that
give a quick access to the most relevant observable of a run. To start
analyzing the ``evo_dir`` project, you need to create a *Simulation*
object associated to it.

.. code:: python

    from phievo.AnalysisTools import Simulation

    sim = Simulation("evo_dir")

From there it is pretty straight forward to explore the architecture of
the results. A simulation contains Seeds which themselves contain
Network. In order not to overload the memory, the Seeds do not contain
exactly the networks but a link to them. As an example, here is how you
would load the best network for generation 326 in the seed number 2:

.. code:: python

    sim = Simulation("evo_dir")
    best_net_2_326 = sim.Seeds[2].get_best_net(326)
    # Equivalent to
    sim = Simulation("evo_dir")
    best_net_2_326 = sim.get_best_net(2,326)

Organization of the results
---------------------------

If you want to understand why the *Simulation* object is organized the
way it is and how to go beyond its possibilities, you will need to have
an idea of how φ-evo stores the results of a simulation.

By default, for every generation *g* only one *Network* is stored using
pickle in a file labelled ``Bests_*g*.net``. When the simulation has
only one fitness objective, this network is the one with the best
fitness in the population. However when the evolution is run using a
multiobjective criterium (like pareto optimisation), the best net is
chosen randomly amongst the network of rank 1.

The former storing method limits the disk space usage. However you might
want to store the whole population either for restarting the algorithm
from a given generation or to analyze every member of the generation. To
add this feature, you can specify a storing period by setting the
``prmt['restart']['freq']`` parameter in the initialization file before
launching the simulation. For example, if you set it to 50, the complete
population will be stored every 50 generations in a python *shelve*
named ``restart_file``.

Other files created:

-  ``data`` is a quick access shelve file to elementary informations
   stored as lists at the following keys:

   -  *generation*: index of the generation
   -  *fitness*: fitness of the best network
   -  *n\_species*: number of species in the best network
   -  *n\_interactions*: number of interaction in the best network

-  ``parameters`` is a copy of the parameter dictionnaries (defined for
   the non default in the initialization file) that were used during the
   simulation.
-  ``log_#.c`` Copy of the input, fitness, history, etc. **c** files
   used for the simulation.
