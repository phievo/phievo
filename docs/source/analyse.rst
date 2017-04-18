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
pickle in a file labelled ``Bests_g.net``. When the simulation has only
one fitness objective, this network is the one with the best fitness in
the population. However when the evolution is run using a multiobjective
criterium (like pareto optimisation), the best net is chosen randomly
amongst the network of rank 1.

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

Analysis Tools
--------------

In this section we will explore the built in functions that are bound to
a *Simulation* object.

custom\_plot
~~~~~~~~~~~~

Plots the seed two observables one against each other. The available
observables are the ones present in the ``data`` file ("generation",
"fitness", "n\_species", "n\_interactions").

.. code:: python

    sim.seeds[1].custom_plot("generation","fitness")
    # Similarly can use the shortcut
    sim.custom_plot(1,"generation","fitness")

plot\_fitness
~~~~~~~~~~~~~

There also exists a method to plot the fitness directly:

.. code:: python

    sim.seeds[1].show_fitness()
    # Similarly can use the shortcut
    sim.show_fitness(1)

get\_best\_net
~~~~~~~~~~~~~~

Get the best net found in a given generation (the function reads the
``Bests_g.net`` file and return the Network object)

.. code:: python

    bestnet_g5_seed3 = sim.seeds[3].get_best_net(5)
    # or
    bestnet_g5_seed3 = sim.get_best_net(3,5)

get\_backup\_net
~~~~~~~~~~~~~~~~

If you want to extract a network from a entirely stored generation, you
can use *get\_backup\_net*. Be careful though, not every population is
stored in the ``restart_file``.

.. code:: python

    net8_g50_seed3 = sim.seeds[3].get_backup_net(50,8)
    # Or
    net8_g50_seed3 = sim.get_backup_net(3,50,8)

stored\_generation\_indexes
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The *stored\_generation\_indexes* is a reminder of which generations are
stored.

.. code:: python

    lost_stored = sim.seeds[1].stored_generation_indexes()
    # Or
    lost_stored = sim.stored_generation_indexes(1)

Running a network's dynamics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By construction φ-evo does not allow yet allow to quickly run the
dynamics of a network. Namely a Network object has no method that
directly returns the derivative from at a given state. Instead φ-evo has
a method to write a **c** file containing the derivative function and
that runs the dynamics on pre-defined inputs. This may seem a bit bulky
but the software was initially written to evaluate the fitness of a
given network and that is better done in **c**.

However the *Simulation* has the method *run\_dynamics* to ease the
access to the results the dynamics.

.. code:: python

    net = sim.get_best_net(3,5)
    dyn_buffer = sim.run_dynamics(net=net,trial=1)

This runs the dynamics that would be run in the evolution algorithm with
the history and input **c** files you provided in the project directory.
You can specify the number of trial you want to run if the dynamics is
stochastic. The buffer returned by the function is dictionary where the
main the "time" and "net" keys give you access to respectively the time
vector and the network used for the run. The other keys are the index of
the trial for which you want to access the data. Note that the buffer is
also stored in the *Simulation* object as *buffer\_data*.

Plotting the results of a dynamics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simulation object allows you to plot the two most obvious result you
would like to see after running a dynamics:

1) The time course of the genes in a given cell with *Plot\_TimeCourse*
2) The evolution of the genes along the system at a given time point
   with *Plot\_Profile*

.. code:: python

    sim.Plot_TimeCourse(trial_index=1,cell=1)
    sim.Plot_Profile(trial_index=1,time=1)

Draw a network's layout
~~~~~~~~~~~~~~~~~~~~~~~

The network object contains a function to draw the layout of its gene
interactions:

.. code:: python

    net = sim.get_best_net(3,5)
    net.draw()

Notebook
--------

To facilitate the use of the former functions, φ-evo as a class
*Notebook* that is used to run them in a `jupyter
notebook <https://jupyter.org>`__.

All the functions described previously can be used directly in a jupyter
notebook but the *Notebook* class increases the usability by handling
the dependencies between widgets. For instance you want the module in
charge of plotting a network's layout to be disabled as long as a Seed
and a Network have not been selected.

A Notebook object serves as a container for all the available modules
you can use in the jupyter notebook. A module contains the material to
handle a cell: its widgets, some update functions and a display function
that displays the widgets in the jupyter notebook. In the end, the user
only needs to run ``myNotebook.myModule.display()`` to create a jupyter
elementary app in a cell. Then the module should be able to handle the
expected inputs from the user.

Creating a custom module
~~~~~~~~~~~~~~~~~~~~~~~~

Every module of contained int the *Notebook* object of the *CellModule*
class. The latter is only a minimal template used to constrain the
minimal requirement a module must have.

-  ``__init__(self,Notebook)`` : The init function take the *Notebook*
   it is contained in as an argument.
-  ``display``: The function must be redefined to display the widgets
   and to handle the relation between them.
-  ``update`` : If the module has dependencies, this function must be
   defined. When dependency is updated, this function is called.

\_\_init\_\_
^^^^^^^^^^^^
