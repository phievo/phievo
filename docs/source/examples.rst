Examples
--------

φ-evo provides a series of examples of project and already run seeds.

Examples of projects
~~~~~~~~~~~~~~~~~~~~

The example of projects are stored in the
`Example <https://github.com/phievo/phievo/tree/master/Examples>`__
directory of the *phievo* package:

-  adaptation  [1]_
-  somites  [2]_
-  StaticHox and StaticHox pareto [3]_
-  `lac\_operon <example-lac-operon.html>`__

Examples of seeds
~~~~~~~~~~~~~~~~~

Because some simulation can take some time to run, we provide the result
seeds we used to generate the figure of the paper:

-  seed\_adaptation
-  seed\_adaptation\_pruning
-  seed\_somite
-  seed\_somite\_pruning
-  seed\_lacOperon
-  seed\_lacOperon\_pruning
-  seed\_hox\_pareto\_light

To download the result of a simulation on your computer, you can use
phievo:

.. code:: python

    import phievo
    phievo.download_example("adaptation")
    # To download the seed also
    phievo.download_example("seed_adaptation")

Hox pareto
~~~~~~~~~~

The complete simulation for the Hox Genes takes a lot of space, only a
portion of the original results is accessible through phievo.

You can manually download the complete simulation
`here <https://mcgill-my.sharepoint.com/personal/adrien_henry_mail_mcgill_ca/_layouts/15/guestaccess.aspx?docid=0f1beb049ce8d4a648261a691f3116cd3&authkey=AUsBUDDWzFpkWDjGIo6n5X4>`__.

References
~~~~~~~~~~

.. [1]
   `François P, Siggia ED. A case study of evolutionary computation of
   biochemical adaptation. Physical Biology.
   2008;5(2):26009. <http://iopscience.iop.org/article/10.1088/1478-3975/5/2/026009/meta;jsessionid=63E2805FAE2CE62F041C2DE212DDB0C1.ip-10-40-1-105>`__

.. [2]
   `François P, Hakim V, Siggia ED. Deriving structure from evolution:
   metazoan segmentation. Molecular Systems Biology. 2007
   Dec;3:9. <http://msb.embopress.org/content/3/1/154.long>`__

.. [3]
   `François P, Siggia ED. Predicting embryonic patterning using mutual
   entropy fitness and in silico evolution. Development (Cambridge,
   England).
   2010;137(14):2385–2395. <http://dev.biologists.org/content/137/14/2385>`__
