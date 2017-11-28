# Immune recognition

The project of immune recognition reproduces the evolution presented in Lalanne et al. 2013 [^1].

It is an advance example of custom network where the new interactions are added to $\varphi$-evo:

- Phosphorylation
- Dephosphorylation
- KPR_Binding
- KPR_Unbinding

## Run evolution

Running this project works as any other project. It is launch with **run_evolution.py** placed in the parent directory:

```python
>>> import phievo
>>> phievo.download_tools() 
>>> phievo.download_example("immune")
```
```bash
$ ./run_evolution.py -m example_immune
```

## Analyse the results

Immune comes with its custom analyse notebook stored in the project directory. Make sure you copy bthis notebook in the parent directory before you start the analysis.

```bash
$ cp  example_immune/Analyse_pMHC.ipynb .
```

You may also need to update the *sys.path* within the notebook first cell if you have decided to use an other name than  *example_immune* for the project.




[^1]: [Lalanne JB, François P. Principles of adaptive sorting revealed by in silico evolution. Physical Review Letters. 2013 May;110(21):218102.](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.110.218102)
