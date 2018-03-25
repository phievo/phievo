# Install φ-evo

φ-evo relies on python>=3.5, pip, and c.

The software has been successfully tested on the three main operating systems(windows,mac OSX, and GNU-linux) but **we recommend using a GNU-linux distribution(ubuntu)** as it has been tested more thoroughly and more regularly on this platform.


### install Anaconda
The _phievo package depends on python>=3.5.
If python is not already installed on your computer, we recommend to install it by using the [anaconda distribution](https://www.continuum.io/downloads).

Among other things, anaconda provides the standard package manager of python _pip_. Before anything, it is good to check that you are working with the most recent version of pip:

```bash
pip install --upgrade pip
```


**Note:** When multiple versions of python are installed on the same computer, you may need to specify the version of python or pip you are using: `python` (`pip`) for python2 and `python3` (`pip3`) for python3. Make sure that the  `which pip` and `which python` return the right pip (and python) installation path. For simplicity we will use `pip` in the following instructions.

**Note:** If you install packages for all the users of your computer, you need to have admidistrator rights and use `sudo` before the pip command. It can happens that your global and your local pip are not the same. To make sure the administrator uses the right pip, run `sudo which pip`. The installation instructions assume you do not need to add `sudo` before `pip`.

### install the package

With pip installed, the installation is straight forward, run:
```bash
pip install https://github.com/phievo/phievo/blob/master/dist/phievo-1.0.zip?raw=true
```

### Install gcc on windows
Windows does not come with the `gcc` compiler installed but the free software foundation provides a minimal distribution of the gnu softwares for windows, it is called [MinGW](http://mingw.org/).

Once you have downloaded `mingw-get-setup.exe`, run it. A selection panel will open. We recommend you to install at least the two following packages(the others are not relevant for φ-evo):
- mingw-developper-toolkit
- mingw32-base

Choose the default directory.

After the installation is finished, update windows `PATH` so that it knows where to look for the `gcc` command. Open a the command prompt and run:

```bash
setx PATH "%path%;C:\MinGW\bin"
```

**Note:** gcc is distributed by other packages such as code blocks or visual basics. In such case, you do not need to install MinGw. Just upload you `PATH` so that windows knows where is the gcc compiler.

### Install gcc on mac osx
OSX does not have the gcc compiler installed by default either. There are different ways to install it. The fastest is probably via [homebrew](https://brew.sh/):

```bash
brew install gcc
```

If `gcc` is not already installed on you system (via macports or Xcode), _homebrew_\'s _gcc_ should be automatically in the system's `PATH`.


### Install pygraphviz

_pygraphviz_ is not included in the default dependencies of _phievo_ because it does not exist natively on windows and we wanted to publish a version that that runs on all the systems. _pygraphviz_ is used only to display network layouts. If it is not installed, _phievo_ will print a warning and use _networkx_ spring layout instead.

On max OSX, you have to use homebrew to install graphvix first :

```bash
brew install graphviz pkg-config
pip install pygraphviz
```

On GNU/linux, installing the dependencies varies depanding on the distribution. We tested the following on debian and ubuntu

```bash
sudo apt-get install graphviz graphviz-dev pkg-config
sudo pip install pygraphviz
```
On other distributions, you want to find the equivalent of _graphviz_, _graphviz-dev_, and _pkg-config_.

We found that sometimes on ubuntu the C linking to the graphviz library does not work properly. The fix is to  be more explicit on the linking for the pip command:

```bash
sudo pip install pygraphviz --install-option="--include-path=/usr/include/graphviz" --install-option="--library-path=/usr/lib/graphviz/"
```

### run_evolution.py script

An extra script ([run_evolution.py](https://raw.githubusercontent.com/phievo/phievo/master/run_evolution.py)) needs to be downloaded with the phievo package to start an evolution. It is stored in the root of the phievo repository.

You can either manually download it or open a python terminal and run
```python
>>> import phievo
>>> phievo.download_tools()
```

The former utility also downloads a jupyter notebook that can be used to analyse the results of a simulation in current directory.

### Analyse notebook

We provide a [jupyter notebook](https://github.com/phievo/phievo/blob/master/Analyse%20Run.ipynb) at the root of the [github repository](https://github.com/phievo/phievo) to help with the analysis of the runs. If you wand to run it, you will need to install several extra python libraries, to help with this, they are writen in [extra.txt](https://raw.githubusercontent.com/phievo/phievo/master/extra.txt).
```bash
pip install -r https://raw.githubusercontent.com/phievo/phievo/master/extra.txt
```

Similarly to the ([run_evolution.py](https://raw.githubusercontent.com/phievo/phievo/master/run_evolution.py)) script, Analyse Run.ipynb is downloaded when you call the `phievo.download_tools()` function.

The jupyter kernel is started with the following command

```bash
jupyter notebook
```
Usually it autmotically opens a new windows in your terminal in which you need to select `Analyse Run.ipynb`. If the windows does not open, it can be open manually by copy-pasting the url printed in your shell after you ran the command in a wer browser.


When using the plotly package, you may find that the plots do dot display well in the notebook (white square), the solution to this problem is to increase the io rate allocated to the notebook by using the `NotebookApp.iopub_data_rate_limit` option when starting jupyter:

```bash
jupyter notebook --NotebookApp.iopub_data_rate_limit=10000000000
```

### Test your installation

To test that everything works properly, we recommend that you run an example simulation. Several examples of simulations are stored in the [github repository](https://github.com/phievo/phievo/tree/master/Examples) Examples directory. You can download all the simulations by cloning the repository with git:

```bash
git clone https://github.com/phievo/phievo.git
```

This will also download phievo's code.

To download a single example there is a built-in tool that can be run in a python shell:

```python
>>> import phievo
# Downloads run_evolution.py and Analyse Run.ipynb in  the current directory
>>> phievo.download_tools() 
# Downloads an example project directory
>>> phievo.download_example("adaptation") 
```

The function `download_example` allows to download one of the following examples:

- adaptation
- somite
- hox
- hox_pareto
- lac_operon
- immune
- seed_adaptation
- seed_adaptation_pruning
- seed_somite
- seed_somite_pruning
- seed_lacOperon
- seed_lacOperon_pruning
- seed_hox_pareto_light

The examples starting with "seed_" keyword also contain the results of the simulations. The results can directly be visualized in the Analyse notebook.

After downloading an example project directory and the *run_evolution.py* script you are all set to start an evolution.

```bash
|-- run_evolution.py
|-- Analyse Run.ipynb
`-- example_adaptation/
	|-- initialization.py
	|-- fitness.c
	|-- init_history.py
	`-- input.c
```


To launch the evolution, simply run

```bash
python run_evolution.py -m example_adaptation
```
**Note:**  You can add the -c option (`./run_evolution.py -cm example_adaptation`) to delete a Seed that was created by a former run and prevents a new run to start. Be careful, a deleted seed cannot be recovered.

If everything works correctly you should see the evolution starting. When an evolution is running it displays regularly updates of its current state in the terminal and a `STOP.txt` file is created at the root of the project. The purpose of the STOP file is to have a quick method to check on the current state of a run when it is launched as a background task. When the *STOP* file is deleted, the run stops.

### Create a new project

To start a new project, the best is to use an existing example as a template and to modify the relevant parameters.

Similarly to the `Analyse notebook`, we also propose the [Project Creator.ipynb](https://github.com/phievo/phievo/blob/master/Project%20Creator.ipynb) notebook to help with the creation of a new project.

```bash
jupyter notebook Project\ Creator.ipynb

```
