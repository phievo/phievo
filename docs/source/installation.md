# Install φ-evo

### install Anaconda
The _phievo package depends on python>=3.5.
In it is not already installed on your computer, we recommand to install it by using the [anaconda distribution](https://www.continuum.io/downloads).

Among other things, anaconda provides the standard package manager of python _pip_. Before anything, it is good to check that you are working with the most recent version of pip:

```bash
	pip3 install --upgrade pip
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
On other distribution, you may want to find the equivalent of _graphviz_, _graphviz-dev_, and _pkg-config_.

We found that sometimes on ubuntu the C linking to the graphviz library does not work properly, to fix this, be more explicit and use the linking for the pip command:

```bash
sudo pip install pygraphviz --install-option="--include-path=/usr/include/graphviz" --install-option="--library-path=/usr/lib/graphviz/"
```

### Analyse notebook

We provide a [jupyter notebook](https://github.com/phievo/phievo/blob/master/Analyse%20Run.ipynb) to help with the analysis of the runs. If you wand to run it, you will need to install several extra python libraries, to help with this, they are writen in [extra.txt](extra.txt).
```bash
pip install -r https://raw.githubusercontent.com/phievo/phievo/master/extra.txt
jupyter nbextension enable --py --sys-prefix widgetsnbextension
```

### Test your installation

TO test that everything works properly, we will run an simulation example.

Copy the project directory `Examples/Somites` and `run_evolution.py` fom  [github](https://github.com/phievo/phievo) on your computer. Then copy `run_evolution.py` at the same place as the  `Somites/` directory.

To launch the evolution, simply run

```bash
	./run_evolution.py -m Somites
```

On windows machine  we recommand that you explicitly tell the system that you are running python (make sure you use the good version).

```bash
	python run_evolution.py -m Somites
```

If everything works correctly you should see the evolution starting and regular terminal print of the population best fitness.

You can also choose to stop the simulation by deleting the `Somites/STOP.txt` file after a few generations. The [jupyter notebook](https://github.com/phievo/phievo/blob/master/Analyse%20Run.ipynb) can then be use to visualize the results.
