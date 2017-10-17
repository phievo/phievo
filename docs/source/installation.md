# Install φ-evo


Since _φ-evo_ is not pusblished on any official repository yet, you need to install it manually.

### install Anaconda
In case python is not already installed on your computer, we recommand to install it via the [anaconda distribution](https://www.continuum.io/downloads).

Among other things, anaconda provides the standard package manager of python _pip_. _pip_ can upgrade itself, it is generally recommanded to do so before starting the installation process.

```bash
	pip3 install --upgrade pip
```

When python2 and python3 are installed on the same computer, it is common to need to specify the version of python or pip you are using: `python` (`pip`) for python2 and `python3` (`pip3`) for python3. Make sure you verify which command correspond to python3 on you computer. In the following instruction you may need replace "pip3" by "pip".

### install the package

From the root of the project, run the following command:
```bash
	sudo pip3 install https://github.com/phievo/phievo/blob/master/dist/phievo-1.0.tar.gz?raw=true
```

## Instructions specific to windows

### Install gcc
Windows does not come with the `gcc` compiler installed but the free software foundation provides a minimal distribution of the gnu softwares for windows, it is called [MinGW](http://mingw.org/).

Once you have downloaded `mingw-get-setup.exe`, run it. A selection panel will open. We recommend you to install at least the two following packages(the others are not relevant for phievo):
- mingw-developper-toolkit
- mingw32-base

Choose the default directory.

After the installation is finished, update windows `PATH` so that it knows where to look for _gcc_ command. Open a the command prompt and run:

```bash
setx PATH "%path%;C:\MinGW\bin"
```

**Note:** If you may use other coding distribution such as code blocks or visual basics that already  contain the _gcc_ compiler. In such case, you do not need to install MinGw. Just upload you `PATH` so that windows knows where is the gcc compiler.

## Instructions specific to OSX

### Install gcc
OSX does not have the gcc compiler installed by default. There are different way to install it. The fastest is probably via [homebrew](https://brew.sh/):

```bash
brew install gcc
```

If _gcc_ is not already installed on you system (via macports or Xcode), _homebrew_'s _gcc_ should be automatically in the system's `PATH`.


## Install pygraphviz

_pygraphviz_ is not istalled by default with _phievo_ because it does not exist natively on windows and we wanted to publish a version that that run on al the systems. _pygraphviz_ is used only to display network layouts. If it is not install _phievo_ will print a warning and use _networkx_ spring layout instead.

On max OSX, you have to use homebrew to install graphvix first :

```bash
brew install graphviz pkg-config
pip3 install pygraphviz
```

On GNU/linux, installing the dependencies varies depanding on the distribution. We tested the following on debian and ubuntu

```bash
sudo apt-get install graphviz graphviz-dev pkg-config
sudo pip3 install pygraphviz
```
On other distribution, you may want to find the equivalent of _graphviz_, _graphviz-dev_, and _pkg-config_.

We found that sometimes on ubuntu the C linking to the graphviz library does not work properly, to fix this, be more explicit on the linking:

```bash
sudo pip3 install pygraphviz --install-option="--include-path=/usr/include/graphviz" --install-option="--library-path=/usr/lib/graphviz/"
```


## Build the documentation

To build the documentation, go to the `docs` directory and run

```bash
make html
```
the documentation is built in `docs/build/html` by default.

Note that you may need to run make html twice because the build process uses the package _numfig_ that first needs to list all the figures before numbering them.

## Analyse notebook

We provide a jupyter notebook to help with the analysis of the runs. If you wand to run it, you will need to install several extra python libraries, to help with this, they are writen in [extra.txt](extra.txt).
```bash
pip3 install -r extra.txt
jupyter nbextension enable --py --sys-prefix widgetsnbextension
```

## Example: Static Hox

Copy the project directory `StaticHox` from `Examples` where you want to run it. Then copy `run_evolution` at the same place as `StaticHox`.

To launch the evolution, simply run

```bash
	./run_evolution.py -m StaticHox
```

On windows machine  we recommand that you explicitly tell the system that you are running python (make sure you use the good version).

```bash
	python run_evolution.py -m StaticHox
```
