# φ-evo

## Install phievo


Since _phievo_ is not pusblished on any offical repository yet, you need to install it manually.

### install Anaconda
I you do not have python already installed on you computer, we recommand to install it via the [anaconda distribution](https://www.continuum.io/downloads).

Among other things, anaconda provides the standard package manager of python _pip_. _pip_ can upgrade itself, it is generally recommanded to do so before starting the installation process.

```bash 
	pip install --upgrade pip
```

When python2 and python3 are installed on the same computer, it is common to need to specify the version of python or pip you are using: `python` (`pip`) for python2 and `python3` (`pip3`) for python3. Make sure you verify which command correspond to python3 on you computer.

### install the package


Simply run the following command:
```bash 
	sudo pip install dist/phievo-1.0.tar.gz
```

## Example: Static Hox

Copy the project directory `StaticHox` from `Examples` where you want to run it. Then copy `run_evolution` at the same place as `StaticHox`.

To launch the evolution, simply run

```bash
	./run_evolution -m StaticHox
```
