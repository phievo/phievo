"""Module in charge of creating various utilities to initialize evolutionary codes: loading specific modules, test files, creating directories,etc. """

from phievo import __silent__,__verbose__
if __verbose__:
    print("Execute initialization_code.py")

import os
import optparse
import random
import shutil
import sys
from importlib import import_module

### working path handling ###
# This part should be slowly erased as we migrate to python3
python_path = os.getcwd()
sys.path.insert(1, python_path)  # make evol modules second item in path after cwd.
ccode_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'CCodes')

multiple_phospho = 0

# ################################################################################
# routines to use imported initialization.py module (inits), extract data and add
# to other imported modules, and return modified modules

def init_evolution(inits, deriv2):
    """import modules mutation and evolution_gillespie and initialize
    """
    # input mutation rates into mutation.py
    import phievo.Networks.mutation as mutation

    mutation.dictionary_ranges.update(inits.dictionary_ranges) #Merge instead of overwriting
    mutation.list_types_output = inits.list_types_output

    # reset delay to units of number of steps
    mutation.dictionary_ranges['CorePromoter.delay'] = int(
        0.5 + mutation.dictionary_ranges['CorePromoter.delay'] / inits.prmt['dt'])

    mutation.dictionary_mutation.update(inits.dictionary_mutation)
    mutation.build_lists(mutation.dictionary_mutation)

    # import parameters/fns into gillespie evolution module

    import phievo.Populations_Types.evolution_gillespie as evo_gis
    evo_gis.prmt.update(inits.prmt)
    evo_gis.compile_and_integrate = deriv2.compile_and_integrate
    mutation.compile_and_integrate = deriv2.compile_and_integrate

    try:
        evo_gis.init_network = inits.init_network
        if __verbose__:
            print('initializing with network from initization file')
    except Exception:
        msg = 'initial network from initialization file could not be found, bye'
        display_error(msg)

    try:
        evo_gis.fitness_treatment = inits.fitness_treatment
    except Exception:
        msg = 'Using default fitness_treatment() in module evolution_gillespie'
        display_error(msg)

    return [mutation, evo_gis]


def check_model_dir(model):
    """input string with name of model directory relative to CWD, check for init* file
    import initialization data and rename 'inits', and return with model-directory.

    Args:
        model (str): name of the model directory.
    Returns:
        list: [model_dir, init_module] with
            - model_dir(str): name of the model directory
            - init_module(dictionary): dictionary with the model's options

    """
    model_dir =  model # works as well as full path
    model_files = os.listdir(model_dir)
    init_candidate = [ff for ff in model_files if ff.startswith('init') & ff.endswith('.py')]

    ##adding to sys.modules does not work, syntax with <> ??  insert ahead of CWD
    sys.path.insert(0, model_dir)  # put model_dir ahead of '' incase duplicate names.
    init_name = init_candidate[0]
    init_name_nopy = init_name.replace('.py', '')
    init_module = import_module(init_name_nopy)
    inits = init_module
    try:
        if __verbose__:
            print('initializing with file=', init_name, 'in model dir=', os.path.abspath(model_dir))
    except IndexError:
        msg = "The program did not find any init file in the %s repository."%model_dir
        display_error(msg)
        raise SystemExit
    # add default configuration for pfile
    if hasattr(init_module,"pfile"):

        for key,ff in init_module.pfile.items():
            ff=ff.replace(".",os.sep)            
            if not os.path.isfile(ff+".py"):
                if os.path.isfile(os.path.join(model_dir,ff)+".py"):
                    inits.pfile[key] = os.path.join(model_dir,ff).replace(os.sep,".")
                else:
                    raise FileNotFoundError("ERROR: A python file cannot be found:\n{} doesn't match a file.".format(ff+".py"))
    if hasattr(init_module,"c_libraries"):
        for key,ff in init_module.c_libraries.items():
            init_module.c_libraries[key].setdefault("dir",model_dir)
            init_module.c_libraries[key].setdefault("o",key+".o")
    if hasattr(init_module,"c_dependencies"):
        for key,ff in init_module.c_dependencies.items():            
            init_module.c_dependencies[key] = os.path.join(model_dir,ff)
            
    for key,ff in init_module.cfile.items():
        if not os.path.isfile(ff):
            if os.path.isfile(os.path.join(model_dir,ff)):
                inits.cfile[key] = os.path.join(model_dir,ff)
            else:
                raise FileNotFoundError("ERROR: The c file cannot be found:\n{} doesn't match a file.".format(ff))

    pfile = {"deriv2" : "phievo.Networks.deriv2",
             "interaction" : "phievo.Networks.interaction",
             "pretty_graph": "phievo.Networks.lovelyGraph",
             "plotdata" : "phievo.Networks.plotdata"}
    cfile = {
        "header" : os.path.join(ccode_dir,'integrator_header.h'),
        "utilities" : os.path.join(ccode_dir,'utilities.c'),
        "geometry" : os.path.join(ccode_dir,'linear_geometry.c'),
        "integrator" : os.path.join(ccode_dir,'euler_integrator.c'),
        "main" : os.path.join(ccode_dir,'main_general.c')
        }

    try:
        pfile.update(init_module.pfile)
    except AttributeError:
        if __verbose__:
            print('Remark: No pfile object in init* file, use default one (see check_model_dir)!')
    try:
        cfile.update(init_module.cfile)
    except AttributeError:  
        raise AttributeError("The cfile is not defined in initialization.")

        
    init_module.pfile = pfile
    init_module.cfile = cfile

    setattr(init_module,"model_dir",model_dir)
    return [model_dir, init_module, model_dir+os.sep+init_name]


def init_networks(inits):
    """
    Uses the inits object(import from the initialization.py) to set the different
    components of the classes Networks. It tells the code where to find the different
    c files or python classes and overwrites the default values if necessaries.

    Args:
        inits: Import of the initialization file.
    Returns:
        deriv2 object
    """
    import phievo.Networks.classes_eds2 as net_class
    from phievo import Networks
    for k in list(inits.__dict__.keys()):
        #copies all attributes, but might restrict to some only in the future
        setattr(net_class, k, inits.__dict__[k])
    net_class.list_unremovable = inits.list_unremovable

    interaction = import_module(inits.pfile["interaction"])
    try:
        #pretty_graph = import_module(inits.pfile["pretty_graph"])
        setattr(Networks,"pretty_graph",inits.pfile["pretty_graph"])
    except (KeyError, AttributeError) as e:
        #pretty_graph = import_module('phievo.Networks.lovelyGraph')
        #setattr(Networks,"pretty_graph",pretty_graph)
        setattr(Networks,"pretty_graph",'phievo.Networks.lovelyGraph')

    deriv2 = import_module('phievo.Networks.deriv2')
    deriv2.cfile.update(inits.cfile)

    try:
        deriv2.c_libraries.update(inits.c_libraries)
    except AttributeError:
        pass

    try:
        deriv2.c_dependencies.update(inits.c_dependencies)
    except AttributeError:
        pass

    # Define default directory for cfile then overwrite with information from inits
    #deriv2.cfile['header'] = os.path.join(ccode_dir,'integrator_header.h')

    if ('langevin_noise' in inits.prmt):
        if (inits.prmt['langevin_noise'] > 0):
            deriv2.noise_flag = 1

    if inits.pfile["deriv2"] and inits.pfile["deriv2"] != "phievo.Networks.deriv2":
        mod_deriv2 = import_module(inits.pfile["deriv2"])
        deriv2 = mod_deriv2.modifier(deriv2)

    return deriv2

def get_network_from_test(test_py, classes_eds2):
    """Imports the file test_py file as a module and loops over its attirbute to find an object of type classes_eds2.Network.

    Args:
        test_py: Python file defining a classes_eds2 network.
        classes_eds2: An instance of the the class_eds2 used to interprete the network
    Returns:
        Network object
    """
    test_file = test_py.replace('.py', '')
    test_file = import_module(test_file)
    ## Test every modules imported from the source file to find one of the type classes_eds2.Network
    for kk in  dir(test_file):
        obj = getattr(test_file,kk)
        if isinstance(obj, classes_eds2.Network):
            return obj

    print('No Network object located in file=', test_py)
    return None

def initialize_test(init_py):
    """from options.init import the inits module, and if not supplied extract prmt dictionary from
    the test file.
    """
    if init_py:
        inits = import_module(init_py.replace('.py', ''))
    else:
        print('No init<*>.py file supplied, can not time step test network')
        return None

def parameters2file(inits, file_name):
    """ from the inits file save various things in named text file
    """
    fh = open(file_name, 'wt')
    fh.write('ccode_dir\n')
    fh.write(ccode_dir + '\n\n')
    fh.write('cfile dictionary\n')
    fh.write(str(inits.cfile) + '\n\n')
    fh.write('prmt dictionary\n')
    fh.write(str(inits.prmt) + '\n\n')
    fh.write('ranges dictionary\n')
    fh.write(str(inits.dictionary_ranges) + '\n\n')
    fh.write('mutation dictionary\n')
    fh.write(str(inits.dictionary_mutation) + '\n\n')
    fh.close()
    print('wrote initialization parameters to file=', file_name)

def make_workplace_dir(parent_dir):
    """define name of dir that will contain C code and exececutables, and creates it
    if not already present
    """

    workplace_dir = os.path.join(parent_dir,'Workplace')
    workplace_dir = os.path.normpath(workplace_dir)
    if ( not os.access(workplace_dir, os.F_OK) ):
        os.mkdir(workplace_dir)
    return workplace_dir # in other modules use dir + filename, inconsistently

def display_error(msg = 'ERROR',filename = 'error.txt'):
    """A pretty function to display the catched errors and log them in a file

    Args:
        msg (str): A commentary on the error you wan't to add
        file (str): the file where errors are loged

    Returns:
        None: print and write in file only
    """
    import sys,time,traceback
    what,args,track = sys.exc_info()[:3]
    what = str(what).split("'")[1]
    file,line = track.tb_frame.f_locals.get('__file__','FNF'),track.tb_lineno
    out = str(msg)+'\n'
    out += 'Type: {what}, says {args}\n'.format(**locals())
    out += str(traceback.extract_tb(track)[0])
    print(out)
    with open(filename,'a') as myfile:
        myfile.write(time.ctime(time.time())+'\n'+out+'-'*50+'\n')
