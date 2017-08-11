from phievo.initialization_code import *
import phievo
from importlib import import_module
import time,random
import os,shutil,glob

### Functions ###
def launch_evolution(options):
    """ The actual top routine, run evolution according to the models directory
    (indicated through -m), and places all output in /model directory
    It is at this level that the programm is spread across the different
    processors if needed and final operation of initialization are done.
    The very genetic algorithm is the evolution method from population class.

    Args:
        options (optparse.Values): a dictionnary like class containing
        the model directory (options.model).

    Returns:
        None
    """

    #initializes deriv2 and mutation for all processors

    [model_dir, inits, init_file] = check_model_dir(options["model"])
    ### Write STOP file
    if options["clear"]:
        clear_project(model_dir,inits)
    if options["network"]:
        def net_func():
            return phievo.read_network(options["network"])
        setattr(inits,"init_network",net_func)
    inits.prmt["stop_file"] = phievo.create_STOP_file(model_dir)
    deriv2 = init_networks(inits)


    [mutation, evolution_gillespie] = init_evolution(inits, deriv2)


    # to distinguish master and slave nodes when running on cluster with pypar
    main_loop = False

    if (inits.prmt['multipro_level'] == 2):
        raise NotImplementedError("pypar belongs to the past. We will release a version working with Multiprocess soon.")
        import pypar #parallel programming module
        if pypar.rank() == 0: main_loop = True
    else:
        main_loop = True
    if main_loop: # this part describe only the master processor 0
        # Recovery from restart file
        if (inits.prmt['restart']['activated'] and inits.prmt['nseed'] > 1):
            if inits.prmt['restart'].get('seed',None) is None:
                ## If no seed is provided, searches the one with the largest index
                seeds = glob.glob(os.path.join(model_dir,"Seed*"))
                if len(seeds) == 0:
                    raise FileExistsError("No seed to start from in {0}.".format(model_dir))

                inits.prmt['restart']['seed'] = max([int(seed.replace(os.path.join(model_dir,"Seed"),"")) for seed in seeds])
            print('WARNING initializing from Seed{0}'.format(inits.prmt['restart']['seed']), 'and exactly continuing prior data')
            #print('Therefore no need to for nseed=', inits.prmt['nseed'], 'to be >1, resetting to 1')
            #inits.prmt['nseed'] = inits.prmt['restart']['seed']
            # 11/2010 EDS noticed bug here, identical restart not working?? problem with rand seed??
            inits.prmt['firstseed'] = inits.prmt['restart']['seed']
        if 'firstseed' in inits.prmt:
            firstseed = inits.prmt['firstseed']
        else:
            firstseed = 0

        ## The following line allows running multiple runs in parallel on the same project
        ## without interfering.
        seeds = list(range(firstseed, firstseed + inits.prmt['nseed']))
        time.sleep(random.random()*10)
        while seeds:

            done_seeds = map(lambda path:os.path.split(path)[-1],glob.glob(os.path.join(model_dir,"Seed*")))
            seed = seeds[0]
            try:
                seeds.remove(seed)
            except ValueError:
                pass
            if seed in done_seeds:
                continue


        #for seed in range(firstseed, firstseed + inits.prmt['nseed']):
            try:
                launch_seed(seed,inits,init_file)
            except KeyboardInterrupt:
                print("\n\tThe run was interrupted by the user.")
                os._exit(0)


    else: # this part describes the slave nodes under pypar
        while True:
            #Wait for a message from the master
            message = pypar.receive(source=0)
            #Evaluate the message and send the result back to the master
            command, msg_locals = message
            locals().update(msg_locals)
            result = eval(command)
            pypar.barrier()  #wait here until job completion of all other processors, do not know why, should try to remove it at some point for performance
            pypar.send(result, 0)

def clear_project(model_dir=None,inits=None,options = None):
    """
    Clears a project's old run files.

    Args:
        model_dir: project directory
        inits: initialization parameters
    """
    if options:
        [model_dir, inits, init_file] = check_model_dir(options["model"])
    if inits.prmt["restart"]["activated"]:
        print("The clear(-c) option can not be activated when prmt[\"restart\"][\"activated\"] is set to True.")
        os._exit(0)
    toClear = glob.glob(os.path.join(model_dir,"__pycache__")) + glob.glob(os.path.join(model_dir,"Seed*")) + glob.glob(os.path.join(model_dir,"Workplace"))
    ## Remove dictionnaries
    for directory in toClear:
        shutil.rmtree(directory, ignore_errors=True)
    ## Remove files
    toClear = glob.glob(os.path.join(model_dir,"STOP.txt")) + glob.glob(os.path.join(model_dir,"Buffer*"))
    for ff in toClear:
        os.remove(ff)

def launch_seed(seed,inits,init_file):
    """
    Launch the evolution for a new seed.

    Args:
        seed: index of the seed to run
        inits: initialization parameter dictionnary
                (obtained from initialization file)
    """
    print('initializing random() with seed=', seed, 'prior to beginning the evolution')
    random.seed(seed)
    namefolder = os.path.join(inits.model_dir,"Seed%i" % seed)

    # Create a directory if needed and check if data already present
    if os.access(namefolder, os.F_OK):
        if (len(os.listdir(namefolder)) > 2) and not inits.prmt["restart"]["activated"]:  #ok to overwrite paramter file, and Bests but not simulation data
            message = 'dir= {0} has data in it, exiting'
            sys.exit(message.format(namefolder))
    else:
        os.mkdir(namefolder)
        # Copy some inits file in the Seed directory
        shutil.copyfile(inits.cfile['fitness'],namefolder+os.sep+'log_fitness.c')
        shutil.copyfile(inits.cfile['input'],namefolder+os.sep+'log_input.c')
        shutil.copyfile(inits.cfile['init_history'],namefolder+os.sep+'log_init_histo.c')
        shutil.copyfile(init_file,namefolder+os.sep+'log_init_file.py')

    parameters2file(inits, os.path.join(namefolder,'parameters'))

    # Population construction for run on several machine with pypar
    if (inits.prmt['multipro_level'] == 2):
        if (inits.prmt['pareto']):
            from phievo.Populations_Types.pareto_population import pareto_parallel_Population
            population = pareto_parallel_Population(namefolder, inits.prmt['npareto_functions'],
                                                    inits.prmt['rshare'])
        else:
            from phievo.Populations_Types.parallel_population import parallel_Population
            population = parallel_Population(namefolder)

    # Population construction for multiprocessor run on one machine
    elif (inits.prmt['multipro_level'] == 1):
        if (inits.prmt['pareto']):
            from phievo.Populations_Types.pareto_population import pareto_thread_Population
            population = pareto_thread_Population(namefolder, inits.prmt['npareto_functions'],
                                                  inits.prmt['rshare'])
        else:
            from phievo.Populations_Types.thread_population import thread_Population
            population = thread_Population(namefolder)

    # Population construction for single processor run
    else:
        if (inits.prmt['pareto']):
            from phievo.Populations_Types.pareto_population import pareto_Population
            population = pareto_Population(namefolder, inits.prmt['npareto_functions'], inits.prmt['rshare'])
        else:
            from phievo.Populations_Types.evolution_gillespie import Population
            population = Population(namefolder)

    # Finaly launch the genetic algorithm
    inits.prmt["workplace_dir"] = make_workplace_dir(os.path.join(inits.model_dir,"Seed{0}".format(seed)))
    population.evolution(inits.prmt)

def test_project(project_path,network=None,return_sim= False):
    """
    Test the project on the initial file.
     - Load the initial network from the initialization file.
     - Genertate the C file containing the ODEs
     - Compile and integrade
     - Print the the network's fitness
     - plot the fitness


    Args:
        project_path: Project to test
        network: path to a test .net file
        return_sim: return The simulation object after running the trials.
    Returns:
        None
    """
    from phievo.AnalysisTools import Simulation
    sim = Simulation(project_path,mode="test")
    if network:
        net = phievo.read_network(network)
    else:
        net = sim.inits.init_network()
    if return_sim:
        sim.run_dynamics(net=net,trial=sim.inits.prmt['ntries'],erase_buffer=True,return_treatment_fitness=False)
        return sim
    data = sim.run_dynamics(net=net,trial=sim.inits.prmt['ntries'],erase_buffer=True,return_treatment_fitness=True)
    cfile = glob.glob(os.path.join(project_path,"Workplace","*.c"))[0]
    print("C file created.")
    print("You can recompile by running:")
    print(str.encode("gcc {cfile} -lm -o executable".format(cfile=cfile)))
    print(data)
    #print("The network's C file and executable are stored in {}".format())
    net.draw()
