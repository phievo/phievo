from phievo.initialization_code import *
from importlib import import_module
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
    [model_dir, inits, init_dir] = check_model_dir(options["model"])
    if options["clear"]:
        for directory in glob.glob(os.path.join(model_dir,"Seed*"))+glob.glob(os.path.join(model_dir,"Workplace*")):
            shutil.rmtree(directory)
        
    [classes_eds2, pretty_graph] = init_classes_eds2(inits)
    workplace_dir = make_workplace_dir(model_dir)
    deriv2 = init_deriv2(inits, workplace_dir, inits.prmt)
    [mutation, evolution_gillespie] = init_evolution(inits, deriv2)

    # to distinguish master and slave nodes when running on cluster with pypar
    main_loop = False
    if (inits.prmt['multipro_level'] == 2):
        import pypar #parallel programming module
        if pypar.rank() == 0: main_loop = True
    else:
        main_loop = True
    if main_loop: # this part describe only the master processor 0
        # Recovery from restart file
        if (inits.prmt['restart']['activated']
            and inits.prmt['restart']['same_seed']
            and inits.prmt['nseed'] > 1):
            print('WARNING initializing from dir=', inits.prmt['restart']['dir'], 'and exactly continuing prior data')
            print('Therefore no need to for nseed=', inits.prmt['nseed'], 'to be >1, resetting to 1')
            inits.prmt['nseed'] = 1
            # 11/2010 EDS noticed bug here, identical restart not working?? problem with rand seed??

        if 'firstseed' in inits.prmt:
            firstseed = inits.prmt['firstseed']
        else:
            firstseed = 0

        for seed in range(firstseed, firstseed + inits.prmt['nseed']):
            print('initializing random() with seed=', seed, 'prior to beginning the evolution')
            random.seed(seed)
            namefolder = model_dir + "Seed%i" % seed

            # Create a directory if needed and check if data already present
            if os.access(namefolder, os.F_OK):
                if (len(os.listdir(namefolder)) > 2):  #ok to overwrite paramter file, and Bests but not simulation data
                    message = 'dir= {0} has data in it, exiting'
                    sys.exit(message.format(namefolder))
            else:
                os.mkdir(namefolder)
                # Copy some inits file in the Seed directory
                shutil.copyfile(inits.cfile['fitness'],namefolder+os.sep+'log_fitness.c')
                shutil.copyfile(inits.cfile['input'],namefolder+os.sep+'log_input.c')
                shutil.copyfile(inits.cfile['init_history'],namefolder+os.sep+'log_init_histo.c')
                shutil.copyfile(init_dir,namefolder+os.sep+'log_init_file.py')

            parameters2file(inits, namefolder + '/parameters')

            # Population construction for run on several machine with pypar
            if (inits.prmt['multipro_level'] == 2):
                if (inits.prmt['pareto'] == 1):
                    from phievo.Populations_Types.pareto_population import pareto_parallel_Population
                    population = pareto_parallel_Population(namefolder, inits.prmt['npareto_functions'],
                                                            inits.prmt['rshare'])
                else:
                    from phievo.Populations_Types.parallel_population import parallel_Population
                    population = parallel_Population(namefolder)

            # Population construction for multiprocessor run on one machine
            elif (inits.prmt['multipro_level'] == 1):
                if (inits.prmt['pareto'] == 1):
                    from phievo.Populations_Types.pareto_population import pareto_thread_Population
                    population = pareto_thread_Population(namefolder, inits.prmt['npareto_functions'],
                                                          inits.prmt['rshare'])
                else:
                    from phievo.Populations_Types.thread_population import thread_Population
                    population = thread_Population(namefolder)

            # Population construction for single processor run
            else:
                if (inits.prmt['pareto'] == 1):
                    from Populations_Types.pareto_population import pareto_Population
                    population = pareto_Population(namefolder, inits.prmt['npareto_functions'], inits.prmt['rshare'])
                else:
                    population = evolution_gillespie.Population(namefolder)

            # Finaly launch the genetic algorithm
            population.evolution()

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

def test_network(options):
    """ Test the behavior of a particular network (indicated by the -t
    option) with respect to a given model (the -m or -i option)

    Args:
        options (optparse.Values): a dictionnary like class containing
        the network (options.test) and model (options.model) or init*.py
        (options.inits) file path.

    Returns:
        None
    """
    plotdata = import_module(name_plotdata)
    
    print(plotdata)
    # Define the model to be used

    if (options["init"]):        
        inits = initialize_test(options["init"])
        test_output_dir = './'
        [classes_eds2, pretty_graph] = init_classes_eds2(inits)
    elif (options["model"]):
        [model_dir, inits, init_dir] = check_model_dir(options["model"])
        test_output_dir = model_dir
        [classes_eds2, pretty_graph] = init_classes_eds2(inits)
    else:
        inits = False

    # Initialize the network and the working directory
    net = get_network_from_test(options["test"], classes_eds2)
    net.write_id()
    workplace_dir = make_workplace_dir(test_output_dir)
    gr = net.draw(workplace_dir+"test_net.pdf")
    print('network diagram= test_net.pdf (all nodes), *dot (species only) drawn for test file=', options["test"])
    print('output sent to dir= ', test_output_dir)

    # Launch the integraion with relevant parameters
    if inits:
        deriv2 = init_deriv2(inits, workplace_dir, inits.prmt)
        plotdata.compile_and_integrate = deriv2.compile_and_integrate
        print('running test file with initialization,')
        my_ntries = inits.prmt['ntries']
        ncell = int(options["ncell"]) if options["ncell"] else -1 #default value
        liszt = [int(n) for n in options["list"].split(',')] if options["list"] else [] #default value
        plotdata.net_test_plot(net, inits.prmt, test_output_dir, 0, my_ntries,ncell,liszt)
    else:
        print("No init* file found, use -m or -i option to specify one")
