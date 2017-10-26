#!/usr/bin/env python3
"""Top level routine to launch the evolutionary algorithm

   the model (the directory) may be specified through the -m option (relative to CWD)
   it have to contain a file called init*.py and will collect all the output
   example: python run_evolution.py -m StaticHox

   a particular file describing a network object may also be tested with -t option
   in this case, an init* file should be provided, through a model with -m options
   or directly with -i option.
   a particular cell may be displayed with -n
   a detailed list of species to display may be ask with -l
   example: python run_evolution.py -m StaticHox -t cobaye.py

   Various subdirectories for data files created by this script if needed
   Type -h or --help for options.

   Need explicitly input the directory with python modules and the Ccode directory, below
   All inputs and their description in init*.py file.
"""
import optparse
import phievo
# Definition of the parser options
pp = optparse.OptionParser()
pp.add_option('--model', '-m', action='store',
              help='name of directory with file called init*.py and which receives all output of evolution simulation')
pp.add_option('--test', '-t', action='store',
              help='name of file in CWD with network to run as test (need -i or -m options to supply initialization file to time .) Output goes to dir defined by -m if used, or CWD. Test.py file created by stat_best_net.py')
pp.add_option('--network', '-n', action='store',
              help='Curtom initial network')
pp.add_option("--clear","-c", action="store_true", dest="clear",default=False)

(options, arg) = pp.parse_args()  # NB arg=[], but required output
options = options.__dict__

############
### MAIN ###
############
if __name__ == "__main__":
    if options["model"]:
        phievo.launch_evolution(options)
    elif options["test"]:
        phievo.test_project(options["test"],network=options["network"])
    elif options["clear"]:
        options["model"] = arg[0]
        phievo.clear_project(options=options)
    else:
        print('must specify either --test file OR --model directory.  Use -h for help')
    
