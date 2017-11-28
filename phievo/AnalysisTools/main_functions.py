import numpy as np
import shelve
import sys,os,glob,pickle,zipfile,re
from urllib.request import urlretrieve
from phievo.AnalysisTools import palette
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from  matplotlib.lines import Line2D
def read_network(filename,verbose=False):
    """Retrieve a whole network from a pickle object named filename

    Args:
        filename (str): the directory where the object is saved

    Returns:
        The stored network
    """
    with open(filename,'rb') as my_file:
        net = pickle.load(my_file)
    if verbose:
        print("Network retrieve from: {}".format(filename))
    return net

def smoothing(array,param):
    """Smoothen an array by averaging over the neighbourhood

    Args:
        array (list): the to be smoothed array
        param (int): the distance of the neighbourhood

    Returns:
        list of same size as array
    """
    length = len(array)
    return [np.mean(array[max(0,i-param):i+param+1]) for i in range(length)]

def load_generation_data(generations,restart_file):
    """
        Searches in the restart file the the informations that has been backed up
        up about the individuals at  a given generations.

        Args:
            generations (list): index of the generations to load_generation_data
            restart_file: path of the restart_file
        Returns:
            dictionary where each key contains the informations about one generation.
    """
    gen_data = {}

    with shelve.open(restart_file) as data:
        restart_generations = sorted([int(xx) for xx in data.dict.keys()])
        for gen in generations:
            if gen not in restart_generations:
                limit_print = 20
                err_str = ""
                err_str += "Generation {0} is not saved in the  restart file.\n".format(gen)
                err_str += "Please choose among the following generations:\n"
                if len(restart_generations)<limit_print:
                    err_str+=", ".join([str(x) for x in restart_generations[:limit_print]])+"\n"
                else:
                    err_str+=", ".join([str(x) for x in restart_generations[:limit_print]])+", etc.\n"
                raise AssertionError(err_str)
            dummy,gen_data[gen] = data[str(gen)]
    return gen_data




            
def download_zip(dir_name,url):
    """
    Download and extract zip file to dir_name.
    """
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    else:
        print("The directory {} already exists, download_example_seed cannot overwrite it.".format(dir_name))
        return 0
    ## Downloading zipfile        
    zip_path = os.path.join(dir_name,dir_name+".zip")
    def dlProgress(count, blockSize, totalSize):
        state = int(count * blockSize * 100 / totalSize)
        if state%2==0:
            print("{}: [".format(dir_name+".zip")+("#"*int(state/2))+(" "*(50-int(state/2)))+"] {}%".format(state),end="\r")
    urlretrieve(url,zip_path,reporthook=dlProgress)
    print("{}: [".format(dir_name+".zip")+("#"*50)+"] 100%",end="\n")
    ## unziping file
    print("Extracting zip file...",end="\r")
    
    zip_ref = zipfile.ZipFile(zip_path, 'r')
    zip_ref.extractall(dir_name)
    zip_ref.close()
    print("Extracting zip file...   done.",end="\n")
    print("Deleting zip file...",end="\r")
    os.remove(zip_path)
    print("Deleting zip file...   done.",end="\n")
    return 1
    
    
def download_tools(run_evolution="run_evolution.py",AnalyseRun="AnalyseRun.ipynb"):
    url_runevo = "https://raw.githubusercontent.com/phievo/phievo/master/run_evolution.py"
    url_jpnb = "https://github.com/phievo/phievo/raw/master/Analyse%20Run.ipynb"
    urlretrieve(url_runevo,run_evolution)
    print("run_evolution.py ... downloaded.")
    urlretrieve(url_jpnb,AnalyseRun)
    print("AnalyseRun.ipynb ... downloaded.")
    
def download_example(example_name,directory=None):
    """
    Download an example seed or project.
    """

    #server_address = "http://www.physics.mcgill.ca/~henrya/seeds_phievo/{}"
    server_examples = "https://github.com/phievo/phievo/blob/master/Examples/{}?raw=true"
    
    existing_examples = {
        "adaptation":"adaptation.zip",
        "somite":"Somites.zip",
        "hox":"StaticHox.zip",
        "hox_pareto":"StaticHox_pareto.zip",
        "lac_operon":"lac_operon.zip",
        "immune":"immune.zip",
    }
    server_seed = "https://github.com/phievo/simulation_examples/blob/master/{}?raw=true"
    existing_seeds = {
        "seed_adaptation":"adaptation.zip",
        "seed_adaptation_pruning":"adaptation_pruning.zip",
        "seed_lacOperon":"lacOperon.zip",
        "seed_lacOperon_pruning":"lacOperon_pruning.zip",
        "seed_somite":"somite.zip",
        "seed_somite_pruning":"somite_pruning.zip",
        "seed_hox_pareto_light":"hox_pareto_light.zip",
    }
    
    with_seed = False
    if "seed" in example_name:
        with_seed = True
        try:            
            zip_name = existing_seeds[example_name]
            example_name = example_name[5:]
            url = server_seed.format(zip_name)
        except KeyError:
            print("Example {} is not available.".format(example_name))
            print("Only the following examples are available:\n\t- "+"\n\t- ".join(list(existing_examples.keys())+list(existing_seeds.keys())))
            return None
    else:
        try:            
            zip_name = existing_examples[example_name]
            url = server_examples.format(zip_name)
        except KeyError:
            print("Example {} is not available.".format(example_name))
            print("Only the following examples are available:\n\t- "+"\n\t- ".join(list(existing_examples.keys())+list(existing_seeds.keys())))
            return None
    if not directory:    
        directory = "example_{}".format(example_name)
    res = download_zip(directory,url)
    if not res:
        return None
    
    if with_seed:
        seed_name = os.path.join(directory,"Seed{}".format(example_name))
        
        os.makedirs(seed_name)
        files = glob.glob(os.path.join(directory,"*"))
        files.remove(seed_name)
        for filename in files :
            try:
                os.rename(filename, filename.replace(directory,seed_name))
            except OSError:
                import pdb;pdb.set_trace()
        print("recovering log files...",end="\r")
        for log_f in glob.glob(os.path.join(seed_name,"log_*")):
            f_name = log_f.split(os.sep)[-1]
            f_name = f_name.replace("log_","")
            os.rename(log_f, os.path.join(directory,f_name))
        with open(os.path.join(directory,"init_file.py"),"r") as init_file:
            init_text = init_file.read()
            init_text = re.sub("(cfile\[[\'\"](\w+)[\'\"]]\s*=\s*).+",r"\1'\2.c'",init_text)
            init_text = re.sub("(pfile\[[\'\"](\w+)[\'\"]]\s*=\s*).+",r"\1'\2.py'",init_text)
        with open(os.path.join(directory,"init_file.py"),"w") as init_file:
            init_file.write(init_text)
        print("recovering log files...   done.",end="\n")
    print("Project saved in {}.".format(directory))
    
 
