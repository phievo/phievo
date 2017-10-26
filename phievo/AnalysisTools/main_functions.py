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


def download_example_seed(seed_name):
    """
    Downloads a seed from the example seed repository.
    """
    #server_address = "http://www.physics.mcgill.ca/~henrya/seeds_phievo/{}"
    server_address = "https://github.com/phievo/simulation_examples/blob/master/{}?raw=true"
    existing_seeds = {
        "adaptation":"adaptation.zip",
        "adaptation_pruning":"adaptation_pruning.zip",
        "lacOperon":"lacOperon.zip",
        "lacOperon_pruning":"lacOperon_pruning.zip",
        "somite":"somite.zip",
        "somite_pruning":"somite_pruning.zip",
        "hox_pareto_light":"hox_pareto_light.zip",
    }
    existing_seeds = {kk:server_address.format(val) for kk,val in existing_seeds.items()}
    try:
        url = existing_seeds[seed_name]
    except KeyError:
        print("Only the following examples are available:\n\t- "+"\n\t- ".join(list(existing_seeds.keys())))
        return None
    directory = "example_{}".format(seed_name)
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        print("The directory {} already exists, download_example_seed cannot overwrite it.".format(directory))
        return None
    ## Downloading zipfile
    def dlProgress(count, blockSize, totalSize):
        state = int(count * blockSize * 100 / totalSize)
        if state%2==0:
            print("{}: [".format(seed_name+".zip")+("#"*int(state/2))+(" "*(50-int(state/2)))+"] {}%".format(state),end="\r")
        
    zip_path = os.path.join(directory,seed_name+".zip")
    urlretrieve(url,zip_path,reporthook=dlProgress)
    print("{}: [".format(seed_name+".zip")+("#"*50)+"] 100%",end="\n")
    ## unziping file
    print("Extracting zip file...",end="\r")
    seed_path = os.path.join(directory,"Seed{}".format(seed_name))
    zip_ref = zipfile.ZipFile(zip_path, 'r')
    zip_ref.extractall(seed_path)
    zip_ref.close()
    print("Extracting zip file...   done.",end="\n")
    print("Deleting zip file...",end="\r")
    os.remove(zip_path)
    print("Deleting zip file...   done.",end="\n")
    print("recovering log files...",end="\r")
    for log_f in glob.glob(os.path.join(seed_path,"log_*")):
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
    
 
