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
        Network: the object having been stored
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
        list: of same size as array
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

def plot_multiGen_front2D(generation_fitness,generation_indexes=None):
    """
        Uses the fitness data for multiple generations to represent the pareto fronts
        of those multiple generations.

        Args:
            generation_fitness: nested dictionnaries:
                                level0 keys: generation
                                level1 keys: rank of the fitness (1,2,etc.)
                                index : index of the fitness doublet (they might be
                                        multiple fitnesses with identical rank).
            generation_indexes: Same dictionnary structure as generation_fitness.
                                Contains the index of each network in its population

    """
    NUM_COLORS = len(generation_fitness)
    shapes = ["o","s","^"]
    color_l = palette.color_generate(NUM_COLORS)
    legend_patches = []
    #plt.legend(handles=[red_patch])
    i = 0
    fig = plt.figure()
    ax = fig.gca()
    for gen in sorted(generation_fitness.keys()):
        gen_dico = generation_fitness[gen]
        legend_patches.append(mpatches.Patch(color=color_l[i], label='Generation {0}'.format(gen)))
        color = color_l[i]
        i +=1
        for rank,points in gen_dico.items():
            F1,F2 = list(zip(*points))
            shape = shapes[rank-1] if rank<3 else shapes[-1]
            ax.scatter(F1,F2,c=color,edgecolor=color,s=50,marker=shape)
            if generation_indexes:
                ind_list = generation_indexes[gen][rank]
                for l in range(len(ind_list)):
                    ax.text(F1[l],F2[l],'%d' % ind_list[l],ha='center', va='bottom')
    ax.set_xlabel('Fitness 1')
    ax.set_ylabel('Fitness 2')
    ax.legend(handles=legend_patches)
    legend_patches.append(Line2D([0], [0], linestyle="none", marker=shapes[0], markersize=10,markerfacecolor="black",label="Rank 1"))
    legend_patches.append(Line2D([0], [0], linestyle="none", marker=shapes[1], markersize=10,markerfacecolor="black",label="Rank 2"))
    legend_patches.append(Line2D([0], [0], linestyle="none", marker=shapes[2], markersize=10,markerfacecolor="black",label="Rank≥3"))
    ax.legend(handles=legend_patches)
    plt.show()
    return fig

def download_example_seed(seed_name):
    """
    Downloads a seed from the seed repository.
    """
    existing_seeds = {
        "adaptation":"https://github.com/phievo/simulation_examples/blob/master/adaptation.zip?raw=true",
        "hox_pareto_light":"https://github.com/phievo/simulation_examples/blob/master/hox_pareto_light.zip?raw=true",
        "lacOperon":"https://github.com/phievo/simulation_examples/blob/master/lacOperon.zip?raw=true",
        "somite":"https://github.com/phievo/simulation_examples/blob/master/somitogenesis.zip?raw=true"
    }
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
    
        
    zip_path = os.path.join(directory,seed_name+".zip")
    urlretrieve(url,zip_path)
    ## unziping file
    seed_path = os.path.join(directory,"Seed{}".format(seed_name))
    zip_ref = zipfile.ZipFile(zip_path, 'r')
    zip_ref.extractall(seed_path)
    zip_ref.close()
    os.remove(zip_path)
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
