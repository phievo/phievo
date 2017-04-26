import numpy as np
import shelve

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
    cm = pylab.get_cmap('gist_rainbow')
    color_l= [colors.rgb2hex(cm(1.*i/NUM_COLORS)) for i in range(NUM_COLORS)]
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
