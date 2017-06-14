"""
This module provide a pareto_Population class to perform a Pareto evolution,
that is, a general frame to evolve Networks according to more than one fitness
function.

See: Warmflash, A., Francois, P., & Siggia, E. D. (2012).
Pareto Evolution of Gene Networks: An Algorithm to Optimize
Multiple Fitness Objectives. Physical Biology, 9(5), 56001.

Coder: A. Warmflash, P. François
"""
print('Execute pareto_population.py')

from phievo.Populations_Types.evolution_gillespie import Population
from phievo.Populations_Types.thread_population import thread_Population
import random
from phievo.Networks import classes_eds2
from math import log,sqrt

#######################
### Local Functions ###
#######################

def single_comparison(x,y):
    """Compare two numbers and return 1 if x>y, -1 if x<y and 0 otherwise"""
    if x>y: return 1
    elif x<y: return -1
    else: return 0

def pcompare(x,y,n_functions): #AW added
    """Perform a pareto comparison of two networks based on
    their different fitness

    Args:
        x,y (:class:`Network <phievo.Networks.mutation.Mutable_Network>`): the object to compare
        n_functions (int): the number of function taken into account

    Returns:
         the comparison of x & y (1 if x>y), 0 indicates that they are pareto equivalent
    """
    # case where some fitnesses are None
    if (None not in x.fitness) and (None in y.fitness):
        print("y is None")
        return -1
    elif (None not in y.fitness) and (None in  x.fitness):
        print("x is None")
        return 1
    elif (None in y.fitness) and (None in x.fitness):
        print("x and y is None")
        return 0
    # general case
    compare = [single_comparison(xf,yf) for xf,yf in zip(x.fitness,y.fitness)]
    SUM,LEN = sum(compare),len(compare)
    if SUM == LEN: return 1
    elif SUM == -LEN: return -1
    else: return 0

def compdist(x,y,n_functions):
    """Compute the distance between the fitness of x and y"""
    dist = 0
    for xf,yf in zip(x.fitness[:n_functions],y.fitness):
        if xf is None or yf is None: return None
        dist += (xf-yf)**2
    return sqrt(dist)

##########################################
### Class pareto_Population Definition ###
##########################################

class pareto_Population(Population):
    """Update the Population to manage a Pareto evolution

    Note that we dynamically change the fitness of the individuals to give
    them a list-like fitness.

    Attributes:
        nfunctions (int): number of functions taken into account by pareto
        rshare (float): parameter for the fitness sharing
    """
    def __init__(self,namefolder,nfunctions,rshare):
        Population.__init__(self,namefolder)  # needed when base class in another file it appears
        self.nfunctions=nfunctions
        self.rshare=rshare

        # update the Networks to behave like pareto guys
        for i in range(self.npopulation):
            self.genus[i].prank = 0
            self.genus[i].fitness = None
            self.genus[i].dlt_fitness = [0]*nfunctions

    def pop_sort(self,verbose = False): # AW added
        """Perform a pareto sorting of the population using the Goldberg algorithm.

        See Van Velhuizen and Lamont. Evol Computation. 8:125 (2000) for details
        To avoid having population dominated by 0,0 function assigns lowest rank
        to networks with this score.

        """
        if verbose:
            print("start_sort")
        random.shuffle(self.genus)
        curr_rank = 1
        to_sort = list(range(self.npopulation))
        zeros = []
        # Find those which are zero for all functions or have at least one None
        for index,individual in enumerate(self.genus):
            if individual.fitness is None:
                zeros.append(index)
                to_sort.remove(index)
                continue

        while to_sort:
            dominated = set()
            for i in to_sort:
                for j in to_sort:
                    if j >= i: pass
                    compij = pcompare(self.genus[i],self.genus[j],self.nfunctions)
                    if compij == 1:
                        dominated.add(i)
                    elif compij == -1:
                        dominated.add(j)

            for i in to_sort:
                if not i in dominated:
                    self.genus[i].prank = curr_rank
                    to_sort.remove(i)
            curr_rank = curr_rank+1

        for i in zeros:
            self.genus[i].prank = curr_rank

        self.genus.sort(key = lambda X: X.prank)
        self.pop_fitness_share()
        self.genus.sort(key = lambda X: X.prank)
        if verbose:
            for ind in self:
                print(ind.prank,ind.fitness)
            print("end_sort")

    def pop_fitness_share(self):
        """Use fitness sharing to increase the diversity of the population.

        That is, it augment the rank of inidividual to close from each other
        to promote diversity in the population.
        The implementation is a variant on the basic fitness sharing algorithm
        in section II of Cioppa et al. IEEE Trans. Evol Comp. 11:453

        """
        k = 0
        while( k < self.npopulation and self.genus[k].prank==1):
            for j in range(0,k):
                dist = compdist(self.genus[k],self.genus[j],self.nfunctions)
                if dist is None:
                    pass
                elif dist < self.rshare:
                    self.genus[k].prank+=(1-float(dist/self.rshare))/self.npopulation
                    self.genus[j].prank+=(1-float(dist/self.rshare))/self.npopulation
            k+=1

    def pop_print_pareto(self,f_pop,f_best):  #AW added
        """Write various information about population in files f_pop and

        Args:
            f_pop (str): short description of all individuals
            f_best (str): complete description of the first rank only
        """
        fp = open(f_pop,'w')
        for i0,ind in enumerate(self.genus):
            print(i0,ind.data_evolution,ind.prank,sep='\t',file=fp)
        fp.close()

        fp = open(f_best,'w')
        for i in range(0,self.npopulation):
            if self.genus[i].prank < 2:
                netnow = classes_eds2.print_Network(self.genus[i])
                fp.write(netnow)
        fp.close()

    def update_fitness(self,nnetwork,integration_result):
        """Update (in place) all the fitnesses and the corresponding dlt_fitness

        Args:
            nnetwork (int): the index of the network in the population
            integration_result (list): the output of compile_and_integrate

        """
        if integration_result:
            current_fitness = [float(integration_result[i]) for i in range(self.nfunctions)]
            self.genus[nnetwork].data_evolution = [data for data in integration_result]
        else: #catches the None fitness
            current_fitness = None
            self.genus[nnetwork].data_evolution = None

        #update dlt_fitness and handle the None fitness option
        if current_fitness and self.genus[nnetwork].fitness:
            self.genus[nnetwork].dlt_fitness = [current_fitness[i] - self.genus[nnetwork].fitness[i] for i in range(self.nfunctions)]
        elif current_fitness:
            self.genus[nnetwork].dlt_fitness = [9999]*self.nfunctions
        elif self.genus[nnetwork].fitness:
            self.genus[nnetwork].dlt_fitness = [-9999]*self.nfunctions
        else:
            self.genus[nnetwork].dlt_fitness = [0]*self.nfunctions
        self.genus[nnetwork].fitness = current_fitness

        #print('update_fitness_from_pareto',current_fitness)


#################################################
### Class pareto_thread_Population Definition ###
#################################################

class pareto_thread_Population(pareto_Population,thread_Population):
    """Update the pareto_Population class to allow threading

    Note, when looking for inherited method, python always choose
    the right most first (here pareto_Population).
    """
    pop_mutate_and_integrate = thread_Population.pop_mutate_and_integrate

if __name__ == "__main__":
    print(pcompare([999,0],[0.5,-0.5],2))
