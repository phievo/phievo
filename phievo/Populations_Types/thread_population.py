"""
Expand the population class of evolution_gillespie to allow threading
"""
from .evolution_gillespie import Population
import threading
import time

class thread_Population(Population):
    """Update the Population class to allow threading
    Only change the pop_mutate_and_integrate method
    """
    def __init__(self,namefolder):
        Population.__init__(self,namefolder)  # needed when base class in another file it appears

    def pop_mutate_and_integrate(self,initial,first_mutated,last_mutated,prmt,net_stat):
        """ Recompute the fitness for half the population and mutate/compute the fitness for the rest.
        Save all the data in net_stat
        
        Args:
            initial (int): index of the first individual in population
            first_mutated (int): index of the first mutated individual in population
            last_mutated (int): index of the last mutated individual in population
            prmt (dict): the inits parameters for integration
            net_stat (NetworkStat): to store the population data
        
        Returns:
            None: in place modification
        """
        list_thread=[]
        self.n_mutations=0
        for nnetwork in range(initial,first_mutated):
            list_thread.append(threading.Thread(None,self.genus_mutate_and_integrate,None,(prmt, nnetwork,0,)))
        for nnetwork in range(first_mutated,last_mutated):
            list_thread.append(threading.Thread(None,self.genus_mutate_and_integrate,None,(prmt, nnetwork,1,)))
        
        #now we send the thread to the procs
        for index in list_thread: index.start()
        for index in list_thread: index.join()

        for individual in self.genus:
            net_stat.add_net(individual)
