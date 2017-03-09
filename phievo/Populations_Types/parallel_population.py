"""
Expand the population class of evolution_gillespie to allow parallelization
"""
from phievo.Populations_Types.evolution_gillespie import Population
import pypar
import gc

class parallel_Population(Population):
    """ Update the Population class to allow parallelization
    Modify the pop_mutate_and_integrate method and add a
    multi_proc_mutate_and_integrate one
    """
    def __init__(self,namefolder):
        Population.__init__(self,namefolder)  # needed when base class in another file it appears
        self.numproc=pypar.size()
        
    def multi_proc_mutate_and_integrate(self,prmt,mutation):
        """subroutine to send jobs to multiple processors using pypar
        
        It send one network - and not the complete population- to each processor
        and explicitly we take care of synchronization.
        
        Args:
            prmt (dict): the inits parameters for integration
            mutation (list): tuple (id,mut) indicating the mutation flag of integration for each network
        
        Returns:
            int: the total number of mutations
        """
        numproc=self.numproc #computes number of available proc for integration : the  proc 0 is used as master proc
        l=len(mutation)
        n_mut=0
        for index_job in range(l):
            nproc=1+index_job%(numproc-1)#computes the proc number where to send the job, we start at 1, proc 0 is master proc
            args={'net':self.genus[index_job],'prmt':prmt,'nnetwork':index_job,'tgeneration':self.tgeneration,'mutation':mutation[index_job]}
            pypar.send(('net.mutate_and_integrate(prmt,nnetwork,tgeneration,mutation)',args),nproc)#send integration job to the selected proc
            if ((index_job+1)%(numproc-1)==0):
                pypar.barrier()# every numproc jobs sent, we wait for all other processors to finish their job to synchronize
                results=[pypar.receive(worker) for worker in range(1,numproc)]#receives results from all processors
                for i in results:
                    n_mut+=i[0]#updates number of mutations
                    self.genus[i[1]]=i[2]#updates mutated network
                    self.update_fitness(i[1],i[3])#updates fitness values
    
         #at that point, we may have jobs still running on some subset of the procs, but we only take the results from the last working processors
        if (l%(numproc-1)>0):
            for i in range(l%(numproc-1),numproc-1):
                pypar.send(('0',{}),i+1)#send dummy jobs to the still processors for synchornization purpose
            pypar.barrier()#synchronize the proc
            results=[pypar.receive(worker) for worker in range(1,numproc)]
            #only take the results we are interested in
            for i in range(l%(numproc-1)):
                n_mut+=results[i][0]
                self.genus[results[i][1]]=results[i][2]      
                self.update_fitness(results[i][1],results[i][3])#updates fitness values
        # Shut down workers
        #for worker in range(1,numproc):
        #    pypar.send(SystemExit(),worker)
        return n_mut


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
        mutation=[]
        #creates table containing the list of networks to mutate
        for nnetwork in range(0,first_mutated):
            mutation.append(0)
        for nnetwork in range(first_mutated,last_mutated):
            mutation.append(1)
        #now we compute the mutated network and the number of mutations
        gc.collect()
        list=gc.get_objects()
        print("Number of garbage collected stuff : %i"%(len(list)))
        print("Number of uncollectibla stuff : %i"%(len(gc.garbage)))
        
        self.n_mutations = self.multi_proc_mutate_and_integrate(prmt,mutation)
        
        for individual in self.genus:
            net_stat.add_net(individual)
