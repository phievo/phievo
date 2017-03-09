"""
This module gather three classes used to produce statistics about the
networks in a population.
"""

class AnyStat(object):
    """Flexible interface to accumulate and print various statistics
    
    To add a new statistic have to change all functions consistently
    
    Attributes:
        label (str): to tag the statistic you want to construct
        data (dict): gather the useful information
    
    Main Methods:
        add_value: add a value to the statistic
        stat2string: normalize things and return a string descriptor
    """
    def __init__(self, label):
        self.label = label
        self.data  = self.init_data()

    def init_data(self):
        """Initialize the data dictionary"""
        data  = { 'sum':0, 'min':1.E100, 'max':-1.E100, 'count':0}
        return data

    def add_value(self, val):
        """Add a value to the statistic
        
        Args:
            val (float): the value to be added
            
        Returns:
            None: in place modification
        """
        if val: #skip the None case
            self.data['sum'] += val
            self.data['min'] = min( val, self.data['min'] )
            self.data['max'] = max( val, self.data['max'] )
            self.data['count'] += 1

    def stat2string(self):
        """Return the useful statistic (count,mean,min and max) as a string
        Args:
            -
        Returns:
            str: statistic output
        """
        self.data['avr'] = float(self.data['sum'])/max(1, self.data['count'])
        message = "{0}: Count= {1[count]}, Avr= {1[avr]:.2e}, Min= {1[min]:.2e}, Max= {1[max]:.2e}"
        return message.format(self.label,self.data)

class NetworkStat(object):
    """Composition of various statistic about a Network population
    
    Attributes:
        cmd_label (dict): {label : command to a Network instance}
        stat_box (dict): {label : corresponding Anystat instance}
    
    Main Methods:
        add_net: add a Network to the whole statistic
        output: return the string description
    """
    def __init__(self, cmd_label):
        self.cmd_label = cmd_label
        self.stat_box = {lbl:AnyStat(lbl) for lbl in cmd_label}

    def add_net(self, net):
        """Extract the information from net and add them to stat_box statistics
        
        Args:
            net (Network):
        
        Returns:
            None: in place modification
        """
        net.write_id()
        #import pdb; pdb.set_trace()
        for lbl, cmd in self.cmd_label.items():
            val = cmd(net)
            if isinstance(val, list): val = len(val)
            self.stat_box[lbl].add_value(val)

    def output(self):
        """Return the full statistic computation"""
        for stat in self.stat_box.values():
            print(stat.stat2string())

def isdifferent(a, b, abs_err, rel_err):
    """Estimate if a and b are similar or not up to absolute and relative errors
    
    Args:
        a,b (float or list): the object to compare
        abs_err (float): absolute difference allowed
        rel_err (float): relative difference allowed
    
    Returns:
        bool: if a and b are very close
    """
    try: #handle the case where a and b are lists
        if type(a) is list: a = sum(a)
        if type(b) is list: b = sum(b)
    except TypeError:
        return False
    
    if(abs(a - b) > abs_err + rel_err*(abs(a) + abs(b))/2 ):
        return True
    else:
        return False

class GenusStat(object):
    """Compute statistic about a population of evolving networks
    
    Attributes:
        abs_err,rel_err (float): error margins to decide when fitness changes are significant
        f_incr,f_decr,f_same (int): number of Networks whom fitness have increase/decrease/unchange
        count_topol (dict): count the number of different topol. in the population
    """
    abs_err = 1.e-2
    rel_err = 1.e-3
    
    def __init__(self):
        self.f_incr = 0
        self.f_decr = 0
        self.f_same = 0
        self.count_topol = {}

    def process_sorted_genus(self, pop):
        """extract stats from the population instance and collect in GenusStat"""
        best_fitness = pop.genus[0].fitness
        if not isdifferent(best_fitness, pop.best_fitness, GenusStat.abs_err, GenusStat.rel_err):
            pop.best_fitness_counter += 1
        pop.best_fitness = best_fitness

        for net in pop.genus:
            # counts of networks that incr/decr fitness after mutation
            try: #handle the case where fitness is a list
                dlt = sum(net.dlt_fitness)
            except TypeError:
                dlt = net.dlt_fitness
            
            err = GenusStat.abs_err + GenusStat.rel_err*abs(dlt)
            if dlt > err:
                self.f_incr += 1
            elif dlt < -err:
                self.f_decr += 1
            else:
                self.f_same += 1

            # enumerate distinct topologies.
            hh = net.hash_topology
            self.count_topol[hh] = self.count_topol.get(hh, 0) + 1

    def output(self):
        """Print the full statistic"""
        max_topol = max( self.count_topol.values() )
        n_topol = len( list(self.count_topol.keys()) )

        print('Number of networks with incr fitness=', self.f_incr, ' decr. fitness=', self.f_decr)
        print('Number of topol distinct networks=', n_topol, ' max number with same topology=', max_topol)
