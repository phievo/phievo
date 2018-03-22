import networkx as nx
from inspect import getframeinfo, stack
"""
This class handles the compatibility of phievo with networkx>=2
"""


class MultiDiGraph(nx.MultiDiGraph):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)

    def list_nodes(self):        
        return list(self.nodes)

    def list_successors(self,node):
        return list(self.successors(node))

    def list_predecessors(self,node):
        return list(self.predecessors(node))
    

 
