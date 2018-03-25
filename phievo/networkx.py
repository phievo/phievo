import networkx as nx
from inspect import getframeinfo, stack
"""
This class handles the compatibility of phievo with networkx>=2
"""


class MultiDiGraph(nx.MultiDiGraph):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)

    if float(nx.__version__)>=2:
        def list_nodes(self):        
            return list(self.nodes)

        def list_successors(self,node):
            return list(self.successors(node))

        def list_predecessors(self,node):
            return list(self.predecessors(node))
    else:
        def list_nodes(self):        
            return self.nodes()

        def list_successors(self,node):
            return self.successors(node)

        def list_predecessors(self,node):
            return self.predecessors(node)

 
