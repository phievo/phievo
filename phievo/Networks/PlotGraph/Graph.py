import matplotlib.pyplot as plt
from . import Components
from .Layout import layout
from math import *
import numpy as np

class Graph:
    """
    Container of a directed graph. It contains mainly two types of objects: nodes and edges.
    """
    def __init__(self):
        """
        Initial set up of the conatiner.
        """
        self.nodes = {}
        self.edges = {}
        self.index_counter = 0
        ## An interaction contains all the edges between two given nodes
        self.interactions = {}
        self.node_size = 1
        self.grouped_interactions = {}
    def add_node(self,*argv,**kwargs):
        """
        Add a node to the graph.

        Args:
            argv (list(str)): Is handled if it contains only one element corresponding to the node label.
            kwargs (dict): The function handles only the keys **size**, **marker** that respectively correspond to the node's area and its shape. It can also deal with **label** if it was not defined in argv. The other keys are passed for latter use by the plotting function.
        Returns
            :class:`Networks.PlotGraph.Components.Node`: The node reference.

        """
        if len(argv) == 1:
            label = argv[0]
        else:

            try:
                label = kwargs["label"]
            except KeyError:
                print("A node cannot be added without a label.")
                return None
        if label in self.nodes.keys():
            print("%s is already in the graph, please choose another label."%(label))
            return None
        marker = kwargs.pop("marker","Circle")
        size = kwargs.pop("size",self.node_size)
        #import pdb;pdb.set_trace()
        try:
            newNode = getattr(Components,marker)(label=label,size=size,**kwargs)
        except AttributeError:
            print("Error: %s does not correspond to a node type."%marker)
            return None
        newNode.index = self.index_counter
        self.index_counter += 1
        self.nodes[label] = newNode
        return self.nodes[label]

    def add_edge(self,*argv,**kwargs):
        """
        Add an edge to the graph.

        Args:
            argv (list(str)): Is handled  if it contains only two elements corresponding to the edge's starting and ending nodes.
            kwargs (dict): The function handles only the keys **style**, **label** that respectively correspond to the edge's style and its label. It can also deal with **nodeFrom** and **nodeFrom** if it was not defined in argv. The other keys are passed for latter use by the plotting function.
        Returns
            :class:`Networks.PlotGraph.Components.Edge`:The edge reference.

        """
        if len(argv) == 2:
            nodeFrom = argv[0]
            nodeTo=argv[1]
        else:
            try:
                nodeFrom =  kwargs["nodeFrom"]
                nodeTo =  kwargs["nodeTo"]
            except KeyError:
                print("An edge cannot be added without a nodeFrom and a nodeTo.")
                return None
        style = kwargs.pop("style","Arrow")
        label = kwargs.pop("label",len(self.edges))


        if nodeFrom not in self.nodes.keys() and nodeFrom not  in self.edges.keys():
            self.add_node(nodeFrom)
        if nodeTo not in self.nodes.keys() and nodeTo not  in self.edges.keys():
            self.add_node(nodeTo)

        if nodeFrom in self.nodes.keys():
            nodeFrom = self.nodes[nodeFrom]
        else:
            nodeFrom = self.edges[nodeFrom]
            nodeTo.setReceiveEdge()

        if nodeTo in self.nodes.keys():
            nodeTo = self.nodes[nodeTo]
        else:
            nodeTo = self.edges[nodeTo]
            nodeTo.setReceiveEdge()
        try:
            newEdge = getattr(Components,style)(nodeFrom=nodeFrom,nodeTo=nodeTo,label=label,**kwargs)
        except AttributeError:
            print("Error: %s does not correspond to a edge style."%style)
            return None

        group_key = str(sorted([nodeFrom.index,nodeTo.index]))


        try:
            self.grouped_interactions[group_key].append(newEdge)
        except KeyError:
            self.grouped_interactions[group_key] = [newEdge]
        self.edges[label] = newEdge

        interactionLabel = "I_"+"-".join([str(xx) for xx in sorted([nodeFrom.index,nodeTo.index])]) # ex I_1-2
        try:
            interaction = self.interactions[interactionLabel]
        except KeyError:
            self.interactions[interactionLabel] = Components.Interaction(nodeFrom.index,nodeTo.index)
            interaction = self.interactions[interactionLabel]
        interaction.add_edge(self.edges[label])

        newEdge.index = self.index_counter
        self.index_counter += 1
        return self.edges[label]

    def set_node_size(self,size):
        """
        Homogenise the node area in the network.

        Args:
            size (float): Relative node area as compare to the default area.
        Returns:
            None
        """
        self.node_size = size

    def node_list(self):
        """ Generate a list of the node labels

        Returns:
            list(str): of the labels for the node contained in the graph
        """

        return [node.label for i,node in self.nodes.items()]

    def edge_list(self):
        """ Generate a list of the node edges

        Returns:
            list((str,str)): Each tuple in the list contains the starting and ending node labels.
        """
        edge_list = []
        for i,edge in self.edges.items():
            A = [edge.nodeFrom.label]
            B = [edge.nodeTo.label]

            if A[0] in self.edges.keys():
                A = [self.edges[A[0]].nodeFrom.label,self.edges[A[0]].nodeTo.label]
            if B[0] in self.edges.keys():
                B = [self.edges[B[0]].nodeFrom.label,self.edges[B[0]].nodeTo.label]
            edge_list += [(a,b) for a in A for b in B]
        return edge_list

    def layout(self,recursion=500):
        """
        Compute a layout for the node and set the node positions.
        """
        radius = sqrt(self.node_size/pi)
        #positions = layout(self.node_list(),self.edge_list(),radius,recursion,layout="graphviz")
        positions = layout(self.node_list(),self.edge_list(),radius=radius,layout="graphviz")
        for label,pos in positions.items():
            self.nodes[label].center = tuple(pos)
        return positions



    def draw(self,file=None,edgeLegend=False,display=True):
        """
        Draw the graph in a matplib framework. The node and edges are generated using patches.

        Args:
            file (str): Optional. When defined, the figure will be saved under the **file** name. Otherwise the program pops up a window with the graph.
        Returns:
            None

        """
        Scale = 2
        self.layout()
        xx = [xx.center[0] for key,xx in self.nodes.items()]
        yy = [yy.center[1] for key,yy in self.nodes.items()]
        maxX = np.max(xx) + 4*sqrt(self.node_size/pi)
        minX = np.min(xx) + 4*sqrt(self.node_size/pi)
        maxY = np.max(yy) + 4*sqrt(self.node_size/pi)
        minY = np.min(yy) + 4*sqrt(self.node_size/pi)
        width = (maxX-minX)
        height = (maxY-minY)

        plt.rc('font', family='sans-serif')
        #fig = plt.figure(figsize=(width,height))
        fig = plt.figure()
        ax = fig.gca()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        for label,node in self.nodes.items():
            patch = node.get_patch()
            ax.add_patch(patch)
            ax.annotate(label, node.center, color='black', weight='bold',fontsize=15, ha='center', va='center')

        new_iter = sorted(self.interactions.items(), key=lambda a: 1-a[1].receiveEdge)

        def add_iter(label,inter):
            """ Add an interaction to the plot """
            patches = inter.get_patches((Scale/30,Scale/30))
            for lab,patch in patches:
                ax.add_patch(patch)
                if edgeLegend:
                    ax.text(lab["pos"][0],lab["pos"][1], lab["text"],ha='center', va='center', color='black',bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))

        for label,inter in new_iter:
            if inter.isAuto: continue
            add_iter(label,inter)

        for label,inter in new_iter:
            if not inter.isAuto: continue
            add_iter(label,inter)


        plt.axis("scaled")
        if file:
            fig.savefig(file)
        elif display:
            plt.show()
        return fig
