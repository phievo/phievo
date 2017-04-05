import random as rd
from matplotlib import pyplot as plt
import numpy as np
import warnings
#dist = lambda A,B,pos: ((pos[A][0]-pos[B][0])**2+(pos[A][1]-pos[B][1])**2)**.5
#reg = lambda x: max(x,-1) if x < 0 else min(x,1)

def layout(node_list,interaction_list,radius=1,layout="graphviz"):
    """
    Use networkx layout function to compute the node centers

    Args:
        node_list (list): List of all the nodes in the nework
        interaction_list (list): List of tuple describing the nodes in interaction
        radius (float): Order of magnitude for a node radius. used to scale the minimal distance.
        layout (str): Use a networkx layout. Choose between:
                       - circular
                       - spring
                       - shell
                       - random
                       - spectral
                       - circular
                       - fruchterman_reingold
                       - pygraphviz

    Return:
        dict: indexed by nodes names and containing their (x,y) position (for use with draw_networkx pos argument typically)
    """
    pos = None
    import networkx as nx
    G = nx.Graph()
    for node in node_list:
        G.add_node(node)
    for edge in interaction_list:
        G.add_edge(edge[0],edge[1])

    if layout=="graphviz":
        try:
            ## Graphviz does not return the positions as numpy arrays
            pos = {kk:np.array(xx) for kk,xx in nx.nx_agraph.graphviz_layout(G).items()}
        except ImportError:
            warnings.warn('pygraphviz is not correctly installed - using spring_layout instead')
            pos = nx.spring_layout(G)
    else:
        pos = getattr(nx,layout+"_layout")

        # except AttributeError:
        #     print("%s is not a compatible layout. Please try with one of the following layouts:
        #     \t - circular
        #     \t - spring
        #     \t - shell
        #     \t - random
        #     \t - spectral
        #     \t - fruchterman_reingold")
            ## Scaling according to the minimal distance between two nodes
    minDist = 100000
    keys = list(pos.keys())
    ## Need to run over the indexes in order not to count twice the same distance
    for iA in range(len(node_list)):
        for iB in range(iA+1,len(node_list)):
            minDist = min(np.linalg.norm(pos[keys[iB]] - pos[keys[iA]]),minDist)
    for iA,arr in pos.items():
        pos[iA][0]*=(7*radius/minDist)
        pos[iA][1]*=(7*radius/minDist)
    return pos
