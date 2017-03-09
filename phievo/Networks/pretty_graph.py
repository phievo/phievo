print("Execute pretty_graph.py")

import pydot
from .classes_eds2 import *
from .interaction import *

# independent of both above routines, goes from Network instance to Dot instance, and latter can be
# written out as .dot file and viewed with display

def pretty_graph(net):
    """ Take a Network instance and return a nice pydot graph with only species as nodes and
    interactions as edges.   Just outline of code here, many display items to fiddle
    Best Doc on pydot is http://dkbza.org/pydot/pydot.html   For attributes see
    http://www.graphviz.org/doc/info
    """

    # compose a short species label from selected types and their values,  Other types inferred from interactions
    #(NB TF might be dispensed with, and its activity recorded in arrow shape

    # complete list of 
    protein_types= ['TF', 'Kinase', 'Phosphatase', 'Ligand', 'Receptor']
   
    def short_label(species):
        label = ' '
        for k in species.types:
            if k == 'TF':
                label += 'TF=' + str(species.activity)
            elif k == 'Output':
                label += ' Out=' + str(species.n_put)
            elif k == 'Input':
                label += ' Input=' + str(species.n_put)
            elif protein_types.count(k):
                label += k
            else:
                pass
        return label

    # Collect display options for types of Nodes, use as eg pydot.Node(name, **attribute_dict['TF'] )
    attribute_dict = {}
    for tt in protein_types:
        attribute_dict[tt] = {}
    attribute_dict['TF']['shape'] = 'box'
    attribute_dict['Kinase']['shape'] = 'house'
    attribute_dict['Phosphatase']['shape'] = 'invhouse'
    attribute_dict['Ligand']['shape'] = 'invtriangle'
    attribute_dict['Receptor']['shape'] = 'invtrapezium'
    attribute_dict['Default'] = {}         # incase none of protein_types found

    # there are many different edge-arrow types available see the graphviz web site above.
    # unclear how to get key with shapes/colores, create subgraph? nodes only?
    
    # define graph with title
    graph = pydot.Dot(graph_type='digraph')
    try:
        graph.set_label('species graph: ' + net.title)
    except:
        graph.set_label('species graph: ')   # for compatability with old code or if change input net->net.graph
    graph.set_fontcolor('red')
    graph.set_labelloc('t')
    
    # create pydot nodes for all species nodes in net and store in dictionary
    map_species = {}
    for nn in net.list_types['Node']:
        if not nn.isinstance('Species'):
            continue
        name = nn.id + short_label(nn)
        # chose the protein type for node display, , order in protein_types list maters
        for kk in protein_types:
            if nn.types.count(kk):
                break
            else:
                kk = 'Default'
        # customize attributes.
        attribute_dict[kk]['style']="filled"
        attribute_dict[kk]['fillcolor']='0 0 1'
        map_species[nn] = pydot.Node(name, **attribute_dict[kk] )
        graph.add_node( map_species[nn] )

   
    # find the species connected by interactions and create pydot.Edges
    # NB pydot has no delete node methods, so can not eliminate just interactions, build new graph
    # and then go back and eliminate the TModules (phys object)
    
    for nn in net.list_types['Node']:
        if isinstance(nn, TModule):  
            for pre in net.graph.predecessors(nn):
                pre_species = net.graph.predecessors(pre)[0]  #only one species input to TFHill
                activity=pre.activity
                pre_dot_name = map_species[pre_species].get_name()
                for post in net.graph.successors(nn):
                    post_species = net.graph.successors(post)[0]  #only one species out from CorePromoter
                    post_dot_name = map_species[post_species].get_name()
                    ee = pydot.Edge(pre_dot_name, post_dot_name, label='Transcription' )
                    if hasattr(net,"fixed_activity_for_TF"):
                        	#if (net.fixed_activity_for_TF==0):
                                    if (activity==1):
                                        ee = pydot.Edge(pre_dot_name, post_dot_name, label='activation' )
                                    else:
                                        ee = pydot.Edge(pre_dot_name, post_dot_name, label='repression' )
                    graph.add_edge( ee )
                    
        elif isinstance(nn, Interaction):
            # these two classes of interaction nodes pre/post TModule and eliminated prev. loop
            if isinstance(nn, TFHill) or isinstance(nn, CorePromoter) :
                continue
            # For interaction with 2 in, 2 out species, draw substrate -> kinase -> product
            if isinstance(nn,Phosphorylation):
                name = nn.label
                [catalyst,listIn,listOut]=net.catal_data(nn)
                ee = pydot.Edge(map_species[listIn[0]].name,map_species[catalyst].get_name(), label=name )
                graph.add_edge( ee )
                ee = pydot.Edge(map_species[catalyst].name,map_species[listOut[0]].get_name(), label=name )
                graph.add_edge( ee )
                continue
            name = nn.label
            for pre in net.graph.predecessors(nn):
                pre_dot_name = map_species[pre].get_name()
                for post in net.graph.successors(nn):
                    post_dot_name = map_species[post].get_name()
                    ee = pydot.Edge(pre_dot_name, post_dot_name, label=name )  # color code interaction types?
                    graph.add_edge( ee )
                    
        else:
            pass
    # end of loop to eliminate interaction nodes.
    
    return graph


