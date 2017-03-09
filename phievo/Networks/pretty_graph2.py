print("Execute pretty_graph2.py")

from phievo.initialization_code import *
from phievo.Networks.classes_eds2 import *
from phievo.Networks.palette import *
from phievo.Networks.interaction import *
import pydot

exec('from '+name_interaction+' import *')
# independent of both above routines, goes from Network instance to Dot instance, and latter can be
# written out as .dot file and viewed with display

version=pydot.__version__
flag=int(version[0])

def compute_name(self):
    if (flag<1):
    	return self.name
    else:
        return self.obj_dict['name']
setattr(pydot.Node,'compute_name',compute_name)#adds this function to the class Node


def pretty_graph(net):
    """ Take a Network instance and return a nice pydot graph with only species as nodes and
    interactions as edges.   Just outline of code here, many display items to fiddle
    Best Doc on pydot is http://dkbza.org/pydot/pydot.html
    For attributes see http://www.graphviz.org/doc/info
    """
    net.__build_list_types__()
    
    size=len(net.list_types['Species'])
    colors=color_generate(size)
    # compose a short species label from selected types and their values,  Other types inferred from interactions
    #(NB TF might be dispensed with, and its activity recorded in arrow shape

    # complete list of
    protein_types= ['Output','Input','Ligand','Receptor']

    def short_label(species):
        label = ' '
        for k in species.types:
            if k == 'Output':
                label += ' Output'
            elif k == 'Input':
                label += ' Input'
            elif protein_types.count(k):
                label += k
            else:
                pass
        return label

    # Collect display options for types of Nodes, use as eg pydot.Node(name, **attribute_dict['TF'] )
    attribute_dict = {}
    for tt in protein_types:
        attribute_dict[tt] = {}
    attribute_dict['Output']['shape'] = 'triangle'
    attribute_dict['Input']['shape'] = 'invtriangle'
    attribute_dict['Ligand']['shape'] = 'house'
    attribute_dict['Receptor']['shape'] = 'invhouse'
    attribute_dict['Default'] = {}         # incase none of protein_types found

    # there are many different edge-arrow types available see the graphviz web site above.
    # unclear how to get key with shapes/colores, create subgraph? nodes only?

    # define graph with title
    graph = pydot.Dot(graph_type='digraph')

    # create pydot nodes for all species nodes in net and store in dictionary
    map_species = {}

    for nn in net.list_types['Node']:
        if not nn.isinstance('Species'):
            continue
        #name = '%i'%nn.int_id() + short_label(nn)
        name = '%i'%nn.int_id()
        # chose the protein type for node display, , order in protein_types list maters
        for kk in protein_types:
            if nn.types.count(kk):
                break
            else:
                kk = 'Default'
        # customize attributes.                    
        attribute_dict[kk]['style']="filled"
        attribute_dict[kk]['fillcolor']=colors[nn.int_id()]
        attribute_dict[kk]['fontsize']=25
        attribute_dict[kk]['width']=0.6
        attribute_dict[kk]['height']=0.8
        attribute_dict[kk]['fontcolor']='#FFFFFF'
        attribute_dict[kk]['fixedsize']='true'
        map_species[nn] = pydot.Node(name, **attribute_dict[kk] )

        graph.add_node( map_species[nn] )


    # find the species connected by interactions and create pydot.Edges
    # NB pydot has no delete node methods, so can not eliminate just interactions, build new graph
    # and then go back and eliminate the TModules (phys object)
    attribute={}
    attribute['repression']={}
    attribute['repression']['color']='0 1 1'
    attribute['activation']={}
    attribute['activation']['color']='0.3 1 0.5'
    attribute['generic']={}
    attribute['generic']['color']='0.6 1 1'
    for nn in net.list_types['Node']:
        if isinstance(nn, TModule):
            for pre in net.graph.predecessors(nn):
                pre_species = net.graph.predecessors(pre)[0]  
                activity=pre.activity
                if (flag<1):
                	pre_dot_name = map_species[pre_species].name
                else:
                	pre_dot_name = map_species[pre_species].obj_dict['name']
                for post in net.graph.successors(nn):
                    post_species = net.graph.successors(post)[0]  
                    if (flag<1):
                    	post_dot_name = map_species[post_species].name
                    else:
                    	post_dot_name = map_species[post_species].obj_dict['name']
                    ee = pydot.Edge(pre_dot_name, post_dot_name, label='Transcription' )
                    if hasattr(net,"fixed_activity_for_TF"):
                                    if (activity==1):
                                        ee = pydot.Edge(pre_dot_name, post_dot_name,**attribute['activation'])
                                    else:
                                        ee = pydot.Edge(pre_dot_name, post_dot_name,**attribute['repression'] )
                                        ee.set_arrowhead('tee')
                    ee.set_arrowsize(2)
                    graph.add_edge(ee)

        elif isinstance(nn, Interaction):
            # these two classes of interaction nodes pre/post TModule and eliminated prev. loop
            if isinstance(nn, TFHill) or isinstance(nn, CorePromoter):
                continue
            # For interaction with 2 in, 2 out species, draw substrate -> kinase -> product
            if isinstance(nn,Phosphorylation):
                name = nn.label
                [catalyst,listIn,listOut]=net.catal_data(nn)
                ee = pydot.Edge(map_species[listIn[0]].compute_name(),map_species[catalyst].compute_name(), label=name )
                graph.add_edge( ee )
                ee = pydot.Edge(map_species[catalyst].compute_name(),map_species[listOut[0]].compute_name(), label=name )
                graph.add_edge( ee )
                continue
            if isinstance(nn,PPI):
                name = 'PPI'
                Complex=net.graph.successors(nn)[0]
                PPI_components=net.graph.predecessors(nn)
                namePPI='PPI%i'%nn.int_id()
                nppi=pydot.Node(namePPI,label='ppi',style='filled',fillcolor='0.6 1 1',radius=0.1,fontcolor='#FFFFFF',shape='circle',fixedsize='true', width=0.3,height=0.3)
                graph.add_node(nppi)
                for k in range(len(PPI_components)):
                	ee = pydot.Edge(map_species[PPI_components[k]].compute_name(),namePPI, **attribute['generic']  )
                	graph.add_edge( ee )
                #ee = pydot.Edge(map_species[PPI_components[1]].compute_name(),namePPI, **attribute['generic']  )
                #graph.add_edge( ee )
                ee = pydot.Edge(namePPI,map_species[Complex].compute_name(), **attribute['generic']  )
                graph.add_edge( ee )
                #ee = pydot.Edge(map_species[PPI_components[0]].name,map_species[Complex].name, label=name,**attribute['generic']  )
                #ee.set_fontcolor('0.6 1 1')
                #graph.add_edge( ee )
                #ee = pydot.Edge(map_species[PPI_components[1]].name,map_species[Complex].name, label=name ,**attribute['generic']  )
                #ee.set_fontcolor('0.6 1 1')
                #graph.add_edge( ee )
                continue
            if isinstance(nn,Degradation):
                name = 'Degradation'
                shreder = net.graph.predecessors(nn)[0]
                degraded = net.graph.successors(nn)[0]
                ee = pydot.Edge(map_species[shreder].compute_name(),map_species[degraded].compute_name(), label=name)
                graph.add_edge( ee )
                continue
            name = nn.label
            print(name)
            for pre in net.graph.predecessors(nn):
                pre_dot_name = map_species[pre].compute_name
                for post in net.graph.successors(nn):
                    post_dot_name = map_species[post].compute_name()
                    ee = pydot.Edge(pre_dot_name, post_dot_name,**attribute['generic'] )  # color code interaction types?
                    graph.add_edge( ee )

        else:
            pass
    # end of loop to eliminate interaction nodes.
    return graph
