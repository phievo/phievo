import pydot
#from phievo.initialization_code import *
import phievo.Networks.classes_eds2 as ceds2
from Immune.interaction_pMHC import *
from phievo.AnalysisTools import palette
# independent of both above routines, goes from Network instance to Dot instance, and latter can be
# written out as .dot file and viewed with display



version=pydot.__version__
flag=eval(version[0])

def compute_name(self):
    if (flag<1):
        return self.name
    else:
        return self.obj_dict['name']
setattr(pydot.Node,'compute_name',compute_name)#adds this function to the class Node


def pretty_graph(net,extended=False,layout=None):
    """ Take a Network instance and return a nice pydot graph with only species as nodes and
    interactions as edges.   Just outline of code here, many display items to fiddle
    Best Doc on pydot is http://dkbza.org/pydot/pydot.html   For attributes see
    http://www.graphviz.org/doc/info
    """
    net.__build_dict_types__()
    size=len(net.dict_types['Species'])+len(net.dict_types['pMHC'])+1
    colors=palette.color_generate(size)
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
    #try:
    #    graph.set_label('species graph: ' + net.title)
    #except:
    #    graph.set_label('species graph: ')   # for compatability with old code or if change input net->net.graph
    #graph.set_fontcolor('red')
    #graph.set_labelloc('t')
    
    # create pydot nodes for all species nodes in net and store in dictionary
    map_species = {}
    for nn in net.dict_types['Node']:
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
    for nn in net.dict_types['Node']:
        if isinstance(nn, ceds2.TModule):  
            for pre in net.graph.predecessors(nn):
                pre_species = net.graph.predecessors(pre)[0]  #only one species input to TFHill
                activity=pre.activity
                if (flag<1):
                    pre_dot_name = map_species[pre_species].name
                else:
                    pre_dot_name = map_species[pre_species].obj_dict['name']
                for post in net.graph.successors(nn):
                    post_species = net.graph.successors(post)[0]  #only one species out from CorePromoter
                    if (flag<1):
                        post_dot_name = map_species[post_species].name
                    else:
                        post_dot_name = map_species[post_species].obj_dict['name']
                    ee = pydot.Edge(pre_dot_name, post_dot_name, label='Transcription' )
                    if hasattr(net,"fixed_activity_for_TF"):
                            #if (net.fixed_activity_for_TF==0):
                                    if (activity==1):
                                        ee = pydot.Edge(pre_dot_name, post_dot_name,**attribute['activation'])
                                        
                                    else:
                                        ee = pydot.Edge(pre_dot_name, post_dot_name,**attribute['repression'] )
                                        ee.set_arrowhead('tee')
                    ee.set_arrowsize(2)
                    graph.add_edge( ee )
                    
        elif isinstance(nn, ceds2.Interaction):
            if isinstance(nn,Initial_Concentration):
                continue
            if isinstance(nn,Simple_Phosphorylation):
                name = nn.label
                [catalyst,listIn,listOut]=net.catal_data(nn)
                ee = pydot.Edge(map_species[listIn[0]].compute_name(),map_species[listOut[0]].compute_name(), label=map_species[catalyst[0]].compute_name())
                graph.add_edge( ee )
                continue
            if isinstance(nn,Simple_Dephosphorylation):
                name = nn.label
                [catalyst,listIn,listOut]=net.catal_data(nn)
                ee = pydot.Edge(map_species[listIn[0]].compute_name(),map_species[listOut[0]].compute_name(), label=map_species[catalyst[0]].compute_name())

                graph.add_edge( ee )
                continue
            if isinstance(nn,KPR_Binding):
                name = 'LR'
                Complex=net.graph.successors(nn)[0]
                PPI_components=net.graph.predecessors(nn)
                #nppi=pydot.Node('PPI',label='PPI',style='filled',fillcolor='0.6 1 1',width=0.6,height=0.8)
                nameLR='LR%i'%nn.int_id()
                nppi=pydot.Node(nameLR,label='LR',style='filled',fillcolor='0.6 1 1',radius=0.1,fontcolor='#FFFFFF',shape='circle',fixedsize='true', width=0.3,height=0.3)
                graph.add_node(nppi)
                for k in range(len(PPI_components)):
                    ee = pydot.Edge(map_species[PPI_components[k]].compute_name(),nameLR, **attribute['generic']  )
                    graph.add_edge( ee )
                ee = pydot.Edge(nameLR,map_species[Complex].compute_name(), **attribute['generic']  )
                graph.add_edge( ee )      
                continue
            if isinstance(nn,KPR_Unbinding):
                name = 'KPRU'
                Complex=net.graph.predecessors(nn)[0]
                [P1,P2]=net.graph.successors(nn)
                if P1.isinstance('Ligand'):
                    L=P1
                    R=P2
                else:
                    R=P1
                    L=P2
                #nppi=pydot.Node('PPI',label='PPI',style='filled',fillcolor='0.6 1 1',width=0.6,height=0.8)
                ee = pydot.Edge(map_species[Complex].compute_name(),map_species[L].compute_name(), label="1/&tau;",**attribute['generic']  )
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
    
    
    
def draw_pMHC(self,file=None,edgeLegend=False,extended=False,display=True,return_graph=False):
    plt.close()
    clear_output()
    graph=Networks.pretty_graph.pretty_graph(self,extended=extended)
    graph.write_png('current_graph.png')
    display(Image(filename='current_graph.png'))
setattr(ceds2,"draw",draw_pMHC)
