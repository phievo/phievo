"""
The lovelyGraph modules contains a set of utilities to plot a network.
It uses the homemade package PlotGraph_.
"""

import phievo.Networks.PlotGraph as PlotGraph
from phievo.initialization_code import *
from phievo.Networks.classes_eds2 import *
from phievo.AnalysisTools import palette
from phievo.Networks.interaction import *

## Node Parameters
node_types= ['Output','Input','Ligand','Receptor',"PPI","TModule"]
global_node_attribute = {}
for tt in node_types:global_node_attribute[tt] = {}
global_node_attribute["Default"] = {}
global_node_attribute['Output']['marker'] = 'TriangleUp'
global_node_attribute['Input']['marker'] = 'TriangleDown'
global_node_attribute['Ligand']['marker'] = 'HouseUp'
global_node_attribute['Receptor']['marker'] = 'HouseDown'
global_node_attribute['TModule']= {'marker':"Square","color":"#FE9A2E"}
global_node_attribute['PPI']= {'marker':"RoundedRectangle","facecolor":"#ffffff","ls":"-","lw":1.5,"edgecolor":"#1c2833","size":0.5}
global_node_attribute['Binding']= {'marker':"RoundedRectangle","facecolor":"#E5E5E5","ls":"-","lw":1.5,"edgecolor":"#1c2833","size":0.5}
#global_node_attribute['PPI']= {}#{'marker':"Square","edgecolor":"#1c2833","size":0.5}
global_node_attribute['Default'] = {'marker':"Circle","color":"#808080"}

for tt in global_node_attribute.keys():
        global_node_attribute[tt].setdefault("size",1)
        global_node_attribute[tt].setdefault("lw",0)


## Edge Parameters
global_edge_attribute={}
global_edge_attribute['repression']={"linewidth":2,"color":'#B00909','arrowstyle': '-|'}
global_edge_attribute['activation']={"linewidth":2,"color":'#09B01C','arrowstyle': '-|>'}
global_edge_attribute['Default']={"linewidth":2,"color":'#000000',"arrowstyle":"->",}
global_edge_attribute['Degradation']={"linewidth":2,"color":'#000000',"arrowstyle":"->","linestyle":':'}
global_edge_attribute['PPI']={"linewidth":2,"color":'#666666',"arrowstyle":"->"}
global_edge_attribute['Binding']={"linewidth":3,"color":'#000000',"arrowstyle":"->"}
global_edge_attribute['Phospho']={"linewidth":2,"color":'#000000',"arrowstyle":"->"}

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

def gettype(node,type_list):

    for tt in type_list:
        if tt in node.types:
            return tt
    return "Default"

def produce_species_name(node_species):
        return '$S_{%i}$'%node_species.int_id()
def produce_TModule_name(node_species):
        return '$TM_{%i}$'%node_species.int_id()
def produce_PPI_name(node_PPI):
        return '$PPI_{%i}$'%node_PPI.int_id()

def produce_TFHill_name(node_reac):
        return '$H_{%i}$'%node_reac.int_id()
def produce_CorePromoter_name(node_reac):
        return '$CP_{%i}$'%node_reac.int_id()
def produce_Degradation_name(node_reac):
        return '$D_{%i}$'%node_reac.int_id()
def produce_Phospho_name(node_reac,cat=False):
        if cat:
                return '$\star_{%i}$'%node_reac.int_id()
        else:
                return '$Ph_{%i}$'%node_reac.int_id()


def pretty_graph(net,extended=True):
        """
        Creates a ready-to-plot graph object from a network.

        Args:
            net (:class:`Mutable_Network <phievo.Networks.mutation.Mutable_Network>`)
        Return:
            returns a :class:`PlotGraph graph <phievo.Networks.PlotGraph.Graph.Graph>`

        """
        net.__build_dict_types__()
        size=len(net.dict_types['Species'])
        colors=palette.color_generate(size)
        graph = PlotGraph.Graph()
        # create pydot nodes for all species nodes in net and store in dictionary
        map_species = {}

        ## Add the nodes to the graph
        for nn in net.dict_types['Node']:
                if not nn.isinstance('Species'):
                        continue
                name = produce_species_name(nn)
                node_type = gettype(nn,node_types)
                node_attributes = dict(global_node_attribute[node_type])
                node_attributes["color"] = colors[nn.int_id()]
                graph.add_node(name,**node_attributes)


        for nn in net.dict_types['Node']:
                if isinstance(nn, TModule):
                        succ = net.graph.successors(net.graph.successors(nn)[0])[0]
                        succ_name = produce_species_name(succ)
                        if extended:
                                ## include TModule to the graph
                                post_species = nn
                                post_name = produce_TModule_name(nn)
                                graph.add_node(post_name,**global_node_attribute['TModule'])
                                param = dict(**global_edge_attribute['Default'])
                                param["label"] = produce_CorePromoter_name(net.graph.successors(nn)[0])
                                ee = graph.add_edge(post_name, succ_name,**param)
                        else:
                                post_species = succ
                                post_name = succ_name
                        for pre in net.graph.predecessors(nn):
                                pre_species = net.graph.predecessors(pre)[0] #only one species input to TFHill
                                pre_name = produce_species_name(pre_species)
                                if hasattr(net,"fixed_activity_for_TF"):
                                        if (pre.activity==1):
                                                param = dict(**global_edge_attribute['activation'])
                                                param["label"] = produce_TFHill_name(pre)
                                                ee = graph.add_edge(pre_name, post_name,**param)
                                        else:
                                                param = dict(**global_edge_attribute['repression'])
                                                param["label"] = produce_TFHill_name(pre)
                                                ee = graph.add_edge(pre_name, post_name,**param )
                                else:
                                        param = dict(**global_edge_attribute['Default'])
                                        param["label"] = produce_TFHill_name(pre)
                                        ee = graph.add_edge(pre_name,post_name,**param)


                elif isinstance(nn,Interaction):
                        if isinstance(nn, TFHill) or isinstance(nn, CorePromoter):
                                continue

                        elif isinstance(nn,Phosphorylation):
                                namePhospho = nn.label
                                [cataList,listIn,listOut]=net.catal_data(nn)
                                param = dict(**global_edge_attribute['Phospho'])
                                param["label"] = produce_Phospho_name(nn)
                                graph.add_edge(produce_species_name(listIn[0]),produce_species_name(listOut[0]), **param)
                                param["arrowstyle"] = "-|>"

                                param["label"] = produce_Phospho_name(nn,cat=True)
                                graph.add_edge(produce_species_name(cataList[0]),produce_Phospho_name(nn),**param)

                        elif isinstance(nn,PPI):
                                namePPI = produce_PPI_name(nn)
                                PPI_components = net.graph.predecessors(nn)
                                PPI_complex=net.graph.successors(nn)[0]
                                graph.add_node(namePPI,**global_node_attribute["PPI"])
                                param = dict(**global_edge_attribute['PPI'])
                                param["label"] = produce_PPI_name(nn)
                                graph.add_edge(namePPI,produce_species_name(PPI_complex),**param)
                                for compo in PPI_components:
                                        param["label"] = ""
                                        param["arrowstyle"]="-"
                                        graph.add_edge(produce_species_name(compo),namePPI,**param)

                        elif isinstance(nn,Degradation):
                                nameDegrad = produce_Degradation_name(nn)
                                shreder = net.graph.predecessors(nn)[0]
                                degraded = net.graph.successors(nn)[0]
                                graph.add_edge(produce_species_name(shreder),produce_species_name(degraded), label=nameDegrad,**global_edge_attribute["Degradation"])

                        elif nn.label == "Simple_Phosphorylation":
                                namePhospho = nn.label
                                [listCata,listIn,listOut]=net.catal_data(nn)
                                param = dict(**global_edge_attribute['Phospho'])
                                param["color"] = "#B42B3C"
                                param["label"] = produce_Phospho_name(nn)
                                graph.add_edge(produce_species_name(listIn[0]),produce_species_name(listOut[0]), **param)
                                param["arrowstyle"] = "-|>"
                                param["label"] = produce_Phospho_name(nn,cat=True)
                                graph.add_edge(produce_species_name(listCata[0]),produce_Phospho_name(nn),**param)

                        elif nn.label == "Simple_Dephosphorylation":
                                namePhospho = nn.label
                                [listCata,listIn,listOut]=net.catal_data(nn)
                                param = dict(**global_edge_attribute['Phospho'])
                                param["color"] = "#2B39B4"
                                param["label"] = produce_Phospho_name(nn)
                                graph.add_edge(produce_species_name(listIn[0]),produce_species_name(listOut[0]), **param)
                                param["arrowstyle"] = "-|>"
                                param["label"] = produce_Phospho_name(nn,cat=True)
                                graph.add_edge(produce_species_name(listCata[0]),produce_Phospho_name(nn),**param)

                        elif nn.label == "KPR_Binding":
                                binding_name = 'LR'
                                binding_components = net.graph.predecessors(nn)
                                binding_complex=net.graph.successors(nn)[0]
                                graph.add_node(binding_name,**global_node_attribute["Binding"])
                                param = dict(**global_edge_attribute['Binding'])
                                param["label"] = binding_name
                                graph.add_edge(binding_name,produce_species_name(binding_complex),**param)
                                for compo in binding_components:
                                        param["label"] = ""
                                        param["arrowstyle"]="-"
                                        graph.add_edge(produce_species_name(compo),binding_name,**param)

                        elif nn.label == "Initial_Concentration":
                                name_conc = nn.id
                                spec = net.graph.successors(nn)[0]
                                graph.add_node(name_conc,**global_node_attribute["PPI"])
                                param = dict(**global_edge_attribute['PPI'])
                                param["label"] = ""
                                param["arrowstyle"]="-"
                                graph.add_edge(produce_species_name(spec),name_conc,**param)

                        elif nn.label == "KPR_Unbinding":
                                continue

                        else:
                                inputs = net.graph.predecessors(nn)
                                outputs = net.graph.successors(nn)
                                for inp in inputs:
                                        for out in outputs:
                                                graph.add_edge(produce_species_name(inp),produce_species_name(out), label=nn.id,**global_edge_attribute["Default"])

        return graph
