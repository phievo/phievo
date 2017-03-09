from classes_eds2 import *
from mutation import *
from interaction import *
data=[]
data.append(0.559337)
data.append(0.002013)
data.append(0.032570)
data.append(3.214800)
data.append(0.220346)
import random
g=random.Random(42)
L=Mutable_Network(g)
L.Cseed=45059
s0=Species()
s0.activity=1
s0.degradation=0.448463581195
s0.removable=False
s0.n_put=0
s0.mutable=1
s0.order=0
s0.types=['Species', 'Degradable', 'TF', 'Input', 'Complexable', 'Kinase']
L.add_Node(s0)
s1=Species()
s1.removable=True
s1.activity=0
s1.mutable=1
s1.degradation=0.57329248605
s1.order=1
s1.types=['Species', 'Degradable', 'Phosphorylable', 'TF', 'Complexable','Output']
s1.n_put=0
L.add_Node(s1)
n2=CorePromoter()
n2.output=['Species']
n2.delay=0
n2.removable=True
n2.mutable=1
n2.input=['TModule']
n2.order=2
L.add_Node(n2)
n3=TModule()
n3.basal=0.0
n3.rate=0.852808475279
n3.removable=True
n3.mutable=1
n3.order=3
L.add_Node(n3)
s2=Species()
s2.degradation=0.00956135556962
s2.removable=True
s2.activity=1
s2.mutable=1
s2.order=5
s2.types=['Species', 'Degradable', 'Complex', 'Phosphorylable', 'TF', 'Kinase']
L.add_Node(s2)
n5=PPI()
n5.output=['Species']
n5.disassociation=0.742660524546
n5.removable=True
n5.mutable=1
n5.input=['Complexable', 'Complexable']
n5.order=6
n5.association=0.946127875968
L.add_Node(n5)
s3=Species()
s3.activity=1
s3.degradation=0.951293189816
s3.removable=True
s3.mutable=1
s3.order=7
s3.types=['Species', 'Degradable', 'Complex', 'Phosphorylable', 'TF']
L.add_Node(s3)
n7=PPI()
n7.mutable=1
n7.output=['Species']
n7.removable=True
n7.disassociation=0.824945598857
n7.input=['Complexable', 'Complexable']
n7.order=8
n7.association=0*0.76722012671
L.add_Node(n7)
L.activator_required=0
L.fixed_activity_for_TF=1
L.graph.add_edge(s0,n5)
L.graph.add_edge(s1,n5)
L.graph.add_edge(s1,n7)
L.graph.add_edge(s1,n7)
L.graph.add_edge(n2,s1)
L.graph.add_edge(n3,n2)
L.graph.add_edge(n5,s2)
L.graph.add_edge(n7,s3)

L.remove_Node(n7)
L.clean_Nodes()
