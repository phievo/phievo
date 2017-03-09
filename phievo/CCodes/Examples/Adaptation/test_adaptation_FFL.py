from classes_eds2 import *
from mutation import *
from interaction import *
data=[]
data.append(0.309183)
data.append(0.309183)
data.append(0.080670)
data.append(0.003967)
data.append(0.047573)
data.append(4.050106)
data.append(0.192666)
import random
g=random.Random(42)
L=Mutable_Network(g)

L.activator_required=0
L.fixed_activity_for_TF=1
s0=Species()
s0.activity=0
s0.n_put=0
s0.degradation=0.867257066362
s0.order=0
s0.types=['Species', 'Degradable', 'TF', 'Input']
L.add_Node(s0)
s1=Species()
s1.activity=1
s1.n_put=0
s1.degradation=0.183929466832
s1.order=1
s1.types=['Species', 'Degradable', 'Phosphorylable', 'TF', 'Complexable', 'Output']
L.add_Node(s1)
n2=CorePromoter()
n2.order=2
n2.delay=16
n2.output=['Species']
n2.input=['TModule']
L.add_Node(n2)
n3=TModule()
n3.rate=0.605289524207
n3.order=3
L.add_Node(n3)
n4=TFHill()
n4.hill=1.04843033105
n4.input=['TF']
n4.threshold=0.53797181977
n4.output=['TModule']
n4.order=4
L.add_Node(n4)
s2=Species()
s2.degradation=0.138432354831
s2.order=5
s2.types=['Species', 'Degradable', 'Phosphorylable', 'Complexable']
L.add_Node(s2)
n6=CorePromoter()
n6.delay=328
n6.output=['Species']
n6.input=['TModule']
n6.order=6
L.add_Node(n6)
n7=TModule()
n7.rate=0.963146138448
n7.order=7
L.add_Node(n7)
s3=Species()
s3.order=8
s3.types=['Species', 'Degradable', 'Phosphorylable', 'Complexable']
s3.degradation=0.347810723862
L.add_Node(s3)
n9=CorePromoter()
n9.order=9
n9.delay=394
n9.output=['Species']
n9.input=['TModule']
L.add_Node(n9)
n10=TModule()
n10.rate=0.911018328448
n10.order=10
L.add_Node(n10)
n11=TFHill()
n11.hill=1.43385606625
n11.input=['TF']
n11.threshold=0.703405821168
n11.output=['TModule']
n11.order=11
L.add_Node(n11)
s4=Species()
s4.activity=1
s4.degradation=0.843611422734
s4.order=12
s4.types=['Species', 'Degradable', 'Complex', 'Phosphorylable', 'TF']
L.add_Node(s4)
n13=PPI()
n13.disassociation=0.309140646394
n13.output=['Species']
n13.input=['Complexable', 'Complexable']
n13.order=13
n13.association=0.895477247988
L.add_Node(n13)
s5=Species()
s5.order=20
s5.types=['Species', 'Degradable', 'Complex', 'Phosphorylable']
s5.degradation=0.596028704786
L.add_Node(s5)
n15=PPI()
n15.disassociation=0.507818898738
n15.output=['Species']
n15.input=['Complexable', 'Complexable']
n15.order=21
n15.association=0.0966735959197
L.add_Node(n15)
L.graph.add_edge(s0,n4)
L.graph.add_edge(s0,n11)
L.graph.add_edge(s1,n13)
L.graph.add_edge(n2,s1)
L.graph.add_edge(n3,n2)
L.graph.add_edge(n4,n3)
L.graph.add_edge(s2,n13)
L.graph.add_edge(n6,s2)
L.graph.add_edge(n7,n6)
L.graph.add_edge(s3,n15)
L.graph.add_edge(s3,n15)
L.graph.add_edge(n9,s3)
L.graph.add_edge(n10,n9)
L.graph.add_edge(n11,n7)
L.graph.add_edge(n13,s4)
L.graph.add_edge(n15,s5)

L.remove_Node(s3)
L.clean_Nodes()
