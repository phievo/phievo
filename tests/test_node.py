"""
Test module for the Node class
phievo.Networks.classes_eds2.Node
"""
import unittest
import phievo

class mock_node():
    '''A pseudo node to test check_consistency'''
    def __init__(self,type_list):
        self.types = type_list
        self.removable = True
    def isinstance(self,name):
        return name in self.types

class mock_interaction(phievo.Networks.classes_eds2.Interaction):
    def __init__(self,list_input,list_output,net):
        net.add_Node(self)
        self.removable  = True
        for input in list_input:
            net.graph.add_edge(input,self)
        for output in list_output:
            net.graph.add_edge(self,output)
    def outputs_to_delete(self,net):
        return net.graph.successors(self)

class TestNode(unittest.TestCase):
    def setUp(self):
        self.node = phievo.Networks.classes_eds2.Node()

    def test_id(self):
        self.assertEqual(self.node.id,'None')
        self.assertEqual(self.node.int_id(),None)

    def test_isinstance(self):
        self.assertTrue(self.node.isinstance('Node'))
        self.assertFalse(self.node.isinstance('Species'))

    def test_list_types(self):
        self.assertEqual(self.node.list_types(),['Node'])

    def test_check_consistency(self):
        A = mock_node(['a'])
        B = mock_node(['a','b'])
        C = mock_node(['c','d'])
        cc = phievo.Networks.classes_eds2.check_consistency
        
        self.assertTrue(cc(['a'],[A]))
        self.assertTrue(cc(['b'],[B]))
        self.assertFalse(cc(['a'],[C,A]))
        self.assertTrue(cc(['a','b'],[B,A]))
        self.assertFalse(cc(['a','b'],[C,A]))
        self.assertTrue(cc(['a','b','d'],[C,B,A]))
        self.assertTrue(cc(['a','b','d'],[C,A,B]))
        self.assertFalse(cc(['a','b','a'],[C,A,B]))

    def test_isremovable(self):
        # Dummy network creation
        self.net = phievo.Networks.classes_eds2.Network()
        self.s1 = self.net.new_Species([['Input',0]])
        self.s2 = self.net.new_Species([['Input',1]])
        self.s3 = self.net.new_Species([['Output',0]])
        self.s4 = self.net.new_Species([['Kinase']])
        self.inter1 = mock_interaction([self.s1,self.s2],[self.s3],self.net)
        self.inter2 = mock_interaction([self.s2],[self.s4],self.net)
        # For now, everything should be removable
        self.assertTrue(self.inter1.isremovable(self.net))
        self.assertTrue(self.inter2.isremovable(self.net))
        # Now, everything should be unremovable
        self.s3.removable = False
        phievo.Networks.classes_eds2.list_unremovable.append('Kinase')
        self.assertFalse(self.inter1.isremovable(self.net))
        self.assertFalse(self.inter2.isremovable(self.net))
        
if __name__ == '__main__':
    unittest.main()
