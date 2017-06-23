"""
Test module for the network class
"""
import unittest
import phievo

class mock_interaction(phievo.Networks.classes_eds2.Interaction):
    def __init__(self,list_input,list_output,net):
        net.add_Node(self)
        for input in list_input:
            net.graph.add_edge(input,self)
        for output in list_output:
            net.graph.add_edge(self,output)

class TestNetwork(unittest.TestCase):
    def setUp(self):
        self.net = phievo.Networks.classes_eds2.Network()
        self.s1 = self.net.new_Species([['Input',0]])
        self.s2 = self.net.new_Species([['Input',1]])
        self.s3 = self.net.new_Species([['Output',0]])
        self.inter1 = mock_interaction([self.s1,self.s2],[self.s3],self.net)
        self.inter2 = mock_interaction([self.s2],[self.s3],self.net)

    def test_add_Node(self):
        self.species1 = phievo.Networks.classes_eds2.Species([['Degradable',0.1]])
        self.assertTrue(self.net.add_Node(self.species1))
        self.assertFalse(self.net.add_Node(self.species1))
        self.assertIn(self.species1,self.net.nodes())

    def test_new_Species(self):
        self.spc = self.net.new_Species([['Degradable',0.1]])
        self.assertIn(self.spc,self.net.nodes())

    def test_number_nodes(self):
        self.net.dict_types = dict(a=[1,1],b=[]) #dummy list_types
        self.assertEqual(self.net.number_nodes('a'),2)
        self.assertEqual(self.net.number_nodes('b'),0)
        self.assertEqual(self.net.number_nodes('c'),0)

    def test_catal_data(self):
        cd = self.net.catal_data
        self.inter3 = mock_interaction([self.s1,self.s2],[self.s3,self.s2],self.net)
        self.inter4 = mock_interaction([self.s1],[self.s1,self.s2],self.net)
        A,B,C = cd(self.inter3)
        self.assertEqual(A,[self.s2])
        self.assertEqual(B,[self.s1])
        self.assertEqual(C,[self.s3])
        A,B,C = cd(self.inter4)
        self.assertEqual(A,[self.s1])
        self.assertEqual(B,[self.s1])
        self.assertEqual(C,[self.s2])
        A,B,C = cd(self.inter1)
        self.assertEqual(A,[])
        self.assertEqual(set(B),{self.s1,self.s2})
        self.assertEqual(C,[self.s3])

    def test_check_existing_binary(self):
        ceb = self.net.check_existing_binary
        self.assertTrue(ceb([self.s1,self.s2],'Interaction'))
        self.assertTrue(ceb([self.s1,self.s2],'mock_interaction'))
        self.assertFalse(ceb([self.s1,self.s3],'mock_interaction'))
        self.assertFalse(ceb([self.s1],'mock_interaction'))
        self.assertTrue(ceb([self.s2],'Interaction'))

    def test_check_existing_link(self):
        cel = self.net.check_existing_link
        self.assertTrue(cel([self.s1,self.s2,self.s3],'mock_interaction'))
        self.assertTrue(cel([self.s2,self.s3],'mock_interaction'))
        self.assertFalse(cel([self.s1,self.s2],'mock_interaction'))
        self.assertFalse(cel([self.s1,self.s3],'mock_interaction'))
        self.assertFalse(cel([self.s1,self.s3],'King_of_the_britton'))
        

if __name__ == '__main__':
    unittest.main()
