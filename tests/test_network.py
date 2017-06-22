"""
Test module for the network class
"""
import unittest
import phievo

class TestNetwork(unittest.TestCase):
    def setUp(self):
        self.net = phievo.Networks.classes_eds2.Network()

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

if __name__ == '__main__':
    unittest.main()
