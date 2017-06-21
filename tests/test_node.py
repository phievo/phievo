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
    def isinstance(self,name):
        return name in self.types

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

if __name__ == '__main__':
    unittest.main()
