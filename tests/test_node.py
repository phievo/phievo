"""
Test module for the Node class
phievo.Networks.classes_eds2.Node
"""
import unittest
import phievo

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

if __name__ == '__main__':
    unittest.main()
