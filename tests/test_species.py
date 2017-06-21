"""
Test module for the Species class
phievo.Networks.classes_eds2.Species
"""
import unittest
import phievo

class TestNode(unittest.TestCase):
    def setUp(self):
        self.degradable = phievo.Networks.classes_eds2.Species([['Degradable',0.1]])
        self.input = phievo.Networks.classes_eds2.Species([['Input',0]])

    def test_isinstance(self):
        self.assertTrue(self.degradable.isinstance('Species'))
        self.assertTrue(self.degradable.isinstance('Degradable'))
        self.assertFalse(self.degradable.isinstance('Kinase'))
        self.assertTrue(self.input.isinstance('Input'))

    def test_list_types(self):
        self.assertIn('Species',self.degradable.list_types())
        self.assertIn('Degradable',self.degradable.list_types())
        self.assertIn('Input',self.input.list_types())

    def test_def_label(self):
        self.assertEqual(self.degradable.def_label(),True)
        self.assertEqual(self.input.def_label(),True)
        self.assertEqual(self.degradable.label,", Species, Degradable, 0.1")
        self.assertEqual(self.input.label,", Species, Input, 0")

if __name__ == '__main__':
    unittest.main()
