"""
Test module for the Species class
phievo.Networks.classes_eds2.Species
"""
import unittest
import phievo

class TestSpecies(unittest.TestCase):
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

    def test_change_type(self):
        self.degradable.change_type('Degradable',[.2])
        self.assertEqual(self.degradable.degradation,.2)
        self.input.change_type('Input',[2])
        self.assertEqual(self.input.n_put,2)
        with self.assertRaises(ValueError):
            self.degradable.change_type('Degradable',[.2,2])

    def test_add_type(self):
        with self.assertRaises(TypeError):
            self.degradable.add_type('Kinase')
        with self.assertRaises(Warning):
            self.degradable.add_type(['Degradable',0.1])
        with self.assertRaises(ValueError):
            self.degradable.add_type(['Shrubbery'])
        self.assertEqual(self.input.add_type(['Degradable',0.15]),True)
        self.assertEqual(self.input.degradation,.15)

    def test_clean_type(self):
        self.degradable.clean_type('Degradable')
        self.assertFalse(self.degradable.isinstance('Degradable'))
        

if __name__ == '__main__':
    unittest.main()
