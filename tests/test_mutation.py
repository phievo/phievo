"""
Test module for the mutation module
"""
import unittest
import phievo
from math import log
from collections import Counter #to test unordered lists
mut = phievo.Networks.mutation #shortcut for the test module

class mock_interaction(phievo.Networks.classes_eds2.Interaction):
    def __init__(self,list_input,list_output,net):
        super().__init__()
        net.add_Node(self)
        self.removable  = True
        self.input = len(list_input)*['Species']
        self.output = len(list_output)*['Species']
        for input in list_input:
            net.graph.add_edge(input,self)
        for output in list_output:
            net.graph.add_edge(self,output)

class mock_rg(object):
    """A false predictable random number generator for test evaluation"""
    def __init__(self,v):
        self.value = v
    def random(self):
        return self.value

class TestFunctions(unittest.TestCase):
    def setUp(self):
        self.net = phievo.Networks.mutation.Mutable_Network()
        self.s1 = self.net.new_Species([['Input',0]])
        self.s2 = self.net.new_Species([['Input',1]])
        self.s3 = self.net.new_Species([['Output',0]])
        self.inter1 = mock_interaction([self.s1,self.s2],[self.s3],self.net)
        self.inter2 = mock_interaction([self.s2],[self.s3],self.net)
    
    def test_build_lists(self):
        dict_mutation = {"mutate_Node('Species')":.1,
                         "mutate_Node('PPI')":.1,
                         "remove_Interaction('TFHill')":.2,
                         "remove_Interaction('PPI')":.2,
                         "random_Interaction('Species')":.05,
                         "random_Interaction('TFHill')":.05}
        l_mutate, l_remove, l_create = mut.build_lists(dict_mutation)
        self.assertTrue(Counter(l_mutate) == Counter(['Species','PPI']))
        self.assertTrue(Counter(l_remove) == Counter(['TFHill','PPI']))
        self.assertTrue(Counter(l_create) == Counter(['TFHill','Species']))

    def test_sample_dictionary_ranges(self):
        SDR = mut.sample_dictionary_ranges
        mut.dictionary_ranges['test_float'] = 1.0
        mut.dictionary_ranges['test_interval'] = (.5,2.5)
        mut.dictionary_ranges['CorePromoter.delay'] = 10
        self.assertEqual(SDR('test_float',mock_rg(.2)),.2)
        self.assertEqual(SDR('test_interval',mock_rg(.3)),.5+.3*2)
        self.assertEqual(SDR('CorePromoter.delay',mock_rg(.367)),3)
        with self.assertRaises(KeyError):
            SDR('not_a_key',mock_rg(0.))

    def test_random_parameters(self):
        self.assertEqual(mut.random_parameters('Species',mock_rg(1.)),[['Degradable', 1.0], ['Phosphorylable'], ['Phospho', 0]])
        self.assertEqual(mut.random_parameters('TF',mock_rg(.8)),[['Degradable', .8], ['Phosphorylable'], ['Phospho', 0],['Diffusible',.8],['TF',1],['Complexable']])
        with self.assertRaises(TypeError):
            mut.random_parameters(['Shruberry'],mock_rg(1.))
    
    def test_rand_modify(self):
        self.inter1.mutable = True
        self.inter1.value = 0.0
        mut.dictionary_ranges['mock_interaction.value'] = 1.
        mut.rand_modify(self.inter1,mock_rg(.2))
        self.assertEqual(self.inter1.value,0.2)
        
        mut.dictionary_ranges['relative_variation'] = 1.
        mut.rand_modify(self.inter1,mock_rg(.2))
        self.assertEqual(self.inter1.value,0.2)
        
        self.inter1.mutable = False
        mut.rand_modify(self.inter1,mock_rg(.3))
        self.assertEqual(self.inter1.value,0.2)

class TestMutableNetwork(unittest.TestCase):
    def setUp(self):
        self.net = phievo.Networks.mutation.Mutable_Network()
        self.s1 = self.net.new_Species([['Input',0]])
        self.s2 = self.net.new_Species([['Input',1]])
        self.s3 = self.net.new_Species([['Output',0]])
        self.inter1 = mock_interaction([self.s1,self.s2],[self.s3],self.net)
        self.inter2 = mock_interaction([self.s2],[self.s3],self.net)
        self.net.write_id()

    def test_random_Species(self):
        s4 = self.net.random_Species('Species')
        self.assertTrue(s4 in self.net.nodes())

    def test_random_Interaction(self):
        # delegate to suited interaction, nothing to test
        pass
    
    def test_remove_Interaction(self):
        self.net.write_id()
        self.assertTrue(self.net.remove_Interaction('mock_interaction'))
        self.assertEqual(len(self.net.dict_types['Interaction']),1)
        self.assertFalse(
        self.net.remove_Interaction('Shruberry'))
        with self.assertRaises(TypeError):
            self.net.remove_Interaction('Species')

    def test_random_remove_output(self):
        self.assertTrue(self.net.random_remove_output())
        self.assertFalse('Output' in self.s3.types)
        self.assertFalse(self.net.random_remove_output())
    
    def test_random_add_output(self):
        self.s4 = self.net.new_Species([['Kinase']])
        self.assertTrue(self.net.random_add_output())
        self.assertTrue('Output' in self.s4.types)
        self.assertFalse(self.net.random_add_output())

    def test_random_duplicate(self):
        pass
    
    def test_mutate_Node(self):
        # Rely on already tested function
        pass
    
    def test_build_mutations(self):
        dict_mutation = {"mutate_Node('Species')":.1,
                         "mutate_Node('PPI')":.1,
                         "remove_Interaction('TFHill')":.2,
                         "remove_Interaction('PPI')":.2,
                         "random_Interaction('Species')":.05,
                         "random_Interaction('TFHill')":.05,
                         "random_remove_output()":1.0,
                         "random_change_output()":1.0}
        self.net.number_Species = lambda :len(self.net.dict_types['Species'])
        self.net.number_TFHill = lambda :0
        mut.dictionary_mutation = dict_mutation
        mut.build_lists(dict_mutation)
        
        new_dict = {"mutate_Node('Species')":.1*3,
                    "random_Interaction('Species')":.05*3,
                    "random_remove_output()":1.0,
                    "random_change_output()":1.0}
        self.assertEqual(self.net.build_mutations(),new_dict)

    def test_compute_next_mutation(self):
        # Modify the set up to run the test
        dict_mutation = {"mutate_Node('Species')":.1,
                         "mutate_Node('PPI')":.1,
                         "remove_Interaction('TFHill')":.2,
                         "remove_Interaction('PPI')":.2,
                         "random_Interaction('Species')":.05,
                         "random_Interaction('TFHill')":.05,
                         "random_remove_output()":1.0,
                         "random_change_output()":1.0}
        mut.dictionary_mutation = dict_mutation
        self.net.number_Species = lambda :len(self.net.dict_types['Species'])
        self.net.number_TFHill = lambda :0
        
        self.net.Random = mock_rg(.5)
        res = self.net.compute_next_mutation()
        self.assertIs(type(res),list)
        t,com = res
        self.assertAlmostEqual(t,-log(.5)/2.45)
        self.assertEqual(com,"random_change_output()")
        
        self.net.Random = mock_rg(.1)
        t,com = self.net.compute_next_mutation()
        self.assertAlmostEqual(t,-log(.1)/2.45)
        self.assertEqual(com,"mutate_Node('Species')")
        
if __name__ == '__main__':
    unittest.main()