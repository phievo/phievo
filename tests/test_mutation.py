"""
Test mutation
"""
import unittest
import phievo
import phievo.Networks.mutation as mutation
import random
class TestMutation:
    def setUp(self):
        seed=0
        g=random.Random(seed)
        self.net = mutation.Mutable_Network(g)
    def test_duplicate(self):
        [tm0, prom0, out0] = self.net.random_gene('TF')
        out0.add_type(['Output',0])
        [tm1, prom1, out1] = self.net.random_gene('TF')
        out1.add_type(['Output',1])
        [tm_in0, prom_in0, in1] = self.net.random_gene('TF')
        in1.add_type(['Input',0])

        [D_tm_in0,D_prom_in0,D_in0] = self.net.duplicate_gene(in0)
        
        [D_tm_out0,D_prom_out0,D_out0] = self.net.duplicate_gene(out0)
        self.net.remove_output_when_duplicate = False
