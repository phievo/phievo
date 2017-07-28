"""
Test module for the deriv2 module
"""
import unittest
import phievo
d2 = phievo.Networks.deriv2 #shortcut for the tested module

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

class Test_Routine_Functions(unittest.TestCase):
    def setUp():
        self.net = phievo.Networks.mutation.Mutable_Network()
        self.s1 = self.net.new_Species([['Input',0]])
        self.s2 = self.net.new_Species([['Input',1]])
        self.s3 = self.net.new_Species([['Output',0]])
        self.inter1 = mock_interaction([self.s1,self.s2],[self.s3],self.net)
        self.inter2 = mock_interaction([self.s2],[self.s3],self.net)
        
    def test_compute_leap(self):
        res = d2.compute_leap(['S1'],['S2','S3'],'RATE')
        res = [line.strip() for line in res.split()]
        hope = ['rate=RATE;', 'increment=rate;', 'dS1-=increment;', 'dS2+=increment;', 'dS3+=increment;']
        self.assertEqual(res,hope)
        
        res = d2.compute_leap([],['S1'],'RATE')
        res = [line.strip() for line in res.split()]
        hope = ['rate=RATE;', 'increment=rate;', 'dS1+=increment;']
        self.assertEqual(res,hope)
        
        d2.noise_flag = True
        res = d2.compute_leap(['S1'],[],'RATE')
        res = [line.strip() for line in res.split()]
        hope = ['rate=RATE;', 'increment=compute_noisy_increment(rate);', 'dS1-=increment;']
        self.assertEqual(res,hope)
        d2.noise_flag = False

if __name__ == '__main__':
    unittest.main()