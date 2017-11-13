"""
Test module for the multiplicative logic
phievo.Networks.classes_eds2.Node
"""
import unittest
import phievo

class TestMultipl(unittest.TestCase):
    def setUp(self):
        self.net = phievo.Networks.mutation.Mutable_Network()
        self.net.fixed_activity_for_TF = False

        self.sp0 = self.net.new_Species([["Input",0],["TF",1]])
        self.sp1 = self.net.new_Species([["Input",1],["TF",1]])
        self.TMod, self.CoreP, self.Output = self.net.new_gene(0.75, 0, [["Output",0],["Degradable",.25]], basal=0.125)
        self.net.new_TFHill(self.sp0, 1, .1, self.TMod, activity=1)
        self.net.new_TFHill(self.sp1, 2, .2, self.TMod, activity=1)
        self.net.write_id()

    def test_show(self):
        print(self.net)

    def test_multipl(self):
        print(phievo.Networks.deriv2.compute_transcription(self.net,self.TMod))

if __name__ == '__main__':
    unittest.main()
