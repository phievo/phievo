"""
Test module for the network class
"""
import unittest
import phievo

class TestNetwork(unittest.TestCase):
    def setUp(self):
        self.net = phievo.Networks.classes_eds2.Network()

    def test_print(self):
        self.assertEqual(self.net.__str__(),"")

if __name__ == '__main__':
    unittest.main()
