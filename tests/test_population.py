import unittest,os,shutil
import phievo
from phievo.Populations_Types.evolution_gillespie import Population
from phievo.initialization_code import init_networks,init_evolution,check_model_dir,make_workplace_dir
class TestPopulation(unittest.TestCase):
    def setUp(self):
        self.test_directory = "tests"+os.sep+"lac_operon"
        shutil.copytree("Examples"+os.sep+"lac_operon",self.test_directory)
    def test_identifier(self):
        [model_dir, inits, init_file] = check_model_dir(self.test_directory)
        inits.prmt.update(dict(npopulation=10,nseed=1,nstep=2,multipro_level=0))
        deriv2 = init_networks(inits)
        [mutation, evolution_gillespie] = init_evolution(inits, deriv2)
        seed_dir = self.test_directory+os.sep+"Seed0"
        os.mkdir(seed_dir)
        pop = Population(seed_dir)
        pop.initialize_identifier()
        #import pdb;pdb.set_trace()
        for i,net  in enumerate(pop.genus):
            self.assertEqual(net.identifier,i)

        inits.prmt["workplace_dir"] = make_workplace_dir(seed_dir)
        #import pdb;pdb.set_trace()

        pop.evolution(inits.prmt)

    def tearDown(self):
        shutil.rmtree(self.test_directory)
