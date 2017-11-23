from phievo.AnalysisTools.main_functions import download_tools,download_example
import unittest
import shutil,os,glob

class TestDownloadFunctions(unittest.TestCase):
    def test_download_tools(self):
        AnalyseRun="test_AnalyseRun.py"
        run_evolution = "test_run_evolution.py"
        download_tools(run_evolution=run_evolution,AnalyseRun=AnalyseRun)
        self.assertTrue(os.path.isfile(run_evolution))
        self.assertTrue(os.path.isfile(AnalyseRun))
        os.remove(run_evolution)
        os.remove(AnalyseRun)
    def test_download_example(self):
        directory = "test_directory_adaptation"
        download_example("adaptation",directory=directory)
        for ff in glob.glob(os.path.join("Examples","adaptation","*")):
            ff=ff.split(os.sep)[-1]
            self.assertTrue(os.path.isfile(os.path.join(directory,ff)))
        shutil.rmtree(directory)
        
