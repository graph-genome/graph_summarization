import unittest
from gfa import GFA

class GFATest(unittest.TestCase):
    """ test class of gfa.py
    """

    def test_gfa(self):
        ### Usage
        self.maxDiff = None
        location_of_xg = "test/xg"
        graph = GFA.load_from_gfa("test/test.gfa")
        graph.save_as_xg("test/test.xg", location_of_xg)
        graph2 = GFA.load_form_xg("test/test.xg", location_of_xg)
        self.assertEqual(len(graph.gfa.to_gfa1_s().split("\n")), len(graph2.gfa.to_gfa1_s().split("\n")))
 
if __name__ == "__main__":
    unittest.main()        

        
