import unittest
from gfa import GFA

class GFATest(unittest.TestCase):
    """ test class of gfa.py
    """

    def test_gfa(self):
        ### Usage
        location_of_xg = "test/xg"
        graph = GFA.load_from_gfa("test/test.gfa")
        graph.save_as_xg("test/test.xg", location_of_xg)
        graph2 = GFA.load_form_xg("test/test.", location_of_xg)
        graph2.save_as_pickle("test/test.pickle")
        graph3 = GFA.load_from_pickle("test/test.pickle")
        self.assertEqual(graph.gfa.to_gfa_1_s, graph3.to_gfa_1_s)


if __name__ == "__main__":
    unittest.main()        

        
