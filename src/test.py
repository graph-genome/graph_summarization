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
        self.assertFalse(self.is_different(graph.gfa, graph2.gfa))
#        self.assertEqual(len(graph.gfa.to_gfa1_s().split("\n")), len(graph2.gfa.to_gfa1_s().split("\n")))

    def is_different(self, gfa1, gfa2):
        different = False
        for s in gfa1.segments:
            s2 = gfa2.segment(s)
            if s2 is None:
                different = True
            if s.diff(s2):
                different = True
                for diff in s.diff(s2):
                    print(diff)
        for s in gfa2.segments:
            s1 = gfa1.segment(s)
            if s1 is None:
                different = True

"""
        for s in gfa1.edges:
            s2 = gfa2.edges[gfa2.edges.index(s)]
            if s2 is None:
                different = True
            if s.diff(s2):
                different = True
                for diff in s.diff(s2):
                    print(diff)
        for s in gfa2.edges:
            print(s)
#            s1 = gfa1.edges[]
            if gfa1.edges.index(s) is None:
                different = True
"""        
        return different

if __name__ == "__main__":
    unittest.main()        

        
