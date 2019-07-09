import unittest
from src.gfa import GFA
from src.graph import Graph, Slice, Node, NoAnchorError, PathOverlapError, NoOverlapError, NodeMissingError

def G(rep):
    """Short hand for Graph construction that returns a slice"""
    if len(rep) > 1:
        raise ValueError("Warning: only the first slice will be returned.", rep)
    return Graph(rep)[0]

class GraphTest(unittest.TestCase):
    # Path 5 is sometimes introduced as a tie breaker for Slice.secondary()
    factory_input = [['ACGT', {1, 2, 3, 4}],
                     ['C', {1, 2, 4}, 'T', {3}],  # SNP
                     ['GGA', {1, 2, 3, 4}],  # anchor
                     ['C', {1, 2, 4}, '', {3}],  # [3] repeated from [1] SNP
                     ['AGTACG', {1, 2, 3}, 'CGTACT', {4}],  # [4] different membership from [3]
                     ['TTG', {1, 2, 3, 4}],  # [5] anchor
                     ['A', {1, 2}, 'C', {4, 5}, 'T', {3}],  # [6] third allele
                     ['GG', {1, 2}, 'TT', {3, 4}],  # [7] equal size nodes
                     ['C',{1, 2, 3, 5}, 'T',{4}], # [8] path slip
                     ['C',{1, 2, 5}, 'T',{3, 4}], # [9] path slip
                     ['C',{1, 2, 3}, 'T',{4}], # [10]path slip
                     ['TATA', {1, 2, 3, 4}]]  # [11] anchor
    def example_graph(self):
        factory_input = [Slice([Node('ACGT', {1,2,3,4})]),
                       Slice([Node('C',{1,2,4}),Node('T', {3})]),
                       Slice([Node('GGA',{1,2,3,4})]),
                       Slice([Node('C',{1,2,4}),Node('', {3})]),
                       Slice([Node('AGTACG',{1,2,3}), Node('CGTACT',{4})]),
                       Slice([Node('TTG',{1,2,3,4})]),
                       Slice([Node('A', {1, 2}), Node('C', {4, 5}), Node('T', {3})]),  # third allele
                       Slice([Node('GG', {1, 2}), Node('TT', {3, 4})]),  # equal size nodes
                       Slice([Node('C', {1, 2, 3, 5}), Node('T', {4})]),
                       Slice([Node('C', {1, 2, 5}), Node('T', {3, 4})]),
                       Slice([Node('C', {1, 2, 3}), Node('T', {4})]),
                       Slice([Node('TATA', {1, 2, 3, 4})])  # anchor
                          ]

        base_graph = Graph.load_from_slices(factory_input)
        return base_graph

    def test_graph_factory(self):
        base_graph = self.example_graph()
        assert base_graph == Graph(self.factory_input), \
            ('\n' + repr(base_graph) + '\n' + str(self.factory_input))
        g_double = Graph(eval(str(base_graph)))
        # WARN: Never compare two string literals: could be order sensitive, one object must be Graph
        #str(g_double) == str(base_graph)
        assert g_double == base_graph, repr(g_double) + '\n' + repr(base_graph)
        assert g_double == self.factory_input
        assert g_double == str(self.factory_input)

    def test_G(self):
        with self.assertRaises(ValueError):
            G([['C', {1, 2, 3, 4}], ['T', {12, 16}]])


class GFATest(unittest.TestCase):
    """ test class of gfa.py
    """

    def test_gfa(self):
        self.maxDiff = None
        location_of_xg = "../test/xg"
        graph = GFA.load_from_gfa("../test/test.gfa")
        graph.save_as_xg("../test/test.xg", location_of_xg)
        graph2 = GFA.load_form_xg("../test/test.xg", location_of_xg)
        self.assertFalse(self.is_different(graph.gfa, graph2.gfa))
#        self.assertEqual(len(graph.gfa.to_gfa1_s().split("\n")), len(graph2.gfa.to_gfa1_s().split("\n")))

    def test_load_gfa_to_graph(self):
        gfa = GFA.load_from_gfa("../test/test.gfa")
        graph = gfa.to_graph
        x = 'x'
        y = 'y'
        z = 'z'
        self.assertEqual(graph, [['CAAATAAG', {x, y, z}], ['A', {y, z}, 'G', {x}], ['C', {x, y, z}], ['TTG', {x, y, z}], ['A', {z}, 'G', {x, y}], ['AAATTTTCTGGAGTTCTAT', {x, y, z}], ['T', {x, y, z}], ['ATAT', {x, y, z}], ['T', {x, y, z}], ['CCAACTCTCTG', {x, y, z}]])

    def test_export_as_gfa(self):
        gfa = GFA.load_from_gfa("../test/test.gfa")
        graph = gfa.to_graph
        new_gfa = GFA.from_graph(graph)
        self.assertFalse(self.is_different(gfa.gfa, new_gfa.gfa))

    def test_load_gfa_to_graph_2(self):
        gfa = GFA.load_from_gfa("../test/test2.gfa")
        graph = gfa.to_graph
        self.assertIsNotNone(graph)

    def test_load_gfa_via_xg(self):
        location_of_xg = "../test/xg"
        graph = GFA.load_from_gfa("../test/test.gfa")
        graph.save_as_xg("../test/test.xg", location_of_xg)
        graph2 = GFA.load_form_xg("../test/test.xg", location_of_xg)
        graph = graph2.to_graph
        x = 'x'
        y = 'y'
        z = 'z'
        self.assertEqual(graph, [['CAAATAAG', {x, y, z}], ['A', {y, z}, 'G', {x}], ['C', {x, y, z}], ['TTG', {x, y, z}],
                                 ['A', {z}, 'G', {x, y}], ['AAATTTTCTGGAGTTCTAT', {x, y, z}], ['T', {x, y, z}],
                                 ['ATAT', {x, y, z}], ['T', {x, y, z}], ['CCAACTCTCTG', {x, y, z}]])

    @staticmethod
    def is_different(gfa1, gfa2):
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
        for p in gfa1.paths:
            p2 = [i for i in gfa2.paths if i.name == p.name]
            if p2 is None:
                different = True
            for i, l in enumerate(p.links):
                l2 = p2[0].links[i]
                if l.line.from_segment.name != l2.line.from_segment.name:
                    print(l, l2)
                    different = True
                if l.line.to_segment.name != l2.line.to_segment.name:
                    print(l, l2)
                    different = True
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

        
