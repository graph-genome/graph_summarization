import unittest
from src.gfa import GFA
from src.graph import Graph, Slice, Node, NoAnchorError, PathOverlapError, NoOverlapError, NodeMissingError, \
    Path


def G(rep):
    """Short hand for Graph construction that returns a slice"""
    if len(rep) > 1:
        raise ValueError("Warning: only the first slice will be returned.", rep)
    return Graph.build(rep)[0]


a, b, c, d, e = 'a', 'b', 'c', 'd', 'e'  # Paths must be created first
class GraphTest(unittest.TestCase):
    """Constructing a node with an existing Path object will modify that Path object (doubly linked)
    which means care must be taken when constructing Graphs.  From factory_input we have an example of
    pure python to Graph.build in one step.  In example_graph, we must first declare the Paths,
    then reference them in order in Node Constructors.  Order matters for Graph identity!"""
    # Path e is sometimes introduced as a tie breaker for Slice.secondary()
    factory_input = [['ACGT', {a, b, c, d}],
                     ['C', {a, b, d}, 'T', {c}],  # SNP
                     ['GGA', {a, b, c, d}],  # anchor
                     ['C', {a, b, d}, '', {c}],  # [3] repeated from [1] SNP
                     ['AGTACG', {a, b, c}, 'CGTACT', {d}],  # [4] different membership from [3]
                     ['TTG', {a, b, c, d}],  # [5] anchor
                     ['A', {a, b}, 'C', {d, e}, 'T', {c}],  # [6] third allele
                     ['GG', {a, b}, 'TT', {c, d}],  # [7] equal size nodes
                     ['C', {a, b, c, e}, 'T', {d}],  # [8] path slip
                     ['C', {a, b, e}, 'T', {c, d}],  # [9] path slip
                     ['C', {a, b, c}, 'T', {d}],  # [10]path slip
                     ['TATA', {a, b, c, d}]]  # [11] anchor

    def example_graph(self):
        # IMPORTANT: Never reuse Paths: Paths must be created fresh for each graph
        a, b, c, d, e = Path('a'), Path('b'), Path('c'), Path('d'), Path('e')
        factory_input = [Slice([Node('ACGT', {a,b,c,d})]),
                       Slice([Node('C',{a,b,d}),Node('T', {c})]),
                       Slice([Node('GGA',{a,b,c,d})]),
                       Slice([Node('C',{a,b,d}),Node('', {c})]),
                       Slice([Node('AGTACG',{a,b,c}), Node('CGTACT',{d})]),
                       Slice([Node('TTG',{a,b,c,d})]),
                       Slice([Node('A', {a, b}), Node('C', {d, e}), Node('T', {c})]),  # third allele
                       Slice([Node('GG', {a, b}), Node('TT', {c, d})]),  # equal size nodes
                       Slice([Node('C', {a, b, c, e}), Node('T', {d})]),
                       Slice([Node('C', {a, b, e}), Node('T', {c, d})]),
                       Slice([Node('C', {a, b, c}), Node('T', {d})]),
                       Slice([Node('TATA', {a, b, c, d})])  # anchor
                          ]

        base_graph = Graph.load_from_slices(factory_input)
        return base_graph

    def test_equalities(self):
        self.assertEqual(Node('A', {}),Node('A', {}))
        self.assertEqual(Node('A', {Path('x')}),Node('A', {Path('x')}))
        self.assertEqual(Node('A', {Path('x'),Path('y')}),Node('A', {Path('x'),Path('y')}))
        self.assertEqual(Slice([Node('ACGT', {Path('a'), Path('b'), Path('c'), Path('d')})]),
                         Slice([Node('ACGT', {Path('a'), Path('b'), Path('c'), Path('d')})]))
        self.assertEqual(Graph.build([['ACGT', {a, b, c, d}]]), Graph.build([['ACGT', {a, b, c, d}]]))

    def test_graph_factory(self):
        base_graph = self.example_graph()
        g1, g2 = Graph.build(self.factory_input), Graph.build(self.factory_input)
        assert g1 == g2, \
            ('\n' + repr(g1) + '\n' + repr(g2))
        g_double = Graph.build(eval(str(base_graph)))
        # WARN: Never compare two string literals: could be order sensitive, one object must be Graph
        #str(g_double) == str(base_graph)
        assert g_double == base_graph, repr(g_double) + '\n' + repr(base_graph)
        assert g_double == self.factory_input
        assert g_double == str(self.factory_input)

    def test_G(self):
        with self.assertRaises(ValueError):
            G([['C', {Path('a'), Path('b')}], ['T', {Path('12'), Path('16')}]])


class GFATest(unittest.TestCase):
    """ test class of gfa.py
    """

    @unittest.expectedFailure
    def test_gfa(self):
        self.maxDiff = None
        location_of_xg = "../test/xg"
        graph = GFA.load_from_gfa("../test/test.gfa")
        graph.save_as_xg("../test/test.xg", location_of_xg)
        graph2 = GFA.load_from_xg("../test/test.xg", location_of_xg)
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

    @unittest.expectedFailure
    def test_load_gfa_via_xg(self):
        location_of_xg = "../test/xg"
        graph = GFA.load_from_gfa("../test/test.gfa")
        graph.save_as_xg("../test/test.xg", location_of_xg)
        graph2 = GFA.load_from_xg("../test/test.xg", location_of_xg)
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
            p2 = [i for i in gfa2.paths if i.accession == p.accession]
            if p2 is None:
                different = True
            for i, l in enumerate(p.links):
                l2 = p2[0].links[i]
                if l.line.from_segment.accession != l2.line.from_segment.accession:
                    print(l, l2)
                    different = True
                if l.line.to_segment.accession != l2.line.to_segment.accession:
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

        
