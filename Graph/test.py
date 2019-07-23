import unittest
import os
from os.path import join
from Graph.gfa import GFA
from Graph.models import Graph, Slice, Node, Path, SlicedGraph
from Graph.sort import DAGify

# Define the working directory
WD = os.path.dirname(__file__)
# as our current setup stores the test data in an extra folder this is a dirty workaround
# hopefully Travis eats this
WD = WD[0:-4]
# Define several test example directories
PATH_TO_TEST_DATA = join(WD, "test_data/")
location_of_xg = "../test_data/xg"


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
        paths = [a, b, c, d, e]
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

        base_graph = SlicedGraph.load_from_slices(factory_input, paths)
        return base_graph

    def test_equalities(self):
        self.assertEqual(Node('A', {}),Node('A', {}))
        self.assertEqual(Node('A', {Path('x')}),Node('A', {Path('x')}))
        self.assertEqual(Node('A', {Path('x'),Path('y')}),Node('A', {Path('x'),Path('y')}))
        self.assertEqual(Slice([Node('ACGT', {Path('a'), Path('b'), Path('c'), Path('d')})]),
                         Slice([Node('ACGT', {Path('a'), Path('b'), Path('c'), Path('d')})]))
        self.assertEqual(SlicedGraph.build([['ACGT', {a, b, c, d}]]), SlicedGraph.build([['ACGT', {a, b, c, d}]]))

    def test_graph_factory(self):
        base_graph = self.example_graph()
        g1, g2 = SlicedGraph.build(self.factory_input), SlicedGraph.build(self.factory_input)
        assert g1 == g2, \
            ('\n' + repr(g1) + '\n' + repr(g2))
        g_double = SlicedGraph.build(eval(str(base_graph)))
        # WARN: Never compare two string literals: could be order sensitive, one object must be Graph
        #str(g_double) == str(base_graph)
        assert g_double == base_graph, repr(g_double) + '\n' + repr(base_graph)
        assert g1 == base_graph, repr(g1) + '\n' + repr(base_graph)
        assert g_double == self.factory_input
        assert g_double == str(self.factory_input)

    def test_G(self):
        with self.assertRaises(ValueError):
            G([['C', {Path('a'), Path('b')}], ['T', {Path('12'), Path('16')}]])

x,y,z,a = 'x', 'y', 'z', 'a'

class DAGifyTest(unittest.TestCase):
    """ test class of sort.py
    """


    def test_dagify(self):
        gfa = GFA.load_from_gfa(PATH_TO_TEST_DATA + "/test.gfa")
        paths = gfa.to_paths
        dagify = DAGify(paths)
        profile = dagify.generate_profiles(0)
        graph = SlicedGraph.load_from_slices(dagify.to_slices(profile), paths)
#        x, y, z = graph.paths['x'], graph.paths['y'], graph.paths['z']

        self.assertEqual([['CAAATAAG', {x,y,z}], ['A', {y,z}, 'G', {x}], ['C', {x,y,z}], ['TTG', {x,y,z}], ['A', {z}, 'G', {x,y}], ['AAATTTTCTGGAGTTCTAT', {x,y,z}], ['T', {x,y,z}], ['ATAT', {x,y,z}], ['T', {x,y,z}], ['CCAACTCTCTG', {x,y,z}]], graph)

    def test_dagify2(self):
        gfa = GFA.load_from_gfa(PATH_TO_TEST_DATA + "/test2.gfa")
        paths = gfa.to_paths
        dagify = DAGify(paths)
        profile = dagify.generate_profiles(0)
        graph = SlicedGraph.load_from_slices(dagify.to_slices(profile), paths)
        x,y,z,a = 'x', 'y', 'z', 'a'
        self.assertEqual([['CAAATAAG', {x, y, z}], ['G', {x}, 'A', {y, z}], ['C', {x, y}, 'T', {z}], ['TTG', {x, y, z}], ['G', {x, y}, 'A', {a, z}], ['AAATTTTCTGGAGTTCTAT', {a, x, y, z}], ['A', {a, z}, 'T', {x, y}], ['ATAT', {x, y, z}], ['T', {x, y, z}], ['CCAACTCTCTG', {x, y, z}]], graph)

    def test_dagify3(self):
        gfa = GFA.load_from_gfa(PATH_TO_TEST_DATA + "/test3.gfa")
        paths = gfa.to_paths
        dagify = DAGify(paths)
        profile, rep_count = dagify.generate_profiles_with_minimizing_replications()
        graph = SlicedGraph.load_from_slices(dagify.to_slices(profile), paths)
        self.assertEqual(rep_count, 1)
        self.assertEqual(graph, [['CAAATAAG', {x, y}], ['CCAACTCTCTG', {y}, 'G', {x}], ['C', {x, y}], ['TTG', {x, y}], ['G', {x, y}], ['AAATTTTCTGGAGTTCTAT', {x, y}], ['T', {x, y}], ['ATAT', {x, y}], ['T', {x, y}], ['CCAACTCTCTG', {x, y}]])

    def test_dagify_altpath(self):
        gfa = GFA.load_from_gfa(PATH_TO_TEST_DATA + "/alternate_paths.gfa")
        paths = gfa.to_paths
        dagify = DAGify(paths)
        profile, rep_count = dagify.generate_profiles_with_minimizing_replications()
        graph = SlicedGraph.load_from_slices(dagify.to_slices(profile), paths)
        self.assertEqual(rep_count, 1)
        self.assertEqual(graph, [['CAAATAAG', {x, y}], ['A', {x}, '', {y}], ['G', {x, y}], ['A', {y}, '', {x}], ['T', {x, y}]])

    def test_dagify_dup(self):
        gfa = GFA.load_from_gfa(PATH_TO_TEST_DATA + "/duplicate.gfa")
        paths = gfa.to_paths
        dagify = DAGify(paths)
        profile, rep_count = dagify.generate_profiles_with_minimizing_replications()
        graph = SlicedGraph.load_from_slices(dagify.to_slices(profile), paths)
        self.assertEqual(rep_count, 2)
        self.assertEqual(graph, [['CAAATAAG', {x, y}], ['', {x}, 'A', {y}], ['', {x}, 'G', {y}], ['A', {x, y}], ['G', {x, y}], ['T', {x, y}]])


    def test_unresolved_repreat(self):
        gfa = GFA.load_from_gfa(PATH_TO_TEST_DATA + "/unresolved_repeat.gfa")
        paths = gfa.to_paths
        dagify = DAGify(paths)
        profile, rep_count = dagify.generate_profiles_with_minimizing_replications()
        graph = SlicedGraph.load_from_slices(dagify.to_slices(profile), paths)
        self.assertEqual([['CAAATAAG', {'x'}, 'T', {'y'}], ['A', {'y', 'x'}], ['G', {'x'}, 'C', {'y'}]], graph)

    @unittest.skip("Inversion is unsupported")
    def test_inversion(self):
        gfa = GFA.load_from_gfa(PATH_TO_TEST_DATA + "/inversion.gfa")
        paths = gfa.to_paths
        dagify = DAGify(paths)
        profile, rep_count = dagify.generate_profiles_with_minimizing_replications()
        graph = SlicedGraph.load_from_slices(dagify.to_slices(profile), paths)
        self.assertEqual(graph, [])

    @unittest.skip("Inversion is unsupported")
    def test_nested_inversion(self):
        gfa = GFA.load_from_gfa(PATH_TO_TEST_DATA + "/nested_inv.gfa")
        paths = gfa.to_paths
        dagify = DAGify(paths)
        profile, rep_count = dagify.generate_profiles_with_minimizing_replications()
        graph = SlicedGraph.load_from_slices(dagify.to_slices(profile), paths)
        self.assertEqual(graph, [])

    @unittest.skip("Inversion is unsupported")
    def test_simple_inversion(self):
        gfa = GFA.load_from_gfa(PATH_TO_TEST_DATA + "/simple_inv.gfa")
        paths = gfa.to_paths
        dagify = DAGify(paths)
        profile, rep_count = dagify.generate_profiles_with_minimizing_replications()
        graph = SlicedGraph.load_from_slices(dagify.to_slices(profile), paths)
        self.assertEqual(graph, [['CAAATAAG', {x,y}], ['AC', {x}, 'AC', {y}], ['G', {x, y}]])




class GFATest(unittest.TestCase):
    """ test class of gfa.py
    """

    @unittest.skipIf(not os.path.isfile(location_of_xg), "XG binary is not found.")
    def test_gfa(self):
        self.maxDiff = None
        graph = GFA.load_from_gfa(PATH_TO_TEST_DATA + "/test.gfa")
        graph.save_as_xg(PATH_TO_TEST_DATA + "/test.xg", location_of_xg)
        graph2 = GFA.load_from_xg(PATH_TO_TEST_DATA + "/test.xg", location_of_xg)
        self.assertFalse(self.is_different(graph.gfa, graph2.gfa))
#        self.assertEqual(len(graph.gfa.to_gfa1_s().split("\n")), len(graph2.gfa.to_gfa1_s().split("\n")))

    def test_load_gfa_to_graph(self):
        graph, gfa = self.make_graph_from_gfa()
        self.assertEqual(len(graph.paths), 3)
        self.assertEqual(len(graph.nodes), 15)

    def test_gfa_to_sliced_graph(self):
        graph, gfa = self.make_graph_from_gfa()
        slices = SlicedGraph.from_graph(graph)
        x = 'x'
        y = 'y'
        z = 'z'
        print(slices)
        self.assertEqual(slices, [['CAAATAAG', {x, y, z}], ['A', {y, z}, 'G', {x}], ['C', {x, y, z}], ['TTG', {x, y, z}], ['A', {z}, 'G', {x, y}], ['AAATTTTCTGGAGTTCTAT', {x, y, z}], ['T', {x, y, z}], ['ATAT', {x, y, z}], ['T', {x, y, z}], ['CCAACTCTCTG', {x, y, z}]])

    def test_gfa_to_sliced_graph_via_dagify(self):
        #TODO: this is currently close but not quite there.
        # Slices must be fully defined in SlicedGraph.compute_slices()
        graph, gfa = self.make_graph_from_gfa()
        slices = SlicedGraph.from_graph(graph)
        x = 'x'
        y = 'y'
        z = 'z'
        print(slices)
        self.assertEqual(slices, [['CAAATAAG', {x, y, z}], ['A', {y, z}, 'G', {x}], ['C', {x, y, z}], ['TTG', {x, y, z}], ['A', {z}, 'G', {x, y}], ['AAATTTTCTGGAGTTCTAT', {x, y, z}], ['T', {x, y, z}], ['ATAT', {x, y, z}], ['T', {x, y, z}], ['CCAACTCTCTG', {x, y, z}]])

    def make_graph_from_gfa(self):
        gfa = GFA.load_from_gfa(PATH_TO_TEST_DATA + "test.gfa")
        graph = gfa.to_graph
        return graph, gfa

    def test_export_as_gfa(self):
        graph, gfa = self.make_graph_from_gfa()
        new_gfa = GFA.from_graph(graph)
        self.assertFalse(self.is_different(gfa.gfa, new_gfa.gfa))

    def test_load_gfa_to_graph_2(self):
        gfa = GFA.load_from_gfa(PATH_TO_TEST_DATA + "test2.gfa")
        graph = gfa.to_graph
        self.assertIsNotNone(graph)

    @unittest.skipIf(not os.path.isfile(location_of_xg), "XG binary is not found.")
    def test_load_gfa_via_xg(self):
        graph = GFA.load_from_gfa(PATH_TO_TEST_DATA + "/test.gfa")
        graph.save_as_xg(PATH_TO_TEST_DATA + "/test.xg", location_of_xg)
        graph2 = GFA.load_from_xg(PATH_TO_TEST_DATA + "/test.xg", location_of_xg)
        graph = graph2.to_graph
        graph = SlicedGraph.from_graph(graph)
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

        
