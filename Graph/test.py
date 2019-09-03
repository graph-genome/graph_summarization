import unittest

from django.test import TestCase
import os
from os.path import join
from Graph.gfa import GFA
from Graph.models import Node, Path, NodeTraversal, ZoomLevel
from Graph.sort import DAGify

# Define the working directory
from Graph.utils import build_graph_from_slices
from vgbrowser.settings import BASE_DIR
PATH_TO_TEST_DATA = join(BASE_DIR, "test_data")
location_of_xg = join(BASE_DIR, "test_data","xg")


a, b, c, d, e = 'a', 'b', 'c', 'd', 'e'  # Paths must be created first
x, y, z = 'x', 'y', 'z'


class GraphTest(TestCase):
    """Constructing a node with an existing Path object will modify that Path object (doubly linked)
    which means care must be taken when constructing Graphs.  From factory_input we have an example of
    pure python to Graph.build in one step.  In example_graph, we must first declare the Paths,
    then reference them in order in Node Constructors.  Order matters for Graph identity!"""

    def test_example_graph(self):
        example = GFA.load_from_gfa(os.path.join(PATH_TO_TEST_DATA, 'factory_input.gfa'))
        return example.to_graph()

    def test_equalities(self):
        self.assertEqual(Node('A', {}),Node('A', {}))
        self.assertEqual(Node('A', {Path('x')}),Node('A', {Path('x')}))
        self.assertEqual(Node('A', {Path('x'),Path('y')}),Node('A', {Path('x'), Path('y')}))
        # self.assertEqual(SlicedGraph.build([['ACGT', {a, b, c, d}]]), SlicedGraph.build([['ACGT', {a, b, c, d}]]))

    def test_graph_factory(self):
        original_test = [['ACGT', {a, b, c, d}],
                        ['C',   {a, b, d}, 'T', {c}],  # SNP
                        ['GGA', {a, b, c, d}],  # anchor
                        ['C',   {a, b, d}],  # [3] repeated from [1] SNP
                        ['AGTACG', {a, b, c}, 'CGTACT', {d}],  # [4] different membership from [3]
                        ['TTG',    {a, b, c, d}],  # [5] anchor
                        ['A', {a, b}, 'C', {d, e}, 'T', {c}],  # [6] third allele
                        ['GG', {a, b}, 'TT', {c, d}],  # [7] equal size nodes
                        ['C', {a, b, c, e}, 'T', {d}],  # [8] path slip
                        ['C', {a, b, e}, 'T', {c, d}],  # [9] path slip
                        ['C', {a, b, c}, 'T', {d}],  # [10]path slip
                        ['TATA', {a, b, c, d}]]  # [11] anchor
        g1, g2 = build_graph_from_slices(original_test), build_graph_from_slices(original_test)
        assert g1 == g2, \
            ('\n' + repr(g1) + '\n' + repr(g2))
        g_from_GFA = self.test_example_graph()  # comes from matching
        assert g1 == g_from_GFA, repr(g1) + '\n' + repr(g_from_GFA)
        # TODO: stricter equality checks: currently just counting nodes

    def test_summary_storage(self):
        graph = self.test_example_graph()

        zoom1 = ZoomLevel.objects.create(graph=graph, zoom=1, blank_layer=False)
        path1 = zoom1.paths.filter(accession='a').first()
        new_node = Node.objects.create(seq='ACGTCGGA', name='2*2', zoom=zoom1)
        base_nodes = graph.nucleotide_level.nodes
        base_nodes.filter(name__in=['1', '2', '4']).update(summarized_by=new_node)
        assert new_node.children.count() == 3
        zoom1.nodes.filter(name__in=['1', '2', '4']).delete()  # delete copies made redundant in next layer
        for node in base_nodes.filter(name__in=['1','2','4']):
            assert node.summarized_by_id == new_node.id
        print(list(zoom1.nodes.all()))
        self.assertEqual(NodeTraversal.objects.get(order=0, path=path1).downstream().downstream().node.name, '4')
        self.assertEqual(graph.zoomlevel_set.count(), 2)
        self.assertTrue(path1.summary_child),  # "Path should be linked to its child."
        self.assertEqual(zoom1.paths.count(), 5)
        # ZoomLevel
        zoom0 = graph.nucleotide_level
        self.assertEqual(len(zoom1), len(zoom0) - 2)
        self.assertEqual(zoom1.node_ids(),set(range(23, 42)),)#{23,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41}
        self.assertEqual(zoom0.node_ids(), set(range(1,21)))
        names = [x.name for x in zoom1.nodes]
        self.assertEqual(names, ['5', '6', '2*2'])
        sequences = [x.seq for x in zoom1.nodes]
        self.assertEqual(sequences, ['C', 'AGTACG', 'ACGTCGGA'])


    # def test_zoom_level(self):
    #     graph = self.test_example_graph()
    #     node_ids

@unittest.skip  # DAGify has not been converted to databases yet.
class DAGifyTest(TestCase):
    """ test class of sort.py
    """
    # def tearDown(self) -> None:
    #     # Cascade delete all test DB entries
    #     GraphGenome.objects.get_queryset(name__contains=os.path.join(PATH_TO_TEST_DATA)).delete()

    def test_dagify(self):
        gfa = GFA.load_from_gfa(join(PATH_TO_TEST_DATA, "test.gfa"))
        graph = gfa.to_graph()
        dagify = DAGify(graph)
        profile = dagify.generate_profiles(0)

        # self.assertEqual([['CAAATAAG', {x,y,z}], ['A', {y,z}, 'G', {x}], ['C', {x,y,z}], ['TTG', {x,y,z}], ['A', {z}, 'G', {x,y}], ['AAATTTTCTGGAGTTCTAT', {x,y,z}], ['T', {x,y,z}], ['ATAT', {x,y,z}], ['T', {x,y,z}], ['CCAACTCTCTG', {x,y,z}]], graph)

    def test_dagify2(self):
        gfa = GFA.load_from_gfa(join(PATH_TO_TEST_DATA, "test2.gfa"))
        graph = gfa.to_graph()
        dagify = DAGify(graph)
        profile = dagify.generate_profiles(0)
        # graph = SlicedGraph.load_from_slices(dagify.to_slices(profile), paths)
        # x,y,z,a = 'x', 'y', 'z', 'a'
        # self.assertEqual([['CAAATAAG', {x, y, z}], ['G', {x}, 'A', {y, z}], ['C', {x, y}, 'T', {z}], ['TTG', {x, y, z}], ['G', {x, y}, 'A', {a, z}], ['AAATTTTCTGGAGTTCTAT', {a, x, y, z}], ['A', {a, z}, 'T', {x, y}], ['ATAT', {x, y, z}], ['T', {x, y, z}], ['CCAACTCTCTG', {x, y, z}]], graph)

    def test_dagify3(self):
        gfa = GFA.load_from_gfa(join(PATH_TO_TEST_DATA, "test3.gfa"))
        graph = gfa.to_graph()
        dagify = DAGify(graph)
        profile, rep_count = dagify.generate_profiles_with_minimizing_replications()
        self.assertEqual(rep_count, 1)
        # graph = SlicedGraph.load_from_slices(dagify.to_slices(profile), paths)
        # self.assertEqual(graph, [['CAAATAAG', {x, y}], ['CCAACTCTCTG', {y}, 'G', {x}], ['C', {x, y}], ['TTG', {x, y}], ['G', {x, y}], ['AAATTTTCTGGAGTTCTAT', {x, y}], ['T', {x, y}], ['ATAT', {x, y}], ['T', {x, y}], ['CCAACTCTCTG', {x, y}]])

    def test_dagify_altpath(self):
        gfa = GFA.load_from_gfa(join(PATH_TO_TEST_DATA, "alternate_paths.gfa"))
        graph = gfa.to_graph()
        dagify = DAGify(graph)
        profile, rep_count = dagify.generate_profiles_with_minimizing_replications()
        self.assertEqual(rep_count, 1)
        # graph = SlicedGraph.load_from_slices(dagify.to_slices(profile), paths)
        # self.assertEqual(graph, [['CAAATAAG', {x, y}], ['A', {x}, '', {y}], ['G', {x, y}], ['A', {y}, '', {x}], ['T', {x, y}]])

    def test_dagify_dup(self):
        gfa = GFA.load_from_gfa(join(PATH_TO_TEST_DATA, "duplicate.gfa"))
        graph = gfa.to_graph()
        dagify = DAGify(graph)
        profile, rep_count = dagify.generate_profiles_with_minimizing_replications()
        self.assertEqual(rep_count, 2)
        # graph = SlicedGraph.load_from_slices(dagify.to_slices(profile), paths)
        # self.assertEqual(graph, [['CAAATAAG', {x, y}], ['', {x}, 'A', {y}], ['', {x}, 'G', {y}], ['A', {x, y}], ['G', {x, y}], ['T', {x, y}]])


    def test_unresolved_repreat(self):
        gfa = GFA.load_from_gfa(join(PATH_TO_TEST_DATA, "unresolved_repeat.gfa"))
        graph = gfa.to_graph()
        dagify = DAGify(graph)
        profile, rep_count = dagify.generate_profiles_with_minimizing_replications()
        # graph = SlicedGraph.load_from_slices(dagify.to_slices(profile), paths)
        # self.assertEqual([['CAAATAAG', {'x'}, 'T', {'y'}], ['A', {'y', 'x'}], ['G', {'x'}, 'C', {'y'}]], graph)

    @unittest.skip("Inversion is unsupported")
    def test_inversion(self):
        gfa = GFA.load_from_gfa(join(PATH_TO_TEST_DATA, "inversion.gfa"))
        graph = gfa.to_graph()
        dagify = DAGify(graph)
        profile, rep_count = dagify.generate_profiles_with_minimizing_replications()
        # graph = SlicedGraph.load_from_slices(dagify.to_slices(profile), paths)
        # self.assertEqual(graph, [])

    @unittest.skip("Inversion is unsupported")
    def test_nested_inversion(self):
        gfa = GFA.load_from_gfa(join(PATH_TO_TEST_DATA, "nested_inv.gfa"))
        graph = gfa.to_graph()
        dagify = DAGify(graph)
        profile, rep_count = dagify.generate_profiles_with_minimizing_replications()
        # graph = SlicedGraph.load_from_slices(dagify.to_slices(profile), paths)
        # self.assertEqual(graph, [])

    @unittest.skip("Inversion is unsupported")
    def test_simple_inversion(self):
        gfa = GFA.load_from_gfa(join(PATH_TO_TEST_DATA, "simple_inv.gfa"))
        graph = gfa.to_graph()
        dagify = DAGify(graph)
        profile, rep_count = dagify.generate_profiles_with_minimizing_replications()
        # graph = SlicedGraph.load_from_slices(dagify.to_slices(profile), paths)
        # self.assertEqual(graph, [['CAAATAAG', {x,y}], ['AC', {x}, 'AC', {y}], ['G', {x, y}]])



class GFATest(TestCase):
    """ test class of gfa.py
    """

    @unittest.skipIf(not os.path.isfile(location_of_xg), "XG binary is not found.")
    def test_gfa(self):
        self.maxDiff = None
        graph = GFA.load_from_gfa(join(PATH_TO_TEST_DATA, "test.gfa"))
        graph.save_as_xg(join(PATH_TO_TEST_DATA, "test.xg"), location_of_xg)
        graph2 = GFA.load_from_xg(join(PATH_TO_TEST_DATA, "test.xg"), location_of_xg)
        self.assertFalse(self.is_different(graph.gfa, graph2.gfa))
#        self.assertEqual(len(graph.gfa.to_gfa1_s().split("\n")), len(graph2.gfa.to_gfa1_s().split("\n")))

    def test_load_gfa_to_graph(self):
        graph, gfa = self.make_graph_from_gfa()
        self.assertEqual(graph.paths.count(), 3)
        self.assertEqual(graph.nucleotide_level.nodes.count(), 15)

    def make_graph_from_gfa(self):
        gfa = GFA.load_from_gfa(join(PATH_TO_TEST_DATA, "test.gfa"))
        graph = gfa.to_graph()
        return graph, gfa

    def test_export_as_gfa(self):
        graph, gfa = self.make_graph_from_gfa()
        new_gfa = GFA.from_graph(graph)
        self.assertFalse(self.is_different(gfa.gfa, new_gfa.gfa))

    def test_load_gfa_to_graph_2(self):
        gfa = GFA.load_from_gfa(join(PATH_TO_TEST_DATA, "test2.gfa"))
        graph = gfa.to_graph()
        self.assertIsNotNone(graph)

    @unittest.skipIf(not os.path.isfile(location_of_xg), "XG binary is not found.")
    def test_load_gfa_via_xg(self):
        graph = GFA.load_from_gfa(join(PATH_TO_TEST_DATA, "test.gfa"))
        graph.save_as_xg(join(PATH_TO_TEST_DATA, "test.xg"), location_of_xg)
        graph2 = GFA.load_from_xg(join(PATH_TO_TEST_DATA, "test.xg"), location_of_xg)
        graph = graph2.to_graph()
        x,y,z = 'x','y','z'
        self.assertEqual(graph, build_graph_from_slices([['CAAATAAG', {x, y, z}], ['A', {y, z}, 'G', {x}],
                                                         ['C', {x, y, z}], ['TTG', {x, y, z}],
                                                         ['A', {z}, 'G', {x, y}], ['AAATTTTCTGGAGTTCTAT', {x, y, z}], ['T', {x, y, z}],
                                                         ['ATAT', {x, y, z}], ['T', {x, y, z}], ['CCAACTCTCTG', {x, y, z}]]))

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

        
