from unittest import skip, skipIf
from django.test import TestCase
import Graph.utils
import Graph.views
from Graph.models import GraphGenome, Path, ZoomLevel, NodeTraversal
from vgbrowser.settings import BASE_DIR
import os
# Create your tests here.
# from HaploBlocker.models import Node, Path, Edge
from HaploBlocker.haplonetwork import Node, split_one_group
from HaploBlocker.haplonetwork import read_data, build_all_slices, build_paths, nodes_from_unique_signatures, \
    populate_transitions, simple_merge, neglect_nodes, split_groups


class ModelTest(TestCase):
    def test_creation(self):
        g = GraphGenome.objects.create(name='test_creation')
        p = Path.objects.create(accession='watermelon', graph=g)
        n = Node.objects.create(seq='ACGT', graph=g)
        NodeTraversal.objects.create(node=n, path=p, strand='+')
        NodeTraversal.objects.create(node=Node.objects.create(seq='AAAA', name='-9', graph=g),
                                     path=p, strand='+')
        NodeTraversal.objects.create(node=Node.objects.create(seq='TTT', name='-8', graph=g),
                                     path=p, strand='+')
        print(NodeTraversal.objects.all())


class HaploTest(TestCase):
    @classmethod
    def setUpClass(self) -> None:
        """Reads the input data file once.  Tests that need a fresh graph must
        call create_graph()"""
        print("Setting up HaploTest:", os.getcwd())
        self.alleles, self.individuals = read_data(os.path.join(BASE_DIR, "test_data/KE_chromo10.txt"))
        print("Finished reading SNP file")
        super(HaploTest, self).setUpClass()

    def create_graph(self, graph_name):
        """Tests that need a fresh graph must call create_graph() FIRST!
        Graph summarization works by side effecting Node objects.  Tests can not run independently
        with order dependent side effects.  This method is slow, so don't use it unless you
        need it.
        :param graph_name: """
        print("Creating test Graph")
        graph = GraphGenome.objects.create(name=graph_name)
        slices = build_all_slices(self.alleles, self.individuals, graph)
        self.paths = build_paths(self.individuals, slices, graph)
        print("Finished test graph.")
        return graph


    def test_read(self):
        assert len(self.alleles) == 32767
        assert len(self.individuals[1]) == 32767
        assert len(self.individuals) == 501


    def test_build_individuals(self):
        self.internal_build_individuals(self.alleles, self.individuals)


    def internal_build_individuals(self, alleles, individuals):
        graph = GraphGenome.objects.create(name='internal_build_individuals')
        unique_signatures = build_all_slices(alleles, individuals, graph)
        peek = repr(list(unique_signatures[21].values()))
        assert peek == '[N0:21-21(00002202022222220202), N1:21-21(00202202022222220202), N2:21-21(00022200000000000000), N3:21-21(00000000000000000000), N4:21-21(00002200000000000000), N5:21-21(00022202022220020002), N6:21-21(02000000000000000000), N7:21-21(00002202222220020022)]', peek
        simplified_individuals = build_paths(individuals, unique_signatures, graph)
        traverses = simplified_individuals[500].nodes.filter(order__lt=100)  # [:100]
        nodes = [t.node for t in traverses.prefetch_related('node').all()]
        peek = repr(nodes)
        self.maxDiff = None  # tells the debugger to show the whole thing
        self.assertEqual(peek, '[N2:0-0(00000000000000000000), N2:1-1(00000000000000000000), N2:2-2(00000000000000000000), N2:3-3(00000000000000000000), N2:4-4(00000000000000000000), N2:5-5(00000000000000000000), N3:6-6(00000000000000000000), N3:7-7(00000000000000000000), N3:8-8(00000000000000000000), N2:9-9(00000000000000000000), N0:10-10(00000000000000000000), N1:11-11(00000000000000000000), N2:12-12(00000000000000000000), N2:13-13(00000000000000000000), N2:14-14(00000000000000000000), N2:15-15(00000000000000000000), N3:16-16(00000000000000000000), N3:17-17(00000000000000000000), N4:18-18(00000000000000000000), N3:19-19(00000000000000000000), N5:20-20(00000000000000000000), N3:21-21(00000000000000000000), N3:22-22(00000000000000000000), N10:23-23(00200000000000000000), N4:24-24(00002200222220002000), N3:25-25(02000222220002020222), N4:26-26(20022000002002220002), N3:27-27(22222202020222220000), N1:28-28(00000000000000000000), N1:29-29(00000000000000000022), N4:30-30(00002222222000002200), N3:31-31(00022222202000000000), N21:32-32(00000020202200022020), N1:33-33(02202220022020222000), N1:34-34(00020002000202222002), N1:35-35(22220002220022200022), N1:36-36(22222200000000000000), N1:37-37(00202002222220000200), N1:38-38(00000200000202022200), N1:39-39(02202000202202220000), N1:40-40(00020222200020000020), N1:41-41(20220020022200022200), N1:42-42(00000000000000000000), N1:43-43(00000000000000000000), N1:44-44(00000000000000000000), N1:45-45(00000000000000000000), N1:46-46(00000002220220020020), N1:47-47(00202220222220222202), N1:48-48(00000000000000000002), N1:49-49(20002200000002220022), N1:50-50(22020002002020202022), N1:51-51(02202222220222202000), N1:52-52(20000020000000000000), N1:53-53(00000000000000000000), N1:54-54(00000000000000000000), N1:55-55(00220220000200000220), N1:56-56(20000000202022022020), N1:57-57(20222022022202222220), N1:58-58(22022202222222020200), N1:59-59(22202200202220202220), N1:60-60(22020022220200022022), N1:61-61(20202220000220000220), N1:62-62(00022002000000000000), N1:63-63(00000220000000000000), N1:64-64(00000000000220200000), N1:65-65(00022020200000020022), N1:66-66(20020222222020200020), N1:67-67(00000000000000202222), N1:68-68(22222222000202222202), N1:69-69(22022222020020000022), N1:70-70(00002002220022222200), N1:71-71(22002020020202000000), N1:72-72(00022202000202220020), N1:73-73(22000000000000200020), N1:74-74(22220222220200202202), N1:75-75(00022202222200000000), N1:76-76(00000220220200200022), N1:77-77(02200202020020200000), N0:78-78(00002000000000000000), N0:79-79(00000000000000000000), N1:80-80(00000000000022220000), N1:81-81(00000000000000000000), N1:82-82(00022220200202202202), N1:83-83(20202222200202202202), N1:84-84(00000020000000000000), N1:85-85(00222022020000000002), N1:86-86(22020222020222222000), N1:87-87(00022222002020222022), N1:88-88(00002222000000000200), N1:89-89(00000000000000220022), N1:90-90(22020202200020222220), N1:91-91(00002000002220002222), N1:92-92(22200000000000000000), N1:93-93(00000000000000000000), N1:94-94(00202022200202222222), N1:95-95(22222202202020222222), N1:96-96(00222220200202222020), N1:97-97(22002202220222222022), N0:98-98(20222222222222020220), N0:99-99(20222222220222222002)]')
        self.assertEqual(len(simplified_individuals),501)
        self.assertEqual(simplified_individuals[60].nodes.count(), 1638)


    def test_get_unique_signatures(self):
        graph = GraphGenome.objects.create(name='test_get_unique_signatures')
        unique_blocks = nodes_from_unique_signatures(self.individuals, 0, graph)
        assert len(unique_blocks) == 4
        peek = repr(list(unique_blocks.values()))
        assert peek == \
               '[N0:0-0(02002020002000220000), N1:0-0(00220202220222002222), ' \
               'N2:0-0(00000000000000000000), N3:0-0(20220202220222002222)]', peek

    @skip
    def test_no_duplicate_nodes(self):
        graph = self.create_graph('test')
        unique_nodes = set()
        duplicates_found = 0
        for locus in self.paths:
            for node in locus:
                # assert isinstance(node, Node)
                if node in unique_nodes:  # If two nodes have the same __hash__ they'll be "in"
                    print(node.details())
                    duplicates_found += 1
                else:
                    unique_nodes.add(node)
        assert duplicates_found == 0, f"Found {duplicates_found} duplicated nodes in the graph"


    def _test_simple_merge(self, graph: GraphGenome, zoom_level: int) -> ZoomLevel:
        # these tests could be made independent of test_workflow, but it would be slower
        assert graph.highest_zoom_level() == zoom_level
        starting_level = ZoomLevel.objects.get(graph=graph, zoom=zoom_level)
        self.assertEqual(len(starting_level), 7180)
        next_level = simple_merge(starting_level)
        #Test every Path has a representative in this ZoomLevel
        self.assertEqual(Path.objects.filter(graph=graph, zoom=zoom_level + 1).count(),
                         Path.objects.filter(graph=graph, zoom=zoom_level + 0).count())
        self.assertEqual(NodeTraversal.objects.filter(graph=graph, zoom=zoom_level+1).count(), 3690) #*501?
        return next_level

    @skip
    def test_simple_merge(self):
        graph = self.create_graph('test')
        summary1 = self._test_simple_merge(graph, 0)


    def _test_neglect_nodes(self, all_nodes, zoom_level):
        summary2 = neglect_nodes(all_nodes)
        assert len(summary2) == 2854
        unchanged = neglect_nodes(summary2, 0)
        assert len([n for n in unchanged if len(n.specimens) == 0]) == 0
        return summary2


    def test_split_one_group(self):
        """Self contained example test to look at outputs of one group
                     ['C', {a, b, d}, 'T', {c}],  # SNP
                     ['GGA', {a, b, c, d}],  # anchor
                     ['C', {a, b, d}, '', {c}],  # [3] repeated from [1] SNP
        """
        nodes = [
            ['90', {1, 2, 3, 4}], #[5]
            ['91', {1, 2, 4}, '92', {3}],
            ['93', {1, 2, 3, 4}],  # [2] anchor
            ['94', {1, 2, 4}, '95', {3}], #[3] [4]
            ['96', {1, 2, 3, 4}] #[6]
        ]
        g = Graph.utils.build_graph_from_slices(nodes, '9')
        first, anchor, third = g.node('91'), g.node('93'), g.node('94')
        new_node = split_one_group(first, anchor, third)  # no mentions of minorities [1] or [4]
        print(new_node.details())
        assert new_node in g.node('90').downstream and g.node('92') in g.node('90').downstream
        assert g.node('91') not in g.node('90').downstream
        assert g.node('90') in new_node.upstream and g.node('96') in new_node.downstream
        assert new_node in g.node('96').upstream and g.node('95') in g.node('96').upstream
        assert g.node('94') not in g.node('96').upstream

    def _test_split_groups(self, graph, zoom_level):
        summary3 = split_groups(graph)
        assert len(summary3) > 10
        return summary3


    def test_workflow(self):
        graph = self.create_graph('test')
        summary1 = self._test_simple_merge(graph, 0)
        summary2 = self._test_neglect_nodes(graph, 1)
        summary3 = self._test_split_groups(graph, 2)
        assert len(summary1) > len(summary2) > len(summary3), "Each summarization should result in less nodes"
        summary4 = simple_merge(summary3, 3)
        bad = summary3[2]
        print(bad.details())

        # test_signatures = build_all_slices(alleles, individuals)
        # test_individuals = build_paths(individuals, test_signatures)
        # populate_transitions(test_individuals)  # no return val
        #
        # test1 = test_simple_merge(test_signatures)
        # test2 = neglect_nodes(test1)
        # test3 = split_groups(test2)




