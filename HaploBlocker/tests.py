from django.test import TestCase
import unittest
import os
# Create your tests here.
from HaploBlocker.models import Node, Path, Edge
from HaploBlocker.haplonetwork import read_data, get_all_signatures, build_individuals, get_unique_signatures, \
    populate_transitions, simple_merge, neglect_nodes, split_groups

#
# class ModelTest(TestCase):
#     def test_creation(self):
#         p = Path.objects.create(accession='watermelon')
#         n = Node.objects.create(seq='ACGT')
#         Edge.objects.create(node=n, path=p)
#         Edge.objects.create(node=Node.objects.create(seq='AAAA'), path=p)
#         Edge.objects.create(node=Node.objects.create(seq='TTT'), path=p)
#         print(Edge.objects.all())


class HaploTest(unittest.TestCase):

    def setUp(self) -> None:
        print(os.getcwd())
        self.alleles, self.individuals = read_data("test_data/KE_chromo10.txt")
        self.unique_signatures = get_all_signatures(self.alleles, self.individuals)
        self.simplified_individuals = build_individuals(self.individuals, self.unique_signatures)
        # G = build_graph(simplified_individuals)
        populate_transitions(self.simplified_individuals)


    def test_master(self):
        pass


    def test_read(self):
        assert len(self.alleles) == 32767
        assert len(self.individuals[1]) == 32767
        assert len(self.individuals) == 501


    def test_build_individuals(self):
        self.internal_build_individuals(self.alleles, self.individuals)


    def internal_build_individuals(self, alleles, individuals):
        unique_signatures = get_all_signatures(alleles, individuals)
        assert repr(unique_signatures[
                        21]) == '{(0, 0, 0, 0, 2, 2, 0, 2, 0, 2, 2, 2, 2, 2, 2, 2, 0, 2, 0, 2): N0(21, 21), (0, 0, 2, 0, 2, 2, 0, 2, 0, 2, 2, 2, 2, 2, 2, 2, 0, 2, 0, 2): N1(21, 21), (0, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0): N2(21, 21), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0): N3(21, 21), (0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0): N4(21, 21), (0, 0, 0, 2, 2, 2, 0, 2, 0, 2, 2, 2, 2, 0, 0, 2, 0, 0, 0, 2): N5(21, 21), (0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0): N6(21, 21), (0, 0, 0, 0, 2, 2, 0, 2, 2, 2, 2, 2, 2, 0, 0, 2, 0, 0, 2, 2): N7(21, 21)}'
        simplified_individuals = build_individuals(individuals, unique_signatures)
        assert repr(simplified_individuals[500][
                    :100]) == '[N2(0, 0), N2(1, 1), N2(2, 2), N2(3, 3), N2(4, 4), N2(5, 5), N3(6, 6), N3(7, 7), N3(8, 8), N2(9, 9), N0(10, 10), N1(11, 11), N2(12, 12), N2(13, 13), N2(14, 14), N2(15, 15), N3(16, 16), N3(17, 17), N4(18, 18), N3(19, 19), N5(20, 20), N3(21, 21), N3(22, 22), N10(23, 23), N4(24, 24), N3(25, 25), N4(26, 26), N3(27, 27), N1(28, 28), N1(29, 29), N4(30, 30), N3(31, 31), N21(32, 32), N1(33, 33), N1(34, 34), N1(35, 35), N1(36, 36), N1(37, 37), N1(38, 38), N1(39, 39), N1(40, 40), N1(41, 41), N1(42, 42), N1(43, 43), N1(44, 44), N1(45, 45), N1(46, 46), N1(47, 47), N1(48, 48), N1(49, 49), N1(50, 50), N1(51, 51), N1(52, 52), N1(53, 53), N1(54, 54), N1(55, 55), N1(56, 56), N1(57, 57), N1(58, 58), N1(59, 59), N1(60, 60), N1(61, 61), N1(62, 62), N1(63, 63), N1(64, 64), N1(65, 65), N1(66, 66), N1(67, 67), N1(68, 68), N1(69, 69), N1(70, 70), N1(71, 71), N1(72, 72), N1(73, 73), N1(74, 74), N1(75, 75), N1(76, 76), N1(77, 77), N0(78, 78), N0(79, 79), N1(80, 80), N1(81, 81), N1(82, 82), N1(83, 83), N1(84, 84), N1(85, 85), N1(86, 86), N1(87, 87), N1(88, 88), N1(89, 89), N1(90, 90), N1(91, 91), N1(92, 92), N1(93, 93), N1(94, 94), N1(95, 95), N1(96, 96), N1(97, 97), N0(98, 98), N0(99, 99)]'
        assert len(simplified_individuals) == 501 and len(simplified_individuals[60]) == 1638


    def test_get_unique_signatures(self):
        unique_blocks = get_unique_signatures(self.individuals, 0)
        assert len(unique_blocks) == 4
        assert unique_blocks.__repr__() == '{(0, 2, 0, 0, 2, 0, 2, 0, 0, 0, 2, 0, 0, 0, 2, 2, 0, 0, 0, 0): N0(0, 0), ' \
                                           '(0, 0, 2, 2, 0, 2, 0, 2, 2, 2, 0, 2, 2, 2, 0, 0, 2, 2, 2, 2): N1(0, 0), ' \
                                           '(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0): N2(0, 0), ' \
                                           '(2, 0, 2, 2, 0, 2, 0, 2, 2, 2, 0, 2, 2, 2, 0, 0, 2, 2, 2, 2): N3(0, 0)}'

    @unittest.skip
    def test_no_duplicate_nodes(self):
        unique_nodes = set()
        duplicates_found = 0
        for locus in self.simplified_individuals:
            for node in locus:
                # assert isinstance(node, Node)
                if node in unique_nodes:  # If two nodes have the same __hash__ they'll be "in"
                    print(node.details())
                    duplicates_found += 1
                else:
                    unique_nodes.add(node)
        assert duplicates_found == 0, f"Found {duplicates_found} duplicated nodes in the graph"


    def _test_simple_merge(self, all_nodes):
        assert len(all_nodes) == 7180
        summary1 = simple_merge(all_nodes)
        assert len(summary1) == 3690
        return summary1

    def _test_neglect_nodes(self, all_nodes):
        summary2 = neglect_nodes(all_nodes)
        assert len(summary2) == 2854
        unchanged = neglect_nodes(summary2, 0)
        assert len([n for n in unchanged if len(n.specimens) == 0]) == 0
        return summary2


    def _test_split_groups(self, all_nodes):
        summary3 = split_groups(all_nodes)
        assert summary3
        return summary3


    def test_workflow(self):
        all_nodes = [node for window in self.unique_signatures for node in window.values()]  # think about referencing and deletion
        summary1 = self._test_simple_merge(all_nodes)
        summary2 = self._test_neglect_nodes(summary1)
        summary3 = self._test_split_groups(summary2)

        # test_signatures = get_all_signatures(alleles, individuals)
        # test_individuals = build_individuals(individuals, test_signatures)
        # populate_transitions(test_individuals)  # no return val
        #
        # test1 = test_simple_merge(test_signatures)
        # test2 = neglect_nodes(test1)
        # test3 = split_groups(test2)




