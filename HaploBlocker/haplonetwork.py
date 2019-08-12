from typing import List

import numpy as np
from collections import defaultdict
import networkx as nx
import os

BLOCK_SIZE = 20
FILTER_THRESHOLD = 4


def first(iterable):
    return next(iter(iterable))


class Point:
    def __init__(self, snp, bp=0):
        self.snp, self.bp = snp, bp

    @property
    def window(self):
        return self.snp // BLOCK_SIZE


class Node:
    def __init__(self, ident, start, end, specimens=None, upstream=None, downstream=None):
        self.ident = ident
        self.start = start #Point()
        self.end = end #Point()
        self.specimens = set() if specimens is None else specimens
        # {NOTHING_NODE:501, Node: 38,  Node: 201, Node: 3}
        self.upstream = defaultdict(lambda: 0) if not upstream else upstream
        # {Node: 38,  Node: 201, Node: 3}
        self.downstream = defaultdict(lambda: 0) if not downstream else downstream
        assert self.start is not None and self.end is not None, self.details()
        assert self.end.snp is not None or (self.end.snp is None and self.start.snp is None), self.details()

    def __len__(self):
        return len(self.specimens)

    def __repr__(self):
        return "N%s(%s, %s)" % (str(self.ident), str(self.start.snp), str(self.end.snp))

    def __hash__(self):
        return hash(self.ident + 1) * hash(self.start.snp) * hash(self.end.snp)

    def details(self):
        return f"""Node{self.ident}: {self.start.snp} - {self.end.snp}
        specimens: {self.specimens}"""
        # upstream: {dict((key, value) for key, value in self.upstream.items())}
        # downstream: {dict((key, value) for key, value in self.downstream.items())}

    def is_nothing(self):
        """Useful in Node class definition to check for NOTHING_NODE"""
        return self.ident == -1 and self.start.snp is None and self.end.snp is None

    def validate(self):
        assert self.specimens, "Specimens are empty" + self.details()
        for n in self.upstream:
            for node, weight in n.downstream.items():
                if not node.is_nothing():
                    assert weight > 0, n.details()
        return True

    def is_beginning(self) -> bool:
        return self.start.snp == 0

    def is_end(self) -> bool:
        return len(self.downstream) == 1 and first(self.downstream).is_nothing()


NOTHING_NODE = Node(-1, Point(None), Point(None))


def read_data(file_path):
    """Individuals are rows, not columns"""
    loci = []
    with open(file_path) as ke:
        for line in ke.readlines():
            loci.append(tuple(int(x) for x in line.split()))

    individuals = np.array(loci).T.tolist()
    return loci, individuals


def signature(individual, start_locus):
    return tuple(individual[start_locus: start_locus + BLOCK_SIZE])


def get_unique_signatures(individuals, start_locus):
    unique_blocks = {}
    for individual in individuals:
        sig = signature(individual, start_locus)
        if sig not in unique_blocks:
            unique_blocks[sig] = Node(len(unique_blocks), Point(start_locus // BLOCK_SIZE, start_locus),
                                      Point(start_locus // BLOCK_SIZE, start_locus + BLOCK_SIZE))  # TODO: -1?
    return unique_blocks


def get_all_signatures(alleles, individuals):
    unique_signatures = []
    for locus_start in range(0, len(alleles) - BLOCK_SIZE, BLOCK_SIZE):  # discards remainder
        sig = get_unique_signatures(individuals, locus_start)
        unique_signatures.append(sig)
    return unique_signatures


def build_individuals(individuals, unique_signatures):
    simplified_individuals = []
    for i_specimen, specimen in enumerate(individuals):
        my_simplification = []
        for w, window in enumerate(unique_signatures):  # the length of the genome
            sig = signature(specimen, w * BLOCK_SIZE)
    #         print(sig, unique_signatures[w][sig])
    #         print(i_specimen, window)
            my_simplification.append(unique_signatures[w][sig])
        simplified_individuals.append(my_simplification)
    return simplified_individuals


def populate_transitions(simplified_individuals):
    """ simplified_individuals is a list of loci which contain a list of Nodes which each contain specimen
    build nodes:  [0] first 4 are the 4 starting signatures in window 0.
    Nodes represent a collection of individuals with the same signature at that locus
    For each node list which individuals are present at that node
    List transition rates from one node to all other upstream and downstream
    :param simplified_individuals:
    :return:
    """
    for i, indiv in enumerate(simplified_individuals):
        # look what variants are present
        for x, node in enumerate(indiv):
            node.specimens.add(i)
            if x + 1 < len(indiv):
                node.downstream[indiv[x + 1]] += 1
            else:
                node.downstream[NOTHING_NODE] += 1
            if x - 1 >= 0:
                node.upstream[indiv[x - 1]] += 1
            else:
                node.upstream[NOTHING_NODE] += 1


def update_transition(node):
    if node is not NOTHING_NODE:
        running = node.upstream.keys()
        node.upstream = defaultdict(lambda: 0)
        for n in running:
            if n is not NOTHING_NODE:
                node.upstream[n] = len(node.specimens.intersection(n.specimens))

        running = node.downstream.keys()
        node.downstream = defaultdict(lambda: 0)
        for n in running:
            if n is not NOTHING_NODE:
                node.downstream[n] = len(node.specimens.intersection(n.specimens))

        accounted_upstream = sum(node.upstream.values()) - node.upstream[NOTHING_NODE]
        node.upstream[NOTHING_NODE] = len(node.specimens) - accounted_upstream
        accounted_downstream = sum(node.downstream.values()) - node.downstream[NOTHING_NODE]
        node.downstream[NOTHING_NODE] = len(node.specimens) - accounted_downstream
    return node


def simple_merge(full_graph):
    """ Side effects full_graph by merging any consecutive nodes that have
    identical specimens and removing the redundant node from full_graph.
    :param full_graph:
    :return: full_graph modified
    """
    n = 0
    while n < len(full_graph):  # size of global_nodes changes, necessitating this weird loop
        node = full_graph[n]
    #     print(node, type(node))
        if len(node.downstream) == 1:
            next_node = first(node.downstream.keys())
            if len(node.specimens) == len(next_node.specimens):
                #Torsten deletes nodeA and modifies next_node
                next_node.upstream = node.upstream
                next_node.start = node.start
                #prepare to delete node by removing references
                for parent in node.upstream.keys():
                    if parent != NOTHING_NODE:
                        count = parent.downstream[node]
                        del parent.downstream[node]  # updating pointer
                        parent.downstream[next_node] = count
                full_graph.remove(node)  #delete node
                # zoom_stack[0].append(merged)
                n -= 1
        n += 1
    return full_graph


def delete_node(node, cutoff):
    """Changes references to this node to add to references to NOTHING_NODE"""
    if cutoff < 1:
        return  # if cutoff is 0, then don't touch upstream and downstream
    for parent, count in node.upstream.items():
        parent.downstream[NOTHING_NODE] += parent.downstream[node]
        del parent.downstream[node]
    for descendant, count in node.downstream.items():
        descendant.upstream[NOTHING_NODE] += descendant.upstream[node]
        del descendant.upstream[node]


def neglect_nodes(all_nodes, deletion_cutoff=FILTER_THRESHOLD):
    nodes_to_delete = set()
    #     filtered_nodes = copy(all_nodes)
    #     filtered_nodes.remove(1)
    #     assert len(all_nodes) != len(filtered_nodes)
    for node in all_nodes:
        if len(node.specimens) <= deletion_cutoff:
            delete_node(node, deletion_cutoff)  # TODO: check if this will orphan
            nodes_to_delete.add(node)
    filtered_nodes = [x for x in all_nodes if x not in nodes_to_delete]
    # TODO: remove orphaned haplotypes in a node that transition to and from zero within a 10 window length
    return filtered_nodes




def split_one_group(prev_node, anchor, next_node):
    """ Called when up.specimens == down.specimens"""
    # Comment: That is actually the case we want to split up to obtain longer blocks later
    # Extension of full windows will take care of potential loss of information later
    my_specimens = anchor.specimens
    if prev_node is not NOTHING_NODE:  # normal case
        my_specimens = my_specimens.intersection(prev_node.specimens)
    if next_node is not NOTHING_NODE:  # normal case
        my_specimens = my_specimens.intersection(next_node.specimens)
    if prev_node is NOTHING_NODE and next_node is NOTHING_NODE:  # exceptional: both are nothing node
        my_specimens = anchor.specimens
        for n in anchor.downstream.keys():
            if n is not NOTHING_NODE:  # don't remove empty set
                my_specimens -= n.specimens
        for n in anchor.upstream.keys():
            if n is not NOTHING_NODE:  # don't remove empty set
                my_specimens -= n.specimens

    my_start, my_end = prev_node.start, next_node.end
    my_upstream, my_downstream = prev_node.upstream, next_node.downstream
    if NOTHING_NODE is prev_node:  # Rare case
        my_start = anchor.start
        my_upstream = anchor.upstream
    if NOTHING_NODE is next_node:  # Rare case
        my_end = anchor.end
        my_downstream = anchor.downstream

    # TODO: what about case where more content is joining downstream?
    new = Node(777, my_start, my_end, my_specimens, my_upstream, my_downstream)

    print(new.details())
    print(new.upstream.keys())
    print(new.upstream.values())
    print(sum(new.upstream.values()))

    # Update upstream/downstream
    running = new.upstream.keys()

    ## n.upstream/downstream contains the same key multiple times?!
    ## My quick fix was to delete all upstream/downstream and just recalculate everything...
    new.upstream = defaultdict(lambda: 0)
    for n in running:
        if n != NOTHING_NODE:
            new.upstream[n] = len(new.specimens.intersection(n.specimens))
            n.downstream[new] = new.upstream[n]
            n.downstream[prev_node] = n.downstream[prev_node] - n.downstream[new]
            if n.downstream[prev_node] == 0:
                del n.downstream[prev_node]

    running = new.downstream.keys()
    new.downstream = defaultdict(lambda: 0)
    for n in running:
        if n != NOTHING_NODE:
            new.downstream[n] = len(new.specimens.intersection(n.specimens))
            n.upstream[new] = new.downstream[n]
            n.upstream[next_node] = n.upstream[next_node] - n.upstream[new]
            if n.upstream[next_node] == 0:
                del n.upstream[next_node]

    print(new.details())
    print(new.upstream.keys())
    print(new.upstream.values())
    print(sum(new.upstream.values()))

    accounted_upstream = sum(new.upstream.values()) - new.upstream[NOTHING_NODE]
    # print(f'upstream {sum(new.upstream.values())} downstream {sum(new.downstream.values())}')
    new.upstream[NOTHING_NODE] = len(new.specimens) - accounted_upstream
    accounted_downstream = sum(new.downstream.values()) - new.downstream[NOTHING_NODE]
    new.downstream[NOTHING_NODE] = len(new.specimens) - accounted_downstream

    assert all([count > -1 for count in new.upstream.values()]), new.details()
    assert all([count > -1 for count in new.downstream.values()]), new.details()
    # Update Specimens in prev_node, anchor, next_node
    anchor.specimens -= new.specimens
    if prev_node != NOTHING_NODE:
        prev_node.specimens -= new.specimens
        update_transition(prev_node)

    if next_node != NOTHING_NODE:
        next_node.specimens -= new.specimens
        update_transition(next_node)

    new.validate()
    return new


def split_groups(all_nodes: List[Node]):
    """This is called crossmerge in the R code"""
    length = len(all_nodes)  # size of global_nodes changes, necessitating this weird loop
    for n in range(length):
        node = all_nodes[n]
        # check if all transitition upstream match with one of my downstream nodes
        # if set(node.upstream.values()) == set(node.downstream.values()): WHY?
        if len(node.specimens) > 0:
            # Matchup upstream and downstream with specimen identities
            for up in tuple(node.upstream.keys()):
                for down in tuple(node.downstream.keys()):

                    set1 = up.specimens
                    set2 = down.specimens
                    if up == NOTHING_NODE:
                        set1 = node.specimens
                        for index in tuple(node.upstream.keys()):
                            set1.intersection(index.specimens)  # =- does not work for empty sets
                    if down == NOTHING_NODE:
                        set2 = node.specimens
                        for index in tuple(node.downstream.keys()):
                            set2.intersection(index.specimens)  # =- does not work for empty sets

                    if set1 == set2 and len(set1) > 0:
                        new_node = split_one_group(up, node, down)
                        all_nodes.append(new_node)

    filtered = neglect_nodes(all_nodes, 0)
    return filtered

