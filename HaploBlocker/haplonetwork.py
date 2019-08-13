from typing import List
import numpy as np
from collections import defaultdict
from copy import copy

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
        upstream: {dict((key, value) for key, value in self.upstream.items())}
        downstream: {dict((key, value) for key, value in self.downstream.items())}
        {len(self.specimens)} specimens: {self.specimens}"""

    def is_nothing(self):
        """Useful in Node class definition to check for NOTHING_NODE"""
        return self.ident == -1 and self.start.snp is None and self.end.snp is None

    def validate(self):
        if not self.specimens:
            assert self.specimens, "Specimens are empty" + self.details()
        for n in self.upstream:
            for node, weight in n.downstream.items():
                if not node.is_nothing() and weight < 0:
                    print(n.details())
                    assert weight > -1, node.details()
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
    """Only transition values for nodes already listed in upstream and downstream will be calculated."""
    if node is not NOTHING_NODE:
        update_stream_transitions(node, 'upstream')
        update_stream_transitions(node, 'downstream')

    return node


def update_stream_transitions(node, stream):
    """This will updated either upstream or downstream transition counts based on the
    the value of 'stream'.  This is a meta-programming function that requires the exact
    name of the class field 'upstream' or 'downstream' to work."""
    g = getattr  #
    running = g(node, stream).keys()
    setattr(node, stream, defaultdict(lambda: 0))
    for n in running:
        if n is not NOTHING_NODE:
            g(node, stream)[n] = len(node.specimens.intersection(n.specimens))
    accounted_upstream = sum(g(node, stream).values()) - g(node, stream)[NOTHING_NODE]
    g(node, stream)[NOTHING_NODE] = len(node.specimens) - accounted_upstream
    assert all([count > -1 for count in g(node, stream).values()]), node.details()
    # Cleans up old keys including NOTHING_NODE
    to_be_deleted = {key for key, count in g(node, stream).items() if count == 0}
    for key in to_be_deleted:
        g(node, stream).pop(key, None)


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
                    if parent is not NOTHING_NODE:
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
    my_specimens = copy(anchor.specimens)
    if prev_node is not NOTHING_NODE:  # normal case
        my_specimens = my_specimens.intersection(prev_node.specimens)
    if next_node is not NOTHING_NODE:  # normal case
        my_specimens = my_specimens.intersection(next_node.specimens)
    if prev_node is NOTHING_NODE and next_node is NOTHING_NODE:  # exceptional: both are nothing node
        my_specimens = copy(anchor.specimens)
        # TODO: why are we removing all specimens that transition to nothing?
        for n in anchor.downstream.keys():
            if n is NOTHING_NODE:  # remove dead leads
                my_specimens -= n.specimens
        for n in anchor.upstream.keys():
            if n is NOTHING_NODE:  # remove dead leads
                my_specimens -= n.specimens

    my_start, my_end = prev_node.start, next_node.end
    my_upstream, my_downstream = copy(prev_node.upstream), copy(next_node.downstream)
    if NOTHING_NODE is prev_node:  # Rare case
        my_start = anchor.start
        my_upstream = copy(anchor.upstream)
    if NOTHING_NODE is next_node:  # Rare case
        my_end = anchor.end
        my_downstream = copy(anchor.downstream)

    # TODO: what about case where more content is joining downstream?
    new = Node(777, my_start, my_end, my_specimens, my_upstream, my_downstream)

    # Update Specimens in prev_node, anchor, next_node
    anchor.specimens -= new.specimens
    prev_node.specimens -= new.specimens
    next_node.specimens -= new.specimens

    # Update upstream/downstream
    update_neighbor_pointers(new)
    suspects = {new, prev_node, anchor, next_node}.union(set(new.upstream.keys()), set(new.downstream.keys()))
    for n in suspects:
        update_transition(n)
    new.validate()
    #TODO: Delete nodes with zero specimens from the Graph?
    return new


def update_neighbor_pointers(new_node):
    """Ensure that my new upstream pointers have matching downstream pointers in neighbors,
    and vice versa.  This does not set correct transition rates, it only makes the nodes connected."""
    for n in new_node.upstream.keys():
        if n is not NOTHING_NODE:
            n.downstream[new_node] = 1
    for n in new_node.downstream.keys():
        if n is not NOTHING_NODE:
            n.upstream[new_node] = 1


def split_groups(all_nodes: List[Node]):
    """This is called crossmerge in the R code"""
    length = len(all_nodes)  # size of global_nodes changes, necessitating this weird loop
    for n in range(length):
        node = all_nodes[n]
        # check if all transition upstream match with one of my downstream nodes
        # if set(node.upstream.values()) == set(node.downstream.values()): WHY?
        if len(node.specimens) > 0:
            # Matchup upstream and downstream with specimen identities
            for up in tuple(node.upstream.keys()):
                for down in tuple(node.downstream.keys()):
                    set1 = copy(up.specimens)
                    set2 = copy(down.specimens)
                    if up is NOTHING_NODE:
                        set1 = copy(node.specimens)
                        for index in tuple(node.upstream.keys()):
                            if index is not NOTHING_NODE:
                                set1.difference_update(index.specimens)
                    if down is NOTHING_NODE:
                        set2 = copy(node.specimens)
                        for index in tuple(node.downstream.keys()):
                            if index is not NOTHING_NODE:
                                set2.difference_update(index.specimens)

                    if set1 == set2 and len(set1) > 0:
                        new_node = split_one_group(up, node, down)
                        all_nodes.append(new_node)

    filtered = neglect_nodes(all_nodes, 0)
    return filtered

