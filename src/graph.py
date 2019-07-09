from typing import List, Iterable
from itertools import zip_longest
import pickle
import sys

from src.utils import keydefaultdict


class NoAnchorError(ValueError):
    pass
class PathOverlapError(ValueError):
    pass
class NoOverlapError(PathOverlapError):
    pass
class NodeMissingError(ValueError):
    pass


class Node:
    def __init__(self, seq: str, paths: 'Iterable[Path]'):
        assert isinstance(seq, str), seq
        self.seq = seq
        self.paths = set()  # Set[PathIndex]
        for p in paths:
            self.append_path(p)

    def __len__(self):
        return len(self.paths)

    def __repr__(self):
        """Paths representation is sorted because set ordering is not guaranteed."""
        return repr(self.seq) + \
        ', {' + ', '.join(str(i) for i in sorted(list(self.paths))) + '}'

    def __eq__(self, other):
        if not isinstance(other, Node):
            print("Warn: comparing Node and ", type(other), other)
            return False
        return self.seq == other.seq and self.paths == other.paths

    def __hash__(self):
        return hash(self.seq)

    def append_path(self, path):
        """Instead: Use Path.append_node if possible"""
        assert isinstance(path, Path), path
        self.paths.add(PathIndex(path, len(path.nodes)))  # not parallelizable
        path.nodes.append(NodeTraversal(self))

    def to_gfa(self, segment_id: int):
        return '\t'.join(['S', str(segment_id), self.seq])

    # Typing is picky about order of declaration, but strings bypass this PEP484
    def merge_minor(self, minor_allele: 'Node') -> 'Node':
        m = Node(self.seq, self.paths.union(minor_allele.paths))
        # TODO: penalize paths with nucleotide mismatch
        return m

    def intersection(self, downstream: 'Node') -> 'Node':
        m = Node(self.seq + downstream.seq,
                 self.paths.intersection(downstream.paths))
        return m

    def union(self, downstream: 'Node') -> 'Node':
        return Node(self.seq + downstream.seq,
                 self.paths.union(downstream.paths))

class Slice:
    def __init__(self, nodes: Iterable[Node]):
        self.nodes = set(nodes)

    def add_node(self, node: Node):
        self.nodes.add(node)

    def alternatives(self, main):
        return self.nodes.difference({main})

    def bystanders(self, first, second):
        return self.nodes.difference({first, second})

    def __len__(self):
        return len(self.nodes)

    def __repr__(self):
        # return '{' + ', '.join(str(i) for i in sorted(list(self.nodes))) + '}'
        return list(self.nodes).__repr__()  # '['+ ','.join(self.paths)+']'

    def __eq__(self, other):
        if isinstance(other, Slice):
            # all(a==b for a,b in zip_longest(self.nodes,other.nodes)) # order dependent
            if not self.nodes == other.nodes:
                print(self.nodes, other.nodes, sep='\n')
            return self.nodes == other.nodes
        else:
            print("Warn: comparing Slice and ", type(other), other)
            return False

    def __iter__(self):
        return iter(self.nodes)

    def primary(self):
        return max(self.nodes, key=len)  # When they're the same size, take the other
    biggest = primary  # alias method

    def secondary(self):
        if len(self.nodes) < 2:
            raise NodeMissingError("Secondary requested when there is no alternative", self.nodes)
        biggest = self.primary()
        return max((x for x in self.nodes if x != biggest), key=len)  # When they're the same size, take the next one

    def smallest(self):
        if len(self.nodes) < 2:
            raise NodeMissingError("Smallest node requested when there is no alternative", self.nodes)
        biggest = self.primary()
        return min((x for x in self.nodes if x != biggest),
                   key=len)  # when they're the same size it will take the last listed

    version = 1.0

class Path:
    """Paths represent the linear order of on particular individual (accession) as its genome
    was sequenced.  A path visits a series of nodes and the ordered concatenation of the node
    sequences is the accession's genome.  Create Paths first from accession names, then append
    them to Nodes to link together."""
    def __init__(self, accession: str):
        self.accession = accession  # one path per accessions
        self.nodes = [] # List[NodeTraversal]
        self.position_checkpoints = {}  # TODO: currently not used

    def __getitem__(self, path_index):
        return self.nodes[path_index]

    def __repr__(self):
        """Warning: the representation strings are very sensitive to whitespace"""
        return "'" + self.accession + "'"

    def __eq__(self, other):
        return self.accession == other.accession

    def __hash__(self):
        return hash(self.accession)

    def append_node(self, node: Node, strand: str):
        """This is the preferred way to build a graph in a truly non-linear way.
        NodeTraversal is appended to Path (order dependent) and PathIndex is added to Node (order independent)."""
        self.nodes.append(NodeTraversal(node, strand))
        node.paths.add(PathIndex(self, len(self.nodes)-1))  # already appended node
        return node

    def to_gfa(self):
        return '\t'.join(['P', self.accession, "+,".join([x.node.name + x.strand for x in self.nodes]) + "+", ",".join(['*' for x in self.nodes])])


class PathIndex:
    """Link from a Node to the place in the path where the Node is referenced.  A Node can appear
    in a Path multiple times.  Index indicates which instance it is."""
    def __init__(self, path: Path, index: int):
        self.path = path
        self.index = index

    def __repr__(self):
        return repr(self.path.accession)

    def __eq__(self, other):
        if self.path.accession == other.path.accession and self.index == other.index:
            return True
        else:
            return False

    def __lt__(self, other):
        return self.path.accession < other.path.accession

    def __hash__(self):
        return hash(self.path.accession) * (self.index if self.index else 1)


class NodeTraversal:
    """Link from a Path to a Node it is currently traversing.  Includes strand"""
    def __init__(self, node: Node, strand: str = '+'):
        self.node = node
        self.strand = strand  # TODO: make this required

    def __repr__(self):
        return self.node.seq


class Graph:
    def __init__(self, paths: List = None):
        """Factory for generating graphs from a representation"""
        self.slices = []
        # This can create orphan Nodes with no traversals
        self.nodes = keydefaultdict(lambda key: Node(key, []))  # node id = Node object
        if all(isinstance(x, str) for x in paths):
            self.paths = {x: Path(x) for x in paths}
        elif all(isinstance(x, Path) for x in paths):
            self.paths = {path.name: path for path in paths}
        else:
            self.paths = {}
        #TODO: calculate slices?

    @staticmethod
    def build(cmd):
        """This factory uses existing slice declarations to build a graph with Paths populated in the order
        that they are mentioned in the slices.  Currently, this is + only and does not support non-linear
        orderings.  Use Path.append_node() to build non-linear graphs."""
        path_dict = keydefaultdict(lambda key: Path(key))  # construct blank path if new
        slices = []
        if isinstance(cmd, str):
            cmd = eval(cmd)
        for sl in cmd:
            current_slice = []
            if isinstance(sl, Slice):
                slices.append(sl)
            else:
                if isinstance(sl[0], Node):  # already Nodes, don't need to build
                    current_slice = sl
                else:
                    try:
                        for i in range(0, len(sl), 2):
                            paths = [path_dict[key] for key in sl[i + 1]]
                            current_slice.append(Node(sl[i], paths))
                    except IndexError:
                        raise IndexError("Expecting two terms: ", sl[0])  # sl[i:i+2])

                slices.append(Slice(current_slice))
        return Graph.load_from_slices(slices)

    @classmethod
    def load_from_slices(cls, slices):
        graph = cls([])
        graph.slices = slices
        return graph

    def __repr__(self):
        """Warning: the representation strings are very sensitive to whitespace"""
        return self.slices.__repr__()

    def __getitem__(self, i):
        return self.slices[i]

    def __eq__(self, representation):
        if isinstance(representation, Graph):
            return all(slice_a == slice_b for slice_a, slice_b in zip_longest(self.slices, representation.slices))
        return self == Graph.build(representation)  # build a graph then compare it

    def load_from_pickle(self, file: str):
        self = pickle.load(file)

    def load_from_xg(self, file: str, xg_bin: str):
        from src.gfa import GFA
        gfa = GFA.load_from_xg(file, xg_bin)
        self = gfa.to_graph()

    def save_as_pickle(self, file):
        pickle.dump(self, file)

    def save_as_xg(self, file: str, xg_bin: str):
        from src.gfa import GFA
        gfa = GFA.from_graph(self)
        gfa.save_as_xg(file, xg_bin)

    def append_node_to_path(self, name, strand, path_name):
        """This is the preferred way to build a graph in a truly non-linear way.
        Nodes will be created if necessary.
        NodeTraversal is appended to Path (order dependent) and PathIndex is added to Node
        (order independent)."""
        if name not in self.nodes:  # hasn't been created yet, need to retrieve from dictionary of guid
            if isinstance(name, str):
                self.nodes[name] = Node(name, [])
            else:
                raise ValueError("Provide the id of the node, not", name)
        self.paths[path_name].append_node(self.nodes[name], strand)


if __name__ == "__main__":
    location_of_xg = sys.argv[0]

    ### Usage  # Unfinished
    # graph = Graph.load_from_xg('../test/test.xg', "../test/xg")
    # graph.save_as_pickle()

