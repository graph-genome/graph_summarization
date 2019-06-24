from typing import Callable, Iterator, Union, Optional, List, Iterable, NamedTuple
from itertools import zip_longest
import pickle
import sys

class Node:
    def __init__(self, seq: str, paths: Iterable[int]):
        assert isinstance(seq, str), seq
        assert not isinstance(paths, str) and isinstance(paths, Iterable), paths
        self.seq = seq
        self.paths = set(paths)

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

def merge(self, smaller: Node) -> Node:
    m = Node(self.seq, self.paths.union(smaller.paths))
    # TODO: penalize paths with nucleotide mismatch
    return m
Node.merge = merge  # Typing is picky about order of declaration


class Slice:
    def __init__(self, nodes: Iterable[Node]):
        self.nodes = set(nodes)

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
        biggest = self.primary()
        return max((x for x in self.nodes if x != biggest), key=len)  # When they're the same size, take the next one

    def smallest(self):
        biggest = self.primary()
        return min((x for x in self.nodes if x != biggest),
                   key=len)  # when they're the same size it will take the last listed

    version = 1.0


class Graph:
    def __init__(self, cmd: List):
        """Factory for generating graphs from a representation"""
        self.slices = []
        if isinstance(cmd, str):
            cmd = eval(cmd)
        for sl in cmd:
            current_slice = []
            if isinstance(sl[0], Node):  # already Nodes, don't need to build
                current_slice = sl
            else:
                try:
                    for i in range(0, len(sl), 2):
                        current_slice.append(Node(sl[i], sl[i + 1]))
                except IndexError:
                    raise IndexError("Expecting two terms: ", sl[0])  # sl[i:i+2])

            self.slices.append(Slice(current_slice))

    def __repr__(self):
        """Warning: the representation strings are very sensitive to whitespace"""
        return self.slices.__repr__()

    def __getitem__(self, i):
        return self.slices[i]

    def __eq__(self, representation):
        if isinstance(representation, Graph):
            return all(slice_a == slice_b for slice_a, slice_b in zip_longest(self.slices, representation.slices))
        return self == Graph(representation)  # build a graph then compare it

    def load_from_pickle(self, file: str):
        self = pickle.load(file)

    def load_form_xg(self, file: str, xg_bin: str):
        raise NotImplementedError()

    def save_as_pickle(self, file):
        pickle.dump(self, file)

    def save_as_xg(self):
        raise NotImplementedError()



class Path:
    """TODO: Paths have not been implemented yet."""
    def __init__(self, name: str, nodes: List[Node]):
        self.name = name
        self.nodes = nodes

# class PathIndex(NamedTuple):
#     node: Node
#     index: int

if __name__ == "__main__":
    location_of_xg = sys.argv[0]

    ### Usage
    graph = Graph.load_form_xg
    graph.save_as_pickle()

