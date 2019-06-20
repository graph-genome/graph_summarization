from typing import List, NamedTuple
import pickle

class PathIndex(NamedTuple):
    node: Node
    index: int

class Graph:
    version = 1.0
    def __init__(self, nodes: list(Nodes), paths: list(Paths)):
        self.nodes = nodes
        self.paths = paths

    def load_from_pickle(self, file: str):
        self = pickle.load(file)

    def load_form_xg(self, file: str, xg_bin: str):
        raise NotImplementedError()

    def dump():
        raise NotImplementedError()

    def save_as_pickle():
        pickle.dump(self, file)

    def save_as_xg():
        raise NotImplementedError()

class Node:
    def __init__(self, id: int, sequence: str, paths: list(PathIndex)):
        self.id = id
        self.sequence = sequence
        self.paths = paths
        #self.strand = strand
    
    def alternative_nodes(self):
        #WIP

class Path:
    def __init__(self, name: str, nodes: list(Node)):
        self.name = name
        self.nodes = paths

if __name__ == "__main__":
    location_of_xg = sys.argv[0]

    ### Usage
    graph = Graph.load_form_xg
    graph.save_as_pickle()

