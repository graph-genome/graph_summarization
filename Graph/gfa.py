from collections import defaultdict
from typing import List, NamedTuple
from itertools import tee
import gfapy
import subprocess
import io
import os
import tempfile
from Graph.models import *


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


class TopologicalSort:
    def __init__(self):
        self.graph = defaultdict(list)  # dictionary containing adjacency List
        self.nodes = {}
        self.V = 0

    # function to add an edge to graph
    def add_edge(self, u, v):
        self.graph[u].append(v)
        self.nodes[u] = 1
        self.nodes[v] = 1

    # A recursive function used by topologicalSort
    def topologicalSortUtil(self, v, visited, stack):

        # Mark the current node as visited.
        visited[v] = True

        # Recur for all the vertices adjacent to this vertex
        for i in self.graph[v]:
            if visited[i] == False:
                self.topologicalSortUtil(i, visited, stack)

        # Push current vertex to stack which stores result
        stack.insert(0, v)

    # The function to do Topological Sort. It uses recursive
    def topologicalSort(self):
        # Mark all the vertices as not visited
        # visited = [False] * len(self.nodes.keys())
        stack = []
        visited = defaultdict(lambda: False)

        # Call the recursive helper function to store Topological
        # Sort starting from all vertices one by one
        for i, v in enumerate(self.nodes.keys()):
            if visited[v] == False:
                self.topologicalSortUtil(v, visited, stack)

        # Print contents of stack
        return stack


class GFA:
    def __init__(self, gfa: gfapy.Gfa, source_path: str):
        self.gfa = gfa
        self.source_path = source_path

    #    @classmethod
    #    def load_from_pickle(cls, file: str):
    #        graph = pickle.load(open(file, 'rb'))
    #        return graph

    @classmethod
    def load_from_xg(cls, file: str, xg_bin: str):
        gfa = gfapy.Gfa()
        process = subprocess.Popen([xg_bin, "-i", file, "--gfa-out"], stdout=subprocess.PIPE)
        with io.open(process.stdout.fileno(), closefd=False) as stream:
            [gfa.add_line(line.rstrip()) for line in stream if line != ""]
        process.wait()
        if process.returncode != 0:
            raise OSError()
        graph = cls(gfa, file)
        process.stdout.close()
        return graph

    @classmethod
    def load_from_gfa(cls, file: str):
        gfa = gfapy.Gfa.from_file(file)
        graph = cls(gfa, file)
        return graph

    #    def save_as_pickle(self, outfile: str):
    #        with open(outfile, 'wb') as pickle_file:
    #            pickle.dump(self.gfa, pickle_file, protocol=2)

    def save_as_xg(self, file: str, xg_bin: str):
        with tempfile.NamedTemporaryFile(mode="w+t", delete=False) as f:
            f.write(self.gfa.to_gfa1_s())

        process = subprocess.check_output([xg_bin, "-o", file, "-g", f.name])
        os.remove(f.name)
        return process

    def save_as_gfa(self, file: str):
        self.gfa.to_file(file)

    @classmethod
    def from_graph(cls, graph: Graph):
        """Constructs the lines of a GFA file listing paths, then sequence nodes in arbitrary order."""
        gfa = gfapy.Gfa()
        for path in graph.paths:
            node_series = ",".join([traverse.node.id + traverse.strand for traverse in path.nodes])
            gfa.add_line('\t'.join(['P', path.accession, node_series, ",".join(['*' for _ in path.nodes])]))
        for node in graph.nodes.values(): # in no particular order
            gfa.add_line('\t'.join(['S', str(node.id), node.seq]))
        return cls(gfa, "from Graph")

    def to_paths(self) -> List[Path]:
        # create parent object for this genome
        gdb = GraphGenome.objects.get_or_create(name=self.source_path)[0]
        for segment in self.gfa.segments:
            node_id = segment.name
            Node.objects.get_or_create(seq=segment.sequence, name=node_id, graph=gdb)

        paths = []
        for path in self.gfa.paths:
            p = Path(path.name, graph=gdb).save()
            p.append_gfa_nodes(path.segment_names)
            paths.append(p)
        # path_names = [path.name for path in self.gfa.paths]
        # list(Path.objects.get(name__in=path_names))
        return paths

    def to_graph(self) -> GraphGenome:
        paths = self.to_paths()
        if paths:
            return paths[0].graph
        else:
            return None

        # Extract all paths into graph
        # path_names = [p.name for p in self.gfa.paths]
        # graph = Graph(path_names)  # Paths can be empty at start
        # for path in self.gfa.paths:
        #     for path_index, node in enumerate(path.segment_names):
        #         graph.append_node_to_path(node.name, node.orient, path.name, path_index)
        # for segment in self.gfa.segments:
        #     graph.nodes[segment.name].seq = segment.sequence
        # graph.paths = self.to_paths()
        # return graph


'''
class XGWrapper:
    @staticmethod
    def save(gfa):
        pass
    
    @staticmethod
    def load(gfa):
        pass

class GraphStack:
    def __init__(graphs: List[Graph]):
        self.graphs = graphs
'''

if __name__ == "__main__":
    location_of_xg = sys.argv[0]
