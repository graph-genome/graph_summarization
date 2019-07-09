from collections import defaultdict
from typing import List, NamedTuple
from itertools import tee
import gfapy
import pickle
import subprocess
import io
from IPython import embed
import os
import tempfile
from src.graph import *


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
    def __init__(self, gfa: gfapy.Gfa):
        self.gfa = gfa

    #    @classmethod
    #    def load_from_pickle(cls, file: str):
    #        graph = pickle.load(open(file, 'rb'))
    #        return graph

    @classmethod
    def load_form_xg(cls, file: str, xg_bin: str):
        gfa = gfapy.Gfa()
        process = subprocess.Popen([xg_bin, "-i", file, "--gfa-out"], stdout=subprocess.PIPE)
        with io.open(process.stdout.fileno(), closefd=False) as stream:
            [gfa.add_line(line.rstrip()) for line in stream if line != ""]
        process.wait()
        if process.returncode != 0:
            raise OSError()
        graph = cls(gfa)
        process.stdout.close()
        return graph

    @classmethod
    def load_from_gfa(cls, file: str):
        gfa = gfapy.Gfa.from_file(file)
        graph = cls(gfa)
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
        gfa = gfapy.Gfa()
        path_list = defaultdict(list)
        segment_id = 0
        for slice in graph.slices:
            for node in slice.nodes:
                segment_id += 1
                gfa.add_line('\t'.join(['S', str(segment_id), node.seq]))
                for path in node.paths:
                    path_list[path].append(segment_id)
        for path_key in path_list:
            path_values = [str(x) for x in path_list[path_key]]
            gfa.add_line('\t'.join(['P', path_key, "+,".join(path_values)+"+", ",".join(['*' for _ in path_values])]))
        return cls(gfa)

    @property
    def to_graph(self):
        topological_sort_helper = TopologicalSort()
        path_dict = defaultdict(list)
        node_hash = {}

        # Extract all paths into graph
        for path in self.gfa.paths:
            for node in path.segment_names:
                path_dict[node.name + node.orient].append(path.name)
            for node_pair in pairwise(path.segment_names):
                topological_sort_helper.add_edge(
                    node_pair[0].name + node_pair[0].orient,
                    node_pair[1].name + node_pair[1].orient)

        # Extract all nodes in the graph.
        for segment in self.gfa.segments:
            node_id = segment.name + "+"
            node = Node(segment.sequence, path_dict[node_id])
            node_hash[node_id] = node

            node_id = segment.name + "-"
            node = Node(segment.sequence, path_dict[node_id])
            node_hash[node_id] = node

        node_stack = topological_sort_helper.topologicalSort()

        # Cluster nodes as multiple slices according to the result of the topological sort.
        factory_input = []
        current_slice = Slice([])
        for node in node_stack:
            if len(path_dict[node]) == len(self.gfa.paths):
                if len(current_slice.nodes) > 0:
                    factory_input.append(current_slice)
                factory_input.append(Slice([node_hash[node]]))
                current_slice = Slice([])
            else:
                all_set = set()
                for items in [x.paths for x in current_slice.nodes]:
                    all_set = all_set | items
                if set(path_dict[node]) & all_set != set():
                    if len(current_slice.nodes) > 0:
                        current_slice.add_node(Node("", set([x.name for x in self.gfa.paths]) - all_set))
                        factory_input.append(current_slice)
                    current_slice = Slice([node_hash[node]])
                else:
                    current_slice.add_node(node_hash[node])

        base_graph = Graph.load_from_slices(factory_input)
        return base_graph


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
    graph = GFA.load_from_gfa(sys.argv[0])
    print(graph)
