from collections import defaultdict
from typing import List, NamedTuple
from itertools import tee
import gfapy
import pickle
import subprocess
import io
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
#        self.V = vertices  # No. of vertices

    # function to add an edge to graph
    def add_edge(self, u, v):
        self.graph[u].append(v)
        self.nodes[u] = 1
        self.nodes[v] = 1
        self.V = len(self.nodes.keys())
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

    # topologicalSortUtil()
    def topologicalSort(self):
        # Mark all the vertices as not visited
        visited = [False] * self.V
        stack = []

        # Call the recursive helper function to store Topological
        # Sort starting from all vertices one by one
        for i in range(self.V):
            if visited[i] == False:
                self.topologicalSortUtil(i, visited, stack)

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
            [gfa.add_line(line) for line in stream]
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

    def save_as_gfa(self, file: str):
        self.gfa.to_file(file)

    @property
    def to_graph(self):
        topological_sort_helper = TopologicalSort()
        path_dict = defaultdict(list)
        node_hash = {}

        # Extract all paths into graph
        for path in self.gfa.paths:
            ### TODO: Replace it as a path class object.
            for node in path.segment_names:
                print(node, node.id, path.name)

                path_dict[node].append(path.name)
            for node_pair in pairwise(path.segment_names):
                topological_sort_helper.add_edge(node_pair[0].id, node_pair[1].id)

        # Extract all nodes in the graph.
        for segment in self.gfa.segments:
            print(segment)
            print(path_dict[segment])
            node = Node(segment.sequence, path_dict[segment])
            node_hash[segment.id] = node

        node_stack = topological_sort_helper.topologicalSort()

        print(path_dict)

# Cluster nodes as multiple slices according to the result of topological sort.
        factory_input = []
        current_slice = Slice()
        for node in node_stack:
            # If completely matched with
            if path_dict(node).length == self.gfa.paths:
                factory_input.push(current_slice)
                factory_input.push(node_hash[node])
                current_slice = Slice()
            else:
                if path_dict(node) in [x.paths for x in current_slice.nodes]:
                    factory_input.push(current_slice)
                    current_slice = Slice([node_hash[node]])
                else:
                    current_slice.add_slice(node_hash[node])

        # Convert paths to slices in the graph.
        for node in node_stack:
            # [Slice([Node('ACGT', {1,2,3,4})]),
            #               Slice([Node('C',{1,2,4}),Node('T', {3})]),
            #               Slice([Node('GGA',{1,2,3,4})]),
            #               Slice([Node('C',{1,2,4}),Node('', {3})]),
            #               Slice([Node('AGTACG',{1,2,3}), Node('CGTACT',{4})]),
            #               Slice([Node('TTG',{1,2,3,4})]) ]
            pass

        base_graph = Graph(factory_input)
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
    location_of_xg = sys.argv[0]
