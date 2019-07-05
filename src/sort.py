from collections import defaultdict
import networkx as nx


class DAGify:
    def __init__(self, graph: nx.MultiDiGraph):
        self.graph = graph

    def remove_feedback_arc(self):
        """

        :return: nothing
        """
        graph = self.graph.copy()
        sources = []
        sinks = []

        while graph.nodes:
            out_degrees = list(graph.out_degree(graph.nodes()))
            in_degrees = list(graph.in_degree(graph.nodes()))

            for out_edge in out_degrees:
                if out_edge[1] == 0:
                    sinks.append(out_edge[0])
                    graph.remove_node(out_edge[0])
                    out_degrees.remove(out_edge)
                    for in_edge in in_degrees:
                        if in_edge[0] == out_edge[0]:
                            in_degrees.remove(in_edge)

            for in_edge in in_degrees:
                if in_edge[1] == 0:
                    sources.append(in_edge[0])
                    graph.remove_node(in_edge[0])
                    in_degrees.remove(in_edge)
                    for out_edge in out_degrees:
                        if out_edge[0] == in_edge[0]:
                            out_degrees.remove(out_edge)

            delta_max = 0
            delete_index = -1

            for out_edge in out_degree:
                for in_edge in in_degree:
                    if out_edge[0] == in_edge[0]:
                        a = out_edge[0]
                        b = (out_edge[1] - in_edge[1])
                        if b >= delta_max:
                            delta_max = b
                            delete_index = b

            for node in graph.nodes:
                if graph.out_degree(node) - graph.in_degree(node) > delta_max:
                    delta_max = graph.out_degree(node) - graph.in_degree(node)
                    delete_index = node

            if delete_index != -1:
                sources.append(delete_index)
                graph.remove_node(delete_index)

        sinks.reverse()
        sources.extend(sinks)

        remove_list = []

        for edge in self.graph.edges():
            source = edge[0]
            sink = edge[1]
            if sources.index(source) >= sources.index(sink):
                remove_list.append(edge)

        self.graph.remove(remove_list)




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
