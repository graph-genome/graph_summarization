from typing import List

from django.shortcuts import render

# View contains the endpoints on the server for the browser to fetch data
from Graph.models import GraphGenome, Path, Node


def build_from_test_slices(cmd: List, graph_name='test_data') -> GraphGenome:
    """This factory uses test data shorthand for linear graph slices to build
    a database GraphGenome with all the necessary Paths and Nodes.  Path order populated in the order
    that they are mentioned in the slices.  Currently, this is + only and does not support non-linear
    orderings.  Use Path.append_node() to build non-linear graphs."""
    if isinstance(cmd, str):
        cmd = eval(cmd)
    # preemptively grab all the path names from every odd list entry
    graph = GraphGenome.objects.get_or_create(name=graph_name)[0]  # + str(datetime.now())
    node_count = 0
    paths = {key for sl in cmd for i in range(0, len(sl), 2) for key in sl[i + 1]}
    path_objs = {}
    for path_name in paths:
        path_objs[path_name] = Path.objects.get_or_create(graph=graph, accession=path_name)[0]
    for sl in cmd:
        try:
            for i in range(0, len(sl), 2):
                paths_mentioned = [path_objs[key] for key in sl[i + 1]]
                node, is_new = Node.objects.get_or_create(seq=sl[i], name=graph.name + str(node_count), graph=graph)
                node_count += 1
                for path in paths_mentioned:
                    path.append_node(node, '+')
        except IndexError:
            raise IndexError("Expecting two terms: ", sl[0])  # sl[i:i+2])

    return graph