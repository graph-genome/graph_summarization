from collections import defaultdict
from typing import List

from Graph.models import GraphGenome, Path, Node


class keydefaultdict(defaultdict):
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError( key )
        else:
            ret = self[key] = self.default_factory(key)
            return ret


def build_graph_from_slices(cmd: List, graph_name='test_data') -> GraphGenome:
    """This factory uses test data shorthand for linear graph slices to build
    a database GraphGenome with all the necessary Paths and Nodes.  Path order populated in the order
    that they are mentioned in the slices.  Currently, this is + only and does not support non-linear
    orderings.  Use Path.append_node() to build non-linear graphs."""
    if isinstance(cmd, str):
        cmd = eval(cmd)
    # preemptive grab all the path names from every odd list entry
    graph = GraphGenome.objects.create(name=graph_name)
    node_count = 0
    paths = {key for sl in cmd for i in range(0, len(sl), 2) for key in sl[i + 1]}
    path_objs = {}
    for path_name in paths:
        if graph.paths.filter(accession=path_name).count() == 0:
            path_objs[path_name] = Path.objects.create(accession=path_name, zoom=graph.nucleotide_level)
    for sl in cmd:
        try:
            for i in range(0, len(sl), 2):
                paths_mentioned = [path_objs[key] for key in sl[i + 1]]
                node, is_new = Node.objects.get_or_create(
                    seq=sl[i],
                    name=graph.name + str(node_count),
                    zoom=graph.nucleotide_level)
                node_count += 1
                for path in paths_mentioned:
                    path.append_node(node, '+')
        except IndexError:
            raise IndexError("Expecting two terms: ", sl[0])  # sl[i:i+2])

    return graph