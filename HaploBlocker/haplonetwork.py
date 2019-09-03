"""
HaploBlocker handles Graph Summarization explained in:
Graph Summarization: https://github.com/graph-genome/vgbrowser/issues/3
HaploBlocker: https://github.com/graph-genome/vgbrowser/issues/19
"""
from typing import List, Iterable
import numpy as np
from collections import defaultdict
from copy import copy
from Graph.models import Node, Path, ZoomLevel, NodeTraversal

BLOCK_SIZE = 20
FILTER_THRESHOLD = 4


def first(iterable):
    return next(iter(iterable))

def read_data(file_path):
    """Reads one of Torsten's SNP files.  In the file, Individuals are columns, not rows.
    This method returns both loci (all 501 individual alleles at one loci) and the Transposed
    individuals (all 32,000 loci for one individual at a time)."""
    loci = []
    with open(file_path) as ke:
        for line in ke.readlines():
            loci.append(tuple(int(x) for x in line.split()))

    individuals = np.array(loci).T.tolist()
    return loci, individuals


def signature(individual, start_locus):
    return tuple(individual[start_locus: start_locus + BLOCK_SIZE])


def nodes_from_unique_signatures(individuals, start_locus, current_graph):
    """A signature is a series of BLOCK_SIZE SNPs inside of a locus.  We want to know how many
    unique signatures are present inside of one locus.  A Node is created for each unique
    signature found.
    EX: Signature(1000001011013100200)

    This method does not currently have error tolerance like R Haploblocker"""
    unique_blocks = {}
    for individual in individuals:
        sig = signature(individual, start_locus)
        if sig not in unique_blocks:
            unique_blocks[sig] = Node.objects.create(  # saves to Database
                name=f'{len(unique_blocks)}:{start_locus // BLOCK_SIZE}-{start_locus // BLOCK_SIZE}',
                                      seq=''.join(str(x) for x in sig),
                                      graph=current_graph)
    return unique_blocks


def build_all_slices(alleles, individuals, current_graph):
    """Each item in this list is a slice, representing all the possible states for one locus.
    Inside a slice is a set of Nodes, one for each unique 'signature' or sequence state.
    Paths that all have the same state in this slice all reference the same Node object."""
    slices = []
    for slice_start in range(0, len(alleles) - BLOCK_SIZE, BLOCK_SIZE):  # discards remainder
        nodes = nodes_from_unique_signatures(individuals, slice_start, current_graph)
        slices.append(nodes)
    return slices


def build_paths(individuals, unique_signatures, graph):
    """Describes an individual as a Path (list of Nodes) that the individual visits (NodeTraversals).
    accessions is a list of loci which contain a list of Nodes which each contain specimen
    build nodes:  [0] first 4 are the 4 starting signatures in window 0.
    Nodes represent a collection of individuals with the same signature at that locus
    For each node list which individuals are present at that node
    :param graph: """
    # TODO: It may be more performant to merge build_all_slices and build_paths so that global lists are never stored
    print(f"Building paths from {len(individuals)} individuals and {len(unique_signatures)} loci")
    accessions = []
    for i_specimen, specimen in enumerate(individuals):
        my_path = Path.objects.create(accession=str(i_specimen), graph=graph)
        my_sigs = [unique_signatures[w][signature(specimen, w * BLOCK_SIZE)] for w in range(len(unique_signatures))]
        traverses = [NodeTraversal(node=sig, path=my_path, strand='+', order=i) for i, sig in enumerate(my_sigs)]
        NodeTraversal.objects.bulk_create(traverses, 100)
        accessions.append(my_path)
    print(f"Done building {len(accessions)}Paths")
    return accessions


def populate_transitions(simplified_individuals):
    """
    List transition rates from one node to all other upstream and downstream.
    This method populates Node.specimens and begins the process of side-effecting Nodes.
    To rebuild a fresh Graph copy, you must start at build_all_slices()
    :param simplified_individuals:
    """
    for i, indiv in enumerate(simplified_individuals):
        # look what variants are present
        for x, node in enumerate(indiv):
            node.specimens.add(i)
            if x + 1 < len(indiv):
                node.downstream[indiv[x + 1]] += 1
            # else:
            #     node.downstream[Node.NOTHING] += 1
            if x - 1 >= 0:
                node.upstream[indiv[x - 1]] += 1
            # else:
            #     node.upstream[Node.NOTHING] += 1


def update_transition(node):
    """Only transition values for nodes already listed in upstream and downstream will be calculated."""
    if not node.is_nothing():
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
        if not n.is_nothing():
            g(node, stream)[n] = len(node.specimens.intersection(n.specimens))
    accounted_upstream = sum(g(node, stream).values()) - g(node, stream)[Node.NOTHING]
    g(node, stream)[Node.NOTHING] = len(node.specimens) - accounted_upstream
    assert all([count > -1 for count in g(node, stream).values()]), node.details()
    # Cleans up old keys including Node.NOTHING
    to_be_deleted = {key for key, count in g(node, stream).items() if count == 0}
    for key in to_be_deleted:
        g(node, stream).pop(key, None)


def simple_merge(current_level: ZoomLevel) -> ZoomLevel:
    """ Side effects full_graph by merging any consecutive nodes that have
    identical specimens and removing the redundant my_node from full_graph.
    :param current_level: Graph that will be read and edited
    :return: full_graph modified
    """
    #TODO: Paths start fully populated with redundant NodeTraversals.  Editing NodeTraversals,
    # moves to newly created Nodes.  Global bool for whether or not a particular path was modified.
    next_level, zoom = prep_next_summary_layer(current_level)
    # TODO: Iterate an optimized query or Remove this nonsense comment
    # Node.objects.filter()
    # NodeTraversal.objects.filter(node_id)
    # traverses = Node.nodetraversal_set.filter(path_zoom=zoom)\
    #     .distinct()\
    #     .filter(count=1)
    # downstream_ids = set(t.downstream_id() for t in traverses)
    # a.path_set == b.path_set
    # Node.objects.filter(path_set == )
    for my_node in current_level.nodes():
        # only one Node Downstream, no matter the number of specimens
        if len(my_node.downstream_ids(zoom)) == 1:
            d = my_node.nodetraversal_set.first().downstream()
            if d:
                next_node = d.node  # fetched from DB
                if my_node.nodetraversal_set.count() == next_node.nodetraversal_set.count():  # Not a complete guarantee...
                    # Torsten deletes my_node and modifies next_node
                    merged_node = Node.objects.create(name=f'{my_node.name}*{next_level.zoom}',
                                                      graph=current_level.graph)
                    for x in [my_node, next_node]:
                        x.summarized_by = merged_node
                        x.save()

                    # edit existing traversals
                    path_ids = next_level.paths.values_list('id', flat=True)
                    NodeTraversal.objects.\
                        filter(node=next_node, path_id__in=path_ids).\
                        update(node_id=merged_node.id)
                    # next_node.nodetraversal_set.filter(zoom=zoom).bulk_update(node_id=merged_node.id)

                    # delete my_node and all associates
                    query = NodeTraversal.objects.filter(node=my_node, path_id__in=path_ids)
                    query._raw_delete(query.db)  # https://www.nickang.com/fastest-delete-django/
                    # TODO: merged_node.start = my_node.start
    return next_level


def prep_next_summary_layer(current_level):
    zoom = current_level.zoom
    assert current_level.graph.highest_zoom_level() == zoom, \
        "You should only be summarizing the topmost layer"
    next_level = ZoomLevel.objects.create(graph=current_level.graph, zoom=zoom + 1)
    return next_level, zoom


def delete_node(node: Node, cutoff: int, layer: ZoomLevel):
    """Changes references to this node to add to references to Node.NOTHING"""
    if cutoff < 1:
        return  # if cutoff is 0, then don't touch upstream and downstream
    node.traverses(layer).delete()
    node.validate()


def neglect_nodes(zoom_level : ZoomLevel, deletion_cutoff=FILTER_THRESHOLD):
    """Deletes nodes if they have too few specimens supporting them defined by
    :param deletion_cutoff
    :returns a new list of nodes lacking the pruned nodes in all_nodes"""

    # next_level, zoom = prep_next_summary_layer(current_level)

    for node in zoom_level.nodes():  # TODO optimize distinct count
        if len(node.specimens(zoom_level)) <= deletion_cutoff:
            delete_node(node, deletion_cutoff, zoom_level)


def split_one_group(prev_node, anchor, next_node, zoom_level: ZoomLevel):
    """ Called when up.specimens == down.specimens
    Comment: That is actually the case we want to split up to obtain longer blocks later
    Extension of full windows will take care of potential loss of information later"""

    my_specimens = anchor.specimens(zoom_level.zoom)  # list of path_ids
    my_specimens = my_specimens.intersection(prev_node.specimens(zoom_level.zoom))
    my_specimens = my_specimens.intersection(next_node.specimens(zoom_level.zoom))
    new_node = Node.objects.create(graph=zoom_level.graph, name=f'{anchor.name}:{zoom_level.zoom}')
    for a in (prev_node, anchor, next_node):
        a.summarized_by = new_node
        a.save()

    NodeTraversal.objects.filter(path_id__in=my_specimens, node_id=anchor.id).update(node_id=new_node.id)
    NodeTraversal.objects.filter(path_id__in=my_specimens, node_id=prev_node.id).delete()
    NodeTraversal.objects.filter(path_id__in=my_specimens, node_id=next_node.id).delete()
    # TODO: if this is slow use query._raw_delete

    new_node.validate()
    return new_node


def update_neighbor_pointers(new_node):
    """Ensure that my new upstream pointers have matching downstream pointers in neighbors,
    and vice versa.  This does not set correct transition rates, it only makes the nodes connected."""
    for n in new_node.upstream.keys():
        if not n.is_nothing():
            n.downstream[new_node] = 1
    for n in new_node.downstream.keys():
        if not n.is_nothing():
            n.upstream[new_node] = 1


def split_groups(zoom_level: ZoomLevel):
    """When two haplotypes have a locus in common with no variation, then the graph represents
    this with a single anchor node flanked by 2 haplotypes on either side.  This means 5 Nodes
    are present where 2 would suffice.  Split groups splits the anchor node and gives pieces
    to each haplotype; reducing 5 Nodes to 2.

    :Returns new graph with less nodes, all_nodes is still modified, but length doesn't change
    Note: This is called crossmerge in the R code.
    TODO: Ideally, the database would retain some record of how many nucleotides are shared between
    the two new haplotype nodes."""

    for node in zoom_level.nodes():
        # check if all transition upstream match with one of my downstream nodes
        if len(node.specimens(zoom_level)) > 0:
            # Matchup upstream and downstream with specimen identities
            for up in node.upstream(zoom_level):
                set1 = up.specimens(zoom_level)
                if len(set1):
                    for down in node.downstream(zoom_level):
                        set2 = down.specimens(zoom_level)
                        if set1 == set2 and len(set2) > 0:
                            new_node = split_one_group(up, node, down, zoom_level)
