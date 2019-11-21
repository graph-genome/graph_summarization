from typing import List, Iterable, Set, Optional
from django.db import models
from django.db.models import QuerySet, Max, Min

from Utils.models import CustomSaveManager


# GraphGenome specific error classes for more informative error catching
class NoAnchorError(ValueError):
    pass
class PathOverlapError(ValueError):
    pass
class NoOverlapError(PathOverlapError):
    pass
class NodeMissingError(ValueError):
    pass


class ZoomLevelManager(models.Manager):
    def create(self, graph, zoom, blank_layer=False) -> 'ZoomLevel':
        """Creates a full collection of Paths for the new ZoomLevel that matches the
        previous level's path names, but are now blank."""
        me = super(ZoomLevelManager, self).create(graph=graph, zoom=zoom)  # now saved to DB
        if zoom == 0 or blank_layer:
            return me  # We've done all the necessary work
        # Copy previous level in entirety
        previous_level = graph.zoomlevel_set.get(zoom=zoom - 1)
        Node.objects.bulk_create([Node(name=n.name, zoom=me) for n in previous_level.nodes_xrange()], 100)
        # TODO: This loop can be sped up by bulk_create and bulk_update
        for path in previous_level.path_xrange():
            name = path.name
            p = Path(accession=name, zoom=me)  # new Path for a new level
            p.save()
            # TODO: this is REALLY SLOW AND WASTEFUL!
            # Makes a full copy of every traversal in the Path so new copies can be edited
            copies = [NodeTraversal(path=p, node=traverse.node, strand=traverse.strand, order=traverse.order)
                      for traverse in path.nodetraversal_set.all()]
            NodeTraversal.objects.bulk_create(copies, 100)
        return me


class ZoomLevel(models.Model):
    """Each Graph Genome has multiple Zoom levels.  A zoom level is a collection of Paths
    at that zoom level.  When a Path is unmodified from a summarization step, it can reuse
    a previous Path without duplicating it.
    Zoom was originally a single integer in Path, but it turned out a helper class for queries
    would result in more robust code and allow Path reuse."""
    graph = models.ForeignKey('GraphGenome', on_delete=models.CASCADE)
    zoom = models.IntegerField(default=0)  # Zoom level starts at 0 for nucleotide level and moves up
    # Nodes can be shared between zoom levels.
    objects = ZoomLevelManager()

    class Meta:
        unique_together = ['graph', 'zoom']  # one ZoomLevel at each integer

    @property
    def paths(self) -> QuerySet:
        return self.path_set

    def __len__(self):
        """The number of unique NodeTraversals in this level"""
        return self.node_set.count()

    def __repr__(self):
        """Warning: the representation strings are very sensitive to whitespace"""
        return f"Graph: {self.graph.name}\n{self.paths.count()} paths  {self.node_set.count()} layers."

    def __eq__(self, other: 'ZoomLevel'):
        if isinstance(other, ZoomLevel):
            return other.node_set.count() == self.node_set.count() and \
                   other.paths.count() == self.paths.count()  # other.name == self.name and \
        return False

    def node_ids(self) -> Set[int]: # TODO: optimize
        # path_ids = self.paths.values_list('id', flat=True)
        # nodes = set(NodeTraversal.objects.filter(path_id__in=path_ids).values_list('node_id', flat=True))
        return set(self.node_set.values_list('id', flat=True))

    @property
    def nodes(self) -> QuerySet:
        """Getter only.  Shortcut for DB."""
        return self.node_set.all()

    def node(self, node_name):
        return Node.objects.get(name=node_name, zoom=self)

    def path_xrange(self):
        """Lazy evaluation of Paths to ensure that we don't overrun SQL"""
        return self.paths.all().iterator(chunk_size=100)

    def nodes_xrange(self):
        """Lazy evaluation of Nodes to ensure that we don't overrun SQL"""
        return self.paths.all().iterator(chunk_size=100)


class GraphManager(models.Manager):
    """custom create() methods for consistency and convenience"""
    def create(self, **kwargs):
        """When a graph is created.  It should also automatically create the first ZoomLevel."""
        graph = super(GraphManager, self).create(**kwargs)
        ZoomLevel.objects.create(graph=graph, zoom=0)  # nucleotide_level created
        return graph


class GraphGenome(models.Model):
    """Graph Genomes in general are defined as a collection of unordered Nodes which contain sequence,
    and one Path per individual.  A Path visits Nodes in any order, on either strand.  This directed
    graph then contains the relationships of all individuals to each other through shared sequence.
    GraphGenomes in this database contain an extra concept not found in GFA or VG: Summarization.
    GraphGenome has multiple ZoomLevels which each contain a full set of Paths.  Paths at higher
    zoom levels are shorter and visit summary nodes that explain or discard trivial variation to
    allow user to focus on larger graph structure."""
    name = models.CharField(max_length=1000)
    # Zoom level is defined in Paths only.
    # Nodes can be shared between zoom levels.

    objects = GraphManager()

    @property
    def nucleotide_level(self) -> 'ZoomLevel':
        return self.zoomlevel_set.filter(zoom=0).first()

    @property
    def paths(self) -> QuerySet:
        """Getter only.  Shortcut for DB."""
        return self.nucleotide_level.paths

    # @property
    # def nodes(self):
    #     """Getter only.  Shortcut for DB."""
    #     return self.nucleotide_level.nodes

    def __eq__(self, other: 'GraphGenome'):
        if isinstance(other, GraphGenome):
            if other.zoomlevel_set.count() == self.zoomlevel_set.count():
                for a,b in zip(other.zoomlevel_set.all(),self.zoomlevel_set.all()):
                    if a != b:
                        return False
                return True
            return False
        return False

    def __repr__(self):
        """Warning: the representation strings are very sensitive to whitespace"""
        return f"Graph: {self.name}\n{self.paths.count()} paths  {self.nucleotide_level.node_set.count()} nodes at base in {self.zoomlevel_set.count()} layers."

    @classmethod
    def load_from_xg(cls, file: str, xg_bin: str) -> 'GraphGenome':
        """XG is a graph format used by VG (variation graph).  This method builds a
        database GraphGenome to exactly mirror the contents of an XG file."""
        from Graph.gfa import GFA
        gfa = GFA.load_from_xg(file, xg_bin)
        return gfa.to_graph()

    def save_as_xg(self, file: str, xg_bin: str):
        """XG is a graph format used by VG (variation graph).  This method exports
        a database GraphGenome as an XG file."""
        from Graph.gfa import GFA
        gfa = GFA.from_graph(self)
        gfa.save_as_xg(file, xg_bin)


    def highest_zoom_level(self):
        return self.zoomlevel_set.all().order_by('-zoom').first().zoom


class Node(models.Model):
    """Nodes are the fundamental content carriers for sequence.  A Path
    can traverse nodes in any order, on both + and - strands.
    Nodes can be utilized by Paths in multiple zoom levels.
    Only first order Nodes (zoom=0) have sequences, but all have names which can be used
    to fetch the node from GraphGenome.node()."""
    seq = models.CharField(max_length=255, blank=True)
    name = models.CharField(max_length=15)
    zoom = models.ForeignKey(ZoomLevel, on_delete=models.CASCADE)
    summarized_by = models.ForeignKey('Node', null=True, blank=True, on_delete=models.SET_NULL,
                                      related_name='children')
    # Break points for haploblocks - Erik Garrison - service for coordinates
    # Start and stop positions for a node

    class Meta:
        unique_together = ['zoom', 'name']

    def __len__(self):
        return self.nodetraversal_set.count()

    # def __repr__(self):
    #     """Paths representation is sorted because set ordering is not guaranteed."""
    #     return repr(self.seq) + \
    #     ', {' + ', '.join(str(i) for i in list(self.paths)) + '}'

    # def __eq__(self, other):
    #     if not isinstance(other, Node):
    #         print("Warn: comparing Node and ", type(other), other)
    #         return False
    #     return self.seq == other.seq and self.paths == other.paths  # and self.id == other.id

    def __hash__(self):
        return (hash(self.seq) + 1) * hash(self.name)

    def to_gfa(self, segment_id: int):
        return '\t'.join(['S', str(segment_id), self.seq])

    def specimens(self) -> Set[int]:
        return set(self.traverses().values_list('path_id', flat=True))

    def traverses(self) -> QuerySet:
        return self.nodetraversal_set.all()

    def upstream_ids(self) -> Set[int]:
        """Returns the node ids that are upstream of this node."""
        traverses = self.traverses()  #.value_list('node_id', flat=True)
        # Node.objects.filter(id__in=traverses).values_list('id', flat=True)
        return set(t.upstream_id() for t in traverses)

    def upstream(self) -> Set['Node']:
        """Returns the node ids that are upstream of this node."""
        traverses = self.traverses()  #.value_list('node_id', flat=True)
        nodes = set()
        for t in traverses:
            n = t.upstream()  # upstream may be None
            if n:
                nodes.add(n.node)
        return nodes

    def downstream_ids(self) -> Set[int]:
        """Returns the node ids that are downstream of this node."""
        traverses = self.traverses()
        return set(t.downstream_id() for t in traverses)

    def downstream(self) -> Set['Node']:
        """Returns the node ids that are upstream of this node."""
        traverses = self.traverses()  #.value_list('node_id', flat=True)
        return set(t.downstream().node for t in traverses if t.downstream())

    def __repr__(self):
        return "N%s(%s)" % (str(self.name), self.seq)

    def details(self, zoom=0):
        return f"""Node{self.name}: {self.seq}
        upstream: {self.upstream_ids()}
        downstream: {self.downstream_ids()}
        {len(self.specimens())} specimens: {self.specimens()}"""

    def is_nothing(self):
        """Useful in Node class definition to check for Node.NOTHING"""
        return self.name == '-1'

    def validate(self):
        """Returns true if the Node has specimens and does not have any negative
        transition values, raises an AssertionError otherwise."""
        if not self.nodetraversal_set.count():
            assert False, "Orphaned: No NodeTraversals are referencing this Node." + self.details()
        return True


# Node.NOTHING is essential "Not Applicable" when used to track transition rates between nodes.
# Node.NOTHING is an important concept to Haploblocker, used to track upstream and downstream
# that transitions to an unknown or untracked state.  As neglect_nodes removes minority
# allele nodes, there will be specimens downstream that "come from" Node.NOTHING, meaning their
# full history is no longer tracked.  Node.NOTHING is a regular exception case for missing data,
# the ends of chromosomes, and the gaps between haplotype blocks.
Node.NOTHING = Node(-1)


class Path(models.Model):
    """Paths represent the linear order of on particular individual (accession) as its genome
    was sequenced.  A path visits a series of nodes and the ordered concatenation of the node
    sequences is the accession's genome.  Create Paths first from accession names, then append
    them to Nodes to link together."""
    accession = models.CharField(max_length=1000)  # one path per accession
    zoom = models.ForeignKey(ZoomLevel, on_delete=models.CASCADE)

    class Meta:  # we need a database check where each accession name only occurs once per zoom level
        unique_together = ['accession', 'zoom']

    def __getitem__(self, path_index):
        return self.nodes.filter(order=path_index)

    def __repr__(self):
        """Warning: the representation strings are very sensitive to whitespace"""
        return "'" + self.accession + "'"

    def __eq__(self, other):
        return self.accession == other.accession

    def __hash__(self):
        return hash(self.accession)

    @property
    def nodes(self) -> QuerySet:
        return NodeTraversal.objects.filter(path=self).order_by('order')

    @property
    def graph(self) -> GraphGenome:
        return self.zoom.graph

    @property
    def name(self):
        return self.accession

    @property
    def summary_child(self):
        if self.zoom.zoom == 0:
            return None
        previous_zoom = self.zoom.graph.zoomlevel_set.get(zoom=self.zoom.zoom - 1)
        return previous_zoom.path_set.get(accession=self.accession)

    def append_gfa_nodes(self, nodes):
        assert hasattr(nodes[0], 'orient') and hasattr(nodes[0], 'name'), 'Expecting gfapy.Gfa.path'
        for node in nodes:
            # TODO: could be sped up with bulk_create after checking and setting order
            NodeTraversal.objects.create(node=Node.objects.get(name=node.name),
                                         path=self, strand=node.orient)

    def append_node(self, node: Node, strand: str):
        """This is the preferred way to build a graph in a truly non-linear way.
        NodeTraversal.order is guaranteed to be contiguous, making upstream() and downstream() work properly.
        NodeTraversal is appended to Path (order dependent) and PathIndex is added to Node (order independent)."""
        NodeTraversal.objects.create(node=node, path=self, strand=strand)  # calculates order


    def to_gfa(self):
        return '\t'.join(['P', self.accession, "+,".join([x.node.name + x.strand for x in self.nodes]) + "+", ",".join(['*' for x in self.nodes])])


class NodeTraversalManager(models.Manager):
    def create(self, node, path, strand, order=None):
        """Checks the largest 'order' value in the current path and increments by 1.
        IMPORTANT NOTE: save() does not get called if you do NodeTraverseal.objects.create
        or get_or_create"""
        if order is None:
            last_traversal = path.nodetraversal_set.all().order_by('-order').first()
            order = 0 if not last_traversal else last_traversal.order + 1
        return super(NodeTraversalManager, self).create(node=node, path=path, strand=strand, order=order)


class NodeTraversal(models.Model):
    """Link from a Path to a Node it is currently traversing.  Includes strand"""
    node = models.ForeignKey(Node, db_index=True, on_delete=models.CASCADE)
    path = models.ForeignKey(Path, db_index=True, on_delete=models.CASCADE, help_text='')
    strand = models.CharField(choices=[('+', '+'),('-', '-')], default='+', max_length=1)
    # order is set automatically in the CustomSaveManager
    order = models.IntegerField(help_text='Defines the order a path lists traversals. '
                                          'The scale of order is not preserved between zoom levels.')

    objects = NodeTraversalManager()

    class Meta:
        unique_together = ['path', 'order']

    def __repr__(self):
        if self.strand == '+':
            return self.node.seq
        else:
            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
            return "".join(complement.get(base, base) for base in reversed(self.node.seq))

    def __eq__(self, other):
        return self.node.id == other.node.id and self.strand == other.strand

    def fetch_neighbor(self, target_index):
        query = NodeTraversal.objects.filter \
            (path=self.path, order=target_index).values_list('node__id', flat=True)
        if query:
            return query[0]
        return -1

    def upstream_id(self):
        target_index = self.order - 1
        return self.neighbor_id(target_index, True)

    def downstream_id(self):
        target_index = self.order + 1
        return self.neighbor_id(target_index, False)

    def neighbor_id(self, target_index: int, less_than: bool):
        """Faster query that just retruns the node_id, not the NodeTraversal."""
        try:  # This query version can tolerate non-contiguous Path sequences
            if less_than:
                return NodeTraversal.objects.\
                    filter(path=self.path, order__lte=target_index).\
                    order_by('-order').values_list('node_id', flat=True).first()
            return NodeTraversal.objects.\
                filter(path=self.path, order__gte=target_index).\
                order_by('order').values_list('node_id', flat=True).first()
        except (NodeTraversal.DoesNotExist, IndexError):
            return None

    def neighbor(self, target_index: int, less_than: bool) -> Optional['NodeTraversal']:
        try:  # This query version can tolerate non-contiguous Path sequences
            if less_than:
                return NodeTraversal.objects. \
                    filter(path=self.path, order__lte=target_index). \
                    order_by('-order').first()
            return NodeTraversal.objects. \
                filter(path=self.path, order__gte=target_index). \
                order_by('order').first()
        except (NodeTraversal.DoesNotExist, IndexError):
            return None

    def upstream(self) -> 'NodeTraversal':
        """Slower queries that return the neighboring NodeTraversal.  This can be chained."""
        return self.neighbor(self.order - 1, True)

    def downstream(self) -> 'NodeTraversal':
        """Slower queries that return the neighboring NodeTraversal.  This can be chained."""
        return self.neighbor(self.order + 1, False)
