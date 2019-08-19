from typing import List, Iterable
from django.db import models
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


class GraphGenome(models.Model):
    name = models.CharField(max_length=1000)

    @property
    def paths(self):
        """Getter only.  Shortcut for DB."""
        return self.path_set.all()

    @property
    def nodes(self):
        """Getter only.  Shortcut for DB."""
        return self.node_set.all()

    def __repr__(self):
        """Warning: the representation strings are very sensitive to whitespace"""
        return f"Graph: {self.name}\n{self.path_set.count()} paths  {self.node_set.count()} nodes."

    def __eq__(self, other):
        if isinstance(other, GraphGenome):
            return other.node_set.count() == self.node_set.count() and \
                   other.path_set.count() == self.path_set.count()  # other.name == self.name and \
        return False

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

    def append_node_to_path(self, node_id, strand, path_name) -> None:
        """This is the preferred way to build a graph in a truly non-linear way.
        Nodes will be created if necessary.
        NodeTraversal is appended to Path (order dependent) and PathIndex is added to Node
        (order independent)."""
        if node_id not in self.nodes:  # hasn't been created yet, need to retrieve from dictionary of guid
            if isinstance(node_id, str):
                self.nodes[node_id] = Node('', [], node_id)
            else:
                raise ValueError("Provide the id of the node, not", node_id)
        Path.objects.get(name=path_name).append_node(Node.objects.get(name=node_id), strand)


class Node(models.Model):
    seq = models.CharField(max_length=255, blank=True)
    name = models.CharField(primary_key=True, max_length=15)
    graph = models.ForeignKey(GraphGenome, on_delete=models.CASCADE)

    class Meta:
        unique_together = ['graph', 'name']

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


class Path(models.Model):
    """Paths represent the linear order of on particular individual (accession) as its genome
    was sequenced.  A path visits a series of nodes and the ordered concatenation of the node
    sequences is the accession's genome.  Create Paths first from accession names, then append
    them to Nodes to link together."""
    accession = models.CharField(max_length=1000)  # one path per accession
    graph = models.ForeignKey(GraphGenome, on_delete=models.CASCADE)

    class Meta:
        unique_together = ['graph', 'accession']

    def __getitem__(self, path_index):
        return self.nodes[path_index]

    def __repr__(self):
        """Warning: the representation strings are very sensitive to whitespace"""
        return "'" + self.accession + "'"

    def __eq__(self, other):
        return self.accession == other.accession

    def __hash__(self):
        return hash(self.accession)

    @property
    def nodes(self) -> Iterable['NodeTraversal']:
        return NodeTraversal.objects.filter(path=self).order_by('order').all()

    def append_gfa_nodes(self, nodes):
        assert hasattr(nodes[0], 'orient') and hasattr(nodes[0], 'name'), 'Expecting gfapy.Gfa.path'
        for node in nodes:
            NodeTraversal(node=Node.objects.get(name=node.name),
                          path=self, strand=node.orient).save()

    def append_node(self, node: Node, strand: str):
        """This is the preferred way to build a graph in a truly non-linear way.
        NodeTraversal is appended to Path (order dependent) and PathIndex is added to Node (order independent)."""
        NodeTraversal(node=node, path=self, strand=strand).save()

    # @classmethod
    # def build(cls, name: str, seq_of_nodes: List[str]):
    #     node = Node.objects.create(seq)
    #     for p in paths:
    #         NodeTraversal.objects.create(node, path)

    def name(self):
        return self.accession

    def to_gfa(self):
        return '\t'.join(['P', self.accession, "+,".join([x.node.name + x.strand for x in self.nodes]) + "+", ",".join(['*' for x in self.nodes])])




class NodeTraversal(models.Model):
    """Link from a Path to a Node it is currently traversing.  Includes strand"""
    node = models.ForeignKey(Node, db_index=True, on_delete=models.CASCADE)
    path = models.ForeignKey(Path, db_index=True, on_delete=models.CASCADE, help_text='')
    strand = models.CharField(choices=[('+', '+'),('-', '-')], default='+', max_length=1)
    order = models.IntegerField(help_text='Defines the order a path lists traversals')  # set automatically

    objects = CustomSaveManager()

    def __repr__(self):
        if self.strand == '+':
            return self.node.seq
        else:
            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
            return "".join(complement.get(base, base) for base in reversed(self.node.seq))

    def __eq__(self, other):
        return self.node.id == other.node.id and self.strand == other.strand

    def save(self, **kwargs):
        """Checks the largest 'order' value in the current path and increments by 1.
        IMPORTANT NOTE: save() does not get called if you do NodeTraverseal.objects.create
        or get_or_create"""
        if self.order is None:
            last_traversal = self.path.nodetraversal_set.all().order_by('-order').first()
            self.order = 0 if not last_traversal else last_traversal.order + 1
        super(NodeTraversal, self).save(**kwargs)

