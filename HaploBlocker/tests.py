from django.test import TestCase

# Create your tests here.
from HaploBlocker.models import Node, Path, Edge


class ModelTest(TestCase):
    def test_creation(self):
        p = Path.objects.create(accession='watermelon')
        n = Node.objects.create(seq='ACGT')
        Edge.objects.create(node=n, path=p)
        Edge.objects.create(node=Node.objects.create(seq='AAAA'), path=p)
        Edge.objects.create(node=Node.objects.create(seq='TTT'), path=p)
        print(Edge.objects.all())