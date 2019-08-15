from django.db import models

# Create your models here.
from django.db import models


class Node(models.Model):
    seq = models.CharField(max_length=1000)


class Path(models.Model):
    accession = models.CharField(unique=True, max_length=1000)


class Edge(models.Model):
    node = models.ForeignKey(Node, on_delete=models.CASCADE)
    path = models.ForeignKey(Path, on_delete=models.CASCADE)

    class Meta:
        ordering = ['path','id']