from django.shortcuts import render

from django.http import HttpResponse



def index(request):
    return HttpResponse()  #"Index of Edges.\n" + str(Edge.objects.all()))

# Create your views here.
