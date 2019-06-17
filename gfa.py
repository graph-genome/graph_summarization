from typing import List, NamedTuple
import gfapy
import pickle

class Graph:
    def __init__(gfa: gfapy.Gfa):
        self.gfa = gfa
    
    def load_from_pickle(self, file: str):
        self.gfa = pickle.load(file)

    def load_form_xg(self, file: str, xg_bin: str):
        raise NotImplementedError()
    
    def load_from_gfa(self, file: str)
        self.gfa = gfapy.Gfa.from_file(file) 

    def save_as_pickle(self):
        pickle.dump(self.gfa, file)

    def save_as_xg(self):
        raise NotImplementedError()

    def save_as_gfa(self, file:str):
        self.gfa.to_file(file) 

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

if __name__ == "__main__":
    location_of_xg = sys.argv[0]

    ### Usage
    graph = Graph.load_form_xg("test/test.xg", location_of_xg)
    graph.save_as_pickle("test/test.pickle")

